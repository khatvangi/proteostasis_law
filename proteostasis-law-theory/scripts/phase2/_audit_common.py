"""shared plumbing for the phase 2 audit of phase 1 experiment C.

nothing in this module changes the phase 1 model, scripts, configs or results.
it only reads them back and re-derives quantities independently.

the one rule the audit follows: a phase 2 number is only trustworthy if it was
recomputed from `samples.tsv` (or from the model itself), never copied out of
`summary.json`.  `summary.json` is the object under audit, so it can never be
the source of the audit's own numbers.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "scripts"))

from proteostasis.model import Params  # noqa: E402

#: the phase 1 run under audit. overridable so the audit can be pointed at a
#: re-run without editing code.
DEFAULT_RUN = REPO_ROOT / "results" / "phase1" / "run_20260731T162946-0500"


def runRoot(override: Optional[str] = None) -> Path:
    if override:
        return Path(override)
    return Path(os.environ.get("PHASE1_RUN", str(DEFAULT_RUN)))


def auditRoot(override: Optional[str] = None) -> Path:
    """the phase 2 output directory; created by the caller, not here."""
    if override:
        p = Path(override)
    else:
        p = Path(os.environ["PHASE2_AUDIT_DIR"])
    p.mkdir(parents=True, exist_ok=True)
    return p


def sha256File(path: os.PathLike | str) -> str:
    """streaming sha256. the phase 1 helper is reimplemented here on purpose:
    the audit must not inherit a bug from the code it is auditing."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def loadSamples(run: Path) -> pd.DataFrame:
    df = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    if "error" in df:
        df["error"] = df["error"].fillna("")
    return df


def loadProvenance(run: Path) -> Dict[str, Any]:
    with open(run / "C" / "provenance.json") as fh:
        return json.load(fh)


# ---------------------------------------------------------------------------
# parameter reconstruction
# ---------------------------------------------------------------------------

#: the sampled names, exactly as experiment C's config declares them. `chi` is
#: not a `Params` field: it is expanded into the conserved pools.
SAMPLED = ("nu", "chi", "rho_U", "rho_A", "alpha_n", "alpha_g", "alpha_d", "m",
           "kappa_ref", "kappa_u", "kappa_a", "kappa_dis", "kappa_cu",
           "kappa_ca", "kappa_du", "kappa_da")


def paramsFromRow(row: pd.Series, base_params: Dict[str, float],
                  rescue_total: float = 1.0) -> Params:
    """rebuild the exact `Params` a phase 1 draw used, from the stored p_* columns.

    this reproduces `sweeps.paramsFromSample` but is driven by the persisted
    table rather than by re-drawing the latin hypercube, so a change in the
    sampler cannot silently shift which parameter set is being audited.
    """
    p = Params(**base_params)
    kw = {k: float(row[f"p_{k}"]) for k in SAMPLED
          if k in Params.__dataclass_fields__ and f"p_{k}" in row.index}
    p = p.with_(**kw)
    chi = float(row["p_chi"])
    p = p.with_(c_tot=chi * rescue_total, d_tot=(1.0 - chi) * rescue_total)
    return p.validate()


def checkReconstruction(df: pd.DataFrame, base_params: Dict[str, float],
                        rescue_total: float = 1.0) -> Dict[str, Any]:
    """independent proof that `paramsFromRow` rebuilds the right parameter set.

    experiment C stored `removal_ceiling = c_tot + (rho_U+rho_A)*d_tot` for every
    draw.  recomputing it from the reconstructed `Params` and comparing is a
    cheap end-to-end check on the whole reconstruction path: allocation split,
    rescue total, and the two rho values all have to be right for it to match.
    """
    worst = 0.0
    n = 0
    for _, row in df.iterrows():
        if not np.isfinite(row.get("removal_ceiling", np.nan)):
            continue
        p = paramsFromRow(row, base_params, rescue_total)
        got = p.c_tot + (p.rho_U + p.rho_A) * p.d_tot
        worst = max(worst, abs(got - float(row["removal_ceiling"]))
                    / max(abs(float(row["removal_ceiling"])), 1e-300))
        n += 1
    return dict(n_checked=int(n), max_rel_error=float(worst))


# ---------------------------------------------------------------------------
# small io helpers
# ---------------------------------------------------------------------------


def writeJson(obj: Any, path: os.PathLike | str) -> str:
    def default(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (np.bool_,)):
            return bool(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, Path):
            return str(o)
        raise TypeError(f"not serializable: {type(o)}")
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True, default=default)
        fh.write("\n")
    return sha256File(path)


def writeTable(df: pd.DataFrame, path: os.PathLike | str) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.reindex(sorted(df.columns), axis=1).to_csv(
        path, sep="\t", index=False, float_format="%.12g", lineterminator="\n")
    return sha256File(path)


def parallelMap(fn, items: Sequence[Any], n_workers: int, chunksize: int = 1):
    """order-preserving parallel map; serial when a single worker is requested."""
    items = list(items)
    if n_workers <= 1 or len(items) <= 1:
        return [fn(x) for x in items]
    import multiprocessing as mp
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=n_workers) as pool:
        return pool.map(fn, items, chunksize=chunksize)
