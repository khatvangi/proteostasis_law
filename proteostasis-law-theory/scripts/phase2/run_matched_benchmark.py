#!/usr/bin/env python3
"""matched boron/nitrogen equivalence benchmark -- the 2x2 factorial.

    model form  x  root protocol
    ----------     -------------
    boron (finite sequestration, at each epsilon on the ladder)
    free  (epsilon = 0 face, nitrogen's convention, same binding constants)

    x  P_BORON    dense 9x9 log multistart + Radau attractor confirmation
       P_NITROGEN four fixed linear guesses + DOP853 return check

every cell is scored on the SAME deterministic 7-D latin hypercube
(seed 20260731, n = 2000) with the remaining sixteen coordinates pinned, so a
cell-to-cell difference can only come from the two factors.

the four cells named in the audit:

    cell 1  boron, epsilon = 1e-6, P_NITROGEN   |  1 vs 2 = code equivalence
    cell 2  free,  epsilon = 1e-6, P_NITROGEN   |
    cell 3  boron, epsilon = 1e-6, P_BORON      |  1 vs 3 = solver effect alone
    cell 4  boron, epsilon = 2.0,  P_BORON      |  3 -> 4 = model form alone

admissibility is NOT folded into the label.  every attractor is scored under all
four defensible criteria -- (burden u+a < H) and (componentwise u < H, a < H),
each in TOTAL and in FREE coordinates -- and all four are written to the row.
boron's `burden_threshold` and nitrogen's `threshold_u`/`threshold_a` are both
1.0, so H = 1 is shared; the criteria still differ in shape.

this script produces no scientific conclusion.  it produces labelled rows.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import platform
import subprocess
import sys
import time
from multiprocessing import Pool
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
import scipy

_HERE = Path(__file__).resolve().parent
if str(_HERE.parent) not in sys.path:
    sys.path.insert(0, str(_HERE.parent))

from phase2 import protocols                                        # noqa: E402
from phase2.lhs import (DEFAULT_N, DEFAULT_SEED, describeDesign,     # noqa: E402
                        parametersForEpsilon, sampleMatrix, sampleHash)
from phase2.mapping import EPSILON_LADDER, boundSubstrateFraction    # noqa: E402
from phase2.models import makeAdapter                                # noqa: E402

THRESHOLD_H = 1.0

ROW_FIELDS = [
    "cell", "model_form", "epsilon", "protocol", "sample", "label", "status",
    "j_u", "n_load", "nucleation", "growth", "ref_capacity",
    "deg_u_capacity", "deg_a_capacity",
    "u_eq", "a_eq", "u_free", "a_free", "residual",
    "eig_real_max", "eig_real_min", "eig_imag_absmax",
    "n_roots", "n_stable",
    "adm_burden_total", "adm_comp_total", "adm_burden_free", "adm_comp_free",
    "continuation_evals", "seconds",
]

# globals for the worker processes; set once by `_initWorker`
_G: Dict = {}


def cellName(model_form: str, epsilon: float, protocol: str) -> str:
    return f"{model_form}_eps{epsilon:g}_{protocol}"


def _admissibility(u: float, a: float, uf: float, af: float, h: float = THRESHOLD_H):
    return {
        "adm_burden_total": bool(u + a < h),
        "adm_comp_total": bool(u < h and a < h),
        "adm_burden_free": bool(uf + af < h),
        "adm_comp_free": bool(uf < h and af < h),
    }


def _initWorker(matrix, model_form, epsilon, protocol, t_end):
    _G.update(matrix=matrix, model_form=model_form, epsilon=epsilon,
              protocol=protocol, t_end=t_end)


def _one(i: int) -> Dict:
    """classify sample i in the cell this worker was initialised for."""
    t0 = time.perf_counter()
    z = _G["matrix"][i]
    q = parametersForEpsilon(z, _G["epsilon"])
    row = {
        "cell": cellName(_G["model_form"], _G["epsilon"], _G["protocol"]),
        "model_form": _G["model_form"], "epsilon": _G["epsilon"],
        "protocol": _G["protocol"], "sample": i,
        "j_u": q.j_u, "n_load": q.n_load, "nucleation": q.nucleation,
        "growth": q.growth, "ref_capacity": q.ref_capacity,
        "deg_u_capacity": q.deg_u_capacity, "deg_a_capacity": q.deg_a_capacity,
    }
    try:
        adapter = makeAdapter(_G["model_form"], q, _G["epsilon"])
        if _G["protocol"] == "boron":
            res = protocols.classifyBoron(adapter, t_end=_G["t_end"])
        else:
            res = protocols.classifyNitrogen(adapter)
    except Exception as exc:                                        # noqa: BLE001
        row.update(label=protocols.LABEL_FAIL,
                   status=f"{type(exc).__name__}: {exc}",
                   u_eq=np.nan, a_eq=np.nan, u_free=np.nan, a_free=np.nan,
                   residual=np.nan, eig_real_max=np.nan, eig_real_min=np.nan,
                   eig_imag_absmax=np.nan, n_roots=0, n_stable=0,
                   adm_burden_total=False, adm_comp_total=False,
                   adm_burden_free=False, adm_comp_free=False,
                   continuation_evals=0, seconds=time.perf_counter() - t0)
        return row

    r = res.get("root")
    row.update(label=res["label"], status=res["status"],
               n_roots=res["n_roots"], n_stable=res["n_stable"])
    if r is None:
        row.update(u_eq=np.nan, a_eq=np.nan, u_free=np.nan, a_free=np.nan,
                   residual=np.nan, eig_real_max=np.nan, eig_real_min=np.nan,
                   eig_imag_absmax=np.nan, adm_burden_total=False,
                   adm_comp_total=False, adm_burden_free=False, adm_comp_free=False)
    else:
        row.update(u_eq=r["u"], a_eq=r["a"], u_free=r["u_free"], a_free=r["a_free"],
                   residual=r["residual"], eig_real_max=r["eig_real_max"],
                   eig_real_min=r["eig_real_min"], eig_imag_absmax=r["eig_imag_absmax"])
        adm = _admissibility(r["u"], r["a"], r["u_free"], r["a_free"])
        # an attractor that was not confirmed is not admissible under any criterion
        if res["label"] != protocols.LABEL_STABLE:
            adm = {k: False for k in adm}
        row.update(adm)
    # non-zero only for the boron arm under the nitrogen protocol, where hybr's
    # linear-coordinate line search probes negative burdens; see
    # `phase2/boron_continuation.py` for why that is not a code-path violation
    row["continuation_evals"] = int(getattr(adapter, "continuation_evaluations", 0))
    row["seconds"] = time.perf_counter() - t0
    return row


def runCell(matrix: np.ndarray, model_form: str, epsilon: float, protocol: str,
            outdir: Path, workers: int, t_end: float) -> Dict:
    """one cell of the factorial -> one TSV plus a counts dict."""
    name = cellName(model_form, epsilon, protocol)
    path = outdir / f"{name}.tsv"
    n = len(matrix)
    t0 = time.time()
    args = (matrix, model_form, epsilon, protocol, t_end)
    if workers > 1:
        with Pool(workers, initializer=_initWorker, initargs=args) as pool:
            rows = list(pool.imap(_one, range(n), chunksize=4))   # order preserving
    else:
        _initWorker(*args)
        rows = [_one(i) for i in range(n)]

    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=ROW_FIELDS, delimiter="\t")
        w.writeheader()
        for row in rows:
            w.writerow(row)

    counts: Dict[str, int] = {}
    for row in rows:
        counts[row["label"]] = counts.get(row["label"], 0) + 1
    adm_counts = {
        k: int(sum(bool(r[k]) for r in rows))
        for k in ("adm_burden_total", "adm_comp_total", "adm_burden_free", "adm_comp_free")
    }
    return {
        "cell": name, "model_form": model_form, "epsilon": epsilon,
        "protocol": protocol, "n": n, "labels": counts,
        "admissible_counts": adm_counts,
        "continuation_evals_total": int(sum(r.get("continuation_evals", 0) for r in rows)),
        "rows_using_continuation": int(sum(1 for r in rows if r.get("continuation_evals", 0))),
        "bound_substrate_fraction_at_low_burden": boundSubstrateFraction(
            parametersForEpsilon(matrix[0], epsilon), epsilon),
        "wall_seconds": time.time() - t0,
        "tsv": path.name,
        "tsv_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def gitHead(repo: Path) -> str:
    try:
        return subprocess.run(["git", "-C", str(repo), "rev-parse", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
    except Exception:                                               # noqa: BLE001
        return "unavailable"


def sourceHashes() -> Dict[str, str]:
    out = {}
    for f in sorted(_HERE.glob("*.py")):
        out[f"scripts/phase2/{f.name}"] = hashlib.sha256(f.read_bytes()).hexdigest()
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--n-samples", type=int, default=DEFAULT_N)
    ap.add_argument("--seed", type=int, default=DEFAULT_SEED)
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--model-forms", default="boron,free",
                    help="comma separated subset of {boron,free}")
    ap.add_argument("--protocols", default="nitrogen,boron",
                    help="comma separated subset of {nitrogen,boron}")
    ap.add_argument("--epsilons", default=",".join(f"{e:g}" for e in EPSILON_LADDER))
    ap.add_argument("--attractor-t-end", type=float, default=protocols.BORON_T_END)
    ap.add_argument("--label", default="matched",
                    help="free-text label recorded in the manifest")
    args = ap.parse_args(argv)

    forms = [s.strip() for s in args.model_forms.split(",") if s.strip()]
    protos = [s.strip() for s in args.protocols.split(",") if s.strip()]
    epsilons = [float(s) for s in args.epsilons.split(",") if s.strip()]
    for f in forms:
        if f not in ("boron", "free"):
            ap.error(f"unknown model form {f!r}")
    for p in protos:
        if p not in protocols.PROTOCOLS:
            ap.error(f"unknown protocol {p!r}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    matrix = sampleMatrix(args.n_samples, args.seed)

    manifest = {
        "benchmark": "phase2A_matched_equivalence",
        "label": args.label,
        "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "host": platform.node(),
        "repo_head": gitHead(_HERE.parents[1]),
        "design": describeDesign(args.n_samples, args.seed),
        "factorial": {"model_forms": forms, "protocols": protos, "epsilons": epsilons},
        "attractor_t_end": args.attractor_t_end,
        "threshold_H": THRESHOLD_H,
        "environment": {
            "python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__, "executable": sys.executable,
            "cpu_count": os.cpu_count(), "workers": args.workers,
            "thread_env": {k: os.environ.get(k) for k in
                           ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                            "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS")},
        },
        "source_sha256": sourceHashes(),
        "cells": [],
        "note": ("labels are operational. 'no_bounded_attractor_operational' means "
                 "no attractor was found within this finite protocol; it is not a "
                 "mathematical nonexistence result. no percentage in this output "
                 "may be interpreted until every matched cell has completed."),
    }
    (args.outdir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    for form in forms:
        for eps in epsilons:
            for proto in protos:
                summary = runCell(matrix, form, eps, proto, args.outdir,
                                  args.workers, args.attractor_t_end)
                manifest["cells"].append(summary)
                print(f"[{time.strftime('%H:%M:%S')}] {summary['cell']}: "
                      f"{summary['labels']}  ({summary['wall_seconds']:.1f}s)",
                      flush=True)
                # rewrite after every cell so a killed job still has a usable record
                manifest["finished_utc"] = time.strftime(
                    "%Y-%m-%dT%H:%M:%SZ", time.gmtime())
                (args.outdir / "manifest.json").write_text(
                    json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    manifest["complete"] = True
    (args.outdir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"done: {len(manifest['cells'])} cells -> {args.outdir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
