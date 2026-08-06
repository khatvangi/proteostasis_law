"""parameter sampling, factorial grids and deterministic parallel evaluation.

determinism rules used throughout:
  * every sample set is generated from an explicit integer seed;
  * `parallelMap` preserves input order (`Pool.map`, not `imap_unordered`), so
    the output table does not depend on worker scheduling;
  * worker functions must be module-level and take a single picklable argument.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable, Dict, Iterable, List, Sequence, Tuple

import multiprocessing as mp
import numpy as np
from scipy.stats import qmc

from .model import Params


@dataclass(frozen=True)
class ParamRange:
    """sampling range for one parameter.

    `scale` is 'log' for rate constants and affinities (they span orders of
    magnitude and a linear prior would concentrate all mass at the top) and
    'linear' for bounded shape parameters such as the nucleation order.
    """
    name: str
    lo: float
    hi: float
    scale: str = "log"

    def draw(self, unit: np.ndarray) -> np.ndarray:
        if self.scale == "log":
            return 10.0 ** (np.log10(self.lo) + unit * (np.log10(self.hi) - np.log10(self.lo)))
        if self.scale == "linear":
            return self.lo + unit * (self.hi - self.lo)
        raise ValueError(f"unknown scale '{self.scale}'")


#: default global ranges for experiment C. deliberately wide: the point is to
#: find where the qualitative predictions FAIL, not to flatter them.
GLOBAL_RANGES: Tuple[ParamRange, ...] = (
    ParamRange("nu", 1e-2, 1e1, "log"),
    ParamRange("chi", 0.05, 0.95, "linear"),       # rescue allocation, mapped to c_tot/d_tot
    ParamRange("rho_U", 1e-1, 1e1, "log"),
    ParamRange("rho_A", 1e-2, 1e1, "log"),
    ParamRange("alpha_n", 1e-3, 1e1, "log"),
    ParamRange("alpha_g", 1e-2, 1e2, "log"),
    ParamRange("alpha_d", 1e-2, 1e1, "log"),
    ParamRange("m", 1.5, 4.0, "linear"),
    ParamRange("kappa_ref", 1e-2, 1e1, "log"),
    ParamRange("kappa_u", 1e-2, 1e1, "log"),
    ParamRange("kappa_a", 1e-2, 1e1, "log"),
    ParamRange("kappa_dis", 1e-2, 1e1, "log"),
    ParamRange("kappa_cu", 1e-2, 1e1, "log"),
    ParamRange("kappa_ca", 1e-2, 1e1, "log"),
    ParamRange("kappa_du", 1e-2, 1e1, "log"),
    ParamRange("kappa_da", 1e-2, 1e1, "log"),
)


def rangesFromConfig(spec: Sequence[Dict[str, Any]] | None) -> Tuple[ParamRange, ...]:
    if not spec:
        return GLOBAL_RANGES
    return tuple(ParamRange(**s) for s in spec)


def latinHypercube(ranges: Sequence[ParamRange], n: int, seed: int) -> List[Dict[str, float]]:
    """scrambled latin hypercube sample; reproducible for a given (n, seed)."""
    sampler = qmc.LatinHypercube(d=len(ranges), seed=seed)
    unit = sampler.random(n)
    return [
        {r.name: float(r.draw(np.array([unit[i, k]]))[0]) for k, r in enumerate(ranges)}
        for i in range(n)
    ]


def paramsFromSample(sample: Dict[str, float], base: Params,
                     rescue_total: float = 1.0) -> Params:
    """build a `Params` from an LHS draw.

    `chi` (rescue allocation) is expanded into the conserved pools c_tot and
    d_tot at fixed total rescue, so allocation and total capacity stay
    independently controllable.
    """
    kw = {k: v for k, v in sample.items() if k in Params.__dataclass_fields__}
    p = base.with_(**kw)
    if "chi" in sample:
        chi = float(sample["chi"])
        p = p.with_(c_tot=chi * rescue_total, d_tot=(1.0 - chi) * rescue_total)
    return p.validate()


def factorialGrid(axes: Dict[str, Sequence[float]]) -> List[Dict[str, float]]:
    """full factorial over named axes, in a deterministic (sorted-key) order."""
    names = sorted(axes)
    mesh = np.meshgrid(*[np.asarray(axes[n], dtype=float) for n in names], indexing="ij")
    flat = [m.ravel() for m in mesh]
    return [{n: float(flat[k][i]) for k, n in enumerate(names)}
            for i in range(flat[0].size)]


def parallelMap(fn: Callable[[Any], Any], items: Sequence[Any], n_workers: int,
                chunksize: int | None = None, progress_every: int = 0) -> List[Any]:
    """order-preserving parallel map with a serial fallback for n_workers <= 1."""
    items = list(items)
    if n_workers <= 1 or len(items) <= 1:
        out = []
        for i, it in enumerate(items):
            out.append(fn(it))
            if progress_every and (i + 1) % progress_every == 0:
                print(f"  [{i+1}/{len(items)}]", flush=True)
        return out
    if chunksize is None:
        chunksize = max(1, len(items) // (n_workers * 8))
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=n_workers) as pool:
        return pool.map(fn, items, chunksize=chunksize)


def resolveWorkers(requested: int | None) -> int:
    """respect the SLURM allocation rather than the physical core count."""
    import os
    if requested and requested > 0:
        cap = requested
    else:
        cap = int(os.environ.get("SLURM_CPUS_PER_TASK", 0)) or (os.cpu_count() or 1)
    return max(1, min(cap, int(os.environ.get("SLURM_CPUS_PER_TASK", 0)) or cap))
