"""how far the coded rate laws sit from a mechanistically exact closure.

Task A1 asks whether the rate laws in `scripts/proteostasis/model.py` follow
from an explicit reaction network or are phenomenological closures. The binding
layer is exact (see theory/RATE_LAWS.md, Part 1). The CATALYTIC layer is not,
and this script measures by how much, at the states the paper actually reports:
the 2767 re-solved fold states of the kinetic box.

Three approximations are separable and each has an observable:

1. NON-COMPETITION. `v_ref` and `v_dis` both draw on the same free chaperone
   `cf`, and the code gives each its own independent Michaelis factor. A single
   machine processing both substrates would share one denominator, so the two
   catalytic occupancies would sum to at most 1. The coded form allows
   `s_ref + s_dis` up to 2. The excess over 1 is the size of the error.
   Same question for the protease pool with `s_u + s_a`.

2. NEGLECTED SUBSTRATE. If the productive complexes are counted inside `cf`
   and `df` (which is what gives the Michaelis factor its saturation), then the
   substrate they hold ought to appear in the totals `u` and `a`. It does not.
   The neglected fraction is the classical sQSSA error, here
   `(cf*s_ref + df*s_u)/u` for the soluble pool.

3. FREE ENZYME VERSUS BOUND COMPLEX. The code writes the catalytic flux against
   free enzyme; the standard scheme writes it against the bound complex. That
   comparison is Task A3 and is not measured here.

Nothing here is a pass/fail gate. It supplies the numbers that decide whether
A1 lands as a derivation or as a declared closure.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import fold_theorem as FT  # noqa: E402

OUT = REPO_ROOT / "data" / "computed" / "closure_audit.json"


def closureDiagnostics(u: float, a: float, p: M.Params) -> dict:
    """the three approximation measures at one state."""
    uf, af, cf, df = M.solveFreePools(u, a, p)

    s_ref = uf / (p.kappa_ref + uf)
    s_dis = af / (p.kappa_dis + af)
    s_u = uf / (p.kappa_u + uf)
    s_a = af / (p.kappa_a + af)

    # 1. non-competition: a shared machine would cap each sum at 1
    chap_sum = s_ref + s_dis
    prot_sum = s_u + s_a

    # 2. neglected substrate held in productive complexes, as a fraction of the
    #    total pool. the catalytic constants are absorbed into the rate
    #    coefficients, so the complex concentrations are the flux prefactors.
    held_u = cf * s_ref + df * s_u
    held_a = cf * s_dis + df * s_a

    return {
        "s_ref": s_ref, "s_dis": s_dis, "s_u": s_u, "s_a": s_a,
        "chap_occupancy_sum": chap_sum,
        "prot_occupancy_sum": prot_sum,
        "neglected_u_frac": held_u / max(u, 1e-300),
        "neglected_a_frac": held_a / max(a, 1e-300),
        "u": u, "a": a,
    }


def _worker(item):
    """one network: re-solve the fold, then measure the closure there."""
    import genericity as GEN
    _, r = item
    try:
        p = FT.paramsFromSampleRow(r)
        u0, a0 = FT.foldStateFromSampleRow(r)
    except (M.ModelError, ValueError, KeyError):
        return None
    s = GEN.polishFold(p, u0, a0) or FT.foldSolve(p)
    if s is None:
        return None
    _, u, a = s
    try:
        return closureDiagnostics(u, a, p)
    except (M.ModelError, ValueError, KeyError):
        return None


def audit(workers: int | None = None) -> dict:
    run = FT.phase1RunDir()
    c = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    c = c[pd.to_numeric(c["fold_burden"], errors="coerce").notna()]
    n_admit = len(c)

    if workers is None:
        workers = max(1, min(56, (os.cpu_count() or 2) - 8))
    items = [(i, r) for i, r in c.iterrows()]
    with mp.get_context("fork").Pool(processes=workers) as pool:
        rows = pool.map(_worker, items, chunksize=8)

    D = pd.DataFrame([r for r in rows if r is not None])

    def q(col):
        v = D[col].to_numpy(float)
        v = v[np.isfinite(v)]
        return {"median": float(np.median(v)), "p99": float(np.quantile(v, 0.99)),
                "max": float(v.max()), "n": int(v.size)}

    out = {
        "n_admit_fold": int(n_admit),
        "n_solved": int(len(D)),
        "chap_occupancy_sum": q("chap_occupancy_sum"),
        "prot_occupancy_sum": q("prot_occupancy_sum"),
        "neglected_u_frac": q("neglected_u_frac"),
        "neglected_a_frac": q("neglected_a_frac"),
        # the counts that decide whether the approximation is benign
        "n_chap_sum_over_1": int((D["chap_occupancy_sum"] > 1.0).sum()),
        "n_prot_sum_over_1": int((D["prot_occupancy_sum"] > 1.0).sum()),
        "n_neglected_u_over_10pct": int((D["neglected_u_frac"] > 0.10).sum()),
        "n_neglected_a_over_10pct": int((D["neglected_a_frac"] > 0.10).sum()),
    }
    return out


def main() -> int:
    o = audit()
    OUT.write_text(json.dumps(o, indent=2) + "\n")
    print(json.dumps(o, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
