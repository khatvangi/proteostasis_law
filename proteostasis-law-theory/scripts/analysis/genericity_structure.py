"""(G4) is not an independent condition. It follows from (G1).

Task B3. The manuscript through v5 says "(G4) is not automatic ... no structural
feature of the model forces it", and reports Table 1's two near-failures of (G1)
and (G4) as a coincidence about two networks. It is a theorem.

THE ARGUMENT. Because `j` enters `du/dt` additively and nothing else,
`F_j = (1, 0)` exactly, so (G4) reduces to `w_1 != 0` for the left null vector
`w`. Mass balance gives `du/dt = j - R - G`, so

    J = [ -R_u - G_u   -R_a - G_a ]
        [   G_u          G_a      ]

and at a fold the gradients are parallel, `grad R = lambda grad G`. Then
`w = (1, 1+lambda)` satisfies `w J = -(grad R - lambda grad G) = 0`. Hence
`w_1 = 1 != 0` whenever `w` exists at all, and normalised,

    |w_1| / ||w||  =  ( 1 + (1+lambda)^2 )^{-1/2}.                        (*)

Transversality is therefore automatic given (G1), and its margin is SMALL
exactly when `|lambda| -> infinity`, which is exactly when `grad G -> 0` at fixed
`grad R`. The co-occurrence of the (G1) and (G4) near-failures is forced.

WHAT THIS SCRIPT CHECKS. (*) is a genuine check rather than a restatement,
because the two sides are computed by independent routes: the left-hand side
from a singular value decomposition of the ANALYTIC Jacobian, the right-hand side
from central-difference gradients of `R` and `G`. Nothing connects them except
the theorem.

Task B7 rides along, because it is the same state loop: Table 1 reports
`|tr J|`, which hides the sign, and the sign is the classification. At a
two-dimensional fold the eigenvalues are `0` and `tr J`, so `tr J < 0` is a
stable node colliding with a saddle -- the fold terminates a stable branch --
while `tr J > 0` is an unstable node colliding with a saddle, and the fold is
not a stability boundary at all.
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

OUT_JSON = REPO_ROOT / "data" / "computed" / "genericity_structure.json"
OUT_TSV = REPO_ROOT / "data" / "computed" / "genericity_structure.tsv"


def structureAt(u: float, a: float, p: M.Params) -> dict | None:
    """the two sides of (*), plus the signed trace, at one solved fold state."""
    try:
        J = M.jacobian(u, a, p)
    except (M.ModelError, np.linalg.LinAlgError):
        return None
    Gu, Ga = FT._centralGradient(FT.aggregateG, u, a, p)
    Ru, Ra = FT._centralGradient(FT.removalR, u, a, p)
    n_G2 = Gu * Gu + Ga * Ga
    if n_G2 <= 0.0:
        return None

    # route 1: eigenvector solve on the analytic Jacobian
    _, sv, vt = np.linalg.svd(J.T)
    w = vt[-1, :]
    lhs = float(abs(w[0]) / max(float(np.linalg.norm(w)), 1e-300))

    # route 2: the gradient ratio, from central differences
    lam = float((Ru * Gu + Ra * Ga) / n_G2)
    rhs = float(1.0 / np.sqrt(1.0 + (1.0 + lam) ** 2))

    tr = float(J[0, 0] + J[1, 1])
    return {
        "u": u, "a": a,
        "lambda": lam,
        "transversality_svd": lhs,
        "transversality_closed_form": rhs,
        "abs_diff": abs(lhs - rhs),
        "rel_diff": abs(lhs - rhs) / max(rhs, 1e-300),
        "grad_G": float(np.hypot(Gu, Ga)),
        "grad_R": float(np.hypot(Ru, Ra)),
        "tr_J": tr,
        "det_J": float(np.linalg.det(J)),
        "sv_ratio": float(sv[-1] / max(sv[0], 1e-300)),
    }


def _worker(item):
    import genericity as GEN
    i, r = item
    try:
        p = FT.paramsFromSampleRow(r)
        u0, a0 = FT.foldStateFromSampleRow(r)
    except (M.ModelError, ValueError, KeyError):
        return None
    s = GEN.polishFold(p, u0, a0) or FT.foldSolve(p)
    if s is None:
        return None
    out = structureAt(float(s[1]), float(s[2]), p)
    if out is not None:
        out["name"] = f"draw{i}"
    return out


def _loadGridWorker(item):
    name, p, u, a = item
    out = structureAt(u, a, p)
    if out is not None:
        out["name"] = name
    return out


def _pool(items, fn, workers=None):
    if workers is None:
        workers = max(1, min(56, (os.cpu_count() or 2) - 8))
    with mp.get_context("fork").Pool(processes=workers) as pool:
        return pool.map(fn, items, chunksize=8)


def run() -> dict:
    import genericity as GEN
    run_dir = FT.phase1RunDir()
    c = pd.read_csv(run_dir / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    c = c[pd.to_numeric(c["fold_burden"], errors="coerce").notna()]

    pops = {
        "kinetic_box": pd.DataFrame(
            [r for r in _pool([(i, r) for i, r in c.iterrows()], _worker)
             if r is not None]),
        "load_grid": pd.DataFrame(
            [r for r in _pool(GEN.loadGridStates(run_dir), _loadGridWorker)
             if r is not None]),
    }

    out = {}
    for pop, D in pops.items():
        D["population"] = pop
        pos = D["tr_J"] > 0.0
        out[pop] = {
            "n": int(len(D)),
            # B3: the closed form against the eigenvector solve
            "closed_form_abs_diff_median": float(D["abs_diff"].median()),
            "closed_form_abs_diff_p99": float(D["abs_diff"].quantile(0.99)),
            "closed_form_abs_diff_max": float(D["abs_diff"].max()),
            "closed_form_rel_diff_median": float(D["rel_diff"].median()),
            "closed_form_rel_diff_p99": float(D["rel_diff"].quantile(0.99)),
            "closed_form_rel_diff_max": float(D["rel_diff"].max()),
            "n_rel_diff_over_1e6": int((D["rel_diff"] > 1e-6).sum()),
            "transversality_min": float(D["transversality_svd"].min()),
            "lambda_absmax": float(D["lambda"].abs().max()),
            "grad_G_min": float(D["grad_G"].min()),
            # B7: the sign of the trace, which |tr J| hides
            "tr_J_min": float(D["tr_J"].min()),
            "tr_J_max": float(D["tr_J"].max()),
            "tr_J_median": float(D["tr_J"].median()),
            "n_tr_J_positive": int(pos.sum()),
            "frac_tr_J_positive": float(pos.mean()),
            "abs_tr_J_min": float(D["tr_J"].abs().min()),
        }
    pd.concat(pops.values()).to_csv(OUT_TSV, sep="\t", index=False)
    return out


def main() -> int:
    o = run()
    OUT_JSON.write_text(json.dumps(o, indent=2) + "\n")
    print(json.dumps(o, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
