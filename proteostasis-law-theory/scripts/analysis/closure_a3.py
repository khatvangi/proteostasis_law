"""Task A3: are the quantitative results a property of the model or of the closure?

`theory/RATE_LAWS.md` establishes that the catalytic layer of the model is a
stipulated closure rather than a derivation. The structural results are
unaffected -- Theorem 1 needs only H1-H3, which hold for any closure -- but
`phi`, the utilisation fractions, the ceiling factor and the Hopf incidence are
all computed FROM the rate laws.

This script recomputes them under the mechanistically standard alternative in
which catalytic flux is proportional to the bound complex:

    michaelis (published):  v_ref = c_f * u_f/(K_ref + u_f)
    complex   (this test):  v_ref = c_f * u_f/K_CU  =  [CU]

The second removes both stipulations at once: there is no second saturation, and
cycles sharing a pool compete for it automatically, because the complexes they
form are counted in the same conservation law. The removal CEILING is unchanged,
`c_tot + (rho_U + rho_A) d_tot`, since every complex concentration is bounded by
its own pool total -- so `phi = j_crit/ceiling` means the same thing in both and
the comparison is well posed.

Same kinetic box, same draws, same seeding. The two pipelines differ in exactly
one field of `Params`.

A NOTE ON WHAT IS COMPARED. `s_ref`, `s_u`, `s_a` are Michaelis factors and have
no counterpart under the complex closure. What both closures possess is the
UTILISATION of nominal capacity, `v/V_max`, which is what the phrase "running at
x% of V_max" should mean and is the quantity `phi` actually aggregates. Under
the michaelis closure the two differ by the free-pool fraction:
`util_ref = (c_f/c_tot) * s_ref <= s_ref`. Both are reported.

Stage 2, the Hopf incidence, re-traces the low-burden branch under each closure
and is expensive; it runs only with --hopf.
"""

from __future__ import annotations

import argparse
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

OUT = REPO_ROOT / "data" / "computed" / "closure_a3.json"
CLOSURES = ("michaelis", "complex")


def utilisation(u: float, a: float, p: M.Params) -> dict:
    """where the fold sits, in quantities both closures possess."""
    f = M.fluxes(u, a, p)
    uf, af, cf, df = f["uf"], f["af"], f["cf"], f["df"]
    ceiling = M.removalCeiling(p)
    R = f["refold"] + f["degrade_u"] + f["degrade_a"]
    out = {
        "phi": R / ceiling,
        "util_ref": f["refold"] / p.c_tot,
        "util_u": f["degrade_u"] / (p.rho_U * p.d_tot),
        "util_a": f["degrade_a"] / (p.rho_A * p.d_tot),
        "cf_frac": cf / p.c_tot,
        "df_frac": df / p.d_tot,
        "u": u, "a": a,
    }
    if p.closure == "michaelis":
        out["s_ref"] = uf / (p.kappa_ref + uf)
        out["s_u"] = uf / (p.kappa_u + uf)
        out["s_a"] = af / (p.kappa_a + af)
    return out


def _foldWorker(item):
    """one draw, both closures, seeded identically from the recorded state."""
    import genericity as GEN
    i, r = item
    try:
        base = FT.paramsFromSampleRow(r)
        u0, a0 = FT.foldStateFromSampleRow(r)
    except (M.ModelError, ValueError, KeyError):
        return None
    row = {"name": f"draw{i}"}
    for clos in CLOSURES:
        p = base.with_(closure=clos).validate()
        try:
            s = GEN.polishFold(p, u0, a0) or FT.foldSolve(p)
        except (M.ModelError, ValueError, OverflowError, np.linalg.LinAlgError):
            s = None
        if s is None:
            row[clos] = None
            continue
        _, u, a = s
        try:
            row[clos] = utilisation(u, a, p)
        except (M.ModelError, ValueError, KeyError):
            row[clos] = None
    return row


def _pool(items, fn, workers=None):
    if workers is None:
        workers = max(1, min(56, (os.cpu_count() or 2) - 8))
    with mp.get_context("fork").Pool(processes=workers) as pool:
        return pool.map(fn, items, chunksize=8)


def _stats(D: pd.DataFrame, cols) -> dict:
    o = {}
    for c in cols:
        if c not in D.columns:
            continue
        v = D[c].to_numpy(float)
        v = v[np.isfinite(v)]
        if not v.size:
            continue
        o[c] = {"median": float(np.median(v)),
                "p5": float(np.quantile(v, 0.05)),
                "p95": float(np.quantile(v, 0.95)),
                "n": int(v.size)}
    return o


def foldComparison() -> dict:
    run = FT.phase1RunDir()
    c = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    c = c[pd.to_numeric(c["fold_burden"], errors="coerce").notna()]
    n_admit = len(c)

    rows = [r for r in _pool([(i, r) for i, r in c.iterrows()], _foldWorker)
            if r is not None]

    out = {"n_draws_admitting_fold_under_michaelis": int(n_admit),
           "n_rows": len(rows)}
    frames = {}
    for clos in CLOSURES:
        D = pd.DataFrame([r[clos] for r in rows if r.get(clos) is not None])
        frames[clos] = D
        cols = ["phi", "util_ref", "util_u", "util_a", "cf_frac", "df_frac",
                "s_ref", "s_u", "s_a"]
        out[clos] = {"n_solved": int(len(D)), **_stats(D, cols)}
        out[clos]["inverse_phi_at_median"] = float(1.0 / np.median(D["phi"]))

    # the paired comparison: only networks solved under BOTH closures
    both = [r for r in rows
            if r.get("michaelis") is not None and r.get("complex") is not None]
    if both:
        pm = np.array([r["michaelis"]["phi"] for r in both], float)
        pc = np.array([r["complex"]["phi"] for r in both], float)
        ratio = pc / pm
        out["paired"] = {
            "n": len(both),
            "phi_median_michaelis": float(np.median(pm)),
            "phi_median_complex": float(np.median(pc)),
            "median_ratio_complex_over_michaelis": float(np.median(ratio)),
            "p5_ratio": float(np.quantile(ratio, 0.05)),
            "p95_ratio": float(np.quantile(ratio, 0.95)),
            "n_ratio_outside_half_to_two": int(
                ((ratio < 0.5) | (ratio > 2.0)).sum()),
        }
        # the decision the work order names: does the headline move by >2x?
        m = out["paired"]["phi_median_michaelis"]
        k = out["paired"]["phi_median_complex"]
        out["headline_moves_more_than_twofold"] = bool(
            max(m / k, k / m) > 2.0)
    return out


# ---------------------------------------------------------------------------
# stage 2: does the Hopf incidence survive the closure?
# ---------------------------------------------------------------------------


def _hopfWorker(item):
    import hopf_check as HC
    name, p, u, a = item
    try:
        s = HC.branchSummary(p, u, a, n=150)
    except Exception:
        s = None
    if s is None:
        return {"name": name, "traced": False}
    s["name"] = name
    s["traced"] = True
    return s


def _solveWorker(item):
    """re-solve one draw's fold under one closure. top level so it pickles."""
    import genericity as GEN
    i, r, clos = item
    try:
        p = FT.paramsFromSampleRow(r).with_(closure=clos).validate()
        u0, a0 = FT.foldStateFromSampleRow(r)
    except (M.ModelError, ValueError, KeyError):
        return None
    try:
        s = GEN.polishFold(p, u0, a0) or FT.foldSolve(p)
    except (M.ModelError, ValueError, OverflowError, np.linalg.LinAlgError):
        return None
    if s is None:
        return None
    return (f"draw{i}", p, float(s[1]), float(s[2]))


def hopfComparison() -> dict:
    """crossing incidence under each closure, at each closure's own fold states."""
    run = FT.phase1RunDir()
    c = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    c = c[pd.to_numeric(c["fold_burden"], errors="coerce").notna()]

    out = {}
    for clos in CLOSURES:
        states = [s for s in
                  _pool([(i, r, clos) for i, r in c.iterrows()], _solveWorker)
                  if s is not None]

        S = pd.DataFrame(_pool(states, _hopfWorker))
        S.to_csv(REPO_ROOT / "data" / "computed" / f"closure_a3_hopf_{clos}.tsv",
                 sep="\t", index=False)
        T = S[S["traced"] == True]  # noqa: E712
        cross = T[T["tr_max"] >= 0.0]
        clean = cross[cross["fold_is_j_max"] == 1]
        out[clos] = {
            "n_fold_states": len(states),
            "n_traced": int(len(T)),
            "n_trace_failed": int(len(S) - len(T)),
            "n_cross": int(len(cross)),
            "n_cross_clean": int(len(clean)),
            "incidence_pct": float(100.0 * len(cross) / max(len(T), 1)),
            "incidence_clean_pct": float(100.0 * len(clean) / max(len(T), 1)),
            "tr_max_median": float(T["tr_max"].median()) if len(T) else None,
        }
        if len(clean):
            f = clean["j_at_first_cross"] / clean["j_crit"]
            out[clos]["first_cross_frac_median"] = float(f.median())
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--hopf", action="store_true",
                    help="also re-trace the branch under each closure (slow)")
    args = ap.parse_args()

    o = foldComparison()
    if args.hopf:
        o["hopf"] = hopfComparison()
    OUT.write_text(json.dumps(o, indent=2) + "\n")
    print(json.dumps(o, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
