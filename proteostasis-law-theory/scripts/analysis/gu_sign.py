"""B6: where is `G_u > 0` structural, and where is it only a competition?

`theory/LEMMA0_BINDING.md` (vi) proves: if `G_af >= 0` then `G_u > 0`. Three of
the four chain-rule products are nonnegative unconditionally, by the M-matrix
property of the binding Jacobian; the fourth carries the sign of

    G_af = alpha_g uf - alpha_d cf kappa_dis/(kappa_dis+af)^2
                      - rho_A  df kappa_a  /(kappa_a  +af)^2

This measures the two things the manuscript needs and cannot assert:

1. HOW OFTEN the sufficient condition holds, and what `G_u` actually is where it
   fails. The standing claim rests on 40 kinetic draws with a minimum of
   6.67e-12 -- a number small enough to be an `a -> 0` boundary artifact rather
   than a margin, which is why the burden is reported beside it here.

2. WHETHER the `G_af < 0` states are the crossing networks of Section 7. The
   oscillatory region is characterised by fast, sharply saturating aggregate
   clearance against growth-dominated aggregation, and that is a statement about
   which term dominates inside `G_af`. If the two sets coincide, Sections 3.3 and
   7 are one object seen from two sides and the paper should say so once.
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

COMPUTED = REPO_ROOT / "data" / "computed"
OUT = COMPUTED / "gu_sign.json"
OUT_TSV = COMPUTED / "gu_sign.tsv"


def signsAt(u: float, a: float, p: M.Params) -> dict | None:
    """the four chain-rule products, G_af, and G_u, at one state."""
    try:
        uf, af, cf, df = M.solveFreePools(u, a, p)
    except (M.ModelError, ValueError):
        return None

    G_uf = p.alpha_n * p.m * (uf ** (p.m - 1.0)) + p.alpha_g * af
    # the three terms of G_af, kept apart. Section 7 describes the oscillatory
    # region as sharply saturating clearance against growth-dominated
    # aggregation, which is a claim about WHICH of these dominates -- and the
    # ensemble ratios alone do not settle it, since rho_A is 7.5x HIGHER there
    # while kappa_a is 11.7x lower and the two enter the clearance term as a
    # product. Measured, not inferred.
    t_growth = p.alpha_g * uf
    t_disagg = -p.alpha_d * cf * p.kappa_dis / (p.kappa_dis + af) ** 2
    t_clear = -p.rho_A * df * p.kappa_a / (p.kappa_a + af) ** 2
    G_af = t_growth + t_disagg + t_clear
    G_cf = -p.alpha_d * af / (p.kappa_dis + af)
    G_df = -p.rho_A * af / (p.kappa_a + af)

    # the closure sensitivities, by central difference on the solved pools
    h = 1e-6 * max(u, 1e-12)
    try:
        ufp, afp, cfp, dfp = M.solveFreePools(u + h, a, p)
        ufm, afm, cfm, dfm = M.solveFreePools(max(u - h, 1e-300), a, p)
    except (M.ModelError, ValueError):
        return None
    d = 2.0 * h
    duf, daf = (ufp - ufm) / d, (afp - afm) / d
    dcf, ddf = (cfp - cfm) / d, (dfp - dfm) / d

    try:
        Gu, _ = FT._centralGradient(FT.aggregateG, u, a, p)
    except (M.ModelError, OverflowError):
        return None

    return {
        "u": u, "a": a, "uf": uf, "af": af,
        "G_af": float(G_af), "G_u": float(Gu),
        "t_growth": float(t_growth), "t_disagg": float(t_disagg),
        "t_clear": float(t_clear),
        "p_uf": float(G_uf * duf), "p_af": float(G_af * daf),
        "p_cf": float(G_cf * dcf), "p_df": float(G_df * ddf),
        "duf_du": float(duf), "daf_du": float(daf),
        "dcf_du": float(dcf), "ddf_du": float(ddf),
        "chain_sum": float(G_uf * duf + G_af * daf + G_cf * dcf + G_df * ddf),
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
    out = signsAt(float(s[1]), float(s[2]), p)
    if out is not None:
        out["name"] = f"draw{i}"
    return out


def _gridWorker(item):
    name, p, u, a = item
    out = signsAt(u, a, p)
    if out is not None:
        out["name"] = name
    return out


def run(workers: int | None = None) -> dict:
    import genericity as GEN
    run_dir = FT.phase1RunDir()
    c = pd.read_csv(run_dir / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    c = c[pd.to_numeric(c["fold_burden"], errors="coerce").notna()]

    if workers is None:
        workers = max(1, min(56, (os.cpu_count() or 2) - 8))
    with mp.get_context("fork").Pool(processes=workers) as pool:
        kin = [r for r in pool.map(_worker, [(i, r) for i, r in c.iterrows()],
                                   chunksize=8) if r is not None]
        grid = [r for r in pool.map(_gridWorker, GEN.loadGridStates(run_dir),
                                    chunksize=4) if r is not None]

    K = pd.DataFrame(kin).assign(population="kinetic_box")
    L = pd.DataFrame(grid).assign(population="load_grid")
    D = pd.concat([K, L], ignore_index=True)
    D.to_csv(OUT_TSV, sep="\t", index=False)

    out = {}
    for lab, S in (("kinetic_box", K), ("load_grid", L)):
        neg = S["G_af"] < 0.0
        out[lab] = {
            "n": int(len(S)),
            "n_G_af_negative": int(neg.sum()),
            "frac_G_af_negative": float(neg.mean()),
            "G_u_min": float(S["G_u"].min()),
            "n_G_u_nonpositive": int((S["G_u"] <= 0.0).sum()),
            # the honest form: the minimum WHERE THE PROOF DOES NOT REACH
            "G_u_min_where_G_af_negative": float(S[neg]["G_u"].min())
            if neg.any() else None,
            "G_u_min_where_G_af_nonneg": float(S[~neg]["G_u"].min())
            if (~neg).any() else None,
            # is the smallest G_u an a -> 0 artifact rather than a margin?
            "a_at_G_u_min": float(S.loc[S["G_u"].idxmin(), "a"]),
            "chain_vs_central_max_rel": float(
                ((S["chain_sum"] - S["G_u"]).abs()
                 / S["G_u"].abs().clip(lower=1e-300)).max()),
        }

    # the overlap: G_af sign against Section 7's groups
    try:
        Z = pd.read_csv(COMPUTED / "hopf_zero_counts.tsv", sep="\t")
        grp = dict(zip(Z["name"], Z["group"]))
        K2 = K.assign(group=K["name"].map(grp).fillna("non_crossing"))
        tab = (K2.assign(G_af_neg=K2["G_af"] < 0.0)
                 .groupby(["group", "G_af_neg"]).size().unstack(fill_value=0))
        out["overlap_with_section7"] = {
            str(g): {str(k): int(v) for k, v in row.items()}
            for g, row in tab.iterrows()}
        for g in ("terminal", "window", "non_crossing"):
            sub = K2[K2["group"] == g]
            if len(sub):
                out.setdefault("frac_G_af_negative_by_group", {})[g] = \
                    float((sub["G_af"] < 0.0).mean())
                # which term carries the sign: the two negative ones are
                # compared against the one positive one on its own scale
                s = sub["t_growth"].abs().clip(lower=1e-300)
                out.setdefault("G_af_terms_by_group", {})[g] = {
                    "n": int(len(sub)),
                    "disagg_over_growth_median": float(
                        (sub["t_disagg"].abs() / s).median()),
                    "clear_over_growth_median": float(
                        (sub["t_clear"].abs() / s).median()),
                    "n_clear_is_larger_negative": int(
                        (sub["t_clear"].abs() > sub["t_disagg"].abs()).sum()),
                }
    except FileNotFoundError:
        out["overlap_with_section7"] = "hopf_zero_counts.tsv absent"
    return out


def main() -> int:
    o = run()
    OUT.write_text(json.dumps(o, indent=2) + "\n")
    print(json.dumps(o, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
