#!/usr/bin/env python3
"""
does the framework's distinguishing prediction hold in its own model?

the paper predicts that perturbations at different burden stages interact
SUPRAADDITIVELY, because the burden terms are causally coupled: independent terms
would combine additively, coupled ones need not. that prediction is untested
experimentally. it can, however, be tested *in the model* -- and if the model
turned out additive, the framework's central claim would be wrong on its own
terms, which is worth knowing before anyone designs an experiment.

design
------
a 2x2 factorial on the two-pool ODE:

    E  raise the per-codon error rate      f_codon *= error_factor      (B_error up)
    C  knock down chaperone capacity       C_tot_uM /= capacity_factor  (C_buffer down)

the readout is the viability margin: how far the operating point sits below the
collapse threshold, in log10 units,

    margin = log10( min( P_dagger / P*,  A_max / A* ) )

taking whichever pool binds. NOTE that both the operating point and the threshold
move under a capacity knockdown -- lowering C_tot reduces v_fold, which lowers
the J-curve and so lowers P_dagger as well. that shared dependence is precisely
the coupling the prediction is about.

damage is margin lost, D = margin_baseline - margin_perturbed, and the
interaction is

    I = D_both - (D_E + D_C)

I > 0 is supraadditive. we also record the sharpest possible form: cases where
each single perturbation leaves the system viable but the combination has no
stable state at all.

three sweeps
------------
1. interaction vs starting margin -- swept over baseline f_codon, to test whether
   the interaction switches on only as the margin closes (prediction 2)
2. an effect-size grid over (error_factor, capacity_factor) at the observed rate,
   so an experiment can be sized
3. the collapse frontier: which single-viable combinations are jointly lethal

parameters are the upstream literature-anchored baseline; nothing is fitted.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
COMP = ROOT / "data" / "computed"
sys.path.insert(0, str(Path(__file__).resolve().parent / "vendor"))

import two_pool_ode as m  # noqa: E402  (vendored, unmodified)

OBSERVED_F = 1e-4          # observed E. coli per-codon mistranslation rate


def J_from_f(f, p):
    """inverse of two_pool_ode.f_codon_from_J."""
    return f * p.N_prot * (1.0 - p.S_avg) * p.p_baseline / p.T_gen_s


def state(f_codon, error_factor=1.0, capacity_factor=1.0):
    """
    operating point, threshold, and viability margin under a perturbation.

    error_factor    > 1 raises the per-codon error rate  (more B_error)
    capacity_factor > 1 knocks down chaperone capacity   (less C_buffer)
    """
    p = m.Params()
    p.C_tot_uM = p.C_tot_uM / capacity_factor

    f = f_codon * error_factor
    J = J_from_f(f, p)

    P_dag, J_crit, mech, _ = m.saddle_node_operational(m.J_curve_two, m.A_qs, p)
    P_star, A_star = m.steady_state(J, p)

    collapsed = (not np.isfinite(P_star)) or (J >= J_crit)
    if collapsed:
        return {"f_codon": f, "error_factor": error_factor,
                "capacity_factor": capacity_factor, "J": J, "J_crit": J_crit,
                "P_star": np.nan, "A_star": np.nan, "P_dagger": P_dag,
                "margin_P": np.nan, "margin_A": np.nan, "margin": np.nan,
                "binding_pool": "none", "mechanism": mech, "collapsed": True}

    margin_P = np.log10(P_dag / P_star)
    margin_A = np.log10(p.A_max / A_star)
    return {"f_codon": f, "error_factor": error_factor,
            "capacity_factor": capacity_factor, "J": J, "J_crit": J_crit,
            "P_star": P_star, "A_star": A_star, "P_dagger": P_dag,
            "margin_P": margin_P, "margin_A": margin_A,
            "margin": min(margin_P, margin_A),
            "binding_pool": "P" if margin_P <= margin_A else "A",
            "mechanism": mech, "collapsed": False}


def factorial(f_base, e, c):
    """the 2x2 and its interaction term."""
    base = state(f_base, 1.0, 1.0)
    only_e = state(f_base, e, 1.0)
    only_c = state(f_base, 1.0, c)
    both = state(f_base, e, c)

    row = {
        "f_base": f_base, "error_factor": e, "capacity_factor": c,
        "margin_baseline": base["margin"],
        "margin_error_only": only_e["margin"],
        "margin_capacity_only": only_c["margin"],
        "margin_both": both["margin"],
        "collapsed_error_only": only_e["collapsed"],
        "collapsed_capacity_only": only_c["collapsed"],
        "collapsed_both": both["collapsed"],
        "binding_pool_baseline": base["binding_pool"],
    }

    # qualitative supraadditivity: both singles survive, the combination does not
    row["qualitative_supraadditive"] = bool(
        both["collapsed"] and not only_e["collapsed"]
        and not only_c["collapsed"] and not base["collapsed"])

    if base["collapsed"] or only_e["collapsed"] or only_c["collapsed"]:
        row.update({"D_error": np.nan, "D_capacity": np.nan, "D_both": np.nan,
                    "interaction": np.nan, "interaction_pct_of_additive": np.nan})
        return row

    D_e = base["margin"] - only_e["margin"]
    D_c = base["margin"] - only_c["margin"]
    if both["collapsed"]:
        # the combination has no stable state; the margin loss is at least the
        # whole baseline margin. record the lower bound, flagged above.
        D_b = base["margin"]
        row.update({"D_error": D_e, "D_capacity": D_c, "D_both": D_b,
                    "interaction": D_b - (D_e + D_c),
                    "interaction_pct_of_additive":
                        100.0 * (D_b - (D_e + D_c)) / (D_e + D_c)
                        if (D_e + D_c) else np.nan})
        return row

    D_b = base["margin"] - both["margin"]
    additive = D_e + D_c
    row.update({"D_error": D_e, "D_capacity": D_c, "D_both": D_b,
                "interaction": D_b - additive,
                "interaction_pct_of_additive":
                    100.0 * (D_b - additive) / additive if additive else np.nan})
    return row


# --------------------------------------------------------------------------
def sweep_margin(e=3.0, c=3.0):
    """sweep 1 -- interaction as a function of how much margin is left."""
    fs = [1e-4, 2e-4, 5e-4, 1e-3, 2e-3, 5e-3]
    rows = [factorial(f, e, c) for f in fs]
    return pd.DataFrame(rows)


def sweep_effect_sizes(f_base=OBSERVED_F):
    """sweep 2 -- effect-size grid at the observed rate."""
    factors = [1.5, 2.0, 3.0, 5.0, 10.0, 20.0]
    rows = [factorial(f_base, e, c) for e in factors for c in factors]
    return pd.DataFrame(rows)


def collapse_frontier(f_base=OBSERVED_F):
    """sweep 3 -- single-viable but jointly lethal combinations."""
    grid = np.geomspace(1.5, 200, 26)
    rows = []
    for e in grid:
        for c in grid:
            r = factorial(f_base, float(e), float(c))
            rows.append({"error_factor": e, "capacity_factor": c,
                         "qualitative_supraadditive": r["qualitative_supraadditive"],
                         "collapsed_both": r["collapsed_both"],
                         "collapsed_error_only": r["collapsed_error_only"],
                         "collapsed_capacity_only": r["collapsed_capacity_only"],
                         "interaction": r["interaction"]})
    return pd.DataFrame(rows)


def main():
    COMP.mkdir(parents=True, exist_ok=True)

    base = state(OBSERVED_F)
    print("=" * 78)
    print("baseline at the observed E. coli rate")
    print("=" * 78)
    print(f"  P* = {base['P_star']:.4e}   P_dagger = {base['P_dagger']:.4e}   "
          f"margin_P = {base['margin_P']:.2f} log10 (x{10**base['margin_P']:.0f})")
    print(f"  A* = {base['A_star']:.4e}   A_max    = 0.25        "
          f"margin_A = {base['margin_A']:.2f} log10 (x{10**base['margin_A']:.0f})")
    print(f"  binding pool: {base['binding_pool']}   "
          f"closure mechanism: {base['mechanism']}")

    s1 = sweep_margin()
    s1.to_csv(COMP / "supraadditivity_margin_sweep.tsv", sep="\t", index=False)
    print("\n" + "=" * 78)
    print("sweep 1 -- interaction vs remaining margin  (error x3, capacity /3)")
    print("=" * 78)
    print(f"  {'f_base':>8} {'margin':>7} {'D_err':>7} {'D_cap':>7} "
          f"{'D_both':>7} {'additive':>9} {'interaction':>12} {'% of add.':>10}")
    for _, r in s1.iterrows():
        add = r.D_error + r.D_capacity
        flag = "  <- both collapse" if r.collapsed_both else ""
        print(f"  {r.f_base:8.0e} {r.margin_baseline:7.2f} {r.D_error:7.3f} "
              f"{r.D_capacity:7.3f} {r.D_both:7.3f} {add:9.3f} "
              f"{r.interaction:+12.3f} {r.interaction_pct_of_additive:+9.1f}%{flag}")

    s2 = sweep_effect_sizes()
    s2.to_csv(COMP / "supraadditivity_effect_grid.tsv", sep="\t", index=False)
    print("\n" + "=" * 78)
    print(f"sweep 2 -- interaction (log10 margin units) at f = {OBSERVED_F:.0e}")
    print("=" * 78)
    piv = s2.pivot_table(index="error_factor", columns="capacity_factor",
                         values="interaction")
    print(piv.round(3).to_string())
    pct = s2.pivot_table(index="error_factor", columns="capacity_factor",
                         values="interaction_pct_of_additive")
    print("\n  as % of the additive expectation:")
    print(pct.round(1).to_string())

    s3 = collapse_frontier()
    s3.to_csv(COMP / "supraadditivity_collapse_frontier.tsv", sep="\t", index=False)
    n_qual = int(s3.qualitative_supraadditive.sum())
    print("\n" + "=" * 78)
    print("sweep 3 -- single-viable but jointly lethal")
    print("=" * 78)
    print(f"  grid points where each single perturbation is survivable but the "
          f"combination collapses: {n_qual} / {len(s3)}")
    if n_qual:
        q = s3[s3.qualitative_supraadditive]
        print(f"  mildest such combination: error x{q.error_factor.min():.1f} "
              f"with capacity /{q[q.error_factor == q.error_factor.min()].capacity_factor.min():.1f}")

    summary = {
        "baseline_margin_log10": base["margin"],
        "baseline_margin_fold": 10 ** base["margin"],
        "binding_pool": base["binding_pool"],
        "perturbations": {"error": "f_codon x factor (raises B_error)",
                          "capacity": "C_tot_uM / factor (lowers C_buffer)"},
        "margin_definition": "log10(min(P_dagger/P*, A_max/A*))",
        "sweep1_error_factor": 3.0,
        "sweep1_capacity_factor": 3.0,
        "interaction_at_observed_rate":
            float(s1.iloc[0].interaction),
        "interaction_pct_at_observed_rate":
            float(s1.iloc[0].interaction_pct_of_additive),
        "interaction_is_supraadditive_at_observed_rate":
            bool(s1.iloc[0].interaction > 0),
        "interaction_grows_as_margin_closes": bool(
            s1.dropna(subset=["interaction"]).interaction.is_monotonic_increasing
            or s1.dropna(subset=["interaction"]).interaction.iloc[-1]
            > s1.dropna(subset=["interaction"]).interaction.iloc[0]),
        "max_interaction_in_effect_grid": float(s2.interaction.max()),
        "n_qualitatively_supraadditive_grid_points": n_qual,
        "n_grid_points": int(len(s3)),
        "caveat": "this tests whether the MODEL predicts supraadditivity, not "
                  "whether cells do. it cannot validate the framework; it checks "
                  "internal coherence and sizes the effect an experiment must "
                  "detect. parameters are the upstream literature-anchored "
                  "baseline and are not fitted.",
    }
    (COMP / "supraadditivity_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nwrote 3 tables + supraadditivity_summary.json to {COMP}")


if __name__ == "__main__":
    main()
