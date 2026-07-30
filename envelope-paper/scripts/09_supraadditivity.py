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

    E  raise the per-codon error rate      f_codon   *= error_factor     (B_error up)
    C  knock down rescue throughput        k_obs_max /= capacity_factor  (C_buffer down)

the readout is the viability margin: how far the operating point sits below the
collapse threshold, in log10 units,

    margin = log10( min( P_dagger / P*,  A_max / A* ) )

taking whichever pool binds. NOTE that both the operating point and the threshold
move under a capacity knockdown -- lowering rescue reduces v_fold, which lowers
the J-curve and so lowers P_dagger as well. that shared dependence is precisely
the coupling the prediction is about.

on the choice of capacity knob: the first version of this analysis used C_tot_uM
and found the capacity arm almost inert. that turned out to be an artifact. at the
observed operating point the folding arm is 97.9% saturated (free chaperone
47.5 uM against 0.052 uM misfolded protein), so a 3-fold C_tot knockdown leaves
v_fold at 96%. k_obs_max is the throughput knob and acts proportionally. both are
reported below, because the discrepancy between them is itself a finding about the
model's parameterization.

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
4. knob comparison: k_obs_max vs C_tot, exposing the saturation artifact

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

# ---------------------------------------------------------------------------
# EVALUATION POINT. this was hardcoded to 1e-4 -- the BOTTOM of the quoted
# observed window -- which is exactly the evaluation point Result 3 rejects as
# inflating the headroom sixfold. Result 5 was therefore anchored at a margin of
# 2.20 log10 while the paper's own internally consistent margin is 1.39, and the
# effect it reported (+0.2% over additive) was ~38x smaller than the effect at
# its own evaluation point (+7.4%). the error ran AGAINST the paper: it argued a
# real result away.
#
# the evaluation point is now read from the burden analysis rather than written
# here, so it cannot drift from Result 1 again. the window bottom is retained as
# an explicitly labelled comparison, not as the anchor.
# ---------------------------------------------------------------------------
F_WINDOW_BOTTOM = 1e-4     # bottom of the quoted observed window; comparison only


def consistent_f():
    """
    the usage-weighted mean per-codon error rate derived in scripts/06.

    read from disk on purpose: hardcoding it is how Result 5 came to be anchored
    at a premise Result 3 spends a page rejecting.
    """
    b = json.loads((COMP / "translation_burden.json").read_text())
    return float(b["usage_weighted_mean_mu_per_codon"])


def J_from_f(f, p):
    """inverse of two_pool_ode.f_codon_from_J."""
    return f * p.N_prot * (1.0 - p.S_avg) * p.p_baseline / p.T_gen_s


def state(f_codon, error_factor=1.0, capacity_factor=1.0,
          capacity_knob="k_obs_max"):
    """
    operating point, threshold, and viability margin under a perturbation.

    error_factor    > 1 raises the per-codon error rate  (more B_error)
    capacity_factor > 1 knocks down chaperone capacity   (less C_buffer)

    WHICH CAPACITY KNOB MATTERS. the rescue term is
    v_fold = k_obs_max · c_free/(c_free + K_d), with
    c_free = C_tot/(1 + M/K_d) and M the misfolded pool. at the observed
    operating point M = 0.052 uM against C_tot = 50 uM and K_d = 1 uM, so
    c_free = 47.5 uM and v_fold sits at 97.9% of k_obs_max -- the folding arm is
    almost fully saturated. lowering C_tot therefore does almost nothing (a
    3-fold knockdown leaves v_fold at 96%), while lowering k_obs_max acts
    proportionally (3-fold -> 33%).

    so `k_obs_max` is the default: it represents a loss of rescue THROUGHPUT,
    which is what a chaperone mutant or an inhibitor does. `C_tot` is retained as
    a comparison because it shows how much of the capacity arm's apparent
    insensitivity is an artifact of the titration term being saturated in this
    parameterization.
    """
    p = m.Params()
    if capacity_knob == "k_obs_max":
        p.k_obs_max = p.k_obs_max / capacity_factor
    elif capacity_knob == "C_tot":
        p.C_tot_uM = p.C_tot_uM / capacity_factor
    else:
        raise ValueError(f"unknown capacity_knob {capacity_knob!r}")

    f = f_codon * error_factor
    J = J_from_f(f, p)

    P_dag, J_crit, mech, _ = m.saddle_node_operational(m.J_curve_two, m.A_qs, p)
    P_star, A_star = m.steady_state(J, p)

    collapsed = (not np.isfinite(P_star)) or (J >= J_crit)
    if collapsed:
        return {"f_codon": f, "error_factor": error_factor,
                "capacity_factor": capacity_factor, "capacity_knob": capacity_knob,
                "J": J, "J_crit": J_crit,
                "P_star": np.nan, "A_star": np.nan, "P_dagger": P_dag,
                "margin_P": np.nan, "margin_A": np.nan, "margin": np.nan,
                "binding_pool": "none", "mechanism": mech, "collapsed": True}

    margin_P = np.log10(P_dag / P_star)
    margin_A = np.log10(p.A_max / A_star)
    return {"f_codon": f, "error_factor": error_factor,
            "capacity_factor": capacity_factor, "capacity_knob": capacity_knob,
            "J": J, "J_crit": J_crit,
            "P_star": P_star, "A_star": A_star, "P_dagger": P_dag,
            "margin_P": margin_P, "margin_A": margin_A,
            "margin": min(margin_P, margin_A),
            "binding_pool": "P" if margin_P <= margin_A else "A",
            "mechanism": mech, "collapsed": False}


def factorial(f_base, e, c, knob="k_obs_max"):
    """the 2x2 and its interaction term."""
    base = state(f_base, 1.0, 1.0, knob)
    only_e = state(f_base, e, 1.0, knob)
    only_c = state(f_base, 1.0, c, knob)
    both = state(f_base, e, c, knob)

    row = {
        "f_base": f_base, "error_factor": e, "capacity_factor": c,
        "capacity_knob": knob,
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
        # the combination has no stable state, so the margin loss is unbounded --
        # the system does not merely lose margin, it loses the branch.
        #
        # DO NOT compute a numeric interaction here. an earlier version of this
        # script substituted D_both = margin_baseline as a "lower bound", which
        # made the interaction come out NEGATIVE whenever D_e + D_c already
        # exceeded the baseline margin -- reporting apparent subadditivity for
        # cases that are in fact the strongest possible supraadditivity. the
        # honest record is the qualitative flag plus the two single-perturbation
        # damages.
        row.update({"D_error": D_e, "D_capacity": D_c, "D_both": np.inf,
                    "interaction": np.nan,
                    "interaction_pct_of_additive": np.nan})
        return row

    D_b = base["margin"] - both["margin"]
    additive = D_e + D_c
    row.update({"D_error": D_e, "D_capacity": D_c, "D_both": D_b,
                "interaction": D_b - additive,
                "interaction_pct_of_additive":
                    100.0 * (D_b - additive) / additive if additive else np.nan})
    return row


# --------------------------------------------------------------------------
def sweep_margin(e=3.0, c=3.0, knob="k_obs_max"):
    """
    sweep 1 -- interaction as a function of how much margin is left.

    the grid spans the observed window and beyond, and explicitly includes both
    the window bottom (for comparison with earlier drafts) and the internally
    consistent evaluation point.
    """
    fc = consistent_f()
    fs = sorted({F_WINDOW_BOTTOM, 2e-4, fc, 1e-3, 2e-3, 5e-3})
    rows = []
    for f in fs:
        r = factorial(f, e, c, knob)
        r["evaluation_point"] = ("window_bottom" if f == F_WINDOW_BOTTOM else
                                 "usage_weighted_mu" if f == fc else "")
        rows.append(r)
    return pd.DataFrame(rows)


def sweep_effect_sizes(f_base=None):
    """sweep 2 -- effect-size grid at the internally consistent evaluation point."""
    f_base = consistent_f() if f_base is None else f_base
    factors = [1.5, 2.0, 3.0, 5.0, 10.0, 20.0]
    rows = [factorial(f_base, e, c) for e in factors for c in factors]
    return pd.DataFrame(rows)


def collapse_frontier(f_base=None):
    """sweep 3 -- single-viable but jointly lethal combinations."""
    f_base = consistent_f() if f_base is None else f_base
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

    F = consistent_f()
    base = state(F)
    wb = state(F_WINDOW_BOTTOM)
    print("=" * 78)
    print(f"baseline at the internally consistent evaluation point "
          f"(f = {F:.3e} /codon)")
    print("=" * 78)
    print(f"  P* = {base['P_star']:.4e}   P_dagger = {base['P_dagger']:.4e}   "
          f"margin_P = {base['margin_P']:.2f} log10 (x{10**base['margin_P']:.0f})")
    print(f"  A* = {base['A_star']:.4e}   A_max    = 0.25        "
          f"margin_A = {base['margin_A']:.2f} log10 (x{10**base['margin_A']:.0f})")
    print(f"  binding pool: {base['binding_pool']}   "
          f"closure mechanism: {base['mechanism']}")
    print(f"\n  for comparison, at the window bottom (f = {F_WINDOW_BOTTOM:.0e}, "
          f"what earlier drafts used):")
    print(f"    margin {wb['margin']:.2f} log10 vs {base['margin']:.2f} here; "
          f"P* {wb['P_star']*300:.3f} uM vs {base['P_star']*300:.3f} uM")

    # saturation diagnostic: how much room each capacity knob actually has
    p0 = m.Params()
    cf = m.c_free(base["P_star"], p0)
    sat = m.v_fold(base["P_star"], p0) / p0.k_obs_max
    print(f"\n  chaperone titration at this operating point:")
    print(f"    misfolded protein M = {base['P_star'] * p0.Prot_tot_uM:.4f} uM  vs  "
          f"free chaperone = {cf:.2f} uM  (C_tot = {p0.C_tot_uM:.0f}, "
          f"K_d = {p0.K_d_uM:.1f})")
    print(f"    v_fold / k_obs_max = {sat:.3f}  -- the folding arm is "
          f"{100*sat:.1f}% saturated, so C_tot is a nearly inert knob")

    s1 = sweep_margin()
    s1.to_csv(COMP / "supraadditivity_margin_sweep.tsv", sep="\t", index=False)
    cons = s1[s1.evaluation_point == "usage_weighted_mu"].iloc[0]
    wbrow = s1[s1.evaluation_point == "window_bottom"].iloc[0]
    print("\n" + "=" * 78)
    print("sweep 1 -- interaction vs remaining margin  (error x3, capacity /3)")
    print("=" * 78)
    print(f"  {'f_base':>8} {'margin':>7} {'D_err':>7} {'D_cap':>7} "
          f"{'D_both':>7} {'additive':>9} {'interaction':>12} {'% of add.':>10}")
    for _, r in s1.iterrows():
        if not np.isfinite(r.D_error):
            print(f"  {r.f_base:8.0e} {r.margin_baseline:7.2f} "
                  f"{'a single perturbation already collapses it':>52}")
            continue
        add = r.D_error + r.D_capacity
        if r.collapsed_both:
            print(f"  {r.f_base:8.0e} {r.margin_baseline:7.2f} {r.D_error:7.3f} "
                  f"{r.D_capacity:7.3f} {'COLLAPSE':>7} {add:9.3f} "
                  f"{'unbounded':>12} {'-- singles both survive':>10}")
            continue
        print(f"  {r.f_base:8.0e} {r.margin_baseline:7.2f} {r.D_error:7.3f} "
              f"{r.D_capacity:7.3f} {r.D_both:7.3f} {add:9.3f} "
              f"{r.interaction:+12.3f} {r.interaction_pct_of_additive:+9.1f}%")

    s2 = sweep_effect_sizes()
    s2.to_csv(COMP / "supraadditivity_effect_grid.tsv", sep="\t", index=False)
    print("\n" + "=" * 78)
    print(f"sweep 2 -- interaction (log10 margin units) at f = {F:.3e}")
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

    # sweep 4 -- knob comparison
    s4 = pd.concat([sweep_margin(knob="k_obs_max").assign(knob="k_obs_max"),
                    sweep_margin(knob="C_tot").assign(knob="C_tot")],
                   ignore_index=True)
    s4.to_csv(COMP / "supraadditivity_knob_comparison.tsv", sep="\t", index=False)
    print("\n" + "=" * 78)
    print("sweep 4 -- capacity knob comparison (error x3, capacity /3)")
    print("=" * 78)
    print(f"  {'knob':<11}{'f_base':>9}{'D_cap':>8}{'D_err':>8}"
          f"{'interaction':>13}{'% of additive':>15}")
    for _, r in s4.iterrows():
        if not np.isfinite(r.interaction):
            if r.collapsed_both and np.isfinite(r.D_error):
                print(f"  {r.knob:<11}{r.f_base:9.0e}{r.D_capacity:8.3f}"
                      f"{r.D_error:8.3f}{'COLLAPSE':>13}{'unbounded':>15}")
            continue
        print(f"  {r.knob:<11}{r.f_base:9.0e}{r.D_capacity:8.3f}{r.D_error:8.3f}"
              f"{r.interaction:+13.3f}{r.interaction_pct_of_additive:+14.1f}%")

    # the experimentally actionable case: modest, individually survivable
    # perturbations that are jointly lethal once the margin is compressed
    actionable = s1[s1.qualitative_supraadditive]
    print(f"\n  at the internally consistent evaluation point: interaction "
          f"{cons.interaction_pct_of_additive:+.2f}% of additive "
          f"(margin {cons.margin_baseline:.2f} log10)")
    print(f"  at the window bottom, the same perturbation pair gives "
          f"{wbrow.interaction_pct_of_additive:+.2f}% "
          f"(margin {wbrow.margin_baseline:.2f}) -- a "
          f"{cons.interaction_pct_of_additive / wbrow.interaction_pct_of_additive:.0f}x "
          f"understatement")
    if len(actionable):
        a = actionable.iloc[0]
        print(f"\n  first margin at which error x{a.error_factor:.0f} and "
              f"capacity /{a.capacity_factor:.0f} are individually survivable "
              f"but jointly lethal:")
        print(f"    f_base = {a.f_base:.0e}, starting margin "
              f"{a.margin_baseline:.2f} log10 (x{10**a.margin_baseline:.0f})")

    summary = {
        "capacity_knob": "k_obs_max (rescue throughput)",
        "qualitative_supraadditivity_first_seen_at_f":
            float(actionable.iloc[0].f_base) if len(actionable) else None,
        "qualitative_supraadditivity_first_seen_at_margin_log10":
            float(actionable.iloc[0].margin_baseline) if len(actionable) else None,
        "folding_arm_saturation_at_observed_rate": float(sat),
        "free_chaperone_uM_at_observed_rate": float(cf),
        "misfolded_protein_uM_at_observed_rate":
            float(base["P_star"] * p0.Prot_tot_uM),
        "baseline_margin_log10": base["margin"],
        "baseline_margin_fold": 10 ** base["margin"],
        "binding_pool": base["binding_pool"],
        "perturbations": {"error": "f_codon x factor (raises B_error)",
                          "capacity": "C_tot_uM / factor (lowers C_buffer)"},
        "margin_definition": "log10(min(P_dagger/P*, A_max/A*))",
        "sweep1_error_factor": 3.0,
        "sweep1_capacity_factor": 3.0,
        "evaluation_point_f_codon": F,
        "evaluation_point_source": "usage-weighted mean per-codon error rate from "
                                   "scripts/06 (translation_burden.json)",
        "baseline_margin_log10_at_window_bottom": float(wb["margin"]),
        "interaction_at_evaluation_point": float(cons.interaction),
        "interaction_pct_at_evaluation_point":
            float(cons.interaction_pct_of_additive),
        "interaction_pct_at_window_bottom": float(wbrow.interaction_pct_of_additive),
        "interaction_is_supraadditive_at_evaluation_point":
            bool(cons.interaction > 0),
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
