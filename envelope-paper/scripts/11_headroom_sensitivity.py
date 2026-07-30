#!/usr/bin/env python3
"""
how far inside the envelope does E. coli actually sit?

this script exists because the headline number was not robust. earlier drafts
reported x158 headroom, obtained by evaluating the two-pool model at
f = 1e-4 /codon -- the BOTTOM of the quoted observed window [1e-4, 1e-3]. two
things are wrong with that.

1. WRONG EVALUATION POINT. the paper computes its own usage-weighted mean
   per-codon error rate from the Landerer data and the E. coli codon usage:
   6.33e-4 /codon (scripts/06). evaluating the headroom at 1e-4 while deriving
   6.33e-4 from the same data is inconsistent, and it is inconsistent in the
   favourable direction -- exactly the error the earlier draft made when it quoted
   the arithmetic bound at its deterministic reference point (1.19e-3) rather than
   its Monte Carlo median (8.35e-3).

2. UNTESTED CHAPERONE ANCHORING. at the observed operating point the model's
   folding arm is 97.9% saturated: 47.5 uM free chaperone against 0.052 uM
   misfolded protein (scripts/09). the capacity evidence the paper itself cites --
   Hsp70 buffering of production costs, and metastable-protein interference --
   describes a network running NEAR capacity, not in ~900-fold excess. so the
   published anchoring is not the only defensible one, and the headroom should be
   reported across the alternatives.

this script crosses the two axes and reports the range. it does not pick a
flattering cell.
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


# ---- axis 1: where to evaluate the observed error rate --------------------
def evaluation_points():
    b = json.loads((COMP / "translation_burden.json").read_text())
    return [
        ("window_bottom", 1e-4,
         "bottom of the quoted observed window; what earlier drafts used"),
        ("usage_weighted_mu", b["usage_weighted_mean_mu_per_codon"],
         "this paper's own mu data weighted by E. coli codon usage -- the "
         "internally consistent choice"),
        ("unweighted_mean_mu", b["unweighted_mean_mu_per_codon"],
         "unweighted mean over the 61 sense codons"),
        ("window_top", 1e-3, "top of the quoted observed window"),
    ]


# ---- axis 2: how to anchor the chaperone arm ------------------------------
# each entry is (label, C_tot_uM, K_d_uM, justification)
ANCHORINGS = [
    ("as_published", 50.0, 1.0,
     "upstream baseline; leaves the folding arm ~98% saturated"),
    ("weaker_binding", 50.0, 10.0,
     "same pool, weaker chaperone-substrate affinity"),
    ("smaller_pool", 5.0, 1.0,
     "10-fold smaller effective chaperone pool"),
    ("near_capacity", 2.0, 1.0,
     "free chaperone a few-fold above K_d -- the regime the cited capacity "
     "evidence describes"),
    ("c_free_at_Kd", 1.0, 1.0,
     "free chaperone comparable to K_d: the network genuinely at capacity"),
    ("Kd_at_C_tot", 50.0, 50.0,
     "affinity comparable to pool size; another route to an unsaturated arm"),
]


def evaluate(f_codon, C_tot, K_d):
    p = m.Params()
    p.C_tot_uM, p.K_d_uM = C_tot, K_d
    P_dag, J_crit, mech, _ = m.saddle_node_operational(m.J_curve_two, m.A_qs, p)
    J = f_codon * p.N_prot * (1.0 - p.S_avg) * p.p_baseline / p.T_gen_s
    P_star, A_star = m.steady_state(J, p)
    if not np.isfinite(P_star):
        return None
    return {
        "P_star": P_star, "A_star": A_star, "P_dagger": P_dag,
        "headroom_P": P_dag / P_star, "headroom_A": p.A_max / A_star,
        "margin_log10": float(np.log10(min(P_dag / P_star, p.A_max / A_star))),
        "c_free_uM": float(m.c_free(P_star, p)),
        "folding_arm_saturation": float(m.v_fold(P_star, p) / p.k_obs_max),
        "mechanism": mech,
    }


def main():
    rows = []
    for ep_key, f, ep_note in evaluation_points():
        for anc_key, Ct, Kd, anc_note in ANCHORINGS:
            r = evaluate(f, Ct, Kd)
            base = {"evaluation_point": ep_key, "f_codon": f,
                    "anchoring": anc_key, "C_tot_uM": Ct, "K_d_uM": Kd}
            if r is None:
                rows.append({**base, "headroom_P": np.nan, "headroom_A": np.nan,
                             "margin_log10": np.nan, "collapsed": True,
                             "c_free_uM": np.nan, "folding_arm_saturation": np.nan})
            else:
                rows.append({**base, **r, "collapsed": False})
    df = pd.DataFrame(rows)
    df.to_csv(COMP / "headroom_sensitivity.tsv", sep="\t", index=False)

    print("=" * 84)
    print("headroom in the misfolded-monomer pool, across both axes")
    print("=" * 84)
    piv = df.pivot_table(index="evaluation_point", columns="anchoring",
                         values="headroom_P")
    order_e = [k for k, _, _ in evaluation_points()]
    order_a = [k for k, _, _, _ in ANCHORINGS]
    piv = piv.reindex(index=order_e, columns=order_a)
    print(piv.round(1).to_string())

    published = df[(df.evaluation_point == "window_bottom")
                   & (df.anchoring == "as_published")].iloc[0]
    central = df[(df.evaluation_point == "usage_weighted_mu")
                 & (df.anchoring == "as_published")].iloc[0]
    ok = df[~df.collapsed]
    lo, hi = ok.headroom_P.min(), ok.headroom_P.max()
    # the internally consistent evaluation point, across all anchorings
    uw = ok[ok.evaluation_point == "usage_weighted_mu"]

    print(f"\n  as previously reported   (window bottom, as published): "
          f"x{published.headroom_P:.0f}   margin {published.margin_log10:.2f} log10")
    print(f"  internally consistent    (usage-weighted mu, as published): "
          f"x{central.headroom_P:.1f}   margin {central.margin_log10:.2f} log10")
    print(f"  usage-weighted mu, across all six anchorings: "
          f"x{uw.headroom_P.min():.1f} to x{uw.headroom_P.max():.1f}")
    print(f"  full grid range: x{lo:.1f} to x{hi:.0f}  "
          f"({np.log10(lo):.2f} to {np.log10(hi):.2f} log10)")
    print(f"\n  so 'roughly two orders of magnitude inside' holds only at the "
          f"favourable\n  corner of the grid. the defensible statement is "
          f"roughly one order, and at\n  the internally consistent evaluation "
          f"point about {np.log10(central.headroom_P):.1f} orders.")

    # where does this put E. coli relative to the supraadditivity onset?
    sup = json.loads((COMP / "supraadditivity_summary.json").read_text())
    onset = sup["qualitative_supraadditivity_first_seen_at_margin_log10"]
    print(f"\n  joint collapse of 3-fold perturbations begins at margin "
          f"{onset:.2f} log10.")
    print(f"  at the internally consistent evaluation point the margin is "
          f"{central.margin_log10:.2f} log10,")
    print(f"  i.e. within {central.margin_log10 - onset:+.2f} log10 of that "
          f"onset -- the prediction is reachable,")
    print(f"  not remote.")

    summary = {
        "previously_reported_headroom_P": float(published.headroom_P),
        "previously_reported_evaluation_point": "f = 1e-4 (window bottom)",
        "why_that_was_wrong": "evaluated at the bottom of the observed window "
                              "while deriving 6.33e-4 from the same data, and "
                              "with the folding arm 98% saturated",
        "internally_consistent_headroom_P": float(central.headroom_P),
        "internally_consistent_headroom_A": float(central.headroom_A),
        "internally_consistent_margin_log10": float(central.margin_log10),
        "internally_consistent_evaluation_point_f": float(central.f_codon),
        "usage_weighted_mu_range_across_anchorings":
            [float(uw.headroom_P.min()), float(uw.headroom_P.max())],
        "full_grid_range_headroom_P": [float(lo), float(hi)],
        "full_grid_range_orders_of_magnitude":
            [float(np.log10(lo)), float(np.log10(hi))],
        "supraadditivity_onset_margin_log10": onset,
        "distance_from_onset_log10":
            float(central.margin_log10 - onset),
        "n_collapsed_cells": int(df.collapsed.sum()),
        "n_cells": int(len(df)),
    }
    (COMP / "headroom_sensitivity_summary.json").write_text(
        json.dumps(summary, indent=2))
    print(f"\nwrote headroom_sensitivity.tsv and "
          f"headroom_sensitivity_summary.json to {COMP}")


if __name__ == "__main__":
    main()
