#!/usr/bin/env python3
"""
the chaperone arm's saturation is a hidden assumption. make it a parameter.

what the previous pass got wrong
--------------------------------
scripts/09 found the folding arm 97.9% saturated at the operating point and I
attributed that to the parameterization contradicting the capacity evidence the
paper cites. reading the upstream anchors (proteostasis-P1/LITERATURE_ANCHORS.md)
shows that is not right. both values are properly sourced and used within range:

    C_tot     30-80 uM, baseline 50    GroEL ~30 uM in exponentially growing
                                       E. coli (Lorimer 1996); DnaK 30-50 uM
    K_d       0.06-2 uM, baseline 1    Pierpaoli et al. 1997 EMBO J, DnaK-substrate
                                       peptide dissociation constants
    k_obs_max 3e-3-8.4e-2 /s, base 1e-2  Pierpaoli 1997, DnaK R-state release

so the saturation is not a bad number. it is structural. the rescue term is

    c_free = C_tot / (1 + M/K_d),   M = P * Prot_tot

and at steady state the damaged pool M is 0.05 uM against K_d = 1 uM, so
c_free ~ C_tot and the arm sits at ~98% of k_obs_max. the model hands the WHOLE
chaperone pool to the damaged-protein pool, because it does not represent the
ordinary nascent-chain folding load that occupies most chaperone capacity in a
growing cell at all. C_tot in this model is therefore not total cellular
chaperone; it is chaperone *available to the damaged pool*, and setting it to the
total silently assumes that availability is 100%.

what this script does
---------------------
introduces that assumption as an explicit parameter -- the fraction of the
chaperone pool already committed elsewhere,

    C_avail = C_tot * (1 - theta)

and sweeps theta together with C_tot and K_d over their DOCUMENTED ranges. this
replaces the six ad-hoc anchorings of scripts/11 with one interpretable axis and
makes the follow-up measurement explicit: constrain theta and the headroom range
collapses to a point.

theta is not measured here. nothing below should be read as an estimate of it.
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

# documented literature ranges (proteostasis-P1/LITERATURE_ANCHORS.md)
C_TOT_RANGE = [30.0, 50.0, 80.0]      # uM   Lorimer 1996; DnaK levels
K_D_RANGE = [0.06, 1.0, 2.0]          # uM   Pierpaoli et al. 1997
THETA = [0.0, 0.5, 0.8, 0.9, 0.95, 0.98, 0.99]

SATURATED = 0.5   # below this the folding arm is no longer saturated (c_free < K_d)


def evaluate(f_codon, C_tot, K_d, theta):
    p = m.Params()
    p.C_tot_uM = C_tot * (1.0 - theta)
    p.K_d_uM = K_d
    if p.C_tot_uM <= 0:
        return None
    P_dag, J_crit, mech, _ = m.saddle_node_operational(m.J_curve_two, m.A_qs, p)
    J = f_codon * p.N_prot * (1.0 - p.S_avg) * p.p_baseline / p.T_gen_s
    P_star, A_star = m.steady_state(J, p)
    if not np.isfinite(P_star):
        return None
    return {
        "C_avail_uM": p.C_tot_uM,
        "c_free_uM": float(m.c_free(P_star, p)),
        "folding_arm_saturation": float(m.v_fold(P_star, p) / p.k_obs_max),
        "P_star": P_star, "P_dagger": P_dag,
        "headroom_P": P_dag / P_star,
        "headroom_A": p.A_max / A_star,
        "margin_log10": float(np.log10(min(P_dag / P_star, p.A_max / A_star))),
        "mechanism": mech,
    }


def main():
    burden = json.loads((COMP / "translation_burden.json").read_text())
    f = burden["usage_weighted_mean_mu_per_codon"]
    sup = json.loads((COMP / "supraadditivity_summary.json").read_text())
    onset = sup["qualitative_supraadditivity_first_seen_at_margin_log10"]

    rows = []
    for Ct in C_TOT_RANGE:
        for Kd in K_D_RANGE:
            for th in THETA:
                r = evaluate(f, Ct, Kd, th)
                base = {"C_tot_uM": Ct, "K_d_uM": Kd, "theta": th, "f_codon": f}
                rows.append({**base, **r, "collapsed": False} if r else
                            {**base, "collapsed": True, "headroom_P": np.nan,
                             "margin_log10": np.nan,
                             "folding_arm_saturation": np.nan,
                             "c_free_uM": np.nan, "C_avail_uM": Ct * (1 - th)})
    df = pd.DataFrame(rows)
    df.to_csv(COMP / "chaperone_availability.tsv", sep="\t", index=False)

    print("=" * 80)
    print(f"headroom at the usage-weighted mean error rate ({f:.2e} /codon)")
    print("theta = fraction of the chaperone pool committed to other load")
    print("=" * 80)
    piv = df[df.K_d_uM == 1.0].pivot_table(index="C_tot_uM", columns="theta",
                                           values="headroom_P")
    print("\n  at K_d = 1 uM (baseline):")
    print(piv.round(1).to_string())
    sat = df[df.K_d_uM == 1.0].pivot_table(index="C_tot_uM", columns="theta",
                                           values="folding_arm_saturation")
    print("\n  folding-arm saturation (v_fold / k_obs_max):")
    print(sat.round(3).to_string())

    ok = df[~df.collapsed]
    base = ok[(ok.C_tot_uM == 50.0) & (ok.K_d_uM == 1.0)]
    theta0 = base[base.theta == 0.0].iloc[0]

    # where does the arm stop being saturated, and where does the margin reach
    # the supraadditivity onset?
    unsat = base[base.folding_arm_saturation < SATURATED]
    theta_unsat = float(unsat.theta.min()) if len(unsat) else None
    crossed = base[base.margin_log10 <= onset]
    theta_onset = float(crossed.theta.min()) if len(crossed) else None

    print(f"\n  at theta = 0 (the implicit assumption): "
          f"saturation {theta0.folding_arm_saturation:.3f}, "
          f"headroom x{theta0.headroom_P:.1f}, margin "
          f"{theta0.margin_log10:.2f} log10")
    if theta_unsat is not None:
        u = base[base.theta == theta_unsat].iloc[0]
        print(f"  the arm stops being saturated at theta >= {theta_unsat:.2f} "
              f"(c_free {u.c_free_uM:.2f} uM < K_d): headroom x{u.headroom_P:.1f}, "
              f"margin {u.margin_log10:.2f}")
    if theta_onset is not None:
        c = base[base.theta == theta_onset].iloc[0]
        print(f"  the margin reaches the supraadditivity onset "
              f"({onset:.2f}) at theta >= {theta_onset:.2f}: "
              f"headroom x{c.headroom_P:.1f}")
    else:
        print(f"  the margin never reaches the supraadditivity onset "
              f"({onset:.2f}) over the theta grid")

    print(f"\n  headroom across the full documented grid: "
          f"x{ok.headroom_P.min():.1f} to x{ok.headroom_P.max():.1f}")
    print(f"  ({int(df.collapsed.sum())} of {len(df)} cells have no stable state)")

    summary = {
        "evaluation_point_f_codon": f,
        "parameter_provenance": {
            "C_tot_uM": "30-80, baseline 50 (Lorimer 1996 GroEL ~30 uM; DnaK "
                        "30-50 uM) -- proteostasis-P1/LITERATURE_ANCHORS.md",
            "K_d_uM": "0.06-2, baseline 1 (Pierpaoli et al. 1997 EMBO J)",
            "theta": "NOT measured; introduced here to expose the model's "
                     "implicit assumption that 100% of the chaperone pool is "
                     "available to the damaged-protein pool",
        },
        "correction_to_earlier_claim":
            "scripts/09 attributed the 98%-saturated folding arm to a "
            "parameterization inconsistent with the paper's own capacity "
            "evidence. That was wrong: C_tot and K_d are properly sourced and "
            "used in range. The saturation is structural -- the model omits the "
            "nascent-chain folding load, so it hands the entire chaperone pool "
            "to the damaged pool.",
        "at_theta_zero": {
            "folding_arm_saturation": float(theta0.folding_arm_saturation),
            "headroom_P": float(theta0.headroom_P),
            "margin_log10": float(theta0.margin_log10),
        },
        "theta_at_which_arm_unsaturates": theta_unsat,
        "theta_at_which_margin_reaches_supraadditivity_onset": theta_onset,
        "supraadditivity_onset_margin_log10": onset,
        "headroom_range_over_documented_grid":
            [float(ok.headroom_P.min()), float(ok.headroom_P.max())],
        "n_collapsed": int(df.collapsed.sum()),
        "n_cells": int(len(df)),
        "what_would_pin_this_down":
            "a measurement of chaperone occupancy by nascent-chain folding in "
            "exponentially growing E. coli would fix theta and collapse the "
            "headroom range to a point",
    }
    (COMP / "chaperone_availability_summary.json").write_text(
        json.dumps(summary, indent=2))
    print(f"\nwrote chaperone_availability.tsv and "
          f"chaperone_availability_summary.json to {COMP}")


if __name__ == "__main__":
    main()
