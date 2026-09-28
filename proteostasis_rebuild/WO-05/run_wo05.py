#!/usr/bin/env python3
"""WO-05 gates G5.1-G5.5. writes wo05_results.json. legacy files are read-only."""
import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import error_semantics as es  # noqa: E402

LAW = Path("/storage/kiran-stuff/proteostasis_law")
EP = LAW / "envelope-paper"
P1 = Path("/storage/kiran-stuff/proteostasis-P1")
RAWD = EP / "data" / "raw"
sys.path.insert(0, str(EP / "scripts" / "vendor"))
import two_pool_ode as L  # noqa: E402  (read-only import)

# primary-source definition (PMC10939442, Methods, "Calculating Error Detection Rates")
LANDERER_DEF = ("For each codon, an error detection rate was calculated as the total number "
                "of times (number of Peptide Spectrum Match [PSM]) peptides covering a position "
                "with this codon and carrying a substitution at that position were detected "
                "divided by the total number of times peptides, modified or unmodified, covering "
                "a position with this codon were detected.")

# G5.2 every legacy error -> flux / threshold mapping. (file, line, anchor, class, input, note)
MAPPINGS = [
    (EP / "scripts/11_headroom_sensitivity.py", 81, "(1.0 - p.S_avg)", "double-discount",
     "MS mu (usage-weighted Data_S2 mean)", "the x25 headroom"),
    (EP / "scripts/12_chaperone_availability.py", 72, "(1.0 - p.S_avg)", "double-discount",
     "MS mu (translation_burden.json)", "theta sweep"),
    (EP / "scripts/09_supraadditivity.py", 99, "(1.0 - p.S_avg)", "double-discount",
     "MS mu (translation_burden.json)", "supraadditivity; margins also compared in f-space "
     "against f_codon_from_J, a raw-level threshold"),
    (EP / "scripts/vendor/two_pool_ode.py", 262, "(1.0 - p.S_avg)", "correct",
     "J_crit -> raw-error threshold", "correct as a raw threshold; comparing it to MS mu "
     "understates the margin by (1-S)"),
    (EP / "scripts/vendor/two_pool_ode.py", 513, "(1.0 - p.S_avg)", "ambiguous",
     "f_codon = 1e-4 'observed E. coli rate'", "the kind of the 1e-4 literature window "
     "is not stated"),
    (P1 / "paired_mc.py", 146, "(1.0 - p.S_avg)", "ambiguous",
     "f_obs = 1e-4", "same undocumented window value"),
    (P1 / "paired_mc.py", 148, "f_eff_obs = f_obs * (1.0 - S) * p_m", "ambiguous",
     "f_obs = 1e-4", "same undocumented window value"),
    (P1 / "arithmetic_stress_test.py", 72, "denom = (1.0 - S_syn) * p_misfold", "correct",
     "threshold on raw error", "exact raw threshold"),
    (P1 / "arithmetic_stress_test.py", 80, "denom = (1.0 - S_syn) * p_misfold", "correct",
     "threshold on raw error", "large-N approximation"),
    (P1 / "arithmetic_stress_test.py", 202, "(1-S)·p_m = 1.0", "ambiguous",
     "Part A point estimate", "factor forced to 1.0: neither a raw nor a substitution "
     "threshold under the stated model; source of 1.19e-3"),
    (P1 / "essential_bound.py", 209, "denom = (1.0 - S_syn) * p_misfold", "correct",
     "threshold on raw error", ""),
    (P1 / "essential_bound.py", 216, "f_eff = f_codon * (1.0 - S_syn) * p_misfold", "correct",
     "raw error -> effective damage", "correct if f_codon is raw"),
    (P1 / "essential_bound.py", 241, "denom = (1.0 - S_syn) * p_misfold", "correct",
     "threshold on raw error", ""),
    (P1 / "figures/fig2_arithmetic.py", 114, "((1.0 - Ssyn) * pmis)", "correct",
     "threshold on raw error", "figure anchor"),
    (EP / "scripts/06_translation_burden.py", 38, "mubar = float((w * m.mu).sum())", "correct",
     "MS mu used directly as f_sub in 1-exp(-mu N)", "no (1-S); a lower bound, "
     "since MS undercounts"),
]


def g51():
    ec = pd.read_excel(RAWD / "Data_S2_error_detection_rate.xlsx", sheet_name="E. coli")
    n = (ec["sd"] / ec["se"]) ** 2
    mu = pd.read_csv(RAWD / "codon_error_rates_ecoli.tsv", sep="\t")
    m = mu.merge(ec, left_on="codon", right_on="Codon")
    return {
        "landerer_definition_quote": LANDERER_DEF,
        "landerer_source": "Landerer, Poehls, Toth-Petroczy 2024 MBE 41:msae048, PMC10939442",
        "legacy_mu_equals_DataS2_mean": bool(np.allclose(m.mu, m["mean"], rtol=0, atol=0)),
        "n_datasets_per_codon_min": float(n.min()), "n_datasets_per_codon_max": float(n.max()),
        "n_max_dev_from_integer": float(np.max(np.abs(n - n.round()))),
        "n_datasets_total": 80,
        "inference": "n < 80 for every codon, so Data_S2 mean is over datasets with >= 1 "
                     "detected substitution at that codon (conditional mean)",
        "legacy_docstring_misdescribes": "02_axis_structure / MANUSCRIPT:202 call mu the "
                                         "'mean across detected substitutions'; it is a mean "
                                         "across datasets of a destination-summed rate",
    }


def g52():
    rows = []
    for f, line, anchor, cls, inp, note in MAPPINGS:
        text = f.read_text().splitlines()
        ok = anchor in text[line - 1]
        rows.append({"file": str(f), "line": line, "anchor": anchor, "anchor_found": ok,
                     "class": cls, "input": inp, "note": note})
    return rows


def g54():
    N, P, S, pm = 300.0, 0.70, 0.30, 0.30
    return {"threshold_raw_exact": es.threshold_raw(N, P, S, pm),
            "threshold_raw_largeN": -math.log(P) / (N * (1 - S) * pm),
            "threshold_sub_exact": es.threshold_sub(N, P, pm),
            "legacy_quoted_minus_lnP_over_N": es.threshold_legacy_quoted(N, P),
            "unfiltered_exact_1_minus_P_pow": 1.0 - P ** (1.0 / N),
            "explanation": "1.19e-3 = -ln(0.7)/300: arithmetic_stress_test.py Part A sets "
                           "(1-S)p_m = 1.0 (line 202), dropping both factors"}


def usage_weighted(values):
    usage = pd.read_csv(RAWD / "global_codon_usage_ecoli.tsv", sep="\t")
    m = usage.merge(values, on="codon")
    w = m["count"] / m["count"].sum()
    return float((w * m.val).sum())


ANCHORINGS = [("as_published", 50.0, 1.0), ("weaker_binding", 50.0, 10.0),
              ("smaller_pool", 5.0, 1.0), ("near_capacity", 2.0, 1.0),
              ("c_free_at_Kd", 1.0, 1.0), ("Kd_at_C_tot", 50.0, 50.0)]


def headroom(J, Ct, Kd):
    p = L.Params()
    p.C_tot_uM, p.K_d_uM = Ct, Kd
    P_dag, _, mech, _ = L.saddle_node_operational(L.J_curve_two, L.A_qs, p)
    P_star, _ = L.steady_state(J, p)
    return (P_dag / P_star if np.isfinite(P_star) else float("nan")), mech


def g55():
    mu = pd.read_csv(RAWD / "codon_error_rates_ecoli.tsv", sep="\t").rename(columns={"mu": "val"})
    mubar = usage_weighted(mu)
    p = L.Params()
    r = es.ErrorRate(mubar, es.MS)
    J_leg = es.flux_legacy(mubar, p.N_prot, p.p_baseline, p.T_gen_s, p.S_avg)
    J_cor = es.flux(r, p.N_prot, p.p_baseline, p.T_gen_s, p.S_avg)
    rows = []
    for name, Ct, Kd in ANCHORINGS:
        h_leg, mech = headroom(J_leg, Ct, Kd)
        h_cor, _ = headroom(J_cor, Ct, Kd)
        rows.append({"anchoring": name, "mechanism": mech,
                     "headroom_P_legacy": h_leg, "headroom_P_corrected": h_cor})
    # disclosed sensitivity, not a replacement: zero-inclusive mean, ASSUMING every
    # codon was covered in all 80 datasets (not verifiable from Data_S2)
    ec = pd.read_excel(RAWD / "Data_S2_error_detection_rate.xlsx", sheet_name="E. coli")
    n = ((ec["sd"] / ec["se"]) ** 2).round()
    zi = pd.DataFrame({"codon": ec.Codon, "val": ec["mean"] * n / 80.0})
    mubar_zi = usage_weighted(zi)
    J_zi = es.flux(es.ErrorRate(mubar_zi, es.MS), p.N_prot, p.p_baseline, p.T_gen_s, p.S_avg)
    h_zi, _ = headroom(J_zi, 50.0, 1.0)
    return {"mubar_usage_weighted": mubar, "J_legacy": J_leg, "J_corrected": J_cor,
            "J_ratio": J_cor / J_leg, "by_anchoring": rows,
            "sensitivity_zero_inclusive_mean": {
                "mubar": mubar_zi, "headroom_P_corrected_as_published": h_zi,
                "status": "ASSUMPTION: all 80 datasets cover every codon; not a headline"},
            "caveats": ["MS undercounts (I/L invisible, rare events missed): the corrected "
                        "flux is a lower bound, so the corrected headroom is an UPPER bound",
                        "P_dagger at as_published is the imposed A_max gate (WO-04 G4.2), "
                        "not a fold",
                        "Landerer datasets pool growth conditions; phase matching is WO-07"]}


def run():
    r = {"G5.1": g51(), "G5.2": g52(), "G5.4": g54(), "G5.5": g55()}
    pub = [x for x in r["G5.5"]["by_anchoring"] if x["anchoring"] == "as_published"][0]
    r["pass"] = {
        "G5.1": bool(r["G5.1"]["legacy_mu_equals_DataS2_mean"]
                     and r["G5.1"]["n_max_dev_from_integer"] < 1e-3
                     and r["G5.1"]["n_datasets_per_codon_max"] < 80),
        "G5.2": all(x["anchor_found"] and x["class"] in ("correct", "double-discount", "ambiguous")
                    for x in r["G5.2"]),
        "G5.4": abs(r["G5.4"]["threshold_raw_exact"] - 5.6581428505e-3) < 1e-12
                and abs(r["G5.4"]["legacy_quoted_minus_lnP_over_N"] - 1.19e-3) / 1.19e-3 < 0.01,
        "G5.5": bool(abs(pub["headroom_P_legacy"] - 24.817281807890474) < 1e-6
                     and np.isfinite(pub["headroom_P_corrected"])
                     and pub["headroom_P_corrected"] < pub["headroom_P_legacy"]),
    }
    return r


if __name__ == "__main__":
    r = run()
    print(json.dumps(r, indent=1, default=float))
    (HERE / "wo05_results.json").write_text(json.dumps(r, indent=2, default=float))
    sys.exit(0 if all(r["pass"].values()) else 1)
