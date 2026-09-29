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
import landerer_s4 as L  # noqa: E402

LAW = Path("/storage/kiran-stuff/proteostasis_law")
EP = LAW / "envelope-paper"
P1 = Path("/storage/kiran-stuff/proteostasis-P1")
RAWD = EP / "data" / "raw"
sys.path.insert(0, str(EP / "scripts" / "vendor"))
import two_pool_ode as m  # noqa: E402  (legacy, read-only)

UPSTREAM = ("deTEL main @ c3593d46 (git.mpi-cbg.de/tothpetroczylab/detel, accessed 2026-09-29), "
            "eTEL/workflow/global_report.py get_all_codon_count: log10(detection_rate), "
            "inf -> NaN, dropna -> zero-detection codon/dataset cells are removed")

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
    (P1 / "two_pool_ode.py", 262, "(1.0 - p.S_avg)", "correct",
     "J_crit -> raw-error threshold", ""),
    (P1 / "two_pool_ode.py", 513, "(1.0 - p.S_avg)", "ambiguous",
     "f_codon = 1e-4 'observed E. coli rate'", "kind of 1e-4 window is unstated"),
    (P1 / "two_pool_ode.backup_uniformN.py", 236, "(1.0 - p.S_avg)", "correct",
     "J_crit -> raw-error threshold", "backup implementation"),
    (P1 / "two_pool_ode.backup_uniformN.py", 483, "(1.0 - p.S_avg)", "ambiguous",
     "f_codon = 1e-4 crosscheck", "kind of 1e-4 window is unstated"),
    (P1 / "paired_mc.py", 146, "(1.0 - p.S_avg)", "ambiguous",
     "f_obs = 1e-4", "same undocumented window value"),
    (P1 / "paired_mc.py", 148, "f_eff_obs = f_obs * (1.0 - S) * p_m", "ambiguous",
     "f_obs = 1e-4", "same undocumented window value"),
    (P1 / "arithmetic_stress_test.py", 72, "denom = (1.0 - S_syn) * p_misfold", "correct",
     "threshold on raw error", "exact raw threshold"),
    (P1 / "arithmetic_stress_test.py", 80, "denom = (1.0 - S_syn) * p_misfold", "correct",
     "threshold on raw error", "large-N approximation"),
    (P1 / "arithmetic_stress_test.py", 113, "denom = (1.0 - S_syn) * p_misfold", "correct",
     "threshold on raw error", ""),
    (P1 / "arithmetic_stress_test.py", 242,
     "f_cod = f_eff / ((1.0 - base.S_syn) * base.p_misfold)", "correct",
     "effective damage -> raw error", ""),
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
     "threshold on raw error", "same threshold form at lines 109, 121, 123"),
    (EP / "scripts/06_translation_burden.py", 37, "mubar = float((w * m.mu).sum())", "correct",
     "MS mu used directly as f_sub in 1-exp(-mu N)", "no (1-S); a lower bound, "
     "since MS undercounts"),
]


def g51():
    ec = L.load_s2("ecoli")
    n = (ec["sd"] / ec["se"]) ** 2
    mu = pd.read_csv(RAWD / "codon_error_rates_ecoli.tsv", sep="\t")
    m = mu.merge(ec, left_on="codon", right_on="Codon")
    cc, se = L.load_s4("ecoli")
    rec = L.reconstruct_s2(cc, ec)
    return {
        "landerer_definition_quote": LANDERER_DEF,
        "landerer_source": "Landerer, Poehls, Toth-Petroczy 2024 MBE 41:msae048, "
                           "doi:10.1093/molbev/msae048, PMC10939442",
        "supplement_zip": str(L.SUPPL_ZIP), "supplement_zip_sha256": L.sha256(L.SUPPL_ZIP),
        "legacy_data_S2_identical_to_supplement": L.sha256(RAWD / "Data_S2_error_detection_rate.xlsx")
        == L.sha256_member("Data_S2_error_detection_rate.xlsx"),
        "legacy_mu_equals_DataS2_mean": bool(np.allclose(m.mu, m["mean"], rtol=0, atol=0)),
        "n_datasets_per_codon_min": float(n.min()), "n_datasets_per_codon_max": float(n.max()),
        "n_max_dev_from_integer": float(np.max(np.abs(n - n.round()))),
        "n_datasets_total": int(cc.ds.nunique()),
        "cell_census": L.cell_census(cc, 80),
        "numerator_mismatches_codon_counts_vs_substitution_rows": L.numerator_check(cc, se),
        "total_substitution_psms": int(cc.error_count.sum()),
        "s2_reconstruction": rec,
        "n_meaning": "n = (sd/se)^2 = number of datasets with >= 1 detected substitution at "
                     "the codon (verified against Data_S4). it is NOT the number of covering "
                     "datasets: 3015 covered zero-error cells are excluded, only 19 cells lack "
                     "coverage.",
        "upstream_mechanism": UPSTREAM,
        "legacy_docstring_misdescribes": "02_axis_structure / MANUSCRIPT:202 call mu the "
                                         "'mean across detected substitutions'; it is a mean "
                                         "across detection-positive datasets of a "
                                         "destination-summed per-dataset rate",
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


# ---------------------------------------------------------------- G5.5
S_STIK = {"value": 1.82e-3, "se": 5.92e-5,
          "source": "Stikeleather, Ali, Ho, Licknack, Lynch 2026 NAR 54(13) gkag674, "
                    "doi:10.1093/nar/gkag674, PMID 42406629, PMC13335486",
          "definition": "total detected substitutions / total sites sampled (MS, I/L merged, "
                        "chemical-artefact substitutions removed); wild type, 3 replicates",
          "condition": "Xac-derived E. coli, LB (Miller), 37 C, overnight, harvested in "
                       "STATIONARY phase",
          "verification": "value 1.82e-3 and stationary-phase protocol read in PMC full text "
                          "2026-09-29; SE exponent lost in PMC text extraction (mantissa 5.92 "
                          "matches), exponent -5 taken from the task statement: UNVERIFIED"}


def headroom(J, p):
    P_dag, J_crit, mech, _ = m.saddle_node_operational(m.J_curve_two, m.A_qs, p)
    P_star, A_star = m.steady_state(J, p)
    return {"J": J, "P_star": P_star, "P_dagger": P_dag, "J_crit": J_crit, "mechanism": mech,
            "headroom_P": P_dag / P_star, "headroom_A": p.A_max / A_star}


def variants(f, p):
    """the same MS input through four mappings, isolating each correction."""
    ms = es.ErrorRate(f, es.MS)
    kw = dict(N_prot=p.N_prot, p_misfold=p.p_baseline, T_gen=p.T_gen_s, S=p.S_avg)
    J = {"legacy_(1-S)_1/T": es.flux_legacy(f, **kw),
         "no_(1-S)_1/T": es.flux(ms, **kw),
         "(1-S)_ln2/T": es.flux_legacy(f, **kw) * math.log(2.0),
         "no_(1-S)_ln2/T": es.flux(ms, balanced_growth=True, **kw)}
    return {k: headroom(v, p) for k, v in J.items()}


def g55():
    p = m.Params()                      # as_published anchoring: C_tot 50 uM, K_d 1 uM
    cc, _ = L.load_s4("ecoli")
    w = L.usage_weights()
    r = L.per_codon_rates(cc)
    ds = L.per_dataset_usage_weighted(cc, w)
    inputs = {
        "conditional_mean_DataS2": (L.usage_weighted(r.conditional_mean, w),
            "legacy mu: usage-weighted Data_S2 means; each codon mean excludes covered "
            "datasets with zero detections -> selected on detection, biased upward"),
        "zero_inclusive_mean": (L.usage_weighted(r.zero_inclusive_mean, w),
            "PRIMARY eTEL aggregate: per codon, mean of per-dataset detection rates over "
            "every dataset covering it (zeros kept); usage-weighted"),
        "pooled_ratio": (L.usage_weighted(r.pooled_ratio, w),
            "per codon sum(error PSM)/sum(covering PSM) across datasets; dominated by the "
            "deepest datasets; usage-weighted"),
        "dataset_median": (float(ds.usage_weighted.median()),
            "median over the 80 datasets of each dataset's own usage-weighted rate"),
    }
    out = {"model": "legacy vendored two_pool_ode (read-only), as_published anchoring; "
                    "legacy-model-conditional (Phi inflow, no dilution, A_max gate: WO-01/02/04)",
           "quantity_label": "MS error-DETECTION rate, amino-acid-substitution level, a lower "
                             "bound on f_sub (I/L invisible, rare events missed); NOT raw "
                             "decoding error; an eTEL aggregate over 80 heterogeneous PRIDE "
                             "datasets (strains, media, phases), NOT an E. coli physiological "
                             "operating point",
           "per_dataset_spread": {"min": float(ds.usage_weighted.min()),
                                  "q25": float(ds.usage_weighted.quantile(.25)),
                                  "median": float(ds.usage_weighted.median()),
                                  "q75": float(ds.usage_weighted.quantile(.75)),
                                  "max": float(ds.usage_weighted.max()),
                                  "spearman_rate_vs_psm_depth": float(
                                      ds[["usage_weighted", "psm"]].corr("spearman").iloc[0, 1])},
           "global_psm_pooled_all_codons": float(cc.error_count.sum() / cc.base_count.sum()),
           "inputs": {}}
    for k, (f, note) in inputs.items():
        v = variants(f, p)
        out["inputs"][k] = {"f": f, "note": note, "headroom": v}
    leg = out["inputs"]["conditional_mean_DataS2"]["headroom"]
    base = leg["legacy_(1-S)_1/T"]["headroom_P"]
    out["old_x"] = base
    out["effect_on_old_x"] = {
        "remove_(1-S)_only": leg["no_(1-S)_1/T"]["headroom_P"],
        "ln2_only": leg["(1-S)_ln2/T"]["headroom_P"],
        "both": leg["no_(1-S)_ln2/T"]["headroom_P"],
        "zero_inclusive_input_both": out["inputs"]["zero_inclusive_mean"]["headroom"]
                                     ["no_(1-S)_ln2/T"]["headroom_P"],
        "pooled_input_both": out["inputs"]["pooled_ratio"]["headroom"]
                             ["no_(1-S)_ln2/T"]["headroom_P"],
    }
    # critical substitution-level rate under the corrected mapping (no (1-S), ln2/T):
    # the f at which J reaches the legacy J_crit. J_crit is the imposed A_max gate
    # ('aggregation_death'), not a fold (WO-04).
    J_crit = leg["legacy_(1-S)_1/T"]["J_crit"]
    f_crit = J_crit * p.T_gen_s / (p.N_prot * p.p_baseline * math.log(2.0))
    out["f_crit_substitution_corrected"] = f_crit
    out["f_crit_raw_legacy"] = m.f_codon_from_J(J_crit, p)
    out["linear_f_margin_corrected"] = {k: f_crit / v[0] for k, v in inputs.items()}
    # corrected headroom across the per-dataset spread: how much the input choice matters
    out["headroom_P_corrected_across_datasets"] = {
        k: headroom(es.flux(es.ErrorRate(v, es.MS), p.N_prot, p.p_baseline, p.T_gen_s,
                            p.S_avg, balanced_growth=True), p)["headroom_P"]
        for k, v in out["per_dataset_spread"].items() if k != "spearman_rate_vs_psm_depth"}
    # stationary-phase estimate: reported, compared, NOT turned into an exponential headroom
    st = dict(S_STIK)
    st["ratio_to_zero_inclusive_eTEL"] = S_STIK["value"] / inputs["zero_inclusive_mean"][0]
    st["ratio_to_global_psm_pooled_eTEL"] = S_STIK["value"] / out["global_psm_pooled_all_codons"]
    st["MISMATCHED_what_if_f_crit_over_value"] = f_crit / S_STIK["value"]
    st["headroom"] = ("NOT COMPUTED: the legacy parameters (T_gen = 3600 s growth synthesis, "
                      "chaperone pool) describe exponential growth; a stationary-phase error "
                      "rate needs a stationary bundle (WO-07)")
    out["stikeleather_2026"] = st
    out["status"] = "COMPUTED"
    out["computed_headroom"] = True
    return out


def run():
    r = {"G5.1": g51(), "G5.2": g52(), "G5.4": g54(), "G5.5": g55()}
    r["pass"] = {
        "G5.1": bool(r["G5.1"]["legacy_mu_equals_DataS2_mean"]
                     and r["G5.1"]["legacy_data_S2_identical_to_supplement"]
                     and r["G5.1"]["s2_reconstruction"]["detected_only_error_count_gt_0"]
                        ["max_abs_dev_mean"] < 1e-8
                     and r["G5.1"]["s2_reconstruction"]["detected_only_error_count_gt_0"]
                        ["n_equals_sd_over_se_sq"]
                     and r["G5.1"]["numerator_mismatches_codon_counts_vs_substitution_rows"] == 0),
        "G5.2": all(x["anchor_found"] and x["class"] in ("correct", "double-discount", "ambiguous")
                    for x in r["G5.2"]),
        "G5.4": abs(r["G5.4"]["threshold_raw_exact"] - 5.6581428505e-3) < 1e-12
                and abs(r["G5.4"]["legacy_quoted_minus_lnP_over_N"] - 1.19e-3) / 1.19e-3 < 0.01,
        "G5.5": bool(r["G5.5"]["computed_headroom"]
                     and abs(r["G5.5"]["old_x"] - 24.817) < 5e-3
                     and all(v["headroom"]["no_(1-S)_ln2/T"]["J"] > 0
                             for v in r["G5.5"]["inputs"].values())),
    }
    return r


if __name__ == "__main__":
    r = run()
    print(json.dumps(r, indent=1, default=float))
    (HERE / "wo05_results.json").write_text(json.dumps(r, indent=2, default=float))
    if "BLOCKED" in r["pass"].values():
        sys.exit(2)
    sys.exit(0 if all(r["pass"].values()) else 1)
