"""WO-07: write bundles.tsv, the EXPONENTIAL and STATIONARY parameter bundles.

every value is read from a saved file (WO-06 audit table / audit_derived.json /
schmidt_pools.json, WO-07 derived.json). nothing is typed in except literal
in-vitro constants, and check_bundles.py confirms each of those occurs in the
cited WO-06 audit row.

one row per (bundle, param, candidate). per (bundle, param) there is either a
MATCHED row or an explicit UNMEASURED row (value NA). MISMATCHED candidates may
sit next to an UNMEASURED row; they never stand in for it.

anchors (the condition each bundle is matched TO):
  EXPONENTIAL  E. coli BW25113, M9 glucose, 37 C, balanced exponential growth
               (Schmidt 2016 glucose: the one condition whose total protein is
               an independent measurement)
  STATIONARY   E. coli Xac (the paper's wild type), LB (Miller), 37 C, overnight
               stationary (Stikeleather 2026: the stationary error input defines it)

two error quantities are kept apart (reviewer finding): what MS measures on a
harvested culture is the substitution FREQUENCY per codon in the STANDING
proteome (substitution_freq_standing_f); what the burden flux needs is the
error rate per codon SYNTHESIZED (synthesis_error_rate_per_codon). they coincide
only in balanced growth with no selective degradation; in a stationary harvest
most chains were made before growth stopped. the second is UNMEASURED in both.
"""
import csv
import json
import sys
from pathlib import Path

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
W6 = ROOT / "WO-06"

EXP, STAT = "EXPONENTIAL", "STATIONARY"
ANCHOR = {EXP: {"strain": "BW25113", "medium": "M9 glucose", "temperature": "37 C"},
          STAT: {"strain": "Xac", "medium": "LB (Miller)", "temperature": "37 C"}}
EXP_COND, STAT_COND, STAT_COND3 = "Glucose", "Stationary phase 1 day", "Stationary phase 3 days"
SCHMIDT_STAT_MEDIUM = "M9 glucose, then starved"   # 'stationary phase after a glucose culture'

COLUMNS = ["bundle", "param", "machine", "candidate", "value", "unit", "condition_tag",
           "value_phase", "organism", "strain", "medium", "temperature", "source",
           "source_key", "measurement_type", "wo06_status", "match_flag", "mismatch_axes",
           "effect_ids", "note"]

# every param each bundle must account for (MATCHED or UNMEASURED)
MACHINES = [("pool_DnaK", "DnaK", "dnaK", "uM DnaK monomer"),
            ("pool_GroEL14", "GroEL", "groL", "uM GroEL14 (tetradecamer)"),
            ("pool_GroES7", "GroES", "groS", "uM GroES7 (heptamer)"),
            ("pool_ClpB6", "ClpB", "clpB", "uM ClpB6 (hexamer)"),
            ("pool_DnaJ2", "DnaJ", "dnaJ", "uM DnaJ2 (dimer)"),
            ("pool_GrpE2", "GrpE", "grpE", "uM GrpE2 (dimer)"),
            ("pool_TF", "TF", "tig", "uM trigger factor monomer"),
            ("pool_HtpG2", "HtpG", "htpG", "uM HtpG2 (dimer)")]
UNMEASURED_BOTH = {
    "synthesis_error_rate_per_codon": ("1/codon synthesized", "-", "MS measures substitution "
                                       "frequency in the standing proteome; the per-synthesis "
                                       "rate needs the synthesis epoch and selective "
                                       "degradation of erroneous chains, neither measured"),
    "GroEL_cycle_rate": ("1/s", "GroEL", "no GroEL cycle or folding rate in any saved record"),
    "ClpB_disaggregation_rate_k_dis": ("1/s", "ClpB", "L11: heat-shock recovery with induced "
                                       "chaperones; the 95%/2 h input is not located"),
    "misfolded_degradation_k_d": ("1/s", "-", "L01c Goldberg 1972: E. coli degrades abnormal "
                                  "proteins faster, no rate constant located; L01a/b are yeast"),
    "aggregation_k_a": ("1/(uM s)", "-", "L05a NOT_SUPPORTED, L05b UNVERIFIED: no second-order "
                        "constant"),
    "misfold_probability_p_misfold": ("1", "-", "L09 MISCITED; the 0.1-0.5 range is loss of "
                                      "function from random mutations, not misfolding per "
                                      "mistranslation"),
    "cycle_partition_phi": ("1", "DnaK", "no saved record measures the productive fraction "
                            "of chaperone cycles"),
    "unfolding_k_mis": ("1/s", "-", "no saved record"),
    "aggregate_degradation_k_dA": ("1/s", "-", "no saved record"),
    "legacy_A_max": ("1 (fraction aggregated)", "-", "L13a-e: 2 MISCITED, 2 UNMATCHED; ASSUMED"),
}
PARAMS = (["growth_rate_mu", "synthesis_rate_s_P", "total_protein_P_T"]
          + [p for p, *_ in MACHINES]
          + ["copy_weighted_length_N", "codon_usage_weights", "substitution_freq_standing_f",
             "DnaK_cycle_rate_k_cat", "DnaK_affinity_K"] + list(UNMEASURED_BOTH))


def g(x):
    return f"{x:.6g}"


def load_sources():
    sch = json.loads((W6 / "schmidt_pools.json").read_text())
    per = {c["condition"]: c for c in sch["per_condition"]}
    der = json.loads((HERE / "derived.json").read_text())
    ad = json.loads((W6 / "audit_derived.json").read_text())
    return per, der, ad


def row(bundle, param, value, unit, value_phase, organism, strain, medium, temperature,
        source, source_key, mtype, wo06, flag, axes="-", effects="-", note="-",
        machine="-", candidate="primary"):
    return {"bundle": bundle, "param": param, "machine": machine, "candidate": candidate,
            "value": value, "unit": unit,
            "condition_tag": f"phase={value_phase};strain={strain};medium={medium};T={temperature}",
            "value_phase": value_phase, "organism": organism, "strain": strain,
            "medium": medium, "temperature": temperature, "source": source,
            "source_key": source_key, "measurement_type": mtype, "wo06_status": wo06,
            "match_flag": flag, "mismatch_axes": axes, "effect_ids": effects, "note": note}


def unmeasured(bundle, param, unit, why, machine="-"):
    return row(bundle, param, "NA", unit, "NONE", "Escherichia coli", "NA", "NA", "NA",
               "NONE", "none", "NONE", "-", "UNMEASURED", note=why, machine=machine,
               candidate="anchor")


SCHMIDT = "Schmidt 2016 Nat Biotechnol 34:104, PMID:26641532, Tables S6/S23"
ECO = "Escherichia coli"


def schmidt_row(bundle, param, cond, field, unit, flag, axes="-", effects="-", note="-",
                machine="-", scale=1.0, candidate="primary"):
    per, _, _ = load_sources()
    ph = STAT if per[cond]["phase"] == "stationary" else EXP
    medium = SCHMIDT_STAT_MEDIUM if ph == STAT else ANCHOR[EXP]["medium"]
    return row(bundle, param, g(per[cond][field] * scale), unit, ph, ECO, "BW25113", medium,
               "37 C", SCHMIDT, f"schmidt:{cond}:{field}", "DERIVED_FROM_MEASUREMENT",
               "VERIFIED", flag, axes, effects, note, machine, candidate)


def invitro_rows(bundle, sfx):
    """DnaK kinetic constants: in vitro in every bundle, so MISMATCHED in both."""
    return [
        unmeasured(bundle, "DnaK_cycle_rate_k_cat", "1/s", "no in-vivo DnaK cycle rate in any "
                   "saved record", "DnaK"),
        row(bundle, "DnaK_cycle_rate_k_cat", "0.04", "1/s", "IN_VITRO",
            "Escherichia coli (purified DnaK, DnaJ, GrpE)", "NA (in vitro)",
            "buffer ('conditions approximating those in the cell')", "UNSPECIFIED_IN_RECORD",
            "Pierpaoli 1997 J Mol Biol 269:757, PMID:9223639", "audit:L02b", "MEASURED",
            "MISMATCHED_CONDITION", "MISMATCHED",
            "IN_VITRO;TEMPERATURE_UNSPECIFIED;QUANTITY", f"{sfx}_KCAT",
            "rate-limiting DnaJ-triggered T->R step with peptides; a cycle step, not a "
            "per-client cycle rate; DnaJ:DnaK stoichiometry in vivo is ~1:27 (S04)",
            "DnaK", "T_to_R_step"),
        row(bundle, "DnaK_cycle_rate_k_cat", "0.003..0.084", "1/s", "IN_VITRO",
            "Escherichia coli (purified DnaK)", "NA (in vitro)", "buffer pH 7.0", "25 C",
            "Pierpaoli 1998 Biochemistry 37:16741, PMID:9843444", "audit:L02", "MEASURED",
            "MISCITED", "MISMATCHED", "IN_VITRO;TEMPERATURE;QUANTITY", f"{sfx}_KCAT",
            "legacy k_obs_max range: a peptide BINDING rate of nucleotide-free DnaK; WO-06 "
            "MISCITED refers to the legacy citation, the values are in the 1998 paper",
            "DnaK", "legacy_k_obs_max"),
        unmeasured(bundle, "DnaK_affinity_K", "uM", "the model needs the cycle's K_M "
                   "(WO-03), not a nucleotide-state K_d; no in-vivo value", "DnaK"),
        row(bundle, "DnaK_affinity_K", "0.06..2", "uM", "IN_VITRO",
            "Escherichia coli (purified DnaK)", "NA (in vitro)", "buffer pH 7.0", "25 C",
            "Pierpaoli 1998 Biochemistry 37:16741, PMID:9843444", "audit:L04", "MEASURED",
            "MISCITED", "MISMATCHED", "IN_VITRO;TEMPERATURE;QUANTITY", f"{sfx}_K",
            "nucleotide-free R-state K_d, short peptides", "DnaK", "R_state_Kd"),
        row(bundle, "DnaK_affinity_K", "2.2..107", "uM", "IN_VITRO",
            "Escherichia coli (purified DnaK)", "NA (in vitro)", "buffer pH 7.0", "25 C",
            "Pierpaoli 1998 Biochemistry 37:16741, PMID:9843444", "audit:L04", "MEASURED",
            "MISCITED", "MISMATCHED", "IN_VITRO;TEMPERATURE;QUANTITY", f"{sfx}_K",
            "DnaK-ATP T-state K_d, short peptides (dominant state with ATP present)",
            "DnaK", "T_state_Kd"),
    ]


def exponential(per, der, ad):
    b = EXP
    gl = per[EXP_COND]
    rows = [
        schmidt_row(b, "growth_rate_mu", EXP_COND, "growth_per_h", "1/s", "MATCHED",
                    scale=1 / 3600, note="0.58 /h; S06 VERIFIED"),
        row(b, "synthesis_rate_s_P", g(gl["growth_per_h"] / 3600 * gl["total_protein_mM_wholecell"]
                                       * 1e3), "uM chains/s", EXP, ECO, "BW25113", "M9 glucose",
            "37 C", SCHMIDT, "schmidt:Glucose:growth_per_h*total_protein_mM_wholecell",
            "DERIVED_FROM_MEASUREMENT", "VERIFIED", "MATCHED",
            note="s_P = mu * P_T: balanced-growth closure, turnover of native protein "
                 "neglected; valid only in exponential growth"),
        schmidt_row(b, "total_protein_P_T", EXP_COND, "total_protein_mM_wholecell", "uM chains",
                    "MATCHED", scale=1e3, note="L06d; whole-cell volume (periplasm included); "
                    "the only independent total in Schmidt"),
        row(b, "total_protein_P_T", g(ad["prot_tot"]["bnid104726_mM"] * 1e3), "uM chains", EXP,
            ECO, "B/r", "aerobic glucose minimal", "37 C", "BioNumbers BNID 104726 (Neidhardt 1996)",
            "audit:L06b", "ESTIMATED", "VERIFIED", "MISMATCHED",
            "STRAIN;MEDIUM;GROWTH_RATE;VOLUME_ASSUMED;QUANTITY", "X_PT_BNID",
            "calculated: 2.35e6 proteins/cell over an ASSUMED ~1 um3 volume, 40 min mass "
            "doubling (mu ~1.04 /h vs anchor 0.58 /h)",
            candidate="BNID_104726"),
        row(b, "total_protein_P_T", "..".join(g(x * 1e3) for x in
                                            ad["prot_tot"]["milo2013_mM_from_2to4e6_per_um3"]),
            "uM chains", "AMBIGUOUS", "bacteria, yeast, mammalian cells (general)", "NA",
            "NA", "NA", "Milo 2013 BioEssays, PMID:24114984", "audit:L06c", "ESTIMATED",
            "MISMATCHED_CONDITION", "MISMATCHED", "ORGANISM;QUANTITY", "X_PT_MILO",
            "generic cross-organism benchmark", candidate="Milo_2013"),
    ]
    for param, mach, gene, unit in MACHINES:
        rows.append(schmidt_row(b, param, EXP_COND, f"{gene}_uM_oligomer", unit, "MATCHED",
                                machine=mach, note="exponential 20-condition span in "
                                "schmidt_pools.json summary; not summed with other machines"))
    ne, ng = der["N_synthesis_weighted_exp_glucose"], der["N_genome_unweighted"]
    rows += [
        row(b, "copy_weighted_length_N", g(ne["weighted_mean_codons"]), "codons (copy-weighted "
            "mean)", EXP, ECO, "BW25113", "M9 glucose", "37 C",
            SCHMIDT + " x NC_000913.3 CDS", "derived:N_synthesis_weighted_exp_glucose",
            "DERIVED_FROM_MEASUREMENT", "VERIFIED", "MATCHED",
            note="copy weights = standing proteome; = synthesis weights only under the "
                 "balanced-growth, no-turnover closure; full distribution used in outputs; CDS "
                 "from MG1655 (BW25113 is a K-12 derivative)"),
        row(b, "copy_weighted_length_N", g(ng["weighted_mean_codons"]), "codons (genome mean)",
            "PHASE_INVARIANT", ECO, "K-12 MG1655", "NA (genome-encoded)",
            "NA (genome-encoded)", "UniProt UP000000625 (legacy length table)", "audit:L08b",
            "MEASURED", "VERIFIED", "MISMATCHED", "QUANTITY", "X_N_GENOME",
            "every gene counted once; not what is synthesized", candidate="genome_unweighted"),
        row(b, "codon_usage_weights", "vector:derived.json#codon_usage_translated_exp_glucose",
            "fraction of translated sense codons", EXP, ECO, "BW25113", "M9 glucose", "37 C",
            SCHMIDT + " x NC_000913.3 CDS", "derived:codon_usage_translated_exp_glucose",
            "DERIVED_FROM_MEASUREMENT", "VERIFIED", "MATCHED",
            note="copies x codon counts (balanced growth)"),
        row(b, "codon_usage_weights", "vector:global_codon_usage_ecoli.tsv", "fraction of "
            "genomic sense codons", "PHASE_INVARIANT", ECO, "K-12 MG1655",
            "NA (genome-encoded)", "NA (genome-encoded)", "envelope-paper global_codon_usage_ecoli.tsv",
            "audit:E05", "MEASURED", "MISMATCHED_CONDITION", "MISMATCHED", "QUANTITY",
            "X_CODON_GENOMIC", "legacy weights (WO-05)", candidate="genomic"),
        unmeasured(b, "substitution_freq_standing_f", "1/codon (MS-detected substitution)",
                   "no exponential-phase-annotated E. coli error rate in saved records; the "
                   "80 Landerer datasets carry no phase annotation (PRIDE metadata not saved)"),
    ]
    et = der["etel_translated_weights_exp_glucose"]
    for k, extra in [("conditional_mean_DataS2", ";SELECTION"), ("zero_inclusive_mean", ""),
                     ("pooled_ratio", ""), ("dataset_median", "")]:
        rows.append(row(b, "substitution_freq_standing_f", g(et[k]), "1/codon (MS-detected "
                        "substitution)", "MIXED", ECO, "MIXED (80 PRIDE datasets)",
                        "MIXED (80 PRIDE datasets)", "MIXED (80 PRIDE datasets)",
                        "Landerer 2024 MBE 41:msae048 Data_S4, PMID:38421032",
                        f"derived:etel_translated_weights_exp_glucose.{k}",
                        "DERIVED_FROM_MEASUREMENT", "MISMATCHED_CONDITION", "MISMATCHED",
                        "PHASE_MIXED;STRAIN;MEDIUM;TEMPERATURE" + extra, "X_F_ETEL",
                        "eTEL aggregate re-weighted by glucose translated codon usage; a lower "
                        "bound (MS; I/L invisible, rare events missed)", candidate=f"eTEL_{k}"))
    rows += invitro_rows(b, "X")
    return rows


def stationary(per, der, ad):
    b = STAT
    axes_s = "STRAIN;MEDIUM;TIME_IN_PHASE"
    norm = ";NORMALIZATION_FROM_EXPONENTIAL"
    rows = [
        unmeasured(b, "growth_rate_mu", "1/s", "no growth rate reported for the anchor culture"),
        schmidt_row(b, "growth_rate_mu", STAT_COND, "growth_per_h", "1/s", "MISMATCHED",
                    axes_s, "S_MU", "Schmidt tabulates -0.01 +- 0.003 /h: net growth ~0, "
                    "dilution is not a sink", scale=1 / 3600, candidate="Schmidt_1d"),
        unmeasured(b, "synthesis_rate_s_P", "uM chains/s", "no stationary synthesis rate in "
                   "any saved record; s_P = mu*P_T is a balanced-growth closure and fails here"),
        unmeasured(b, "total_protein_P_T", "uM chains", "no total protein for the anchor"),
        schmidt_row(b, "total_protein_P_T", STAT_COND, "total_protein_mM_wholecell",
                    "uM chains", "MISMATCHED", axes_s + norm, "S_MEDIUM_P_T;S_NORM",
                    "Schmidt scaled every non-glucose mass/cell by the GLUCOSE-exponential "
                    "volumetric protein concentration; 3 d: "
                    + g(per[STAT_COND3]["total_protein_mM_wholecell"] * 1e3),
                    scale=1e3, candidate="Schmidt_1d"),
    ]
    for param, mach, gene, unit in MACHINES:
        rows.append(unmeasured(b, param, unit, "no pool measured in the anchor culture", mach))
        rows.append(schmidt_row(b, param, STAT_COND, f"{gene}_uM_oligomer", unit, "MISMATCHED",
                                axes_s + norm, f"S_MEDIUM_{mach};S_NORM",
                                "3 d: " + g(per[STAT_COND3][f"{gene}_uM_oligomer"]) + "; absolute "
                                "uM inherits the exponential normalization, ratios to P_T do not",
                                mach, candidate="Schmidt_1d"))
    ng, ns = der["N_genome_unweighted"], der["N_abundance_weighted_stationary_1d"]
    rows += [
        unmeasured(b, "copy_weighted_length_N", "codons", "no proteome of the anchor culture "
                   "(Xac, LB, overnight) in saved records"),
        row(b, "copy_weighted_length_N", g(ng["weighted_mean_codons"]), "codons (genome mean)",
            "PHASE_INVARIANT", ECO, "K-12 MG1655", "NA (genome-encoded)",
            "NA (genome-encoded)", "UniProt UP000000625 (legacy length table)", "audit:L08b",
            "MEASURED", "VERIFIED", "MISMATCHED", "QUANTITY", "S_N",
            "every gene counted once", candidate="genome_unweighted"),
        row(b, "copy_weighted_length_N", g(ns["weighted_mean_codons"]), "codons (copy-weighted "
            "mean)", STAT, ECO, "BW25113", SCHMIDT_STAT_MEDIUM, "37 C",
            SCHMIDT + " x NC_000913.3 CDS", "derived:N_abundance_weighted_stationary_1d",
            "DERIVED_FROM_MEASUREMENT", "VERIFIED", "MISMATCHED",
            axes_s, "S_N", "standing proteome of a different culture; copies are not "
            "synthesis weights in non-growing cells", candidate="copy_weighted_1d"),
        unmeasured(b, "codon_usage_weights", "fraction of translated codons", "not needed for "
                   "the anchor input: Stikeleather f is pooled over all sites sampled"),
        row(b, "substitution_freq_standing_f", "0.00182", "1/codon (MS-detected substitution)", STAT,
            ECO, ANCHOR[STAT]["strain"], ANCHOR[STAT]["medium"], "37 C",
            "Stikeleather 2026 NAR 54(13) gkag674, PMID:42406629", "audit:E03", "MEASURED",
            "VERIFIED", "MATCHED", note="SE 5.92e-5; total detected substitutions / sites "
            "sampled in the STANDING proteome harvested at stationary phase, I/L merged; most "
            "chains were synthesized before growth stopped, so this is NOT a stationary "
            "per-synthesis error rate (that is UNMEASURED); genotype Xac [ara, "
            "delta(lac-proAB), gyrA, rpoB, argE(am)]; never used in the EXPONENTIAL bundle"),
    ]
    rows += invitro_rows(b, "S")
    return rows


def unmeasured_tail(bundle):
    return [unmeasured(bundle, p, u, why, mach) for p, (u, mach, why) in UNMEASURED_BOTH.items()]


def build():
    per, der, ad = load_sources()
    rows = exponential(per, der, ad) + unmeasured_tail(EXP)
    rows += stationary(per, der, ad) + unmeasured_tail(STAT)
    return rows


def write(rows, path=HERE / "bundles.tsv"):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    rows = build()
    write(rows)
    print(f"{len(rows)} rows ->", HERE / "bundles.tsv")
