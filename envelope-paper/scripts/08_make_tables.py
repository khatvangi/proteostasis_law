#!/usr/bin/env python3
"""
generate the paper's tables from data/computed/.

emits, into ../tables/:
  Table1_burden_terms.tsv          burden terms and their flux operationalization
  Table2_bounds.tsv                bounds on the tolerable per-codon error rate
  Table3_axis_tests.tsv            permutation tests on both axes, both mu stats
  Table4_metal_site_backgrounds.tsv  the removed metal-site result, both backgrounds
  TableS1_codon_coordinates.tsv    per-codon (mu, nu) and standardized coordinates
  TableS2_delta_per_aa.tsv         per-amino-acid operational spread
  TableS3_tai_validation.tsv       the nu-axis validation, both candidate vectors
  TableS4_crossspecies.tsv         the removed cross-species result
  TABLES.md                        all of the above, formatted for reading

the TSVs carry full precision for deposition; TABLES.md carries the rounded,
typeset values. numbers in TABLES.md are formatted the same way as in the
manuscript so that tests/ can check the two against each other -- drift between
a generated table and the manuscript text is exactly the failure mode that sank
the previous draft.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
COMP = ROOT / "data" / "computed"
TAB = ROOT / "tables"

SUP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def sci(x, sig=3):
    """format as '8.35 × 10⁻³', matching the manuscript's notation."""
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "—"
    if x == 0:
        return "0"
    exp = int(np.floor(np.log10(abs(x))))
    mant = x / 10 ** exp
    mant_s = f"{mant:.{sig - 1}f}"
    if float(mant_s) == 10.0:            # rounding pushed the mantissa over
        mant_s, exp = f"{1.0:.{sig-1}f}", exp + 1
    return f"{mant_s} × 10{str(exp).translate(SUP)}"


def load_json(name):
    return json.loads((COMP / name).read_text())


# --------------------------------------------------------------------------
def table1():
    """the burden decomposition. static by design -- it defines the vocabulary."""
    rows = [
        ("B_error", "decoding ambiguity / mistranslation burden",
         "rate of production of mistranslated chains "
         "(proteome-scale misincorporation frequency × synthesis flux)"),
        ("B_fold", "folding-failure burden",
         "rate at which nascent or mature chains fail to reach native state "
         "(folding-yield deficit × flux)"),
        ("B_agg", "aggregation / toxic load",
         "rate of irreversible conversion of misfolded monomer into aggregate"),
        ("B_qc", "rescue and clearance demand",
         "chaperone and protease occupancy imposed by the above"),
        ("C_buffer", "effective proteostasis network capacity",
         "finite throughput of the rescue-and-clearance network"),
        ("B_total", "summed burden; viability requires B_total ≤ C_buffer",
         "sum of the four burden fluxes"),
    ]
    return pd.DataFrame(rows, columns=["term", "meaning", "flux_operationalization"])


def table2():
    """bounds on the maximum tolerable per-codon mistranslation rate."""
    b = load_json("bounds_summary.json")
    h = load_json("headroom_sensitivity_summary.json")
    lo, hi = b["observed_window_per_codon"]
    rows = [
        {"quantity": "Observed E. coli mistranslation rate",
         "value_per_codon": np.nan, "ci_low": lo, "ci_high": hi,
         "basis": "literature consensus"},
        {"quantity": "Arithmetic bound, deterministic reference point",
         "value_per_codon": b["arithmetic_deterministic"],
         "ci_low": np.nan, "ci_high": np.nan,
         "basis": "N=300, P_correct=0.7, p_misfold=0.3, S_syn=0.3"},
        {"quantity": "Arithmetic bound, paired MC median",
         "value_per_codon": b["arithmetic_paired_median"],
         "ci_low": b["arithmetic_paired_ci95"][0],
         "ci_high": b["arithmetic_paired_ci95"][1],
         "basis": "5,000 paired draws"},
        {"quantity": "Two-pool ODE bound, deterministic",
         "value_per_codon": b["ode_deterministic"],
         "ci_low": np.nan, "ci_high": np.nan,
         "basis": "aggregation-death closure"},
        {"quantity": "Two-pool ODE bound, paired MC median",
         "value_per_codon": b["ode_paired_median"],
         "ci_low": b["ode_paired_ci95"][0],
         "ci_high": b["ode_paired_ci95"][1],
         "basis": "5,000 paired draws"},
    ]
    df = pd.DataFrame(rows)

    stats = pd.DataFrame([
        {"statistic": "ODE / arithmetic at the paired median",
         "value": b["ode_over_arith_at_median"], "unit": "ratio"},
        {"statistic": "Median paired ratio r = f_arith / f_ODE",
         "value": b["paired_median_ratio_arith_over_ode"], "unit": "ratio"},
        {"statistic": "P(arithmetic is the tighter bound)",
         "value": b["paired_P_arith_tighter"], "unit": "probability"},
        {"statistic": "Draws closed by aggregation-death, not monomer runaway",
         "value": b["mechanism_frac_aggregation_death"], "unit": "fraction"},
        {"statistic": "Headroom, misfolded-monomer pool, at the usage-weighted "
                      "mean error rate (this paper's estimate)",
         "value": h["internally_consistent_headroom_P"], "unit": "fold"},
        {"statistic": "Headroom, aggregated pool, at the same rate",
         "value": h["internally_consistent_headroom_A"], "unit": "fold"},
        {"statistic": "Headroom at f = 10\u207b\u2074 (window bottom; quoted in "
                      "earlier drafts, not used here)",
         "value": b["headroom_P_at_window_bottom"], "unit": "fold"},
    ])
    return df, stats


def table3():
    """permutation tests on both axes, under both mu summary statistics."""
    frames = []
    for stat, fname in (("mean", "axis_tests.tsv"),
                        ("median", "axis_tests_median.tsv")):
        t = pd.read_csv(COMP / fname, sep="\t")
        t["mu_stat"] = stat
        frames.append(t)
    df = pd.concat(frames, ignore_index=True)
    df = df[["axis", "mu_stat", "null", "observed", "null_mean", "null_sd",
             "z", "p_one_sided", "direction"]]

    # the nu axis does not depend on how mu is summarized, so the median rows for
    # nu are byte-identical duplicates of the mean rows. keeping both would read
    # as two independent tests. drop them and mark nu as applying to both.
    dup = (df.axis == "nu") & (df.mu_stat == "median")
    df = df[~dup].copy()
    # "both", not "n/a": the nu result holds under either mu statistic, and
    # pandas parses "n/a" as NaN on read-back
    df.loc[df.axis == "nu", "mu_stat"] = "both"

    # "clustered"/"spread" describes the sign, which is meaningless when the test
    # is non-significant. say so instead.
    df.loc[df.p_one_sided > 0.05, "direction"] = "no signal"

    order = {"mu": 0, "nu": 1, "2D": 2}
    df = (df.assign(_o=df.axis.map(order),
                    _s=(df.mu_stat == "median").astype(int),
                    _n=(df.null == "full").astype(int))
            .sort_values(["_o", "_s", "_n"])
            .drop(columns=["_o", "_s", "_n"]).reset_index(drop=True))
    return df


def table4():
    """the removed metal-site result under both backgrounds."""
    df = pd.read_csv(COMP / "removed_metal_site_test.tsv", sep="\t")
    aa_name = {"D": "Asp", "C": "Cys", "E": "Glu", "H": "His"}
    df.insert(1, "amino_acid_name", df.amino_acid.map(aa_name))
    return df.sort_values("amino_acid_name").reset_index(drop=True)


def table5():
    """the supraadditivity test: interaction vs remaining margin."""
    df = pd.read_csv(COMP / "supraadditivity_margin_sweep.tsv", sep="\t")
    df["additive_expectation"] = df.D_error + df.D_capacity
    return df[["f_base", "margin_baseline", "error_factor", "capacity_factor",
               "D_error", "D_capacity", "additive_expectation", "D_both",
               "interaction", "interaction_pct_of_additive",
               "collapsed_both", "qualitative_supraadditive"]]


def table_s5():
    return pd.read_csv(COMP / "supraadditivity_effect_grid.tsv", sep="\t")[
        ["error_factor", "capacity_factor", "margin_baseline", "D_error",
         "D_capacity", "D_both", "interaction", "interaction_pct_of_additive",
         "collapsed_both", "qualitative_supraadditive"]]


def table_s6():
    return pd.read_csv(COMP / "supraadditivity_knob_comparison.tsv", sep="\t")[
        ["knob", "f_base", "margin_baseline", "D_capacity", "D_error",
         "interaction", "interaction_pct_of_additive", "collapsed_both"]]


def table6():
    """headroom across evaluation point and chaperone anchoring."""
    df = pd.read_csv(COMP / "headroom_sensitivity.tsv", sep="\t")
    return df[["evaluation_point", "f_codon", "anchoring", "C_tot_uM", "K_d_uM",
               "c_free_uM", "folding_arm_saturation", "P_star", "P_dagger",
               "headroom_P", "headroom_A", "margin_log10", "collapsed"]]


def table7():
    """headroom vs chaperone availability theta, over documented C_tot / K_d."""
    return pd.read_csv(COMP / "chaperone_availability.tsv", sep="\t")[
        ["C_tot_uM", "K_d_uM", "theta", "C_avail_uM", "c_free_uM",
         "folding_arm_saturation", "headroom_P", "margin_log10", "collapsed"]]


def table_s1():
    df = pd.read_csv(COMP / "codon_axes.tsv", sep="\t")
    return df[["codon", "aa", "degeneracy", "mu", "log_mu", "mu_z",
               "nu_tai", "nu_z"]].sort_values(["aa", "codon"]).reset_index(drop=True)


def table_s2():
    m = pd.read_csv(COMP / "delta_per_aa.tsv", sep="\t")
    d = pd.read_csv(COMP / "delta_per_aa_median.tsv", sep="\t")
    out = m.merge(d[["aa", "delta_mu"]], on="aa",
                  suffixes=("", "_mu_median"))
    return out.rename(columns={"delta_mu_mu_median": "delta_mu_median_stat"})


def table_s3():
    rep = load_json("tai_validation_report.json")
    return pd.DataFrame(rep)


def table_s4():
    s = load_json("removed_crossspecies_test.json")
    rows = []
    pub = s["published_rho_2D"]
    for comp in pub:
        a, b = comp.split("_vs_")
        rows.append({
            "comparison": f"{a} vs {b}",
            "rho_published_2D_shared_mu": pub[comp],
            "rho_nu_only_vectors_as_used":
                s["rho_nu_only_vectors_as_used"][comp]["rho"],
            "p_nu_only_vectors_as_used":
                s["rho_nu_only_vectors_as_used"][comp]["p"],
            "rho_nu_only_independent_tai":
                s["rho_nu_only_independent_canonical_tai"][comp]["rho"],
            "p_nu_only_independent_tai":
                s["rho_nu_only_independent_canonical_tai"][comp]["p"],
        })
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
def md_table(df, cols=None, fmt=None):
    """render a dataframe as a markdown table."""
    cols = cols or list(df.columns)
    fmt = fmt or {}
    head = "| " + " | ".join(cols) + " |"
    rule = "|" + "|".join("---" for _ in cols) + "|"
    lines = [head, rule]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            f = fmt.get(c)
            if f:
                cells.append(f(v))
            elif isinstance(v, float):
                cells.append("—" if not np.isfinite(v) else f"{v:g}")
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main():
    TAB.mkdir(parents=True, exist_ok=True)

    t1 = table1()
    t2, t2s = table2()
    t5, s5, s6 = table5(), table_s5(), table_s6()
    t6, t7 = table6(), table7()
    t3 = table3()
    t4 = table4()
    s1, s2, s3, s4 = table_s1(), table_s2(), table_s3(), table_s4()

    written = {}
    for name, df in (("Table1_burden_terms", t1),
                     ("Table2_bounds", t2),
                     ("Table2_bounds_statistics", t2s),
                     ("Table3_axis_tests", t3),
                     ("Table4_metal_site_backgrounds", t4),
                     ("TableS1_codon_coordinates", s1),
                     ("TableS2_delta_per_aa", s2),
                     ("TableS3_tai_validation", s3),
                     ("TableS4_crossspecies", s4),
                     ("Table5_supraadditivity", t5),
                     ("TableS5_supraadditivity_grid", s5),
                     ("TableS6_capacity_knob_comparison", s6),
                     ("Table6_headroom_sensitivity", t6),
                     ("Table7_chaperone_availability", t7)):
        p = TAB / f"{name}.tsv"
        df.to_csv(p, sep="\t", index=False)
        written[name] = len(df)

    # ---------------- TABLES.md ----------------
    # unicode minus, matching the manuscript, so tests can compare directly
    z3 = lambda v: f"{v:+.2f}".replace("-", "\u2212")
    p4 = lambda v: f"{v:.4f}"
    d4 = lambda v: f"{v:.4f}"

    def rng(r):
        if not np.isfinite(r.ci_low):
            return "—"
        return f"{sci(r.ci_low)} – {sci(r.ci_high)}"

    t2v = t2.copy()
    t2v["value_fmt"] = t2v.value_per_codon.map(
        lambda v: "—" if not np.isfinite(v) else sci(v))
    t2v["ci_fmt"] = [rng(r) for _, r in t2v.iterrows()]

    def stat_fmt(r):
        if r.unit == "fold":
            return f"×{r.value:,.0f}"
        if r.unit in ("probability", "fraction"):
            return f"{r.value:.4f}"
        # a ratio labelled r is just a number; only the ODE/arithmetic comparison
        # reads naturally with a multiplication sign
        if r.statistic.startswith("ODE /"):
            return f"{r.value:.2f}×"
        return f"{r.value:.3f}"

    t2sv = t2s.copy()
    t2sv["value_fmt"] = [stat_fmt(r) for _, r in t2sv.iterrows()]

    parts = [
        "# Tables",
        "",
        "Generated by `scripts/08_make_tables.py` from `data/computed/`. "
        "Full-precision versions are the `.tsv` files in this directory.",
        "",
        "---",
        "",
        "## Table 1. Burden terms and their flux operationalization",
        "",
        "Terms are fluxes through stages, not disjoint pools, so a single "
        "misfolded protein is not counted three times.",
        "",
        md_table(t1, ["term", "meaning", "flux_operationalization"]),
        "",
        "---",
        "",
        "## Table 2. Bounds on the maximum tolerable per-codon mistranslation rate",
        "",
        "All two-bound comparisons use the **paired** Monte Carlo, in which the "
        "parameters shared by both bounds are drawn once per sample "
        "(5,000 draws, seed 17). Statistics derived from separately parameterized "
        "marginal runs are not valid and are not reported.",
        "",
        md_table(t2v, ["quantity", "value_fmt", "ci_fmt", "basis"]).replace(
            "| quantity | value_fmt | ci_fmt | basis |",
            "| Quantity | Value (/codon) | 95% CI | Basis |"),
        "",
        "### Table 2b. Derived comparison statistics",
        "",
        md_table(t2sv, ["statistic", "value_fmt"]).replace(
            "| statistic | value_fmt |", "| Statistic | Value |"),
        "",
        "---",
        "",
        "## Table 3. Permutation tests on the operational axes",
        "",
        "Δ is the mean pairwise distance among synonymous codons in standardized "
        "coordinates, averaged over the 18 multi-codon amino acids (59 codons). "
        "10,000 permutations per test, seed 42; one-sided empirical p in the "
        "observed direction with a +1 correction.",
        "",
        md_table(t3, ["axis", "mu_stat", "null", "observed", "null_mean",
                      "null_sd", "z", "p_one_sided", "direction"],
                 fmt={"observed": d4, "null_mean": d4, "null_sd": d4,
                      "z": z3, "p_one_sided": p4}).replace(
            "| axis | mu_stat | null | observed | null_mean | null_sd | z | p_one_sided | direction |",
            "| Axis | μ statistic | Null model | Observed Δ | Null mean | Null SD | z | p | Direction |"),
        "",
        "---",
        "",
        "## Table 4. A removed result: metal-binding-site codon enrichment",
        "",
        "**This result is excluded from the paper.** It is reported here because "
        "the negative is informative. The published version used a genome-wide "
        "background; metalloproteins are enriched for abundant enzymes and "
        "abundant genes have stronger codon bias, so that comparison confounds "
        "site-level selection with gene-level expression. Against non-metal "
        "positions of the *same* genes, nothing survives.",
        "",
        md_table(t4, ["amino_acid_name", "enriched_codon", "alt_codon",
                      "n_sites_enriched", "n_sites_alt",
                      "genomewide_or", "genomewide_p",
                      "within_gene_or", "within_gene_p",
                      "published_or", "published_p"],
                 fmt={"genomewide_or": lambda v: f"{v:.3f}",
                      "genomewide_p": lambda v: f"{v:.4f}",
                      "within_gene_or": lambda v: f"{v:.3f}",
                      "within_gene_p": lambda v: f"{v:.4f}",
                      "published_or": lambda v: f"{v:.3f}",
                      "published_p": lambda v: f"{v:.4f}"}).replace(
            "| amino_acid_name | enriched_codon | alt_codon | n_sites_enriched | "
            "n_sites_alt | genomewide_or | genomewide_p | within_gene_or | "
            "within_gene_p | published_or | published_p |",
            "| Amino acid | Enriched | Alternative | n sites (enr) | n sites (alt) | "
            "OR (genome-wide) | p | **OR (within-gene)** | **p** | OR (published) | p |"),
        "",
        "---",
        "",
        "## Table 5. The distinguishing prediction, tested in the model",
        "",
        "2x2 factorial on the two-pool ODE: error rate x3 (raising `B_error`) "
        "against rescue throughput /3 (lowering `C_buffer`). The readout is the "
        "viability margin, log10(min(P_dagger/P*, A_max/A*)). Damage is margin "
        "lost; the interaction is observed joint damage minus the additive "
        "expectation, so positive means supraadditive. Where the combination has "
        "no stable state the margin loss is unbounded and no numeric interaction "
        "is defined -- those rows are marked, not assigned a value.",
        "",
        md_table(t5, ["f_base", "margin_baseline", "D_error", "D_capacity",
                      "additive_expectation", "D_both", "interaction",
                      "interaction_pct_of_additive", "collapsed_both"],
                 fmt={"f_base": lambda v: sci(v),
                      "margin_baseline": lambda v: f"{v:.2f}",
                      "D_error": d4, "D_capacity": d4,
                      "additive_expectation": d4,
                      "D_both": lambda v: ("unbounded" if not np.isfinite(v)
                                           else f"{v:.4f}"),
                      "interaction": lambda v: ("—" if not np.isfinite(v)
                                                else f"{v:+.4f}"),
                      "interaction_pct_of_additive":
                          lambda v: ("—" if not np.isfinite(v) else f"{v:+.1f}%"),
                      "collapsed_both": lambda v: "yes" if v else ""}).replace(
            "| f_base | margin_baseline | D_error | D_capacity | "
            "additive_expectation | D_both | interaction | "
            "interaction_pct_of_additive | collapsed_both |",
            "| Baseline f (/codon) | Starting margin (log₁₀) | D error ×3 | "
            "D capacity ÷3 | Additive | **Observed, both** | Interaction | "
            "% above additive | Joint collapse |"),
        "",
        "---",
        "",
        "## Table 6. How far inside the envelope, across two explicit choices",
        "",
        "Headroom to collapse depends on where in the observed window the model is "
        "evaluated and on how the chaperone arm is anchored. The internally "
        "consistent cell is **usage_weighted_mu x as_published** (x25); the x158 "
        "quoted in earlier drafts is **window_bottom x as_published**, the most "
        "favourable corner. `folding_arm_saturation` shows why the anchoring "
        "matters: at the published values the rescue arm is 97.9% saturated.",
        "",
        md_table(t6, ["evaluation_point", "f_codon", "anchoring",
                      "C_tot_uM", "K_d_uM", "folding_arm_saturation",
                      "headroom_P", "headroom_A", "margin_log10"],
                 fmt={"f_codon": lambda v: sci(v),
                      "C_tot_uM": lambda v: f"{v:g}", "K_d_uM": lambda v: f"{v:g}",
                      "folding_arm_saturation":
                          lambda v: "—" if not np.isfinite(v) else f"{v:.3f}",
                      "headroom_P":
                          lambda v: "collapsed" if not np.isfinite(v) else f"×{v:.1f}",
                      "headroom_A":
                          lambda v: "—" if not np.isfinite(v) else f"×{v:.0f}",
                      "margin_log10":
                          lambda v: "—" if not np.isfinite(v) else f"{v:.2f}"}).replace(
            "| evaluation_point | f_codon | anchoring | C_tot_uM | K_d_uM | "
            "folding_arm_saturation | headroom_P | headroom_A | margin_log10 |",
            "| Evaluated at | f (/codon) | Chaperone anchoring | C_tot (µM) | "
            "K_d (µM) | Folding arm saturation | **Headroom, P pool** | "
            "Headroom, A pool | Margin (log₁₀) |"),
        "",
        "---",
        "",
        "## Table 7. Headroom against chaperone availability",
        "",
        "The model gives the whole chaperone pool to the damaged-protein pool "
        "because it does not represent nascent-chain folding. `theta` makes that "
        "assumption explicit: `C_avail = C_tot(1 - theta)`. `C_tot` (30-80 µM) and "
        "`K_d` (0.06-2 µM) are swept over their **documented** ranges; **theta is "
        "not measured** and is not estimated here. Shown at K_d = 1 µM; the full "
        "grid is in the TSV. theta = 0 is the assumption earlier drafts made "
        "implicitly.",
        "",
        md_table(t7[(t7.K_d_uM == 1.0)],
                 ["C_tot_uM", "theta", "C_avail_uM", "c_free_uM",
                  "folding_arm_saturation", "headroom_P", "margin_log10"],
                 fmt={"C_tot_uM": lambda v: f"{v:g}",
                      "theta": lambda v: f"{v:g}",
                      "C_avail_uM": lambda v: f"{v:.2f}",
                      "c_free_uM":
                          lambda v: "—" if not np.isfinite(v) else f"{v:.2f}",
                      "folding_arm_saturation":
                          lambda v: "—" if not np.isfinite(v) else f"{v:.3f}",
                      "headroom_P":
                          lambda v: "collapsed" if not np.isfinite(v) else f"×{v:.1f}",
                      "margin_log10":
                          lambda v: "—" if not np.isfinite(v) else f"{v:.2f}"}).replace(
            "| C_tot_uM | theta | C_avail_uM | c_free_uM | "
            "folding_arm_saturation | headroom_P | margin_log10 |",
            "| C_tot (µM) | θ | C_avail (µM) | Free chaperone (µM) | "
            "Folding arm saturation | **Headroom** | Margin (log₁₀) |"),
        "",
        "---",
        "",
        "# Supplementary tables",
        "",
        "## Table S1. Per-codon operational coordinates",
        "",
        f"All {len(s1)} analysed sense codons. μ from Landerer et al. Data S2 "
        "(per-codon mean); ν is the validated E. coli tAI. Met and Trp are "
        "excluded as single-codon amino acids.",
        "",
        md_table(s1, ["codon", "aa", "degeneracy", "mu", "nu_tai", "mu_z", "nu_z"],
                 fmt={"mu": lambda v: sci(v), "nu_tai": lambda v: f"{v:.4f}",
                      "mu_z": z3, "nu_z": z3}).replace(
            "| codon | aa | degeneracy | mu | nu_tai | mu_z | nu_z |",
            "| Codon | AA | Degeneracy | μ (/codon) | ν (tAI) | μ̃ | ν̃ |"),
        "",
        "## Table S2. Operational spread per amino acid",
        "",
        md_table(s2, ["aa", "degeneracy", "delta_mu", "delta_nu", "delta_2D",
                      "delta_mu_median_stat"],
                 fmt={"delta_mu": d4, "delta_nu": d4, "delta_2D": d4,
                      "delta_mu_median_stat": d4}).replace(
            "| aa | degeneracy | delta_mu | delta_nu | delta_2D | delta_mu_median_stat |",
            "| AA | Degeneracy | Δ (μ) | Δ (ν) | Δ (μ,ν) | Δ (μ, median stat) |"),
        "",
        "## Table S3. Validation of the translational-supply (ν) axis",
        "",
        "Translational selection predicts that highly expressed genes shift "
        "synonymous usage toward well-served codons. The per-amino-acid sign test "
        "is the discriminating criterion, because Δ is computed strictly within "
        "amino acids; the pooled correlation is ~+0.31 for both vectors and is "
        "therefore not diagnostic on its own.",
        "",
        md_table(s3, ["vector", "rho_tai_vs_delta_usage", "p_delta",
                      "rho_tai_vs_usage_ribosomal",
                      "n_aa_where_high_tai_synonym_gains",
                      "p_sign_test", "verdict"]).replace(
            "| vector | rho_tai_vs_delta_usage | p_delta | "
            "rho_tai_vs_usage_ribosomal | n_aa_where_high_tai_synonym_gains | "
            "p_sign_test | verdict |",
            "| Candidate ν vector | ρ(tAI, usage shift) | p | ρ(tAI, usage in "
            "ribosomal genes) | AAs where higher-tAI synonym gains | Sign test p | Verdict |"),
        "",
        "## Table S4. A removed result: cross-species conservation of Δ",
        "",
        "**This result is excluded from the paper.** The published ρ used Δ in "
        "combined (μ, ν) space with **E. coli μ for all three species**, so the "
        "correlated vectors shared an identical coordinate. The ν-only columns "
        "isolate the genuinely species-varying part. The independent-tAI column is "
        "indicative only: it uses the prior script's own hardcoded tRNA counts, "
        "which are unverified.",
        "",
        md_table(s4, ["comparison", "rho_published_2D_shared_mu",
                      "rho_nu_only_vectors_as_used", "p_nu_only_vectors_as_used",
                      "rho_nu_only_independent_tai", "p_nu_only_independent_tai"],
                 fmt={"p_nu_only_vectors_as_used": lambda v: f"{v:.3g}",
                      "p_nu_only_independent_tai": lambda v: f"{v:.3g}"}).replace(
            "| comparison | rho_published_2D_shared_mu | "
            "rho_nu_only_vectors_as_used | p_nu_only_vectors_as_used | "
            "rho_nu_only_independent_tai | p_nu_only_independent_tai |",
            "| Comparison | ρ published (2D, shared μ) | ρ (ν only, vectors as "
            "used) | p | ρ (ν only, independent tAI) | p |"),
        "",
        "## Table S5. Supraadditivity effect-size grid at the observed rate",
        "",
        f"All {len(s5)} perturbation combinations at f = 1e-4. Blank interaction "
        "means the combination collapses, where the loss is unbounded.",
        "",
        md_table(s5, ["error_factor", "capacity_factor", "D_error", "D_capacity",
                      "D_both", "interaction_pct_of_additive", "collapsed_both"],
                 fmt={"D_error": d4, "D_capacity": d4,
                      "D_both": lambda v: ("unbounded" if not np.isfinite(v)
                                           else f"{v:.4f}"),
                      "interaction_pct_of_additive":
                          lambda v: ("—" if not np.isfinite(v) else f"{v:+.2f}%"),
                      "collapsed_both": lambda v: "yes" if v else ""}).replace(
            "| error_factor | capacity_factor | D_error | D_capacity | D_both | "
            "interaction_pct_of_additive | collapsed_both |",
            "| Error × | Capacity ÷ | D error | D capacity | Observed, both | "
            "% above additive | Joint collapse |"),
        "",
        "## Table S6. The two capacity knobs are not equivalent",
        "",
        "At the observed operating point the folding arm is 97.9% saturated "
        "(47.5 µM free chaperone against 0.052 µM misfolded protein), so shrinking "
        "the chaperone pool `C_tot` barely changes the rescue rate while cutting "
        "throughput `k_obs_max` acts proportionally. This is a property of the "
        "model's parameterization, not of the framework — see Limitations.",
        "",
        md_table(s6, ["knob", "f_base", "margin_baseline", "D_capacity",
                      "D_error", "interaction_pct_of_additive", "collapsed_both"],
                 fmt={"f_base": lambda v: sci(v),
                      "margin_baseline": lambda v: f"{v:.2f}",
                      "D_capacity": d4, "D_error": d4,
                      "interaction_pct_of_additive":
                          lambda v: ("—" if not np.isfinite(v) else f"{v:+.1f}%"),
                      "collapsed_both": lambda v: "yes" if v else ""}).replace(
            "| knob | f_base | margin_baseline | D_capacity | D_error | "
            "interaction_pct_of_additive | collapsed_both |",
            "| Capacity knob | Baseline f | Starting margin | D capacity ÷3 | "
            "D error ×3 | % above additive | Joint collapse |"),
        "",
    ]
    (TAB / "TABLES.md").write_text("\n".join(parts) + "\n")

    print(f"wrote {len(written) + 1} files to {TAB}")
    for name, n in written.items():
        print(f"  {name + '.tsv':42s} {n:>4} rows")
    print(f"  {'TABLES.md':42s}")


if __name__ == "__main__":
    main()
