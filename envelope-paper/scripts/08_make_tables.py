#!/usr/bin/env python3
"""
generate the paper's tables from data/computed/.

table numbering follows manuscript_v2_draft.md, in which the two excluded
analyses are no longer paper tables:

  Table 1   burden terms and their flux operationalization
  Table 2   bounds on the tolerable per-codon error rate (+ 2b, derived stats)
  Table 3   headroom against chaperone availability theta
  Table 4   supraadditivity against starting margin
  Table 5   permutation tests on the operational axes
  Table S1  per-codon (mu, nu) coordinates
  Table S2  per-amino-acid operational spread
  Table S3  validation of the nu axis
  Table S4  supraadditivity effect-size grid at the operating point
  Table S5  the two capacity knobs are not equivalent
  Table S6  what the axis tests could have detected
  Table S7  leave-one-codon-out jackknife on the mu axis

plus three files that are NOT paper tables and are named so they cannot be
mistaken for one: the two excluded analyses, and the retired anchoring grid that
Table 3 supersedes.

the functions below are named for their CONTENT, not their table number. the
number appears exactly once, in OUTPUTS, so renumbering cannot leave a function
called table7() emitting Table 3.

the TSVs carry full precision for deposition; TABLES.md carries the rounded,
typeset values. numbers in TABLES.md are formatted the same way as in the
manuscript so that tests/ can check the two against each other -- drift between
a generated table and the manuscript text is exactly the failure mode that sank
the previous draft.

every descriptive number in a TABLES.md caption is interpolated from the data,
not typed. captions went stale twice: they still quoted the window-bottom
saturation triple (97.9% / 47.5 uM / 0.052 uM) and "f = 1e-4" after the
manuscript had corrected both.
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


# ---------------------------------------------------------------- main tables
def burden_terms():
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


def bounds():
    """bounds on the maximum tolerable per-codon mistranslation rate."""
    b = load_json("bounds_summary.json")
    h = load_json("headroom_sensitivity_summary.json")
    lo, hi = b["observed_window_per_codon"]
    rows = [
        {"quantity": "Observed E. coli mistranslation rate",
         "value_per_codon": np.nan, "ci_low": lo, "ci_high": hi,
         "basis": "literature consensus"},
        {"quantity": "Combinatorial bound, deterministic reference point",
         "value_per_codon": b["arithmetic_deterministic"],
         "ci_low": np.nan, "ci_high": np.nan,
         "basis": "N=300, P_correct=0.7, p_misfold=0.3, S_syn=0.3"},
        {"quantity": "Combinatorial bound, paired MC median",
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
        {"statistic": "ODE / combinatorial at the paired median",
         "value": b["ode_over_arith_at_median"], "unit": "ratio"},
        {"statistic": "Median paired ratio r = f_arith / f_ODE",
         "value": b["paired_median_ratio_arith_over_ode"], "unit": "ratio"},
        {"statistic": "P(combinatorial is the tighter bound)",
         "value": b["paired_P_arith_tighter"], "unit": "probability"},
        {"statistic": "Draws closed by aggregation-death, not monomer runaway",
         "value": b["mechanism_frac_aggregation_death"], "unit": "fraction"},
        {"statistic": "Headroom, misfolded-monomer pool, at the usage-weighted "
                      "mean error rate (this paper's estimate)",
         "value": h["internally_consistent_headroom_P"], "unit": "fold"},
        {"statistic": "Headroom, aggregated pool, at the same rate",
         "value": h["internally_consistent_headroom_A"], "unit": "fold"},
        {"statistic": "Headroom at f = 10⁻⁴ (window bottom; quoted in "
                      "earlier drafts, not used here)",
         "value": b["headroom_P_at_window_bottom"], "unit": "fold"},
    ])
    return df, stats


def chaperone_availability():
    """headroom vs chaperone availability theta, over documented C_tot / K_d."""
    return pd.read_csv(COMP / "chaperone_availability.tsv", sep="\t")[
        ["C_tot_uM", "K_d_uM", "theta", "C_avail_uM", "c_free_uM",
         "folding_arm_saturation", "headroom_P", "margin_log10", "collapsed"]]


def supraadditivity_margin():
    """the supraadditivity test: interaction vs remaining margin."""
    df = pd.read_csv(COMP / "supraadditivity_margin_sweep.tsv", sep="\t")
    df["additive_expectation"] = df.D_error + df.D_capacity
    return df[["f_base", "margin_baseline", "error_factor", "capacity_factor",
               "D_error", "D_capacity", "additive_expectation", "D_both",
               "interaction", "interaction_pct_of_additive",
               "collapsed_both", "qualitative_supraadditive"]]


def axis_tests():
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


# ------------------------------------------------------- supplementary tables
def codon_coordinates():
    df = pd.read_csv(COMP / "codon_axes.tsv", sep="\t")
    return df[["codon", "aa", "degeneracy", "mu", "log_mu", "mu_z",
               "nu_tai", "nu_z"]].sort_values(["aa", "codon"]).reset_index(drop=True)


def delta_per_aa():
    m = pd.read_csv(COMP / "delta_per_aa.tsv", sep="\t")
    d = pd.read_csv(COMP / "delta_per_aa_median.tsv", sep="\t")
    out = m.merge(d[["aa", "delta_mu"]], on="aa", suffixes=("", "_mu_median"))
    return out.rename(columns={"delta_mu_mu_median": "delta_mu_median_stat"})


def tai_validation():
    return pd.DataFrame(load_json("tai_validation_report.json"))


def supraadditivity_grid():
    return pd.read_csv(COMP / "supraadditivity_effect_grid.tsv", sep="\t")[
        ["error_factor", "capacity_factor", "margin_baseline", "D_error",
         "D_capacity", "D_both", "interaction", "interaction_pct_of_additive",
         "collapsed_both", "qualitative_supraadditive"]]


def knob_comparison():
    return pd.read_csv(COMP / "supraadditivity_knob_comparison.tsv", sep="\t")[
        ["knob", "f_base", "margin_baseline", "D_capacity", "D_error",
         "interaction", "interaction_pct_of_additive", "collapsed_both"]]


def axis_power():
    """minimum detectable effect and power on each axis."""
    return pd.read_csv(COMP / "nu_power_sweep.tsv", sep="\t")


def mu_jackknife():
    """leave-one-codon-out jackknife on the mu axis."""
    return pd.read_csv(COMP / "mu_jackknife.tsv", sep="\t")


# ------------------------------------------- not paper tables: excluded work
def excluded_metal_sites():
    """the excluded metal-site result under both backgrounds."""
    df = pd.read_csv(COMP / "removed_metal_site_test.tsv", sep="\t")
    aa_name = {"D": "Asp", "C": "Cys", "E": "Glu", "H": "His"}
    df.insert(1, "amino_acid_name", df.amino_acid.map(aa_name))
    return df.sort_values("amino_acid_name").reset_index(drop=True)


def excluded_crossspecies():
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


def retired_anchoring_grid():
    """
    RETIRED as a main table. this crossed the evaluation point against six ad-hoc
    (C_tot, K_d) anchorings, but that axis is not independent of Table 3's theta:
    `usage_weighted_mu x c_free_at_Kd` (C_tot = 1) IS `C_tot = 50, theta = 0.98`
    (saturation 0.346, x4.6, margin 0.67), and the grid reached its low end only by
    using C_tot = 1-2 uM, outside the documented 30-80 uM range. It is one
    availability parameter written two ways. Kept as a repository cross-check
    only; Table 3 is the principled version.
    """
    df = pd.read_csv(COMP / "headroom_sensitivity.tsv", sep="\t")
    return df[["evaluation_point", "f_codon", "anchoring", "C_tot_uM", "K_d_uM",
               "c_free_uM", "folding_arm_saturation", "P_star", "P_dagger",
               "headroom_P", "headroom_A", "margin_log10", "collapsed"]]


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


# the one place a table number is written down
OUTPUTS = (
    ("Table1_burden_terms", burden_terms),
    ("Table2_bounds", None),                       # filled in main(), returns 2 frames
    ("Table2_bounds_statistics", None),
    ("Table3_chaperone_availability", chaperone_availability),
    ("Table4_supraadditivity", supraadditivity_margin),
    ("Table5_axis_tests", axis_tests),
    ("TableS1_codon_coordinates", codon_coordinates),
    ("TableS2_delta_per_aa", delta_per_aa),
    ("TableS3_tai_validation", tai_validation),
    ("TableS4_supraadditivity_grid", supraadditivity_grid),
    ("TableS5_capacity_knob_comparison", knob_comparison),
    ("TableS6_axis_power", axis_power),
    ("TableS7_mu_jackknife", mu_jackknife),
    ("Excluded_metal_site_backgrounds", excluded_metal_sites),
    ("Excluded_crossspecies", excluded_crossspecies),
    ("Retired_anchoring_grid", retired_anchoring_grid),
)

# files the earlier numbering wrote. they would otherwise sit in tables/ looking
# current, and a reader has no way to tell a stale Table 7 from a live one.
SUPERSEDED = ("Table3_axis_tests", "Table4_metal_site_backgrounds",
              "Table5_supraadditivity", "Table7_chaperone_availability",
              "TableS4_crossspecies", "TableS5_supraadditivity_grid",
              "TableS6_capacity_knob_comparison", "TableS7_axis_power",
              "TableS8_mu_jackknife", "TableS9_retired_anchoring_grid")


def main():
    TAB.mkdir(parents=True, exist_ok=True)

    t2, t2s = bounds()
    frames = {}
    for name, fn in OUTPUTS:
        frames[name] = t2 if name == "Table2_bounds" else (
            t2s if name == "Table2_bounds_statistics" else fn())

    for name, df in frames.items():
        df.to_csv(TAB / f"{name}.tsv", sep="\t", index=False)

    removed = []
    for name in SUPERSEDED:
        p = TAB / f"{name}.tsv"
        if p.exists():
            p.unlink()
            removed.append(name)

    t1 = frames["Table1_burden_terms"]
    t3 = frames["Table3_chaperone_availability"]
    t4 = frames["Table4_supraadditivity"]
    t5 = frames["Table5_axis_tests"]
    s1, s2, s3 = (frames["TableS1_codon_coordinates"],
                  frames["TableS2_delta_per_aa"],
                  frames["TableS3_tai_validation"])
    s4, s5 = (frames["TableS4_supraadditivity_grid"],
              frames["TableS5_capacity_knob_comparison"])
    s6, s7 = frames["TableS6_axis_power"], frames["TableS7_mu_jackknife"]
    xm, xc = (frames["Excluded_metal_site_backgrounds"],
              frames["Excluded_crossspecies"])

    # ---------------- TABLES.md ----------------
    # unicode minus, matching the manuscript, so tests can compare directly
    z3 = lambda v: f"{v:+.2f}".replace("-", "−")
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
        # a ratio labelled r is just a number; only the ODE/combinatorial
        # comparison reads naturally with a multiplication sign
        if r.statistic.startswith("ODE /"):
            return f"{r.value:.2f}×"
        return f"{r.value:.3f}"

    t2sv = t2s.copy()
    t2sv["value_fmt"] = [stat_fmt(r) for _, r in t2sv.iterrows()]

    # caption values, read from the data rather than typed. the previous version
    # of this file described the knob table with the window-bottom saturation
    # triple and the grid as "f = 1e-4", both retracted in the manuscript.
    sup = load_json("supraadditivity_summary.json")
    f_eval = sup["evaluation_point_f_codon"]
    sat_pct = 100 * sup["folding_arm_saturation_at_observed_rate"]
    c_free = sup["free_chaperone_uM_at_observed_rate"]
    p_mis = sup["misfolded_protein_uM_at_observed_rate"]
    ca = load_json("chaperone_availability_summary.json")
    hr_lo, hr_hi = ca["headroom_range_over_documented_grid"]
    pw = load_json("nu_power_summary.json")["minimum_detectable_effect"]
    jk = load_json("mu_jackknife_summary.json")
    mu_max = load_json("translation_burden.json")["mu_max"]

    parts = [
        "# Tables",
        "",
        "Generated by `scripts/08_make_tables.py` from `data/computed/`. "
        "Full-precision versions are the `.tsv` files in this directory. "
        "Numbering follows the manuscript; the three files at the end are **not** "
        "paper tables.",
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
        "## Table 3. Headroom against chaperone availability",
        "",
        "The model gives the whole chaperone pool to the damaged-protein pool "
        "because it does not represent nascent-chain folding. θ makes that "
        "assumption explicit: `C_avail = C_tot(1 − θ)`. `C_tot` (30 to 80 µM) and "
        "`K_d` (0.06 to 2 µM) are swept over their **documented** ranges; **θ is "
        "not measured** and is not estimated here. Across the full grid the "
        f"headroom spans ×{hr_lo:.1f} to ×{hr_hi:.0f}. Shown at K_d = 1 µM; the "
        "full grid is in `tables/Table3_chaperone_availability.tsv`. θ = 0 is the "
        "assumption earlier drafts made implicitly. This table supersedes the "
        "retired anchoring grid in `tables/TABLES.md`.",
        "",
        md_table(t3[(t3.K_d_uM == 1.0)],
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
        "## Table 4. Supraadditivity against starting margin",
        "",
        "2x2 factorial on the two-pool ODE: error rate x3 (raising `B_error`) "
        "against rescue throughput /3 (lowering `C_buffer`). The readout is the "
        "viability margin, log10(min(P_dagger/P*, A_max/A*)). Damage is margin "
        "lost; the interaction is observed joint damage minus the additive "
        "expectation, so positive means supraadditive. Where the combination has "
        "no stable state the margin loss is unbounded and no numeric interaction "
        "is defined -- those rows are marked, not assigned a value. The paper's "
        f"own operating point is f = {sci(f_eval)} per codon.",
        "",
        md_table(t4, ["f_base", "margin_baseline", "D_error", "D_capacity",
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
        "## Table 5. Permutation tests on the operational axes",
        "",
        "Δ is the mean pairwise distance among synonymous codons in standardized "
        "coordinates, averaged over the 18 multi-codon amino acids (59 codons). "
        "10,000 permutations per test, seed 42; one-sided empirical p in the "
        "observed direction with a +1 correction. The combined (μ, ν) rows are "
        "listed for completeness and are not interpreted: Δ in that space is a "
        "weighted mixture of a structured and an unstructured axis.",
        "",
        md_table(t5, ["axis", "mu_stat", "null", "observed", "null_mean",
                      "null_sd", "z", "p_one_sided", "direction"],
                 fmt={"observed": d4, "null_mean": d4, "null_sd": d4,
                      "z": z3, "p_one_sided": p4}).replace(
            "| axis | mu_stat | null | observed | null_mean | null_sd | z | p_one_sided | direction |",
            "| Axis | μ statistic | Null model | Observed Δ | Null mean | Null SD | z | p | Direction |"),
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
        "## Table S4. Supraadditivity effect-size grid at the operating point",
        "",
        f"All {len(s4)} perturbation combinations at the paper's own evaluation "
        f"point, f = {sci(f_eval)} per codon. Blank interaction means the "
        "combination collapses, where the loss is unbounded rather than large; "
        f"{int(s4.qualitative_supraadditive.sum())} of {len(s4)} combinations are "
        "survivable alone but lethal together.",
        "",
        md_table(s4, ["error_factor", "capacity_factor", "D_error", "D_capacity",
                      "D_both", "interaction_pct_of_additive", "collapsed_both",
                      "qualitative_supraadditive"],
                 fmt={"D_error": d4, "D_capacity": d4,
                      "D_both": lambda v: ("unbounded" if not np.isfinite(v)
                                           else f"{v:.4f}"),
                      "interaction_pct_of_additive":
                          lambda v: ("—" if not np.isfinite(v) else f"{v:+.2f}%"),
                      "collapsed_both": lambda v: "yes" if v else "",
                      "qualitative_supraadditive":
                          lambda v: "SL" if v else ""}).replace(
            "| error_factor | capacity_factor | D_error | D_capacity | D_both | "
            "interaction_pct_of_additive | collapsed_both | "
            "qualitative_supraadditive |",
            "| Error × | Capacity ÷ | D error | D capacity | Observed, both | "
            "% above additive | Joint collapse | Synthetic lethal |"),
        "",
        "## Table S5. The two capacity knobs are not equivalent",
        "",
        f"At the operating point the folding arm is {sat_pct:.1f}% saturated "
        f"({c_free:.1f} µM free chaperone against {p_mis:.3f} µM misfolded "
        f"protein, a {c_free / p_mis:.0f}-fold excess), so shrinking the "
        "chaperone pool `C_tot` barely changes the rescue rate while cutting "
        "throughput `k_obs_max` acts proportionally. This is a property of the "
        "model's parameterization, not of the framework — the model omits the "
        "nascent-chain folding load, so it hands the whole pool to the damaged "
        "pool. See Limitations.",
        "",
        md_table(s5, ["knob", "f_base", "margin_baseline", "D_capacity",
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
        "## Table S6. What the axis tests could have detected",
        "",
        "Within-amino-acid deviations are shrunk by a factor `s` (s = 1 the observed "
        "code, s = 0 all synonyms identical) and the permutation test re-run. "
        "`pct_below_null` makes the two axes comparable. The ν axis rejects only "
        f"once spread is tightened by "
        f"{pw['nu']['mde_delta_reduction_pct']:.0f}% "
        f"({pw['nu']['mde_pct_below_null']:.1f}% below null); μ's own clustering "
        f"sits {pw['mu']['observed_pct_below_null']:.1f}% below its null, so an "
        "effect of μ's size would have been seen on ν. See also "
        "`nu_power_curve.tsv` for the subset-model power curve.",
        "",
        md_table(s6[lambda d: d.s >= 0.5],
                 ["axis", "s", "delta", "null_mean", "z", "p", "reject",
                  "delta_reduction_pct", "pct_below_null"],
                 fmt={"s": lambda v: f"{v:.2f}", "delta": d4,
                      "null_mean": d4, "z": lambda v: f"{v:+.2f}",
                      "p": lambda v: f"{v:.4f}",
                      "reject": lambda v: "yes" if v else "",
                      "delta_reduction_pct": lambda v: f"{v:.0f}%",
                      "pct_below_null": lambda v: f"{v:.1f}%"}).replace(
            "| axis | s | delta | null_mean | z | p | reject | "
            "delta_reduction_pct | pct_below_null |",
            "| Axis | s | Δ | Null mean | z | p | Rejects | Δ reduction | "
            "% below null |"),
        "",
        "## Table S7. Leave-one-codon-out jackknife on the μ axis",
        "",
        "Each row drops one codon and re-runs the within-degeneracy permutation "
        f"test (each subset restandardized). The clustering survives all "
        f"{jk['jackknife']['n_total']} single deletions, including CCC — the "
        f"{sci(mu_max, 2)} maximum that sets the "   # 2 sig figs, as the prose has it
        f"{jk['mu_span_fold_with_CCC']:.0f}-fold span on its own, without which "
        f"the span falls to {jk['mu_span_fold_without_CCC']:.0f}-fold. Ordered "
        "from strongest to weakest z.",
        "",
        md_table(s7, ["dropped_codon", "aa", "mu", "z", "p", "n_codons"],
                 fmt={"mu": lambda v: sci(v), "z": lambda v: f"{v:+.2f}",
                      "p": lambda v: f"{v:.4f}"}).replace(
            "| dropped_codon | aa | mu | z | p | n_codons |",
            "| Codon dropped | AA | its μ | z | p | n codons |"),
        "",
        "---",
        "",
        "# Not paper tables",
        "",
        "The three sections below are kept for the record and are **not** results "
        "of this paper. The first two are the analyses excluded during "
        "verification; the third is a superseded version of Table 3.",
        "",
        "## Excluded analysis A: metal-binding-site codon enrichment",
        "",
        "**Excluded from the paper.** Reported here because the negative is "
        "informative. The earlier version used a genome-wide background; "
        "metalloproteins are enriched for abundant enzymes and abundant genes have "
        "stronger codon bias, so that comparison confounds site-level selection "
        "with gene-level expression. Against non-metal positions of the *same* "
        "genes, nothing survives.",
        "",
        md_table(xm, ["amino_acid_name", "enriched_codon", "alt_codon",
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
        "## Excluded analysis B: cross-species conservation of Δ",
        "",
        "**Excluded from the paper.** The earlier ρ used Δ in combined (μ, ν) "
        "space with **E. coli μ for all three species**, so the correlated vectors "
        "shared an identical coordinate. The ν-only columns isolate the genuinely "
        "species-varying part. The independent-tAI column is indicative only: it "
        "uses the prior script's own hardcoded tRNA counts, which are unverified.",
        "",
        md_table(xc, ["comparison", "rho_published_2D_shared_mu",
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
        "## Retired: the evaluation-point × anchoring grid",
        "",
        "**Superseded by Table 3.** This crossed the evaluation point against six "
        "ad-hoc (C_tot, K_d) anchorings, but that axis is not independent of "
        "Table 3's θ, and it reached its low end only by using C_tot = 1–2 µM, "
        "outside the documented 30–80 µM range. It is one availability parameter "
        "written two ways. Full grid in `Retired_anchoring_grid.tsv`.",
        "",
    ]
    (TAB / "TABLES.md").write_text("\n".join(parts) + "\n")

    print(f"wrote {len(frames) + 1} files to {TAB}")
    for name, df in frames.items():
        print(f"  {name + '.tsv':42s} {len(df):>4} rows")
    print(f"  {'TABLES.md':42s}")
    if removed:
        print("removed superseded files from the old numbering:")
        for name in removed:
            print(f"  {name}.tsv")


if __name__ == "__main__":
    main()
