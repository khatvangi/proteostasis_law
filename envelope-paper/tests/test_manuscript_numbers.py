#!/usr/bin/env python3
"""
every quantitative claim in MANUSCRIPT.md is asserted here against a computed
artifact in data/computed/ or a raw input in data/raw/.

this exists because of how the previous version of this project failed. its
numbers traced to a markdown summary, which traced to another markdown summary.
its verification pass compared the draft against precomputed TSVs but never
compared those TSVs against generating code, so it certified a corrupt supply
axis, a superseded Monte Carlo run, and a statistic the upstream project had
explicitly retracted -- 17/20 "pass".

so the rule enforced here is: a number may appear in the manuscript only if it
can be recomputed from raw inputs by a script in scripts/. run:

    python -m unittest discover -s tests -v

from the paper root. scripts/01 through scripts/15 must have been run first.

a second rule was added after the v2 restructuring: consistency between two
places is itself asserted. per-number provenance is satisfied by a document whose
sections are each sourced and mutually contradictory -- which is exactly what
happened when Result 5's script hardcoded the evaluation point Result 3 rejects.
"""
import json
import re
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
COMP = ROOT / "data" / "computed"
MSDIR = ROOT / "manuscript"
MS = MSDIR / "MANUSCRIPT.md"
PAPER = MSDIR / "PAPER.md"
FIGS = ROOT / "figures"
SCRIPTS = ROOT / "scripts"
TAB = ROOT / "tables"

# section headings, in one place, so renaming a section is a one-line change
# rather than a hunt through scoped assertions
H_ABSTRACT = "## Abstract"
H_BURDEN = "### Burden is measurable and buffering capacity is finite"
H_BOUNDS = "### Two independent bounds agree on where the envelope's edge lies"
H_HEADROOM = "### *E. coli* operates about one order of magnitude inside the envelope"
H_SUPRA = ("### Burden and capacity perturbations interact supraadditively "
           "where *E. coli* sits")
H_AXES = ("### The measured error axis is organized at the amino-acid level; "
          "the supply axis is not")
H_EXCLUDED = "### Two analyses excluded during verification"
H_LIMITS = "### Limitations"
H_PREDICTIONS = "### Predictions"


def manuscript():
    return MS.read_text()


def section(heading, text=None):
    """
    the text of one section, from its heading to the next heading of the same or
    higher level.

    scoping matters. an assertion that a value appears SOMEWHERE in the
    manuscript passes even when the abstract and the results disagree -- which is
    precisely how the previous draft carried a retracted statistic in one section
    and the corrected one in another. headline claims are asserted against the
    section that makes them.
    """
    t = manuscript() if text is None else text
    i = t.index(heading)
    level = len(heading) - len(heading.lstrip("#"))
    rest = t[i + len(heading):]
    # a marker like "\n## " matches only level-2 headings, since a level-3
    # heading reads "\n### " and fails on the space
    end = len(rest)
    for lv in range(1, level + 1):
        j = rest.find("\n" + "#" * lv + " ")
        if j != -1:
            end = min(end, j)
    return heading + rest[:end]


def flat(text=None):
    """
    manuscript text with runs of whitespace collapsed.

    prose is hard-wrapped, so a phrase assertion can straddle a newline and an
    indent. normalize before matching phrases, so the tests constrain the CLAIM
    rather than the line breaks -- otherwise reflowing a paragraph breaks a test
    and the temptation is to reword the paper to fit the test.
    """
    return " ".join((manuscript() if text is None else text).split())


def load_json(p):
    return json.loads(Path(p).read_text())


class ManuscriptContainsClaim(unittest.TestCase):
    """the manuscript must literally state the values the analysis produced."""

    def setUp(self):
        self.text = manuscript()
        self.flat = flat()

    def assertStated(self, needle, why):
        self.assertIn(needle, self.text,
                      f"manuscript does not state {needle!r} ({why})")

    # ---------- burden magnitude ----------
    def test_mu_span_and_range(self):
        b = load_json(COMP / "translation_burden.json")
        self.assertAlmostEqual(b["mu_fold_span"], 613, delta=1)
        self.assertStated("613-fold", "per-codon mu span")
        self.assertStated("3.3 × 10⁻⁵", "mu minimum")
        self.assertStated("2.0 × 10⁻²", "mu maximum")

    def test_usage_weighted_mu(self):
        b = load_json(COMP / "translation_burden.json")
        self.assertAlmostEqual(b["usage_weighted_mean_mu_per_codon"], 6.334e-4,
                               delta=5e-6)
        self.assertStated("6.3 × 10⁻⁴", "usage-weighted mean mu")

    def test_fraction_of_proteins_carrying_an_error(self):
        b = load_json(COMP / "translation_burden.json")
        lo = b["frac_proteins_with_ge1_error_at_median_length"]
        hi = b["frac_proteins_with_ge1_error_at_mean_length"]
        self.assertTrue(0.15 < lo < 0.165, lo)
        self.assertTrue(0.17 < hi < 0.185, hi)
        self.assertIn("16 to 18% of synthesized proteins", self.flat)
        # and the claim the earlier draft made must NOT reappear
        self.assertNotIn("one-fifth to one-quarter", self.text,
                         "the unsupported 20-25% claim is back")

    def test_proteome_description(self):
        b = load_json(COMP / "translation_burden.json")
        self.assertEqual(b["proteome_n"], 4403)
        self.assertEqual(b["median_length_aa"], 271)
        self.assertStated("n = 4,403", "proteome size")
        self.assertStated("median 271 aa", "median protein length")

    # ---------- bounds ----------
    def test_deterministic_bounds(self):
        s = load_json(COMP / "bounds_summary.json")
        self.assertAlmostEqual(s["arithmetic_deterministic"], 1.188e-3, delta=1e-5)
        self.assertAlmostEqual(s["ode_deterministic"], 1.000e-2, delta=1e-4)
        self.assertStated("1.19 × 10⁻³", "combinatorial deterministic bound")
        self.assertStated("1.00 × 10⁻²", "ODE deterministic bound")

    def test_paired_medians_and_intervals(self):
        s = load_json(COMP / "bounds_summary.json")
        self.assertAlmostEqual(s["arithmetic_paired_median"], 8.348e-3, delta=1e-5)
        self.assertAlmostEqual(s["ode_paired_median"], 2.589e-2, delta=1e-4)
        self.assertStated("8.35 × 10⁻³", "combinatorial paired median")
        self.assertStated("2.59 × 10⁻²", "ODE paired median")
        self.assertStated("1.38 × 10⁻³", "combinatorial CI low")
        self.assertStated("8.79 × 10⁻²", "combinatorial CI high")
        self.assertStated("1.25 × 10⁻³", "ODE CI low")
        self.assertStated("5.87 × 10⁻¹", "ODE CI high")

    def test_paired_ratio_and_probability(self):
        s = load_json(COMP / "bounds_summary.json")
        self.assertAlmostEqual(s["paired_P_arith_tighter"], 0.7684, delta=1e-3)
        self.assertAlmostEqual(s["paired_median_ratio_arith_over_ode"], 0.3337,
                               delta=1e-3)
        self.assertAlmostEqual(s["ode_over_arith_at_median"], 3.101, delta=1e-2)
        # the prose quotes the ratio of medians and the median ratio the way round
        # a reader expects (ODE over combinatorial); the artifact stores the
        # reciprocal, so this asserts the relationship, not a restated literal
        self.assertAlmostEqual(1 / s["paired_median_ratio_arith_over_ode"], 3.0,
                               delta=0.05)
        self.assertIn("differ by a factor of 3.1", self.flat)
        self.assertIn("median of the paired ratio is 3.0", self.flat)
        self.assertIn("tighter in 76.8% of draws", self.flat)

    def test_retracted_statistic_is_absent(self):
        """
        PAIRED_MC_TASK.md states that P(arith < ODE) = 0.654 is not a valid
        paired statistic. the earlier draft reported it as 0.65. it must not
        appear here, and neither may the superseded backup ODE median of 1.7e-2.
        """
        for bad, why in (
            ("0.654", "retracted paired statistic"),
            ("= 0.65", "retracted paired statistic, rounded"),
            ("1.7 × 10⁻²", "superseded backup_uniformN ODE MC median"),
        ):
            self.assertNotIn(bad, self.text, f"manuscript contains {bad!r}: {why}")

    def test_headroom_is_reported_at_the_consistent_evaluation_point(self):
        """
        the headline headroom must come from the usage-weighted mean error rate
        this paper derives, NOT from the bottom of the observed window. earlier
        drafts quoted x158 by evaluating at 1e-4 while deriving 6.3e-4 from the
        same data -- an inconsistency in the favourable direction.
        """
        h = load_json(COMP / "headroom_sensitivity_summary.json")
        self.assertAlmostEqual(h["internally_consistent_headroom_P"], 24.8, delta=0.2)
        self.assertAlmostEqual(h["previously_reported_headroom_P"], 158.06, delta=0.1)
        self.assertStated("×25", "headroom at the usage-weighted mean rate")
        # scoped: the ABSTRACT must carry the corrected figure, not the old one
        abstract = flat(section(H_ABSTRACT))
        self.assertIn("×25 headroom", abstract,
                      "the abstract does not state the corrected headroom")
        self.assertNotIn("×158 headroom", abstract,
                         "the abstract states the superseded headroom as fact")
        # the inflated figure may appear only as the value being corrected
        for para in self.text.split("\n\n"):
            if "×158" in para:
                self.assertTrue(
                    any(w in para.lower() for w in
                        ("earlier draft", "window", "quoted", "previously")),
                    "×158 appears without being marked as the superseded figure:\n"
                    + para[:200])

    def test_headroom_range_is_reported_not_a_point(self):
        c = load_json(COMP / "chaperone_availability_summary.json")
        lo, hi = c["headroom_range_over_documented_grid"]
        self.assertAlmostEqual(lo, 1.9, delta=0.15)
        self.assertAlmostEqual(hi, 25.3, delta=0.2)
        self.assertStated("×1.9", "low end of the availability range")
        # scoped: the abstract and the headroom section must carry the range
        for name, heading in (("abstract", H_ABSTRACT),
                              ("the headroom section", H_HEADROOM)):
            sec = flat(section(heading))
            self.assertIn("×1.9", sec,
                          f"{name} reports a point headroom, not the range")
        self.assertIn("roughly one order of magnitude inside", self.flat)
        self.assertNotIn("roughly two orders of magnitude inside", self.flat,
                         "the superseded claim is back")

    def test_the_anchoring_grid_is_retired_not_silently_kept(self):
        """
        the old anchoring grid duplicated the theta sweep and reached its low end
        outside the documented C_tot range. it must not be cited as a table of
        this paper, and the old Table 6/Table 7 numbering must be gone.
        """
        for stale in ("Table 6", "Table 7", "Table S8", "Table S9"):
            self.assertNotIn(stale, flat(),
                             f"the manuscript still cites {stale} from the old "
                             "numbering")
        self.assertTrue((TAB / "Retired_anchoring_grid.tsv").exists())
        src = (SCRIPTS / "08_make_tables.py").read_text()
        self.assertIn("RETIRED as a main table", src)

    def test_margin_is_near_the_supraadditivity_onset(self):
        h = load_json(COMP / "headroom_sensitivity_summary.json")
        s = load_json(COMP / "supraadditivity_summary.json")
        self.assertAlmostEqual(h["internally_consistent_margin_log10"], 1.39,
                               delta=0.02)
        # the gap between the operating margin and the onset is the reason the
        # paper says the margin "matters"; assert it from the artifacts rather
        # than from a restated literal
        onset = s["qualitative_supraadditivity_first_seen_at_margin_log10"]
        self.assertAlmostEqual(
            h["internally_consistent_margin_log10"] - onset, 0.20, delta=0.02)
        self.assertIn("1.39 log₁₀", self.flat)
        self.assertIn("1.19 log₁₀", self.flat)

    def test_aggregation_death_dominates(self):
        s = load_json(COMP / "bounds_summary.json")
        self.assertAlmostEqual(s["mechanism_frac_aggregation_death"], 0.9994,
                               delta=1e-4)
        self.assertStated("99.94%", "fraction of draws closed by aggregation-death")

    # ---------- axis structure ----------
    def test_mu_axis_clustering(self):
        t = pd.read_csv(COMP / "axis_tests.tsv", sep="\t")
        wd = t[(t.axis == "mu") & (t.null == "within_degeneracy")].iloc[0]
        fu = t[(t.axis == "mu") & (t.null == "full")].iloc[0]
        self.assertAlmostEqual(wd.z, -3.56, delta=0.02)
        self.assertAlmostEqual(fu.z, -2.87, delta=0.02)
        self.assertLess(wd.p_one_sided, 0.005)
        self.assertStated("z = −3.56", "mu within-degeneracy z")
        self.assertStated("z = −2.87", "mu full-shuffle z")

    def test_nu_axis_is_null(self):
        t = pd.read_csv(COMP / "axis_tests.tsv", sep="\t")
        wd = t[(t.axis == "nu") & (t.null == "within_degeneracy")].iloc[0]
        fu = t[(t.axis == "nu") & (t.null == "full")].iloc[0]
        self.assertAlmostEqual(wd.z, -0.06, delta=0.02)
        self.assertAlmostEqual(fu.z, +0.47, delta=0.02)
        self.assertGreater(wd.p_one_sided, 0.05, "nu must be non-significant")
        self.assertGreater(fu.p_one_sided, 0.05, "nu must be non-significant")
        self.assertStated("z = −0.06", "nu within-degeneracy z")
        self.assertStated("z = +0.47", "nu full-shuffle z")

    def test_mu_median_sensitivity(self):
        t = pd.read_csv(COMP / "axis_tests_median.tsv", sep="\t")
        wd = t[(t.axis == "mu") & (t.null == "within_degeneracy")].iloc[0]
        self.assertAlmostEqual(wd.z, -2.59, delta=0.02)
        self.assertLess(wd.p_one_sided, 0.05,
                        "the mu result must survive the median sensitivity")
        self.assertStated("z = −2.59", "mu z under median summarization")

    def test_variance_decomposition(self):
        v = load_json(COMP / "mu_variance_decomposition.json")
        self.assertAlmostEqual(v["eta_squared_log_mu_between_aa"], 0.556, delta=0.002)
        self.assertEqual(v["n_amino_acids"], 18)
        self.assertEqual(v["n_codons"], 59)
        self.assertStated("η² = 0.556", "between-amino-acid variance fraction")
        self.assertStated("59 codons in 18 amino acids", "analysed codon set")

    def test_the_2D_space_is_not_interpreted(self):
        """
        Delta in combined (mu, nu) space mixes a structured and an unstructured
        axis, so its value tracks the mixing weights. the table lists it; the
        paper must decline to read anything into it.
        """
        sec = flat(section(H_AXES))
        self.assertIn("We do not interpret it", sec)
        self.assertIn("weighted mixture", sec)

    # ---------- Fig 1 parameters ----------
    def test_ode_figure_parameters_are_declared_illustrative(self):
        self.assertStated("ρ = 4.0, χ = 0.15", "Fig 1 parameters")
        self.assertStated("λ_crit = 4.80", "Fig 1 fold location")
        self.assertTrue(
            re.search(r"[Pp]arameters are illustrative and not\s+estimated from data",
                      self.text),
            "the manuscript must state that Fig 1 parameters are not fitted")


class ExcludedAnalysesAreDocumented(unittest.TestCase):
    """
    the two dropped analyses must stay reported, with reproducible numbers. the
    v2 restructuring removed this section; it is restored because a reader of the
    repository can see scripts/07 and the Excluded_* tables and is entitled to
    know why they are not results.
    """

    def setUp(self):
        self.text = manuscript()
        self.sec = flat(section(H_EXCLUDED))

    def test_metal_site_collapses_under_within_gene_background(self):
        df = pd.read_csv(COMP / "removed_metal_site_test.tsv", sep="\t")
        self.assertEqual(len(df), 4)
        # every genome-wide test significant, no within-gene test significant
        self.assertTrue((df.genomewide_p < 0.05).all(),
                        "genome-wide background should reproduce the old result")
        self.assertTrue((df.within_gene_p > 0.05).all(),
                        "within-gene background should kill the old result")
        his = df[df.amino_acid == "H"].iloc[0]
        self.assertLess(his.within_gene_or, 1.15)
        self.assertGreater(his.within_gene_p, 0.5)
        for s in ("OR 1.19, p = 0.17", "OR 1.22, p = 0.14",
                  "OR 1.26, p = 0.14", "OR 1.07, p = 0.62"):
            self.assertIn(s, self.sec, f"within-gene result {s!r} not stated")

    def test_site_position_concordance_is_reported(self):
        d = load_json(COMP / "removed_metal_site_diagnostics.json")
        self.assertAlmostEqual(d["position_concordance_frac"], 0.599, delta=0.01)
        # the codon lookup itself is fine; the position mapping is the bug
        self.assertGreater(d["annotation_codon_matches_cds_codon_frac"], 0.95)
        self.assertIn("~60%", self.sec)

    def test_crossspecies_artifacts_are_reported(self):
        s = load_json(COMP / "removed_crossspecies_test.json")
        shared = s["published_delta_A_reproduced_using_ecoli_mu_for_all_species"]
        self.assertTrue(shared["exact"],
                        "shared-mu reproduction of the published table must be exact")
        ind = s["tai_vector_independence"]
        self.assertEqual(ind["ecoli_vs_scerevisiae_bit_identical"], 44)
        self.assertEqual(ind["n_compared"], 60)
        self.assertIn("44 of the 60", self.sec)
        self.assertIn("E. coli μ for all three species", self.sec)

    def test_they_are_not_presented_as_results_of_this_paper(self):
        self.assertIn("failed verification and are excluded", self.sec)
        # and they must not be cited as numbered tables of the paper
        self.assertNotIn("Table 4)", self.sec)
        self.assertNotIn("Table S4)", self.sec)


class TaiAxisIsValidated(unittest.TestCase):
    """regression guard: the corrupt supply vector must never come back."""

    def test_validation_report_discriminates(self):
        rep = load_json(COMP / "tai_validation_report.json")
        by = {r["vector"]: r for r in rep}
        ref = [v for k, v in by.items() if k.startswith("reference")][0]
        bad = [v for k, v in by.items() if k.startswith("CORRUPT")][0]
        self.assertEqual(ref["verdict"], "PASS")
        self.assertEqual(bad["verdict"], "FAIL",
                         "the known-corrupt tAI vector must fail validation")

    def test_validated_vector_is_biologically_sane(self):
        nu = pd.read_csv(COMP / "nu_tai_ecoli_validated.tsv", sep="\t")
        nu = dict(zip(nu.codon, nu.nu_tai))
        # CTG is E. coli's dominant Leu codon and must not be near the floor;
        # this is the single check the corrupt vector failed most visibly
        self.assertGreater(nu["CTG"], 0.5, "CTG supply weight is implausibly low")
        self.assertGreater(nu["CTG"], nu["CTA"], "CTG must outrank rare CTA")
        # AAA is the preferred Lys codon, AAG the minor one
        self.assertGreater(nu["AAA"], nu["AAG"])
        # GAA is the preferred Glu codon
        self.assertGreater(nu["GAA"], nu["GAG"])

    def test_corrupt_vector_is_retained_for_the_regression_test(self):
        p = RAW / "ecoli_tai_CORRUPT_for_regression_test.tsv"
        self.assertTrue(p.exists(),
                        "the corrupt vector must stay on disk so the guard runs")
        bad = pd.read_csv(p, sep="\t")
        bad = dict(zip(bad[bad.columns[0]], bad[bad.columns[1]]))
        self.assertLess(bad["CTG"], 0.1,
                        "this file is supposed to be the corrupt one")

    def test_the_supply_proxy_is_labelled_a_proxy(self):
        """tAI is gene copy number, not a measured elongation rate."""
        self.assertIn("supply proxy derived from tRNA gene copy number",
                      flat(section(H_LIMITS)))


class CitationIntegrity(unittest.TestCase):
    """
    the earlier draft had 20 of 36 references uncited and pointed its
    error-minimization citation at a chaperone protocol. numbering is checked
    mechanically here. v2 switched to parenthetical citations, so the parser
    accepts only parentheticals made entirely of numbers and separators --
    "(1, 2)" and "(23-26)" are citations, "(5,000 draws, seed 17)" is not.
    """

    def setUp(self):
        self.body, self.refs = manuscript().split("## References")

    def cited(self):
        out = set()
        for m in re.finditer(r"\((\d+(?:\s*[,–—-]\s*\d+)*)\)", self.body):
            inner = m.group(1)
            if re.search(r"[–—-]", inner):
                a, b = (int(x) for x in re.split(r"\s*[–—-]\s*", inner))
                out.update(range(a, b + 1))
            else:
                out.update(int(x) for x in inner.split(","))
        return out

    def listed(self):
        return {int(m.group(1))
                for m in re.finditer(r"^(\d+)\.\s", self.refs, re.M)}

    def test_the_citation_parser_sees_citations(self):
        """guard the guard: a parser that matches nothing passes vacuously."""
        self.assertGreater(len(self.cited()), 20,
                           "the citation parser found almost nothing, so the "
                           "tests below would pass vacuously")

    def test_citation_style_is_uniform(self):
        """
        v1 used [1,2] and v2 uses (1, 2). a leftover bracket citation is invisible
        to a reader but means the parser above silently ignores it, so a reference
        cited only in bracket form would read as uncited.
        """
        stray = re.findall(r"\[\d[\d,\s–—-]*\]", self.body)
        self.assertEqual(stray, [],
                         f"bracket-style citations left in the body: {stray}")

    def test_reference_numbering_is_contiguous(self):
        listed = self.listed()
        self.assertEqual(sorted(listed), list(range(1, len(listed) + 1)),
                         "reference numbers must run 1..N with no gaps")

    def test_every_reference_is_cited(self):
        uncited = sorted(self.listed() - self.cited())
        self.assertEqual(uncited, [], f"references never cited in text: {uncited}")

    def test_every_citation_resolves(self):
        dangling = sorted(self.cited() - self.listed())
        self.assertEqual(dangling, [],
                         f"citations with no reference entry: {dangling}")

    def test_all_references_carry_an_identifier(self):
        missing = [line.strip()[:70]
                   for line in self.refs.splitlines()
                   if re.match(r"^\d+\.\s", line.strip())
                   and not re.search(r"PMID|PMC|doi", line, re.I)]
        self.assertEqual(missing, [], f"references without PMID/PMC/DOI: {missing}")

    def test_no_unresolved_placeholders(self):
        for bad in ("VERIFY", "in prep", "preprint, 2025", "[repository URL]"):
            self.assertNotIn(bad, manuscript(),
                             f"unresolved placeholder {bad!r} in manuscript")


class ArtifactsExist(unittest.TestCase):
    """figures must exist, be generated by a script, and be embedded."""

    STEMS = ("Fig1_envelope", "Fig2_axis_structure", "Fig3_bounds",
             "Fig4_supraadditivity")

    def test_figures_present(self):
        for stem in self.STEMS:
            for ext in ("png", "pdf"):
                p = FIGS / f"{stem}.{ext}"
                self.assertTrue(p.exists(), f"missing figure {p.name}")

    def test_every_figure_has_a_generating_script(self):
        src = "\n".join(p.read_text() for p in SCRIPTS.glob("*.py"))
        for stem in self.STEMS:
            self.assertIn(stem, src, f"no script in scripts/ generates {stem}")

    def test_every_figure_is_embedded_with_a_legend(self):
        text = manuscript()
        for i, stem in enumerate(self.STEMS, start=1):
            self.assertIn(f"](../figures/{stem}.png)", text,
                          f"{stem} is not embedded in the manuscript")
        for n in (1, 2, 3, 4):
            self.assertIn(f"**Fig. {n}.", text, f"Fig. {n} has no legend")

    def test_computed_outputs_present(self):
        for name in ("nu_tai_ecoli_validated.tsv", "tai_validation_report.json",
                     "axis_tests.tsv", "axis_tests_median.tsv",
                     "mu_variance_decomposition.json", "translation_burden.json",
                     "bounds_summary.json", "removed_metal_site_test.tsv",
                     "removed_crossspecies_test.json",
                     "supraadditivity_summary.json",
                     "supraadditivity_margin_sweep.tsv",
                     "supraadditivity_effect_grid.tsv",
                     "supraadditivity_collapse_frontier.tsv",
                     "supraadditivity_knob_comparison.tsv",
                     "headroom_sensitivity.tsv",
                     "headroom_sensitivity_summary.json",
                     "chaperone_availability.tsv",
                     "chaperone_availability_summary.json",
                     "nu_power_sweep.tsv", "nu_power_curve.tsv",
                     "nu_power_summary.json", "mu_jackknife.tsv",
                     "mu_depth_sensitivity.tsv", "mu_jackknife_summary.json"):
            self.assertTrue((COMP / name).exists(), f"missing artifact {name}")

    def test_no_triplet_origin_claim(self):
        """the project's standing constraint: this is not a code-origin paper."""
        text = manuscript()
        self.assertTrue(
            "no code-origin claim" in text or "no claim about the origin" in text,
            "the manuscript must explicitly disclaim a code-origin reading")


class Supraadditivity(unittest.TestCase):
    """
    the framework's distinguishing prediction, tested in its own model. the point
    of these assertions is that the result must stay anchored at the paper's own
    operating point -- evaluating it anywhere else changes the answer by 38x.
    """

    def setUp(self):
        self.text = manuscript()
        self.flat = flat()
        self.s = load_json(COMP / "supraadditivity_summary.json")

    def test_interaction_is_supraadditive_never_sub(self):
        t = pd.read_csv(COMP / "supraadditivity_effect_grid.tsv", sep="\t")
        defined = t[t.interaction.notna()]
        self.assertGreater(len(defined), 0)
        self.assertTrue((defined.interaction >= -1e-9).all(),
                        "the model must not produce subadditive interactions")

    def test_evaluated_at_the_internally_consistent_point(self):
        """
        THE defect this class exists to prevent. scripts/09 hardcoded f = 1e-4 --
        the window bottom, the evaluation point the headroom section spends a page
        rejecting -- so the factorial was anchored at margin 2.20 while the paper's
        own margin is 1.39, and reported the interaction as +0.2% instead of
        +7.4%. The error ran AGAINST the paper: it argued a real result away. And
        this suite pinned it in place by asserting the string "+0.2% of the
        additive expectation".
        """
        b = load_json(COMP / "translation_burden.json")
        self.assertAlmostEqual(self.s["evaluation_point_f_codon"],
                               b["usage_weighted_mean_mu_per_codon"], places=12,
                               msg="the factorial is not evaluated at the point "
                                   "the burden analysis establishes")
        h = load_json(COMP / "headroom_sensitivity_summary.json")
        self.assertAlmostEqual(self.s["baseline_margin_log10"],
                               h["internally_consistent_margin_log10"], delta=0.01,
                               msg="the factorial and the headroom analysis "
                                   "disagree on the margin")

    def test_the_evaluation_point_is_not_hardcoded(self):
        """a literal 1e-4 anchor is how this went wrong; keep it read-from-disk."""
        src = (SCRIPTS / "09_supraadditivity.py").read_text()
        self.assertIn("def consistent_f()", src)
        self.assertIn("translation_burden.json", src)
        self.assertNotIn("OBSERVED_F = 1e-4", src,
                         "the evaluation point is hardcoded again")

    def test_effect_size_at_the_evaluation_point(self):
        pct = self.s["interaction_pct_at_evaluation_point"]
        self.assertAlmostEqual(pct, 7.44, delta=0.1)
        self.assertGreater(pct, 5.0,
                           "the effect at E. coli's own margin is not negligible")
        self.assertIn("7.4% at *E. coli*'s own margin", self.flat)
        self.assertIn("9.6%", self.flat)
        for bad in ("+0.2% of the additive expectation", "not a detectable effect",
                    "far too small to measure"):
            self.assertNotIn(bad, self.flat,
                             f"superseded conclusion {bad!r} is back")

    def test_largest_viable_excess_matches_its_grid_cell(self):
        """9.6% is a specific cell (error x2, capacity /5); check it is that one."""
        g = pd.read_csv(COMP / "supraadditivity_effect_grid.tsv", sep="\t")
        viable = g[g.interaction.notna()]
        top = viable.loc[viable.interaction_pct_of_additive.idxmax()]
        self.assertAlmostEqual(top.interaction_pct_of_additive, 9.59, delta=0.05)
        self.assertEqual((top.error_factor, top.capacity_factor), (2.0, 5.0))
        self.assertIn("2-fold error increase with a 5-fold knockdown", self.flat)

    def test_window_bottom_is_kept_only_as_a_comparison(self):
        self.assertAlmostEqual(self.s["interaction_pct_at_window_bottom"], 0.196,
                               delta=0.01)
        self.assertAlmostEqual(self.s["baseline_margin_log10_at_window_bottom"],
                               2.199, delta=0.01)
        # where 0.2% appears it must be marked as the window bottom, not E. coli
        for para in self.text.split("\n\n"):
            if "0.2%" in para:
                self.assertTrue(
                    any(w in para.lower() for w in
                        ("window bottom", "earlier draft")),
                    "0.2% appears unmarked as the window bottom:\n" + para[:220])

    def test_interaction_grows_as_margin_closes(self):
        self.assertTrue(self.s["interaction_grows_as_margin_closes"])
        t = pd.read_csv(COMP / "supraadditivity_margin_sweep.tsv", sep="\t")
        d = t[t.interaction.notna()].sort_values("margin_baseline", ascending=False)
        self.assertTrue(d.interaction_pct_of_additive.is_monotonic_increasing,
                        "interaction should rise monotonically as margin closes")

    def test_margin_ladder_matches_computed(self):
        """the inline ladder in the prose must match the margin sweep."""
        t = pd.read_csv(COMP / "supraadditivity_margin_sweep.tsv", sep="\t")
        for margin, pct in ((2.20, 0.2), (1.90, 0.6), (1.39, 7.4)):
            row = t.iloc[(t.margin_baseline - margin).abs().argmin()]
            self.assertAlmostEqual(row.interaction_pct_of_additive, pct, delta=0.15,
                                   msg=f"margin {margin}: computed "
                                       f"{row.interaction_pct_of_additive:+.2f}%")
            self.assertIn(f"{margin:.2f}", self.flat)

    def test_qualitative_synthetic_lethality(self):
        margin = self.s["qualitative_supraadditivity_first_seen_at_margin_log10"]
        self.assertAlmostEqual(margin, 1.195, delta=0.01)
        self.assertIn("1.19", self.flat)
        self.assertEqual(self.s["n_qualitatively_supraadditive_grid_points"], 114)
        self.assertIn("114 of 676", self.flat)
        # the headline count is over the 36-cell grid at E. coli's own margin
        g = pd.read_csv(COMP / "supraadditivity_effect_grid.tsv", sep="\t")
        self.assertEqual(len(g), 36)
        self.assertEqual(int(g.qualitative_supraadditive.sum()), 12)
        self.assertIn("12 of the 36", self.flat)
        self.assertIn("does not require a sensitized background", self.flat)

    def test_the_mildest_lethal_pair_matches_the_grid(self):
        g = pd.read_csv(COMP / "supraadditivity_effect_grid.tsv", sep="\t")
        sl = g[g.qualitative_supraadditive]
        mildest = sl.loc[(sl.error_factor * sl.capacity_factor).idxmin()]
        self.assertEqual((mildest.error_factor, mildest.capacity_factor),
                         (1.5, 10.0))
        self.assertIn("1.5-fold error increase with a 10-fold rescue knockdown",
                      self.flat)

    def test_collapse_cells_have_no_numeric_interaction(self):
        """
        substituting a lower bound there yields spurious NEGATIVE interactions.
        the tables must leave those cells undefined.
        """
        for f in ("supraadditivity_margin_sweep.tsv",
                  "supraadditivity_effect_grid.tsv"):
            t = pd.read_csv(COMP / f, sep="\t")
            bad = t[t.collapsed_both & t.interaction.notna()]
            self.assertEqual(len(bad), 0,
                             f"{f}: {len(bad)} collapsed cells carry a numeric "
                             "interaction, which cannot be meaningful")

    def test_saturation_numbers_are_at_the_evaluation_point(self):
        """
        the "97.9% / 47.5 uM / 0.052 uM" triple was the window-bottom operating
        point and contradicted the theta = 0 row of the availability sweep. all
        three must now be the evaluation-point values, and must agree with Table 3.
        """
        sat = self.s["folding_arm_saturation_at_observed_rate"]
        self.assertAlmostEqual(sat, 0.974, delta=0.002)
        self.assertAlmostEqual(self.s["free_chaperone_uM_at_observed_rate"], 37.56,
                               delta=0.1)
        self.assertAlmostEqual(self.s["misfolded_protein_uM_at_observed_rate"],
                               0.331, delta=0.005)
        self.assertIn("97.4% saturated", self.flat)
        self.assertIn("37.6 µM free", self.flat)
        self.assertIn("0.331 µM", self.flat)
        # the excess is a derived quantity; assert it against the two it derives
        # from rather than trusting the prose
        excess = (self.s["free_chaperone_uM_at_observed_rate"]
                  / self.s["misfolded_protein_uM_at_observed_rate"])
        self.assertAlmostEqual(excess, 113, delta=1)
        self.assertIn("113-fold excess", self.flat)
        for bad in ("97.9% saturated", "47.5 µM free", "0.052 µM misfolded",
                    "900-fold"):
            self.assertNotIn(bad, self.flat,
                             f"window-bottom figure {bad!r} is back")

        # and it must agree with Table 3's theta = 0, C_tot = 50 row
        t3 = pd.read_csv(TAB / "Table3_chaperone_availability.tsv", sep="\t")
        row = t3[(t3.theta == 0.0) & (t3.C_tot_uM == 50.0)
                 & (t3.K_d_uM == 1.0)].iloc[0]
        self.assertAlmostEqual(row.folding_arm_saturation, sat, delta=0.002,
                               msg="the factorial and Table 3 disagree on saturation")
        self.assertAlmostEqual(
            row.c_free_uM, self.s["free_chaperone_uM_at_observed_rate"], delta=0.1,
            msg="the factorial and Table 3 disagree on free chaperone")
        # it must be named as a limitation, not just a methods note -- but as the
        # STRUCTURAL omission it is, not as the parameter criticism an earlier
        # pass wrongly made (see ChaperoneAvailability)
        lim = flat(section(H_LIMITS))
        self.assertIn("omits the nascent-chain folding load", lim,
                      "the saturation problem must appear in Limitations, framed "
                      "as the structural omission it is")
        self.assertIn("parameters are not at fault", lim)

    def test_the_two_knobs_are_distinguished_with_their_computed_costs(self):
        k = pd.read_csv(COMP / "supraadditivity_knob_comparison.tsv", sep="\t")
        at = k[k.evaluation_point == "usage_weighted_mu"]
        c_tot = at[at.knob == "C_tot"].iloc[0].D_capacity
        k_obs = at[at.knob == "k_obs_max"].iloc[0].D_capacity
        self.assertAlmostEqual(c_tot, 0.022, delta=0.001)
        self.assertAlmostEqual(k_obs, 0.460, delta=0.001)
        self.assertLess(c_tot * 10, k_obs, "the knobs should differ by ~20x")
        self.assertIn("costs only 0.022 log₁₀", self.flat)
        self.assertIn("costs 0.460", self.flat)

    def test_model_only_caveat_is_stated(self):
        self.assertIn("not whether cells do", self.flat)
        self.assertIn("cannot validate the framework", self.flat)

    def test_predictions_were_rewritten_not_left_stale(self):
        """prediction 1 must no longer promise a wild-type additivity deviation."""
        pred = flat(section(H_PREDICTIONS))
        self.assertIn("Synthetic lethality in wild type", pred)
        self.assertIn("no sensitization is required", pred)
        self.assertNotIn("synthetic lethality in a sensitized background", self.flat,
                         "the superseded sensitized-background framing is back")
        self.assertNotIn("Combining a synonymous change that raises folding "
                         "burden with a reduction in chaperone capacity",
                         self.flat)

    def test_vendored_model_is_unmodified(self):
        a = (SCRIPTS / "vendor" / "two_pool_ode.py").read_bytes()
        b = Path("/storage/kiran-stuff/proteostasis-P1/two_pool_ode.py")
        if b.exists():
            self.assertEqual(a, b.read_bytes(),
                             "the vendored model has drifted from upstream")


class ChaperoneAvailability(unittest.TestCase):
    """
    the folding arm's saturation is a structural assumption, not a bad parameter.
    these tests keep that distinction -- and keep theta labelled as unmeasured.
    """

    def setUp(self):
        self.text = manuscript()
        self.flat = flat()
        self.s = load_json(COMP / "chaperone_availability_summary.json")

    def test_parameters_are_stated_to_be_in_sourced_ranges(self):
        prov = self.s["parameter_provenance"]
        self.assertIn("Lorimer", prov["C_tot_uM"])
        self.assertIn("Pierpaoli", prov["K_d_uM"])
        # and the manuscript must say the parameters are not at fault
        self.assertIn("used within their sourced ranges", self.flat)
        self.assertIn("30 to 80 µM", self.flat)
        self.assertIn("0.06 to 2 µM", self.flat)

    def test_the_earlier_overstatement_is_gone(self):
        """
        an earlier pass claimed the parameterization contradicted the paper's own
        capacity evidence. it does not -- the omission of nascent-chain folding
        does. that correction must not silently revert.
        """
        for bad in ("understates how tightly", "sits awkwardly against"):
            self.assertNotIn(bad, self.flat,
                             f"the superseded overstatement {bad!r} is back")
        self.assertIn("It is structural.", self.text)
        self.assertIn("does not represent the ordinary nascent-chain folding load",
                      self.flat)

    def test_theta_is_labelled_unmeasured(self):
        self.assertIn("NOT measured", self.s["parameter_provenance"]["theta"])
        self.assertIn("θ is not measured here", self.flat)
        self.assertIn("θ is unmeasured", flat(section(H_LIMITS)))

    def test_theta_thresholds(self):
        self.assertAlmostEqual(self.s["theta_at_which_arm_unsaturates"], 0.98,
                               delta=1e-6)
        self.assertAlmostEqual(
            self.s["theta_at_which_margin_reaches_supraadditivity_onset"], 0.90,
            delta=1e-6)
        self.assertIn("θ ≈ 0.98", self.flat)
        self.assertIn("θ ≈ 0.90", self.flat)

    def test_theta_zero_agrees_with_the_headroom_analysis(self):
        h = load_json(COMP / "headroom_sensitivity_summary.json")
        z = self.s["at_theta_zero"]
        self.assertAlmostEqual(z["headroom_P"],
                               h["internally_consistent_headroom_P"], delta=0.1)
        self.assertAlmostEqual(z["folding_arm_saturation"], 0.974, delta=0.002)

    def test_headroom_range_over_documented_grid(self):
        lo, hi = self.s["headroom_range_over_documented_grid"]
        self.assertAlmostEqual(lo, 1.9, delta=0.15)
        self.assertAlmostEqual(hi, 25.3, delta=0.2)
        self.assertIn("×1.9 to ×25", self.flat)
        # Limitations must carry the range, not a point
        self.assertIn("×1.9–×25", flat(section(H_LIMITS)))

    def test_the_pinning_measurement_is_named(self):
        self.assertIn("occupancy", self.s["what_would_pin_this_down"])
        self.assertIn("chaperone occupancy by nascent-chain folding", self.flat)


class AxisPowerAndRobustness(unittest.TestCase):
    """
    the two objections a reviewer raises about the axis result: an unpowered nu
    null, and a mu signal that may be leverage or mass-spec detectability rather
    than biology.
    """

    def setUp(self):
        self.flat = flat()
        self.pw = load_json(COMP / "nu_power_summary.json")
        self.jk = load_json(COMP / "mu_jackknife_summary.json")

    # ---------- nu resolution and power ----------
    def test_resolution_asymmetry_is_reported(self):
        t = self.pw["tie_structure"]
        self.assertEqual(t["nu"]["n_distinct"], 21)
        self.assertEqual(t["log_mu"]["n_distinct"], 59)
        self.assertEqual(t["nu"]["largest_tied_group"], 13)
        self.assertIn("21 of 59 codons carry distinct ν values", self.flat)
        self.assertIn("largest tied group spanning 13", self.flat)

    def test_mu_effect_exceeds_the_nu_detection_floor(self):
        c = self.pw["minimum_detectable_effect"]["comparison"]
        self.assertTrue(c["mu_effect_exceeds_nu_floor"],
                        "if this fails, the axis contrast is a resolution "
                        "artifact and must be withdrawn")
        self.assertAlmostEqual(c["mu_observed_pct_below_null"], 37.4, delta=0.5)
        self.assertAlmostEqual(c["nu_detection_floor_pct_below_null"], 19.6,
                               delta=0.5)
        self.assertAlmostEqual(c["margin_percentage_points"], 17.8, delta=0.6)
        # and the margin must be the difference of the two the paper quotes
        self.assertAlmostEqual(c["mu_observed_pct_below_null"]
                               - c["nu_detection_floor_pct_below_null"],
                               c["margin_percentage_points"], places=6)
        self.assertIn("37.4% below its null mean", self.flat)
        self.assertIn("19.6% below its", self.flat)
        self.assertIn("17.8 percentage points", self.flat)

    def test_the_subset_power_limit_is_conceded(self):
        """the floor assumes a uniform effect; the paper must say so."""
        self.assertIsNone(self.pw["s_at_80_percent_power"],
                          "80% power was reached; update the concession")
        self.assertIn("uniformly across amino acids", self.flat)
        self.assertIn("is not excluded", self.flat)
        # the blanket claim must not be asserted anywhere
        self.assertNotIn("no structure on the supply axis at all", self.flat.lower())
        self.assertIn("narrower than a blanket absence", self.flat)
        # a figure title is a claim too, and it is the one place a blanket
        # absence can slip past prose review
        fig = (SCRIPTS / "04_fig2_axis.py").read_text()
        self.assertNotIn('"No structure on supply"', fig,
                         "Fig 2b asserts a blanket absence the paper disclaims")

    def test_nu_p_values_are_labelled_descriptive(self):
        self.assertIn("descriptive rather than evidential", self.flat)
        self.assertIn("would be ≈0.68", self.flat)

    # ---------- mu leverage and detectability ----------
    def test_jackknife_survives_every_deletion(self):
        j = self.jk["jackknife"]
        self.assertTrue(j["all_remain_significant"])
        self.assertEqual(j["n_significant"], 59)
        self.assertEqual(j["n_total"], 59)
        self.assertIn("59 of 59 deletions", self.flat)
        self.assertAlmostEqual(j["z_min"], -3.86, delta=0.01)
        self.assertAlmostEqual(j["z_max"], -2.62, delta=0.01)
        self.assertIn("z between −3.86 and −2.62", self.flat)

    def test_CCC_does_not_carry_the_result(self):
        j = self.jk["jackknife"]
        self.assertLess(j["p_without_CCC"], 0.05)
        self.assertAlmostEqual(j["z_without_CCC"], -3.67, delta=0.1)
        self.assertIn("z = −3.67", self.flat)
        # and the span must be reported both ways
        self.assertAlmostEqual(self.jk["mu_span_fold_without_CCC"], 286, delta=2)
        self.assertIn("286-fold", self.flat)

    def test_detectability_confound_is_disclosed_with_numbers(self):
        d = self.jk["detectability"]
        self.assertAlmostEqual(d["spearman_depth_vs_log_mu"], -0.366, delta=0.01)
        self.assertAlmostEqual(d["eta2_between_aa_sampling_depth"], 0.560,
                               delta=0.005)
        # the damaging comparison: depth carries as much AA structure as mu
        self.assertAlmostEqual(d["eta2_between_aa_sampling_depth"],
                               d["eta2_between_aa_log_mu"], delta=0.02)
        self.assertIn("η² = 0.560 between amino acids for depth against 0.556",
                      self.flat)
        self.assertIn("cannot be separated from amino-acid-level structure in "
                      "detectability", self.flat)
        self.assertIn("1 to 16 detected substitutions (median 10)", self.flat)
        self.assertIn("Spearman ρ = −0.37, p = 0.004", self.flat)
        # and it must be conceded in Limitations, not only in Results
        self.assertIn("detectability", flat(section(H_LIMITS)))

    def test_clustering_survives_dropping_thin_codons(self):
        thin = {r["drop_below_percentile"]: r
                for r in self.jk["detectability"]["depth_sensitivity"]}
        self.assertLess(thin[25]["p"], 0.05)
        self.assertIn("z = −2.85, p = 0.007", self.flat)


class TablesAgreeWithData(unittest.TestCase):
    """
    the generated tables must match data/computed/, the manuscript must match the
    generated tables, and the old numbering must be gone from disk. drift between
    a table and the prose that describes it is the failure mode that sank the
    previous draft.
    """

    MAIN = {"Table1_burden_terms": "Table 1",
            "Table2_bounds": "Table 2",
            "Table3_chaperone_availability": "Table 3",
            "Table4_supraadditivity": "Table 4",
            "Table5_axis_tests": "Table 5"}
    SUPP = ["TableS1_codon_coordinates", "TableS2_delta_per_aa",
            "TableS3_tai_validation", "TableS4_supraadditivity_grid",
            "TableS5_capacity_knob_comparison", "TableS6_axis_power",
            "TableS7_mu_jackknife"]
    NOT_PAPER = ["Excluded_metal_site_backgrounds", "Excluded_crossspecies",
                 "Retired_anchoring_grid"]
    # written by the pre-v2 numbering; a stale Table 7 sitting in tables/ is
    # indistinguishable from a live one
    SUPERSEDED = ["Table3_axis_tests", "Table4_metal_site_backgrounds",
                  "Table5_supraadditivity", "Table7_chaperone_availability",
                  "TableS4_crossspecies", "TableS5_supraadditivity_grid",
                  "TableS6_capacity_knob_comparison", "TableS7_axis_power",
                  "TableS8_mu_jackknife", "TableS9_retired_anchoring_grid"]

    def test_all_tables_written(self):
        for name in (list(self.MAIN) + ["Table2_bounds_statistics"]
                     + self.SUPP + self.NOT_PAPER):
            self.assertTrue((TAB / f"{name}.tsv").exists(),
                            f"missing tables/{name}.tsv")
        self.assertTrue((TAB / "TABLES.md").exists(), "missing tables/TABLES.md")

    def test_the_old_numbering_is_removed_from_disk(self):
        left = [n for n in self.SUPERSEDED if (TAB / f"{n}.tsv").exists()]
        self.assertEqual(left, [],
                         f"files from the superseded numbering are still in "
                         f"tables/: {left}")

    def test_table2_matches_bounds_summary(self):
        b = load_json(COMP / "bounds_summary.json")
        t = pd.read_csv(TAB / "Table2_bounds.tsv", sep="\t")
        got = dict(zip(t.quantity, t.value_per_codon))
        self.assertAlmostEqual(got["Combinatorial bound, paired MC median"],
                               b["arithmetic_paired_median"], places=12)
        self.assertAlmostEqual(got["Two-pool ODE bound, paired MC median"],
                               b["ode_paired_median"], places=12)

        s = pd.read_csv(TAB / "Table2_bounds_statistics.tsv", sep="\t")
        stats = dict(zip(s.statistic, s.value))
        self.assertAlmostEqual(stats["P(combinatorial is the tighter bound)"],
                               b["paired_P_arith_tighter"], places=12)

    def test_table3_matches_the_availability_sweep(self):
        t = pd.read_csv(TAB / "Table3_chaperone_availability.tsv", sep="\t")
        src = pd.read_csv(COMP / "chaperone_availability.tsv", sep="\t")
        self.assertEqual(len(t), len(src))
        for col in ("theta", "headroom_P", "folding_arm_saturation"):
            pd.testing.assert_series_equal(t[col], src[col], check_names=False)

    def test_table5_matches_axis_tests(self):
        t = pd.read_csv(TAB / "Table5_axis_tests.tsv", sep="\t")
        src = pd.read_csv(COMP / "axis_tests.tsv", sep="\t")
        for _, r in src.iterrows():
            m = t[(t.axis == r.axis) & (t.null == r.null)
                  & (t.mu_stat.isin(["mean", "both"]))]
            self.assertEqual(len(m), 1,
                             f"Table 5 has {len(m)} rows for {r.axis}/{r.null}")
            self.assertAlmostEqual(float(m.iloc[0].z), float(r.z), places=10)

    def test_table5_does_not_duplicate_the_nu_rows(self):
        """
        nu does not depend on the mu summary statistic, so listing it twice would
        present one test as two.
        """
        t = pd.read_csv(TAB / "Table5_axis_tests.tsv", sep="\t")
        nu = t[t.axis == "nu"]
        self.assertEqual(len(nu), 2, "nu should appear once per null model")
        self.assertTrue((nu.mu_stat == "both").all(),
                        "nu rows must not be attributed to a mu statistic")

    def test_table5_labels_null_results_as_null(self):
        t = pd.read_csv(TAB / "Table5_axis_tests.tsv", sep="\t")
        ns = t[t.p_one_sided > 0.05]
        self.assertGreater(len(ns), 0, "expected at least one null row")
        self.assertTrue((ns.direction == "no signal").all(),
                        "non-significant rows must not be labelled "
                        "clustered/spread on the strength of their sign")

    def test_excluded_table_matches_the_recomputed_metal_site_test(self):
        t = pd.read_csv(TAB / "Excluded_metal_site_backgrounds.tsv", sep="\t")
        src = pd.read_csv(COMP / "removed_metal_site_test.tsv", sep="\t")
        self.assertEqual(len(t), len(src))
        for aa in src.amino_acid:
            a = t[t.amino_acid == aa].iloc[0]
            b = src[src.amino_acid == aa].iloc[0]
            self.assertAlmostEqual(a.within_gene_or, b.within_gene_or, places=10)
            self.assertAlmostEqual(a.within_gene_p, b.within_gene_p, places=10)

    def test_supplementary_table_sizes(self):
        s1 = pd.read_csv(TAB / "TableS1_codon_coordinates.tsv", sep="\t")
        s2 = pd.read_csv(TAB / "TableS2_delta_per_aa.tsv", sep="\t")
        v = load_json(COMP / "mu_variance_decomposition.json")
        self.assertEqual(len(s1), v["n_codons"], "S1 must list every analysed codon")
        self.assertEqual(len(s2), v["n_amino_acids"],
                         "S2 must list every multi-codon amino acid")

    def test_manuscript_matches_generated_tables(self):
        """
        the values the manuscript states in prose must be the ones the table
        files contain, formatted identically.
        """
        text = manuscript()
        md = (TAB / "TABLES.md").read_text()
        for needle in ("8.35 × 10⁻³", "2.59 × 10⁻²", "1.19 × 10⁻³",
                       "1.00 × 10⁻²", "z = −3.56"):
            bare = needle.replace("z = ", "")
            self.assertIn(bare, md,
                          f"{bare!r} is in the manuscript but not in TABLES.md")
            self.assertIn(needle, text,
                          f"{needle!r} is in TABLES.md but not in the manuscript")

    def test_table_captions_are_not_stale(self):
        """
        the captions quoted the window-bottom saturation triple and "f = 1e-4"
        for a full commit after the manuscript had corrected both. they are now
        interpolated from the data; this asserts they stay that way.
        """
        md = (TAB / "TABLES.md").read_text()
        for bad in ("97.9% saturated", "47.5 µM free", "0.052 µM misfolded",
                    "combinations at f = 1e-4"):
            self.assertNotIn(bad, md, f"a stale caption value is back: {bad!r}")
        s = load_json(COMP / "supraadditivity_summary.json")
        self.assertIn("97.4% saturated", md)
        self.assertIn("6.33 × 10⁻⁴", md)
        self.assertAlmostEqual(s["evaluation_point_f_codon"], 6.334e-4, delta=1e-6)
        src = (SCRIPTS / "08_make_tables.py").read_text()
        self.assertIn("sat_pct", src,
                      "caption values must be interpolated, not typed")

    def test_manuscript_cites_the_tables(self):
        text = manuscript()
        for label in list(self.MAIN.values()) + [
                "Table S1", "Table S2", "Table S3", "Table S4", "Table S5",
                "Table S6", "Table S7"]:
            self.assertIn(label, text, f"manuscript never refers to {label}")

    def test_not_paper_tables_are_marked_as_such(self):
        """a reader must not mistake an excluded analysis for a finding."""
        md = (TAB / "TABLES.md").read_text()
        for heading in ("## Excluded analysis A:", "## Excluded analysis B:"):
            i = md.index(heading)
            self.assertIn("Excluded from the paper", md[i:i + 700],
                          f"{heading} does not say the result is excluded")
        i = md.index("## Retired: the evaluation-point")
        self.assertIn("Superseded by Table 3", md[i:i + 400])


class AssembledPaper(unittest.TestCase):
    """
    scripts/15_build_paper.py inlines the generated tables into the prose. the
    point of the placeholder mechanism is that a table body is never typed into
    the manuscript, so these tests check the seam.
    """

    def setUp(self):
        self.ms = manuscript()
        self.placeholders = re.findall(r"<!-- TABLE:(.+?) -->", self.ms)

    def test_every_main_table_has_a_placeholder(self):
        self.assertEqual(sorted(self.placeholders),
                         ["Table 1", "Table 2", "Table 3", "Table 4", "Table 5"])

    def test_the_manuscript_does_not_inline_table_bodies_by_hand(self):
        """
        a hand-typed table row is a copy of a generated number, which is how the
        previous draft went stale. the prose may contain no markdown table rows.
        """
        rows = [ln for ln in self.ms.splitlines()
                if ln.startswith("|") and "---" not in ln]
        self.assertEqual(rows, [],
                         f"MANUSCRIPT.md contains {len(rows)} hand-written table "
                         "rows; use a <!-- TABLE:... --> placeholder instead")

    def test_every_placeholder_resolves(self):
        md = (TAB / "TABLES.md").read_text()
        for label in self.placeholders:
            self.assertIn(f"## {label}.", md,
                          f"no block in TABLES.md for placeholder {label!r}")

    def test_paper_is_assembled_and_complete(self):
        self.assertTrue(PAPER.exists(),
                        "manuscript/PAPER.md missing; run scripts/15_build_paper.py")
        paper = PAPER.read_text()
        self.assertNotIn("<!-- TABLE:", paper,
                         "the assembled paper still has an unresolved placeholder")
        for label in self.placeholders:
            self.assertIn(f"**{label}.", paper,
                          f"{label} is not inlined in the assembled paper")
        for stem in ("Fig1_envelope", "Fig2_axis_structure", "Fig3_bounds",
                     "Fig4_supraadditivity"):
            self.assertIn(stem, paper, f"{stem} missing from the assembled paper")
        self.assertIn("## References", paper)

    def test_the_assembled_tables_are_the_generated_ones(self):
        """spot-check a value that only exists in the generated Table 3."""
        paper = PAPER.read_text()
        t3 = pd.read_csv(TAB / "Table3_chaperone_availability.tsv", sep="\t")
        row = t3[(t3.theta == 0.0) & (t3.C_tot_uM == 50.0)
                 & (t3.K_d_uM == 1.0)].iloc[0]
        self.assertIn(f"×{row.headroom_P:.1f}", paper)
        self.assertIn(f"{row.c_free_uM:.2f}", paper)

    def test_the_build_script_fails_closed(self):
        src = (SCRIPTS / "15_build_paper.py").read_text()
        self.assertIn("sys.exit", src,
                      "the build script must fail on a missing table or figure")


if __name__ == "__main__":
    unittest.main(verbosity=2)
