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

from the paper root. scripts/01 through scripts/07 must have been run first.
"""
import json
import re
import unittest
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
COMP = ROOT / "data" / "computed"
MS = ROOT / "manuscript" / "MANUSCRIPT.md"
FIGS = ROOT / "figures"
SCRIPTS = ROOT / "scripts"
TAB = ROOT / "tables"


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

    def assertStated(self, needle, why):
        self.assertIn(needle, self.text,
                      f"manuscript does not state {needle!r} ({why})")

    # ---------- Result 1: burden magnitude ----------
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
        self.assertStated("16–18%", "fraction of proteins with >=1 error")
        # and the claim the earlier draft made must NOT reappear
        self.assertNotIn("one-fifth to one-quarter", self.text,
                         "the unsupported 20-25% claim is back")

    def test_proteome_description(self):
        b = load_json(COMP / "translation_burden.json")
        self.assertEqual(b["proteome_n"], 4403)
        self.assertEqual(b["median_length_aa"], 271)
        self.assertStated("n = 4,403", "proteome size")
        self.assertStated("median 271 aa", "median protein length")

    # ---------- Result 2/3: bounds ----------
    def test_deterministic_bounds(self):
        s = load_json(COMP / "bounds_summary.json")
        self.assertAlmostEqual(s["arithmetic_deterministic"], 1.188e-3, delta=1e-5)
        self.assertAlmostEqual(s["ode_deterministic"], 1.000e-2, delta=1e-4)
        self.assertStated("1.19 × 10⁻³", "arithmetic deterministic bound")
        self.assertStated("1.00 × 10⁻²", "ODE deterministic bound")

    def test_paired_medians_and_intervals(self):
        s = load_json(COMP / "bounds_summary.json")
        self.assertAlmostEqual(s["arithmetic_paired_median"], 8.348e-3, delta=1e-5)
        self.assertAlmostEqual(s["ode_paired_median"], 2.589e-2, delta=1e-4)
        self.assertStated("8.35 × 10⁻³", "arithmetic paired median")
        self.assertStated("2.59 × 10⁻²", "ODE paired median")
        self.assertStated("1.38 × 10⁻³", "arithmetic CI low")
        self.assertStated("8.79 × 10⁻²", "arithmetic CI high")
        self.assertStated("1.25 × 10⁻³", "ODE CI low")
        self.assertStated("5.87 × 10⁻¹", "ODE CI high")

    def test_paired_ratio_and_probability(self):
        s = load_json(COMP / "bounds_summary.json")
        self.assertAlmostEqual(s["paired_P_arith_tighter"], 0.7684, delta=1e-3)
        self.assertAlmostEqual(s["paired_median_ratio_arith_over_ode"], 0.3337,
                               delta=1e-3)
        self.assertAlmostEqual(s["ode_over_arith_at_median"], 3.101, delta=1e-2)
        self.assertStated("0.768", "valid paired P(arithmetic tighter)")
        self.assertStated("0.334", "median paired ratio")
        self.assertStated("3.10×", "ODE/arithmetic at the median")

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
        abstract = flat(section("## Abstract"))
        self.assertIn("×25 headroom", abstract,
                      "the abstract does not state the corrected headroom")
        self.assertNotIn("×158 headroom", abstract,
                         "the abstract states the superseded headroom as fact")
        # the inflated figure may appear only as the value being corrected
        for para in self.text.split("\n\n"):
            if "×158" in para:
                self.assertTrue(
                    any(w in para.lower() for w in
                        ("earlier draft", "window bottom", "quoted", "favourable",
                         "as published", "previously")),
                    "×158 appears without being marked as the superseded figure:\n"
                    + para[:200])

    def test_headroom_range_is_reported_not_a_point(self):
        h = load_json(COMP / "headroom_sensitivity_summary.json")
        lo, hi = h["usage_weighted_mu_range_across_anchorings"]
        self.assertAlmostEqual(lo, 4.6, delta=0.2)
        self.assertAlmostEqual(hi, 24.8, delta=0.2)
        self.assertStated("×4.6", "low end of the chaperone-anchoring range")
        # scoped: both the abstract and Result 3 must carry the range, not a point
        for name, heading in (("abstract", "## Abstract"),
                              ("Result 3", "### Result 3")):
            sec = flat(section(heading))
            self.assertIn("×4.6", sec,
                          f"{name} reports a point headroom, not the range")
        self.assertIn("roughly one order of magnitude inside", flat())
        self.assertNotIn("roughly two orders of magnitude inside", flat(),
                         "the superseded claim is back")

    def test_result3_inline_table_matches_computed_headroom(self):
        """
        Result 3 quotes a headroom grid inline. every cell must match
        data/computed/headroom_sensitivity.tsv -- an inline table is just as able
        to go stale as a prose number.
        """
        hs = pd.read_csv(COMP / "headroom_sensitivity.tsv", sep="\t")
        ok = hs[~hs.collapsed]
        sec = section("### Result 3")
        rows = [ln for ln in sec.splitlines()
                if ln.startswith("| ") and "×" in ln and "Evaluated at" not in ln]
        self.assertEqual(len(rows), 4,
                         f"expected 4 headroom rows in Result 3, found {len(rows)}")

        # match on the parenthetical, not the exponent: "10⁻⁴" is a substring of
        # "6.3 × 10⁻⁴" and silently matched the wrong row
        key = {"window bottom": "window_bottom",
               "usage-weighted": "usage_weighted_mu",
               "unweighted mean": "unweighted_mean_mu",
               "window top": "window_top"}
        for line in rows:
            cells = [c.strip().replace("*", "") for c in line.strip("|").split("|")]
            label, published, rng = cells[0], cells[1], cells[2]
            ep = next(v for k, v in key.items() if k in label)
            g = ok[ok.evaluation_point == ep]
            pub = g[g.anchoring == "as_published"].iloc[0].headroom_P
            self.assertAlmostEqual(
                float(published.lstrip("×")), pub, delta=0.6,
                msg=f"{ep}: manuscript says {published}, computed ×{pub:.1f}")
            lo_s, hi_s = [x.strip().lstrip("×") for x in rng.split("–")]
            self.assertAlmostEqual(float(lo_s), g.headroom_P.min(), delta=0.6,
                                   msg=f"{ep}: range low disagrees")
            self.assertAlmostEqual(float(hi_s), g.headroom_P.max(), delta=0.6,
                                   msg=f"{ep}: range high disagrees")

    def test_margin_is_near_the_supraadditivity_onset(self):
        h = load_json(COMP / "headroom_sensitivity_summary.json")
        self.assertAlmostEqual(h["internally_consistent_margin_log10"], 1.39,
                               delta=0.02)
        self.assertAlmostEqual(h["distance_from_onset_log10"], 0.20, delta=0.02)
        self.assertIn("1.39 log₁₀", flat())
        self.assertIn("0.20", flat())

    def test_aggregation_death_dominates(self):
        s = load_json(COMP / "bounds_summary.json")
        self.assertAlmostEqual(s["mechanism_frac_aggregation_death"], 0.9994,
                               delta=1e-4)
        self.assertStated("99.94%", "fraction of draws closed by aggregation-death")

    # ---------- Result 4: axis structure ----------
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

    # ---------- Fig 1 parameters ----------
    def test_ode_figure_parameters_are_declared_illustrative(self):
        self.assertStated("ρ = 4.0, χ = 0.15", "Fig 1 parameters")
        self.assertStated("λ_crit = 4.80", "Fig 1 fold location")
        self.assertTrue(
            re.search(r"[Pp]arameters are illustrative and not\s+estimated from data",
                      self.text),
            "the manuscript must state that Fig 1 parameters are not fitted")


class RemovedResultsAreDocumented(unittest.TestCase):
    """the two dropped results must be reported, with reproducible numbers."""

    def setUp(self):
        self.text = manuscript()

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
            self.assertIn(s, self.text, f"within-gene result {s!r} not stated")

    def test_site_position_concordance_is_reported(self):
        d = load_json(COMP / "removed_metal_site_diagnostics.json")
        self.assertAlmostEqual(d["position_concordance_frac"], 0.599, delta=0.01)
        # the codon lookup itself is fine; the position mapping is the bug
        self.assertGreater(d["annotation_codon_matches_cds_codon_frac"], 0.95)
        self.assertIn("~60%", self.text)

    def test_crossspecies_artifacts_are_reported(self):
        s = load_json(COMP / "removed_crossspecies_test.json")
        shared = s["published_delta_A_reproduced_using_ecoli_mu_for_all_species"]
        self.assertTrue(shared["exact"],
                        "shared-mu reproduction of the published table must be exact")
        ind = s["tai_vector_independence"]
        self.assertEqual(ind["ecoli_vs_scerevisiae_bit_identical"], 44)
        self.assertEqual(ind["n_compared"], 60)
        self.assertIn("44 of the 60", self.text)
        self.assertIn("E. coli μ for all three species", self.text)


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


class CitationIntegrity(unittest.TestCase):
    """
    the earlier draft had 20 of 36 references uncited and pointed its
    error-minimization citation at a chaperone protocol. numbering is checked
    mechanically here.
    """

    def setUp(self):
        self.body, self.refs = manuscript().split("## References")

    def cited(self):
        out = set()
        for m in re.finditer(r"\[([0-9,\s\-]+)\]", self.body):
            for part in m.group(1).split(","):
                part = part.strip()
                if part.isdigit():
                    out.add(int(part))
                elif re.fullmatch(r"\d+\s*-\s*\d+", part):
                    a, b = (int(x) for x in part.split("-"))
                    out.update(range(a, b + 1))
        return out

    def listed(self):
        return {int(m.group(1))
                for m in re.finditer(r"^(\d+)\.\s", self.refs, re.M)}

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
    """figures must exist and have a generating script; no orphan images."""

    def test_figures_present(self):
        for stem in ("Fig1_envelope", "Fig2_axis_structure", "Fig3_bounds"):
            for ext in ("png", "pdf"):
                p = FIGS / f"{stem}.{ext}"
                self.assertTrue(p.exists(), f"missing figure {p.name}")

    def test_every_figure_has_a_generating_script(self):
        src = "\n".join(p.read_text() for p in SCRIPTS.glob("*.py"))
        for stem in ("Fig1_envelope", "Fig2_axis_structure", "Fig3_bounds",
                     "Fig4_supraadditivity"):
            self.assertIn(stem, src,
                          f"no script in scripts/ generates {stem}")

    def test_figures_referenced_in_manuscript(self):
        text = manuscript()
        for stem in ("Fig1_envelope", "Fig2_axis_structure", "Fig3_bounds",
                     "Fig4_supraadditivity"):
            self.assertIn(stem, text, f"{stem} has no figure legend")

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
                     "supraadditivity_knob_comparison.tsv"):
            self.assertTrue((COMP / name).exists(), f"missing artifact {name}")

    def test_no_triplet_origin_claim(self):
        """the project's standing constraint: this is not a code-origin paper."""
        text = manuscript()
        self.assertTrue(
            "no code-origin claim" in text or "no claim about the origin" in text,
            "the manuscript must explicitly disclaim a code-origin reading")


class Supraadditivity(unittest.TestCase):
    """
    Result 5. the framework's distinguishing prediction, tested in its own model.
    the point of these assertions is that the result must stay HONEST: the effect
    is real in sign but negligible at wild-type margin, and the paper has to keep
    saying so.
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

    def test_effect_is_negligible_at_wild_type_margin(self):
        pct = self.s["interaction_pct_at_observed_rate"]
        self.assertLess(pct, 1.0, "if this grew, the claim below would change")
        self.assertGreater(pct, 0.0)
        self.assertIn("+0.2% of the additive expectation", self.flat)
        self.assertIn("not a detectable effect", self.flat)

    def test_interaction_grows_as_margin_closes(self):
        self.assertTrue(self.s["interaction_grows_as_margin_closes"])
        t = pd.read_csv(COMP / "supraadditivity_margin_sweep.tsv", sep="\t")
        d = t[t.interaction.notna()].sort_values("margin_baseline", ascending=False)
        self.assertTrue(d.interaction_pct_of_additive.is_monotonic_increasing,
                        "interaction should rise monotonically as margin closes")
        for stated in ("+0.6% at 1.90", "+4.2% at 1.50"):
            self.assertIn(stated, self.flat)

    def test_qualitative_synthetic_lethality(self):
        margin = self.s["qualitative_supraadditivity_first_seen_at_margin_log10"]
        self.assertAlmostEqual(margin, 1.195, delta=0.01)
        self.assertIn("1.19 log₁₀", self.flat)
        self.assertIn("×16 headroom", self.flat)
        self.assertEqual(self.s["n_qualitatively_supraadditive_grid_points"], 338)
        self.assertIn("338 of 676", self.flat)

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

    def test_saturation_limitation_is_disclosed(self):
        sat = self.s["folding_arm_saturation_at_observed_rate"]
        self.assertAlmostEqual(sat, 0.979, delta=0.002)
        self.assertIn("97.9% saturated", self.flat)
        self.assertIn("47.5 µM free", self.flat)
        # and it must be named as a limitation, not just a methods note
        lim = flat(self.text[self.text.index("### Limitations"):])
        self.assertIn("understates how tightly", lim,
                      "the saturation problem must appear in Limitations")

    def test_model_only_caveat_is_stated(self):
        self.assertIn("not whether cells do", self.flat)
        self.assertIn("cannot validate the framework", self.flat)

    def test_predictions_were_rewritten_not_left_stale(self):
        """prediction 1 must no longer promise a wild-type additivity deviation."""
        self.assertIn("synthetic lethality in a sensitized background", self.flat)
        self.assertNotIn("Combining a synonymous change that raises folding "
                         "burden with a reduction in chaperone capacity",
                         self.flat)

    def test_vendored_model_is_unmodified(self):
        a = (SCRIPTS / "vendor" / "two_pool_ode.py").read_bytes()
        b = Path("/storage/kiran-stuff/proteostasis-P1/two_pool_ode.py")
        if b.exists():
            self.assertEqual(a, b.read_bytes(),
                             "the vendored model has drifted from upstream")


class TablesAgreeWithData(unittest.TestCase):
    """
    the generated tables must match data/computed/, and the manuscript must
    match the generated tables. drift between a table and the prose that
    describes it is the failure mode that sank the previous draft.
    """

    EXPECTED = [
        "Table1_burden_terms", "Table2_bounds", "Table2_bounds_statistics",
        "Table3_axis_tests", "Table4_metal_site_backgrounds",
        "TableS1_codon_coordinates", "TableS2_delta_per_aa",
        "TableS3_tai_validation", "TableS4_crossspecies",
        "Table5_supraadditivity", "TableS5_supraadditivity_grid",
        "TableS6_capacity_knob_comparison",
    ]

    def test_all_tables_written(self):
        for name in self.EXPECTED:
            self.assertTrue((TAB / f"{name}.tsv").exists(),
                            f"missing tables/{name}.tsv")
        self.assertTrue((TAB / "TABLES.md").exists(), "missing tables/TABLES.md")

    def test_table2_matches_bounds_summary(self):
        b = load_json(COMP / "bounds_summary.json")
        t = pd.read_csv(TAB / "Table2_bounds.tsv", sep="\t")
        got = dict(zip(t.quantity, t.value_per_codon))
        self.assertAlmostEqual(
            got["Arithmetic bound, paired MC median"],
            b["arithmetic_paired_median"], places=12)
        self.assertAlmostEqual(
            got["Two-pool ODE bound, paired MC median"],
            b["ode_paired_median"], places=12)

        s = pd.read_csv(TAB / "Table2_bounds_statistics.tsv", sep="\t")
        stats = dict(zip(s.statistic, s.value))
        self.assertAlmostEqual(stats["P(arithmetic is the tighter bound)"],
                               b["paired_P_arith_tighter"], places=12)

    def test_table3_matches_axis_tests(self):
        t = pd.read_csv(TAB / "Table3_axis_tests.tsv", sep="\t")
        src = pd.read_csv(COMP / "axis_tests.tsv", sep="\t")
        for _, r in src.iterrows():
            m = t[(t.axis == r.axis) & (t.null == r.null)
                  & (t.mu_stat.isin(["mean", "both"]))]
            self.assertEqual(len(m), 1,
                             f"Table 3 has {len(m)} rows for {r.axis}/{r.null}")
            self.assertAlmostEqual(float(m.iloc[0].z), float(r.z), places=10)

    def test_table3_does_not_duplicate_the_nu_rows(self):
        """
        nu does not depend on the mu summary statistic, so listing it twice would
        present one test as two.
        """
        t = pd.read_csv(TAB / "Table3_axis_tests.tsv", sep="\t")
        nu = t[t.axis == "nu"]
        self.assertEqual(len(nu), 2, "nu should appear once per null model")
        self.assertTrue((nu.mu_stat == "both").all(),
                        "nu rows must not be attributed to a mu statistic")

    def test_table3_labels_null_results_as_null(self):
        t = pd.read_csv(TAB / "Table3_axis_tests.tsv", sep="\t")
        ns = t[t.p_one_sided > 0.05]
        self.assertGreater(len(ns), 0, "expected at least one null row")
        self.assertTrue((ns.direction == "no signal").all(),
                        "non-significant rows must not be labelled "
                        "clustered/spread on the strength of their sign")

    def test_table4_matches_removed_metal_site_test(self):
        t = pd.read_csv(TAB / "Table4_metal_site_backgrounds.tsv", sep="\t")
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
                       "1.00 × 10⁻²", "z = −3.56", "0.768"):
            bare = needle.replace("z = ", "")
            self.assertIn(bare, md,
                          f"{bare!r} is in the manuscript but not in TABLES.md")
            self.assertIn(needle, text,
                          f"{needle!r} is in TABLES.md but not in the manuscript")

    def test_manuscript_cites_the_tables(self):
        text = manuscript()
        for label in ("Table 1", "Table 2", "Table 3", "Table 4", "Table 5",
                      "Table S1", "Table S2", "Table S3", "Table S4",
                      "Table S5", "Table S6"):
            self.assertIn(label, text, f"manuscript never refers to {label}")

    def test_removed_result_tables_are_marked_as_removed(self):
        """a reader must not mistake Table 4 or S4 for a finding of this paper."""
        md = (TAB / "TABLES.md").read_text()
        for heading in ("## Table 4.", "## Table S4."):
            i = md.index(heading)
            block = md[i:i + 700]
            self.assertIn("excluded from the paper", block,
                          f"{heading} does not say the result is excluded")


if __name__ == "__main__":
    unittest.main(verbosity=2)
