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


def manuscript():
    return MS.read_text()


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

    def test_headroom(self):
        s = load_json(COMP / "bounds_summary.json")
        self.assertAlmostEqual(s["headroom_P"], 158.06, delta=0.1)
        self.assertAlmostEqual(s["headroom_A"], 11091, delta=1)
        self.assertStated("×158", "misfolded-monomer headroom")
        self.assertStated("×1.1 × 10⁴", "aggregated-pool headroom")

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
        for stem in ("Fig1_envelope", "Fig2_axis_structure", "Fig3_bounds"):
            self.assertIn(stem, src,
                          f"no script in scripts/ generates {stem}")

    def test_figures_referenced_in_manuscript(self):
        text = manuscript()
        for stem in ("Fig1_envelope", "Fig2_axis_structure", "Fig3_bounds"):
            self.assertIn(stem, text, f"{stem} has no figure legend")

    def test_computed_outputs_present(self):
        for name in ("nu_tai_ecoli_validated.tsv", "tai_validation_report.json",
                     "axis_tests.tsv", "axis_tests_median.tsv",
                     "mu_variance_decomposition.json", "translation_burden.json",
                     "bounds_summary.json", "removed_metal_site_test.tsv",
                     "removed_crossspecies_test.json"):
            self.assertTrue((COMP / name).exists(), f"missing artifact {name}")

    def test_no_triplet_origin_claim(self):
        """the project's standing constraint: this is not a code-origin paper."""
        text = manuscript()
        self.assertTrue(
            "no code-origin claim" in text or "no claim about the origin" in text,
            "the manuscript must explicitly disclaim a code-origin reading")


if __name__ == "__main__":
    unittest.main(verbosity=2)
