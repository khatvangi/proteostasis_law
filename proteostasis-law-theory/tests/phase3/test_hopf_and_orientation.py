"""phase 3: the manuscript's Block 1 numbers, asserted against their generators.

Every number Section 3, Remark 3, Corollary 1 and Section 3.2 assert is
recomputed here from `data/computed/`, never from a serialised summary of a
summary. The rule these enforce is D046's: solving `det J = 0` returns candidates
without distinguishing them, so orientation and branch identity must be COMPUTED
rather than inferred from a solver having converged.

The counts are asserted, not printed. A count nobody asserts against is a
comment -- and three of the four Block 1 artefacts were caught by exactly the
counts a tidier report would have omitted.
"""

import re
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parents[2]
_COMPUTED = _REPO_ROOT / "data" / "computed"
_PAPER = _REPO_ROOT / "manuscript" / "bmb_v4.md"


def _load(name):
    p = _COMPUTED / name
    if not p.exists():
        raise unittest.SkipTest(f"{name} not built; run the phase 3 generators")
    return pd.read_csv(p, sep="\t")


class TestGenericityConditions(unittest.TestCase):
    """Theorem 1's converse holds where the paper says it holds."""

    def setUp(self):
        self.text = _PAPER.read_text()

    def testLoadGridMarginsMatchTheManuscript(self):
        D = _load("genericity_load_grid.tsv")
        ok = D[D["ok"] == True]  # noqa: E712
        self.assertEqual(len(ok), 325)
        for col, quoted in (("grad_G", 0.106), ("tr_J", 0.303),
                            ("d2R_signed", 0.0929), ("transversality", 0.341)):
            got = float(ok[col].dropna().abs().min())
            self.assertAlmostEqual(got, quoted, delta=0.0006,
                                   msg=f"{col} margin moved: {got:.4g}")
            self.assertIn(f"{quoted:g}".rstrip("0").rstrip("."), self.text)

    def testKineticBoxHasExactlyTwoViolationsAndTheyAreAtTheDomainEdge(self):
        D = _load("genericity_kinetic_box.tsv")
        ok = D[D["ok"] == True]  # noqa: E712
        bad = ok[(ok["grad_G"].abs() <= 1e-6) | (ok["transversality"].abs() <= 1e-6)]
        self.assertEqual(len(bad), 2)
        # both sit at the edge of the domain, which is why they are reported
        # rather than treated as counterexamples
        self.assertLess(float(bad["a"].max()), 1e-10)
        self.assertIn("2765 of 2767", self.text)
        self.assertIn("117 of 2884", self.text)

    def testTransversalityIsTheConditionWithNoStructuralGuarantee(self):
        """(G4) reduces to w_1 != 0 because dF/dj = (1,0) exactly."""
        seg = self.text[self.text.index("(G4)"):][:1200]
        self.assertIn("(1, 0)", seg + self.text[:0] or seg)
        self.assertIn("no structural guarantee", self.text)


class TestFoldOrientation(unittest.TestCase):
    """Remark 3 is a classification: folds come in two orientations."""

    def setUp(self):
        self.text = _PAPER.read_text()

    def testTheLoadGridIsEntirelyCollapseOriented(self):
        D = _load("genericity_load_grid.tsv")
        d2 = D[D["ok"] == True]["d2R_signed"].dropna()  # noqa: E712
        self.assertEqual(int((d2 < 0).sum()), 325)
        self.assertEqual(int((d2 > 0).sum()), 0)
        self.assertIn("325 of 325 load-grid folds", self.text)

    def testTwentySixKineticBoxFoldsAreBirthOriented(self):
        D = _load("genericity_kinetic_box.tsv")
        d2 = D[D["ok"] == True]["d2R_signed"].dropna()  # noqa: E712
        self.assertEqual(int((d2 > 0).sum()), 26)
        self.assertEqual(int((d2 < 0).sum()), 2739)
        self.assertIn("2739 of 2765", self.text)

    def testEveryBirthOrientedNetworkAlsoCarriesACollapsePointAboveIt(self):
        """the claim that makes them hysteresis loops rather than counterexamples."""
        R = _load("fold_orientation.tsv")
        pos = R[(R["group"] == "d2R>0") & (R["ok"] == True)]  # noqa: E712
        self.assertEqual(len(pos), 26)
        self.assertEqual(int(pos["has_both"].sum()), 26,
                         "a birth-oriented fold with no collapse point above it "
                         "would be a counterexample, not a hysteresis loop")
        self.assertEqual(int((pos["j_at_min"] < pos["j_at_max"]).sum()), 26)
        self.assertIn("all 26 of those networks carry a collapse-oriented", self.text)

    def testTheTwoOrientationRoutesAgree(self):
        """arclength d2R/ds2 and the shape of j along the branch are independent."""
        R = _load("fold_orientation.tsv")
        pos = R[(R["group"] == "d2R>0") & (R["ok"] == True)]  # noqa: E712
        ctl = R[(R["group"] == "d2R<0") & (R["ok"] == True)]  # noqa: E712
        self.assertEqual((int(pos["orient_agree"].sum()),
                          int(pos["orient_n"].sum())), (54, 54))
        self.assertEqual((int(ctl["orient_agree"].sum()),
                          int(ctl["orient_n"].sum())), (162, 167))
        self.assertEqual(int(ctl["has_both"].sum()), 7)


class TestMultiplicityScopesCorollaryOne(unittest.TestCase):
    """Corollary 1 claims a saving with a named residual, not an unqualified one."""

    def setUp(self):
        self.text = _PAPER.read_text()

    def testNinePercentOfNetworksCarryMoreThanOneCandidate(self):
        B = _load("branch_kinetic_box.tsv")
        B = B[B["traced"] == True]  # noqa: E712
        n_multi = int((B["n_singular"] > 1).sum())
        self.assertEqual((n_multi, len(B)), (252, 2765))
        self.assertAlmostEqual(100.0 * n_multi / len(B), 9.1, delta=0.05)
        self.assertEqual(int(B["n_singular"].max()), 5)
        self.assertIn("252 of 2765", self.text)
        self.assertIn("up to five", self.text)

    def testCorollaryOneSaysCandidatesNotTheFold(self):
        seg = self.text[self.text.index("**Corollary 1"):][:1400]
        self.assertIn("candidates", seg)
        self.assertIn("9.1%", seg)


class TestHopfPrecedesTheFold(unittest.TestCase):
    """Section 3.2. The result, and the control that gives it meaning."""

    def setUp(self):
        self.text = _PAPER.read_text()
        self.S = _load("hopf_refined_kinetic_box.tsv")
        self.L = _load("hopf_refined_load_grid.tsv")
        self.I = _load("hopf_integration_kinetic_box.tsv")
        T = self.S[self.S["traced"] == True]  # noqa: E712
        self.clean = T[(T["tr_max"] >= 0.0) & (T["fold_is_j_max"] == 1)]
        self.traced = T

    def testLoadGridNeverCrossesAndTheManuscriptSaysWhatThatProves(self):
        L = self.L[self.L["traced"] == True]  # noqa: E712
        self.assertEqual(len(L), 325)
        self.assertEqual(int((L["tr_max"] >= 0.0).sum()), 0)
        self.assertAlmostEqual(float(L["tr_max"].max()), -0.243, delta=0.0006)
        # the load grid is ONE POINT in kinetic space, and the manuscript must
        # not let its cleanliness read as a bound on the oscillatory corner
        seg = self.text[self.text.index("### 3.2"):]
        seg = seg[:seg.index("## 4.")]
        self.assertIn("not that the corner is small", seg)
        self.assertIn("holds kinetics fixed and sweeps load", seg)

    def testOneHundredAndFourNetworksCrossBeforeTheFold(self):
        self.assertEqual((len(self.clean), len(self.traced)), (104, 2766))
        self.assertIn("104 of 2766", self.text)
        frac = self.clean["j_at_first_cross"] / self.clean["j_crit"]
        self.assertAlmostEqual(float(frac.median()), 0.83, delta=0.006)
        self.assertAlmostEqual(float(frac.min()), 0.12, delta=0.006)

    def testEveryCrossingHasPositiveDeterminantSoItIsAHopfByDefinition(self):
        det = self.clean["det_at_first_cross"]
        self.assertEqual(int((det <= 0).sum()), 0)
        self.assertGreaterEqual(float(det.min()), 1.4e-6)
        self.assertIn("Hopf points by definition", self.text)

    def testTheBranchPointsAreExactEquilibria(self):
        """reprojection is what makes tr J at these points mean anything."""
        self.assertLess(float(self.traced["G_absmax"].max()), 1e-13)
        self.assertLess(float(self.L[self.L["traced"] == True]["G_absmax"].max()),  # noqa: E712
                        1e-16)

    def testTheControlBlockDoesNotMoveAtAll(self):
        """without this the integration test is unfalsifiable."""
        ok = self.I[self.I["tested"] == True]  # noqa: E712
        cross = ok[ok["name"].str.startswith("+")]
        ctrl = ok[ok["name"].str.startswith("-")]
        self.assertEqual((len(cross), len(ctrl)), (104, 205))
        self.assertEqual(int(ctrl["escaped"].sum()), 0)
        self.assertEqual(int(ctrl["grew"].sum()), 0)
        self.assertAlmostEqual(float(ctrl["ratio_max"].median()), 1.0, delta=0.01)
        # and the two distributions must not overlap
        self.assertGreater(float(cross["ratio_max"].min()),
                           float(ctrl["ratio_max"].max()))
        self.assertEqual(int(cross["grew"].sum()), 104)
        self.assertEqual(int(cross["escaped"].sum()), 102)

    def testIntegrationRecoversTheEigenvalueItWasNotGiven(self):
        ok = self.I[self.I["tested"] == True]  # noqa: E712
        cr = ok[ok["name"].str.startswith("+")]
        fit = cr[cr["slope"].notna()]
        rel = (fit["slope"] - fit["lambda_max"]).abs() / fit["lambda_max"].abs()
        self.assertEqual(int((rel < 0.05).sum()), 90)
        per = cr[cr["period_measured"].notna()]
        pr = (per["period_measured"] - per["period_predicted"]).abs() \
            / per["period_predicted"]
        self.assertEqual((int((pr < 0.05).sum()), len(per)), (47, 49))
        self.assertLess(float(pr.median()), 1.1e-3)

    def testThePredictedPeriodIsPiOverOmegaNotTwoPi(self):
        """|delta|^2 oscillates at 2*omega; the first version was off by 2x."""
        ok = self.I[self.I["tested"] == True]  # noqa: E712
        per = ok[ok["period_measured"].notna() & ok["omega"].notna()]
        pred = np.pi / per["omega"]
        self.assertTrue(np.allclose(per["period_predicted"], pred, rtol=1e-9))
        self.assertIn("`π/ω`", self.text)

    def testTheParameterCornerIsWhereTheManuscriptSaysItIs(self):
        g = pd.read_csv(_COMPUTED / "hopf_parameter_corner.tsv", sep="\t",
                        index_col=0)
        self.assertAlmostEqual(1.0 / float(g.loc["p_kappa_a", "ratio"]), 11.7,
                               delta=0.1)
        self.assertAlmostEqual(float(g.loc["p_rho_A", "ratio"]), 7.5, delta=0.05)
        self.assertAlmostEqual(1.0 / float(g.loc["p_alpha_n", "ratio"]), 5.5,
                               delta=0.05)
        for quoted in ("11.7 times lower", "7.5 times higher",
                       "5.5 times lower", "3.7 to 33.6"):
            self.assertIn(quoted, self.text)

    def testTheIncidenceRateIsNotOfferedAsAPrediction(self):
        seg = self.text[self.text.index("### 3.2"):]
        seg = seg[:seg.index("## 4.")]
        self.assertIn("3.76%", seg)
        self.assertIn("property of the sampling box", seg)
        self.assertIn("do not offer 3.76% as a prediction", seg)


class TestTitleAndTerminology(unittest.TestCase):
    """the paper says fold, not collapse threshold (D046)."""

    def setUp(self):
        self.text = _PAPER.read_text()

    def testTheTitleIsTheFoldCondition(self):
        first = self.text.splitlines()[0]
        self.assertEqual(
            first,
            "# An Exact Fold Condition for Mass-Balanced Models of "
            "Protein Quality Control")

    def testCollapseThresholdIsNotUsedAsALABELForTheFold(self):
        """pin the PROPERTY, not the token.

        A first version of this test banned the string outright and failed on two
        sentences that use the term correctly -- Remark 3 saying a birth-oriented
        fold is NOT a collapse threshold, and the discussion saying only one
        orientation IS one. Both need the term to say what they say. What must go
        is the term used as a NAME for the fold, which is where it appears in
        titles, keywords, headings and the lead of a figure caption.
        """
        lines = self.text.splitlines()
        naming = [ln for ln in lines
                  if (ln.startswith("#") or ln.startswith("**Keywords")
                      or re.match(r"\*\*Fig\. \d+\*\*", ln))
                  and ("collapse threshold" in ln.lower()
                       or "collapse boundary" in ln.lower())]
        self.assertEqual(naming, [], f"the fold is being NAMED a collapse "
                                     f"threshold in: {naming}")
        # and the phenomenon is still discussed, so the word itself must survive
        self.assertIn("proteostasis collapse", self.text.lower())
        # the two legitimate uses distinguish rather than name; keep them honest
        self.assertIn("not a collapse threshold at all", self.text)

    def testTheMetadataTitlesTrackTheManuscript(self):
        title = ("An Exact Fold Condition for Mass-Balanced Models of "
                 "Protein Quality Control")
        for f in ("CITATION.cff", ".zenodo.json"):
            self.assertIn(title, (_REPO_ROOT / f).read_text(), f"{f} title stale")


class TestThreeDegeneracies(unittest.TestCase):
    def testSectionThreeOneNamesAllThree(self):
        text = _PAPER.read_text()
        seg = text[text.index("**Neither numerator."):][:1400]
        for token in ("`det H = 0`", "`det J = 0`", "`tr J = 0`", "Hopf"):
            self.assertIn(token, seg)


if __name__ == "__main__":
    unittest.main()
