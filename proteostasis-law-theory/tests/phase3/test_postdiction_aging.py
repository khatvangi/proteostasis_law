"""phase 3: the aging/rejuvenation post-diction, which FAILED.

these tests pin a negative result so it cannot drift into a positive one. the
tempting error is to quote the constant-dilution regime -- which is bistable and
therefore looks like a success -- while ignoring that the same regime predicts
zero reproductive loss, contradicting the very measurement being explained.
"""

import sys
import unittest
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (_REPO_ROOT / "scripts", _REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import dilution as D  # noqa: E402
import boundary_structure as BS  # noqa: E402
import postdiction_aging as PA  # noqa: E402


class TestObservationProvenance(unittest.TestCase):
    def testCitationIsRecordedWithTheNumber(self):
        self.assertEqual(PA.LINDNER2008["doi"], "10.1073/pnas.0708931105")
        self.assertEqual(PA.LINDNER2008["pmid"], "18287048")
        self.assertEqual(PA.LINDNER2008["pmcid"], "PMC2268587")
        # the measured quantity is DERIVED from two source numbers, not quoted
        self.assertAlmostEqual(PA.LINDNER2008["aging_effect"], 0.0395)
        self.assertEqual(PA.LINDNER2008["aggregate_share"], (0.30, 0.40))
        lo, hi = PA.aggregateAttributableLoss()
        self.assertAlmostEqual(lo, 0.01035, places=6)
        self.assertAlmostEqual(hi, 0.01780, places=6)

    def testTheAbstractsThirtyPercentIsAShareNotTheEffect(self):
        """reading '>30% OF the loss' as '>30% loss' was wrong by ~30x."""
        _, hi = PA.aggregateAttributableLoss()
        self.assertLess(hi, 0.02)


class TestConstantDilutionIsBistableButPredictsNoGrowthCost(unittest.TestCase):
    """the trap: bistable, so it looks like the post-diction works."""

    def testItIsBistable(self):
        df = PA.scanInfluxes(M.Params().validate(), D.Growth(0.04))
        self.assertFalse(df.empty)
        self.assertGreater(df["aggregate_ratio"].min(), 5.0)

    def testAndItPredictsExactlyZeroReproductiveLoss(self):
        """k_mu = inf means growth cannot respond to burden, by construction."""
        g = D.Growth(0.04)
        for u, a in ((0.05, 0.1), (0.1, 4.0)):
            self.assertEqual(g.rate(u, a), g.mu0)

    def testRejuvenationThresholdIsAFractionNotAllOrNothing(self):
        df = PA.scanInfluxes(M.Params().validate(), D.Growth(0.04))
        self.assertTrue((df["shed_fraction_needed"] > 0.0).all())
        self.assertTrue((df["shed_fraction_needed"] < 1.0).all())


class TestUnderFeedbackTheModelFails(unittest.TestCase):
    """the honest result: bistability is not generic once growth responds."""

    @classmethod
    def setUpClass(cls):
        cls.fb = PA.feedbackScan(M.Params().validate())

    def testMostSettingsAreMonostable(self):
        self.assertFalse(self.fb.empty)
        self.assertLess(self.fb["bistable"].mean(), 0.75,
                        "bistability must not be reported as generic under feedback")

    def testWhereBistableTheGrowthCostExceedsTheObserved(self):
        bi = self.fb[self.fb["bistable"]]
        if bi.empty:
            self.skipTest("no bistable feedback setting converged")
        self.assertGreater(bi["reproductive_loss"].median(),
                           PA.aggregateAttributableLoss()[1],
                           "predicted loss should exceed the measured band")


class TestLinearArrestGivesNoRejuvenationAtAll(unittest.TestCase):
    def testNoBistableWindow(self):
        df = PA.scanInfluxes(M.Params().validate(), BS.LinearGrowth(0.3, 1.5625))
        self.assertTrue(df.empty)


class TestVerdictIsRecordedAsNegative(unittest.TestCase):
    def testModuleStatesTheFailureAndThePointer(self):
        src = (_REPO_ROOT / "scripts" / "phase3" / "postdiction_aging.py").read_text()
        self.assertIn("does not post-dict", src.lower())
        self.assertIn("sequestration", src.lower())


if __name__ == "__main__":
    unittest.main()
