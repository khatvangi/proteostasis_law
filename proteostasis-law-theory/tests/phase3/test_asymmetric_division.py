"""phase 3: the framing test -- was bistability ever required?

the result is a PASS, and passes are the outcome with no track record in this
project, so these pin the things that make it a pass rather than an artefact:
that the f = 0.5 control is EXACTLY zero rather than merely small, that the band
is imported and not re-derived, that the effect is monotone in f, and that the
model carries no second attractor at all.
"""

import sys
import unittest
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (_REPO_ROOT / "scripts", _REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import dilution as D  # noqa: E402
import postdiction_aging as PA  # noqa: E402
import asymmetric_division as A  # noqa: E402


class TestTheBandIsImportedNotRederived(unittest.TestCase):
    def testEdgesComeFromTheFixedSource(self):
        self.assertEqual((A.BAND_LO, A.BAND_HI), PA.aggregateAttributableLoss())

    def testControlFractionIsAHalf(self):
        self.assertEqual(A.CONTROL_F, 0.5)


class TestTheControlIsExactlyZeroNotMerelySmall(unittest.TestCase):
    """symmetric partitioning in half the volume leaves concentration unchanged,
    which is what the -mu.x dilution term already encodes."""

    def testSymmetricSplitIsTheIdentityOnConcentration(self):
        for a in (0.0, 1e-3, 0.25, 7.5):
            self.assertEqual(2.0 * A.CONTROL_F * a, a)

    def testAgingEffectVanishesAtTheControl(self):
        p0 = M.Params().validate()
        g = A.calibratedHyperbolic(0.1, 0.02)
        fold = D.foldSolveDil(p0, g)
        self.assertIsNotNone(fold)
        p = p0.with_(j=0.5 * fold[0]).validate()
        r = A.lineageContrast(p, g, A.CONTROL_F)
        self.assertIsNotNone(r)
        self.assertEqual(r["aging_effect"], 0.0)
        self.assertEqual(r["gr_old"], r["gr_new"])


class TestTheEffectRisesWithAsymmetry(unittest.TestCase):
    def testMonotoneInF(self):
        p0 = M.Params().validate()
        g = A.calibratedHyperbolic(0.1, 0.02)
        fold = D.foldSolveDil(p0, g)
        p = p0.with_(j=0.5 * fold[0]).validate()
        vals = []
        for f in (0.5, 0.6, 0.7, 0.8, 0.9):
            r = A.lineageContrast(p, g, f)
            self.assertIsNotNone(r)
            vals.append(r["aging_effect"])
        self.assertTrue(all(b > a for a, b in zip(vals, vals[1:])),
                        f"aging effect must rise with f, got {vals}")

    def testTheOldPoleDaughterCarriesMoreAggregateAndGrowsSlower(self):
        p0 = M.Params().validate()
        g = A.calibratedHyperbolic(0.1, 0.02)
        fold = D.foldSolveDil(p0, g)
        p = p0.with_(j=0.5 * fold[0]).validate()
        r = A.lineageContrast(p, g, 0.8)
        self.assertGreater(r["a_old_start"], r["a_new_start"])
        self.assertLess(r["gr_old"], r["gr_new"])


class TestNoSecondAttractorIsUsed(unittest.TestCase):
    """the whole point: this is the MONOSTABLE model, with no mechanism added."""

    def testTheGenerationMapUsesTheUnmodifiedDilutedField(self):
        p = M.Params().validate()
        g = A.calibratedHyperbolic(0.1, 0.02)
        out = A.oneGeneration((1e-3, 1e-5), p, g)
        self.assertIsNotNone(out)
        # interdivision time is an OUTPUT, defined by volume doubling
        self.assertAlmostEqual(out["growth_rate"] * out["T"], np.log(2.0), places=9)

    def testHeavierBurdenLengthensTheGeneration(self):
        p = M.Params().validate()
        g = A.calibratedHyperbolic(0.1, 0.02)
        light = A.oneGeneration((1e-3, 1e-5), p, g)
        heavy = A.oneGeneration((1e-3, 0.5), p, g)
        self.assertGreater(heavy["T"], light["T"])


class TestVerdictRespectsRuleSix(unittest.TestCase):
    def testControlCannotCarryThePass(self):
        import pandas as pd
        df = pd.DataFrame([
            {"f": 0.5, "aging_effect": 0.012, "settled": True},
            {"f": 0.9, "aging_effect": 0.400, "settled": True},
        ])
        v = A.verdict(df)
        self.assertEqual(v["in_band_control"], 1)
        self.assertEqual(v["in_band_mechanism_on"], 0)
        self.assertFalse(v["mechanism_passes"])


class TestTheWithdrawalIsRecorded(unittest.TestCase):
    def testD031WithdrawsD026sSurvivingClaim(self):
        doc = (_REPO_ROOT / "DECISIONS.md").read_text()
        i = doc.index("## D031 —")
        entry = doc[i:]
        self.assertIn("WITHDRAWN", entry)
        self.assertIn("continuously renewed perturbation, not an attractor", entry)

    def testTheAuditFlagIsSettledAsLawSaturationNotSampling(self):
        doc = (_REPO_ROOT / "DECISIONS.md").read_text()
        entry = doc[doc.index("## D031 —"):]
        self.assertIn("independent", entry)
        self.assertIn("clamp, not a measurement", entry)

    def testLinearArrestReallyDoesReturnExactlyZeroPastKmu(self):
        """the mechanism behind the identical 1.000 values."""
        import boundary_structure as B
        g = B.LinearGrowth(0.05, 0.5)
        self.assertEqual(g.rate(0.0, 17.03), 0.0)
        self.assertEqual(g.rate(0.0, 1.71), 0.0)
        self.assertGreater(g.rate(0.0, 0.1), 0.0)


if __name__ == "__main__":
    unittest.main()
