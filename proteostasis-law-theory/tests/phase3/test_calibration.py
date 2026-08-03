"""phase 3 calibration against a measured growth-burden relation.

the anchor is Geiler-Samerotte et al. 2011 PNAS, doi:10.1073/pnas.1017570108
(PMID 21187411), retrieved via PubMed: a 3.2 % growth-rate reduction at less than
0.1 % of total cellular protein misfolded, in yeast.

these tests pin the conversion arithmetic and the two substantive outcomes: a
boundary survives everywhere in the calibrated range, and bistability does NOT --
it depends on the high-burden shape of the growth law, which the measurement does
not constrain.
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
import boundary_structure as B  # noqa: E402
import calibration as C  # noqa: E402


class TestAnchorArithmetic(unittest.TestCase):
    def testSlopeIsDerivedFromTheRecordedMeasurement(self):
        self.assertAlmostEqual(
            C.SLOPE_PER_PROTEOME_FRACTION,
            C.GS2011["growth_rate_reduction"] / C.GS2011["misfolded_proteome_fraction"],
            places=12)
        self.assertAlmostEqual(C.SLOPE_PER_PROTEOME_FRACTION, 32.0, places=9)

    def testArrestFractionIsTheSlopeReciprocal(self):
        self.assertAlmostEqual(C.ARREST_PROTEOME_FRACTION,
                               1.0 / C.SLOPE_PER_PROTEOME_FRACTION, places=15)
        self.assertAlmostEqual(C.ARREST_PROTEOME_FRACTION, 0.03125, places=9)

    def testProvenanceIsRecordedWithTheNumber(self):
        """a measured constant must carry its source in the artefact itself."""
        self.assertEqual(C.GS2011["doi"], "10.1073/pnas.1017570108")
        self.assertEqual(C.GS2011["pmid"], "21187411")
        self.assertIn("cerevisiae", C.GS2011["organism"])

    def testConversionIsMonotoneAndRejectsNonFractions(self):
        ks = [C.kMuFromProteomeShare(x) for x in (0.005, 0.01, 0.02, 0.05, 0.1)]
        self.assertTrue(all(a > b for a, b in zip(ks, ks[1:])), ks)
        for bad in (0.0, 1.0, -0.1, 2.0):
            with self.assertRaises(ValueError):
                C.kMuFromProteomeShare(bad)


class TestCalibratedRangeBracketsThePriorGuesses(unittest.TestCase):
    """the values chosen before any data was consulted were 0.5 and 2.0."""

    def testGuessesLieInsideTheCalibratedRange(self):
        lo = C.kMuFromProteomeShare(0.10)
        hi = C.kMuFromProteomeShare(0.005)
        for guess in (0.5, 2.0):
            self.assertGreater(guess, lo)
            self.assertLess(guess, hi)


class TestBoundarySurvivesCalibration(unittest.TestCase):
    def testBoundaryExistsAcrossTheCalibratedRange(self):
        S = C.sweepCalibration(p_qcs=(0.01, 0.02, 0.05), mu0s=(0.01, 0.1, 0.3))
        self.assertEqual(len(S), 9)
        self.assertTrue(S["boundary"].all(),
                        "a calibrated setting lost the collapse boundary")

    def testEnzymaticConditionStaysInANarrowBand(self):
        S = C.sweepCalibration(p_qcs=(0.01, 0.02, 0.05), mu0s=(0.01, 0.1, 0.3))
        ok = S[S["boundary"]]
        self.assertGreater(ok["phi_enz"].min(), 0.02)
        self.assertLess(ok["phi_enz"].max(), 0.30)


class TestBistabilityIsFormDependentNotMeasured(unittest.TestCase):
    """the substantive negative result from calibration."""

    def testHyperbolicFormGivesAGenuineSecondAttractor(self):
        p = M.Params().validate()
        h = B.hysteresisSweep(p, D.Growth(0.04),
                              [0.10, 0.155, 0.17, 0.19, 0.196, 0.21])
        self.assertTrue(h["all_settled"], "the D013 case must settle on both branches")
        self.assertIsNotNone(h["window"])

    def testMeasuredLinearFormDoesNot(self):
        r = C.bistabilityUnderCalibration(p_qc=0.02, mu0=0.05, n_j=6)
        self.assertTrue(r["boundary"])
        self.assertEqual(r["down_settled"], 0,
                         "the linear law should give a runaway, not an attractor")
        self.assertIsNone(r["window"])

    def testTheReasonIsCompleteVersusAsymptoticArrest(self):
        """linear arrest switches dilution off entirely; hyperbolic never does."""
        big = 10.0
        self.assertEqual(B.LinearGrowth(mu0=0.3, k_mu=1.5625).rate(big, big), 0.0)
        self.assertGreater(D.Growth(mu0=0.3, k_mu=1.5625).rate(big, big), 0.0)


if __name__ == "__main__":
    unittest.main()
