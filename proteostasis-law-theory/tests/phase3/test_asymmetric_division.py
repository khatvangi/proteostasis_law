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


class TestTheParameterFreeForm(unittest.TestCase):
    """D032: the match is an accommodation until the knobs are removed."""

    def testPqcCancelsExactly(self):
        """dB/k_mu == 32 * dB_as_proteome_fraction, identically."""
        import calibration as C
        for p_qc in (0.005, 0.02, 0.10):
            k_mu = C.kMuFromProteomeShare(p_qc)
            for dB in (1e-4, 1e-2, 0.5):
                self.assertAlmostEqual(
                    dB / k_mu, A.SLOPE_PER_PROTEOME_FRACTION * dB * p_qc,
                    places=12)

    def testFIsPinnedAtOneByTheIndivisibleFocus(self):
        self.assertEqual(A.MEASURED_F, 1.0)

    def testTheRequirementIsSubTenthOfAPercentOfProteome(self):
        lo, hi = A.requiredAggregateFraction(damping=0.4386)
        self.assertLess(hi, 1e-3)
        self.assertGreater(lo, 1e-4)

    def testTheRequirementSitsFarBelowTheOnlyMeasuredAggregateFraction(self):
        """Tomoyasu 2001: rpoH-null 5-10%; wild type UNDETECTED, so a bound."""
        self.assertIsNone(A.TOMOYASU2001["wild_type_30C"])
        lo, hi = A.requiredAggregateFraction(damping=0.4386)
        rpoh_lo, _ = A.TOMOYASU2001["rpoH_null_30C"]
        self.assertGreater(rpoh_lo / hi, 50.0)

    def testDecisionEntryCallsTheMatchAnAccommodation(self):
        doc = (_REPO_ROOT / "DECISIONS.md").read_text()
        entry = doc[doc.index("## D032 —"):]
        self.assertIn("accommodation", entry)
        self.assertIn("FOUR", entry)
        # and keeps the structural conclusion, which needs no fit
        self.assertIn("does not depend on any of", entry)

    def testD012CarriesTheWithdrawnEmpiricalMotivation(self):
        doc = (_REPO_ROOT / "DECISIONS.md").read_text()
        i = doc.index("## D012 —")
        entry = doc[i:doc.index("## D013 —")]
        self.assertIn("Empirical motivation withdrawn", entry)

    def testD029SeverityIsMarkedAsALowerBound(self):
        doc = (_REPO_ROOT / "DECISIONS.md").read_text()
        i = doc.index("## D029 —")
        entry = doc[i:doc.index("## D030 —")]
        self.assertIn("UNDERSTATE", entry)
        self.assertIn("lower bound on the miss", entry)


if __name__ == "__main__":
    unittest.main()


class TestBetaTwoCompartmentPartition(unittest.TestCase):
    """D033: the focus is not the whole of `a`, so f is determined, not pinned."""

    def testBetaOneReducesToD031ExactlyAtZero(self):
        self.assertEqual(A.fEffFromBeta(1.0), A.MEASURED_F)
        self.assertEqual(A.partitionMultipliers(1.0), (2.0, 0.0))

    def testBetaZeroReducesToTheControlExactlyAtZero(self):
        self.assertEqual(A.fEffFromBeta(0.0), A.CONTROL_F)
        self.assertEqual(A.partitionMultipliers(0.0), (1.0, 1.0))

    def testTheTwoCompartmentRuleIsTheScalarRuleReparameterised(self):
        for b in (0.0, 0.145, 0.25, 0.5, 0.75, 1.0):
            old, new = A.partitionMultipliers(b)
            f = A.fEffFromBeta(b)
            self.assertAlmostEqual(old, 2.0 * f, places=12)
            self.assertAlmostEqual(new, 2.0 * (1.0 - f), places=12)

    def testTheControlStillGivesExactlyZeroUnderTheBetaForm(self):
        p0 = M.Params().validate()
        g = A.calibratedHyperbolic(0.1, 0.02)
        fold = D.foldSolveDil(p0, g)
        p = p0.with_(j=0.5 * fold[0]).validate()
        r = A.lineageContrast(p, g, A.fEffFromBeta(0.0))
        self.assertEqual(r["aging_effect"], 0.0)

    def testRequirementIsMonotoneDecreasingInBeta(self):
        """lower focus share needs MORE aggregate. this is the reportable direction."""
        prev = None
        for b in (1.0, 0.75, 0.5, 0.25, 0.145):
            lo, hi = A.requiredAggregateFractionBeta(b, 0.35)
            if prev is not None:
                self.assertGreater(lo, prev[0])
                self.assertGreater(hi, prev[1])
            prev = (lo, hi)

    def testRequirementScalesAsInverseBeta(self):
        a = A.requiredAggregateFractionBeta(1.0, 0.35)
        b = A.requiredAggregateFractionBeta(0.5, 0.35)
        self.assertAlmostEqual(b[0] / a[0], 2.0, places=9)
        self.assertAlmostEqual(b[1] / a[1], 2.0, places=9)

    def testZeroBetaGivesAnUnboundedRequirement(self):
        lo, hi = A.requiredAggregateFractionBeta(0.0, 0.35)
        self.assertEqual(lo, float("inf"))
        self.assertEqual(hi, float("inf"))

    def testDecisionEntryRefusesTheParameterFreeClaim(self):
        doc = (_REPO_ROOT / "DECISIONS.md").read_text()
        entry = doc[doc.index("## D033 —"):]
        self.assertIn("not parameter-free", entry)
        self.assertIn("beta`-indexed interval", entry)
        # the weak bound must carry both caveats
        self.assertIn("our arithmetic, not their", entry)
        self.assertIn("wrong condition", entry)

    def testTheSectionTextDoesNotClaimParameterFreeness(self):
        p = _REPO_ROOT / "manuscript" / "SECTION_8_4_LINEAGE_PREDICTION.md"
        txt = p.read_text()
        self.assertIn("We do not claim this", txt)
        self.assertIn("parameter-free", txt)
        self.assertIn("closes this", txt)
