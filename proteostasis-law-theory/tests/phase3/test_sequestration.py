"""phase 3: spatial sequestration, the mechanism a failed post-diction named.

these tests are written against the PREREGISTRATION (D028), not against the
results. the structural half -- reduction at k_seq = 0, mass balance, the D024
identity, and the band being the one that was declared -- is what stops the
verdict from being reachable by editing the criterion after the run.
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
import boundary_structure as B  # noqa: E402
import fold_theorem as FT  # noqa: E402
import sequestration as S  # noqa: E402

_P = M.Params().validate()
_STATES = ((0.12, 0.08, 0.05), (0.30, 0.20, 0.90), (0.02, 0.01, 0.004))
_LAWS = (None, D.Growth(0.1, 0.5), B.LinearGrowth(0.1, 0.5))


class TestReductionToTheTwoStateModel(unittest.TestCase):
    """k_seq = 0 must be the frozen model EXACTLY, as the dilution work asserts."""

    def testUndilutedReductionIsExact(self):
        z = S.Sequestration(0.0, 0.0, 1.0).validate()
        du, da, das = S.rhs3s(0.12, 0.08, 0.0, _P, z, None)
        ref = M.rhs(0.12, 0.08, _P)
        self.assertEqual(du, ref[0])
        self.assertEqual(da, ref[1])
        self.assertEqual(das, 0.0)

    def testDilutedReductionIsExact(self):
        z = S.Sequestration(0.0, 0.0, 1.0).validate()
        g = D.Growth(0.1, 0.5)
        du, da, _ = S.rhs3s(0.12, 0.08, 0.0, _P, z, g)
        ref = D.rhsDil(0.12, 0.08, _P, g)
        self.assertEqual(du, ref[0])
        self.assertEqual(da, ref[1])

    def testZeroRatesGiveZeroFluxes(self):
        self.assertEqual(S.Sequestration(0.0, 0.0, 1.0).validate().fluxes(3.0, 4.0),
                         (0.0, 0.0))


class TestTheMechanismIsKineticNotBookkeeping(unittest.TestCase):
    def testMassBalanceIsUntouchedBySequestration(self):
        """du + da_r + da_s = j - R exactly: the k_seq/k_rel pair must cancel."""
        for seq in S.SEQ_GRID:
            for u, ar, asq in _STATES:
                f = S.rhs3s(u, ar, asq, _P, seq.validate(), None)
                want = _P.j - FT.removalR(u, ar, _P)
                self.assertAlmostEqual(sum(f), want, places=12)

    def testMassBalanceHoldsUnderDilutionToo(self):
        g = D.Growth(0.1, 0.5)
        for seq in S.SEQ_GRID:
            for u, ar, asq in _STATES:
                f = S.rhs3s(u, ar, asq, _P, seq.validate(), g)
                want = _P.j - S.removalR3s(u, ar, asq, _P, g)
                self.assertAlmostEqual(sum(f), want, places=12)

    def testSequesteredAggregateDoesNotBindMachineryOrNucleate(self):
        """raising a_s alone must not change any reactive flux."""
        seq = S.Sequestration(0.0, 0.0, 1.0).validate()   # no transfer either way
        a = S.rhs3s(0.12, 0.08, 0.0, _P, seq, None)
        b = S.rhs3s(0.12, 0.08, 7.5, _P, seq, None)
        self.assertEqual(a[0], b[0])
        self.assertEqual(a[1], b[1])

    def testReleaseReturnsSequesteredAggregateToTheReactivePool(self):
        seq = S.Sequestration(0.0, 0.5, 1.0).validate()
        _, dar, das = S.rhs3s(0.12, 0.08, 2.0, _P, seq, None)
        self.assertAlmostEqual(dar - M.rhs(0.12, 0.08, _P)[1], 1.0, places=12)
        self.assertAlmostEqual(das, -1.0, places=12)


class TestD024AppliesAndIsVerifiedNotAssumed(unittest.TestCase):
    """D028 requires the identity be checked BEFORE any number is interpreted."""

    def testIdentityHoldsOnTheExtendedSystem(self):
        errs = []
        for seq in S.SEQ_GRID:
            for g in _LAWS:
                for st in _STATES:
                    for costs in (True, False):
                        errs.append(S.identity3s(*st, _P, seq.validate(), g,
                                                 costs)["rel_err"])
        errs = np.array(errs)
        self.assertGreater(len(errs), 100)
        self.assertLess(float(np.median(errs)), 1e-9)
        self.assertLess(float(errs.max()), 1e-6)


class TestThePreregisteredBandIsTheOneDeclared(unittest.TestCase):
    """the criterion must not be reachable by editing it after the run."""

    def testEdgesMatchD028(self):
        self.assertEqual((S.BAND_LO, S.BAND_HI), (0.30, 0.60))

    def testDecisionEntryStatesBothEdgesAndTheConflict(self):
        doc = (_REPO_ROOT / "DECISIONS.md").read_text()
        i = doc.index("## D028 —")
        entry = doc[i: doc.find("\n## ", i + 5) if doc.find("\n## ", i + 5) > 0
                    else len(doc)]
        self.assertIn("0.30 <= (1 - mu_high/mu_low) <= 0.60", entry)
        self.assertIn("Upper edge 0.60", entry)
        self.assertIn("Declared conflict of interest with a past verdict", entry)

    def testBandIsClosedOnBothSides(self):
        self.assertFalse(S.inBand(0.29))
        self.assertTrue(S.inBand(0.30))
        self.assertTrue(S.inBand(0.60))
        self.assertFalse(S.inBand(0.61))
        self.assertFalse(S.inBand(0.95))

    def testConstantDilutionIsDisqualifiedInAdvance(self):
        """protocol rule 4: a match under a rejected regime is a failure."""
        self.assertIn("constant", S.DISQUALIFIED)


class TestVerdictCannotBeCarriedByADisqualifiedRegime(unittest.TestCase):
    def testPassRequiresAQualifiedLaw(self):
        import pandas as pd
        df = pd.DataFrame([
            {"law": "constant", "bistable": True, "in_band": True},
            {"law": "hyperbolic", "bistable": False, "in_band": False},
        ])
        v = S.verdict(df)
        self.assertFalse(v["passes"])
        self.assertEqual(v["disqualified_in_band"], 1)
        self.assertEqual(v["qualified_in_band"], 0)


class TestValidation(unittest.TestCase):
    def testNegativeRatesAndSublinearOrderAreRejected(self):
        with self.assertRaises(M.ModelError):
            S.Sequestration(-1.0, 0.0, 1.0).validate()
        with self.assertRaises(M.ModelError):
            S.Sequestration(1.0, 0.0, 0.5).validate()


class TestProtocolIsRecorded(unittest.TestCase):
    def testProtocolFileStatesTheFiveRulesAndTheRetroApplication(self):
        doc = (_REPO_ROOT / "notes" / "POSTDICTION_PROTOCOL.md").read_text()
        self.assertIn("Find the measured NUMBER in the paper", doc)
        self.assertIn("D.Growth(0.04)", doc)
        self.assertIn("State the upper edge", doc)
        self.assertRegex(doc, r"counts as a FAILURE")


if __name__ == "__main__":
    unittest.main()
