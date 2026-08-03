"""phase 3, antecedent check A1: capacity self-damage and the determinant identity.

the result these pin is a SURVIVAL, and survivals are the easy thing to fake by
running a check that cannot fail. so they also pin the reasons the check is weak:
that `eps = 0` is exact, that the residual is roundoff rather than identity error
(measured by its slope in the differencing step), and that the fold-finding
SHORTCUT dies even though the theorem does not.
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
import fold_theorem as FT  # noqa: E402
import capacity_self_damage as CS  # noqa: E402

_HAS_RUN = (FT.phase1RunDir() / "C" / "samples.tsv").is_file()


class TestZeroEpsilonIsTheFrozenModelExactly(unittest.TestCase):
    def testFactorIsExactlyOne(self):
        sd = CS.SelfDamage(0.0, CS.INFLUX).validate()
        self.assertEqual(sd.factor(0.3, 0.1, 0.2), 1.0)

    def testParamsAreReturnedUnchanged(self):
        p = M.Params().validate()
        for mode in (CS.INFLUX, CS.BURDEN):
            pe = CS.effectiveParams(p, 0.1, 0.2, CS.SelfDamage(0.0, mode).validate())
            self.assertEqual(pe.c_tot, p.c_tot)
            self.assertEqual(pe.d_tot, p.d_tot)

    def testFieldsMatchTheFrozenModelAtZero(self):
        p = M.Params().validate()
        sd = CS.SelfDamage(0.0, CS.INFLUX).validate()
        self.assertEqual(CS.removalRsd(0.1, 0.05, p, sd), FT.removalR(0.1, 0.05, p))
        self.assertEqual(CS.aggregateGsd(0.1, 0.05, p, sd), FT.aggregateG(0.1, 0.05, p))


class TestCapacityFallsAndStaysPositive(unittest.TestCase):
    def testFactorIsInUnitIntervalAndMonotoneDecreasing(self):
        p = M.Params(j=0.2).validate()
        prev = 1.0
        for eps in CS.EPS_LADDER:
            f = CS.SelfDamage(eps, CS.INFLUX).validate().factor(p.j, 0.1, 0.1)
            self.assertTrue(0.0 < f <= 1.0)
            self.assertLess(f, prev)
            prev = f

    def testBurdenModeRespondsToStateAndInfluxModeDoesNot(self):
        p = M.Params(j=0.2).validate()
        inf = CS.SelfDamage(1.0, CS.INFLUX).validate()
        bur = CS.SelfDamage(1.0, CS.BURDEN).validate()
        self.assertEqual(inf.factor(p.j, 0.1, 0.1), inf.factor(p.j, 9.0, 9.0))
        self.assertLess(bur.factor(p.j, 9.0, 9.0), bur.factor(p.j, 0.1, 0.1))


@unittest.skipUnless(_HAS_RUN, "phase 1 run artefacts absent (results/ is gitignored)")
class TestTheIdentitySurvivesSelfDamage(unittest.TestCase):
    """the headline: det J = -(grad R x grad G) holds over four decades of eps."""

    def testInfluxModeStaysAtTheFloor(self):
        df = CS.identityLadder(k=12, mode=CS.INFLUX)
        self.assertGreater(len(df), 4)
        self.assertLess(df["rel_err_grad_median"].max(), 1e-10)

    def testBurdenModeStaysAtTheFloorEvenThoughCapacityIsStateDependent(self):
        """this is the mode that could break it: capacity enters both gradients."""
        df = CS.identityLadder(k=12, mode=CS.BURDEN)
        self.assertLess(df["rel_err_grad_median"].max(), 1e-10)
        # and the coupling is genuinely severe at the top of the ladder
        self.assertLess(df["capacity_factor"].min(), 0.05)

    def testTheResidualIsRoundoffNotIdentityError(self):
        """slope ~ -1 in h. identity error would be flat in h; truncation +2.

        this is the test that stops 'the identity survives' from resting on a
        check that merely happened to be run at a forgiving step size.
        """
        d = CS.refinementCheck(k=12, mode=CS.BURDEN, eps=100.0,
                               steps=(1e-2, 1e-3, 1e-4, 1e-5, 1e-6))
        self.assertLess(d.attrs["slope_in_h"], -0.5)
        self.assertGreater(d.attrs["slope_in_h"], -1.5)


class TestTheShortcutDiesEvenThoughTheTheoremDoesNot(unittest.TestCase):
    def testFrozenFoldSolvesTwoEquationsAndSelfDamagedSolvesThree(self):
        p = M.Params().validate()
        base = FT.foldSolve(p)
        self.assertIsNotNone(base)
        out = CS.foldSolveSd(p, CS.SelfDamage(0.0, CS.INFLUX).validate())
        self.assertIsNotNone(out)
        # at eps = 0 the three-equation solve must reproduce the two-equation one
        for got, want in zip(out, base):
            self.assertAlmostEqual(got, want, delta=1e-6 * max(1.0, abs(want)))

    def testSelfConsistencyIsImposedNotEvaluated(self):
        """R = j must hold AT the solved fold, which is the third equation."""
        p = M.Params().validate()
        sd = CS.SelfDamage(1.0, CS.INFLUX).validate()
        out = CS.foldSolveSd(p, sd)
        if out is None:
            self.skipTest("no fold at this setting")
        j, u, a = out
        pe = CS.effectiveParams(p.with_(j=j), u, a, sd)
        self.assertAlmostEqual(FT.removalR(u, a, pe) / j, 1.0, places=6)
        self.assertLess(abs(FT.aggregateG(u, a, pe)), 1e-7)


class TestTheAnalyticCeilingBecomesSquareRoot(unittest.TestCase):
    def testItReducesToTheFrozenCeilingAtZero(self):
        p = M.Params().validate()
        sd = CS.SelfDamage(0.0, CS.INFLUX).validate()
        self.assertEqual(CS.influxCeilingSelfDamaged(p, sd), M.removalCeiling(p))

    def testItSolvesItsDefiningInequalityWithEquality(self):
        p = M.Params().validate()
        C0 = M.removalCeiling(p)
        for eps in CS.EPS_LADDER:
            j = CS.influxCeilingSelfDamaged(p, CS.SelfDamage(eps, CS.INFLUX).validate())
            self.assertAlmostEqual(j * (1.0 + eps * j), C0, places=10)

    def testItApproachesSqrtCoverEpsForLargeEps(self):
        p = M.Params().validate()
        C0 = M.removalCeiling(p)
        j = CS.influxCeilingSelfDamaged(p, CS.SelfDamage(1e6, CS.INFLUX).validate())
        self.assertAlmostEqual(j / np.sqrt(C0 / 1e6), 1.0, places=2)

    def testBurdenModeFallsBackToTheFrozenCeiling(self):
        p = M.Params().validate()
        sd = CS.SelfDamage(10.0, CS.BURDEN).validate()
        self.assertEqual(CS.influxCeilingSelfDamaged(p, sd), M.removalCeiling(p))


class TestValidation(unittest.TestCase):
    def testNegativeEpsAndUnknownModeAreRejected(self):
        with self.assertRaises(M.ModelError):
            CS.SelfDamage(-1.0, CS.INFLUX).validate()
        with self.assertRaises(M.ModelError):
            CS.SelfDamage(1.0, "nonsense").validate()


if __name__ == "__main__":
    unittest.main()
