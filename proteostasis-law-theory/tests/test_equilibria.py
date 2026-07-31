"""equilibrium, stability and fold-detection tests.

the load-bearing check here is that the two independent equilibrium finders --
blind multi-start root solving and warm-started continuation -- agree. they
share the rhs but nothing else, so a disagreement localises a numerical
failure rather than a modelling one.
"""

import unittest

import numpy as np

import _context  # noqa: F401
from proteostasis import (Params, findEquilibria, lowestStableEquilibrium, findFold,
                          traceBranch, solveEquilibriumFrom, jacobian, rhsVector,
                          removalCeiling, simulate)
from proteostasis.equilibria import residualScale, stableEquilibriumAt

BASE = Params(j=0.02)


class TestEquilibriumSolving(unittest.TestCase):

    def testResidualsAreSmall(self):
        for j in (0.002, 0.02, 0.08, 0.14):
            for eq in findEquilibria(BASE.with_(j=j), n_grid=7):
                f = rhsVector((eq.u, eq.a), BASE.with_(j=j))
                self.assertLess(float(np.max(np.abs(f))) / residualScale(BASE.with_(j=j)),
                                1e-8, f"large residual at j={j}")

    def testStabilityMatchesJacobianEigenvalues(self):
        for eq in findEquilibria(BASE, n_grid=7):
            eig = np.linalg.eigvals(jacobian(eq.u, eq.a, BASE))
            self.assertEqual(eq.stable, bool(np.max(eig.real) < 0.0))

    def testEquilibriaAreStrictlyPositive(self):
        """with j > 0 no equilibrium may sit on the nonnegative boundary."""
        for eq in findEquilibria(BASE, n_grid=7):
            self.assertGreater(eq.u, 0.0)
            self.assertGreater(eq.a, 0.0)

    def testSaddleNodePairStructure(self):
        """below the fold: exactly one stable and one unstable state."""
        eqs = findEquilibria(BASE.with_(j=0.05), n_grid=9, hi=1e5)
        self.assertEqual(len(eqs), 2, f"expected a pair, got {len(eqs)}")
        self.assertEqual(sum(e.stable for e in eqs), 1)
        self.assertLess(eqs[0].burden, eqs[1].burden)
        self.assertTrue(eqs[0].stable)

    def testNoEquilibriumAboveRemovalCeiling(self):
        p = BASE.with_(j=1.5 * removalCeiling(BASE))
        self.assertEqual(findEquilibria(p, n_grid=9), [])


class TestContinuationAgreesWithBlindSearch(unittest.TestCase):

    def testAgreementAcrossInflux(self):
        j_grid = [0.001, 0.005, 0.02, 0.05, 0.09, 0.12, 0.15, 0.2, 0.4]
        for bp in traceBranch(BASE, "j", j_grid):
            blind = [e for e in findEquilibria(BASE.with_(j=bp.value), n_grid=9)
                     if e.stable]
            self.assertEqual(bool(blind), bp.eq is not None,
                             f"existence disagreement at j={bp.value}")
            if blind and bp.eq is not None:
                self.assertLess(abs(blind[0].u - bp.eq.u) / max(blind[0].u, 1e-30), 1e-4)
                self.assertLess(abs(blind[0].a - bp.eq.a) / max(blind[0].a, 1e-30), 1e-4)


class TestFoldDetection(unittest.TestCase):

    def setUp(self):
        self.fold = findFold(BASE, "j", 1e-4, 2.0)

    def testFoldIsFoundAndBracketed(self):
        self.assertTrue(self.fold.found, self.fold.reason)
        lo, hi = self.fold.bracket
        self.assertLess(lo, hi)
        self.assertLessEqual((hi - lo) / hi, 1e-5)

    def testStateExistsBelowAndNotAboveTheBracket(self):
        """the bracket must be verified by the INDEPENDENT blind finder.

        this is the check that caught continuation reporting a fold three
        orders of magnitude too low: a solver stall looks exactly like a fold
        unless a second method is asked.
        """
        lo, hi = self.fold.bracket
        self.assertIsNotNone(lowestStableEquilibrium(BASE.with_(j=lo * 0.999), n_grid=9))
        self.assertIsNone(lowestStableEquilibrium(BASE.with_(j=hi * 1.001), n_grid=9))

    def testFoldLiesBelowTheAnalyticRemovalCeiling(self):
        self.assertLess(self.fold.fold_value, removalCeiling(BASE))

    def testEquilibriumCountDropsTwoToZeroAcrossTheFold(self):
        lo, hi = self.fold.bracket
        self.assertEqual(len(findEquilibria(BASE.with_(j=lo * 0.99), n_grid=9, hi=1e5)), 2)
        self.assertEqual(len(findEquilibria(BASE.with_(j=hi * 1.01), n_grid=9, hi=1e5)), 0)

    def testNascentOccupancyShrinksTheFeasibleWindow(self):
        """nu carries no damage influx, so any shrinkage is capacity competition."""
        folds = [findFold(BASE.with_(nu=nu), "j", 1e-4, 2.0).fold_value
                 for nu in (0.0, 1.0, 4.0, 10.0)]
        self.assertTrue(all(v is not None for v in folds))
        self.assertTrue(all(b < a for a, b in zip(folds, folds[1:])),
                        f"j_fold not decreasing in nu: {folds}")

    def testNoFoldReportedWhenNoStateExistsAtAll(self):
        p = BASE.with_(j=1e-4, alpha_g=1e4, alpha_n=1e4)
        fr = findFold(p, "j", 1.0, 2.0)
        self.assertFalse(fr.found)
        self.assertIn("no stable low branch", fr.reason)


class TestReachability(unittest.TestCase):

    def testHealthyInitialConditionReachesTheLowBranch(self):
        eq = lowestStableEquilibrium(BASE)
        tr = simulate(BASE, 0.0, 0.0, t_end=5e4)
        self.assertEqual(tr.status, "converged")
        self.assertLess(abs(tr.final_u - eq.u) / eq.u, 1e-4)
        self.assertLess(abs(tr.final_a - eq.a) / eq.a, 1e-4)

    def testAboveTheUnstableStateTheSystemEscapes(self):
        """the saddle's unstable manifold separates recovery from escape."""
        eqs = findEquilibria(BASE, n_grid=9, hi=1e5)
        low = [e for e in eqs if e.stable][0]
        upper = [e for e in eqs if not e.stable][0]
        values, vectors = np.linalg.eig(jacobian(upper.u, upper.a, BASE))
        unstable = np.flatnonzero(values.real > 0.0)
        self.assertEqual(len(unstable), 1, f"not a planar saddle: {values}")
        direction = vectors[:, unstable[0]].real
        direction /= np.linalg.norm(direction)
        step = 0.01 * upper.burden

        sides = []
        for sign in (-1.0, 1.0):
            y0 = np.array([upper.u, upper.a]) + sign * step * direction
            self.assertTrue(np.all(y0 > 0.0))
            burden_rate = float(np.sum(rhsVector(y0, BASE)))
            sides.append((burden_rate, y0))
        high_rate, high_y0 = max(sides, key=lambda side: side[0])
        low_rate, low_y0 = min(sides, key=lambda side: side[0])
        self.assertGreater(high_rate, 0.0)
        self.assertLess(low_rate, 0.0)

        # A burden of 100 is already twenty times beyond the saddle.  The
        # influx-limited asymptotic growth cannot reach the generic 1e4 guard
        # within this test's finite horizon.
        high_tr = simulate(BASE, *high_y0, t_end=5e4, blowup=100.0)
        low_tr = simulate(BASE, *low_y0, t_end=5e4, blowup=100.0)
        self.assertEqual(high_tr.status, "blowup")
        self.assertEqual(low_tr.status, "converged")
        self.assertLess(abs(low_tr.final_u - low.u) / low.u, 1e-4)
        self.assertLess(abs(low_tr.final_a - low.a) / low.a, 1e-4)


if __name__ == "__main__":
    unittest.main()
