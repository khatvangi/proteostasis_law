"""phase 3 follow-ups: uniqueness, the dilution-proof decomposition, bistability.

all model-level; no phase 1 artefacts needed except the threshold sweep, which
is not tested here because it is a distribution rather than an invariant.
imports are set inline for the reason documented in `test_fold_theorem.py`.
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


class TestDecompositionSurvivesDivision(unittest.TestCase):
    """j_crit = C_enz . phi_enz / (1 - delta), exactly.

    this is the replacement for the old margin, whose denominator stops being an
    upper bound once cells divide.
    """

    def testIdentityIsExact(self):
        p = M.Params().validate()
        for mu in (0.0, 0.02, 0.04, 0.06):
            d = B.boundaryDecomposition(p, D.Growth(mu0=mu))
            self.assertIsNotNone(d, f"no boundary at mu={mu}")
            self.assertLess(d["identity_rel_err"], 1e-12, f"mu={mu}: {d}")

    def testBothComponentsAreDimensionlessFractions(self):
        p = M.Params().validate()
        for mu in (0.0, 0.04, 0.08):
            d = B.boundaryDecomposition(p, D.Growth(mu0=mu))
            self.assertIsNotNone(d)
            for k in ("phi_enz", "delta"):
                self.assertGreaterEqual(d[k], 0.0)
                self.assertLess(d[k], 1.0)

    def testZeroDilutionGivesZeroDilutionShare(self):
        p = M.Params().validate()
        d = B.boundaryDecomposition(p, D.Growth(mu0=0.0))
        self.assertEqual(d["delta"], 0.0)
        # and then the decomposition reduces to the undiluted margin
        self.assertAlmostEqual(d["phi_enz"], d["j_crit"] / M.removalCeiling(p),
                               places=12)

    def testDilutionShareRisesWithDilution(self):
        p = M.Params().validate()
        deltas = []
        for mu in (0.0, 0.02, 0.04, 0.06, 0.08):
            d = B.boundaryDecomposition(p, D.Growth(mu0=mu))
            self.assertIsNotNone(d)
            deltas.append(d["delta"])
        self.assertTrue(all(x < y for x, y in zip(deltas, deltas[1:])), deltas)


class TestLinearGrowthLaw(unittest.TestCase):
    """the proteome-partitioning shape, as a second FORM (not a calibration)."""

    def testReducesToConstantWhenNoFeedback(self):
        g = B.LinearGrowth(mu0=0.3, k_mu=float("inf"))
        self.assertEqual(g.rate(0.0, 0.0), 0.3)
        self.assertEqual(g.rate(5.0, 5.0), 0.3)

    def testReachesZeroAtFiniteBurden(self):
        """unlike the hyperbolic form, growth actually arrests."""
        g = B.LinearGrowth(mu0=0.3, k_mu=2.0)
        self.assertAlmostEqual(g.rate(1.0, 1.0), 0.0, places=15)
        self.assertEqual(g.rate(5.0, 5.0), 0.0)
        self.assertGreater(D.Growth(mu0=0.3, k_mu=2.0).rate(5.0, 5.0), 0.0)

    def testGradientMatchesCentralDifference(self):
        g = B.LinearGrowth(mu0=0.3, k_mu=2.0)
        u, a, h = 0.3, 0.2, 1e-7
        num = (g.rate(u + h, a) - g.rate(u - h, a)) / (2 * h)
        self.assertAlmostEqual(g.rateGradient(u, a)[0], num, places=7)

    def testBoundarySurvivesUnderEitherForm(self):
        """the qualitative conclusion must not depend on the functional form."""
        p = M.Params().validate()
        for mu0 in (0.1, 0.3, 0.6):
            self.assertIsNotNone(D.foldSolveDil(p, D.Growth(mu0, 0.5)),
                                 f"hyperbolic law lost the boundary at {mu0}")
            self.assertIsNotNone(D.foldSolveDil(p, B.LinearGrowth(mu0, 2.0)),
                                 f"linear law lost the boundary at {mu0}")


class TestUniqueness(unittest.TestCase):
    def testUndilutedCriticalPointIsUnique(self):
        """with no dilution the nullcline closes and carries exactly one."""
        p = M.Params().validate()
        r = B.criticalPointsOnNullcline(p, D.Growth(0.0))
        self.assertGreater(r["n_lower"], 50)
        # a properly closed curve has comparable branch lengths
        self.assertAlmostEqual(r["n_lower"] / r["n_upper"], 1.0, delta=0.2)
        self.assertEqual(r["n_sign_changes"], 1)

    def testDilutionAdmitsASecondDistinctCriticalPoint(self):
        """uniqueness FAILS once cells divide -- both points are genuine."""
        p = M.Params().validate()
        g = D.Growth(0.04)
        pts = []
        for seed in ((0.45, 0.35), (0.16, 1.95)):
            got = B._solveSeeded(p, g, seed)
            self.assertIsNotNone(got, f"seed {seed} did not converge")
            u, a = got
            self.assertLess(abs(D.aggregateGDil(u, a, p, g)), 1e-9)
            self.assertLess(abs(float(np.linalg.det(D.jacobianDil(u, a, p, g)))),
                            1e-7)
            pts.append((D.removalTotal(u, a, p, g), u, a))
        self.assertGreater(abs(pts[0][1] - pts[1][1]) / pts[0][1], 1e-3,
                           "the two seeds converged to the same point")


class TestDilutionMakesTheSystemBistable(unittest.TestCase):
    """losing the low branch is a jump to a bounded state, not divergence."""

    @classmethod
    def setUpClass(cls):
        cls.p = M.Params().validate()
        cls.g = D.Growth(0.04)
        cls.js = [0.10, 0.155, 0.17, 0.19, 0.196, 0.21]
        cls.h = B.hysteresisSweep(cls.p, cls.g, cls.js)

    def testHighBurdenStateIsBounded(self):
        u, a = self.h["up"][0.21]
        self.assertTrue(np.isfinite(u) and np.isfinite(a))
        self.assertLess(u + a, 1e3)
        self.assertGreater(a, 1.0, "expected a high-aggregate state")

    def testSweepsDisagreeOverAnInterval(self):
        self.assertTrue(self.h["hysteretic_j"], "no hysteresis detected")

    def testWindowIsBracketedByTheTwoCriticalPoints(self):
        lo = B._solveSeeded(self.p, self.g, (0.16, 1.95))
        hi = B._solveSeeded(self.p, self.g, (0.45, 0.35))
        j_lo = D.removalTotal(lo[0], lo[1], self.p, self.g)
        j_hi = D.removalTotal(hi[0], hi[1], self.p, self.g)
        w0, w1 = self.h["window"]
        self.assertGreaterEqual(w0, min(j_lo, j_hi))
        self.assertLessEqual(w1, max(j_lo, j_hi))


if __name__ == "__main__":
    unittest.main()
