"""unit tests for the phase 2 audit machinery itself.

these run against the model, not against stored results, so they hold on a clean
checkout with no `results/` directory. their job is to make sure the audit's own
tools are trustworthy before anything is concluded from them.
"""

import unittest

import numpy as np

import _context  # noqa: F401
from _roots import (clusterRoots, continueBranch, denseMultiStart, matchRoot,
                    nullclineRoots, certifyRoot, Root)
from proteostasis import Params

#: the phase 1 default parameter set at an influx well inside its stable range.
SINGLE = Params(j=0.02)


def _root(u, a, stable=True, res=0.0):
    return Root(u=u, a=a, stable=stable, eig_real_max=-1.0 if stable else 1.0,
                eig_real_min=-2.0, eig_imag_max=0.0, residual_rel=res,
                residual_cancel=res, method="test")


class TestClustering(unittest.TestCase):

    def testDistinctRootsSurviveTightTolerance(self):
        rs = [_root(1.0, 1.0), _root(10.0, 5.0)]
        self.assertEqual(len(clusterRoots(rs, 1e-5)), 2)

    def testNearIdenticalRootsMerge(self):
        """two solver copies of one root must collapse, and the merge must be
        counted -- n_hits is how the audit distinguishes a well-found root from
        a lucky single hit."""
        rs = [_root(1.0, 1.0), _root(1.0 + 1e-12, 1.0 - 1e-12)]
        merged = clusterRoots(rs, 1e-5)
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0].n_hits, 2)

    def testMergingRequiresBothCoordinates(self):
        """agreement in u alone is not agreement: the phase 1 dedupe rule is a
        conjunction, and an audit that loosened it to a disjunction would erase
        genuinely distinct roots."""
        rs = [_root(1.0, 1.0), _root(1.0, 50.0)]
        self.assertEqual(len(clusterRoots(rs, 1e-3)), 2)

    def testLooserToleranceNeverIncreasesCount(self):
        rs = [_root(1.0, 1.0), _root(1.02, 1.02), _root(80.0, 3.0)]
        counts = [len(clusterRoots(rs, t)) for t in (1e-9, 1e-3, 1e-1)]
        self.assertEqual(counts, sorted(counts, reverse=True))

    def testMatchRoot(self):
        rs = [_root(1.0, 2.0)]
        self.assertIsNotNone(matchRoot(rs, 1.0 + 1e-6, 2.0 - 1e-6, rtol=1e-3))
        self.assertIsNone(matchRoot(rs, 1.5, 2.0, rtol=1e-3))


class TestIndependentFindersAgree(unittest.TestCase):
    """the three finders must agree on a case with a known, simple answer.

    if they disagree here, no disagreement they report on a candidate can be
    attributed to the candidate.
    """

    @classmethod
    def setUpClass(cls):
        cls.hybr = clusterRoots(denseMultiStart(SINGLE, n_grid=15, method="hybr"), 1e-5)
        cls.lm = clusterRoots(denseMultiStart(SINGLE, n_grid=15, method="lm"), 1e-5)
        cls.nl = clusterRoots(nullclineRoots(SINGLE, n_u=200, n_a=400)[0], 1e-5)

    def testAllThreeFindTheSameStableCount(self):
        n = [sum(r.stable for r in rs) for rs in (self.hybr, self.lm, self.nl)]
        self.assertEqual(len(set(n)), 1, f"finders disagree on stable count: {n}")

    def testEveryHybrRootIsFoundByTheBracketingMethod(self):
        """the newton-free method is the check against invented roots."""
        for r in self.hybr:
            self.assertIsNotNone(matchRoot(self.nl, r.u, r.a, rtol=1e-3),
                                 f"hybr root u={r.u:.6g} a={r.a:.6g} not bracketed")

    def testDefaultParametersAreMonostable(self):
        """phase 1's base parameter set is a single-attractor case; if this ever
        changes, every 'rare' statement in the audit needs rereading."""
        self.assertEqual(sum(r.stable for r in self.hybr), 1)

    def testRootsAreStrictlyPositiveAndCertify(self):
        for r in self.hybr:
            c = certifyRoot(r.u, r.a, SINGLE)
            self.assertTrue(c["positive"])
            self.assertTrue(c["closure_unique"])
            self.assertTrue(c["free_pools_bounded"])
            self.assertLess(c["residual_cancel"], 1e-9)
            self.assertLess(c["pool_balance_max"], 1e-10)
            self.assertLess(abs(c["mass_balance_residual"]), 1e-10)
            self.assertIsNotNone(c["jac_rel_error"])
            self.assertLess(c["jac_rel_error"], 1e-4,
                            "analytic jacobian disagrees with richardson "
                            "extrapolation; stability classification is unsafe")


class TestDeterminism(unittest.TestCase):

    def testDenseSearchIsDeterministic(self):
        a = clusterRoots(denseMultiStart(SINGLE, n_grid=11, method="hybr"), 1e-5)
        b = clusterRoots(denseMultiStart(SINGLE, n_grid=11, method="hybr"), 1e-5)
        self.assertEqual(len(a), len(b))
        for x, y in zip(a, b):
            self.assertEqual(x.u, y.u)
            self.assertEqual(x.a, y.a)


class TestContinuation(unittest.TestCase):

    def testContinuationTracksTheKnownBranch(self):
        """from a stable equilibrium, arclength continuation in j must stay on a
        branch of genuine equilibria (residual ~ 0 at every traced point)."""
        eq = [r for r in clusterRoots(denseMultiStart(SINGLE, n_grid=11), 1e-5)
              if r.stable][0]
        c = continueBranch(SINGLE, eq.u, eq.a, SINGLE.j, max_steps=200)
        self.assertGreater(c["n_points"], 20)
        self.assertLess(max(q["residual"] for q in c["points"]), 1e-8)

    def testContinuationRevealsAFoldOnAMonostableCase(self):
        """even a monostable parameter set has a fold in j -- that is what
        experiment C's j_fold is. continuation must find at least one."""
        eq = [r for r in clusterRoots(denseMultiStart(SINGLE, n_grid=11), 1e-5)
              if r.stable][0]
        c = continueBranch(SINGLE, eq.u, eq.a, SINGLE.j, max_steps=2000)
        self.assertGreaterEqual(c["n_folds"], 1)


if __name__ == "__main__":
    unittest.main()
