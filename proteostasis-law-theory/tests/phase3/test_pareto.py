"""phase 3 Pareto layer: the trade-off surface the framework asserted but never built.

model-level only; no artefacts needed. imports are set inline for the reason
documented in `test_fold_theorem.py`.
"""

import sys
import unittest
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (_REPO_ROOT / "scripts", _REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import pareto as P  # noqa: E402


class TestEconomyIsAGenuineTradeOff(unittest.TestCase):
    """both strategy coordinates must cost throughput, or there is no trade-off."""

    def testThroughputFallsInBothCoordinates(self):
        e = P.Economy()
        base = e.throughput(0.3, 0.05)
        self.assertLess(e.throughput(0.6, 0.05), base)
        self.assertLess(e.throughput(0.3, 0.10), base)

    def testAccuracyBuysLowerErrorRate(self):
        e = P.Economy()
        eps = [e.errorRate(a) for a in (0.0, 0.25, 0.5, 0.75, 1.0)]
        self.assertTrue(all(x > y for x, y in zip(eps, eps[1:])), eps)

    def testInfluxFallsWithAccuracy(self):
        """the exponential must dominate the linear throughput term."""
        e = P.Economy()
        js = [e.influx(a, 0.05) for a in (0.0, 0.25, 0.5, 0.75, 1.0)]
        self.assertTrue(all(x > y for x, y in zip(js, js[1:])), js)


class TestConstraintIsTheDerivedOne(unittest.TestCase):
    def testCriticalInfluxRisesWithQCInvestment(self):
        e = P.Economy()
        jc = [P.criticalInflux(R, e) for R in (0.02, 0.04, 0.08, 0.12)]
        self.assertTrue(all(x is not None for x in jc), jc)
        self.assertTrue(all(x < y for x, y in zip(jc, jc[1:])), jc)

    def testCriticalInfluxComesFromTheFoldTheorem(self):
        """j_crit here must equal foldSolve at the same rescue pools."""
        from fold_theorem import foldSolve
        e = P.Economy()
        R, chi = 0.05, 0.6
        tot = e.rescueTotal(R)
        p = M.Params().with_(c_tot=chi * tot, d_tot=(1 - chi) * tot).validate()
        self.assertAlmostEqual(P.criticalInflux(R, e, chi), foldSolve(p)[0],
                               places=12)


class TestOptimumBindsTheConstraint(unittest.TestCase):
    """the substantive result: throughput is limited by proteostasis, exactly."""

    @classmethod
    def setUpClass(cls):
        cls.exact = P.optimumOnBoundary(n_R=60)

    def testAnOptimumExists(self):
        self.assertIsNotNone(self.exact)
        self.assertGreater(self.exact["throughput"], 0.0)

    def testItSitsExactlyOnTheFeasibilityBoundary(self):
        self.assertAlmostEqual(self.exact["margin_ratio"], 1.0, places=6)

    def testItBeatsTheGridOptimum(self):
        """a grid can only under-report the optimum, never exceed it."""
        grid = P.optimum(P.strategyGrid(n_alpha=12, n_R=10))
        self.assertIsNotNone(grid)
        self.assertGreaterEqual(self.exact["throughput"] + 1e-12,
                                grid["throughput"])

    def testBoundaryTracerIsConsistentWithFeasibility(self):
        e = P.Economy()
        for R in (0.03, 0.05, 0.08):
            a = P.minimumAccuracyFor(R, e)
            self.assertIsNotNone(a)
            # just above the boundary is feasible, just below is not
            self.assertTrue(P.evaluate(a + 1e-4, R, e)["feasible"])
            if a > 1e-3:
                self.assertFalse(P.evaluate(a - 1e-3, R, e)["feasible"])


class TestFrontSpansAMarginRange(unittest.TestCase):
    """not every viable strategy sits near the boundary -- only the fast end."""

    def testMarginVariesAlongTheFront(self):
        df = P.strategyGrid(n_alpha=16, n_R=14)
        front = P.paretoFront(df)
        self.assertGreater(len(front), 3)
        self.assertLess(front["margin_ratio"].min(), 0.6)
        self.assertGreater(front["margin_ratio"].max(), 0.85)

    def testFrontIsNonDominated(self):
        df = P.strategyGrid(n_alpha=12, n_R=10)
        front = P.paretoFront(df)
        for _, r in front.iterrows():
            better = ((front["throughput"] >= r["throughput"])
                      & (front["accuracy"] >= r["accuracy"])
                      & ((front["throughput"] > r["throughput"])
                         | (front["accuracy"] > r["accuracy"])))
            self.assertFalse(bool(better.any()), "front contains a dominated point")


if __name__ == "__main__":
    unittest.main()
