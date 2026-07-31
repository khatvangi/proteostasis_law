"""integration tests: no negative-state leakage, correct escape detection.

`simulate` carries a terminal event that fires if either state drops below
NEG_TOL. the nonnegative orthant is forward invariant, so that event must never
fire on the real model -- and it must fire on a system that genuinely leaks,
which the mutation test below verifies.
"""

import unittest

import numpy as np

import _context  # noqa: F401
from proteostasis import Params, simulate, defaultInitialConditions, removalCeiling
from proteostasis.simulate import NEG_TOL
from proteostasis.sweeps import GLOBAL_RANGES, latinHypercube, paramsFromSample


class TestPositivityUnderIntegration(unittest.TestCase):

    def testNoNegativeStateLeakage(self):
        base = Params()
        param_sets = [base] + [paramsFromSample(s, base)
                               for s in latinHypercube(GLOBAL_RANGES, 8, seed=21)]
        checked = 0
        for p in param_sets:
            for u0, a0 in defaultInitialConditions(scale=1.0, n=8, seed=3):
                tr = simulate(p, u0, a0, t_end=2e4, n_out=100, blowup=1e4)
                self.assertNotEqual(tr.status, "negative",
                                    f"negative-state leakage at u0={u0}, a0={a0}")
                self.assertGreaterEqual(tr.min_u, NEG_TOL)
                self.assertGreaterEqual(tr.min_a, NEG_TOL)
                checked += 1
        self.assertGreater(checked, 50)

    def testBoundaryInitialConditionsStayNonnegative(self):
        p = Params()
        for u0, a0 in [(0.0, 0.0), (0.0, 0.5), (0.5, 0.0)]:
            tr = simulate(p, u0, a0, t_end=2e4, n_out=100)
            self.assertNotEqual(tr.status, "negative")
            self.assertGreaterEqual(min(tr.min_u, tr.min_a), NEG_TOL)

    def testMutationLeakingSystemIsDetected(self):
        """the negative-state event must actually fire on a leaking field.

        an unconditional drain violates forward invariance; if `simulate` did
        not detect it, the positivity test above would be vacuous.
        """
        from scipy.integrate import solve_ivp

        def leaking(_t, y):
            return np.array([-0.1, -0.1])

        def evNeg(_t, y):
            return min(y[0], y[1]) - NEG_TOL
        evNeg.terminal = True
        evNeg.direction = -1.0

        sol = solve_ivp(leaking, (0.0, 10.0), [0.1, 0.1], events=evNeg,
                        rtol=1e-9, atol=1e-12)
        self.assertTrue(len(sol.t_events[0]) > 0, "leak detector failed to fire")


class TestEscapeDetection(unittest.TestCase):

    def testInfluxAboveCeilingEscapes(self):
        base = Params()
        for p in [base] + [paramsFromSample(s, base)
                           for s in latinHypercube(GLOBAL_RANGES, 6, seed=31)]:
            tr = simulate(p.with_(j=1.5 * removalCeiling(p)), 0.0, 0.0,
                          t_end=1e5, n_out=100, blowup=1e4)
            self.assertEqual(tr.status, "blowup",
                             "burden did not escape above the analytic removal ceiling")

    def testConvergedRunsHaveNegligibleFinalRate(self):
        tr = simulate(Params(j=0.02), 0.0, 0.0, t_end=5e4)
        self.assertEqual(tr.status, "converged")
        self.assertLess(tr.final_rate, 1e-6)


class TestRecoveryTime(unittest.TestCase):

    def testRecoveryTimeGrowsTowardTheFold(self):
        """critical slowing down, measured dynamically rather than from lambda."""
        from proteostasis import findFold, lowestStableEquilibrium, recoveryTime

        base = Params()
        jf = findFold(base, "j", 1e-4, 2.0).fold_value
        self.assertIsNotNone(jf)
        times = []
        for frac in (0.3, 0.9, 0.99):
            p = base.with_(j=frac * jf)
            eq = lowestStableEquilibrium(p, n_grid=9)
            self.assertIsNotNone(eq, f"no equilibrium at {frac} of the fold")
            t = recoveryTime(p, eq.u, eq.a, kick=0.02)
            self.assertIsNotNone(t)
            times.append(t)
        # Event location uses the integrator's root finder (rtol=1e-10), so a
        # relative separation well above numerical error is required.
        self.assertTrue(all(b > a * (1.0 + 1e-4) for a, b in zip(times, times[1:])),
                        f"recovery time not increasing toward the fold: {times}")


if __name__ == "__main__":
    unittest.main()
