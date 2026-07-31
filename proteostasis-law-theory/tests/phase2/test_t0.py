"""T0 -- the gate.  nothing downstream of this is meaningful if it fails.

the pass criteria here are the ones stated in
`theory/MATCHED_IMPLEMENTATION_PROTOCOL.md` section 5, derived a priori from the
O(epsilon) expansion rather than fitted to the observed numbers.  the slope
tolerances in particular are bounded by |r| <= 2 on the second-order coefficient,
which follows from su - 1 = epsilon*sigma with sigma <= 2.
"""

import unittest

import numpy as np

import _context  # noqa: F401
from phase2 import t0_equivalence as t0


class TestT0(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.report = t0.runT0()

    def testT0Passes(self):
        failed = [k for k, v in self.report["checks"].items() if not v]
        self.assertEqual(failed, [], t0.failureReport(self.report))

    def testAnchorRhsMeetsTheRequiredTolerance(self):
        """the task's hard requirement: < 1e-5 relative at epsilon = 1e-6."""
        anchor = next(c for c in self.report["cells"] if c["epsilon"] == 1e-6)
        self.assertLess(anchor["rhs_rel_max"], 1e-5)

    def testAnchorJacobianMeetsTheRequiredTolerance(self):
        anchor = next(c for c in self.report["cells"] if c["epsilon"] == 1e-6)
        self.assertLess(anchor["jac_rel_max"], 1e-5)

    def testDiscrepancyIsFirstOrderNotZerothOrSecond(self):
        """an exponent of 0 or 2 would falsify the mapping, in opposite ways."""
        s = self.report["scaling"]
        for key in ("rhs_slope_loglog", "jac_slope_loglog"):
            self.assertAlmostEqual(s[key], 1.0, delta=t0.SLOPE_TOL_FULL, msg=key)
        for key in ("rhs_slope_asymptotic", "jac_slope_asymptotic"):
            self.assertAlmostEqual(s[key], 1.0, delta=t0.SLOPE_TOL_ASYMPTOTIC, msg=key)

    def testLogLogRelationIsLinear(self):
        s = self.report["scaling"]
        self.assertGreaterEqual(s["rhs_r2"], t0.R2_MIN)
        self.assertGreaterEqual(s["jac_r2"], t0.R2_MIN)

    def testDiscrepancyIsMonotoneInEpsilon(self):
        s = self.report["scaling"]
        self.assertTrue(np.all(np.diff(s["rhs_rel_max"]) > 0))
        self.assertTrue(np.all(np.diff(s["jac_rel_max"]) > 0))

    def testGridIsTheTwelveStatesOfExperimentA(self):
        for c in self.report["cells"]:
            self.assertEqual(len(c["per_state"]), 12)

    def testDeterministic(self):
        again = t0.runT0()
        for a, b in zip(self.report["cells"], again["cells"]):
            self.assertEqual(a["rhs_rel_max"], b["rhs_rel_max"])
            self.assertEqual(a["jac_rel_max"], b["jac_rel_max"])

    def testMappedParametersMatchTheAuditTable(self):
        """the audit's section 5.2 numbers, recomputed from the config."""
        anchor = next(c for c in self.report["cells"] if c["epsilon"] == 1e-6)
        q = anchor["nitrogen_params"]
        expect = {"ref_capacity": 0.6, "deg_u_capacity": 0.4, "deg_a_capacity": 0.2,
                  "disaggregation": 0.18, "nucleation": 0.5, "growth": 1.0,
                  "ref_k": 0.5, "deg_u_k": 0.5, "deg_a_k": 0.5,
                  "disaggregation_k": 0.5, "n_load": 0.5, "j_u": 0.02}
        for k, v in expect.items():
            self.assertAlmostEqual(q[k], v, places=12, msg=k)
        self.assertAlmostEqual(q["c_u"], 1e-6 / 0.6, places=18)
        self.assertAlmostEqual(q["d_u"], 1e-6 / 0.4, places=18)


class TestFailureReporting(unittest.TestCase):
    """a failing T0 must name the mismatched terms, not just say 'failed'."""

    def testFailureReportNamesStatesAndFluxes(self):
        broken = dict(self.__class__._report())
        text = t0.failureReport(broken)
        self.assertIn("T0 FAILED", text)
        self.assertIn("F_boron", text)
        self.assertIn("worst state", text)

    @staticmethod
    def _report():
        r = t0.runT0()
        # force the report into its failed shape without touching the models
        r["checks"] = {k: False for k in r["checks"]}
        for c in r["cells"]:
            c["rhs_rel_max"] = 1.0
        return r


if __name__ == "__main__":
    unittest.main()
