"""the boron <-> nitrogen mapping is an identity, not a resemblance.

these tests assert the algebra of `theory/MATCHED_IMPLEMENTATION_PROTOCOL.md`
section 3 directly, so that a typo in one coefficient cannot survive by being
numerically small.
"""

import unittest

import numpy as np

import _context  # noqa: F401
from phase2.lhs import parametersForEpsilon, sampleMatrix
from phase2.mapping import (D_TOT_GAUGE, EPSILON_LADDER, MU_GAUGE,
                            bindingConstants, boronToNitrogen,
                            boundSubstrateFraction, nitrogenToBoron,
                            sequestrationEpsilon)
from phase2.models import BoronAdapter, FreeLimitAdapter
from proteostasis.model import Params

MATRIX = sampleMatrix(24, 20260731)


class TestLadderAlgebra(unittest.TestCase):

    def testBindingConstantsAreExactReciprocals(self):
        """c_u = 1/kappa_cu is what makes both arms share one occupancy term."""
        for eps in EPSILON_LADDER:
            for c_tot, d_tot in ((0.6, 0.4), (1.0, 1.0), (3.7, 1.0)):
                kb = bindingConstants(eps, c_tot, d_tot)
                self.assertAlmostEqual(kb["c_u"] * kb["kappa_cu"], 1.0, delta=1e-14)
                self.assertAlmostEqual(kb["c_a"] * kb["kappa_ca"], 1.0, delta=1e-14)
                self.assertAlmostEqual(kb["d_u"] * kb["kappa_du"], 1.0, delta=1e-14)
                self.assertAlmostEqual(kb["d_a"] * kb["kappa_da"], 1.0, delta=1e-14)

    def testEpsilonZeroIsRefusedAsBoronParams(self):
        with self.assertRaises(ValueError):
            bindingConstants(0.0, 0.6, 0.4)

    def testLadderReproducesItsOwnEpsilon(self):
        """the audit's epsilon read back off the constructed Params."""
        for eps in EPSILON_LADDER:
            q = parametersForEpsilon(MATRIX[3], eps)
            p = Params(**nitrogenToBoron(q, eps)).validate()
            self.assertAlmostEqual(sequestrationEpsilon(p) / eps, 1.0, places=12)

    def testKappaCaEqualsKappaCuOnTheLadder(self):
        """su == sa is what collapses the discrepancy to a single scalar."""
        for eps in EPSILON_LADDER:
            p = Params(**nitrogenToBoron(parametersForEpsilon(MATRIX[7], eps), eps))
            self.assertEqual(p.kappa_ca, p.kappa_cu)
            self.assertEqual(p.kappa_da, p.kappa_du)


class TestRoundTrip(unittest.TestCase):

    def testBoronToNitrogenToBoronIsIdentity(self):
        base = dict(j=0.02, nu=0.5, c_tot=0.6, d_tot=0.4, rho_U=1.0, rho_A=0.5,
                    alpha_n=0.5, alpha_g=1.0, alpha_d=0.3, m=2.0, kappa_ref=0.5,
                    kappa_u=0.5, kappa_a=0.5, kappa_dis=0.5)
        for eps in EPSILON_LADDER:
            kb = bindingConstants(eps, base["c_tot"], base["d_tot"])
            p = Params(**base, **{k: kb[k] for k in
                                  ("kappa_cu", "kappa_ca", "kappa_du", "kappa_da")})
            back = nitrogenToBoron(boronToNitrogen(p, eps), eps)
            # d_tot is a gauge freedom the nitrogen groups do not fix, so the
            # round trip lands in the D_TOT_GAUGE gauge; every RATE must survive.
            self.assertAlmostEqual(back["d_tot"], D_TOT_GAUGE, places=15)
            self.assertAlmostEqual(back["c_tot"], p.c_tot, places=15)
            self.assertAlmostEqual(back["rho_U"] * back["d_tot"],
                                   p.rho_U * p.d_tot, places=15)
            self.assertAlmostEqual(back["rho_A"] * back["d_tot"],
                                   p.rho_A * p.d_tot, places=15)
            self.assertAlmostEqual(back["alpha_d"] * back["c_tot"],
                                   p.alpha_d * p.c_tot, places=15)
            for k in ("j", "nu", "alpha_n", "alpha_g", "m",
                      "kappa_ref", "kappa_u", "kappa_a", "kappa_dis"):
                self.assertAlmostEqual(back[k], getattr(p, k), places=15)

    def testMuGaugeIsUnity(self):
        """T0 and the benchmark both assume boron and nitrogen clocks agree."""
        self.assertEqual(MU_GAUGE, 1.0)


class TestBenchmarkConsistency(unittest.TestCase):

    def testAdaptersSeeTheSameBindingConstants(self):
        """the free arm's c_u and the boron arm's kappa_cu must be reciprocal.

        if they were not, the two arms would differ in resource occupancy as
        well as in sequestration, and the factorial would be confounded.
        """
        for eps in (1e-6, 1e-2, 2.0):
            for i in (0, 5, 11, 23):
                q = parametersForEpsilon(MATRIX[i], eps)
                b = BoronAdapter(q, eps)
                f = FreeLimitAdapter(q, eps)
                self.assertAlmostEqual(f.params.c_u * b.params.kappa_cu, 1.0, delta=1e-13)
                self.assertAlmostEqual(f.params.d_u * b.params.kappa_du, 1.0, delta=1e-13)

    def testRemovalCeilingsAgree(self):
        """both arms must apply the SAME relative residual scale."""
        for eps in (1e-6, 0.3, 2.0):
            for i in (0, 9, 20):
                q = parametersForEpsilon(MATRIX[i], eps)
                b, f = BoronAdapter(q, eps), FreeLimitAdapter(q, eps)
                self.assertAlmostEqual(b.removalCeiling(), f.removalCeiling(), places=12)
                self.assertAlmostEqual(b.residualScale(), f.residualScale(), places=12)


class TestBoundFraction(unittest.TestCase):

    def testMonotoneAndBounded(self):
        q = parametersForEpsilon(MATRIX[0], 1.0)
        vals = [boundSubstrateFraction(q, e) for e in EPSILON_LADDER]
        self.assertTrue(all(0.0 <= v < 1.0 for v in vals))
        self.assertTrue(np.all(np.diff(vals) > 0.0))

    def testVanishesAtTheAnchor(self):
        q = parametersForEpsilon(MATRIX[0], 1e-6)
        self.assertLess(boundSubstrateFraction(q, 1e-6), 1e-5)


if __name__ == "__main__":
    unittest.main()
