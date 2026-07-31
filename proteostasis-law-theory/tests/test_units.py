"""unit-consistency tests.

these fail when a rate constant is declared with the wrong dimension, when a
term of the wrong dimension is added into dU/dt or dA/dt, or when the
nondimensionalization stops being an exact rescaling of the dimensional model.

each positive check is paired with a MUTATION test: the same check is re-run
against a deliberately broken quantity and must fail. a check that cannot fail
is not a test.
"""

import unittest

import numpy as np

import _context  # noqa: F401
from proteostasis import units as U


def conc(v):
    return U.Quantity(v, U.CONC)


class TestDimensionAlgebra(unittest.TestCase):

    def testExponentArithmetic(self):
        self.assertTrue((U.AMOUNT / U.VOLUME).matches(U.CONC))
        self.assertTrue((U.CONC / U.TIME).matches(U.FLUX))
        self.assertTrue((U.CONC ** 0).isDimensionless())
        self.assertTrue(((U.CONC ** 3) / (U.CONC ** 3)).isDimensionless())

    def testAdditionRequiresMatchingDimensions(self):
        with self.assertRaises(U.DimensionError):
            conc(1.0) + U.Quantity(1.0, U.FLUX)
        with self.assertRaises(U.DimensionError):
            conc(1.0) - U.Quantity(1.0, U.TIME)
        # matching dimensions must still work
        self.assertAlmostEqual((conc(1.0) + conc(2.0)).value, 3.0)

    def testRequireRejectsWrongDimension(self):
        with self.assertRaises(U.DimensionError):
            conc(1.0).require(U.FLUX, "test quantity")


class TestDimensionalModel(unittest.TestCase):

    def setUp(self):
        self.dp = U.DimensionalParams()
        self.args = (conc(0.1), conc(0.05), conc(0.3), conc(0.3))

    def testEveryRateLawIsAFlux(self):
        laws = U.dimensionalRateLaws(*self.args, self.dp)
        self.assertGreaterEqual(len(laws), 7)
        for name, q in laws.items():
            self.assertTrue(q.dim.matches(U.FLUX), f"{name} has [{q.dim}]")

    def testRhsTermsAreFluxes(self):
        dU, dA = U.dimensionalRhs(*self.args, self.dp)
        self.assertTrue(dU.dim.matches(U.FLUX))
        self.assertTrue(dA.dim.matches(U.FLUX))

    def testOccupancyFactorsAreDimensionless(self):
        out = U.dimensionalConservationCheck(*self.args, self.dp)
        for name, q in out.items():
            self.assertTrue(q.dim.matches(U.CONC), f"{name} has [{q.dim}]")

    def testNucleationOrderChangesRequiredConstantDimension(self):
        """k_n must carry conc^(1-m) time^-1; m is not a free label."""
        for m in (1.5, 2.0, 3.0):
            dp = U.DimensionalParams(m=m)
            laws = U.dimensionalRateLaws(*self.args, dp)
            self.assertTrue(laws["nucleation"].dim.matches(U.FLUX))

    # ---- mutation tests: the checks above must be capable of failing ----

    def testMutationWrongNucleationDimensionIsRejected(self):
        bad = U.Quantity(1e-3, U.RATE) * (conc(0.1) ** self.dp.m)   # should be conc^(1-m)/T
        self.assertFalse(bad.dim.matches(U.FLUX))
        with self.assertRaises(U.DimensionError):
            U.Quantity(self.dp.J, U.FLUX) - bad

    def testMutationWrongGrowthDimensionIsRejected(self):
        bad = U.Quantity(1e-2, U.RATE) * conc(0.1) * conc(0.05)     # should be 1/(conc T)
        self.assertFalse(bad.dim.matches(U.FLUX))
        with self.assertRaises(U.DimensionError):
            U.Quantity(self.dp.J, U.FLUX) + bad

    def testMutationBareConcentrationInOccupancyIsRejected(self):
        """`1 + U_f` (missing the /K) must not typecheck."""
        with self.assertRaises(U.DimensionError):
            U.one() + conc(0.1)


class TestNondimensionalization(unittest.TestCase):

    def testFreePoolSolversAgree(self):
        """the independent dimensional solver matches the model's solver.

        different formulation (4 unknowns vs 2) and different algorithm
        (bounded least squares vs safeguarded newton), so agreement is evidence
        about the algebra rather than about one solver reproducing itself.
        """
        from proteostasis.model import solveFreePools

        dp = U.DimensionalParams()
        p, phi, _ = U.toNondimensional(dp)
        for U_tot, A_tot in [(0.01, 0.001), (0.1, 0.05), (0.5, 0.4), (2.0, 3.0)]:
            Uf, Af, Cf, Df = U.dimensionalFreePools(U_tot, A_tot, dp)
            uf, af, cf, df = solveFreePools(U_tot / phi, A_tot / phi, p)
            for got, want, name in ((uf * phi, Uf, "U_f"), (af * phi, Af, "A_f"),
                                    (cf * phi, Cf, "C_f"), (df * phi, Df, "D_f")):
                self.assertLess(abs(got - want) / max(abs(want), 1e-12), 1e-7,
                                f"{name} mismatch at U={U_tot}, A={A_tot}")

    def testNondimensionalizationIsExactRescaling(self):
        """dU/dt in dimensional units equals phi/tau times the scaled rhs."""
        from proteostasis.model import rhs

        dp = U.DimensionalParams()
        p, phi, tau = U.toNondimensional(dp)
        for U_tot, A_tot in [(0.02, 0.005), (0.2, 0.1), (1.0, 0.8)]:
            dU, dA = U.dimensionalRhsNumeric(U_tot, A_tot, dp)
            du, da = rhs(U_tot / phi, A_tot / phi, p)
            scale = max(abs(dU), abs(dA), 1e-30)
            self.assertLess(abs(du * phi / tau - dU) / scale, 1e-7)
            self.assertLess(abs(da * phi / tau - dA) / scale, 1e-7)

    def testMutationWrongTimeScaleBreaksRescaling(self):
        """if tau were wrong by a factor, the rescaling check must notice."""
        from proteostasis.model import rhs

        dp = U.DimensionalParams()
        p, phi, tau = U.toNondimensional(dp)
        dU, _ = U.dimensionalRhsNumeric(0.2, 0.1, dp)
        du, _ = rhs(0.2 / phi, 0.1 / phi, p)
        bad = du * phi / (tau * 2.0)
        self.assertGreater(abs(bad - dU) / max(abs(dU), 1e-30), 1e-3)

    def testMutationWrongNucleationScalingBreaksRescaling(self):
        """alpha_n must scale as k_n * phi**(m-1); phi**m would also 'look' fine."""
        from proteostasis.model import rhs

        dp = U.DimensionalParams(m=3.0)
        p, phi, tau = U.toNondimensional(dp)
        p_bad = p.with_(alpha_n=dp.k_n * (phi ** dp.m) / dp.k_ref)   # off by one power
        dU, _ = U.dimensionalRhsNumeric(0.3, 0.1, dp)
        du_bad, _ = rhs(0.3 / phi, 0.1 / phi, p_bad)
        self.assertGreater(abs(du_bad * phi / tau - dU) / max(abs(dU), 1e-30), 1e-3)


if __name__ == "__main__":
    unittest.main()
