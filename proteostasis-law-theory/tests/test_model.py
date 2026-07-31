"""model-level tests: conservation, positivity, jacobian correctness.

these are the checks that must fail on negative-state leakage, on a broken
conserved-pool closure, and on an incorrect jacobian. mutation tests confirm
each one has the power to fail.
"""

import unittest

import numpy as np

import _context  # noqa: F401
from proteostasis import (Params, ModelError, rhs, fluxes, jacobian, numericalJacobian,
                          solveFreePools, solveFreePoolsCertified, poolBalanceResiduals,
                          massBalanceResidual, removalCeiling, allocationParams)
from proteostasis.sweeps import GLOBAL_RANGES, latinHypercube, paramsFromSample

STATES = [(0.0, 0.0), (0.0, 0.5), (0.5, 0.0), (1e-9, 1e-9), (0.01, 0.001),
          (0.2, 0.05), (1.0, 1.0), (5.0, 0.3), (0.3, 5.0), (50.0, 80.0)]


def sampleParams(n=12, seed=7):
    base = Params()
    return [base] + [paramsFromSample(s, base)
                     for s in latinHypercube(GLOBAL_RANGES, n, seed)]


class TestParams(unittest.TestCase):

    def testRejectsNonPositivePools(self):
        for kw in ({"c_tot": 0.0}, {"d_tot": -1.0}, {"kappa_cu": 0.0}, {"alpha_g": -2.0}):
            with self.assertRaises(ModelError):
                Params(**kw).validate()

    def testRejectsNegativeLoads(self):
        with self.assertRaises(ModelError):
            Params(j=-1e-6).validate()
        with self.assertRaises(ModelError):
            Params(nu=-0.1).validate()

    def testRejectsNucleationOrderAtOrBelowOne(self):
        """m <= 1 removes the superlinear transfer the fold structure needs."""
        for m in (0.5, 1.0):
            with self.assertRaises(ModelError):
                Params(m=m).validate()

    def testAllocationSplitsAFixedRescueBudget(self):
        p = allocationParams(Params(), 0.25, 1.0)
        self.assertAlmostEqual(p.c_tot, 0.25)
        self.assertAlmostEqual(p.d_tot, 0.75)
        self.assertAlmostEqual(p.rescueTotal, 1.0)
        self.assertAlmostEqual(p.chaperoneShare, 0.25)
        for chi in (0.0, 1.0, -0.1, 1.2):
            with self.assertRaises(ModelError):
                allocationParams(Params(), chi, 1.0)


class TestConservedPools(unittest.TestCase):

    def testFreePoolSolutionIsUnique(self):
        """least and greatest fixed points of the monotone map must coincide.

        by knaster-tarski they bracket EVERY fixed point, so agreement is a
        uniqueness certificate rather than a convergence anecdote.
        """
        for p in sampleParams():
            for u, a in STATES:
                *_, cert = solveFreePoolsCertified(u, a, p)
                self.assertTrue(cert["unique"],
                                f"non-unique free pools at u={u}, a={a}: gap={cert['gap']:.3e}")

    def testClosureReproducesTotals(self):
        for p in sampleParams():
            for u, a in STATES:
                res = poolBalanceResiduals(u, a, p)
                for name, v in res.items():
                    self.assertLess(abs(v), 1e-10, f"{name} residual {v:.3e}")

    def testFreePoolsNeverExceedTotals(self):
        """free resource cannot exceed the conserved pool; free substrate cannot
        exceed the total substrate. a violation means occupancy was double
        counted or dropped."""
        for p in sampleParams():
            for u, a in STATES:
                uf, af, cf, df = solveFreePools(u, a, p)
                self.assertGreaterEqual(cf, 0.0)
                self.assertGreaterEqual(df, 0.0)
                self.assertLessEqual(cf, p.c_tot * (1 + 1e-12))
                self.assertLessEqual(df, p.d_tot * (1 + 1e-12))
                self.assertLessEqual(uf, u * (1 + 1e-12))
                self.assertLessEqual(af, a * (1 + 1e-12))

    def testNascentLoadConsumesCapacityWithoutAddingInflux(self):
        """the mechanism the theory requires, stated as three separate facts.

        ordinary nascent-chain folding is NOT damage influx, yet it must be
        able to move the stability boundary. so: (a) the influx term is
        literally untouched by nu, (b) free chaperone falls, (c) refolding flux
        falls and the soluble burden accumulates faster. if nu ever leaked into
        the influx, (a) would fail; if nu were decoupled from capacity, (b) and
        (c) would fail.
        """
        u, a = 0.2, 0.05
        p0 = Params(nu=0.0)
        p1 = p0.with_(nu=3.0)
        f0, f1 = fluxes(u, a, p0), fluxes(u, a, p1)

        self.assertEqual(f0["influx"], f1["influx"])        # (a) no damage added
        self.assertEqual(f0["influx"], p0.j)
        self.assertLess(f1["cf"], f0["cf"])                 # (b) capacity consumed
        self.assertLess(f1["refold"], f0["refold"])         # (c) rescue weakened
        self.assertGreater(rhs(u, a, p1)[0], rhs(u, a, p0)[0])

    def testNascentLoadIsMonotoneInCapacityConsumption(self):
        u, a = 0.2, 0.05
        cfs = [solveFreePools(u, a, Params(nu=nu))[2] for nu in (0.0, 0.5, 2.0, 8.0)]
        self.assertTrue(all(b < c for c, b in zip(cfs, cfs[1:])),
                        f"free chaperone not monotone decreasing in nu: {cfs}")

    def testRejectsNegativeTotals(self):
        with self.assertRaises(ModelError):
            solveFreePools(-1e-9, 0.0, Params())


class TestMassBalance(unittest.TestCase):

    def testTransferTermsCancel(self):
        for p in sampleParams():
            for u, a in STATES:
                self.assertLess(abs(massBalanceResidual(u, a, p)), 1e-12,
                                f"mass balance broken at u={u}, a={a}")

    def testMutationBrokenStoichiometryIsDetected(self):
        """if nucleation entered A with the wrong coefficient, the check fires."""
        p = Params()
        u, a = 0.3, 0.1
        f = fluxes(u, a, p)
        du = (f["influx"] - f["refold"] - f["degrade_u"] - f["nucleate"]
              - f["grow"] + f["disaggregate"])
        da_bad = 0.9 * f["nucleate"] + f["grow"] - f["disaggregate"] - f["degrade_a"]
        expected = f["influx"] - f["refold"] - f["degrade_u"] - f["degrade_a"]
        scale = max(abs(v) for k, v in f.items()
                    if k in ("influx", "refold", "degrade_u", "degrade_a",
                             "nucleate", "grow", "disaggregate"))
        self.assertGreater(abs((du + da_bad) - expected) / scale, 1e-6)


class TestPositivity(unittest.TestCase):

    def testBoundaryFluxesPointInward(self):
        """forward invariance of the nonnegative orthant, stated exactly.

        at u=0 the rhs must not push u negative; at a=0 likewise. this is the
        algebraic statement whose integrated consequence is 'no negative-state
        leakage', and it fails immediately if a removal-term sign is flipped.
        """
        for p in sampleParams(n=20, seed=11):
            for a in (0.0, 1e-6, 0.1, 2.0):
                du, _ = rhs(0.0, a, p)
                self.assertGreaterEqual(du, -1e-15, f"du/dt<0 at u=0, a={a}")
            for u in (0.0, 1e-6, 0.1, 2.0):
                _, da = rhs(u, 0.0, p)
                self.assertGreaterEqual(da, -1e-15, f"da/dt<0 at a=0, u={u}")

    def testMutationFlippedSignBreaksInvariance(self):
        """a removal term wrongly applied at the boundary must be caught."""
        p = Params(j=0.0)
        f = fluxes(0.0, 0.5, p)
        du_ok = (f["influx"] - f["refold"] - f["degrade_u"] - f["nucleate"]
                 - f["grow"] + f["disaggregate"])
        du_bad = du_ok - 2.0 * f["disaggregate"]      # disaggregation sign flipped
        self.assertGreaterEqual(du_ok, -1e-15)
        self.assertLess(du_bad, -1e-15)


class TestJacobian(unittest.TestCase):

    def testAnalyticMatchesFiniteDifference(self):
        for p in sampleParams(n=16, seed=3):
            for u, a in STATES:
                J = jacobian(u, a, p)
                Jn = numericalJacobian(u, a, p)
                denom = max(float(np.max(np.abs(Jn))), 1e-12)
                err = float(np.max(np.abs(J - Jn))) / denom
                self.assertLess(err, 1e-5, f"jacobian error {err:.3e} at u={u}, a={a}")

    def testJacobianAccountsForResourceSequestration(self):
        """the implicit free-pool dependence must contribute.

        a jacobian that treated cf and df as constants would be wrong. that
        'frozen-resource' jacobian is computed here and must differ measurably
        from the true one wherever sequestration is significant.
        """
        p = Params()
        u, a = 0.8, 0.4
        uf, af, cf, df = solveFreePools(u, a, p)
        self.assertLess(cf, p.c_tot * 0.9, "test state must actually sequester chaperone")
        J = jacobian(u, a, p)
        frozen = _frozenResourceJacobian(u, a, p, cf, df)
        self.assertGreater(float(np.max(np.abs(J - frozen))) / float(np.max(np.abs(J))), 1e-3)

    def testMutationPerturbedJacobianIsRejected(self):
        """confirm the comparison in testAnalyticMatchesFiniteDifference bites."""
        p = Params()
        u, a = 0.2, 0.05
        J = jacobian(u, a, p)
        Jn = numericalJacobian(u, a, p)
        J_bad = J.copy()
        J_bad[0, 1] *= 1.01
        denom = max(float(np.max(np.abs(Jn))), 1e-12)
        self.assertGreater(float(np.max(np.abs(J_bad - Jn))) / denom, 1e-5)


def _frozenResourceJacobian(u, a, p, cf, df):
    """finite-difference jacobian with the free resources held fixed."""
    def f(uu, aa):
        su = 1.0 + cf / p.kappa_cu + df / p.kappa_du
        sa = 1.0 + cf / p.kappa_ca + df / p.kappa_da
        uf, af = uu / su, aa / sa
        du = (p.j - cf * uf / (p.kappa_ref + uf) - p.rho_U * df * uf / (p.kappa_u + uf)
              - p.alpha_n * uf ** p.m - p.alpha_g * uf * af
              + p.alpha_d * cf * af / (p.kappa_dis + af))
        da = (p.alpha_n * uf ** p.m + p.alpha_g * uf * af
              - p.alpha_d * cf * af / (p.kappa_dis + af)
              - p.rho_A * df * af / (p.kappa_a + af))
        return np.array([du, da])
    h = 1e-6
    return np.column_stack([(f(u + h, a) - f(u - h, a)) / (2 * h),
                            (f(u, a + h) - f(u, a - h)) / (2 * h)])


class TestAnalyticBounds(unittest.TestCase):

    def testRemovalCeilingBoundsTotalRemoval(self):
        """no state may remove more than c_tot + (rho_U+rho_A) d_tot."""
        for p in sampleParams(n=16, seed=5):
            ceiling = removalCeiling(p)
            for u, a in STATES:
                f = fluxes(u, a, p)
                removal = f["refold"] + f["degrade_u"] + f["degrade_a"]
                self.assertLessEqual(removal, ceiling * (1 + 1e-12),
                                     f"removal {removal:.6g} exceeds ceiling {ceiling:.6g}")

    def testAboveCeilingBurdenAlwaysGrows(self):
        for p in sampleParams(n=10, seed=13):
            pv = p.with_(j=1.5 * removalCeiling(p))
            for u, a in STATES:
                du, da = rhs(u, a, pv)
                self.assertGreater(du + da, 0.0,
                                   f"burden not growing above the removal ceiling at u={u}, a={a}")


if __name__ == "__main__":
    unittest.main()
