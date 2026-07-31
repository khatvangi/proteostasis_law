"""the analytic continuation must be invisible where it matters.

`scripts/phase2/boron_continuation.py` exists only so that nitrogen's
linear-coordinate root protocol can be run against boron's model at all.  the
claim that this does not compromise cell 1 rests on two properties, both
asserted here:

  1. on the nonnegative orthant the continuation DELEGATES -- every value is
     bit-identical to the shipped `proteostasis.model`, so the accepted root,
     its residual, its jacobian and its eigenvalues come from the phase 1 code
     path unchanged;
  2. past the boundary the continuation is a genuine solution of the same
     binding residual, and its analytic jacobian matches central differences of
     its own right-hand side.

if (1) ever fails, the benchmark is silently running a different model.
"""

import unittest

import numpy as np

import _context  # noqa: F401
from phase2 import boron_continuation as bc
from phase2.lhs import parametersForEpsilon, sampleMatrix
from phase2.mapping import nitrogenToBoron
from proteostasis import model as bm

MATRIX = sampleMatrix(16, 20260731)
NONNEGATIVE = [(0.0, 0.0), (1e-8, 1e-8), (0.01, 0.001), (0.3, 0.2),
               (1.0, 1.0), (7.0, 3.0), (100.0, 100.0)]
NEGATIVE = [(0.001, -1e-4), (-1e-4, 0.002), (-1e-3, -1e-3),
            (0.5, -0.05), (2.0, -0.7)]


def _params(i: int, eps: float):
    return bm.Params(**nitrogenToBoron(parametersForEpsilon(MATRIX[i], eps), eps)).validate()


class TestDelegationIsExact(unittest.TestCase):

    def testFreePoolsBitIdentical(self):
        for eps in (1e-6, 1e-2, 0.3, 2.0):
            for i in (0, 5, 15):
                p = _params(i, eps)
                for u, a in NONNEGATIVE:
                    np.testing.assert_array_equal(
                        np.array(bc.solveFreePoolsExtended(u, a, p)),
                        np.array(bm.solveFreePools(u, a, p)))

    def testRhsBitIdentical(self):
        for eps in (1e-6, 0.3, 2.0):
            for i in (1, 9):
                p = _params(i, eps)
                for u, a in NONNEGATIVE:
                    np.testing.assert_array_equal(
                        bc.rhsVectorExtended((u, a), p), bm.rhsVector((u, a), p))

    def testJacobianBitIdentical(self):
        for eps in (1e-6, 0.3, 2.0):
            for i in (2, 11):
                p = _params(i, eps)
                for u, a in NONNEGATIVE:
                    np.testing.assert_array_equal(
                        bc.jacobianExtended(u, a, p), bm.jacobian(u, a, p))

    def testContinuationFlagOnlyFiresOffTheOrthant(self):
        for u, a in NONNEGATIVE:
            self.assertFalse(bc.continuationUsed(u, a))
        for u, a in NEGATIVE:
            self.assertTrue(bc.continuationUsed(u, a))


class TestContinuationIsCorrect(unittest.TestCase):

    def testBindingResidualIsSatisfied(self):
        """the continued (cf, df) really solves the same pool balances."""
        for eps in (1e-6, 0.1, 2.0):
            for i in (0, 7, 14):
                p = _params(i, eps)
                for u, a in NEGATIVE:
                    uf, af, cf, df = bc.solveFreePoolsExtended(u, a, p)
                    (r1, r2), _, _ = bm._bindingResidual(cf, df, u, a, p)
                    scale = max(1.0, p.c_tot, p.d_tot)
                    self.assertLess(max(abs(r1), abs(r2)) / scale, 1e-8,
                                    f"eps={eps} i={i} state=({u},{a})")

    def testAnalyticJacobianMatchesCentralDifference(self):
        for eps in (1e-6, 0.1, 2.0):
            for i in (3, 12):
                p = _params(i, eps)
                for u, a in NEGATIVE:
                    got = bc.jacobianExtended(u, a, p)
                    num = np.zeros((2, 2))
                    x = np.array([u, a], dtype=float)
                    for k in range(2):
                        h = 1e-7 * max(abs(x[k]), 1e-3)
                        xp, xm = x.copy(), x.copy()
                        xp[k] += h
                        xm[k] -= h
                        num[:, k] = (bc.rhsVectorExtended(xp, p)
                                     - bc.rhsVectorExtended(xm, p)) / (2 * h)
                    scale = max(float(np.max(np.abs(got))), 1e-30)
                    self.assertLess(float(np.max(np.abs(got - num))) / scale, 1e-5,
                                    f"eps={eps} i={i} state=({u},{a})")

    def testFieldIsContinuousAcrossTheBoundary(self):
        """approaching a = 0 from below reproduces the value at a = 0."""
        for eps in (1e-6, 0.3, 2.0):
            p = _params(4, eps)
            at_zero = bc.rhsVectorExtended((0.4, 0.0), p)
            approach = bc.rhsVectorExtended((0.4, -1e-10), p)
            np.testing.assert_allclose(approach, at_zero, rtol=1e-7, atol=1e-12)

    def testNonIntegerNucleationOrderIsRefusedNotFudged(self):
        """uf < 0 with fractional m would make alpha_n*uf**m complex.

        the benchmark pins m = 2.0 so this never fires in practice; the guard is
        here so a future non-integer m fails loudly instead of returning nan.
        note it is uf, not af, that is raised to a power -- a state with a < 0
        but u > 0 is perfectly fine.
        """
        p = _params(0, 1.0).with_(m=2.5)
        self.assertIsNotNone(bc.fluxesExtended(0.1, -0.01, p))     # af < 0 is fine
        with self.assertRaises(bm.ModelError):
            bc.fluxesExtended(-0.01, 0.1, p)
        with self.assertRaises(bm.ModelError):
            bc.jacobianExtended(-0.01, 0.1, p)


class TestClampWasRejectedForCause(unittest.TestCase):
    """the alternative fix, recorded so it is not silently reintroduced.

    clamping every evaluation at max(x, 0) was measured to change a large
    fraction of free-limit results, which is why the continuation exists.  this
    test re-measures on a small sample and fails if the clamp ever became
    harmless -- in which case the design note in `boron_continuation.py` would
    need revisiting rather than being trusted.
    """

    def testClampIsNotANoOp(self):
        from phase2 import protocols
        from phase2.models import FreeLimitAdapter

        class Clamped:
            def __init__(self, inner):
                self.i = inner

            def __getattr__(self, n):
                return getattr(self.i, n)

            def rhsVector(self, x):
                return self.i.rhsVector((max(float(x[0]), 0.0), max(float(x[1]), 0.0)))

            def jacobian(self, x):
                return self.i.jacobian((max(float(x[0]), 0.0), max(float(x[1]), 0.0)))

        differing = 0
        for i in range(16):
            q = parametersForEpsilon(MATRIX[i], 1e-6)
            a = FreeLimitAdapter(q, 1e-6)
            r0 = protocols.classifyNitrogen(a)
            r1 = protocols.classifyNitrogen(Clamped(a))
            same = (r0["label"] == r1["label"] and (r0["root"] is None) == (r1["root"] is None)
                    and (r0["root"] is None
                         or abs(r0["root"]["u"] - r1["root"]["u"]) < 1e-12))
            differing += 0 if same else 1
        self.assertGreater(differing, 0,
                           "clamping is now a no-op; revisit boron_continuation.py")


if __name__ == "__main__":
    unittest.main()
