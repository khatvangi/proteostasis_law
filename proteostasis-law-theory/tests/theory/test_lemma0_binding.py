"""Lemma 0, checked numerically over the kinetic box.

Task A2. The proof in theory/LEMMA0_BINDING.md does not depend on this file --
that is the point of proving it. What the checks add is a guard: if a future
edit breaks the structure the proof relies on (1:1 rapid equilibrium, additive
constant in each denominator, no cooperativity), the M-matrix property is the
first thing that goes, and it goes silently, because `solveFreePools` would
still return a number.

Three claims are checked, each an inequality the proof derives in closed form:

  * off-diagonal entries of the binding Jacobian are nonpositive (Z-matrix);
  * the two contractions against (C_f, D_f) exceed C_f(1+nu) and D_f;
  * det >= 1 + nu.

The uniqueness certificate is checked separately, from below and from above.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (_REPO_ROOT / "scripts", _REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import fold_theorem as FT  # noqa: E402


def _mMatrixMargins(u: float, a: float, p: M.Params):
    """(max off-diagonal, min row contraction margin, det - (1+nu))."""
    _, _, cf, df = M.solveFreePools(u, a, p)
    (_, _), (j11, j12, j21, j22), _ = M._bindingResidual(cf, df, u, a, p)
    det = j11 * j22 - j12 * j21
    row1 = (cf * j11 + df * j12) / cf - (1.0 + p.nu)
    row2 = (cf * j21 + df * j22) / df - 1.0
    return max(j12, j21), min(row1, row2), det - (1.0 + p.nu)


class TestBindingJacobianIsAnMMatrix(unittest.TestCase):

    def testOverAWideParameterAndStateSweep(self):
        rng = np.random.default_rng(20260805)
        base = M.Params().validate()
        off, row, det = -np.inf, np.inf, np.inf
        n = 0
        for _ in range(60):
            p = base.with_(
                nu=float(rng.uniform(0.0, 3.0)),
                c_tot=float(rng.uniform(0.1, 0.9)),
                d_tot=float(rng.uniform(0.1, 0.9)),
                kappa_cu=float(10 ** rng.uniform(-2, 1)),
                kappa_ca=float(10 ** rng.uniform(-2, 1)),
                kappa_du=float(10 ** rng.uniform(-2, 1)),
                kappa_da=float(10 ** rng.uniform(-2, 1)),
            ).validate()
            for u in np.geomspace(1e-5, 50.0, 8):
                for a in np.geomspace(1e-6, 50.0, 8):
                    o, r, d = _mMatrixMargins(u, a, p)
                    off, row, det = max(off, o), min(row, r), min(det, d)
                    n += 1
        self.assertEqual(n, 60 * 64)
        self.assertLessEqual(off, 0.0, "binding Jacobian is not a Z-matrix")
        self.assertGreater(row, -1e-12, "row contraction bound violated")
        self.assertGreater(det, -1e-12, "det J < 1 + nu; M-matrix bound violated")

    def testAtTheSolvedFoldStatesOfTheKineticBox(self):
        """the population the paper reports, not a convenience sample."""
        run = FT.phase1RunDir()
        if not (run / "C" / "samples.tsv").exists():
            self.skipTest("phase 1 run root absent (results/ is gitignored)")
        c = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
        c = c[c["C1_fold_exists"] == True]  # noqa: E712
        c = c[pd.to_numeric(c["fold_burden"], errors="coerce").notna()]
        off, row, det, n = -np.inf, np.inf, np.inf, 0
        for _, r in c.iterrows():
            try:
                p = FT.paramsFromSampleRow(r)
                u, a = FT.foldStateFromSampleRow(r)
                o, rw, d = _mMatrixMargins(u, a, p)
            except (M.ModelError, ValueError, KeyError):
                continue
            off, row, det = max(off, o), min(row, rw), min(det, d)
            n += 1
        self.assertGreater(n, 2800, f"only {n} states evaluated")
        self.assertLessEqual(off, 0.0)
        self.assertGreater(row, -1e-10)
        self.assertGreater(det, -1e-10)


class TestTheFixedPointIsUniqueAndPositive(unittest.TestCase):

    def testMonotoneIterationFromBothEndsAgrees(self):
        """Knaster-Tarski brackets every fixed point; agreement is uniqueness."""
        p = M.Params().validate()
        worst_gap = 0.0
        for u in np.geomspace(1e-5, 50.0, 12):
            for a in np.geomspace(1e-6, 50.0, 12):
                *_, cert = M.solveFreePoolsCertified(u, a, p)
                self.assertTrue(cert["unique"], f"gap {cert['gap']:.3e} at {u},{a}")
                worst_gap = max(worst_gap, cert["gap"])
        self.assertLess(worst_gap, 1e-9)

    def testTheSolutionIsStrictlyPositiveWithTheProvedFloor(self):
        p = M.Params().validate()
        for u, a in ((0.0, 0.0), (1e-8, 1e-8), (50.0, 50.0), (200.0, 0.0)):
            _, _, cf, df = M.solveFreePools(u, a, p)
            floor_c = p.c_tot / (1.0 + p.nu + u / p.kappa_cu + a / p.kappa_ca)
            floor_d = p.d_tot / (1.0 + u / p.kappa_du + a / p.kappa_da)
            self.assertGreater(cf, 0.0)
            self.assertGreater(df, 0.0)
            self.assertGreaterEqual(cf, floor_c * (1.0 - 1e-9))
            self.assertGreaterEqual(df, floor_d * (1.0 - 1e-9))


class TestTheMonotonicityFactsUsedInSectionThreeThree(unittest.TestCase):
    """d(u_f)/du > 0 and d(a_f)/du >= 0 -- the second breaks the old sign argument."""

    def testFreeMonomerRisesAndFreeAggregateAlsoRises(self):
        p = M.Params().validate()
        min_duf, min_daf = np.inf, np.inf
        for u in np.geomspace(1e-3, 20.0, 20):
            for a in np.geomspace(1e-4, 20.0, 20):
                h = 1e-6 * u
                ufp, afp, *_ = M.solveFreePools(u + h, a, p)
                ufm, afm, *_ = M.solveFreePools(u - h, a, p)
                min_duf = min(min_duf, (ufp - ufm) / (2 * h))
                min_daf = min(min_daf, (afp - afm) / (2 * h))
        self.assertGreater(min_duf, 0.0, "d(u_f)/du must be strictly positive")
        self.assertGreaterEqual(min_daf, -1e-12,
                                "d(a_f)/du must be nonnegative")


class TestTheGuardIsDocumentedAsUnreachable(unittest.TestCase):
    def testModelSaysTheSingularityCheckCannotFire(self):
        src = (_REPO_ROOT / "scripts" / "proteostasis" / "model.py").read_text()
        i = src.index("singular binding jacobian")
        self.assertIn("CANNOT FIRE", src[max(0, i - 600): i])


if __name__ == "__main__":
    unittest.main()
