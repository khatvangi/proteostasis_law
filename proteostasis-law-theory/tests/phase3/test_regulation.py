"""phase 3: the theorem beyond two states, and the regulated model.

the load-bearing test here is that det J = -det[grad R; grad G; grad C] holds in
THREE states. if it does, the fold theorem is a structural property rather than
an artefact of the two-state reduction, which is the main objection to it.
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
import regulation as R  # noqa: E402

STATES = [(0.05, 0.01), (0.3, 0.1), (1.0, 0.5)]


class TestReductionToFrozenModel(unittest.TestCase):
    """sigma0 = 0 must reproduce the two-state model bit for bit."""

    def testFieldIsIdentical(self):
        p = M.Params().validate()
        reg = R.Regulator(sigma0=0.0)
        for u, a in STATES:
            du, da = M.rhs(u, a, p)
            du3, da3, dc3 = R.rhs3(u, a, p.c_tot, p, reg)
            self.assertEqual(du, du3)
            self.assertEqual(da, da3)
            self.assertEqual(dc3, 0.0, "unregulated chaperone pool must not move")

    def testRemovalIsIdentical(self):
        p = M.Params().validate()
        for u, a in STATES:
            f = M.fluxes(u, a, p)
            self.assertAlmostEqual(
                R.removalR3(u, a, p.c_tot, p),
                f["refold"] + f["degrade_u"] + f["degrade_a"], places=14)


class TestGeneralisedIdentity(unittest.TestCase):
    """THE result: the theorem is not about two states."""

    def testHoldsInThreeStates(self):
        p = M.Params().validate()
        for reg in (R.Regulator(0.0), R.Regulator(0.6, 0.1, 0.5),
                    R.Regulator(1.5, 0.05, 0.8)):
            for u, a in STATES:
                d = R.determinantIdentity3(u, a, p.c_tot, p, reg)
                self.assertLess(d["rel_err"], 1e-5, f"{reg} at u={u}, a={a}: {d}")

    def testHoldsAwayFromTheNominalChaperonePool(self):
        """c is a STATE now, so the identity must hold off its initial value."""
        p = M.Params().validate()
        reg = R.Regulator(0.6, 0.1, 0.5)
        for c in (0.3, 0.6, 1.2):
            d = R.determinantIdentity3(0.3, 0.1, c, p, reg)
            self.assertLess(d["rel_err"], 1e-5, f"c={c}: {d}")


class TestRegulator(unittest.TestCase):
    def testSynthesisRisesAsFreeChaperoneFalls(self):
        """the sigma-32 titration mechanism, not a burden-sensing phenomenology."""
        reg = R.Regulator(0.6, 0.1, 0.5)
        s = [reg.synthesis(cf) for cf in (0.01, 0.05, 0.2, 1.0)]
        self.assertTrue(all(x > y for x, y in zip(s, s[1:])), s)

    def testDisabledRegulatorSynthesisesNothing(self):
        reg = R.Regulator(sigma0=0.0)
        for cf in (0.0, 0.1, 1.0):
            self.assertEqual(reg.synthesis(cf), 0.0)

    def testSteadyStatePoolIsPositiveAndFinite(self):
        reg = R.Regulator(0.6, 0.1, 0.5)
        # dc/dt = 0 at c = synthesis(cf)/delta
        for cf in (0.05, 0.3):
            c = reg.synthesis(cf) / reg.delta
            self.assertGreater(c, 0.0)
            self.assertTrue(np.isfinite(c))
            self.assertAlmostEqual(reg.dcdt(c, cf), 0.0, places=14)


class TestRegulatedFoldSolves(unittest.TestCase):
    def testSolvedFoldSatisfiesAllThreeEquations(self):
        p = M.Params().validate()
        reg = R.Regulator(0.6, 0.1, 0.5)
        out = R.foldSolve3(p, reg, seed=(0.3, 0.1, p.c_tot))
        if out is None:
            self.skipTest("no regulated fold at the base parameters")
        j, u, a, c = out
        _, da, dc = R.rhs3(u, a, c, p, reg)
        self.assertLess(abs(da), 1e-6)
        self.assertLess(abs(dc), 1e-6)
        det = float(np.linalg.det(R.jacobian3(u, a, c, p, reg)))
        self.assertLess(abs(det), 1e-5)
        self.assertAlmostEqual(j, R.removalR3(u, a, c, p), places=10)


if __name__ == "__main__":
    unittest.main()
