"""phase 3 growth dilution: the theorem survives, the removal ceiling does not.

all tests here are MODEL-LEVEL and run on a clean checkout -- dilution needs no
phase 1 artefacts. imports are set inline for the reason documented in
`test_fold_theorem.py`.
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
import dilution as D  # noqa: E402

STATES = [(0.05, 0.01), (0.3, 0.1), (1.0, 0.5)]
LAWS = [("constant", D.Growth(0.02)), ("feedback", D.Growth(0.2, 0.5))]


class TestReductionToFrozenModel(unittest.TestCase):
    """mu = 0 must reproduce the upstream model bit for bit.

    this is the analogue of the phase 2 T0 epsilon -> 0 test: it licenses reading
    the diluted results as an extension rather than a different model.
    """

    def testExactAtZeroDilution(self):
        p = M.Params().validate()
        g = D.Growth(mu0=0.0)
        for u, a in STATES:
            self.assertEqual(M.rhs(u, a, p), D.rhsDil(u, a, p, g))

    def testJacobianExactAtZeroDilution(self):
        p = M.Params().validate()
        g = D.Growth(mu0=0.0)
        for u, a in STATES:
            self.assertTrue(np.array_equal(M.jacobian(u, a, p),
                                           D.jacobianDil(u, a, p, g)))


class TestDilutedJacobian(unittest.TestCase):
    def testAnalyticMatchesNumerical(self):
        p = M.Params().validate()
        for name, g in LAWS:
            for u, a in STATES:
                diff = np.max(np.abs(D.jacobianDil(u, a, p, g)
                                     - D.numericalJacobianDil(u, a, p, g)))
                self.assertLess(diff, 1e-6, f"{name} law at u={u}, a={a}")

    def testConstantDilutionIsJacobianMinusMuIdentity(self):
        """for constant mu the diluted jacobian must be exactly J - mu.I."""
        p = M.Params().validate()
        mu = 0.037
        g = D.Growth(mu0=mu)
        for u, a in STATES:
            expected = M.jacobian(u, a, p) - mu * np.eye(2)
            self.assertTrue(np.allclose(D.jacobianDil(u, a, p, g), expected,
                                        rtol=0.0, atol=1e-15))


class TestTheoremSurvivesDilution(unittest.TestCase):
    """det J_dil == -(grad R_tot x grad G_dil), for any dilution law."""

    def testIdentityHolds(self):
        p = M.Params().validate()
        for name, g in LAWS:
            for u, a in STATES:
                d = D.determinantIdentityDil(u, a, p, g)
                self.assertLess(d["rel_err"], 1e-5, f"{name} law: {d}")

    def testMassBalanceIncludesDilution(self):
        """du/dt + da/dt = j - R_tot exactly; the transfer still cancels."""
        p = M.Params().validate()
        for name, g in LAWS:
            for u, a in STATES:
                du, da = D.rhsDil(u, a, p, g)
                expected = p.j - D.removalTotal(u, a, p, g)
                scale = max(abs(expected), abs(du + da), 1.0)
                self.assertLess(abs((du + da) - expected) / scale, 1e-12, name)

    def testInfluxStillEntersOnlyTheSolubleEquation(self):
        """the structural fact the whole theorem rests on, under dilution."""
        base = M.Params()
        for name, g in LAWS:
            for u, a in STATES:
                ref = D.aggregateGDil(u, a, base.with_(j=0.0).validate(), g)
                for j in (0.01, 0.5):
                    self.assertEqual(
                        D.aggregateGDil(u, a, base.with_(j=j).validate(), g), ref,
                        f"{name} law: da/dt moved with j")


class TestSolvedFoldsAreGenuine(unittest.TestCase):
    def testSolvedFoldSatisfiesBothEquations(self):
        p = M.Params().validate()
        for name, g in (("no dilution", D.Growth(0.0)),
                        ("constant mu=0.02", D.Growth(0.02)),
                        ("feedback mu0=0.3", D.Growth(0.3, 0.5))):
            out = D.foldSolveDil(p, g)
            self.assertIsNotNone(out, f"{name}: no fold found")
            j, u, a = out
            self.assertLess(abs(D.aggregateGDil(u, a, p, g)), 1e-8, name)
            self.assertLess(abs(float(np.linalg.det(D.jacobianDil(u, a, p, g)))),
                            1e-7, name)
            self.assertAlmostEqual(j, D.removalTotal(u, a, p, g), places=12)

    def testZeroDilutionAgreesWithTheUndilutedSolver(self):
        import fold_theorem as FT
        p = M.Params().validate()
        a_ = FT.foldSolve(p)
        b_ = D.foldSolveDil(p, D.Growth(0.0))
        self.assertIsNotNone(a_)
        self.assertIsNotNone(b_)
        self.assertLess(abs(a_[0] - b_[0]) / a_[0], 1e-6)


class TestDilutionRaisesAndThenAbolishesTheBoundary(unittest.TestCase):
    """the substantive result: constant dilution destroys the fold."""

    def testCriticalInfluxRisesWithDilution(self):
        p = M.Params().validate()
        js = []
        for mu in (0.0, 0.02, 0.04, 0.06):
            out = D.foldSolveDil(p, D.Growth(mu0=mu))
            self.assertIsNotNone(out, f"fold lost early at mu={mu}")
            js.append(out[0])
        self.assertTrue(all(x < y for x, y in zip(js, js[1:])),
                        f"j_crit not increasing with dilution: {js}")

    def testConstantDilutionEventuallyAbolishesTheFold(self):
        p = M.Params().validate()
        self.assertIsNone(D.foldSolveDil(p, D.Growth(mu0=0.3)),
                          "a fold survived unbounded constant dilution")

    def testGrowthFeedbackRestoresTheFold(self):
        """same rate, but growth slows with burden -> the boundary comes back."""
        p = M.Params().validate()
        self.assertIsNone(D.foldSolveDil(p, D.Growth(mu0=0.3)))
        self.assertIsNotNone(D.foldSolveDil(p, D.Growth(mu0=0.3, k_mu=0.5)))

    def testCeilingIsNoLongerAnUpperBoundOnCriticalInflux(self):
        """A8 bounds enzymatic removal only; dilution is unbounded in burden."""
        p = M.Params().validate()
        u, a = 3.0, 4.0
        g = D.Growth(mu0=0.5)
        self.assertGreater(D.removalTotal(u, a, p, g), M.removalCeiling(p))


if __name__ == "__main__":
    unittest.main()
