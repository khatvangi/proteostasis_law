"""phase 3 fold theorem: the identity, the solver, and the phi decomposition.

two kinds of test live here and they fail for different reasons.

MODEL-LEVEL tests pin the theorem itself -- mass balance, the j-independence of
the aggregate nullcline, and the determinant identity `det J = -(grad R x grad G)`
at ARBITRARY states rather than only at folds. these run on a clean checkout and
are the ones that actually pin the mathematics.

ARTEFACT tests assert the recorded phase 1 numbers. they skip when the gitignored
run root is absent, which is the same convention `tests/phase2` uses.
"""

NOTE_ON_IMPORTS = """these tests deliberately do NOT use a `_context.py` shim.

`tests/_context.py` and `tests/phase2/_context.py` already share a module name,
and pytest prepends each test file's own directory to sys.path, so a third one
here would be resolved by collection order rather than by location. the paths
are set explicitly instead, and the run-root helper is taken from the phase 3
module itself, which already owns it.
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
import fold_theorem as FT  # noqa: E402


# a spread of states that are all solvable, deliberately not near any fold
STATES = [(0.01, 0.001), (0.05, 0.01), (0.2, 0.05), (0.5, 0.2),
          (1.0, 0.5), (2.0, 1.0), (0.3, 1.5)]


class TestMassBalance(unittest.TestCase):
    """du/dt + da/dt = j - R exactly; the internal transfer must cancel."""

    def testTransferCancels(self):
        p = M.Params().validate()
        for u, a in STATES:
            du, da = M.rhs(u, a, p)
            expected = p.j - FT.removalR(u, a, p)
            scale = max(abs(expected), abs(du + da), 1.0)
            self.assertLess(abs((du + da) - expected) / scale, 1e-12,
                            f"mass balance failed at u={u}, a={a}")


class TestNullclineIsLoadIndependent(unittest.TestCase):
    """the aggregate nullcline must not move when j moves.

    this is the structural fact the whole theorem rests on: j enters du/dt and
    nowhere else, so {G = 0} is a fixed curve.
    """

    def testGIndependentOfJ(self):
        base = M.Params()
        for u, a in STATES:
            g_ref = FT.aggregateG(u, a, base.with_(j=0.0).validate())
            for j in (0.01, 0.1, 0.5, 1.0):
                g = FT.aggregateG(u, a, base.with_(j=j).validate())
                self.assertEqual(g, g_ref,
                                 f"da/dt changed with j at u={u}, a={a}")

    def testRIndependentOfJ(self):
        """total removal is also j-free; only the balance j - R depends on j."""
        base = M.Params()
        for u, a in STATES:
            r_ref = FT.removalR(u, a, base.with_(j=0.0).validate())
            for j in (0.01, 0.5):
                self.assertEqual(FT.removalR(u, a, base.with_(j=j).validate()), r_ref)


class TestDeterminantIdentity(unittest.TestCase):
    """det J = -(grad R x grad G), at arbitrary states, not only at folds.

    the residual floor is the central-difference error in the gradients, not the
    identity, so the tolerance is set accordingly.
    """

    def testIdentityAtArbitraryStates(self):
        p = M.Params().validate()
        for u, a in STATES:
            d = FT.determinantIdentity(u, a, p)
            self.assertLess(d["rel_err"], 1e-4,
                            f"identity failed at u={u}, a={a}: {d}")

    def testIdentityAcrossParameterVariation(self):
        base = M.Params()
        for nu in (0.0, 0.5, 5.0):
            for chi in (0.2, 0.5, 0.8):
                p = M.allocationParams(base.with_(nu=nu), chi).validate()
                d = FT.determinantIdentity(0.3, 0.1, p)
                self.assertLess(d["rel_err"], 1e-4,
                                f"identity failed at nu={nu}, chi={chi}: {d}")


class TestFoldSolver(unittest.TestCase):
    """the solved fold must satisfy both defining equations."""

    def testSolvedFoldSatisfiesBothEquations(self):
        p = M.Params().validate()
        out = FT.foldSolve(p)
        self.assertIsNotNone(out, "no fold found for the base parameters")
        j_crit, u_s, a_s = out

        self.assertLess(abs(FT.aggregateG(u_s, a_s, p)), 1e-8,
                        "solved fold is not on the aggregate nullcline")
        detJ = float(np.linalg.det(M.jacobian(u_s, a_s, p)))
        self.assertLess(abs(detJ), 1e-7, "solved fold is not a saddle-node")

        # and the critical influx is the removal flux there
        self.assertAlmostEqual(j_crit, FT.removalR(u_s, a_s, p), places=12)

    def testFoldLiesBelowTheRemovalCeiling(self):
        """C2 restated as a consequence rather than as a sampled fraction."""
        p = M.Params().validate()
        out = FT.foldSolve(p)
        self.assertIsNotNone(out)
        self.assertLess(out[0], M.removalCeiling(p))


class TestPhiDecomposition(unittest.TestCase):
    """the counterfactual split must be arithmetically coherent."""

    def testPhiReconstructsRemoval(self):
        p = M.Params().validate()
        for u, a in STATES:
            d = FT.phiDecomposition(u, a, p)
            self.assertAlmostEqual(d["phi"] * M.removalCeiling(p),
                                   FT.removalR(u, a, p), places=12)

    def testCounterfactualsRelaxTheConstraint(self):
        """removing either deficit can only raise phi, never lower it."""
        p = M.Params().validate()
        for u, a in STATES:
            d = FT.phiDecomposition(u, a, p)
            self.assertGreaterEqual(d["phi_if_saturated"], d["phi"] - 1e-15)
            self.assertGreaterEqual(d["phi_if_full_pools"], d["phi"] - 1e-15)
            for k in ("s_ref", "s_u", "s_a", "cf_frac", "df_frac"):
                self.assertGreaterEqual(d[k], 0.0)
                self.assertLessEqual(d[k], 1.0 + 1e-12)


class TestAgainstPhase1Run(unittest.TestCase):
    """artefact-dependent: the recorded phase 1 numbers. skips without results/."""

    @classmethod
    def setUpClass(cls):
        run = FT.phase1RunDir()
        if not run.is_dir():
            raise unittest.SkipTest(f"phase 1 run root absent: {run}")
        cls.v = FT.verifyAgainstRun(run)

    def testIdentityHoldsAtRecordedFolds(self):
        self.assertLess(self.v["identity_rel_err_median"], 1e-5)

    def testParallelismResidualIsBracketTolerance(self):
        """if the residual were a real failure it would not track the eigenvalue."""
        self.assertGreater(self.v["parallelism_tracks_eigenvalue"], 0.95)

    def testSolverReproducesTheContinuationSweep(self):
        self.assertGreaterEqual(self.v["n_solver"], 15)
        self.assertLess(self.v["solver_rel_err_max"], 1e-5)

    def testPhiRebuildsFromFirstPrinciples(self):
        self.assertGreaterEqual(self.v["n_phi"], 2800)
        self.assertLess(self.v["phi_rebuild_err_median"], 1e-10)
        self.assertLess(self.v["phi_rebuild_err_max"], 1e-6)

    def testCollapseHappensFarBelowSaturation(self):
        """the central mechanistic claim: machinery is nowhere near V_max."""
        for k in ("s_ref_median", "s_u_median", "s_a_median"):
            self.assertLess(self.v[k], 0.25, f"{k} is not far below saturation")
        self.assertLess(self.v["s_a_median"], self.v["s_ref_median"])

    def testSaturationDominatesSequestration(self):
        self.assertGreater(self.v["shortfall_share_saturation"],
                           self.v["shortfall_share_sequestration"])

    def testHeadlineValuesAsRecorded(self):
        self.assertAlmostEqual(self.v["phi_median"], 0.0769, places=3)
        self.assertAlmostEqual(self.v["shortfall_share_saturation"], 0.358, places=2)
        self.assertAlmostEqual(self.v["shortfall_share_sequestration"], 0.126, places=2)


if __name__ == "__main__":
    unittest.main()
