"""phase 3: the discriminating prediction -- harmless folding load moves the boundary.

the exponent (D021) is robust but generic to saddle-nodes, so it cannot select
this model. this prediction can: it follows from rescue capacity being SHARED
between ordinary nascent chains and damaged protein, which an independent-
handling model does not have.

these tests pin the mechanism, the direction, and the confound requirement.
"""

import sys
import unittest
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (_REPO_ROOT / "scripts", _REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import dilution as D  # noqa: E402
import calibration as C  # noqa: E402
import gate4_discriminating as G  # noqa: E402


class TestMechanism(unittest.TestCase):
    """nu must consume capacity WITHOUT contributing damage -- both halves."""

    def testNascentLoadLowersFreeChaperone(self):
        base = M.Params().validate()
        prev = None
        for nu in (0.0, 0.5, 2.0, 5.0):
            _, _, cf, _ = M.solveFreePools(0.2, 0.05, base.with_(nu=nu).validate())
            if prev is not None:
                self.assertLess(cf, prev, f"free chaperone did not fall at nu={nu}")
            prev = cf

    def testNascentLoadAddsNoDamageInflux(self):
        """if nu contributed damage the perturbation would not be orthogonal."""
        base = M.Params().validate()
        for nu in (0.0, 1.0, 5.0):
            self.assertEqual(base.with_(nu=nu).validate().j, base.j)


class TestDirectionIsRobust(unittest.TestCase):
    def testCriticalInfluxFallsWithNascentLoad(self):
        for g in (D.Growth(0.0), C.calibratedGrowth(0.05, 0.02)):
            r = G.criticalInfluxVsNascent(M.Params().validate(), g)
            self.assertIsNotNone(r, f"ladder did not converge for {g}")
            self.assertTrue(r["any_drop"], "j_crit did not fall with nascent load")
            self.assertTrue(r["monotone_down"], "shift was not monotone")
            self.assertGreater(r["fold_drop"], 1.0)

    def testHoldsAcrossDifferentKinetics(self):
        base = M.Params()
        variants = [base.with_(alpha_n=0.3), base.with_(kappa_cu=0.5),
                    base.with_(rho_A=0.8), base.with_(m=2.5)]
        seen = 0
        for p in variants:
            r = G.criticalInfluxVsNascent(p.validate(), D.Growth(0.0))
            if r is None:
                continue
            seen += 1
            self.assertTrue(r["any_drop"], f"no drop for {p}")
        self.assertGreaterEqual(seen, 3, "too few variants converged to judge")


class TestLadderChoice(unittest.TestCase):
    def testLadderIsAscendingAndSpansAboutHundredFold(self):
        lad = G.NU_LADDER
        self.assertTrue(all(x < y for x, y in zip(lad, lad[1:])), lad)
        span = lad[-1] / lad[0]
        self.assertGreater(span, 50.0)
        self.assertLess(span, 200.0)

    def testWiderLadderGivesAtLeastAsMuchEffect(self):
        p = M.Params().validate()
        g = D.Growth(0.0)
        narrow = G.criticalInfluxVsNascent(p, g, (0.05, 0.15, 0.5, 1.5))
        wide = G.criticalInfluxVsNascent(p, g, G.NU_LADDER)
        self.assertIsNotNone(narrow)
        self.assertIsNotNone(wide)
        self.assertGreaterEqual(wide["fold_drop"], narrow["fold_drop"] - 1e-9)


class TestConfoundIsRecorded(unittest.TestCase):
    """the design is invalid in batch culture; that must not be lost in editing."""

    def testModuleDeclaresFixedGrowthRequirement(self):
        self.assertTrue(G.REQUIRES_FIXED_GROWTH_RATE)

    def testProposalStatesTheBatchCultureInvalidity(self):
        doc = (_REPO_ROOT / "empirical" / "GATE4_PROPOSAL.md").read_text()
        self.assertIn("chemostat", doc.lower())
        self.assertIn("K5", doc)
        self.assertIn("K6", doc)
        self.assertRegex(doc, r"batch version of this experiment\s+cannot test")


if __name__ == "__main__":
    unittest.main()
