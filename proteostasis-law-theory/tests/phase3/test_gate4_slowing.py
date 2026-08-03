"""phase 3: the critical-slowing exponent that Gate 4's primary outcome rests on.

the saddle-node normal form predicts |lambda| ~ (j_crit - j)^(1/2) exactly, with
no dependence on any rate constant. these tests pin that the fitted exponent is
1/2 and that it survives dilution -- the prerequisite in GATE4_PROPOSAL.md 10.4.

model-level; only the sweep needs phase 1 artefacts and it skips without them.
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
import calibration as C  # noqa: E402
import gate4_slowing as S  # noqa: E402


class TestExponentAtBaseParameters(unittest.TestCase):
    """the headline: slope 1/2, and dilution does not move it."""

    def testUndilutedExponentIsOneHalf(self):
        r = S.slowingExponent(M.Params().validate(), D.Growth(0.0))
        self.assertIsNotNone(r, "no ladder converged without dilution")
        self.assertAlmostEqual(r["slope"], 0.5, delta=0.05)
        self.assertGreater(r["r2"], 0.99)

    def testDilutedExponentIsAlsoOneHalf(self):
        for g in (D.Growth(0.02), C.calibratedGrowth(0.05, 0.02),
                  C.calibratedGrowth(0.10, 0.02)):
            r = S.slowingExponent(M.Params().validate(), g)
            self.assertIsNotNone(r, f"no ladder converged for {g}")
            self.assertAlmostEqual(r["slope"], 0.5, delta=0.05)
            self.assertGreater(r["r2"], 0.99)

    def testRecoveryTimeExponentIsTheNegative(self):
        r = S.slowingExponent(M.Params().validate(), D.Growth(0.0))
        self.assertAlmostEqual(r["tau_exponent"], -r["slope"], places=12)


class TestEquilibriaAreGenuine(unittest.TestCase):
    def testReturnedStateIsAStableRoot(self):
        p = M.Params().validate()
        g = D.Growth(0.02)
        out = D.foldSolveDil(p, g)
        self.assertIsNotNone(out)
        j_crit, u_f, a_f = out
        j = j_crit * (1.0 - 1e-3)
        eq = S.stableEquilibrium(j, p, g, (u_f, a_f))
        self.assertIsNotNone(eq)
        u, a, lead = eq
        du, da = D.rhsDil(u, a, p.with_(j=j).validate(), g)
        self.assertLess(max(abs(du), abs(da)), 1e-8)
        self.assertLess(lead, 0.0, "must be the STABLE branch")

    def testEigenvalueShrinksTowardTheBoundary(self):
        """the substance of critical slowing, independent of the fitted slope."""
        p = M.Params().validate()
        g = D.Growth(0.02)
        j_crit, u_f, a_f = D.foldSolveDil(p, g)
        lams, guess = [], (u_f, a_f)
        for rel in (1e-2, 3e-3, 1e-3, 3e-4):
            eq = S.stableEquilibrium(j_crit * (1 - rel), p, g, guess)
            self.assertIsNotNone(eq, f"no equilibrium at rel={rel}")
            guess = (eq[0], eq[1])
            lams.append(abs(eq[2]))
        self.assertTrue(all(x > y for x, y in zip(lams, lams[1:])), lams)


class TestExponentIsParameterFree(unittest.TestCase):
    """1/2 must not depend on the rate constants -- that is the whole point."""

    def testHoldsAcrossDifferentKinetics(self):
        base = M.Params()
        variants = [base.with_(alpha_n=0.2), base.with_(m=2.5),
                    base.with_(kappa_ref=1.2), base.with_(rho_U=1.6)]
        got = 0
        for p in variants:
            r = S.slowingExponent(p.validate(), D.Growth(0.0))
            if r is None:
                continue
            got += 1
            self.assertAlmostEqual(r["slope"], 0.5, delta=0.08)
        self.assertGreaterEqual(got, 3, "too few variants converged to judge")


if __name__ == "__main__":
    unittest.main()
