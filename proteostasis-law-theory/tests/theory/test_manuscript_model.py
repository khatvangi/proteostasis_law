"""Section 2.1 must be the model that was solved, not a nearby one.

Task A1. Through v5 the manuscript printed rate laws in TOTAL concentrations,
`v_ref = k_ref C_f u/(K_ref + u)`, while `model.py` solves them in FREE
concentrations. A reader implementing the printed §2.1 would have got a
different model, and every quantitative result in Sections 6-8 is computed from
the coded one.

The check here is the one that would have caught it: §2.1 is implemented AGAIN,
from the printed equations alone, by a different algorithm -- plain damped
fixed-point iteration on the closure rather than the safeguarded Newton in
`model.py` -- and the two right-hand sides are compared on a grid. Independence
is the point. Importing anything from `model.py` other than `Params` and the
functions under test would make this a tautology, which is the failure mode
Section 5 warns about for the identity residual.
"""

from __future__ import annotations

import re
import sys
import unittest
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[2]
if str(_REPO_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT / "scripts"))

from proteostasis import model as M  # noqa: E402

_MS = (_REPO_ROOT / "manuscript" / "MANUSCRIPT_BMB_v5.md").read_text()


def _section21() -> str:
    start = _MS.index("### 2.1 State equations")
    return _MS[start: _MS.index("### 2.2 Hypotheses")]


# manuscript symbol -> the `Params` field carrying it. `k_ref` is 1 by the
# nondimensionalisation (time scale tau = 1/k_ref), so it has no field.
SYMBOL_MAP = {
    "j": "j", "ν": "nu", "m": "m",
    "C_tot": "c_tot", "D_tot": "d_tot",
    "k_U": "rho_U", "k_A": "rho_A",
    "k_n": "alpha_n", "k_g": "alpha_g", "k_dis": "alpha_d",
    "K_ref": "kappa_ref", "K_U": "kappa_u",
    "K_dis": "kappa_dis", "K_A": "kappa_a",
    "K_CU": "kappa_cu", "K_CA": "kappa_ca",
    "K_DU": "kappa_du", "K_DA": "kappa_da",
}


def freePools(u: float, a: float, p: M.Params, tol: float = 1e-15,
              maxiter: int = 200_000):
    """the four printed conservation laws, solved by naive fixed-point iteration.

    Deliberately the simplest possible solver: substitute, repeat. It converges
    because the map is monotone and bounded (Lemma 0), which is exactly the
    property the proof in theory/LEMMA0_BINDING.md establishes.
    """
    C_f, D_f = p.c_tot, p.d_tot
    for _ in range(maxiter):
        u_f = u / (1.0 + C_f / p.kappa_cu + D_f / p.kappa_du)
        a_f = a / (1.0 + C_f / p.kappa_ca + D_f / p.kappa_da)
        C_new = p.c_tot / (1.0 + p.nu + u_f / p.kappa_cu + a_f / p.kappa_ca)
        D_new = p.d_tot / (1.0 + u_f / p.kappa_du + a_f / p.kappa_da)
        if max(abs(C_new - C_f), abs(D_new - D_f)) <= tol:
            C_f, D_f = C_new, D_new
            break
        C_f, D_f = C_new, D_new
    u_f = u / (1.0 + C_f / p.kappa_cu + D_f / p.kappa_du)
    a_f = a / (1.0 + C_f / p.kappa_ca + D_f / p.kappa_da)
    return u_f, a_f, C_f, D_f


def manuscriptRhs(u: float, a: float, p: M.Params):
    """du/dt, da/dt from the equations as §2.1 prints them."""
    u_f, a_f, C_f, D_f = freePools(u, a, p)
    v_ref = 1.0 * C_f * u_f / (p.kappa_ref + u_f)      # k_ref = 1
    v_degU = p.rho_U * D_f * u_f / (p.kappa_u + u_f)
    n = p.alpha_n * u_f ** p.m
    g = p.alpha_g * u_f * a_f
    v_dis = p.alpha_d * C_f * a_f / (p.kappa_dis + a_f)
    v_degA = p.rho_A * D_f * a_f / (p.kappa_a + a_f)
    return (p.j - v_ref - v_degU - n - g + v_dis,
            n + g - v_dis - v_degA)


class TestSectionTwoOneIsTheCodedModel(unittest.TestCase):

    def testTheTwoImplementationsAgreeOnAGrid(self):
        """independent solver, independent transcription, same field."""
        p = M.Params().validate()
        worst = 0.0
        n = 0
        for u in np.geomspace(1e-4, 20.0, 26):
            for a in np.geomspace(1e-6, 20.0, 26):
                du_m, da_m = manuscriptRhs(u, a, p)
                du_c, da_c = M.rhs(u, a, p)
                scale = max(abs(du_c), abs(da_c), 1e-12)
                worst = max(worst, abs(du_m - du_c) / scale,
                            abs(da_m - da_c) / scale)
                n += 1
        self.assertEqual(n, 676)
        self.assertLess(worst, 1e-10,
                        f"§2.1 and model.py disagree by {worst:.3e} relative")

    def testAgreementSurvivesAwayFromTheBaseParameters(self):
        """one grid at one parameter set would not have caught total-vs-free."""
        rng = np.random.default_rng(20260805)
        base = M.Params().validate()
        worst = 0.0
        for _ in range(40):
            p = base.with_(
                nu=float(rng.uniform(0.0, 3.0)),
                c_tot=float(rng.uniform(0.2, 0.8)),
                d_tot=float(rng.uniform(0.2, 0.8)),
                rho_U=float(rng.uniform(0.2, 3.0)),
                rho_A=float(rng.uniform(0.1, 2.0)),
                alpha_n=float(rng.uniform(0.05, 2.0)),
                alpha_g=float(rng.uniform(0.1, 4.0)),
                alpha_d=float(rng.uniform(0.05, 1.5)),
                m=float(rng.uniform(1.2, 3.0)),
                kappa_ref=float(rng.uniform(0.05, 2.0)),
                kappa_u=float(rng.uniform(0.05, 2.0)),
                kappa_a=float(rng.uniform(0.05, 2.0)),
                kappa_dis=float(rng.uniform(0.05, 2.0)),
                kappa_cu=float(rng.uniform(0.05, 2.0)),
                kappa_ca=float(rng.uniform(0.05, 2.0)),
                kappa_du=float(rng.uniform(0.05, 2.0)),
                kappa_da=float(rng.uniform(0.05, 2.0)),
            ).validate()
            for u, a in ((0.05, 0.01), (0.4, 0.26), (3.0, 1.0), (10.0, 8.0)):
                du_m, da_m = manuscriptRhs(u, a, p)
                du_c, da_c = M.rhs(u, a, p)
                scale = max(abs(du_c), abs(da_c), 1e-12)
                worst = max(worst, abs(du_m - du_c) / scale,
                            abs(da_m - da_c) / scale)
        self.assertLess(worst, 1e-9, f"disagreement {worst:.3e} off baseline")

    def testTheTotalConcentrationFormIsDetectablyDifferent(self):
        """the guard has teeth: the v5 printed form does NOT reproduce rhs.

        Without this, a test asserting agreement could pass for the wrong
        reason -- e.g. if free and total happened to coincide numerically.
        """
        p = M.Params().validate()
        u, a = 0.4, 0.26
        _, _, C_f, D_f = freePools(u, a, p)
        v_ref_free = C_f * (u / (1.0 + C_f / p.kappa_cu + D_f / p.kappa_du)) / \
            (p.kappa_ref + u / (1.0 + C_f / p.kappa_cu + D_f / p.kappa_du))
        v_ref_total = C_f * u / (p.kappa_ref + u)      # the v5 printed form
        self.assertGreater(abs(v_ref_total - v_ref_free) / v_ref_free, 0.2)


class TestSectionTwoOnePrintsWhatTheCodeCarries(unittest.TestCase):

    def testEverySymbolMapsToAParamsField(self):
        for field in SYMBOL_MAP.values():
            self.assertIn(field, M.Params.__dataclass_fields__,
                          f"§2.1 symbol maps to missing field {field}")

    def testTheStateVariablesAreDeclaredTotal(self):
        sec = _section21()
        self.assertIn("total", sec.lower())
        self.assertRegex(sec, r"\*\*total\*\* pools")

    def testTheClosureIsPrintedInFreeConcentrations(self):
        """all four conservation laws, and no total substrate in a rate law."""
        sec = _section21()
        for lhs in ("u_f = u /(1", "a_f = a /(1",
                    "C_f = C_tot /(1", "D_f = D_tot /(1"):
            self.assertIn(lhs, sec, f"closure line missing: {lhs}")
        # every Michaelis factor must be in a free concentration
        for law in ("K_ref + u_f", "K_U + u_f", "K_dis + a_f", "K_A + a_f"):
            self.assertIn(law, sec, f"rate law not in free concentrations: {law}")
        self.assertNotIn("K_ref + u)", sec, "total-concentration rate law survives")
        self.assertNotIn("K_A + a)", sec, "total-concentration rate law survives")

    def testNucleationAndElongationUseFreeMonomer(self):
        sec = _section21()
        self.assertIn("k_n u_f^m", sec)
        self.assertIn("k_g u_f a_f", sec)

    def testTheStipulationIsStatedWhereTheLawsAre(self):
        """the closure is phenomenological and §2.1 must say so, with numbers."""
        sec = _section21()
        self.assertIn("stipulation", sec)
        self.assertIn("334", sec)
        self.assertIn("2767", sec)


if __name__ == "__main__":
    unittest.main()
