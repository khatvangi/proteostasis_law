import sys
import unittest
from pathlib import Path

import numpy as np
import sympy as sp

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import model as m  # noqa: E402
import proofs  # noqa: E402


class TestConservation(unittest.TestCase):
    def test_symbolic_totals(self):  # G2.1
        self.assertEqual(proofs.g21_conservation(),
                         {"P_T_identity": True, "C_T_identity": True})

    def test_every_internal_flux_balanced(self):
        # each non-source, non-sink flux must net to zero protein and zero chaperone
        sources_sinks = {"synth_native", "synth_nonnative", "synth_chaperone",
                         "degrade_U", "degrade_A"}
        for f, st in m.STOICH.items():
            if f in sources_sinks:
                continue
            self.assertEqual(sum(m.PROTEIN_W[s] * c for s, c in st.items()), 0, f)
            self.assertEqual(sum(m.CHAP_W[s] * c for s, c in st.items()), 0, f)

    def test_negative_control_phi_like_source_breaks_conservation(self):
        # add a legacy-style state-dependent inflow to U with no donor:
        # the P_T identity must then fail
        bad = list(m.RHS_SYM)
        bad[1] = bad[1] + m.k_a * m.U**2 * m.eps
        dPT = sum(m.PROTEIN_W[s] * r for s, r in zip(m.STATE_NAMES, bad))
        want = m.s_P - m.k_d * m.U - m.k_dA * m.A - m.mu * m.P_T
        self.assertNotEqual(sp.simplify(dPT - want), 0)


class TestDimensionsAndInvariance(unittest.TestCase):
    def test_dimensions(self):  # G2.2
        self.assertEqual(proofs.g22_dimensions()["n_bad_terms"], 0)

    def test_orthant_invariance(self):  # G2.3
        self.assertTrue(all(proofs.g23_invariance().values()))

    def test_invariance_negative_control(self):
        # a sink on U that does not vanish at U=0 must be caught
        a, b = sp.symbols("a b", nonnegative=True)
        face = (m.RHS_SYM[1] - m.k_d * m.N).subs(m.U, 0)
        face = face.subs({m.eps: a / (a + b), m.phi: sp.Rational(1, 2)})
        num = sp.expand(sp.cancel(face * (a + b)))
        coeffs = sp.Poly(num, *num.free_symbols).coeffs()
        self.assertFalse(all(c >= 0 for c in coeffs))


class TestNumeric(unittest.TestCase):
    def test_small_integration_sample(self):  # G2.4 (reduced n for speed)
        r = proofs.g24_numeric(n=8, seed=1)
        self.assertEqual(r["n_integration_failures"], 0)
        self.assertLess(r["max_rel_err_P_T"], 1e-6)
        self.assertLess(r["max_rel_err_C_T"], 1e-6)
        self.assertGreater(r["most_negative_state"], -1e-9)

    def test_scenario_params_balanced(self):
        p = m.scenario_params()
        self.assertAlmostEqual(p["s_P"], p["mu"] * p["P_T"])
        self.assertTrue(np.all(np.isfinite(m.pvec(p))))


class TestLegacyPhi(unittest.TestCase):
    def test_phi_has_no_donor(self):  # G2.5
        r = proofs.g25_phi_phantom()
        self.assertFalse(r["donor_pool_found"])
        self.assertTrue(r["unaccounted_equals_J_bare_times_Phi_minus_1"])
        self.assertGreater(r["phantom_fraction_of_inflow_at_P_dagger"],
                           r["phantom_fraction_of_inflow_at_P_star"])


if __name__ == "__main__":
    unittest.main()
