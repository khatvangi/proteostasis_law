"""tests for audit_a6_a7.py. hypothetical dimensionless coefficients only.
run from this directory:  python -m unittest test_audit_a6_a7 -v
"""

import unittest

import numpy as np

import audit_a6_a7 as au


class TestA6Domain(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s = au.deriveA6Symbolic()

    def test_saddle_node_matches_equations_md(self):
        # reproduces A6.2 of EQUATIONS.md before any domain restriction
        self.assertEqual(str(self.s["A_c"]), "(M_tot*k_c + k_a + mu)/(2*beta*k_c)")

    def test_fold_physical_iff_kcM_exceeds_ka_plus_mu(self):
        # A_max - A_c = (k_c M_tot - k_a - mu) / (2 beta k_c)
        self.assertTrue(self.s["A_gap_ok"])

    def test_pole_lies_beyond_titration(self):
        # the convex blow-up of A6.1 at A = a/b sits a0/(k_c beta) past the point where M_free hits 0
        self.assertTrue(self.s["pole_ok"])

    def test_L_sat_never_exceeds_L_crit(self):
        # L_crit - L_sat = (mu+delta)(k_c M_tot - k_a - mu)^2 / (4 beta k_a k_c) >= 0
        self.assertTrue(self.s["L_gap_ok"])

    def test_written_upper_root_unphysical_in_checks_instance(self):
        # checks.py's a=3, b=1, k_a=1, mu+delta=1 at 0.9 L_crit: upper root 1.974 > A_max = 1.75
        q = au.unpackA6(au.P_FOLD)
        L = 0.9 * au.criticalLoads(au.P_FOLD)[0]
        r = np.roots([q["m"] * q["b"], -q["m"] * q["a"], au.P_FOLD["k_a"] * L])
        self.assertGreater(max(r.real), q["A_max"])


class TestA6Saturated(unittest.TestCase):
    def test_root_counts_fold_set(self):
        p = au.P_FOLD
        Lc = au.criticalLoads(p)[0]
        self.assertEqual(len(au.steadyStatesA6(0.9 * Lc, p)), 1)    # checks.py claims 2
        self.assertEqual(len(au.steadyStatesA6(0.99 * Lc, p)), 3)   # bistable window
        self.assertEqual(len(au.steadyStatesA6(1.1 * Lc, p)), 1)    # checks.py claims 0

    def test_no_fold_set_same_checks_numbers(self):
        # same a, b, k_a, mu+delta as checks.py but k_c M_tot < k_a + mu: never bistable
        p = au.P_NOFOLD
        Lc = au.criticalLoads(p)[0]
        for f in np.linspace(0.05, 3.0, 300):
            self.assertEqual(len(au.steadyStatesA6(f * Lc, p)), 1)

    def test_closed_form_matches_grid(self):
        for p in (au.P_FOLD, au.P_NOFOLD):
            Lc = au.criticalLoads(p)[0]
            for f in (0.5, 0.9, 0.99, 1.1, 3.0):
                self.assertEqual(len(au.steadyStatesA6(f * Lc, p)), au.countSignChanges(f * Lc, p))

    def test_above_Lcrit_converges_to_finite_titrated_state(self):
        p = au.P_FOLD
        L = 1.1 * au.criticalLoads(p)[0]
        r = au.integrateA6(L, p)
        A_hi = p["k_a"] * L / ((p["k_a"] + p["mu"]) * (p["mu"] + p["delta"]))
        self.assertAlmostEqual(r["A_final"], A_hi, places=6)
        self.assertGreater(A_hi, au.unpackA6(p)["A_max"])

    def test_written_model_goes_unphysical(self):
        p = au.P_FOLD
        r = au.integrateA6(1.1 * au.criticalLoads(p)[0], p, saturated=False)
        self.assertTrue(r["hit_pole"])
        self.assertLess(r["M_free_final_written"], 0.0)

    def test_hysteresis_jumps(self):
        h = au.hysteresisSweep(au.P_FOLD)
        self.assertLessEqual(abs(h["L_up_jump"] - h["L_crit"]), 2 * h["dL"])
        self.assertLessEqual(abs(h["L_down_jump"] - h["L_sat"]), 2 * h["dL"])

    def test_random_sweep_topology(self):
        r = au.randomSweep(n=5000)
        self.assertEqual(sum(r["violations"].values()), 0, r["violations"])
        self.assertGreater(r["n_bistable"], 0)

    def test_zero_loss_is_the_only_runaway(self):
        p = dict(au.P_FOLD, mu=0.0, delta=0.0)
        r = au.integrateA6(2.0 * au.criticalLoads(au.P_FOLD)[0], p, t_end=200.0)
        self.assertGreater(r["A_final"], 100 * au.unpackA6(p)["A_max"])

    def test_disaggregation_needing_free_chaperone_still_finite(self):
        p = au.P_FOLD
        L = 1.1 * au.criticalLoads(p)[0]
        for d_deg in (0.0, 0.05):
            dfun = lambda Mf, d=d_deg: d + 0.7 * Mf / p["M_tot"]
            r = au.integrateA6(L, p, delta_of_free=dfun)
            self.assertAlmostEqual(r["A_final"], p["k_a"] * L / ((p["k_a"] + p["mu"]) * (p["mu"] + d_deg)), places=5)


class TestA7(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.s = au.deriveA7Symbolic()

    def test_face_derivative_has_moving_boundary_term(self):
        txt = str(self.s["phidot_face"])
        self.assertIn("Derivative(N_req(mu(t)), mu(t))*Derivative(mu(t), t)", txt)
        self.assertIn("k_dN", txt)

    def test_A71_holds_but_state_leaves_V(self):
        r = au.movingQuotaCounterexample()
        self.assertEqual(r["A71_slack_at_t0"], 0.0)
        self.assertLess(r["phidot_at_t0"], 0.0)
        self.assertLess(r["gap_t_0p1"], 0.0)

    def test_quota_at_mu0_includes_native_degradation(self):
        self.assertEqual(str(self.s["quota_mu0_with_kdN"]), "n_0*(k_dN + k_u)")

    def test_growth_interval_closes_at_positive_mu(self):
        self.assertEqual(str(self.s["mu_star"]), "sqrt(k_u)*sqrt(n_0)/sqrt(n_1)")
        s0, n0, n1, ku = 10.0, 1.0, 2.0, 0.1
        yc = (np.sqrt(n0) + np.sqrt(n1 * ku)) ** 2 / s0
        lo, hi = au.feasibleGrowthInterval(1.5 * yc, s0, n0, n1, ku)
        self.assertGreater(lo, 0.0)                     # no viability near mu = 0 with s = s0 mu
        self.assertIsNone(au.feasibleGrowthInterval(0.99 * yc, s0, n0, n1, ku))
        lo2, hi2 = au.feasibleGrowthInterval(1.01 * yc, s0, n0, n1, ku)
        self.assertTrue(lo2 < np.sqrt(ku * n0 / n1) < hi2)


if __name__ == "__main__":
    unittest.main()
