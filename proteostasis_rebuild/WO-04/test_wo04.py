import sys
import unittest
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import bifurcation as bf  # noqa: E402
import run_wo04 as r  # noqa: E402


def scan_example(variant, idx):
    """regenerate the idx-th parameter set of the scan (seeds 100/101/102 for V0/V1/V2)."""
    rng = np.random.default_rng(100 + ("V0", "V1", "V2").index(variant))
    for _ in range(idx + 1):
        p = bf.sample(rng, variant)
    return p


class TestScalar(unittest.TestCase):
    def test_g41_fold_and_cubic(self):  # G4.1
        s = r.g41_scalar()
        self.assertTrue(np.allclose(s["cubic_coeffs"], [-0.3, 0.4, 1.7, 5.0]))
        self.assertEqual(len(s["x_fold"]), 1)
        self.assertAlmostEqual(s["x_fold"][0], 3.890758, places=6)
        self.assertAlmostEqual(s["lambda_fold"], 4.802189, places=6)
        self.assertTrue(s["g2_identity_negative_definite"])
        self.assertLess(max(s["residuals"]), 1e-12)


class TestLegacyMechanism(unittest.TestCase):
    def test_g42_collapse_is_the_gate(self):  # G4.2
        rows = {x["anchoring"]: x for x in r.g42_legacy()}
        a = rows["as_published"]
        self.assertEqual(a["mechanism"], "aggregation_death")
        self.assertAlmostEqual(a["A_at_operational"], 0.25, places=3)
        # the mathematical fold lies beyond the gate
        self.assertGreater(a["P_math_fold"], 2 * a["P_operational"])
        # without the donorless Phi the J-curve has no interior maximum at all
        self.assertTrue(all(x["J_curve_monotone_without_Phi"] for x in rows.values()))


class TestConservativeCount(unittest.TestCase):
    def test_variants_conserve(self):
        okP, okC = bf.conservation_ok()
        self.assertTrue(okP and okC)
        self.assertTrue(bf.reduces_to_wo02())

    def test_g43_analytic(self):  # G4.3 analytical half
        a = bf.analytic_v0()
        for k in ("reduction_unique", "curve_satisfies_other_balances",
                  "dG_dU_decomposition_holds", "G0_formula_holds"):
            self.assertTrue(a[k], k)

    def test_g43_numeric_agrees(self):  # G4.3 numerical half (small deterministic scan)
        s = r.scan("V0", 60, seed=100)
        self.assertEqual(s["count_hist"], {1: 60})
        self.assertEqual(s["n_bad_residual"], 0)

    def test_eps_enters_linearly(self):
        p = scan_example("V1", 16)
        U = np.geomspace(1e-3, 10, 50)
        for e in (0.0, 0.2, 0.9):
            lhs = bf.G(U, {**p, "eps": e})
            rhs = e * bf.eps_coeff(p) + bf.G(U, {**p, "eps": 0.0})
            self.assertLess(np.max(np.abs(lhs - rhs)) / p["s_P"], 1e-12)


class TestStabilityAndFold(unittest.TestCase):
    def test_g44_pattern(self):  # G4.4
        p = scan_example("V1", 16)
        rs, _ = bf.roots(p)
        self.assertEqual(len(rs), 3)
        pat = [bf.classify(bf.full_state(u, p), p)["stable"] for u in rs]
        self.assertEqual(pat, [True, False, True])

    def test_g45_continuation_folds(self):  # G4.5
        c = bf.continue_branch(scan_example("V1", 16), n=2500)
        phys = sorted(f["eps_fold"] for f in c["folds"] if f["physical"])
        self.assertEqual(len(phys), 2)
        self.assertAlmostEqual(phys[0], 0.00391425, places=7)
        self.assertAlmostEqual(phys[1], 0.3349029, places=6)
        for f in c["folds"]:
            self.assertLess(f["min_abs_eig_at_fold"], 1e-4 * f["ref_min_abs_eig_nearby"])
            self.assertLess(f["det_sign_left"] * f["det_sign_right"], 0)

    def test_no_spurious_folds_from_grid_cut(self):
        # regression: cutting the grid to 0<eps<1 before differencing made
        # fake sign changes; each det change must match exactly one fold
        c = bf.continue_branch(scan_example("V1", 38), n=2500)
        self.assertEqual(len(c["det_sign_changes"]), len(c["folds"]))

    def test_v0_has_no_fold(self):
        c = bf.continue_branch(scan_example("V0", 3), n=1200)
        self.assertEqual(c["folds"], [])
        # H(U) strictly decreasing => eps(U) = -H/c strictly increasing
        self.assertTrue(np.all(np.diff(c["eps"]) > 0))


if __name__ == "__main__":
    unittest.main()
