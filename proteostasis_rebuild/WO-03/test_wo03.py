import sys
import unittest
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import chaperone as ch  # noqa: E402
import run_wo03 as r  # noqa: E402


class TestEquilibrium(unittest.TestCase):
    def test_audit_values(self):  # G3.1
        a = r.g31_equilibrium()
        self.assertAlmostEqual(a["at_M_T_50"]["C_f_exact"], 6.58872, places=5)
        self.assertAlmostEqual(a["at_M_T_50"]["fold_ratio_exact_over_legacy"], 1.75382, places=5)
        self.assertLess(a["max_residual_client"], 1e-10)

    def test_legacy_closure_violates_client_balance(self):
        # implied bound chaperone C_T - C_free_legacy exceeds total client
        C_T, M_T, K = 50.0, 10.0, 1.0
        implied_bound = C_T - ch.free_legacy(C_T, M_T, K)
        self.assertGreater(implied_bound, M_T)

    def test_limits(self):
        # K -> 0: B = min(C_T, M_T); M_T -> 0: B -> 0
        self.assertAlmostEqual(ch.bound_exact(50.0, 20.0, 1e-12), 20.0, places=6)
        self.assertAlmostEqual(ch.bound_exact(50.0, 1e-12, 1.0), 0.0, places=9)


class TestCycle(unittest.TestCase):
    def test_K_M_governs_driven_steady_state(self):  # G3.2
        rows = r.g32_cycle()["simple_cycle"]
        for x in rows:
            self.assertLess(x["rel_err_K_M"], 1e-6)
        self.assertGreater(rows[-1]["rel_err_K_d"], 0.3)

    def test_ultra_affinity_needs_drive(self):
        f = r.g32_cycle()["four_state"]
        self.assertTrue(f["driven_tighter_than_both_states"])
        self.assertTrue(f["undriven_within_state_range"])

    def test_four_state_distribution_normalised(self):
        occ, pi = ch.four_state_occupancy(1.0, 10, 100, 0.01, 0.01, 10, 0.01, 0.1)
        self.assertAlmostEqual(pi.sum(), 1.0, places=12)
        self.assertTrue(np.all(pi > -1e-12))


class TestCompetition(unittest.TestCase):
    def test_conservation_and_reduction(self):  # G3.3
        okP, okC = ch.ext_conservation()
        self.assertTrue(okP and okC)
        same, stays = ch.ext_reduces_symbolically()
        self.assertTrue(same and stays)

    def test_theta_increases_with_nascent_load(self):
        rows = r.g33_competition()["theta_sweep"]
        th = [x["theta_BX_over_C_T"] for x in rows]
        self.assertTrue(all(b > a for a, b in zip(th, th[1:])))


class TestQSS(unittest.TestCase):
    def test_tqssa_accuracy(self):  # G3.5
        q = r.g35_qss()
        self.assertLess(q["max_rel_err_U"], 1e-2)
        self.assertLess(q["max_rel_err_A"], 1e-2)


if __name__ == "__main__":
    unittest.main()
