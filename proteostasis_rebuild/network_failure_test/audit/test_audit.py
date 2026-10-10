"""AUDIT 2026-10-09 tests (new work). bounded: runs in about a minute.
run from audit/:  PYTHONDONTWRITEBYTECODE=1 python -m unittest test_audit -v
"""
import itertools
import sys
import unittest
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import cf_model as C  # noqa: E402
import model as M  # noqa: E402
from pools import findPools  # noqa: E402

P = findPools()


class TestCounterfactualHarness(unittest.TestCase):
    def test_orig_mode_reproduces_original_instance(self):
        # W=0 removes the only draw that moved to a separate stream
        for s in range(30):
            a, b = M.makeInstance(s, W=0, pools=P), C.makeInstanceCF(s, W=0, pools=P)
            self.assertEqual(a["routes"], b["routes"])
            for k in ("U", "d", "c"):
                np.testing.assert_array_equal(a[k], b[k])

    def test_nocut_has_no_single_cut_and_no_lethal_single_knockout(self):
        for s in range(15):
            i = C.makeInstanceCF(s, routes="nocut", pools=P)
            self.assertEqual(M.findStructuralCuts(i, 1), [])
            _, W1, _ = M.findKnockdowns(i, "C", 0.0)
            self.assertTrue(np.all(W1 > 0))

    def test_extra_routes_keep_base_instance_paired(self):
        for s in range(15):
            a = C.makeInstanceCF(s, pools=P)
            b = C.makeInstanceCF(s, routes="nocut", agg_free=True, pools=P)
            n = len(a["routes"])
            self.assertEqual(a["routes"], b["routes"][:n])
            np.testing.assert_array_equal(a["U"], b["U"][:n])
            np.testing.assert_array_equal(a["d"], b["d"])
            np.testing.assert_array_equal(a["c"], b["c"])


class TestSingleCutIsStructural(unittest.TestCase):
    def test_cuts_independent_of_capacity_abundance_and_load(self):
        # same seed -> same routes in the original code, whatever W / abundance
        for s in range(20):
            ref = M.findStructuralCuts(M.makeInstance(s, W=0, pools=P), 1)
            for kw in (dict(W=4), dict(W=2, abundance="blind")):
                self.assertEqual(M.findStructuralCuts(M.makeInstance(s, pools=P, **kw), 1), ref)

    def test_single_cut_rate_follows_from_rule_probabilities(self):
        # no cut iff GRO has its KJE alternative (0.2), HTP primary is HtpG
        # alone with a KJE alternative (0.25), DK has a GRO or TF alternative
        # (0.65), DMG has >= 2 proteases (0.784), AGG has SEQ or LON (0.79)
        analytic = 1 - 0.2 * 0.25 * 0.65 * (1 - 0.3 ** 3 - 3 * 0.7 * 0.09) * 0.79
        n = 4000
        hits = sum(bool(C.findCutMachines(M.findRoutesBio(np.random.default_rng(s)))) for s in range(n))
        se = np.sqrt(analytic * (1 - analytic) / n)
        self.assertLess(abs(hits / n - analytic), 4 * se)

    def test_clpb_exclusivity_is_not_the_driver(self):
        n = 3000
        a = b = 0
        for s in range(n):
            R = M.findRoutesBio(np.random.default_rng(s))
            a += bool(C.findCutMachines(R))
            b += bool(C.findCutMachines([r for r in R if r[0] != "AGG"]))
        # dropping the aggregate class altogether changes the rate by < 2 points
        self.assertLess((a - b) / n, 0.02)


class TestTheorems(unittest.TestCase):
    def test_partial_knockdown_cannot_be_lethal_below_kappa(self):
        # lambda*(c') >= lambda*(kappa c) = kappa lambda* whenever c' >= kappa c
        for s in range(6):
            i = M.makeInstance(s, pools=P)
            _, W1, W2 = M.findKnockdowns(i, "C", 0.5)
            self.assertGreaterEqual(min(W2.values()), 0.5 - 1e-9)

    def test_structural_SL_pairs_are_size2_cuts(self):
        for s in range(8):
            i = M.makeInstance(s, pools=P)
            r = C.findPairSplit(i)
            n2 = sum(len(c) == 2 for c in M.findStructuralCuts(i, 2))
            self.assertLessEqual(r["SL_struct_rho0.5"], n2)

    def test_route_preference_moves_saturation_not_viability(self):
        for s in range(6):
            a = C.makeInstanceCF(s, pools=P)
            b = C.makeInstanceCF(s, cost="reverse", pools=P)
            self.assertAlmostEqual(M.solveLamC(a) / M.solveLamC(b), 1.0, places=9)


class TestKnownSolverFailure(unittest.TestCase):
    def test_presolve_false_infeasible_reproduced_and_resolved(self):
        # original failure: alpha=0.3, seed 65, intact lambda*_C. the LP always
        # has the feasible point x=0, lambda=0, so 'infeasible' is a solver error
        i = M.makeInstance(65, W=2, alpha=0.3, pools=P)
        with self.assertRaises(RuntimeError):
            M.solveLamC(i)
        U, cls, d, c = i["U"], i["cls"], i["d"], i["c"]
        R = len(U)
        Aeq = np.zeros((6, R + 1))
        Aeq[cls, np.arange(R)] = 1
        Aeq[:, R] = -d
        obj = np.zeros(R + 1)
        obj[R] = -1
        r = linprog(obj, A_ub=np.hstack([(U / c).T, np.zeros((15, 1))]), b_ub=np.ones(15),
                    A_eq=Aeq, b_eq=np.zeros(6), method="highs", options={"presolve": False})
        self.assertEqual(r.status, 0)
        x = r.x
        self.assertTrue(np.all(U.T @ x[:R] <= c * (1 + 1e-7)))
        self.assertGreater(x[R], 0)


if __name__ == "__main__":
    unittest.main()
