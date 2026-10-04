"""deterministic tests for the network-capacity model.
run: PYTHONDONTWRITEBYTECODE=1 python -m unittest test_network -v
"""
import itertools
import sys
import unittest

import numpy as np

sys.dont_write_bytecode = True

import feedback as fbm
import model as M
from pools import WO07_ACC, findPools, readSchmidtGlucose, readWO07Pools

P = findPools()


def makeBipartite(seed, K=5, Mm=6):
    """single-machine routes, unit usage: the setting of Hall's theorem."""
    rng = np.random.default_rng(seed)
    routes, cls = [], []
    for k in range(K):
        ms = rng.choice(Mm, size=rng.integers(1, 4), replace=False)
        for m in ms:
            routes.append(m)
            cls.append(k)
    U = np.zeros((len(routes), Mm))
    U[np.arange(len(routes)), routes] = 1.0
    cls = np.array(cls)
    return {"U": U, "cls": cls, "d": rng.dirichlet(np.ones(K)), "c": 10 ** rng.uniform(-1, 1, Mm),
            "primary": np.array([np.flatnonzero(cls == k)[0] for k in range(K)])}


class TestPools(unittest.TestCase):
    def test_schmidt_path_reproduces_wo07(self):
        a, b = readSchmidtGlucose(WO07_ACC), readWO07Pools()
        for m in b:
            self.assertAlmostEqual(a[m] / b[m], 1.0, places=5, msg=m)

    def test_all_pools_positive(self):
        self.assertEqual(set(P), set(M.MACHINES))
        self.assertTrue(all(v > 0 for v in P.values()))


class TestFormulation(unittest.TestCase):
    def test_ordering_A_ge_C_ge_B(self):
        for s in range(40):
            for st in ("bio", "unstructured"):
                i = M.makeInstance(s, structure=st, pools=P)
                a, b, c = M.solveLamA(i), M.solveLamB(i), M.solveLamC(i)
                self.assertGreaterEqual(a, c * (1 - 1e-9))
                self.assertGreaterEqual(c, b * (1 - 1e-9))

    def test_all_to_all_network_is_a_pool(self):
        # every class can use every machine at unit cost -> C collapses to A
        rng = np.random.default_rng(0)
        K, Mm = 4, 5
        U = np.tile(np.eye(Mm), (K, 1))
        cls = np.repeat(np.arange(K), Mm)
        i = {"U": U, "cls": cls, "d": rng.dirichlet(np.ones(K)), "c": rng.random(Mm) + 0.1}
        self.assertAlmostEqual(M.solveLamC(i), M.solveLamA(i), places=9)

    def test_single_route_network_is_B(self):
        for s in range(20):
            i = M.makeInstance(s, pools=P)
            keep = i["primary"]
            j = dict(i, U=i["U"][keep], cls=i["cls"][keep], primary=np.arange(len(keep)))
            self.assertAlmostEqual(M.solveLamC(j) / M.solveLamB(j), 1.0, places=7)

    def test_hall_min_cut_equals_lp(self):
        for s in range(30):
            i = makeBipartite(s)
            K = len(i["d"])
            best = np.inf
            for n in range(1, K + 1):
                for S in itertools.combinations(range(K), n):
                    nb = set(np.flatnonzero(i["U"][np.isin(i["cls"], S)].sum(axis=0) > 0))
                    best = min(best, i["c"][list(nb)].sum() / i["d"][list(S)].sum())
            self.assertAlmostEqual(M.solveLamC(i) / best, 1.0, places=7)

    def test_feasible_set_is_an_interval(self):
        for s in range(15):
            i = M.makeInstance(s, pools=P)
            lam = M.solveLamC(i)
            for f in (0.01, 0.3, 0.7, 0.999):
                self.assertIsNotNone(M.solveMinCost(i, f * lam), f)
            self.assertIsNone(M.solveMinCost(i, 1.001 * lam))
            self.assertIsNone(M.solveMinCost(i, 5 * lam))

    def test_strong_duality_certificate(self):
        for s in range(15):
            i = M.makeInstance(s, pools=P)
            lam, x, y = M.solveLamC(i, duals=True)
            self.assertTrue(np.all(y >= -1e-9))
            self.assertAlmostEqual(y @ i["c"] / lam, 1.0, places=6)

    def test_reserve_margin_is_lambda_star(self):
        # max reserve at load lam is 1 - lam/lambda*: not an independent objective
        i = M.makeInstance(3, pools=P)
        lam = M.solveLamC(i)
        self.assertIsNotNone(M.solveMinCost(i, 0.4 * lam, cap_scale=0.4 + 1e-7))
        self.assertIsNone(M.solveMinCost(i, 0.4 * lam, cap_scale=0.4 - 1e-4))


class TestKnockdowns(unittest.TestCase):
    def test_B_pairs_are_weakest_link(self):
        for s in range(10):
            i = M.makeInstance(s, pools=P)
            for kappa in (0.0, 0.5):
                _, W1, W2 = M.findKnockdowns(i, "B", kappa)
                for (a, b), w in W2.items():
                    self.assertAlmostEqual(w, min(W1[a], W1[b]), places=12)
                st = M.findPairStats(W1, W2, 0.5)
                self.assertEqual(st["synthetic_lethal_pairs"], 0)
                self.assertEqual(st["pairs_beyond_weakest_link"], 0)

    def test_A_pairs_closed_form(self):
        i = M.makeInstance(4, pools=P)
        C = i["c"].sum()
        _, W1, W2 = M.findKnockdowns(i, "A", 0.5)
        for (a, b), w in W2.items():
            self.assertAlmostEqual(w, 1 - 0.5 * (i["c"][a] + i["c"][b]) / C, places=12)

    def test_structural_cuts_are_minimal_and_lethal(self):
        for s in range(10):
            i = M.makeInstance(s, pools=P)
            for cut in M.findStructuralCuts(i):
                c = i["c"].copy()
                c[list(cut)] = 0
                self.assertEqual(M.solveLamC(i, c), 0.0)
                for n in range(1, len(cut)):
                    for sub in itertools.combinations(cut, n):
                        c = i["c"].copy()
                        c[list(sub)] = 0
                        self.assertGreater(M.solveLamC(i, c), 0.0)

    def test_C_knockdown_never_beats_intact(self):
        i = M.makeInstance(7, pools=P)
        _, W1, W2 = M.findKnockdowns(i, "C", 0.5)
        self.assertTrue(np.all(W1 <= 1 + 1e-9))
        for (a, b), w in W2.items():
            self.assertLessEqual(w, min(W1[a], W1[b]) + 1e-9)


class TestSaturationAndPareto(unittest.TestCase):
    def test_events_strictly_before_failure(self):
        i = M.makeInstance(2, pools=P)
        ev = M.findSaturationEvents(i, n=50)
        self.assertTrue(all(0 < f < 1 for _, f in ev["order"]))

    def test_pareto_front_monotone(self):
        i = M.makeInstance(5, pools=P)
        pts = M.findParetoFront(i, 0.5, n=10)
        self.assertAlmostEqual(pts[0]["reserve_margin"], 0.5, places=6)
        atp = [p["atp"] for p in pts]
        self.assertTrue(all(b <= a + 1e-7 * max(1, a) for a, b in zip(atp, atp[1:])))


class TestFeedback(unittest.TestCase):
    def test_no_feedback_no_jump_onset_at_lambda_star(self):
        i = M.makeInstance(1, pools=P)
        r = fbm.findCascade(i, {"kind": "none", "gamma": 0, "beta": 0}, n=41)
        step = r["grid"][1] - r["grid"][0]
        self.assertEqual(r["n_jumps"], 0)
        self.assertLess(r["hysteresis"], 1e-6)
        self.assertTrue(1.0 <= r["onset_frac"] < 1.0 + step + 1e-12)

    def test_overload_feedback_cannot_move_viability_threshold(self):
        i = M.makeInstance(1, pools=P)
        fb = {"kind": "titration", "gamma": 5.0, "beta": 5.0}
        r = fbm.findCascade(i, fb, n=41)
        self.assertGreaterEqual(r["onset_frac"], 1.0)

    def test_jump_detector_positive_control(self):
        # hand-built loop: one folding class served only by DnaK, one aggregate
        # class served only by DnaK. unserved folding -> aggregates (gamma),
        # unserved aggregates titrate DnaK (beta). gain gamma*beta > 1 gives a
        # runaway at lambda* (analytically: once z>0, DnaK capacity collapses)
        Mm = len(M.MACHINES)
        U = np.zeros((6, Mm))
        for r in range(6):
            U[r, M.MIDX["DnaK"] if r in (1, 5) else M.MIDX["TF"]] = 1.0
        i = {"U": U, "cls": np.arange(6), "d": np.array([0.01, 0.6, 0.01, 0.01, 0.01, 0.36]),
             "c": np.ones(Mm) * 100.0, "cost": np.ones(6), "primary": np.arange(6)}
        i["c"][M.MIDX["DnaK"]] = 1.0
        r = fbm.findCascade(i, {"kind": "titration", "gamma": 3.0, "beta": 3.0}, n=61)
        self.assertGreaterEqual(r["n_jumps"], 1)


class TestDeterminism(unittest.TestCase):
    def test_same_seed_same_instance(self):
        a, b = M.makeInstance(11, pools=P), M.makeInstance(11, pools=P)
        for k in ("U", "cls", "d", "c", "cost"):
            np.testing.assert_array_equal(a[k], b[k])
        self.assertEqual(M.solveLamC(a), M.solveLamC(b))


if __name__ == "__main__":
    unittest.main()
