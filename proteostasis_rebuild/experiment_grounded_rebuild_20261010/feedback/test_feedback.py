"""algebra checks for feedback.py. every parameter here is a dimensionless placeholder;
nothing in this file is, or may be reported as, a physiological estimate."""
import unittest
import numpy as np
import feedback as fb

BASE = dict(L=1.0, k=2.0, K=1.0, mu=0.1, delta=0.0, b0=0.05, bm=1.0, Ks=1.0, n=2)


def drawParams(rng):
    # random positive placeholders spanning several decades
    lg = lambda lo, hi: 10 ** rng.uniform(lo, hi)
    return dict(L=lg(-1, 1), k=lg(-1, 1), K=lg(-1, 1), mu=lg(-2, 0), delta=lg(-3, 0),
                b0=lg(-3, -1), bm=lg(-1, 1), Ks=lg(-1, 1), n=rng.choice([1, 2, 4]))


class TestSignsAndUniqueness(unittest.TestCase):
    def test_jacobian_sign_pattern_and_stability_without_delay(self):
        rng = np.random.default_rng(0)
        for _ in range(200):
            p = drawParams(rng)
            eq = fb.findEquilibria(p)
            # fixed production, constant growth: exactly one equilibrium (index argument)
            self.assertEqual(len(eq), 1)
            U, C = eq[0]
            J = fb.findJacobian(U, C, p)
            Jn = fb.findJacobianNumeric(U, C, p)
            np.testing.assert_allclose(J, Jn, rtol=1e-5, atol=1e-9)
            self.assertLess(J[0, 0], 0)   # more client -> faster removal
            self.assertLess(J[0, 1], 0)   # more chaperone -> less client
            self.assertGreater(J[1, 0], 0)  # more client -> less free chaperone -> more synthesis
            self.assertLess(J[1, 1], 0)   # more chaperone -> more free -> less synthesis + dilution
            self.assertLess(np.trace(J), 0)
            self.assertGreater(np.linalg.det(J), 0)

    def test_box_is_positively_invariant(self):
        rng = np.random.default_rng(1)
        for _ in range(100):
            p = drawParams(rng)
            Umax = p["L"] / (p["mu"] + p["delta"])
            Cmin, Cmax = fb.findCapMin(p), fb.findCapMax(p)
            for s in np.linspace(0, 1, 21):
                C = Cmin + s * (Cmax - Cmin)
                U = s * Umax
                self.assertGreater(fb.findRates(0.0, C, p)[0], 0)
                self.assertLess(fb.findRates(Umax, C, p)[0], 0)
                self.assertGreater(fb.findRates(U, Cmin, p)[1], 0)
                self.assertLess(fb.findRates(U, Cmax, p)[1], 0)


    def test_divergence_negative_everywhere_in_box(self):
        # bendixson: no closed orbit without delay
        rng = np.random.default_rng(6)
        for _ in range(100):
            p = drawParams(rng)
            U = rng.uniform(1e-3, 1) * p["L"] / (p["mu"] + p["delta"])
            C = fb.findCapMin(p) + rng.uniform(0, 1) * (fb.findCapMax(p) - fb.findCapMin(p))
            self.assertLess(np.trace(fb.findJacobianNumeric(U, C, p)), 0)


class TestConservationAndBounds(unittest.TestCase):
    def test_influx_equals_served_plus_degraded_plus_diluted(self):
        rng = np.random.default_rng(2)
        for _ in range(100):
            p = drawParams(rng)
            U, C = fb.findEquilibria(p)[0]
            self.assertAlmostEqual(sum(fb.findFluxPartition(U, C, p)), p["L"], places=9)

    def test_load_above_max_service_is_never_fully_served(self):
        # whatever the gains, served flux <= k*Cmax, so the excess stays unserved
        rng = np.random.default_rng(3)
        for _ in range(100):
            p = drawParams(rng)
            p["L"] = fb.findCapMax(p) * p["k"] * rng.uniform(1.01, 10)
            U, C = fb.findEquilibria(p)[0]
            served, _, _ = fb.findFluxPartition(U, C, p)
            self.assertLessEqual(C, fb.findCapMax(p) * (1 + 1e-12))
            self.assertLessEqual(served, p["k"] * fb.findCapMax(p))
            self.assertGreaterEqual(p["L"] - served, p["L"] - p["k"] * fb.findCapMax(p) - 1e-12)

    def test_degradation_balances_the_books_by_destroying_clients(self):
        # raising the protease route lowers U* and the chaperone-served flux, while
        # the influx balance still closes: balance is not productive folding
        p = dict(BASE, L=15.0)
        served, deg = [], []
        for d in [0.0, 0.1, 1.0, 10.0]:
            q = dict(p, delta=d)
            U, C = fb.findEquilibria(q)[0]
            s, g, _ = fb.findFluxPartition(U, C, q)
            served.append(s)
            deg.append(g)
        self.assertTrue(all(np.diff(served) < 0))
        self.assertTrue(all(np.diff(deg) > 0))

    def test_growth_coupled_ceiling_is_an_allocation_not_a_rate(self):
        # with synthesis proportional to growth, Cmax = P*(phi0+phim) for any mu
        for mu in [0.01, 0.1, 1.0]:
            p = dict(mode="growth", P=10.0, phi0=0.01, phim=0.05, k=1.0, K=1.0, mu=mu,
                     delta=0.0, Ks=1.0, n=2, L=100.0)
            U, C = fb.findEquilibria(p)[0]
            self.assertLessEqual(C, 0.6 + 1e-12)


class TestSensitivityAndDelay(unittest.TestCase):
    def test_feedback_divides_steady_state_error_by_one_plus_loop_gain(self):
        rng = np.random.default_rng(4)
        for _ in range(50):
            p = drawParams(rng)
            U, C = fb.findEquilibria(p)[0]
            a, b, c, e, mu = fb.findLoopTerms(U, C, p)
            G = fb.findLoopGain(a, b, c, e, mu)
            det = np.linalg.det(fb.findJacobian(U, C, p))
            self.assertAlmostEqual(det / (a * mu * (1 + G)), 1.0, places=9)
            # dF*/dL by finite difference against -(C K/(K+U)^2) / (a (1+G))
            h = 1e-6 * p["L"]
            U2, C2 = fb.findEquilibria(dict(p, L=p["L"] + h))[0]
            dF = (fb.findFree(C2, U2, p["K"]) - fb.findFree(C, U, p["K"])) / h
            pred = -C * p["K"] / (p["K"] + U) ** 2 / (a * (1 + G))
            self.assertLess(dF, 0)
            self.assertAlmostEqual(dF / pred, 1.0, places=3)

    def test_high_gain_loop_destabilises_beyond_critical_delay(self):
        p = dict(BASE, bm=5.0, n=4)
        U, C = fb.findEquilibria(p)[0]
        a, b, c, e, mu = fb.findLoopTerms(U, C, p)
        self.assertGreater(fb.findLoopGain(a, b, c, e, mu), 1)
        w, tau = fb.findCriticalDelay(a, b, c, e, mu)
        s = 1j * w
        resid = (s + a) * (s + mu) + (e * s + e * a + b * c) * np.exp(-s * tau)
        self.assertLess(abs(resid), 1e-9)
        T = 60 * (2 * np.pi / w)
        self.assertLess(fb.findLinearDdeGrowth(a, b, c, e, mu, 0.8 * tau, T), 0)
        self.assertGreater(fb.findLinearDdeGrowth(a, b, c, e, mu, 1.25 * tau, T), 0)

    def test_loop_gain_above_one_always_has_a_finite_critical_delay(self):
        rng = np.random.default_rng(5)
        hits = 0
        for _ in range(300):
            p = drawParams(rng)
            U, C = fb.findEquilibria(p)[0]
            terms = fb.findLoopTerms(U, C, p)
            if fb.findLoopGain(*terms) > 1:
                hits += 1
                self.assertIsNotNone(fb.findCriticalDelay(*terms))
        self.assertGreater(hits, 10)

    def test_weak_loop_without_crossing_is_stable_for_every_delay(self):
        a, b, c, e, mu = 1.0, 0.1, 0.1, 0.05, 0.5   # G < 1 and a^2+mu^2 > e^2
        self.assertLess(fb.findLoopGain(a, b, c, e, mu), 1)
        self.assertIsNone(fb.findCriticalDelay(a, b, c, e, mu))


class TestGrowthToxicity(unittest.TestCase):
    def test_growth_inhibition_can_create_a_collapsed_second_state(self):
        # load-dependent growth arrest removes dilution: a served state, a saddle,
        # and a collapsed state with the chaperone pool already at its ceiling
        p = dict(mode="growth", P=10.0, phi0=0.01, phim=0.05, k=0.1, K=1.0, mu=0.1,
                 Ks=1.0, n=2, Ktox=5.0, ntox=2, L=0.25, delta=0.001, Uscan=1e5)
        eq = fb.findEquilibria(p)
        self.assertEqual(len(eq), 3)
        dets = [np.linalg.det(fb.findJacobianNumeric(U, C, p)) for U, C in eq]
        self.assertGreater(dets[0], 0)
        self.assertLess(dets[1], 0)          # saddle separates the two basins
        self.assertGreater(eq[2][1], 0.99 * fb.findCapMax(p))
        # raising load slightly removes the served state altogether
        self.assertEqual(len(fb.findEquilibria(dict(p, L=0.32, delta=0.0005))), 1)


if __name__ == "__main__":
    unittest.main()
