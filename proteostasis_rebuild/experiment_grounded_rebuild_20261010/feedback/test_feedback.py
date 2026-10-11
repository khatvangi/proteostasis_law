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
        # at equilibria only: served flux <= k*Cmax, so the excess stays unserved.
        # trajectories that start above Cmax can exceed it transiently (TestCeiling)
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
        # the influx balance still closes. "served" is chaperone turnover; the model
        # has no native-product state, so this says nothing about native yield
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
        a, b, c, e, mu = 1.0, 0.1, 0.1, 0.05, 0.5   # G < 1 (which already forces e < mu)
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


def drawLoop(rng):
    # random positive linear loop terms a, b, c, e, mu over two decades
    return tuple(10 ** rng.uniform(-1, 1, size=5))


def findOscillatingHistory(t):
    # nonconstant history on [-tau, 0], so stability is not checked from one special start
    return 1e-3 * np.cos(3.0 * t), 1e-3 * np.sin(1.7 * t)


class TestDelayTheorem(unittest.TestCase):
    """linear, local result: stable at every delay iff G <= 1; G > 1 -> stable iff tau < tau*."""

    def test_crossing_polynomial_signs_follow_the_gain(self):
        rng = np.random.default_rng(10)
        weak = strong = 0
        for _ in range(5000):
            a, b, c, e, mu = drawLoop(rng)
            G = fb.findLoopGain(a, b, c, e, mu)
            A2 = a * a + mu * mu - e * e
            A0 = (a * mu) ** 2 - (e * a + b * c) ** 2
            self.assertAlmostEqual(A0 / ((a * mu) ** 2 * (1 - G * G)), 1.0, places=6)
            if G <= 1:
                weak += 1
                # G <= 1 -> e*a < a*mu (since b*c > 0) -> e < mu -> A2 > 0; and A0 >= 0
                self.assertLess(e, mu)
                self.assertGreater(A2, 0)
                self.assertIsNone(fb.findCriticalDelay(a, b, c, e, mu))
            else:
                strong += 1
                zs = [z for z in np.roots([1.0, A2, A0]) if abs(z.imag) < 1e-12 and z.real > 0]
                self.assertEqual(len(zs), 1)                  # A0 < 0: exactly one positive root
                self.assertGreater(2 * zs[0].real + A2, 0)    # F'(z+) > 0: every crossing destabilises
                w, tau = fb.findCriticalDelay(a, b, c, e, mu)
                s = 1j * w
                f = (s + a) * (s + mu) + (e * s + e * a + b * c) * np.exp(-s * tau)
                self.assertLess(abs(f), 1e-8 * max(1.0, abs((s + a) * (s + mu))))
        self.assertGreater(weak, 500)
        self.assertGreater(strong, 500)

    def test_equality_case_has_no_crossing_and_still_decays(self):
        # G = 1 exactly: A0 = 0, roots z = 0 and z = -A2. z = 0 would be s = 0, but
        # f(0) = a*mu*(1+G) = 2 a mu != 0, so there is no crossing at any finite delay
        a, b, c, e, mu = 1.0, 1.0, 0.5, 0.5, 1.0
        self.assertAlmostEqual(fb.findLoopGain(a, b, c, e, mu), 1.0, places=12)
        self.assertIsNone(fb.findCriticalDelay(a, b, c, e, mu))
        self.assertAlmostEqual(a * mu + e * a + b * c, 2 * a * mu, places=12)
        self.assertLess(fb.findLinearDdeGrowth(a, b, c, e, mu, 1.0, 120, hist=findOscillatingHistory), 0)

    def test_weak_gain_is_stable_at_every_delay(self):
        # G = 0.5: decays from a nonconstant history at short, medium and long delay
        a, b, c, e, mu = 1.0, 1.0, 0.2, 0.3, 1.0
        self.assertAlmostEqual(fb.findLoopGain(a, b, c, e, mu), 0.5, places=12)
        self.assertIsNone(fb.findCriticalDelay(a, b, c, e, mu))
        for tau in [0.5, 5.0, 20.0]:
            g = fb.findLinearDdeGrowth(a, b, c, e, mu, tau, 10 * tau + 60, hist=findOscillatingHistory)
            self.assertLess(g, 0, msg=f"tau={tau}")

    def test_strong_gain_is_unstable_at_every_delay_beyond_tau_star(self):
        # e < mu and e^2 < a^2+mu^2 both hold, yet G = 2.5 > 1, so a finite tau* exists.
        # stable just below tau*, unstable just above and further out (crossings never reverse)
        a, b, c, e, mu = 1.0, 1.0, 2.0, 0.5, 1.0
        self.assertGreater(fb.findLoopGain(a, b, c, e, mu), 1)
        w, tau = fb.findCriticalDelay(a, b, c, e, mu)
        self.assertAlmostEqual(tau, 1.2666, places=3)
        T = 60 * (2 * np.pi / w)
        H = findOscillatingHistory
        self.assertLess(fb.findLinearDdeGrowth(a, b, c, e, mu, 0.9 * tau, T, hist=H), 0)
        for f in [1.1, 2.0, 3.0]:
            self.assertGreater(fb.findLinearDdeGrowth(a, b, c, e, mu, f * tau, T, hist=H), 0,
                               msg=f"{f} tau*")


class TestCeiling(unittest.TestCase):
    """Cmax = (b0+bm)/mu attracts and is forward-invariant; it is not an absolute bound."""
    P = dict(L=5.0, k=2.0, K=1.0, mu=0.5, delta=0.0, b0=0.05, bm=1.0, Ks=1.0, n=2)

    def test_pool_starting_at_or_below_ceiling_stays_below(self):
        p = self.P
        Cmax = fb.findCapMax(p)
        for U0, C0 in [(0.0, Cmax), (10.0, Cmax), (1.0, 0.0), (50.0, 0.5 * Cmax)]:
            _, _, C = fb.findTrajectory(U0, C0, p, 20)
            self.assertLessEqual(C.max(), Cmax * (1 + 1e-9))

    def test_initial_pool_above_ceiling_exceeds_it_and_serves_above_k_cmax(self):
        # counterexample to "the pool can never exceed Cmax": a pre-induced start
        p = self.P
        Cmax = fb.findCapMax(p)
        C0 = 3 * Cmax
        t, U, C = fb.findTrajectory(1.0, C0, p, 20)
        served = p["k"] * C * U / (p["K"] + U)
        self.assertGreater(t[C > Cmax].max(), 3.0)            # above the ceiling for t > 3
        self.assertGreater(served.max(), p["k"] * Cmax)       # served above k*Cmax ...
        self.assertGreater(served.max(), p["L"])              # ... and above the influx
        # the true statement: the excess decays at least as fast as exp(-mu t)
        self.assertTrue(np.all(C - Cmax <= (C0 - Cmax) * np.exp(-p["mu"] * t) + 1e-9))
        self.assertLess(C[-1], Cmax)

    def test_growth_upshift_starts_pool_above_its_new_ceiling(self):
        # nonconstant growth: steady state at mu = 0.25, then mu steps to 0.5
        p = self.P
        U0, C0 = fb.findEquilibria(dict(p, mu=0.25))[0]
        Cmax = fb.findCapMax(p)
        self.assertGreater(C0, Cmax)
        t, _, C = fb.findTrajectory(U0, C0, p, 20)
        self.assertGreater(t[C > Cmax].max(), 1.0)
        self.assertLess(C[-1], Cmax)


class TestModelAssumptions(unittest.TestCase):
    def test_sigma_derivative_at_zero_free_chaperone(self):
        for Ks in [0.5, 2.0]:
            self.assertAlmostEqual(fb.findSigmaActivity(0.0, Ks, 1)[1], -1.0 / Ks, places=12)
            self.assertEqual(fb.findSigmaActivity(0.0, Ks, 2)[1], 0.0)
            self.assertEqual(fb.findSigmaActivity(0.0, Ks, 4)[1], 0.0)
            self.assertEqual(fb.findSigmaActivity(0.0, Ks, 0.5)[1], -np.inf)
        # and the closed form matches finite differences away from zero
        rng = np.random.default_rng(11)
        for _ in range(50):
            F, Ks, n, h = rng.uniform(0.01, 5), rng.uniform(0.1, 5), rng.choice([0.5, 1, 2, 4]), 1e-6
            fd = (fb.findSigmaActivity(F + h, Ks, n)[0] - fb.findSigmaActivity(F - h, Ks, n)[0]) / (2 * h)
            self.assertAlmostEqual(fb.findSigmaActivity(F, Ks, n)[1] / fd, 1.0, places=5)

    def test_collapse_rests_on_synthesis_proportional_to_growth(self):
        # growth mode: synthesis = mu(U)*P*(...) so it vanishes when growth stops ...
        p = dict(mode="growth", P=10.0, phi0=0.01, phim=0.05, k=0.1, K=1.0, mu=0.0,
                 Ks=1.0, n=2, L=0.25, delta=0.001)
        self.assertEqual(fb.findProduction(1.0, 0.1, p), 0.0)
        # ... and at the collapsed state of the bistable case synthesis is ~1000x below
        # the served state, because growth (and with it synthesis) has nearly stopped
        q = dict(p, mu=0.1, Ktox=5.0, ntox=2, Uscan=1e5)
        eq = fb.findEquilibria(q)
        self.assertEqual(len(eq), 3)
        syn = [fb.findProduction(U, fb.findFree(C, U, q["K"]), q) for U, C in eq]
        self.assertLess(syn[2], 1e-2 * syn[0])
        # fixed mode, the opposite assumption: no ceiling as growth stops
        f = dict(BASE)
        caps = [fb.findCapMax(dict(f, mu=m)) for m in [0.5, 0.05, 0.005]]
        self.assertAlmostEqual(caps[2] / caps[0], 100.0, places=9)
        with self.assertRaises(ValueError):
            fb.findProduction(1.0, 0.1, dict(f, Ktox=5.0))
        with self.assertRaises(ValueError):
            fb.findLoopTerms(1.0, 0.5, q)


if __name__ == "__main__":
    unittest.main()
