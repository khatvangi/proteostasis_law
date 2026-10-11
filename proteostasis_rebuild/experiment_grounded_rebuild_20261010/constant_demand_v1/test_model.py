"""tests for model.py.   python -m unittest test_model -v
all parameter values are hypothetical; the tests check logic, not biology."""
import unittest

import numpy as np

import model as m

RNG = np.random.default_rng(20261010)


class TestScalar(unittest.TestCase):
    def test_steady_state_is_exact_and_attracting(self):
        # E2/E3: N(t) -> P/lam from any start, and P/lam is a fixed point
        for _ in range(200):
            P, lam, N0 = RNG.uniform(1, 100), RNG.uniform(0.001, 1), RNG.uniform(0, 1e4)
            Ns = m.steadyLevel(P, lam)
            self.assertAlmostEqual(P - lam * Ns, 0.0, places=9)
            self.assertAlmostEqual(m.levelAt(60 / lam, P, lam, N0), Ns,
                                   delta=1e-12 * max(Ns, N0) + 1e-6)

    def test_condition_is_necessary_and_sufficient(self):
        # E6 against an independent brute force: is there any P on a fine grid
        # of [0, Pmax] whose steady state meets D?  draws avoid the knife edge
        for _ in range(500):
            a, lam, Pmax = RNG.uniform(1, 1e3), RNG.uniform(1e-3, 1), RNG.uniform(1, 1e3)
            D = a * Pmax / lam * RNG.choice([RNG.uniform(0.1, 0.99), RNG.uniform(1.01, 3)])
            grid = np.linspace(0, Pmax, 2001)
            brute = bool(np.any(a * m.steadyLevel(grid, lam) >= D))
            self.assertEqual(brute, m.scalarFeasible(D, a, Pmax, lam))

    def test_required_supply_scales_with_turnover_not_with_service(self):
        # E5: catalysis does not consume the catalyst; cost is lam*Nreq.
        D, a = 1.8e6, 600.0
        self.assertAlmostEqual(m.requiredSupply(D, a, 0.02), 60.0)
        self.assertAlmostEqual(m.requiredSupply(D, a, 0.04), 120.0)  # x2 turnover
        self.assertEqual(m.requiredSupply(D, a, 0.0), 0.0)            # no turnover
        # E9: only if each service event inactivates eps molecules does D enter directly
        self.assertAlmostEqual(m.requiredSupply(D, a, 0.02, eps=1e-5), 60.0 + 18.0)

    def test_route_feedback_breaks_sufficiency_not_necessity(self):
        # E10: P(N) = Pmax - beta N. Q* = a Pmax/(lam+beta) < a Pmax/lam
        a, lam, Pmax, beta = 600.0, 0.02, 95.0, 0.01
        Qstar = a * m.feedbackLevel(Pmax, beta, lam)
        D = 0.5 * (Qstar + a * Pmax / lam)            # between the two
        self.assertTrue(m.scalarFeasible(D, a, Pmax, lam))   # condition passes
        self.assertLess(Qstar, D)                            # yet demand unmet
        for _ in range(200):                                 # still necessary
            b = RNG.uniform(0, 1)
            self.assertLessEqual(a * m.feedbackLevel(Pmax, b, lam), a * Pmax / lam + 1e-9)

    def test_installed_capacity_decays_to_sustainable(self):
        # E8: installed N0 above Nreq does not mean sustainable
        Nreq, lam = 3000.0, 0.02
        Ns = m.steadyLevel(m.realizedRouting(m.ROUTING, m.POOR)['P'], lam)
        t = m.timeToShortfall(5000.0, Nreq, Ns, lam)
        self.assertTrue(np.isfinite(t) and t > 0)
        P = lam * Ns
        self.assertAlmostEqual(m.levelAt(t, P, lam, 5000.0), Nreq, places=6)
        self.assertEqual(m.timeToShortfall(5000.0, Nreq, 3500.0, lam), np.inf)


class TestRouting(unittest.TestCase):
    def test_same_counts_poor_routing_fails_without_saturation(self):
        p = m.ROUTING
        Preq = m.requiredSupply(m.D_SCALAR, m.A_SPEC, m.LAM)
        Pmax = m.potentialSupply(p)            # depends on counts, not on bK, bG
        good, poor = m.realizedRouting(p, m.GOOD), m.realizedRouting(p, m.POOR)
        self.assertTrue(m.scalarFeasible(m.D_SCALAR, m.A_SPEC, Pmax, m.LAM))
        self.assertGreaterEqual(good['P'], Preq)
        self.assertLess(poor['P'], Preq)
        # the shortfall is not a saturation effect: both machines mostly free
        self.assertLess(max(poor['occK'], poor['occG']), 0.10)
        # doubling both installed pools still does not rescue the poor case
        p2 = dict(p, Ktot=2 * p['Ktot'], Gtot=2 * p['Gtot'])
        self.assertLess(m.realizedRouting(p2, m.POOR)['P'], Preq)

    def test_realized_never_exceeds_potential(self):
        # realized fluxes are one feasible point of the LP, so P_real <= Pmax
        for _ in range(100):
            p = dict(m.ROUTING, kf=RNG.uniform(0.05, 2), ka=RNG.uniform(0.05, 2),
                     yK=RNG.uniform(0, 1), yG=RNG.uniform(0.3, 1),
                     Ktot=RNG.uniform(5, 300), Gtot=RNG.uniform(5, 100))
            b = dict(bK=10 ** RNG.uniform(-4, -1), bG=10 ** RNG.uniform(-4, -1))
            r = m.realizedRouting(p, b)
            self.assertLessEqual(r['P'], m.potentialSupply(p) + 1e-7)
            self.assertLessEqual(max(r['occK'], r['occG']), 1 + 1e-9)


class TestVector(unittest.TestCase):
    def test_lp_reduces_to_scalar_condition(self):
        # one client, one function: LP feasibility == D <= a Pmax / lam
        clients, machines = m.scalarClient(m.ROUTING)
        Pmax = m.potentialSupply(m.ROUTING)
        for scale in (0.5, 0.99, 1.01, 2.0):
            D = scale * m.A_SPEC * Pmax / m.LAM
            lp = m.buildLP(clients, machines, [[m.A_SPEC]], [D])
            self.assertEqual(m.feasible(lp)[0], m.scalarFeasible(D, m.A_SPEC, Pmax, m.LAM))

    def test_hall_counting_passes_but_serial_route_is_infeasible(self):
        # client 1 must go K then G; client 2 only G. capacity 1 each, dwell 1 min.
        clients = [dict(name='c1', lam=1.0, length=1.0, routes=[(1.0, {'K': 1.0, 'G': 1.0})]),
                   dict(name='c2', lam=1.0, length=1.0, routes=[(1.0, {'G': 1.0})])]
        machines = {'K': dict(fixed=1.0), 'G': dict(fixed=1.0)}
        need = [1.0, 0.5]                          # required productive flux
        self.assertTrue(m.hallSubsetCheck(need, [{'K', 'G'}, {'G'}], {'K': 1.0, 'G': 1.0}))
        lp = m.buildLP(clients, machines, np.eye(2), need)
        self.assertFalse(m.feasible(lp)[0])        # G would carry 1.5 > 1
        u, w = m.certificate(lp)
        self.assertTrue(m.checkCertificate(lp, u, w))

    def test_vector_example_feasibility_and_certificates(self):
        D0 = np.array([1.8e6, 0.9e6])
        cl, mc, A = m.vectorExample(60000.0, 0.02)
        lp = m.buildLP(cl, mc, A, D0, 60000.0)
        ok, x = m.feasible(lp)
        self.assertTrue(ok)
        self.assertTrue(np.all(lp['F'] @ x >= D0 * (1 - 1e-9)))
        self.assertTrue(np.all(lp['G'] @ x <= lp['h'] + 1e-6))
        self.assertIsNone(m.certificate(lp))       # no certificate when feasible
        lp15 = m.buildLP(cl, mc, A, 1.5 * D0, 60000.0)
        self.assertFalse(m.feasible(lp15)[0])
        self.assertTrue(m.checkCertificate(lp15, *m.certificate(lp15)))

    def test_machine_maintenance_can_decide_feasibility(self):
        # same demand and budget; only the cost of sustaining machines differs
        D0 = np.array([1.8e6, 0.9e6])
        verdict = {}
        for lam_m in (0.02, 0.0):
            cl, mc, A = m.vectorExample(48000.0, lam_m)
            verdict[lam_m] = m.feasible(m.buildLP(cl, mc, A, D0, 48000.0))[0]
        self.assertFalse(verdict[0.02])
        self.assertTrue(verdict[0.0])

    def test_tampered_certificate_is_rejected(self):
        cl, mc, A = m.vectorExample(60000.0, 0.02)
        lp = m.buildLP(cl, mc, A, 1.5 * np.array([1.8e6, 0.9e6]), 60000.0)
        u, w = m.certificate(lp)
        self.assertFalse(m.checkCertificate(lp, np.zeros_like(u), w))


if __name__ == '__main__':
    unittest.main()
