"""tests for laws.py.

    cd experiment_grounded_rebuild_20261010/laws && python -m unittest test_laws -v

every parameter is a dimensionless placeholder. a passing test means a stated
algebraic or numerical property holds for the instances given; it says nothing
about E. coli parameter values.
"""
import unittest
from fractions import Fraction

import numpy as np

import laws as L


def randomLP(rng, nC=3, nM=3, nR=6, serialProb=0.5):
    """random route LP with serial and parallel routes, yields and returns."""
    machines = {f'm{k}': rng.uniform(0.2, 3.0) for k in range(nM)}
    cl = [f'c{k}' for k in range(nC)]
    routes = []
    for k in range(nR):
        c = cl[k % nC]
        nUse = 2 if rng.random() < serialProb else 1
        ms = rng.choice(list(machines), size=nUse, replace=False)
        routes.append({'client': c, 'a': {m: rng.uniform(0.2, 3.0) for m in ms},
                       'rho': rng.uniform(0.2, 1.0), 'b': rng.uniform(0, 0.5)})
    supply = {c: rng.uniform(0.5, 4.0) for c in cl}
    demand = {c: rng.uniform(0.05, 0.8) for c in cl}
    return L.RouteLP(routes, machines, supply, demand)


class TestTimeToDepletion(unittest.TestCase):
    def test_symbolic_solution_and_crossing(self):
        r = L.checkCrossingSymbolic()
        self.assertTrue(r['ode_solution_matches'])
        self.assertTrue(r['crossing_hits_req'])

    def test_failure_time_inside_window_random(self):
        # P(t), lam(t), Nreq(t) wander inside their bounds; the first quota
        # violation must lie in [t_lo, t_up]
        rng = np.random.default_rng(11)
        n = 0
        for _ in range(150):
            lamLo = rng.uniform(0.2, 1.0); lamHi = lamLo * rng.uniform(1.0, 2.0)
            reqLo = rng.uniform(0.5, 1.0); reqHi = reqLo * rng.uniform(1.0, 1.3)
            pHi = rng.uniform(0.1, 0.95) * lamLo * reqLo
            pLo = rng.uniform(0.0, 1.0) * pHi
            n0 = reqHi * rng.uniform(1.1, 3.0)
            w = rng.uniform(0.5, 5.0, 3); ph = rng.uniform(0, 6.3, 3)
            pF = lambda t: pLo + (pHi - pLo) * 0.5 * (1 + np.sin(w[0] * t + ph[0]))
            lF = lambda t: lamLo + (lamHi - lamLo) * 0.5 * (1 + np.sin(w[1] * t + ph[1]))
            rF = lambda t: reqLo + (reqHi - reqLo) * 0.5 * (1 + np.sin(w[2] * t + ph[2]))
            tLo, tUp = L.findFailureWindow(n0, pLo, pHi, lamLo, lamHi, reqLo, reqHi)
            self.assertTrue(np.isfinite(tUp))
            tf, _, _ = L.simulateScalarQuota(n0, pF, lF, rF, 1.2 * tUp + 1.0)
            grid = (1.2 * tUp + 1.0) / 4000
            self.assertGreaterEqual(tf, tLo - grid)
            self.assertLessEqual(tf, tUp + grid)
            n += 1
        self.assertEqual(n, 150)

    def test_upper_bound_is_attained_by_extreme_inputs(self):
        tUp = L.findCrossingTime(2.0, 0.3, 0.5, 1.0)
        tf, _, _ = L.simulateScalarQuota(2.0, lambda t: 0.3, lambda t: 0.5,
                                         lambda t: 1.0, 2 * tUp, nPts=200001)
        self.assertAlmostEqual(tf, tUp, delta=2 * tUp / 200000)

    def test_cotranslational_floor_excludes_failure(self):
        # if sigma*s >= lamHi*reqHi the quota can never fail, whatever machines do
        sigS, lamHi, req = 0.85, 0.8, 1.0
        self.assertGreaterEqual(sigS, lamHi * req)   # protected instance
        self.assertEqual(L.findCrossingTime(1.5, sigS, lamHi, req), np.inf)
        tLo, _ = L.findFailureWindow(1.5, sigS, 2.0, 0.5, lamHi, req, req)
        self.assertEqual(tLo, np.inf)
        # just below the floor the lower bound becomes finite
        tLo2, _ = L.findFailureWindow(1.5, 0.75, 2.0, 0.5, lamHi, req, req)
        self.assertTrue(np.isfinite(tLo2))

    def test_moving_quota_face_term_is_needed(self):
        st, mv = L.findMovingQuotaP()
        self.assertAlmostEqual(st, 1.2)
        self.assertAlmostEqual(mv, 3.3)
        self.assertLess(L.findMovingQuotaRun(st), -1.0)        # static condition fails
        self.assertGreaterEqual(L.findMovingQuotaRun(mv), -1e-9)  # face condition holds


class TestRouteFeasibility(unittest.TestCase):
    def test_serial_route_breaks_hall_union(self):
        lp = L.findSerialToyLP(1.5)
        self.assertFalse(lp.findFeasible()[0])
        self.assertTrue(L.checkHallUnion(lp))     # casual Hall test wrongly passes

    def test_parallel_route_breaks_bottleneck_min(self):
        lp = L.findParallelToyLP(1.5)
        self.assertTrue(lp.findFeasible()[0])
        self.assertFalse(L.checkBottleneckMin(lp))  # serial-min rule wrongly fails

    def test_serial_certificate_exact(self):
        # certificate z = 1 (per native chain), w_A = 1 (per machine-time), others 0.
        # route: w_A a_A + w_B a_B + u(1-b) - z rho = 1 + 0 + 0 - 1 = 0 >= 0
        # value: z d - w_A M_A - w_B M_B - u e = 3/2 - 1 = 1/2 > 0
        z, wA, wB, u = Fraction(1), Fraction(1), Fraction(0), Fraction(0)
        aA, aB, rho, d, MA, MB, e = (Fraction(1), Fraction(1), Fraction(1),
                                     Fraction(3, 2), Fraction(1), Fraction(1), Fraction(10))
        self.assertGreaterEqual(wA * aA + wB * aB + u - z * rho, 0)
        self.assertEqual(z * d - wA * MA - wB * MB - u * e, Fraction(1, 2))
        # the solver finds an equivalent certificate
        y, val = lp_cert = L.findSerialToyLP(1.5).findFarkasCertificate()
        self.assertIsNotNone(y)
        red, gap = L.findSerialToyLP(1.5).checkCertificate(y)
        self.assertTrue(np.all(red >= -1e-12))
        self.assertGreater(gap, 0)

    def test_farkas_alternative_random(self):
        # exactly one of: feasible point, or certificate (checked independently)
        rng = np.random.default_rng(5)
        nFeas = nCert = 0
        for _ in range(300):
            lp = randomLP(rng)
            feas, x = lp.findFeasible()
            y, _ = lp.findFarkasCertificate()
            self.assertNotEqual(feas, y is not None)
            if feas:
                G, h, _ = lp.buildRows()
                self.assertTrue(np.all(G @ x <= h + 1e-8))
                self.assertTrue(np.all(x >= -1e-12))
                nFeas += 1
            else:
                red, gap = lp.checkCertificate(y)
                self.assertTrue(np.all(red >= -1e-9))
                self.assertGreater(gap, 1e-9)
                nCert += 1
        self.assertGreater(nFeas, 20)
        self.assertGreater(nCert, 20)

    def test_pooling_is_a_strict_relaxation(self):
        pl = L.findPoolingToyLP()
        self.assertFalse(pl.findFeasible()[0])
        self.assertTrue(pl.findPooled().findFeasible()[0])
        rng = np.random.default_rng(6)
        for _ in range(200):
            lp = randomLP(rng)
            if lp.findFeasible()[0]:
                self.assertTrue(lp.findPooled().findFeasible()[0])

    def test_gale_cut_condition_exact_for_single_machine_routes(self):
        # special case where a Hall/Gale cut condition is exact
        rng = np.random.default_rng(7)
        for _ in range(200):
            cl = [f'c{k}' for k in range(4)]; ms = [f'm{k}' for k in range(3)]
            M = {m: rng.uniform(0.2, 2.0) for m in ms}
            tau = {m: rng.uniform(0.3, 2.0) for m in ms}
            qual = {c: set(rng.choice(ms, size=rng.integers(1, 3), replace=False))
                    for c in cl}
            demand = {c: rng.uniform(0.1, 1.5) for c in cl}
            routes = [{'client': c, 'a': {m: tau[m]}, 'rho': 1.0}
                      for c in cl for m in sorted(qual[c])]
            lp = L.RouteLP(routes, M, {c: 100.0 for c in cl}, demand)
            self.assertEqual(lp.findFeasible()[0], L.checkGaleSubsets(demand, M, tau, qual))

    def test_time_unit_invariance(self):
        # rescale time by c: a -> c a, fluxes -> flux / c, pools unchanged
        rng = np.random.default_rng(8)
        for _ in range(50):
            lp = randomLP(rng)
            c = rng.uniform(0.01, 100.0)
            r2 = [dict(r, a={m: v * c for m, v in r['a'].items()}) for r in lp.routes]
            lp2 = L.RouteLP(r2, lp.machines, {k: v / c for k, v in lp.supply.items()},
                            {k: v / c for k, v in lp.demand.items()})
            self.assertEqual(lp.findFeasible()[0], lp2.findFeasible()[0])
            t1, t2 = lp.findMargin()[0], lp2.findMargin()[0]
            self.assertAlmostEqual(t1, t2, delta=1e-7 * max(1.0, t1))

    def test_margin_monotone_in_every_pool(self):
        rng = np.random.default_rng(9)
        for _ in range(100):
            lp = randomLP(rng)
            t0 = lp.findMargin()[0]
            m = rng.choice(list(lp.machines))
            M2 = dict(lp.machines); M2[m] *= rng.uniform(1.0, 3.0)
            t1 = L.RouteLP(lp.routes, M2, lp.supply, lp.demand).findMargin()[0]
            self.assertGreaterEqual(t1, t0 - 1e-9)

    def test_realized_steady_state_is_lp_feasible(self):
        # necessity: route fluxes of a kinetic steady state satisfy the LP rows
        rng = np.random.default_rng(10)
        for _ in range(200):
            kf, kl = rng.uniform(0, 1), rng.uniform(0.01, 1)
            kon = rng.uniform(0.1, 10, 2); M = rng.uniform(0.1, 5, 2)
            tau = rng.uniform(0.1, 3, 2); pN = rng.uniform(0, 0.9, 2)
            pU = rng.uniform(0, 1, 2) * (1 - pN)
            e = rng.uniform(0.1, 3)
            st = L.findSteadySingleClient(e, kf, kl, kon, M, tau, pN, pU)
            lp = L.findSingleClientRouteLP(e, st['Y'] * e, kf, kon, M, tau, pN, pU,
                                           st['U'], ['G', 'H'])
            x = np.concatenate([[kf * st['U']], st['J']])
            G, h, _ = lp.buildRows()
            self.assertTrue(np.all(G @ x <= h + 1e-9 * (1 + np.abs(h))))


class TestCapacityVersusRealized(unittest.TestCase):
    def test_sign_law_symbolic(self):
        self.assertTrue(L.checkSignLawSymbolic())

    def test_sign_law_transfers_to_finite_pool_single_machine(self):
        # one machine plus spontaneous routes: sign dY/dM = sign(committed - Y)
        rng = np.random.default_rng(12)
        agree = 0
        for _ in range(500):
            kf, kl = rng.uniform(0, 1), rng.uniform(0.01, 1)
            kon, M, tau = rng.uniform(0.1, 10), rng.uniform(0.1, 5), rng.uniform(0.1, 3)
            pN = rng.uniform(0, 0.95); pU = rng.uniform(0, 0.99) * (1 - pN)
            e = rng.uniform(0.1, 3)
            f = lambda MM: L.findSteadySingleClient(e, kf, kl, [kon], [MM], [tau],
                                                    [pN], [pU])['Y']
            Y = f(M)
            committed = pN / (1 - pU)
            if abs(committed - Y) < 1e-4:
                continue
            dY = (f(M * (1 + 1e-6)) - f(M * (1 - 1e-6)))
            self.assertEqual(np.sign(dY), np.sign(committed - Y))
            agree += 1
        self.assertGreater(agree, 450)

    def test_two_machines_sign_transfer_recorded_not_claimed(self):
        # with two finite pools the transfer is not proved; only count exceptions
        bad, tot = L.findTwoMachineSignExceptions()
        self.assertGreater(tot, 300)
        self.assertLessEqual(bad, tot)

    def test_s1_capacity_up_realized_down(self):
        rows = L.findS1Comparison()
        th = [r['margin_LP'] for r in rows]; Y = [r['Y_realized'] for r in rows]
        self.assertTrue(all(b >= a - 1e-9 for a, b in zip(th, th[1:])))
        self.assertTrue(all(b < a for a, b in zip(Y, Y[1:])))
        self.assertTrue(all(r['committed_H'] < r['Y_realized'] for r in rows))
        self.assertGreater(Y[0] / Y[-1], 3.0)

    def test_quota_failure_at_vanishing_occupancy(self):
        for eps in (0.1, 0.03, 0.01):
            st = L.findLowOccupancyFailure(eps)
            self.assertLess(st['Y'], eps)
            self.assertLess(np.max(st['occupancy']), eps)


class TestInductionRace(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.r = L.findRaceSummary()
        cls.f = L.findRaceSummary(dict(L.RACE, forceU=1.0))

    def test_bound_symbolic(self):
        r = L.checkRaceBoundSymbolic()
        self.assertTrue(r['solves_comparison_ode'])
        self.assertTrue(r['initial_value'])

    def test_machine_ceiling_any_regulation(self):
        self.assertLess(L.findRandomControlCeiling(), 1e-8)

    def test_trajectory_below_bound_and_conserved(self):
        for r in (self.r, self.f):
            self.assertLessEqual(r['max_N_minus_bound'], 1e-8)
            self.assertLess(r['mass_error'], 1e-9)
            self.assertGreaterEqual(r['min_state'], -1e-10)
            self.assertGreaterEqual(r['min_free_machine'], -1e-10)

    def test_static_feasible_yet_guaranteed_transient_violation(self):
        r = self.r
        self.assertGreater(r['steady_bound_full_induction'], r['Nreq'])  # static capacity suffices
        self.assertGreater(r['Nss_post'], r['Nreq'])    # closed-loop equilibrium is viable
        self.assertGreater(r['N0'], r['Nreq'])          # start is viable
        self.assertLess(r['min_bound'], r['Nreq'])      # bound guarantees a violation
        self.assertIsNotNone(r['first_violation_t'])
        self.assertLessEqual(r['first_violation_t'], r['bound_guarantees_violation_by_t'] + 0.011)

    def test_post_shock_equilibrium_locally_stable(self):
        ev, resid = L.findRaceEquilibriumEigen()
        self.assertLess(resid, 1e-9)
        self.assertLess(np.max(ev.real), 0)

    def test_violation_survives_maximal_instant_induction(self):
        self.assertIsNotNone(self.f['first_violation_t'])
        self.assertGreater(self.f['Nss_post'], self.f['Nreq'])


class TestCoRegulation(unittest.TestCase):
    def test_capacity_set_up_closed(self):
        # LP feasible at regulator u implies feasible at u = 1 (pools only grow)
        rng = np.random.default_rng(14)
        for _ in range(100):
            lp = randomLP(rng)
            base = {m: 0.3 * v for m, v in lp.machines.items()}
            gain = {m: rng.uniform(0, 2) * v for m, v in lp.machines.items()}
            feas = [L.RouteLP(lp.routes, {m: base[m] + gain[m] * u for m in base},
                              lp.supply, lp.demand).findFeasible()[0]
                    for u in np.linspace(0, 1, 6)]
            if any(feas):
                self.assertTrue(feas[-1])

    def test_realized_yield_has_interior_optimum(self):
        us = np.linspace(0, 1, 101)
        ys = np.array([L.findCoRegulatedYield(u) for u in us])
        j = int(np.argmax(ys))
        self.assertTrue(0 < j < 100)
        self.assertGreater(ys[j] - max(ys[0], ys[-1]), 0.05)


class TestCodon(unittest.TestCase):
    def test_cotranslational_symbolic(self):
        self.assertTrue(L.checkCotranslationalSymbolic())

    def test_two_client_cross_effect_and_ode_agreement(self):
        a, b = L.findTwoClientSteady(0.0), L.findTwoClientSteady(0.3)
        self.assertLess(a['resid'], 1e-12)
        self.assertGreater(b['Y'][0], a['Y'][0])   # client 1 gains
        self.assertLess(b['Y'][1], a['Y'][1])      # client 2 loses
        self.assertGreater(b['F'][1], a['F'][1])   # through more free H
        for sg, st in ((0.0, a), (0.3, b)):
            od = L.findTwoClientOde(sg)
            self.assertLess(od['massErr'], 1e-8)
            self.assertGreaterEqual(od['minState'], -1e-10)
            self.assertTrue(np.allclose(od['Y'], st['Y'], atol=1e-6))

    def test_margin_sensitivity_binding_vs_slack(self):
        r = L.findMarginSensitivity()
        self.assertAlmostEqual(r['fd_G'], r['pred_G'], delta=1e-5)
        self.assertAlmostEqual(r['fd_K'], 0.0, delta=1e-9)
        self.assertGreater(r['w_G'], 0)
        self.assertEqual(r['w_K'], 0)

    def test_mistranslation_first_order_bound(self):
        rng = np.random.default_rng(15)
        for _ in range(500):
            n = rng.integers(0, 30, 5); eps = 10 ** rng.uniform(-6, -2, 5)
            exact, bound = L.findMistranslationLoadBound(n, eps)
            self.assertLessEqual(exact, bound + 1e-15)


class TestEnergyOptional(unittest.TestCase):
    def test_front_exists_only_with_energy_row(self):
        budgets = [0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 10.0]
        front, free = L.findEnergyFront(budgets)
        self.assertTrue(all(b >= a - 1e-9 for a, b in zip(front, front[1:])))
        self.assertAlmostEqual(front[-1], free, places=7)
        self.assertLess(front[0], 0.5 * free)


if __name__ == '__main__':
    unittest.main()
