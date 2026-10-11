"""tests and negative controls for model.py.  python -m unittest test_model -v"""
import unittest

import numpy as np

import model as m

MEAS = m.loadMeasured()
MU, M0 = MEAS['mu'], MEAS['pools']


def clients(seed, n=1, fH_zero=True):
    rng = np.random.default_rng([99, seed])
    return [m.sampleClient(rng, MEAS, fH_zero) for _ in range(n)]


def odeEnd(cl, M, t=400.0, y0=None):
    sol = m.integrate(cl, M, MU, t / MU, y0)
    return sol, sol.y[:, -1].reshape(len(cl), m.NS)


class TestInputs(unittest.TestCase):
    def test_pools_are_the_measured_values(self):
        # WO-07 MATCHED values, Schmidt 2016 BW25113 glucose exponential
        self.assertAlmostEqual(M0[m.K], 10.8803)
        self.assertAlmostEqual(M0[m.G], 0.949851)
        self.assertAlmostEqual(M0[m.B], 0.0116582)
        self.assertAlmostEqual(M0[m.T], 14.5782)
        self.assertAlmostEqual(M0[m.H], 0.612733)
        self.assertAlmostEqual(MEAS['mu'], 0.000161111)
        # GroES7 exceeds GroEL14, so GroEL14 is the limiting chaperonin unit
        self.assertGreater(MEAS['GroES7'], M0[m.G])


class TestExactness(unittest.TestCase):
    def test_closed_form_matches_flux_yield_and_ode(self):
        for seed in range(3):
            cl = clients(seed, 2, fH_zero=False)
            st = m.solveSteady(cl, M0, MU)
            self.assertTrue(st['ok'])
            _, y = odeEnd(cl, M0)
            for i, r in enumerate(st['clients']):
                self.assertAlmostEqual(r['Y'], r['Y_formula'], places=12)
                for k, col in [('U', 0), ('CG', 3), ('CH', 4), ('A', 6)]:
                    self.assertLess(abs(r[k] - y[i, col]),
                                    1e-6 * max(r[k], 1e-12) + 1e-13)

    def test_route_criterion_theorem(self):
        # at fixed free pools: increasing capture of route x raises Y iff
        # pi_x = a_x/(1-b_x) exceeds Ymarg = P/(Q + 2 k' U)
        checked = 0
        for seed in range(40):
            c = clients(seed, 1, fH_zero=bool(seed % 2))[0]
            F = m.solveSteady([c], M0, MU)['F']
            r = m.clientSteady(c, F, MU)
            for j in range(4):
                gap = r['pi'][j] - r['Ymarg']
                if abs(gap) < 1e-6:
                    continue
                c2 = m.copyClient(c)
                c2['kon'][j] *= 1 + 1e-7
                dY = m.clientSteady(c2, F, MU)['Y'] - r['Y']
                self.assertEqual(np.sign(dY), np.sign(gap))
                checked += 1
        self.assertGreater(checked, 100)

    def test_faster_folding_raises_yield_at_fixed_pools(self):
        for seed in range(30):
            c = clients(seed)[0]
            F = m.solveSteady([c], M0, MU)['F']
            y0 = m.clientSteady(c, F, MU)['Y']
            y1 = m.clientSteady(m.copyClient(c, kf=1.3 * c['kf']), F, MU)['Y']
            self.assertGreater(y1, y0)

    def test_pure_holdase_cannot_raise_yield(self):
        # a route with no productive exit (a_x = 0) never raises Y
        c = clients(5)[0]
        c.update(fT=0.0, kh=0.0)
        F = m.solveSteady([c], M0, MU)['F']
        y0 = m.clientSteady(c, F, MU)['Y']
        c2 = m.copyClient(c)
        c2['kon'][0] *= 3
        self.assertLessEqual(m.clientSteady(c2, F, MU)['Y'], y0)


class TestConservation(unittest.TestCase):
    def test_mass_conservation(self):
        for seed in range(3):
            cl = clients(seed, 2, fH_zero=False)
            sol, y = odeEnd(cl, M0, t=50.0)
            for i, c in enumerate(cl):
                total = y[i].sum()
                self.assertLess(abs(total - c['s'] * sol.t[-1]),
                                1e-8 * c['s'] * sol.t[-1])

    def test_positivity_and_capacity(self):
        for seed in range(3):
            cl = clients(seed, 2, fH_zero=False)
            sol, _ = odeEnd(cl, M0, t=50.0)
            self.assertGreaterEqual(sol.y.min(), -1e-12)
            st = m.solveSteady(cl, M0, MU)
            self.assertTrue(np.all(st['F'] >= 0))
            occ = sum(r['occ'] for r in st['clients'])
            self.assertTrue(np.all(occ <= M0 * (1 + 1e-9)))
            for r in st['clients']:
                for k in ['U', 'CT', 'CK', 'CG', 'CH', 'CV', 'A', 'CB']:
                    self.assertGreaterEqual(r[k], 0.0)
                self.assertTrue(0.0 <= r['Y'] <= 1.0)

    def test_steady_state_independent_of_start(self):
        cl = clients(7, 2, fH_zero=False)
        y0 = np.zeros(2 * m.NS)
        y0[0] = y0[m.NS] = 5.0            # load free U in both clients
        _, a = odeEnd(cl, M0)
        _, b = odeEnd(cl, M0, y0=y0)
        self.assertLess(np.max(np.abs(a[:, :8] - b[:, :8])), 1e-8)


class TestNegativeControls(unittest.TestCase):
    def test_htpg_route_disabled_means_no_htpg_effect(self):
        for seed in range(5):
            cl = clients(seed, 2, fH_zero=True)
            for c in cl:
                c['kon'][m.CAPTURE.index(m.H)] = 0.0
            y = m.solveSteady(cl, M0, MU)['clients'][0]['Y']
            yk = m.solveSteady(*m.perturb(cl, M0, knock=(m.H,)), MU)
            self.assertEqual(y, yk['clients'][0]['Y'])

    def test_groel_rescue_requires_direct_bypass(self):
        for seed in range(10):
            cl = clients(seed)
            ys = [m.solveSteady(*m.perturb(cl, M0, knock=(m.T, m.K),
                                           scale={m.G: f}, serial=True),
                                MU)['clients'][0]['Y'] for f in (1, 5)]
            self.assertEqual(ys[0], ys[1])           # serial: inexpressible
            yb = [m.solveSteady(*m.perturb(cl, M0, knock=(m.T, m.K),
                                           scale={m.G: f}),
                                MU)['clients'][0]['Y'] for f in (1, 5)]
            self.assertNotEqual(yb[0], yb[1])        # bypass: it can move

    def test_same_total_swap(self):
        cl = clients(3, 2, fH_zero=False)
        swap = M0.copy()
        swap[m.K] -= 0.4
        swap[m.G] += 0.4
        r = lambda M, **kw: m.solveSteady(cl, M, MU, **kw)['clients'][0]['Y']
        # identical up to floating-point summation order of the pool total
        self.assertAlmostEqual(r(M0, pooled=True, Mref=M0),
                               r(swap, pooled=True, Mref=M0), places=13)
        self.assertNotAlmostEqual(r(M0), r(swap), places=9)

    def test_pooled_equals_resolved_at_vanishing_load(self):
        cl = [m.copyClient(c, s=c['s'] * 1e-9) for c in clients(4, 2)]
        a = m.solveSteady(cl, M0, MU)['clients'][0]['Y']
        b = m.solveSteady(cl, M0, MU, pooled=True, Mref=M0)['clients'][0]['Y']
        self.assertAlmostEqual(a, b, places=8)

    def test_no_clpb_quota(self):
        # demand is set by the client's synthesis share alone; removing ClpB
        # leaves every client parameter unchanged and aggregation is an
        # outcome, not a constraint that must be met
        cl = clients(2, 2)
        c2, M2 = m.perturb(cl, M0, knock=(m.B,))
        for a, b in zip(cl, c2):
            self.assertEqual(a['s'], b['s'])
            self.assertEqual(a['kf'], b['kf'])
        st = m.solveSteady(c2, M2, MU)
        self.assertTrue(st['ok'])
        for r in st['clients']:
            self.assertGreater(r['Y'], 0.0)
            self.assertGreater(r['A'], 0.0)

    def test_htpG_loss_with_no_productive_exit_raises_single_client_yield(self):
        # this sign is ENCODED by fH = 0, not a prediction
        for seed in range(10):
            cl = clients(seed, 1, fH_zero=True)
            y = m.solveSteady(cl, M0, MU)['clients'][0]['Y']
            yk = m.solveSteady(*m.perturb(cl, M0, knock=(m.H,)),
                               MU)['clients'][0]['Y']
            self.assertGreater(yk, y)


if __name__ == '__main__':
    unittest.main()
