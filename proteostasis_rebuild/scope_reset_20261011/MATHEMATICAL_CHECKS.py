#!/usr/bin/env python3
"""Checks for the bounded reassessment. All numerical witnesses are illustrative.
Run: python3 MATHEMATICAL_CHECKS.py > CHECK_RESULTS.txt 2>&1
Uses the standard library; adds symbolic checks if SymPy is installed.
No biological fitting, files written, or old ensembles rerun.
"""
import random
import sys
import unittest
from fractions import Fraction as Q
try:
    import sympy as sp
except ImportError:
    sp = None


def rhs(x, p):
    U, C, A, N, F = x
    s, q, ku, r, f, kf, du, dc, da, dn, mu, kon, ka, alpha = p
    capture, agg = kon*F*U, ka*U*U
    return (q*s+ku*N+r*C-capture-kf*U-du*U-agg-mu*U,
            capture-(r+f+dc+mu)*C,
            agg-(da+mu)*A,
            (1-q)*s+kf*U+f*C-(ku+dn+mu)*N,
            alpha-mu*F-capture+(r+f+dc)*C)


def cases():
    rng = random.Random(20261011)
    for _ in range(100):
        x = [Q(rng.randrange(21), 7) for _ in range(5)]
        p = [Q(rng.randrange(21), 13) for _ in range(14)]
        p[1] = Q(rng.randrange(14), 13)
        yield x, p


class Accounting(unittest.TestCase):
    def test_01_client_conservation_exact(self):
        for x, p in cases():
            U,C,A,N,F = x
            du,dc,da,dn,mu = p[6:11]
            self.assertEqual(sum(rhs(x,p)[:4]), p[0]-mu*sum(x[:4])-du*U-dc*C-da*A-dn*N)

    def test_02_cumulative_sinks(self):
        for x,p in cases():
            degraded = sum(a*b for a,b in zip(p[6:10], x[:4]))
            diluted = p[10]*sum(x[:4])
            self.assertEqual(sum(rhs(x,p)[:4])+degraded+diluted, p[0])

    def test_03_machine_balance(self):
        for x,p in cases():
            dx = rhs(x,p)
            self.assertEqual(dx[1]+dx[4], p[13]-p[10]*(x[1]+x[4]))

    def test_04_aggregation_transfers_not_creates(self):
        for x,p in cases():
            p0 = p.copy(); p0[12] = 0
            change = [a-b for a,b in zip(rhs(x,p),rhs(x,p0))]
            self.assertEqual(change[0], -p[12]*x[0]**2)
            self.assertEqual(change[2], p[12]*x[0]**2)
            self.assertEqual(sum(change[:4]),0)

    def test_05_counterfactual_outcome_identity(self):
        s=Q(5); pi=[Q(1,2),Q(1,3),Q(1,6)]; p=[Q(1,5),Q(4,5),Q(1,10)]
        self.assertEqual(s*sum(a*b for a,b in zip(pi,p))-s*p[0],
                         s*sum(pi[i]*(p[i]-p[0]) for i in (1,2)))

    def test_06_incremental_workload_can_have_either_sign(self):
        s,eps,w0=Q(3),Q(1,10),Q(2)
        for we,sign in [(Q(4),1),(Q(2),0),(Q(1),-1)]:
            full=s*((1-eps)*w0+eps*we)-s*w0
            self.assertEqual(full,s*eps*(we-w0))
            self.assertEqual((full>0)-(full<0),sign)

    def test_07_changed_environment_requires_full_difference(self):
        eps0,eps1=Q(1,10),Q(1,5)
        w00,we0,w01,we1=Q(1),Q(3),Q(2),Q(4)
        full=(1-eps1)*w01+eps1*we1-((1-eps0)*w00+eps0*we0)
        fixed=(eps1-eps0)*(we0-w00)
        self.assertNotEqual(full,fixed)

    def test_08_dimensions_with_chain_site_stoichiometry(self):
        # Exponents: chain count, service-site count, volume, time.
        add=lambda *ds: tuple(map(sum,zip(*ds)))
        chain=(1,0,-1,0); site=(0,1,-1,0); flux=(1,0,-1,-1)
        synthesis=flux; workload=(-1,1,0,1); tau=(0,0,0,1)
        kon=(0,-1,1,-1); ka=(-1,0,1,-1); rate=(0,0,0,-1)
        stoich=(-1,1,0,0) # one occupied site per bound chain
        self.assertEqual(add(synthesis,workload),site)
        self.assertEqual(add(kon,site,chain),flux)
        self.assertEqual(add(ka,chain,chain),flux)
        self.assertEqual(add(rate,chain),flux)
        self.assertEqual(add(flux,tau,stoich),site)
        self.assertNotEqual(flux,site)
        self.assertNotEqual(chain,flux)

    def test_09_all_zero_faces_point_inward(self):
        for x,p in cases():
            for i in range(5):
                face=x.copy(); face[i]=Q(0)
                self.assertGreaterEqual(rhs(face,p)[i],0)

    def test_10_fixed_pool_boundary(self):
        for x,p in cases():
            x[4]=0; p[13]=p[10]*x[1]
            dx=rhs(x,p)
            self.assertLessEqual(dx[1],0)
            self.assertEqual(dx[1]+dx[4],0)

    def test_11_multiclient_pool_boundary(self):
        M=Q(5); C=[Q(2),Q(3)]; rates=[Q(1,2),Q(3,4)]; mu=Q(1,10)
        self.assertEqual(M-sum(C),0)
        self.assertLess(-sum((rate+mu)*c for rate,c in zip(rates,C)),0)

    def test_12_steady_occupancy_identity(self):
        M=Q(4); U=[Q(2),Q(3)]; kon=[Q(1,3),Q(1,2)]; tau=[Q(2),Q(1)]
        F=M/(1+sum(k*u*t for k,u,t in zip(kon,U,tau)))
        capture=[k*F*u for k,u in zip(kon,U)]
        C=[b*t for b,t in zip(capture,tau)]
        self.assertEqual(F+sum(C),M)
        self.assertLessEqual(sum(C),M)
        self.assertTrue(all(b-c/t==0 for b,c,t in zip(capture,C,tau)))

    def test_13_native_identity_includes_degradation(self):
        incoming,mu,ku,dn=Q(3),Q(1,10),Q(1,5),Q(1,4)
        N=incoming/(mu+ku+dn)
        self.assertEqual(incoming,(mu+ku+dn)*N)
        self.assertNotEqual(incoming,(mu+ku)*N)
        self.assertGreaterEqual(incoming,(mu+ku+dn)*(N/2))

    def test_14_positive_dilution_bounds_mass(self):
        for x,p in cases():
            self.assertLessEqual(sum(rhs(x,p)[:4]),p[0]-p[10]*sum(x[:4]))
        # At P>s/mu the upper comparison field is negative.
        s,mu=Q(2),Q(1,5)
        self.assertLess(s-mu*(s/mu+1),0)

    def test_15_zero_dilution_removal_absent_counterexample(self):
        x=[Q(0)]*5; p=[Q(0)]*14; p[0]=Q(1); p[1]=Q(1)
        self.assertEqual(sum(rhs(x,p)[:4]),1) # P'=1, unbounded if sustained

    def test_16_signal_does_not_force_induction(self):
        # Constant b(S), a legitimate unresponsive closure, gives identical alpha.
        gamma,kt=Q(2),Q(3)
        b=lambda S: Q(4)
        self.assertEqual(kt*b(1)/gamma,kt*b(10)/gamma)
        # Responsive b still gives no production with kt=0.
        self.assertEqual(Q(0)*Q(10)/gamma,0)

    def test_17_more_capture_not_necessarily_more_native(self):
        # Unbound client folds with kf, is lost with ell, or captured for disposal h.
        # Model consequence only: native fraction kf/(kf+ell+h).
        kf,ell=Q(1),Q(1)
        yield_at=lambda h: kf/(kf+ell+h)
        self.assertLess(yield_at(Q(2)),yield_at(Q(0)))

    def test_18_feasibility_not_realized_output(self):
        # Two routes, success 1 and 0, both with adequate pools.
        s=Q(1); quota=Q(4,5); productive_hazard=Q(1); disposal_hazard=Q(3)
        optimum=s; realized=s*productive_hazard/(productive_hazard+disposal_hazard)
        self.assertGreaterEqual(optimum,quota)
        self.assertLess(realized,quota)

    @unittest.skipIf(sp is None, 'SymPy not installed; exact rational checks still run')
    def test_19_symbolic_conservation_and_machine_identity(self):
        U,C,A,N,F=sp.symbols('U C A N F', nonnegative=True)
        p=sp.symbols('s q ku r f kf du dc da dn mu kon ka alpha')
        x=[U,C,A,N,F]; dx=rhs(x,p)
        expected=p[0]-p[10]*(U+C+A+N)-p[6]*U-p[7]*C-p[8]*A-p[9]*N
        self.assertEqual(sp.expand(sum(dx[:4])-expected),0)
        self.assertEqual(sp.expand(dx[1]+dx[4]-p[13]+p[10]*(C+F)),0)

    @unittest.skipIf(sp is None, 'SymPy not installed; exact rational checks still run')
    def test_20_symbolic_steady_occupancy_and_boundary_signs(self):
        b,r,f,d,mu=sp.symbols('b r f d mu', positive=True)
        C=b/(r+f+d+mu)
        self.assertEqual(sp.simplify(b-(r+f+d+mu)*C),0)
        U,C,A,N,F=sp.symbols('U C A N F', nonnegative=True)
        s,ku,r,f,kf,du,dc,da,dn,mu,kon,ka,alpha=sp.symbols(
            's ku r f kf du dc da dn mu kon ka alpha', nonnegative=True)
        q=sp.symbols('q', nonnegative=True)
        dx=rhs([U,C,A,N,F],[s,q,ku,r,f,kf,du,dc,da,dn,mu,kon,ka,alpha])
        for x,derivative in zip([U,C,A,F],[dx[0],dx[1],dx[2],dx[4]]):
            self.assertTrue(sp.expand(derivative.subs(x,0)).is_nonnegative)
        # N-face: replace 1-q with an explicitly nonnegative fraction.
        z=sp.symbols('z', nonnegative=True)
        self.assertTrue(sp.expand(dx[3].subs(N,0).subs(q,1-z)).is_nonnegative)


if __name__ == '__main__':
    print('Translation-associated proteostasis: mathematical checks', flush=True)
    print('All witnesses illustrative; no physiological headroom, fit or wet-lab validation.', flush=True)
    print('Python:',sys.version.split()[0], 'SymPy:', sp.__version__ if sp else 'absent', flush=True)
    suite=unittest.defaultTestLoader.loadTestsFromTestCase(Accounting)
    result=unittest.TextTestRunner(verbosity=2,stream=sys.stdout).run(suite)
    passed=result.testsRun-len(result.failures)-len(result.errors)-len(result.skipped)
    print(f'COUNTS: {result.testsRun} run; {passed} passed; {len(result.failures)} failed; '
          f'{len(result.errors)} errors; {len(result.skipped)} skipped')
    print('SCIENTIFIC STATUS: bookkeeping checked; empirical validation unavailable.')
    sys.exit(0 if result.wasSuccessful() else 1)
