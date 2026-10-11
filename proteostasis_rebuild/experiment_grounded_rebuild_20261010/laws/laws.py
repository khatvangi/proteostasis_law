"""
proteostasis laws: reproducible numerics and symbolic checks for THEOREMS.md.

every parameter in this file is a dimensionless placeholder chosen to exhibit a
sign, a bound or a counterexample. none is an E. coli estimate and nothing here
is fitted to data.

sections (numbering follows FORMAL_MODEL.md)
  L3  time-to-depletion bounds for dN/dt = P - lam N
  L4  route-qualified linear feasibility (LP), Farkas certificate, margin
  L5  realized competing-fates kinetics (single client, finite pools)
  S1  capacity-monotone versus realized-sign separation
  S2  induction race: steady-state feasibility does not give maintenance
  L6  sigma32 co-regulation of several machines
  L8  codon-conditioned perturbations
  L9  optional energy row

run `python laws.py` to write results.json; tests live in test_laws.py.
"""
import itertools
import json
import os

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, linprog

HERE = os.path.dirname(os.path.abspath(__file__))


# ----------------------------------------------------------------------------
# L3. time-to-depletion bounds
# ----------------------------------------------------------------------------

def findNUpper(t, n0, pHi, lamLo):
    """comparison solution of y' = pHi - lamLo*y, y(0)=n0. upper bound on N(t)
    whenever P(t) <= pHi, lam(t) >= lamLo > 0 and N >= 0."""
    return n0 * np.exp(-lamLo * t) + pHi / lamLo * (1.0 - np.exp(-lamLo * t))


def findNLower(t, n0, pLo, lamHi):
    """comparison solution of y' = pLo - lamHi*y. lower bound on N(t)
    whenever P(t) >= pLo >= 0 and lam(t) <= lamHi."""
    return n0 * np.exp(-lamHi * t) + pLo / lamHi * (1.0 - np.exp(-lamHi * t))


def findCrossingTime(n0, p, lam, req):
    """time at which y' = p - lam*y, y(0)=n0 reaches req.
    returns 0 if n0 <= req, inf if the asymptote p/lam >= req (never crosses)."""
    if n0 <= req:
        return 0.0
    ss = p / lam
    if ss >= req:
        return np.inf
    return np.log((n0 - ss) / (req - ss)) / lam


def findFailureWindow(n0, pLo, pHi, lamLo, lamHi, reqLo, reqHi):
    """[t_lo, t_up] bracketing the first time N(t) < Nreq(t), given
    P in [pLo, pHi], lam in [lamLo, lamHi], Nreq in [reqLo, reqHi], n0 > reqHi.
    t_up = inf means failure is not guaranteed; t_lo = inf means it is excluded."""
    tUp = findCrossingTime(n0, pHi, lamLo, reqLo)
    tLo = findCrossingTime(n0, pLo, lamHi, reqHi)
    return tLo, tUp


def checkCrossingSymbolic():
    """sympy: the comparison solution and the closed-form crossing time."""
    t, n0, p, lam, req = sp.symbols('t n0 p lam req', positive=True)
    y = sp.Function('y')
    sol = sp.dsolve(sp.Eq(y(t).diff(t), p - lam * y(t)), y(t), ics={y(0): n0}).rhs
    closed = n0 * sp.exp(-lam * t) + p / lam * (1 - sp.exp(-lam * t))
    tc = sp.log((n0 - p / lam) / (req - p / lam)) / lam
    return {'ode_solution_matches': sp.simplify(sol - closed) == 0,
            'crossing_hits_req': sp.simplify(closed.subs(t, tc) - req) == 0}


def simulateScalarQuota(n0, pFun, lamFun, reqFun, tEnd, nPts=4001):
    """integrate dN/dt = P(t) - lam(t) N and return the first time N < Nreq."""
    sol = solve_ivp(lambda t, y: [pFun(t) - lamFun(t) * y[0]], (0, tEnd), [n0],
                    dense_output=True, rtol=1e-10, atol=1e-12, max_step=tEnd / 2000)
    ts = np.linspace(0, tEnd, nPts)
    ns = sol.sol(ts)[0]
    gap = ns - np.array([reqFun(x) for x in ts])
    bad = np.nonzero(gap < 0)[0]
    return (ts[bad[0]] if len(bad) else np.inf), ts, ns


def findMovingQuotaRun(P, k=0.1, n0=1.0, n1=2.0, mu0=0.5, mu1=1.0, r=1.0, tEnd=10.0):
    """moving quota face. Nreq(mu) = n0 + n1 mu, mu relaxes mu0 -> mu1 at rate r,
    dN/dt = P - (mu + k) N, N(0) = Nreq(mu0). returns min_t [N - Nreq(mu(t))]."""
    def rhs(t, y):
        N, mu = y
        return [P - (mu + k) * N, r * (mu1 - mu)]
    sol = solve_ivp(rhs, (0, tEnd), [n0 + n1 * mu0, mu0], rtol=1e-11, atol=1e-13,
                    dense_output=True)
    ts = np.linspace(0, tEnd, 4001)
    N, mu = sol.sol(ts)
    return float(np.min(N - (n0 + n1 * mu)))


def findMovingQuotaP(k=0.1, n0=1.0, n1=2.0, mu0=0.5, mu1=1.0, r=1.0):
    """(static) P that meets the quota with zero slack at mu0, and the
    (moving) P that meets the face condition (A7.1) at every mu in [mu0, mu1]:
    P >= (mu + k) Nreq(mu) + Nreq'(mu) * r (mu1 - mu)."""
    static = (mu0 + k) * (n0 + n1 * mu0)
    mus = np.linspace(mu0, mu1, 2001)
    moving = float(np.max((mus + k) * (n0 + n1 * mus) + n1 * r * (mu1 - mus)))
    return static, moving


# ----------------------------------------------------------------------------
# L4. route-qualified linear feasibility
# ----------------------------------------------------------------------------

class RouteLP:
    """linear route feasibility.

    routes: list of dicts
        client: client name
        a:      {machine: expected machine-time per route entry}  [time]
        rho:    probability the route ends native                   [-]
        b:      probability the route returns the chain to U        [-] (default 0)
        cap:    optional ceiling on route flux                      [conc/time]
        energy: optional cost per unit route flux                   [energy/conc]
    machines: {machine: available pool M_m}                          [conc]
    supply:   {client: e_i, non-native entry flux}                  [conc/time]
    demand:   {client: d_i, native deficit to be met by routes}     [conc/time]

    rows:  sum_r a_mr x_r <= M_m ;  sum_{r in i} (1-b_r) x_r <= e_i ;
           sum_{r in i} rho_r x_r >= d_i ;  x_r <= cap_r ;  x >= 0.
    """

    def __init__(self, routes, machines, supply, demand, energyBudget=None):
        self.routes = routes
        self.machines = dict(machines)
        self.supply = dict(supply)
        self.demand = dict(demand)
        self.energyBudget = energyBudget
        self.mNames = list(self.machines)
        self.cNames = list(self.supply)

    def buildRows(self, theta=None):
        """inequality rows G x <= h (x = route fluxes, optionally plus theta).
        returns G, h and row labels. if theta is given as 'var', demand rows
        become theta*d_i - sum rho x <= 0 and a last column holds theta."""
        R = len(self.routes)
        ncol = R + (1 if theta == 'var' else 0)
        G, h, lab = [], [], []
        for m in self.mNames:
            G.append([r['a'].get(m, 0.0) for r in self.routes] + [0.0] * (ncol - R))
            h.append(self.machines[m]); lab.append(('capacity', m))
        for c in self.cNames:
            G.append([(1.0 - r.get('b', 0.0)) if r['client'] == c else 0.0
                      for r in self.routes] + [0.0] * (ncol - R))
            h.append(self.supply[c]); lab.append(('supply', c))
        for c in self.cNames:
            row = [-r['rho'] if r['client'] == c else 0.0 for r in self.routes]
            if theta == 'var':
                G.append(row + [self.demand.get(c, 0.0)]); h.append(0.0)
            else:
                G.append(row); h.append(-self.demand.get(c, 0.0))
            lab.append(('demand', c))
        for k, r in enumerate(self.routes):
            if r.get('cap') is not None:
                row = [0.0] * ncol; row[k] = 1.0
                G.append(row); h.append(r['cap']); lab.append(('ceiling', k))
        if self.energyBudget is not None:
            G.append([r.get('energy', 0.0) for r in self.routes] + [0.0] * (ncol - R))
            h.append(self.energyBudget); lab.append(('energy', None))
        return np.array(G, float), np.array(h, float), lab

    def findFeasible(self):
        """phase-one LP. returns (feasible, x)."""
        G, h, _ = self.buildRows()
        res = linprog(np.zeros(G.shape[1]), A_ub=G, b_ub=h, bounds=(0, None),
                      method='highs', options={'presolve': False})
        if res.status == 0:
            x = res.x
            return bool(np.all(G @ x <= h + 1e-9)), x
        return False, None

    def findFarkasCertificate(self):
        """search y >= 0 with G^T y >= 0 and h^T y < 0 (Farkas alternative for
        G x <= h, x >= 0). y is boxed in [0,1] so the LP is bounded.
        returns (y, h^T y) or (None, value)."""
        G, h, lab = self.buildRows()
        res = linprog(h, A_ub=-G.T, b_ub=np.zeros(G.shape[1]),
                      bounds=[(0, 1)] * len(h), method='highs',
                      options={'presolve': False})
        if res.status == 0 and res.fun < -1e-9:
            return dict(zip(lab, res.x)), res.fun
        return None, (res.fun if res.status == 0 else None)

    def checkCertificate(self, y):
        """verify a certificate given in physical form (dict keyed by row label).
        returns (route reduced costs, weighted demand minus weighted resources)."""
        G, h, lab = self.buildRows()
        yv = np.array([y.get(l, 0.0) for l in lab])
        if np.any(yv < 0):
            return None, None
        return G.T @ yv, -(h @ yv)

    def findMargin(self):
        """max-load margin: largest theta with theta*d feasible. returns
        (theta, x, row duals keyed by label). duals are d theta / d h_row."""
        G, h, lab = self.buildRows(theta='var')
        R = len(self.routes)
        cost = np.zeros(R + 1); cost[-1] = -1.0
        res = linprog(cost, A_ub=G, b_ub=h, bounds=[(0, None)] * R + [(0, 1e6)],
                      method='highs', options={'presolve': False})
        if res.status != 0:
            return None, None, None
        duals = dict(zip(lab, -res.ineqlin.marginals))
        return res.x[-1], res.x[:R], duals

    def findPooled(self):
        """global-pooling relaxation: one machine whose pool is sum M_m and whose
        coefficient for a route is sum_m a_mr. returns a new RouteLP."""
        pooled = [dict(r, a={'pool': sum(r['a'].values())}) for r in self.routes]
        return RouteLP(pooled, {'pool': sum(self.machines.values())},
                       self.supply, self.demand, self.energyBudget)


def checkHallUnion(lp):
    """the casual 'Hall' test: for every client subset S, sum of demands <=
    sum over machines in the union neighbourhood of M_m / min a_mr. treats the
    machines a client touches as substitutable. returns True if it passes."""
    cl = lp.cNames
    for k in range(1, len(cl) + 1):
        for S in itertools.combinations(cl, k):
            nbr = {}
            for r in lp.routes:
                if r['client'] in S:
                    for m, a in r['a'].items():
                        if a > 0:
                            nbr[m] = min(nbr.get(m, np.inf), a)
            capacity = sum(lp.machines[m] / a for m, a in nbr.items())
            if sum(lp.demand[c] for c in S) > capacity + 1e-12:
                return False
    return True


def checkBottleneckMin(lp):
    """the opposite casual rule: each client limited by the smallest
    M_m / a_mr among the machines it touches (treats every route as serial)."""
    for c in lp.cNames:
        rates = [lp.machines[m] / a for r in lp.routes if r['client'] == c
                 for m, a in r['a'].items() if a > 0]
        if rates and lp.demand[c] > min(rates) + 1e-12:
            return False
    return True


def checkGaleSubsets(demand, M, tau, qual):
    """Gale/Hall cut condition for single-machine routes with rho = 1 and a
    per-machine service time: for all client subsets S,
    sum_{i in S} d_i <= sum_{m in N(S)} M_m / tau_m."""
    cl = list(demand)
    for k in range(1, len(cl) + 1):
        for S in itertools.combinations(cl, k):
            nbr = set().union(*[qual[c] for c in S])
            if sum(demand[c] for c in S) > sum(M[m] / tau[m] for m in nbr) + 1e-12:
                return False
    return True


def findSerialToyLP(demand=1.5):
    """one client, one serial route through A then B (a_A = a_B = 1), M_A = M_B = 1."""
    routes = [{'client': 'i', 'a': {'A': 1.0, 'B': 1.0}, 'rho': 1.0}]
    return RouteLP(routes, {'A': 1.0, 'B': 1.0}, {'i': 10.0}, {'i': demand})


def findParallelToyLP(demand=1.5):
    """one client, two parallel single-machine routes, M_A = M_B = 1."""
    routes = [{'client': 'i', 'a': {'A': 1.0}, 'rho': 1.0},
              {'client': 'i', 'a': {'B': 1.0}, 'rho': 1.0}]
    return RouteLP(routes, {'A': 1.0, 'B': 1.0}, {'i': 10.0}, {'i': demand})


def findPoolingToyLP():
    """client g can use only machine G (obligate); client k only K.
    G is short, K has plenty: pooled feasible, qualified infeasible."""
    routes = [{'client': 'g', 'a': {'G': 2.0}, 'rho': 1.0},
              {'client': 'k', 'a': {'K': 1.0}, 'rho': 1.0}]
    return RouteLP(routes, {'G': 0.5, 'K': 10.0}, {'g': 5.0, 'k': 5.0},
                   {'g': 1.0, 'k': 1.0})


# ----------------------------------------------------------------------------
# L5. realized competing-fates kinetics, single client, finite pools
# ----------------------------------------------------------------------------

def findYieldFixedPools(kf, kl, theta, pN, pU):
    """per-entry native yield from U with fixed binding hazards theta_m.
    Y = (kf + sum theta pN) / (kf + kl + sum theta (1-pU))."""
    theta, pN, pU = map(np.asarray, (theta, pN, pU))
    return (kf + np.sum(theta * pN)) / (kf + kl + np.sum(theta * (1 - pU)))


def checkSignLawSymbolic():
    """sympy: dY/dtheta_1 = (1-pU1)/Q * (pN1/(1-pU1) - Y) for two machines."""
    kf, kl, t1, t2, n1, n2, u1, u2 = sp.symbols('kf kl th1 th2 pN1 pN2 pU1 pU2',
                                                positive=True)
    Q = kf + kl + t1 * (1 - u1) + t2 * (1 - u2)
    Y = (kf + t1 * n1 + t2 * n2) / Q
    claim = (1 - u1) / Q * (n1 / (1 - u1) - Y)
    return sp.simplify(sp.diff(Y, t1) - claim) == 0


def findSteadySingleClient(e, kf, kl, kon, M, tau, pN, pU):
    """steady state of one client with finite pools (no sequestration).
    bound C_m = kon_m M_m U tau_m / (1 + kon_m tau_m U)  (C = J tau, J = kon (M-C) U)
    entry balance e = (kf+kl) U + sum (1-pU_m) C_m / tau_m is increasing in U,
    so the root is unique. returns dict with U, C, free pools, Y, occupancy."""
    kon, M, tau, pN, pU = map(np.asarray, (kon, M, tau, pN, pU))

    def findBound(U):
        return kon * M * U * tau / (1 + kon * tau * U)

    def findResidual(U):
        return (kf + kl) * U + np.sum((1 - pU) * findBound(U) / tau) - e

    hi = 1.0
    while findResidual(hi) < 0:
        hi *= 2
    U = brentq(findResidual, 0.0, hi, xtol=1e-15, rtol=1e-14)
    C = findBound(U)
    Y = (kf * U + np.sum(pN * C / tau)) / e
    occ = np.divide(C, M, out=np.zeros_like(C), where=M > 0)
    return {'U': U, 'C': C, 'free': M - C, 'Y': Y, 'occupancy': occ,
            'J': C / tau, 'theta': kon * (M - C)}


def findSingleClientRouteLP(e, demand, kf, kon, M, tau, pN, pU, uCap, names):
    """route LP built from the same single-client instance: one spontaneous
    route (rho 1, ceiling kf*uCap), one one-step route per machine
    (a = tau, rho = pN, b = pU, ceiling kon*M*uCap)."""
    routes = [{'client': 'i', 'a': {}, 'rho': 1.0, 'cap': kf * uCap}]
    for k, nm in enumerate(names):
        routes.append({'client': 'i', 'a': {nm: tau[k]}, 'rho': pN[k], 'b': pU[k],
                       'cap': kon[k] * M[k] * uCap})
    return RouteLP(routes, dict(zip(names, M)), {'i': e}, {'i': demand})


# S1 instance: a productive machine G and a disposal-dominated machine H
S1 = dict(e=1.0, kf=0.3, kl=0.2, kon=[1.0, 6.0], tau=[2.0, 1.0],
          pN=[0.9, 0.1], pU=[0.05, 0.0], names=['G', 'H'])


def findS1Comparison(MH_list=(0.0, 0.5, 2.0, 5.0), MG=5.0, demand=0.6, uFactor=2.0):
    """realized yield versus LP margin as the H pool grows."""
    p = S1
    rows = []
    for MH in MH_list:
        M = np.array([MG, MH])
        st = findSteadySingleClient(p['e'], p['kf'], p['kl'], p['kon'], M,
                                    p['tau'], p['pN'], p['pU'])
        # ceiling on the non-native inventory: a fixed tolerated level,
        # identical for every MH, so the LP compares capacity only
        uCap = uFactor
        lp = findSingleClientRouteLP(p['e'], demand, p['kf'], p['kon'], M,
                                     p['tau'], p['pN'], p['pU'], uCap, p['names'])
        theta, _, _ = lp.findMargin()
        rows.append({'MH': MH, 'Y_realized': st['Y'], 'margin_LP': theta,
                     'occupancy_G': st['occupancy'][0],
                     'occupancy_H': st['occupancy'][1] if MH > 0 else 0.0,
                     'committed_H': p['pN'][1] / (1 - p['pU'][1])})
    return rows


def findTwoMachineSignExceptions(seed=13, nDraws=400):
    """with two finite pools, how often does sign(dY/dM_1) differ from
    sign(committed_1 - Y)? the single-pool transfer is proved; this is not."""
    rng = np.random.default_rng(seed)
    bad = tot = 0
    for _ in range(nDraws):
        kf, kl = rng.uniform(0, 1), rng.uniform(0.01, 1)
        kon = rng.uniform(0.1, 10, 2); M = rng.uniform(0.1, 5, 2)
        tau = rng.uniform(0.1, 3, 2); pN = rng.uniform(0, 0.95, 2)
        pU = rng.uniform(0, 0.99, 2) * (1 - pN); e = rng.uniform(0.1, 3)

        def findY(m1):
            return findSteadySingleClient(e, kf, kl, kon, [m1, M[1]], tau, pN, pU)['Y']
        Y = findY(M[0]); committed = pN[0] / (1 - pU[0])
        if abs(committed - Y) < 1e-4:
            continue
        dY = findY(M[0] * (1 + 1e-6)) - findY(M[0] * (1 - 1e-6))
        tot += 1
        bad += int(np.sign(dY) != np.sign(committed - Y))
    return bad, tot


def findLowOccupancyFailure(eps):
    """family showing quota failure with every machine below occupancy eps.
    a disposal machine with huge capture and committed success ~0."""
    kf, kl = 1.0, 0.0
    big = 1.0 / eps ** 2
    st = findSteadySingleClient(1.0, kf, kl, [1.0, big], [1.0 / eps, 1.0 / eps],
                                [1.0, eps ** 3], [1.0, 0.0], [0.0, 0.0])
    return st


# ----------------------------------------------------------------------------
# S2. induction race (closed-loop mechanistic model with sigma32-type control)
# ----------------------------------------------------------------------------

# chosen from a 400-draw scan (seed 1) as a case where the bound itself dips
# below the closed-loop post-shock steady state, then rounded. placeholders.
RACE = dict(s=2.5, sigma=0.2, kon=100.0, kU=0.1, kN=1.5, kA=0.0, ka=0.1, kpU=0.05,
            ku0=0.02, ku1=1.0, kdN=0.05, kdA=0.1, mu=0.035, b0=0.01, b1=0.9,
            K=0.1, n=2.0)
RACE_NREQ = 12.0


def findRaceRhs(t, y, p, ku):
    """states U, C, N, A, D(cum degraded), L(cum diluted), M(total machine).
    machine synthesis b0 + b1 h(F), F = M - C free machine, h = 1/(1+(F/K)^n).
    machine dilution only (gamma = mu), so a bound unit leaves with its client."""
    U, C, N, A, D, L, M = y
    F = max(M - C, 0.0)
    # 'forceU' replaces the sigma32-type law by a fixed regulator activity,
    # e.g. forceU = 1: maximal induction from t = 0 with zero delay
    h = p['forceU'] if 'forceU' in p else 1.0 / (1.0 + (F / p['K']) ** p['n'])
    bind = p['kon'] * F * U
    kout = p['kN'] + p['kU'] + p['kA']
    dU = ((1 - p['sigma']) * p['s'] + ku * N + p['kU'] * C - bind
          - (p['ka'] + p['kpU'] + p['mu']) * U)
    dC = bind - (kout + p['mu']) * C
    dN = p['sigma'] * p['s'] + p['kN'] * C - (ku + p['kdN'] + p['mu']) * N
    dA = p['ka'] * U + p['kA'] * C - (p['kdA'] + p['mu']) * A
    dD = p['kpU'] * U + p['kdN'] * N + p['kdA'] * A
    dL = p['mu'] * (U + C + N + A)
    dM = p['b0'] + p['b1'] * h - p['mu'] * M
    return [dU, dC, dN, dA, dD, dL, dM]


def findRacePreShock(p=RACE, tEnd=4000.0):
    """pre-shock steady state (ku = ku0) by long integration."""
    y0 = [0, 0, 0, 0, 0, 0, p['b0'] / p['mu']]
    sol = solve_ivp(findRaceRhs, (0, tEnd), y0, args=(p, p['ku0']),
                    rtol=1e-11, atol=1e-13, method='LSODA')
    y = sol.y[:, -1].copy()
    y[4] = y[5] = 0.0
    return y


def findRaceEquilibriumEigen(p=RACE, ku=None):
    """eigenvalues of the Jacobian (U, C, N, A, M) at the closed-loop
    equilibrium reached after the shock (central finite differences)."""
    ku = p['ku1'] if ku is None else ku
    sol = solve_ivp(findRaceRhs, (0, 6000.0), findRacePreShock(p), args=(p, ku),
                    rtol=1e-12, atol=1e-14, method='LSODA')
    y = sol.y[:, -1]
    idx = [0, 1, 2, 3, 6]
    J = np.zeros((5, 5))
    for c, k in enumerate(idx):
        h = 1e-6 * max(1.0, abs(y[k]))
        yp, ym = y.copy(), y.copy(); yp[k] += h; ym[k] -= h
        fp = np.array(findRaceRhs(0, yp, p, ku))[idx]
        fm = np.array(findRaceRhs(0, ym, p, ku))[idx]
        J[:, c] = (fp - fm) / (2 * h)
    resid = float(np.max(np.abs(np.array(findRaceRhs(0, y, p, ku))[idx])))
    return np.linalg.eigvals(J), resid


def findRaceUpperBound(t, n0, p0, kappa, M0, Mc, gam, lam):
    """N(t) <= n0 e^{-lam t} + (p0 + kappa Mc)/lam (1 - e^{-lam t})
               + kappa (M0 - Mc)(e^{-gam t} - e^{-lam t})/(lam - gam)."""
    return (n0 * np.exp(-lam * t) + (p0 + kappa * Mc) / lam * (1 - np.exp(-lam * t))
            + kappa * (M0 - Mc) * (np.exp(-gam * t) - np.exp(-lam * t)) / (lam - gam))


def checkRaceBoundSymbolic():
    """sympy: the race bound solves y' = p0 + kappa*Mbar(t) - lam*y."""
    t, n0, p0, ka, M0, Mc, g, lam = sp.symbols('t n0 p0 kappa M0 Mc gamma lam',
                                               positive=True)
    Mbar = Mc + (M0 - Mc) * sp.exp(-g * t)
    yb = (n0 * sp.exp(-lam * t) + (p0 + ka * Mc) / lam * (1 - sp.exp(-lam * t))
          + ka * (M0 - Mc) * (sp.exp(-g * t) - sp.exp(-lam * t)) / (lam - g))
    ode = sp.simplify(yb.diff(t) - (p0 + ka * Mbar - lam * yb))
    return {'solves_comparison_ode': ode == 0,
            'initial_value': sp.simplify(yb.subs(t, 0) - n0) == 0}


def findRaceRun(p=RACE, tEnd=200.0, nPts=20001):
    """post-shock trajectory (ku0 -> ku1 at t = 0) from the pre-shock state.
    the pre-shock state always uses the sigma32-type law. returns trajectories,
    the race bound, quota bookkeeping and conservation."""
    pre = {k: v for k, v in p.items() if k != 'forceU'}
    y0 = findRacePreShock(pre)
    sol = solve_ivp(findRaceRhs, (0, tEnd), y0, args=(p, p['ku1']),
                    rtol=1e-11, atol=1e-13, method='LSODA', dense_output=True)
    ts = np.linspace(0, tEnd, nPts)
    Y = sol.sol(ts)
    U, C, N, A, D, L, M = Y
    lam = p['mu'] + p['ku1'] + p['kdN']
    Mc = (p['b0'] + p['b1']) / p['mu']
    p0 = p['sigma'] * p['s']
    bound = findRaceUpperBound(ts, y0[2], p0, p['kN'], y0[6], Mc, p['mu'], lam)
    # conservation control: U+C+N+A+D+L must equal its start plus s*t
    mass = U + C + N + A + D + L
    massErr = np.max(np.abs(mass - (mass[0] + p['s'] * ts)))
    # long-run post-shock state
    sol2 = solve_ivp(findRaceRhs, (0, 4000.0), y0, args=(p, p['ku1']),
                     rtol=1e-11, atol=1e-13, method='LSODA')
    Nss = sol2.y[2, -1]
    return {'t': ts, 'N': N, 'M': M, 'C': C, 'U': U, 'bound': bound,
            'N0': y0[2], 'M0': y0[6], 'Mc': Mc, 'lam': lam, 'p0': p0,
            'Nss_post': Nss, 'massErr': massErr,
            'minState': float(np.min(Y[:4])), 'minFree': float(np.min(M - C)),
            'steady_bound_full_induction': (p0 + p['kN'] * Mc) / lam}


def findRaceSummary(p=RACE, Nreq=RACE_NREQ):
    """race counterexample in one dict. the quota lies between the transient
    minimum of the bound and the closed-loop post-shock steady state."""
    r = findRaceRun(p)
    below = np.nonzero(r['N'] < Nreq)[0]
    boundBelow = np.nonzero(r['bound'] < Nreq)[0]
    return {'Nreq': Nreq, 'N0': r['N0'], 'Nss_post': r['Nss_post'],
            'steady_bound_full_induction': r['steady_bound_full_induction'],
            'min_N': float(np.min(r['N'])), 'min_bound': float(np.min(r['bound'])),
            'first_violation_t': float(r['t'][below[0]]) if len(below) else None,
            'bound_guarantees_violation_by_t': float(r['t'][boundBelow[0]]) if len(boundBelow) else None,
            'max_N_minus_bound': float(np.max(r['N'] - r['bound'])),
            'mass_error': float(r['massErr']), 'min_state': r['minState'],
            'min_free_machine': r['minFree'],
            'lam': r['lam'], 'mu': p['mu'], 'M0': r['M0'], 'Mc': r['Mc'],
            't_min_bound': float(r['t'][int(np.argmin(r['bound']))]),
            'post_shock_occupancy': float(r['C'][-1] / r['M'][-1])}


def findRandomControlCeiling(seed=0, nRuns=50):
    """machine ceiling M(t) <= Mbar(t) under arbitrary bounded regulation
    u(t) in [0,1] (random piecewise-constant, with arbitrary delays folded in)."""
    rng = np.random.default_rng([2026, seed])
    worst = -np.inf
    for _ in range(nRuns):
        b0, b1, gam = rng.uniform(0.01, 1), rng.uniform(0, 2), rng.uniform(0.05, 1)
        M0 = rng.uniform(0, 3 * (b0 + b1) / gam)
        knots = np.sort(rng.uniform(0, 30, 12)); vals = rng.uniform(0, 1, 13)
        u = lambda t: vals[np.searchsorted(knots, t)]
        sol = solve_ivp(lambda t, y: [b0 + b1 * u(t) - gam * y[0]], (0, 30), [M0],
                        dense_output=True, rtol=1e-10, atol=1e-12, max_step=0.05)
        ts = np.linspace(0, 30, 3001)
        Mc = (b0 + b1) / gam
        Mbar = Mc + (M0 - Mc) * np.exp(-gam * ts)
        worst = max(worst, float(np.max(sol.sol(ts)[0] - Mbar)))
    return worst


# ----------------------------------------------------------------------------
# L6. sigma32 co-regulation: one input, several machines
# ----------------------------------------------------------------------------

def findCoRegulatedYield(u, base=(0.1, 0.0), gain=(2.0, 2.0), kf=0.05, **kw):
    """single client, productive machine G and disposal machine H both driven by
    one regulator activity u in [0,1]: M_m = base_m + gain_m * u.
    returns the realized yield (finite pools)."""
    p = dict(S1); p.update(kw); p['kf'] = kf
    M = np.array([base[0] + gain[0] * u, base[1] + gain[1] * u])
    return findSteadySingleClient(p['e'], p['kf'], p['kl'], p['kon'], M,
                                  p['tau'], p['pN'], p['pU'])['Y']


# ----------------------------------------------------------------------------
# L8. codon-conditioned perturbations
# ----------------------------------------------------------------------------

def checkCotranslationalSymbolic():
    """Y_tot = sigma + (1-sigma) Y_birth  =>  dY_tot/dsigma = 1 - Y_birth."""
    sg, yb = sp.symbols('sigma Y_birth', nonnegative=True)
    Yt = sg + (1 - sg) * yb
    return sp.simplify(sp.diff(Yt, sg) - (1 - yb)) == 0


TWO = dict(s=(1.0, 1.0), kf=(0.05, 0.6), kl=(0.3, 0.2),
           kon=((1.0, 8.0), (1.0, 1.0)), tau=(2.0, 1.0),
           pN=((0.9, 0.0), (0.9, 0.05)), pU=((0.05, 0.0), (0.05, 0.0)),
           M=(3.0, 2.0))


def findTwoClientSteady(sigma1, p=TWO, M=None):
    """two clients sharing machines G (index 0) and H (index 1). per client
    entry e_i = (1-sigma_i) s_i (sigma_2 = 0), per machine free pool solved
    self-consistently: C_mi = kon_mi F_m U_i tau_m, F_m = M_m - sum_i C_mi,
    U_i from e_i = (kf_i+kl_i) U_i + sum_m (1-pU_mi) C_mi / tau_m.
    returns yields Y_i (native flux per synthesized chain) and free pools."""
    M = np.array(p['M'] if M is None else M, float)
    tau = np.array(p['tau'])
    kon = np.array(p['kon']); pN = np.array(p['pN']); pU = np.array(p['pU'])
    e = np.array([(1 - sigma1) * p['s'][0], p['s'][1]])
    kf = np.array(p['kf']); kl = np.array(p['kl'])

    def findU(F):
        # given free pools, each U_i is linear
        Q = kf + kl + np.sum(kon * F[None, :] * (1 - pU), axis=1)
        return e / Q

    def findFree(F):
        U = findU(F)
        C = kon * F[None, :] * U[:, None] * tau[None, :]
        return M - C.sum(axis=0)

    F = M.copy()
    for _ in range(20000):
        Fn = findFree(F)
        Fn = 0.5 * F + 0.5 * np.clip(Fn, 1e-15, M)
        if np.max(np.abs(Fn - F)) < 1e-15:
            F = Fn
            break
        F = Fn
    U = findU(F)
    C = kon * F[None, :] * U[:, None] * tau[None, :]
    nat = kf * U + np.sum(pN * C / tau[None, :], axis=1)
    sigmaV = np.array([sigma1, 0.0])
    Y = (sigmaV * np.array(p['s']) + nat) / np.array(p['s'])
    resid = float(np.max(np.abs(F - (M - C.sum(axis=0)))))
    return {'Y': Y, 'F': F, 'U': U, 'C': C, 'resid': resid}


def findTwoClientOde(sigma1, p=TWO, tEnd=2000.0):
    """same two-client system integrated as an ODE (independent of the fixed
    point) with conservation counters. returns end state yields and mass error."""
    M = np.array(p['M'], float); tau = np.array(p['tau'])
    kon = np.array(p['kon']); pN = np.array(p['pN']); pU = np.array(p['pU'])
    kf = np.array(p['kf']); kl = np.array(p['kl']); s = np.array(p['s'])
    sig = np.array([sigma1, 0.0])
    koff = 1.0 / tau   # total exit rate from a bound state (pN + pU + pD = 1)

    def rhs(t, y):
        U = y[0:2]; C = y[2:6].reshape(2, 2); Ncum = y[6:8]; Lcum = y[8:10]
        F = M - C.sum(axis=0)
        bind = kon * F[None, :] * U[:, None]
        out = C * koff[None, :]
        dU = (1 - sig) * s + np.sum(pU * out, axis=1) - bind.sum(axis=1) - (kf + kl) * U
        dC = bind - out
        dN = sig * s + kf * U + np.sum(pN * out, axis=1)
        dL = kl * U + np.sum((1 - pN - pU) * out, axis=1)
        return np.concatenate([dU, dC.ravel(), dN, dL])

    sol = solve_ivp(rhs, (0, tEnd), np.zeros(10), rtol=1e-11, atol=1e-13,
                    method='LSODA', dense_output=True)
    y = sol.y[:, -1]
    t1, t2 = tEnd - 200.0, tEnd
    y1 = sol.sol(t1)
    Y = (y[6:8] - y1[6:8]) / (t2 - t1) / s
    mass = sol.y[0:2] + sol.y[2:6].reshape(2, 2, -1).sum(axis=1) + sol.y[6:8] + sol.y[8:10]
    massErr = float(np.max(np.abs(mass - s[:, None] * sol.t[None, :])))
    C = y[2:6].reshape(2, 2)
    return {'Y': Y, 'massErr': massErr, 'maxOccupancy': float(np.max(C.sum(axis=0) / M)),
            'minState': float(np.min(sol.y[0:6]))}


def findMistranslationLoadBound(nCodons, eps):
    """fraction of chains carrying at least one error from the listed codon
    classes, and its first-order upper bound sum n_c eps_c."""
    nCodons, eps = np.asarray(nCodons, float), np.asarray(eps, float)
    exact = 1.0 - np.prod((1.0 - eps) ** nCodons)
    return exact, float(np.sum(nCodons * eps))


def findMarginSensitivity(delta=1e-6):
    """margin LP with one binding and one slack machine. compares finite
    differences in a route coefficient with -w_m * x_r * theta-scaled dual."""
    routes = [{'client': 'i', 'a': {'G': 2.0, 'K': 0.5}, 'rho': 0.9},
              {'client': 'j', 'a': {'K': 1.0}, 'rho': 0.8}]
    machines = {'G': 1.0, 'K': 20.0}
    supply = {'i': 5.0, 'j': 5.0}
    demand = {'i': 0.3, 'j': 0.5}
    lp = RouteLP(routes, machines, supply, demand)
    th, x, du = lp.findMargin()
    out = {'theta': th, 'x': x.tolist(), 'w_G': du[('capacity', 'G')],
           'w_K': du[('capacity', 'K')]}
    for m in ('G', 'K'):
        r2 = [dict(r, a=dict(r['a'])) for r in routes]
        r2[0]['a'][m] += delta
        th2, _, _ = RouteLP(r2, machines, supply, demand).findMargin()
        out['fd_' + m] = (th2 - th) / delta
        out['pred_' + m] = -du[('capacity', m)] * x[0]
    return out


# ----------------------------------------------------------------------------
# L9. optional energy row
# ----------------------------------------------------------------------------

def findEnergyFront(budgets):
    """margin as a function of an energy budget. a cheap machine-free route
    that is mostly lossy (its losses cost resynthesis) versus an expensive
    productive machine route. without the energy row the problem has one
    objective (margin); the front exists only once the row is added."""
    routes = [{'client': 'i', 'a': {'G': 1.0}, 'rho': 0.95, 'energy': 7.0},
              {'client': 'i', 'a': {}, 'rho': 0.3, 'cap': 1.0, 'energy': 1.0}]
    out = []
    for E in budgets:
        lp = RouteLP(routes, {'G': 0.6}, {'i': 3.0}, {'i': 0.5}, energyBudget=E)
        th, _, _ = lp.findMargin()
        out.append(th)
    lp0 = RouteLP(routes, {'G': 0.6}, {'i': 3.0}, {'i': 0.5})
    return out, lp0.findMargin()[0]


# ----------------------------------------------------------------------------
# results
# ----------------------------------------------------------------------------

def findAllResults():
    res = {}
    res['L3_symbolic'] = checkCrossingSymbolic()
    tLo, tUp = findFailureWindow(1.0, 0.1, 0.25, 0.5, 0.8, 0.6, 0.6)
    res['L3_window_example'] = {'t_lo': tLo, 't_up': tUp}
    st, mv = findMovingQuotaP()
    res['L1_moving_quota'] = {'P_static': st, 'P_moving': mv,
                              'min_gap_static': findMovingQuotaRun(st),
                              'min_gap_moving': findMovingQuotaRun(mv)}
    res['L4_serial'] = {'lp_feasible': findSerialToyLP().findFeasible()[0],
                        'hall_union_passes': checkHallUnion(findSerialToyLP())}
    res['L4_parallel'] = {'lp_feasible': findParallelToyLP().findFeasible()[0],
                          'bottleneck_min_passes': checkBottleneckMin(findParallelToyLP())}
    pl = findPoolingToyLP()
    res['L4_pooling'] = {'qualified_feasible': pl.findFeasible()[0],
                         'pooled_feasible': pl.findPooled().findFeasible()[0]}
    cert, val = findSerialToyLP().findFarkasCertificate()
    res['L4_serial_certificate'] = {str(k): v for k, v in cert.items()} if cert else None
    res['S1_sign_symbolic'] = checkSignLawSymbolic()
    res['S1_comparison'] = findS1Comparison()
    res['S1_low_occupancy'] = {str(e): {'Y': float(st['Y']),
                                        'max_occupancy': float(np.max(st['occupancy']))}
                               for e, st in ((e, findLowOccupancyFailure(e))
                                             for e in (0.1, 0.03, 0.01))}
    bad, tot = findTwoMachineSignExceptions()
    res['S1_two_machine_sign_exceptions'] = {'exceptions': bad, 'evaluated': tot}
    res['S2_symbolic'] = checkRaceBoundSymbolic()
    res['S2_race'] = findRaceSummary()
    ev, resid = findRaceEquilibriumEigen()
    res['S2_post_shock_equilibrium'] = {'max_real_eigenvalue': float(np.max(ev.real)),
                                        'rhs_residual': resid}
    forced = findRaceSummary(dict(RACE, forceU=1.0))
    res['S2_race_forced_max_induction'] = {k: forced[k] for k in
                                           ('min_N', 'first_violation_t', 'max_N_minus_bound',
                                            'Nss_post', 'mass_error')}
    res['S2_ceiling_worst_excess'] = findRandomControlCeiling()
    us = np.linspace(0, 1, 101)
    ys = [findCoRegulatedYield(u) for u in us]
    res['L6_coregulated'] = {'u_best': float(us[int(np.argmax(ys))]),
                             'Y_at_0': ys[0], 'Y_best': float(max(ys)), 'Y_at_1': ys[-1]}
    res['L8_cotranslational_symbolic'] = checkCotranslationalSymbolic()
    a, b = findTwoClientSteady(0.0), findTwoClientSteady(0.3)
    res['L8_two_client'] = {'Y_sigma1_0': a['Y'].tolist(), 'Y_sigma1_0.3': b['Y'].tolist(),
                            'free_0': a['F'].tolist(), 'free_0.3': b['F'].tolist()}
    res['L8_margin_sensitivity'] = findMarginSensitivity()
    budgets = [0.5, 1.0, 2.0, 3.0, 4.0, 6.0, 10.0]
    front, free = findEnergyFront(budgets)
    res['L9_energy'] = {'budgets': budgets, 'margin': front, 'margin_without_energy_row': free}
    return res


def findJsonSafe(o):
    if isinstance(o, dict):
        return {str(k): findJsonSafe(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [findJsonSafe(v) for v in o]
    if isinstance(o, (np.floating, float)):
        return None if not np.isfinite(o) else float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return findJsonSafe(o.tolist())
    return o


if __name__ == '__main__':
    out = findJsonSafe(findAllResults())
    with open(os.path.join(HERE, 'results.json'), 'w') as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print(json.dumps(out, indent=1, sort_keys=True))
