"""constant-demand steady-state model (see EQUATIONS.md for numbering).

three levels under a fixed, externally specified demand:
  D   required function            [fu/min]
  Q   available functional capacity [fu/min],  Q = a N
  P   supply of functional protein  [molecules/min], which must offset loss lam*N

every parameter value in this file is HYPOTHETICAL. none is an E. coli estimate.
units: time in min, protein and machines in molecules (or complexes) per cell,
function in fu ("function units", e.g. substrate molecules converted).

run:  python model.py            (prints the worked examples)
"""
from itertools import combinations

import numpy as np
from scipy.optimize import linprog, root

TOL = 1e-9


# ---------------------------------------------------------------- scalar model
def steadyLevel(P, lam):
    """N* = P/lam  (E3). P in molecules/min, lam in 1/min."""
    return P / lam


def levelAt(t, P, lam, N0):
    """exact solution of dN/dt = P - lam N with constant P  (E2)."""
    Nstar = P / lam
    return Nstar + (N0 - Nstar) * np.exp(-lam * t)


def requiredSupply(D, a, lam, eps=0.0):
    """minimum functional supply that sustains Q >= D  (E5, E9).
    eps = molecules inactivated per fu delivered; eps = 0 means catalysis
    does not consume the catalyst, so the only cost is turnover lam."""
    return D * (lam / a + eps)


def scalarFeasible(D, a, Pmax, lam):
    """the scalar condition D <= a Pmax / lam  (E6)."""
    return D <= a * Pmax / lam


def timeToShortfall(N0, Nreq, Nstar, lam):
    """time for an installed level N0 > Nreq to decay to Nreq when the
    sustainable level Nstar < Nreq  (E8). inf if it never gets there."""
    if N0 <= Nreq:
        return 0.0
    if Nstar >= Nreq:
        return np.inf
    return np.log((N0 - Nstar) / (Nreq - Nstar)) / lam


def feedbackLevel(Pmax, beta, lam):
    """steady state when supply falls with N: P(N) = Pmax - beta N  (E10)."""
    return Pmax / (lam + beta)


# ------------------------------------------- routing: potential vs realized
# one client, two machines (K = a DnaK-like holder, G = a GroEL-like folder).
# routes a chain can take, with yield (native fraction) and machine cost
# (machine-min per chain = dwell time tau):
#   spont : no machine, yield kf/(kf+ka)
#   K     : K first (dwell tauK); fraction yK native, rest passed to G (serial)
#   G     : G directly (dwell tauG), yield yG
ROUTING = dict(s=100.0,                 # chains synthesized per min
               kf=0.3, ka=0.7,          # spontaneous fold / aggregate commit [1/min]
               yK=0.5, yG=0.9,          # native fraction after K, after G
               tauK=0.05, tauG=0.2,     # dwell per chain [min]
               Ktot=200.0, Gtot=40.0)   # installed machine counts (same in both cases)
GOOD = dict(bK=0.05, bG=0.02)    # capture rate constants [1/(molecule min)]
POOR = dict(bK=0.001, bG=0.002)  # same machines, weak client recognition

LAM, A_SPEC = 0.02, 600.0        # client loss [1/min], activity [fu/(molecule min)]
D_SCALAR = 1.8e6                 # demand [fu/min] -> Nreq 3000, Preq 60/min


def routeTable(p):
    """per-route yield and machine cost, shared by the LP and the kinetics."""
    return {'spont': (p['kf'] / (p['kf'] + p['ka']), {}),
            'K': (p['yK'] + (1 - p['yK']) * p['yG'],
                  {'K': p['tauK'], 'G': (1 - p['yK']) * p['tauG']}),
            'G': (p['yG'], {'G': p['tauG']})}


def scalarClient(p, lam=LAM):
    """routing example as a one-client problem for buildLP (machines installed)."""
    routes = list(routeTable(p).values())
    clients = [dict(name='E', lam=lam, length=1.0, smax=p['s'], routes=routes)]
    machines = {'K': dict(fixed=p['Ktot']), 'G': dict(fixed=p['Gtot'])}
    return clients, machines


def potentialSupply(p):
    """Pmax: best productive flux if routing were freely allocable (LP, E12)."""
    clients, machines = scalarClient(p)
    return maxProductive(clients, machines)


def realizedRouting(p, b):
    """P_real: routing set by kinetic competition for the newly made chain.
    hazards: kf, ka (spontaneous), bK*Kfree, bG*Gfree. free machine counts
    solve the occupancy balance free = tot - tau*flux (Little's law, E13)."""
    y = routeTable(p)

    def fluxes(Kf, Gf):
        hK, hG = b['bK'] * Kf, b['bG'] * Gf
        H = p['kf'] + p['ka'] + hK + hG
        return {'spont': p['s'] * (p['kf'] + p['ka']) / H,
                'K': p['s'] * hK / H, 'G': p['s'] * hG / H}

    def occupancy(J):
        return {m: sum(J[r] * y[r][1].get(m, 0.0) for r in J) for m in ('K', 'G')}

    def resid(z):
        occ = occupancy(fluxes(*z))
        return [z[0] - (p['Ktot'] - occ['K']), z[1] - (p['Gtot'] - occ['G'])]

    sol = root(resid, [p['Ktot'], p['Gtot']], tol=1e-12)
    assert sol.success and max(abs(np.array(resid(sol.x)))) < 1e-8
    J = fluxes(*sol.x)
    occ = occupancy(J)
    return dict(P=sum(J[r] * y[r][0] for r in J), J=J,
                occK=occ['K'] / p['Ktot'], occG=occ['G'] / p['Gtot'])


# ------------------------------------ multi-client fixed demand: LP + Farkas
def buildLP(clients, machines, A=None, D=None, R=None):
    """variables x = [route fluxes J_{i,r} (chains/min) ; maintained machines M_m].

    clients : dicts(name, lam [1/min], length [aa], routes=[(yield, {m: cost})],
              optional smax [chains/min])
    machines: m -> dict(fixed=count)                    installed, cost ignored
                or dict(lam=[1/min], length=[aa])       must be sustained (E16)
    A [fu/(molecule min)], D [fu/min]: requirement A N >= D with N = P/lam (E14)
    R [aa/min]: shared synthesis budget, includes machine maintenance (E17)
    returns Y (productive flux per client), F (= A Lam^-1 Y), G x <= h."""
    cols = [(i, r) for i, c in enumerate(clients) for r in range(len(c['routes']))]
    mvar = [m for m, s in machines.items() if 'fixed' not in s]
    nJ, nx = len(cols), len(cols) + len(mvar)

    Y = np.zeros((len(clients), nx))
    for k, (i, r) in enumerate(cols):
        Y[i, k] = clients[i]['routes'][r][0]

    rows, h = [], []
    for m, spec in machines.items():                   # occupancy (E13)
        row = np.zeros(nx)
        for k, (i, r) in enumerate(cols):
            row[k] = clients[i]['routes'][r][1].get(m, 0.0)
        if 'fixed' in spec:
            h.append(spec['fixed'])
        else:
            row[nJ + mvar.index(m)] = -1.0
            h.append(0.0)
        rows.append(row)
    for i, c in enumerate(clients):                    # per-client synthesis cap
        if 'smax' in c:
            row = np.zeros(nx)
            row[[k for k, (ii, _) in enumerate(cols) if ii == i]] = 1.0
            rows.append(row)
            h.append(c['smax'])
    if R is not None:                                  # shared budget (E17)
        row = np.zeros(nx)
        for k, (i, _) in enumerate(cols):
            row[k] = clients[i]['length']              # failed chains cost too
        for j, m in enumerate(mvar):
            row[nJ + j] = machines[m]['length'] * machines[m]['lam']
        rows.append(row)
        h.append(R)

    lp = dict(Y=Y, G=np.array(rows), h=np.array(h, float), cols=cols, mvar=mvar)
    if A is not None:
        lam = np.array([c['lam'] for c in clients])
        lp['F'] = np.atleast_2d(A) @ (Y / lam[:, None])
        lp['D'] = np.atleast_1d(np.asarray(D, float))
    return lp


def maxProductive(clients, machines, R=None, client=0):
    """largest productive flux of one client under the constraints (Pmax)."""
    lp = buildLP(clients, machines, R=R)
    res = linprog(-lp['Y'][client], A_ub=lp['G'], b_ub=lp['h'],
                  bounds=(0, None), method='highs')
    assert res.status == 0
    return -res.fun


def feasible(lp):
    """does some x >= 0 satisfy G x <= h and F x >= D?  (E15)"""
    res = linprog(np.zeros(lp['G'].shape[1]),
                  A_ub=np.vstack([lp['G'], -lp['F']]),
                  b_ub=np.concatenate([lp['h'], -lp['D']]),
                  bounds=(0, None), method='highs')
    assert res.status in (0, 2), res.message
    return res.status == 0, (res.x if res.status == 0 else None)


def certificate(lp):
    """search for Farkas prices (u >= 0 on resources, w >= 0 on demands) with
    G^T u >= F^T w and h.u < D.w  (E18). returns (u, w) or None."""
    nu, nw = len(lp['h']), len(lp['D'])
    c = np.concatenate([lp['h'], -lp['D']])
    A_ub = np.hstack([-lp['G'].T, lp['F'].T])          # -(G^T u - F^T w) <= 0
    res = linprog(c, A_ub=A_ub, b_ub=np.zeros(A_ub.shape[0]),
                  A_eq=np.ones((1, nu + nw)), b_eq=[1.0],
                  bounds=(0, None), method='highs')
    if res.status != 0 or res.fun >= -TOL:
        return None
    return res.x[:nu], res.x[nu:]


def checkCertificate(lp, u, w):
    """verify a certificate by hand, independent of the solver that found it.
    if it holds, no x >= 0 can meet the demand: for any feasible x,
    D.w <= w.Fx <= u.Gx <= u.h, contradicting h.u < D.w."""
    ok_cols = np.all(lp['G'].T @ u - lp['F'].T @ w >= -TOL)
    return bool(ok_cols and (lp['h'] @ u < lp['D'] @ w - TOL)
                and np.all(u >= -TOL) and np.all(w >= -TOL))


def hallSubsetCheck(need, touches, capacity, tau=1.0):
    """Hall-type counting: for every subset S of clients, total occupancy
    demand <= total capacity of the machines S touches. necessary for some
    networks, NOT sufficient once a route uses more than one machine."""
    n = len(need)
    for k in range(1, n + 1):
        for S in combinations(range(n), k):
            ms = set().union(*(touches[i] for i in S))
            if tau * sum(need[i] for i in S) > sum(capacity[m] for m in ms) + TOL:
                return False
    return True


# ------------------------------------------------- hypothetical vector example
def vectorExample(R=60000.0, lam_m=0.02):
    """two clients sharing a GroEL-like machine G; E2 also uses K in series.
    machines must be sustained against loss lam_m from the same budget R."""
    clients = [
        dict(name='E1', lam=0.02, length=400.0,
             routes=[(0.1, {}), (0.9, {'G': 0.2})]),
        dict(name='E2', lam=0.02, length=300.0,
             routes=[(0.6, {}), (0.95, {'K': 0.05, 'G': 0.1})]),
    ]
    machines = {'K': dict(lam=lam_m, length=640.0),
                'G': dict(lam=lam_m, length=8400.0)}
    A = np.array([[600.0, 0.0], [0.0, 300.0]])
    return clients, machines, A


def main():
    p = ROUTING
    Nreq = D_SCALAR / A_SPEC
    Preq = requiredSupply(D_SCALAR, A_SPEC, LAM)
    Pmax = potentialSupply(p)
    print('HYPOTHETICAL scalar example')
    print(f'  D={D_SCALAR:.3g} fu/min  a={A_SPEC} fu/(molecule min)  lam={LAM}/min')
    print(f'  Nreq={Nreq:.0f} molecules  Preq={Preq:.2f} molecules/min')
    print(f'  Pmax (LP, allocable routing) = {Pmax:.2f}  -> aPmax/lam = {A_SPEC*Pmax/LAM:.3g}')
    for tag, b in (('GOOD', GOOD), ('POOR', POOR)):
        r = realizedRouting(p, b)
        print(f'  {tag}: P_real={r["P"]:.2f}  Q_sus={A_SPEC*r["P"]/LAM:.3g}  '
              f'occK={100*r["occK"]:.1f}%  occG={100*r["occG"]:.1f}%  '
              f'route fluxes ' + ', '.join(f'{k}={v:.1f}' for k, v in r['J'].items()))
    rp = realizedRouting(p, POOR)
    Ns = steadyLevel(rp['P'], LAM)
    print(f'  POOR, start installed N0=5000: N*={Ns:.0f}, falls below Nreq after '
          f'{timeToShortfall(5000, Nreq, Ns, LAM):.1f} min')
    p2 = dict(p, Ktot=2 * p['Ktot'], Gtot=2 * p['Gtot'])
    print(f'  POOR with both pools doubled: P_real={realizedRouting(p2, POOR)["P"]:.2f}')

    print('HYPOTHETICAL vector example (two clients, shared G, maintained machines)')
    D0 = np.array([1.8e6, 0.9e6])
    for R, lam_m, scale in ((60000, 0.02, 1.0), (60000, 0.02, 1.5),
                            (48000, 0.02, 1.0), (48000, 0.0, 1.0)):
        cl, mc, A = vectorExample(R, lam_m)
        lp = buildLP(cl, mc, A, scale * D0, R)
        ok, x = feasible(lp)
        line = f'  R={R} aa/min lam_m={lam_m} D={scale}xD0: feasible={ok}'
        if not ok:
            u, w = certificate(lp)
            line += (f'  certificate verified={checkCertificate(lp, u, w)}  '
                     f'binding prices u={np.round(u, 4).tolist()} w={np.round(w, 4).tolist()}')
        print(line)


if __name__ == '__main__':
    main()
