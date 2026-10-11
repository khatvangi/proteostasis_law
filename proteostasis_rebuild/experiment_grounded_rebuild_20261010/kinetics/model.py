"""
conservative, state-resolved kinetic model of client maintenance by the
E. coli chaperone/protease network, with finite, non-interchangeable
machine occupancy.

units
  concentrations: uM (per cell volume)
  first-order rates: 1/s ; bimolecular capture constants: 1/(uM s)
  aggregation constant ka: 1/(uM s)  (second order in free unfolded client)

per client i, states (all uM):
  U   free non-native client (aggregation-prone)
  CT  client held by trigger factor           (occupies TF)
  CK  client held by DnaK/DnaJ/GrpE (KJE)     (occupies DnaK)
  CG  client inside GroEL/GroES               (occupies GroEL14)
  CH  client held by HtpG                     (occupies HtpG2)
  CV  client engaged by HslUV for degradation (occupies HslUV)
  A   aggregated client, monomer equivalents
  CB  client being extracted by ClpB          (occupies ClpB6)
  NAT native client (productive pool)
  DEG cumulative degraded client   (sink)
  DIL cumulative diluted client    (sink, includes diluted native)

positive flows (all mass-action, all rate constants >= 0):
  source        s            -> U
  spontaneous   kf U         U  -> NAT          (codon-conditioned hazard)
  aggregation   ka U^2       U  -> A
  capture x     kon_x F_x U  U  -> Cx,  x in {T,K,G,H}
                (kon_G = 0 when the direct GroEL bypass is disabled)
  handoff       kh F_G CT    CT -> CG ;  kh F_G CK  CK -> CG
  release       rT CT, rK CK, rG CG, rH CH       -> U
  rescue        fT CT, fK CK, fG CG, fH CH       -> NAT
  disposal      kv F_V CH    CH -> CV ;  dV CV  CV -> DEG
  disaggregate  kb F_B A     A  -> CB ;  rB CB  CB -> U
  dilution      mu * every in-cell state          -> DIL

machine pools: F_x = M_x - sum_i C_x,i  (free sites). total pools are held
constant (balanced growth: machine synthesis balances dilution), so machine
conservation is exact by construction. pools are not interchangeable: a site
of one machine can never serve another machine's route.

client mass conservation: d/dt(sum of all eleven states) = s exactly.

outputs: productive yield Y = (native formation flux)/s at steady state;
aggregate A (uM) and aggregate share of total protein A/P_T.
there is no viability threshold anywhere in this file.
"""
import csv
import hashlib
import os

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root

HERE = os.path.dirname(os.path.abspath(__file__))
REBUILD = os.path.abspath(os.path.join(HERE, '..', '..'))
WO07 = os.path.join(REBUILD, 'WO-07', 'bundles.tsv')
ABUND = os.path.join(REBUILD, 'paper_reassessment_20261010', 'evidence',
                     'machine_abundance.tsv')

MACHINES = ['TF', 'DnaK', 'GroEL', 'HtpG', 'HslUV', 'ClpB']
T, K, G, H, V, B = range(6)
CAPTURE = [T, K, G, H]          # machines that capture free U directly
STATES = ['U', 'CT', 'CK', 'CG', 'CH', 'CV', 'A', 'CB', 'NAT', 'DEG', 'DIL']
NS = len(STATES)


def sha256(path):
    with open(path, 'rb') as fh:
        return hashlib.sha256(fh.read()).hexdigest()


def loadMeasured():
    """read the only measured inputs: machine pools, growth rate, synthesis
    rate and total protein. all from Schmidt 2016 BW25113 glucose exponential,
    as matched in WO-07. HslV/HslU are not in the WO-07 bundle; they come from
    the reassessment re-extraction of the same Schmidt tables."""
    rows = list(csv.DictReader(open(WO07), delimiter='\t'))
    exp = {(r['param'], r['candidate']): r for r in rows
           if r['bundle'] == 'EXPONENTIAL' and r['match_flag'] == 'MATCHED'}

    def val(p):
        return float(exp[(p, 'primary')]['value'])

    ab = {r['machine']: float(r['unit_pool_uM'])
          for r in csv.DictReader(open(ABUND), delimiter='\t')}
    # functional HslVU = one HslV12 + two HslU6; limiting unit sets the pool
    hslvu = min(ab['HslV'], ab['HslU'] / 2.0)
    pools = np.array([val('pool_TF'), val('pool_DnaK'), val('pool_GroEL14'),
                      val('pool_HtpG2'), hslvu, val('pool_ClpB6')])
    return {
        'pools': pools,                       # uM binding sites (1 per unit)
        'mu': val('growth_rate_mu'),          # 1/s
        's_P': val('synthesis_rate_s_P'),     # uM chains / s
        'P_T': val('total_protein_P_T'),      # uM chains
        'GroES7': val('pool_GroES7'),         # checked: exceeds GroEL14
        'hashes': {'WO-07/bundles.tsv': sha256(WO07),
                   'machine_abundance.tsv': sha256(ABUND)},
    }


# ---------------------------------------------------------------------------
# client given free pools: exact steady state (the 1-client analytic core)
# ---------------------------------------------------------------------------
def exitProbabilities(c, F, mu):
    """committed outcome probabilities of each held state.
    a_x: leaves to NAT without revisiting U; b_x: returns to U.
    1 - a_x - b_x: lost (degraded or diluted while not native)."""
    eG = c['fG'] + c['rG'] + mu
    aG, bG = c['fG'] / eG, c['rG'] / eG
    hand = c['kh'] * F[G]
    eT = c['rT'] + c['fT'] + hand + mu
    aT, bT = (c['fT'] + hand * aG) / eT, (c['rT'] + hand * bG) / eT
    eK = c['rK'] + c['fK'] + hand + mu
    aK, bK = (c['fK'] + hand * aG) / eK, (c['rK'] + hand * bG) / eK
    vH = c['kv'] * F[V]
    eH = c['rH'] + c['fH'] + vH + mu
    aH, bH = c['fH'] / eH, c['rH'] / eH
    eA = c['kb'] * F[B] + mu
    eB = c['rB'] + mu
    bA = (c['kb'] * F[B] / eA) * (c['rB'] / eB)
    return dict(a=np.array([aT, aK, aG, aH]), b=np.array([bT, bK, bG, bH]),
                bA=bA, eT=eT, eK=eK, eG=eG, eH=eH, eA=eA, eB=eB, hand=hand,
                vH=vH)


def clientSteady(c, F, mu):
    """exact steady state of one client for fixed free pools F.
    s = U (Q + k' U) with Q = kf + mu + sum theta_x (1-b_x), k' = ka (1-bA),
    P = kf + sum theta_x a_x, and Y = P U / s."""
    e = exitProbabilities(c, F, mu)
    theta = c['kon'] * F[CAPTURE]
    Q = c['kf'] + mu + np.sum(theta * (1.0 - e['b']))
    P = c['kf'] + np.sum(theta * e['a'])
    kp = c['ka'] * (1.0 - e['bA'])
    s = c['s']
    U = 2.0 * s / (Q + np.sqrt(Q * Q + 4.0 * kp * s))   # stable quadratic root
    CT = theta[0] * U / e['eT']
    CK = theta[1] * U / e['eK']
    CG = (theta[2] * U + e['hand'] * (CT + CK)) / e['eG']
    CH = theta[3] * U / e['eH']
    CV = e['vH'] * CH / (c['dV'] + mu)
    A = c['ka'] * U * U / e['eA']
    CB = c['kb'] * F[B] * A / e['eB']
    nat_flux = c['kf'] * U + c['fT'] * CT + c['fK'] * CK + c['fG'] * CG \
        + c['fH'] * CH
    occ = np.array([CT, CK, CG, CH, CV, CB])
    return dict(U=U, CT=CT, CK=CK, CG=CG, CH=CH, CV=CV, A=A, CB=CB, occ=occ,
                Y=nat_flux / s, Y_formula=P * U / s, P=P, Q=Q, kp=kp,
                Ymarg=P / (Q + 2.0 * kp * U), pi=e['a'] / np.maximum(
                    1.0 - e['b'], 1e-300), deg_flux=c['dV'] * CV)


# ---------------------------------------------------------------------------
# finite-capacity closure: F = M - occupancy(F)
# ---------------------------------------------------------------------------
def _occupancy(clients, Feff, mu):
    res = [clientSteady(c, Feff, mu) for c in clients]
    return res, sum(r['occ'] for r in res)


def solveSteady(clients, M, mu, pooled=False, Mref=None):
    """steady state with shared finite pools.
    resolved: each machine has its own pool M_x.
    pooled comparator: one shared pool sum(M); route x sees an effective free
    pool F_pool * Mref_x / sum(Mref) (Mref = wild-type pools), so at empty
    wild type it is identical to the resolved model; only the capacity
    accounting differs. a deleted machine must also have its route disabled
    (done by the caller through kon/kv/kb = 0)."""
    M = np.asarray(M, float)
    if pooled:
        Mref = np.asarray(Mref, float)
        Mpool = M.sum()
        w = Mref / Mref.sum()

        def gap(Fp):
            return Fp + _occupancy(clients, Fp * w, mu)[1].sum() - Mpool
        Fp = 0.0 if Mpool <= 0 else brentq(gap, 0.0, Mpool, xtol=1e-15,
                                           rtol=1e-14)
        res, occ = _occupancy(clients, Fp * w, mu)
        resid = abs(Fp + occ.sum() - Mpool) / max(Mpool, 1e-300)
        return dict(F=Fp * w, clients=res, resid=resid, ok=resid < 1e-9)

    live = M > 0

    def residual(z):
        F = np.zeros(6)
        F[live] = M[live] * np.clip(z, 0.0, 1.0)
        occ = _occupancy(clients, F, mu)[1]
        return (F[live] + occ[live] - M[live]) / M[live]

    # gauss-seidel warm start with brentq on each coordinate (each scalar
    # residual is -M at F_x = 0 and >= 0 at F_x = M, so a root is bracketed)
    F = M.copy()
    for _ in range(30):
        Fold = F.copy()
        for j in np.where(live)[0]:
            def g(x, j=j):
                Ft = F.copy()
                Ft[j] = x
                return x + _occupancy(clients, Ft, mu)[1][j] - M[j]
            F[j] = brentq(g, 0.0, M[j], xtol=1e-16, rtol=1e-14)
        if np.max(np.abs(F - Fold)[live] / M[live]) < 1e-12:
            break
    def polish(z0):
        sol = root(residual, z0, method='hybr', tol=1e-14)
        zs = np.clip(sol.x, 0, 1)
        r0, r1 = np.max(np.abs(residual(z0))), np.max(np.abs(residual(zs)))
        return (zs, r1) if r1 < r0 else (z0, r0)

    z, r = polish(F[live] / M[live])
    if r > 1e-10 and live.any():
        # fallback: relax the full ODE (slow modes such as HtpG -> HslUV
        # stall gauss-seidel), then polish from the relaxed free pools
        sol = integrate(clients, M, mu, 400.0 / mu)
        Yc = sol.y[:, -1].reshape(len(clients), NS)
        occ = np.array([Yc[:, 1].sum(), Yc[:, 2].sum(), Yc[:, 3].sum(),
                        Yc[:, 4].sum(), Yc[:, 5].sum(), Yc[:, 7].sum()])
        z2, r2 = polish(np.clip((M - occ)[live] / M[live], 0, 1))
        if r2 < r:
            z, r = z2, r2
    F = np.zeros(6)
    F[live] = M[live] * z
    res, occ = _occupancy(clients, F, mu)
    resid = float(np.max(np.abs(residual(z)))) if live.any() else 0.0
    return dict(F=F, clients=res, resid=resid, ok=resid < 1e-9)


# ---------------------------------------------------------------------------
# full ODE (used to verify the algebraic steady state, conservation, sign)
# ---------------------------------------------------------------------------
def rhs(t, y, clients, M, mu):
    n = len(clients)
    Y = y.reshape(n, NS)
    occ = np.zeros(6)
    occ[T] = Y[:, 1].sum()
    occ[K] = Y[:, 2].sum()
    occ[G] = Y[:, 3].sum()
    occ[H] = Y[:, 4].sum()
    occ[V] = Y[:, 5].sum()
    occ[B] = Y[:, 7].sum()
    F = M - occ                        # free pools, exact by construction
    d = np.zeros_like(Y)
    for i, c in enumerate(clients):
        U, CT, CK, CG, CH, CV, A, CB = Y[i, :8]
        cap = c['kon'] * F[CAPTURE] * U
        hT, hK = c['kh'] * F[G] * CT, c['kh'] * F[G] * CK
        fold = c['kf'] * U
        agg = c['ka'] * U * U
        disp = c['kv'] * F[V] * CH
        dis = c['kb'] * F[B] * A
        d[i, 0] = (c['s'] - fold - agg - cap.sum() + c['rT'] * CT
                   + c['rK'] * CK + c['rG'] * CG + c['rH'] * CH
                   + c['rB'] * CB - mu * U)
        d[i, 1] = cap[0] - hT - (c['rT'] + c['fT'] + mu) * CT
        d[i, 2] = cap[1] - hK - (c['rK'] + c['fK'] + mu) * CK
        d[i, 3] = cap[2] + hT + hK - (c['rG'] + c['fG'] + mu) * CG
        d[i, 4] = cap[3] - disp - (c['rH'] + c['fH'] + mu) * CH
        d[i, 5] = disp - (c['dV'] + mu) * CV
        d[i, 6] = agg - dis - mu * A
        d[i, 7] = dis - (c['rB'] + mu) * CB
        d[i, 8] = (fold + c['fT'] * CT + c['fK'] * CK + c['fG'] * CG
                   + c['fH'] * CH - mu * Y[i, 8])
        d[i, 9] = c['dV'] * CV
        d[i, 10] = mu * Y[i, :9].sum()
    return d.ravel()


def integrate(clients, M, mu, t_end, y0=None):
    n = len(clients)
    y0 = np.zeros(n * NS) if y0 is None else y0
    return solve_ivp(rhs, (0, t_end), y0, args=(clients, np.asarray(M, float),
                     mu), method='LSODA', rtol=1e-10, atol=1e-14,
                     dense_output=False)


# ---------------------------------------------------------------------------
# parameter sampling: every kinetic multiplier is unknown and sampled
# log-uniformly in units of the measured dilution rate mu. nothing is fitted.
# ---------------------------------------------------------------------------
RANGES = {                       # log10 bounds of dimensionless multipliers
    'kf': (-1, 3),               # kf / mu
    'agg': (-1, 4),              # ka (s/mu) / mu: aggregation hazard at U=s/mu
    'cap': (0, 4),               # kon_x M_x / mu: capture hazard, empty pool
    'exit': (0, 3),              # r_x / mu, f_x / mu, dV / mu, rB / mu
    'hand': (0, 3),              # kh M_G / mu
    'disp': (0, 3),              # kv M_V / mu
    'disagg': (-1, 3),           # kb M_B / mu
    'fH': (-1, 2),               # fH / mu when HtpG has a productive exit
    'phi': (-3, -1),             # client share of total synthesis s / s_P
}


def sampleClient(rng, meas, fH_zero):
    mu, M = meas['mu'], meas['pools']

    def lu(key):
        lo, hi = RANGES[key]
        return 10 ** rng.uniform(lo, hi)
    s = lu('phi') * meas['s_P']
    c = dict(s=s, kf=lu('kf') * mu, ka=lu('agg') * mu * mu / s)
    c['kon'] = np.array([lu('cap') * mu / M[x] for x in CAPTURE])
    for k in ['rT', 'fT', 'rK', 'fK', 'rG', 'fG', 'rH', 'dV', 'rB']:
        c[k] = lu('exit') * mu
    c['fH'] = 0.0 if fH_zero else lu('fH') * mu
    c['kh'] = lu('hand') * mu / M[G]
    c['kv'] = lu('disp') * mu / M[V]
    c['kb'] = lu('disagg') * mu / M[B]
    return c


def copyClient(c, **kw):
    d = {k: (v.copy() if isinstance(v, np.ndarray) else v)
         for k, v in c.items()}
    d.update(kw)
    return d


# ---------------------------------------------------------------------------
# perturbations. a deletion removes the pool AND (for the pooled comparator)
# the route; in the resolved model a zero pool already silences the route.
# ---------------------------------------------------------------------------
def perturb(clients, M, knock=(), scale=None, serial=False):
    M = np.asarray(M, float).copy()
    cl = [copyClient(c) for c in clients]
    for x in knock:
        M[x] = 0.0
    if scale:
        for x, f in scale.items():
            M[x] *= f
    for c in cl:
        for x in knock:
            if x in CAPTURE:
                c['kon'][CAPTURE.index(x)] = 0.0
            if x == G:
                c['kh'] = 0.0
            if x == V:
                c['kv'] = 0.0
            if x == B:
                c['kb'] = 0.0
        if serial:
            c['kon'][CAPTURE.index(G)] = 0.0
    return cl, M


def growthIndex(Y_ess, Y_ref, beta=1.0):
    """separately stated, uncalibrated growth relation: relative growth =
    (Y_ess / Y_ref)^beta with unknown beta > 0. only its sign is ever used;
    no viability threshold is defined."""
    return (Y_ess / Y_ref) ** beta
