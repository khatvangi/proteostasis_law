"""
bounded sweeps over unknown kinetic multipliers (no fitting) and the
negative controls. writes outputs/*.json and outputs/*.tsv.

  python run.py            # 200 one-client + 200 two-client configs
"""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

import model as m

OUT = os.path.join(m.HERE, 'outputs')
SEED = 20261010
N1 = int(os.environ.get("N1", 200))                     # one-client configs
N2 = int(os.environ.get("N2", 200))                     # two-client configs
CAP_GRID = [0.1, 0.3, 1.0, 3.0, 10.0]   # all pools scaled together
GROEL_OVER = 5.0             # overproduction fold (not reported in our records)
R_GAS = 8.314e-3             # kJ/(mol K)
EPS = 1e-9                   # relative change below this counts as zero

MEAS = m.loadMeasured()
MU, M0 = MEAS['mu'], MEAS['pools']


def sgn(new, old):
    d = new - old
    return 0 if abs(d) <= EPS * max(abs(old), 1e-300) else int(np.sign(d))


def steady(cl, M, **kw):
    st = m.solveSteady(cl, M, MU, **kw)
    if not st['ok']:
        raise RuntimeError('steady state not converged: %g' % st['resid'])
    return st


def Y1(cl, M, **kw):
    return steady(cl, M, **kw)['clients'][0]['Y']


def condition(cl, knock=(), scale=None, serial=False, M=None):
    c2, M2 = m.perturb(cl, M0 if M is None else M, knock, scale, serial)
    return steady(c2, M2)['clients'][0]


def cool(c, ea):
    """30 -> 23 C with unknown activation energies (kJ/mol) for folding,
    aggregation, machine kinetics and growth. used only to show whether the
    sign of cooling is identified."""
    q = {k: np.exp(-v / R_GAS * (1 / 296.15 - 1 / 303.15)) for k, v in
         ea.items()}
    d = m.copyClient(c)
    d['kf'] *= q['fold']
    d['ka'] *= q['agg']
    d['kon'] = d['kon'] * q['mach']
    for k in ['rT', 'fT', 'rK', 'fK', 'rG', 'fG', 'rH', 'fH', 'kh', 'kv',
              'dV', 'kb', 'rB']:
        d[k] *= q['mach']
    d['s'] *= q['grow']
    return d, MU * q['grow']


# ---------------------------------------------------------------------------
# part A: one client
# ---------------------------------------------------------------------------
def oneClient(i):
    rng = np.random.default_rng([SEED, 1, i])
    fH_zero = (i % 2 == 0)
    c = m.sampleClient(rng, MEAS, fH_zero)
    ea = {k: rng.uniform(0, 150) for k in ['fold', 'agg', 'mach', 'grow']}
    cl = [c]
    wt = condition(cl)
    out = dict(i=i, fH_zero=fH_zero, Y_wt=wt['Y'], A_wt=wt['A'])

    # theorem check at fixed wild-type free pools: route x helps iff
    # pi_x = a_x/(1-b_x) > Ymarg = P/(Q + 2 k' U)
    st = steady(cl, M0)
    F = st['F']
    agree = 0
    for j in range(4):
        if F[m.CAPTURE[j]] <= 0:
            continue
        c2 = m.copyClient(c)
        c2['kon'][j] *= 1 + 1e-6
        dY = m.clientSteady(c2, F, MU)['Y'] - m.clientSteady(c, F, MU)['Y']
        pred = np.sign(wt['pi'][j] - wt['Ymarg'])
        agree += int(np.sign(dY) == pred or abs(wt['pi'][j] - wt['Ymarg'])
                     < 1e-9)
    out['theorem_agree'] = agree

    # signs in the self-consistent finite-capacity model
    g2 = condition(cl, scale={m.G: 2.0})
    out['GroELx2_dY'] = sgn(g2['Y'], wt['Y'])
    out['GroELx2_dA'] = sgn(g2['A'], wt['A'])
    h0 = condition(cl, knock=(m.H,))
    out['dHtpG_dY'] = sgn(h0['Y'], wt['Y'])
    out['dHtpG_dA'] = sgn(h0['A'], wt['A'])
    for name, x in [('TF', m.T), ('DnaK', m.K), ('GroEL', m.G),
                    ('HslUV', m.V)]:
        out['d%s_dA' % name] = sgn(condition(cl, knock=(x,))['A'], wt['A'])

    # codon-conditioned folding hazard
    cf = [m.copyClient(c, kf=c['kf'] * 1.5)]
    wtf = condition(cf)
    out['kf_dY'] = sgn(wtf['Y'], wt['Y'])
    out['kf_dA'] = sgn(wtf['A'], wt['A'])
    depG = wt['Y'] - condition(cl, knock=(m.G,))['Y']
    depGf = wtf['Y'] - condition(cf, knock=(m.G,))['Y']
    depK = wt['Y'] - condition(cl, knock=(m.K,))['Y']
    depKf = wtf['Y'] - condition(cf, knock=(m.K,))['Y']
    out['depG'], out['depG_kf'] = depG, depGf
    out['kf_dDepG'] = sgn(depGf, depG) if depG > 0 else None
    out['kf_dDepK'] = sgn(depKf, depK) if depK > 0 else None

    # experimental directions (filters, not fits)
    dko = condition(cl, knock=(m.T, m.K))
    dko_g = condition(cl, knock=(m.T, m.K), scale={m.G: GROEL_OVER})
    out['E1a'] = dko['Y'] < wt['Y'] * (1 - EPS)
    out['E1b'] = dko_g['Y'] > dko['Y'] * (1 + EPS)
    out['E1b_aggregate_down'] = dko_g['A'] < dko['A'] * (1 - EPS)
    dkj = condition(cl, knock=(m.K,))
    dkj_h = condition(cl, knock=(m.K, m.H))
    dkj_v = condition(cl, knock=(m.K, m.V))
    out['E2a'] = dkj_h['Y'] > dkj['Y'] * (1 + EPS)
    out['E2b'] = dkj_v['Y'] > dkj['Y'] * (1 + EPS)
    # held out: insoluble (aggregate) in dKJ dHtpG vs dKJ, Fauvet 30 C
    out['D1_dA'] = sgn(dkj_h['A'], dkj['A'])
    out['D1_dA_pp'] = 100 * (dkj_h['A'] - dkj['A']) / MEAS['P_T']
    # pelletable-complex mapping: insoluble = A + CH + CV
    ins = lambda r: r['A'] + r['CH'] + r['CV']
    out['D1_dIns_with_complexes'] = sgn(ins(dkj_h), ins(dkj))
    out['dkj_CH_plus_CV'] = dkj['CH'] + dkj['CV']
    # second held-out direction: Calloni 30 C, dKJ dtig aggregates more
    out['D2_dA'] = sgn(condition(cl, knock=(m.K, m.T))['A'], dkj['A'])

    # serial (no direct GroEL bypass) variant of E1b
    s_dko = condition(cl, knock=(m.T, m.K), serial=True)
    s_dko_g = condition(cl, knock=(m.T, m.K), scale={m.G: GROEL_OVER},
                        serial=True)
    out['serial_E1b_dY'] = s_dko_g['Y'] - s_dko['Y']

    # cooling with unknown activation energies
    for tag, knock in [('wt', ()), ('dko', (m.T, m.K))]:
        c30, M30 = m.perturb(cl, M0, knock)
        c23 = [cool(c30[0], ea)[0]]
        mu23 = cool(c30[0], ea)[1]
        y30 = m.solveSteady(c30, M30, MU)['clients'][0]['Y']
        y23 = m.solveSteady(c23, M30, mu23)['clients'][0]['Y']
        out['cool_%s_dY' % tag] = sgn(y23, y30)

    # ClpB quota control: key signs at ClpB pool x0, x1, x10
    keys = []
    for fB in [1e-12, 1.0, 10.0]:
        base = condition(cl, scale={m.B: fB})
        keys.append((sgn(condition(cl, scale={m.B: fB, m.G: 2.0})['Y'],
                         base['Y']),
                     sgn(condition(cl, knock=(m.H,), scale={m.B: fB})['Y'],
                         base['Y']),
                     sgn(condition([m.copyClient(c, kf=c['kf'] * 1.5)],
                                   scale={m.B: fB})['Y'], base['Y'])))
    out['clpb_signs_invariant'] = len(set(keys)) == 1
    return out


# ---------------------------------------------------------------------------
# part B: two clients sharing finite pools
# ---------------------------------------------------------------------------
def twoClient(i):
    rng = np.random.default_rng([SEED, 2, i])
    cl = [m.sampleClient(rng, MEAS, bool(rng.random() < 0.5))
          for _ in range(2)]
    out = dict(i=i, fH1_zero=cl[0]['fH'] == 0.0, fH2_zero=cl[1]['fH'] == 0.0)
    eff = {'dHtpG': [], 's2x2': [], 'kf2x2': [], 'GroELx2': []}
    for cap in CAP_GRID:
        M = M0 * cap
        y = Y1(cl, M)
        eff['dHtpG'].append(sgn(Y1(*m.perturb(cl, M, knock=(m.H,))), y))
        eff['s2x2'].append(sgn(Y1([cl[0], m.copyClient(cl[1], s=2 *
                                   cl[1]['s'])], M), y))
        eff['kf2x2'].append(sgn(Y1([cl[0], m.copyClient(cl[1], kf=2 *
                                    cl[1]['kf'])], M), y))
        eff['GroELx2'].append(sgn(Y1(*m.perturb(cl, M, scale={m.G: 2.0})),
                                  y))
    out['eff'] = eff
    # held-out D1 with two clients: total aggregate, dKJ dHtpG vs dKJ
    a = lambda st: sum(r['A'] for r in st['clients'])
    out['D1_dA'] = sgn(a(steady(*m.perturb(cl, M0, knock=(m.K, m.H)))),
                       a(steady(*m.perturb(cl, M0, knock=(m.K,)))))
    # the one-client direct effect of HtpG loss on client 1 (client 2 absent)
    out['dHtpG_alone'] = [sgn(Y1(*m.perturb([cl[0]], M0 * cap,
                                             knock=(m.H,))),
                              Y1([cl[0]], M0 * cap)) for cap in CAP_GRID]

    # pooled-capacity comparator at measured pools
    def both(pooled):
        kw = dict(pooled=True, Mref=M0) if pooled else {}
        y = Y1(cl, M0, **kw)
        r = {}
        c2, M2 = m.perturb(cl, M0, knock=(m.H,))
        r['dHtpG'] = sgn(Y1(c2, M2, **kw), y)
        r['GroELx2'] = sgn(Y1(cl, M0 * np.array([1, 1, 2, 1, 1, 1.]), **kw),
                           y)
        c2, M2 = m.perturb(cl, M0, knock=(m.T, m.K), scale={m.G: GROEL_OVER})
        c3, M3 = m.perturb(cl, M0, knock=(m.T, m.K))
        r['dko_GroEL'] = sgn(Y1(c2, M2, **kw), Y1(c3, M3, **kw))
        swap = M0.copy()
        dlt = 0.5 * M0[m.G]                  # move sites DnaK -> GroEL
        swap[m.K] -= dlt
        swap[m.G] += dlt
        r['swap'] = sgn(Y1(cl, swap, **kw), y)
        return r
    out['resolved'] = both(False)
    out['pooled'] = both(True)
    return out


def summarize1(rows):
    n = len(rows)
    S = {'n': n}

    def tally(key, sub=None):
        rr = [r[key] for r in rows if (sub is None or sub(r)) and
              r[key] is not None]
        return {'+': sum(v == 1 for v in rr), '0': sum(v == 0 for v in rr),
                '-': sum(v == -1 for v in rr), 'n': len(rr)}
    S['theorem_agree_routes'] = sum(r['theorem_agree'] for r in rows)
    for k in ['GroELx2_dY', 'GroELx2_dA', 'kf_dY', 'kf_dA', 'kf_dDepG',
              'kf_dDepK', 'dTF_dA', 'dDnaK_dA', 'dGroEL_dA', 'dHslUV_dA',
              'dHtpG_dA', 'D1_dA', 'D1_dIns_with_complexes', 'D2_dA',
              'cool_wt_dY', 'cool_dko_dY']:
        S[k] = tally(k)
    S['dHtpG_dY_fH0'] = tally('dHtpG_dY', lambda r: r['fH_zero'])
    S['dHtpG_dY_fHpos'] = tally('dHtpG_dY', lambda r: not r['fH_zero'])
    for k in ['E1a', 'E1b', 'E2a', 'E2b']:
        S[k] = sum(r[k] for r in rows)
    ok = lambda r: r['E1a'] and r['E1b'] and r['E2a'] and r['E2b']
    S['pass_all_filters'] = sum(ok(r) for r in rows)
    S['pass_fH0'] = sum(ok(r) for r in rows if r['fH_zero'])
    S['pass_fHpos'] = sum(ok(r) for r in rows if not r['fH_zero'])
    for k in ['GroELx2_dY', 'kf_dDepG', 'D1_dA', 'D2_dA', 'cool_dko_dY',
              'cool_wt_dY']:
        S['filtered_' + k] = tally(k, ok)
    S['filtered_E1b_aggregate_down'] = sum(r['E1b_aggregate_down'] for r in
                                           rows if ok(r))
    S['serial_E1b_max_abs_dY'] = max(abs(r['serial_E1b_dY']) for r in rows)
    S['clpb_signs_invariant'] = sum(r['clpb_signs_invariant'] for r in rows)
    pp = [r['D1_dA_pp'] for r in rows]
    S['D1_dA_pp_range'] = [min(pp), float(np.median(pp)), max(pp)]
    S['dkj_CH_plus_CV_max_uM'] = max(r['dkj_CH_plus_CV'] for r in rows)
    return S


def pattern(v):
    nz = set(x for x in v if x != 0)
    if not nz:
        return 'zero'
    if nz == {1}:
        return 'always+'
    if nz == {-1}:
        return 'always-'
    return 'reverses'


def summarize2(rows):
    S = {'n': len(rows)}
    for k in ['dHtpG', 's2x2', 'kf2x2', 'GroELx2']:
        pats = [pattern(r['eff'][k]) for r in rows]
        S[k] = {p: pats.count(p) for p in ['always+', 'always-', 'reverses',
                                            'zero']}
    # HtpG loss lowers client-1 yield although it helps client 1 alone
    flip = [r for r in rows if r['fH1_zero']]
    S['dHtpG_fH1zero_n'] = len(flip)
    S['dHtpG_fH1zero_any_negative'] = sum(-1 in r['eff']['dHtpG']
                                          for r in flip)
    S['dHtpG_alone_any_negative_fH1zero'] = sum(-1 in r['dHtpG_alone']
                                                for r in flip)
    S['D1_dA'] = {v: sum(r['D1_dA'] == v for r in rows) for v in (1, 0, -1)}
    for k in ['dHtpG', 'GroELx2', 'dko_GroEL', 'swap']:
        S['pooled_vs_resolved_' + k] = {
            'disagree': sum(r['resolved'][k] != r['pooled'][k] for r in rows),
            'resolved_nonzero': sum(r['resolved'][k] != 0 for r in rows),
            'pooled_nonzero': sum(r['pooled'][k] != 0 for r in rows)}
    return S


def main():
    os.makedirs(OUT, exist_ok=True)
    procs = min(60, os.cpu_count() or 1)
    with Pool(procs) as pool:
        r1 = pool.map(oneClient, range(N1))
        r2 = pool.map(twoClient, range(N2))
    s1, s2 = summarize1(r1), summarize2(r2)
    meta = dict(seed=SEED, N1=N1, N2=N2, cap_grid=CAP_GRID,
                groel_over=GROEL_OVER, ranges=m.RANGES, eps=EPS,
                pools_uM=dict(zip(m.MACHINES, M0.tolist())), mu=MU,
                s_P=MEAS['s_P'], P_T=MEAS['P_T'], GroES7=MEAS['GroES7'],
                input_sha256=MEAS['hashes'])
    json.dump(dict(meta=meta, one_client=s1, two_client=s2),
              open(os.path.join(OUT, 'summary.json'), 'w'), indent=1,
              default=lambda o: o.tolist() if hasattr(o, 'tolist') else
              (bool(o) if isinstance(o, np.bool_) else str(o)))
    json.dump(dict(one=r1, two=r2), open(os.path.join(OUT, 'rows.json'), 'w'),
              default=lambda o: o.tolist() if hasattr(o, 'tolist') else
              (bool(o) if isinstance(o, np.bool_) else int(o)
               if isinstance(o, np.integer) else str(o)))
    print(json.dumps(dict(one_client=s1, two_client=s2), indent=1,
                     default=str))


if __name__ == '__main__':
    sys.exit(main())
