"""Two direct confrontations of the route architecture with existing experiments.

(1) Vorderwulbecke 2004 (FEBS Lett 559:181-187): GroEL/GroES overproduction restores
    growth of delta-tig delta-dnaK52 E. coli at 30 C. Test whether the imposed
    serial upstream requirement makes that rescue impossible in this model, and
    whether the no_serial variant admits it.
(2) The all-client radial normalization requires productive service of the aggregate
    class too, whose only productive route consumes the very small ClpB pool.
    Quantify how much that depresses the demand scale.
No biological rates are fitted here; these are structural LP facts.
"""
import json
from pathlib import Path
import numpy as np
from model import architecture, draw, evaluate, radial_capacity, pools, MACHINES

HERE = Path(__file__).resolve().parent
MULTS = (1., 2., 5., 10., 100.)

def rescue(n=200, seed=20261010):
    out = {'groel_multipliers': list(MULTS), 'full': {}, 'no_serial': {}, 'radial': {}}
    rows = {k: {m: [] for m in MULTS} for k in ('full', 'no_serial')}
    lam_all, lam_noagg, agg_share = [], [], []
    for i in range(n):
        shares, C, costs, _ = draw(seed, i)
        net = architecture(costs)
        lam = radial_capacity(net, shares, C)
        d = .9 * lam * shares
        # same demand and the same TF/DnaK ablation in both architectures
        ko = C.copy(); ko[[0, 1]] = 0.
        for control, key in (('full', 'full'), ('no_serial', 'no_serial')):
            alt = architecture(costs, control)
            for m in MULTS:
                c = ko.copy(); c[4] *= m; c[5] *= m   # GroEL and GroES together
                rows[key][m].append(evaluate(alt, d, c, 'architecture')['essential_fraction'])
        # radial scale with and without requiring productive service of the aggregate class
        keep = np.array([0, 1, 2, 4])
        s2 = shares.copy(); s2[3] = 0.
        lam_all.append(lam); lam_noagg.append(radial_capacity(net, s2, C)); agg_share.append(shares[3])
    for key in ('full', 'no_serial'):
        out[key] = {f'{m:g}x': dict(median_essential_fraction=float(np.median(rows[key][m])),
                                    trials_with_any_essential_flux=int(sum(v > 1e-9 for v in rows[key][m])),
                                    trials=len(rows[key][m])) for m in MULTS}
    lam_all = np.array(lam_all); lam_noagg = np.array(lam_noagg)
    ratio = lam_noagg / lam_all
    out['radial'] = dict(trials=n,
        lambda_all_classes_q05_50_95=np.quantile(lam_all, [.05, .5, .95]).tolist(),
        lambda_without_aggregate_q05_50_95=np.quantile(lam_noagg, [.05, .5, .95]).tolist(),
        ratio_q05_50_95=np.quantile(ratio, [.05, .5, .95]).tolist(),
        trials_ratio_above_2=int((ratio > 2).sum()),
        aggregate_share_q05_50_95=np.quantile(agg_share, [.05, .5, .95]).tolist(),
        clpb_pool_uM=float(pools()[MACHINES.index('ClpB')]),
        dnak_pool_uM=float(pools()[MACHINES.index('DnaK')]))
    (HERE / 'confront.json').write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))

if __name__ == '__main__':
    rescue()


def kappa_sweep(n=200, seed=20261010, kappas=(1., 2., 5., 10., 20., 50.)):
    """Compatible alternative: a GroEL-only route for obligate clients that costs
    kappa times the normal GroEL/GroES service. Required behaviour, from existing
    experiments: at baseline GroEL capacity the TF+DnaK double loss still fails
    (Genevaux 2004 synthetic lethality at 30 C), while GroEL/GroES overproduction
    rescues it (Vorderwulbecke 2004). Report which kappa reproduces both.
    """
    import numpy as np
    from model import Network, CLASSES, MACHINES
    res = {}
    for kappa in kappas:
        fail_base, rescued_10x = 0, 0
        for i in range(n):
            shares, C, costs, _ = draw(seed, i)
            net = architecture(costs)
            d = .9 * radial_capacity(net, shares, C) * shares
            cols, hs, prod, sinks, atps, names = list(net.A.T), list(net.H.T), list(net.productive), list(net.sink), list(net.atp), list(net.names)
            for k in (0, 1):   # the two obligate classes
                a = np.zeros(len(MACHINES)); a[[4, 5]] = kappa * costs[[4, 5], k]
                h = np.zeros(len(CLASSES)); h[k] = 1
                cols.append(a); hs.append(h); prod.append(1.); sinks.append('rescue')
                atps.append(float(np.array([0,1,0,0,7,0,2,4,0.]) @ a))
                names.append(CLASSES[k] + ':GroEL>GroES_direct:rescue')
            alt = Network(np.array(cols).T, np.array(hs).T, np.array(prod, float), np.array(sinks), np.array(atps), names)
            ko = C.copy(); ko[[0, 1]] = 0.
            base = evaluate(alt, d, ko, 'architecture')['essential_fraction']
            over = ko.copy(); over[4] *= 10; over[5] *= 10
            up = evaluate(alt, d, over, 'architecture')['essential_fraction']
            fail_base += base < .999
            rescued_10x += up > .999
        res[f'kappa_{kappa:g}'] = dict(trials=n, trials_double_KO_still_short=int(fail_base),
                                       trials_rescued_by_10x_GroEL=int(rescued_10x))
    p = HERE / 'confront.json'
    data = json.loads(p.read_text()); data['kappa_sweep'] = res
    p.write_text(json.dumps(data, indent=2) + '\n')
    print(json.dumps(res, indent=2))
