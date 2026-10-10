"""Finite-capacity route LP; no legacy model imports or stochastic route edges."""
from dataclasses import dataclass
import csv
from pathlib import Path
import numpy as np
from scipy.optimize import linprog

MACHINES = ('TF', 'DnaK', 'DnaJ', 'GrpE', 'GroEL', 'GroES', 'ClpB', 'PROTEASE', 'STORAGE')
CLASSES = ('essential_obligate', 'other_obligate', 'overlap', 'aggregate', 'other')
PARAMS = ('pool_TF', 'pool_DnaK', 'pool_DnaJ2', 'pool_GrpE2', 'pool_GroEL14', 'pool_GroES7', 'pool_ClpB6')
ROOT = Path(__file__).resolve().parents[2]
SOLVER_OPTIONS = {'primal_feasibility_tolerance': 1e-10, 'dual_feasibility_tolerance': 1e-10, 'ipm_optimality_tolerance': 1e-10}

@dataclass
class Network:
    A: np.ndarray
    H: np.ndarray
    productive: np.ndarray
    sink: np.ndarray
    atp: np.ndarray
    names: list


def pools():
    with (ROOT / 'WO-07/bundles.tsv').open() as f:
        rows = {r['param']: float(r['value']) for r in csv.DictReader(f, delimiter='\t')
                if r['bundle'] == 'EXPONENTIAL' and r['match_flag'] == 'MATCHED'
                and r['param'] in PARAMS}
    return np.array([rows[p] for p in PARAMS])


def architecture(costs=None, control='full'):
    """Fixed illustrative class routes, informed but NOT quantitatively measured.

    remove_obligate adds a finite KJE bypass ONLY to the two obligate classes.
    remove_upstream adds a finite downstream GroEL route to overlap clients.
    no_serial removes upstream stages before GroEL in the obligate classes.
    """
    if costs is None:
        costs = np.ones((len(MACHINES), len(CLASSES)))
    tf, kje, gro, recovery = (0,), (1, 2, 3), (4, 5), (6, 1, 2, 3)
    routes = [[tf + gro, kje + gro], [tf + gro, kje + gro],
              [tf, kje], [recovery], [tf, kje, gro]]
    if control == 'remove_obligate':
        routes[0].append(kje); routes[1].append(kje)
    elif control == 'remove_upstream':
        routes[2].append(gro)
    elif control == 'no_serial':
        routes[0] = [gro]; routes[1] = [gro]
    elif control != 'full':
        raise ValueError(control)
    cols, hs, prod, sinks, atps, names = [], [], [], [], [], []
    # Illustrative ATP per stage, not fitted hydrolysis measurements.
    stage_atp = np.array([0, 1, 0, 0, 7, 0, 2, 4, 0.])
    for k, rr in enumerate(routes):
        for route, sink in [(r, 'rescue') for r in rr] + [((7,), 'degradation'), ((8,), 'sequestration')]:
            a = np.zeros(len(MACHINES)); a[list(route)] = costs[list(route), k]
            h = np.zeros(len(CLASSES)); h[k] = 1
            cols.append(a); hs.append(h); prod.append(sink == 'rescue'); sinks.append(sink)
            atps.append(float(stage_atp @ a))
            names.append(CLASSES[k] + ':' + '>'.join(MACHINES[m] for m in route) + ':' + sink)
    return Network(np.array(cols).T, np.array(hs).T, np.array(prod, dtype=float),
                   np.array(sinks), np.array(atps), names)


def constraints(net, C, kind):
    """Nested relaxations: full subset independent subset scalar.

    Independent clients each get a separate COPY of each capacity; a global
    scalar work bound is retained. Copies are an explicit optimistic null.
    """
    if kind == 'scalar':
        return net.A.sum(axis=0, keepdims=True), np.array([C.sum()])
    if kind == 'independent':
        blocks = [net.A * h for h in net.H]
        return np.vstack([net.A.sum(axis=0), *blocks]), np.r_[C.sum(), np.tile(C, len(net.H))]
    if kind == 'architecture':
        return net.A, C
    raise ValueError(kind)


def solve(net, d, C, kind, reward, extra=()):
    """Max reward*x subject to Hx<=d and resources; u=d-Hx exactly.

    Returns a checked primal/dual certificate for every solve. Extra rows
    implement exact lexicographic objectives (up to solver tolerance).
    """
    B, b = constraints(net, C, kind)
    B = np.vstack([B, net.H]); b = np.r_[b, d]
    if extra:
        B = np.vstack([B, *[r for r, _ in extra]])
        b = np.r_[b, [v for _, v in extra]]
    res = linprog(-np.asarray(reward), A_ub=B, b_ub=b, bounds=(0, None), method='highs', options=SOLVER_OPTIONS)
    if not res.success:
        raise RuntimeError(res.message)
    x = res.x; prices = -res.ineqlin.marginals
    primal = float(reward @ x); dual = float(b @ prices)
    error = max(float(np.max(B @ x - b)), float(np.max(reward - B.T @ prices)),
                float(max(0, -prices.min())), abs(primal - dual), float(max(0, -x.min())))
    if error > 2e-9 * max(1, abs(primal)):
        raise AssertionError(('certificate', error))
    return x, primal, prices, error


def evaluate(net, d, C, kind):
    essential = net.H[0] * net.productive
    x, ev, _, e1 = solve(net, d, C, kind, essential)
    fix = [(-essential, -ev)]
    x, folded, prices, e2 = solve(net, d, C, kind, net.productive, fix)
    fix.append((-net.productive, -folded))
    x, handled, _, e3 = solve(net, d, C, kind, np.ones(len(x)), fix)
    fix.append((-np.ones(len(x)), -handled))
    x, _, _, e4 = solve(net, d, C, kind, -net.atp, fix)
    u = d - net.H @ x
    if np.min(u) < -2e-9 or (kind == 'architecture' and np.max(net.A @ x - C) > 2e-9):
        raise AssertionError('conservation violation')
    return dict(essential_fraction=ev/d[0] if d[0] else 1., folded=folded,
                folded_fraction=folded/d.sum(), degradation=float(x[net.sink == 'degradation'].sum()),
                sequestration=float(x[net.sink == 'sequestration'].sum()),
                unserved=max(0., float(u.sum())), atp=float(net.atp @ x),
                essential_shortfall=float(ev < d[0] - 1e-8),
                certificate_error=max(e1,e2,e3,e4))


def radial_capacity(net, shares, C):
    """Maximum lambda delivering lambda*shares productively for EVERY class."""
    A = net.A[:, net.productive.astype(bool)]
    H = net.H[:, net.productive.astype(bool)]
    n = A.shape[1]
    res = linprog(np.r_[np.zeros(n), -1.], A_ub=np.c_[A, np.zeros(len(C))], b_ub=C,
                  A_eq=np.c_[H, -shares], b_eq=np.zeros(len(shares)), bounds=(0,None), method='highs', options=SOLVER_OPTIONS)
    if not res.success:
        raise RuntimeError(res.message)
    return float(res.x[-1])


def draw(seed, trial, width=1.):
    # Trial stream independent of execution order, model variant and trial count.
    rng = np.random.default_rng(np.random.SeedSequence([seed, trial]))
    obligate = 10**rng.uniform(np.log10(.005), np.log10(.35))
    essential_share = rng.uniform(.02, .4)  # flux fraction within 85, NOT 13/85
    aggregate = 10**rng.uniform(-3, -1)
    overlap = rng.uniform(.05, .45)
    shares = np.array([obligate*essential_share, obligate*(1-essential_share), overlap,
                       aggregate, 1-obligate-overlap-aggregate])
    multipliers = 10**rng.uniform(-width, width, 7)
    raw = pools()*multipliers
    C = np.r_[raw/raw.sum(), .05, .05]  # sink capacities synthetic; same all models
    costs = 10**rng.uniform(np.log10(.5), np.log10(2.), (9,5))
    return shares, C, costs, multipliers


def perturbations(C):
    out = {'baseline': C.copy()}
    # Exact equal absolute work loss without capacity creation or compensation.
    delta = .8 * min(C[1], C[4], C[6])
    for name, ix in [('DnaK',1),('GroEL',4),('ClpB',6)]:
        c=C.copy(); c[ix]-=delta; out['matched_'+name]=c
    for a,b in [(1,4),(1,6),(4,6),(0,1)]:
        delta_pair = .4 * min(C[a],C[b])
        for label, indices in [('a',[a]),('b',[b]),('ab',[a,b])]:
            c=C.copy()
            for ix in indices: c[ix]-=delta_pair
            out[f'pair_{MACHINES[a]}_{MACHINES[b]}_{label}']=c
        # Same total removal as ab, concentrated on one endpoint.
        for label, ix in [('aa',a),('bb',b)]:
            c=C.copy(); c[ix]-=2*delta_pair
            out[f'pair_{MACHINES[a]}_{MACHINES[b]}_{label}']=c
    for name, ids in [('GroEL',[4]),('TF',[0]),('DnaK',[1]),('TF_DnaK',[0,1])]:
        c=C.copy(); c[ids]=0; out['KO_'+name]=c
    return out
