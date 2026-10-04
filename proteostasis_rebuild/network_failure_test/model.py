"""constrained-allocation model of proteostasis network capacity.

STRUCTURAL / SENSITIVITY model. not calibrated, not a biological prediction.

objects
  machines m = 1..M with capacity c_m = a_m * s_m
      a_m  measured glucose-exponential pool (pools.py), used only as a
           relative size
      s_m  per-unit throughput multiplier: UNMEASURED, drawn log-uniform over
           a declared width W (decades) and swept
  client classes k with demand lambda * d_k (sum d_k = 1, lambda = global load)
  routes r: a class can be resolved by any of its routes; route r consumes
      U[r, m] units of machine m's capacity per unit of client flux. a route may
      use several machines (DnaK + DnaJ + GrpE; ClpB + KJE; ClpX + ClpP)

viability at load lambda (problem i, feasibility)
      exists x >= 0:  sum_{r in k} x_r = lambda d_k  (every class fully served)
                      U^T x <= c                     (no machine over capacity)
maximal tolerable load (problem ii, single objective)
      lambda* = max lambda subject to the same constraints
              = max_x min_m c_m / usage_m(x)     (maximin over allocations)
the feasible set is a convex polyhedron containing (x=0, lambda=0) and scaling
x by t in [0,1] maps feasible to feasible, so {lambda: feasible} = [0, lambda*]:
(i) and (ii) are the same problem. LP duality gives the bottleneck certificate
      lambda* = min_{y>=0} c.y  s.t.  sum_k d_k min_{r in k} U[r].y >= 1
which for single-machine routes is the Hall / min-cut formula
      lambda* = min_{S subset of classes} c(N(S)) / d(S).

architectures, same machines, same capacities, same demand
  A  pooled bucket: the machine constraints are summed into one,
     sum_m usage_m <= sum_m c_m  (a relaxation of C, so lambda*_A >= lambda*_C)
  B  independent, non-overlapping: each class uses only its primary route and a
     shared machine is split statically in proportion to nominal load; no
     rerouting (a restriction of C, so lambda*_C >= lambda*_B).
     lambda*_B = min_m c_m / L_m,  L_m = sum_k d_k U[primary_k, m]
  C  overlapping network with rerouting: the full LP above
"""
import itertools

import numpy as np
from scipy.optimize import linprog

from pools import findPools

MACHINES = ["TF", "DnaK", "DnaJ", "GrpE", "GroEL", "GroES", "HtpG", "ClpB",
            "IbpA", "IbpB", "Lon", "ClpX", "ClpP", "HslU", "HslV"]
CLASSES = ["NASC", "DK", "GRO", "HTP", "DMG", "AGG"]
FOLDING = ["NASC", "DK", "GRO", "HTP"]
MIDX = {m: i for i, m in enumerate(MACHINES)}

# functional modules (machine sets). cofactors are separate machines so that
# a scarce cofactor (DnaJ2 at 0.2 uM) can be the binding constraint
MODULES = {
    "TF": ("TF",), "KJE": ("DnaK", "DnaJ", "GrpE"), "GRO": ("GroEL", "GroES"),
    "HTPG": ("HtpG",), "CLPB": ("ClpB",), "SEQ": ("IbpA", "IbpB"),
    "LON": ("Lon",), "CLPXP": ("ClpX", "ClpP"), "HSL": ("HslU", "HslV"),
}
PROTEASES = ["LON", "CLPXP", "HSL"]


# ---------------------------------------------------------------- instances

def findRoutesBio(rng):
    """random compatibility constrained to the declared structure.

    rules (assumptions, not data):
      DnaK broad: KJE resolves NASC, DK and HTP always; GRO rarely (p=0.2)
      GroEL subset: GRO always resolves GRO, DK with p=0.5, NASC with p=0.2
      TF: NASC always, DK with p=0.3
      HtpG: HTP always (with KJE attached with p=0.5)
      ClpB only on the aggregate/rescue branch: every ClpB route is an AGG
        route and needs KJE (bichaperone disaggregation)
      degradation is a separate sink: DMG uses proteases only (each p=0.7,
        at least one); AGG can be disaggregated then degraded (each p=0.3) or
        degraded directly by Lon (p=0.3)
      sequestration is a separate sink: AGG -> IbpA+IbpB with p=0.7
    the first route listed per class is its primary (canonical) route.
    returns list of (class, tuple of modules).
    """
    R = []
    p = rng.random
    R += [("NASC", ("TF",)), ("NASC", ("KJE",))]
    if p() < 0.2: R.append(("NASC", ("GRO",)))
    R.append(("DK", ("KJE",)))
    if p() < 0.5: R.append(("DK", ("GRO",)))
    if p() < 0.3: R.append(("DK", ("TF",)))
    R.append(("GRO", ("GRO",)))
    if p() < 0.2: R.append(("GRO", ("KJE",)))
    R.append(("HTP", ("HTPG", "KJE") if p() < 0.5 else ("HTPG",)))
    if p() < 0.5: R.append(("HTP", ("KJE",)))
    prots = [q for q in PROTEASES if p() < 0.7] or [PROTEASES[rng.integers(3)]]
    R += [("DMG", (q,)) for q in prots]
    R.append(("AGG", ("CLPB", "KJE")))
    if p() < 0.7: R.append(("AGG", ("SEQ",)))
    R += [("AGG", ("CLPB", "KJE", q)) for q in PROTEASES if p() < 0.3]
    if p() < 0.3: R.append(("AGG", ("LON",)))
    return R


def findRoutesUnstructured(rng):
    """adversarial null: 1-4 routes per class, each a random set of 1-3
    machines, no biological rules at all."""
    R = []
    for k in CLASSES:
        for _ in range(rng.integers(1, 5)):
            ms = rng.choice(MACHINES, size=rng.integers(1, 4), replace=False)
            R.append((k, tuple(f"@{m}" for m in ms)))
    return R


def makeInstance(seed, W=2.0, wU=1.0, alpha=1.0, abundance="measured",
                 structure="bio", pools=None):
    """one random network.
    W          width (decades) of the log-uniform capacity multiplier s_m
    wU         width (decades) of the log-uniform client-specific usage U[r,m]
    alpha      Dirichlet concentration of the demand mix d
    abundance  'measured' (a_m = glucose pool) or 'blind' (a_m = 1)
    """
    rng = np.random.default_rng(seed)
    routes = findRoutesBio(rng) if structure == "bio" else findRoutesUnstructured(rng)
    M, R = len(MACHINES), len(routes)
    U = np.zeros((R, M))
    cls = np.zeros(R, dtype=int)
    for i, (k, mods) in enumerate(routes):
        cls[i] = CLASSES.index(k)
        for mod in mods:
            for m in (MODULES[mod] if not mod.startswith("@") else (mod[1:],)):
                U[i, MIDX[m]] = 10 ** rng.uniform(-wU / 2, wU / 2)
    d = rng.dirichlet(alpha * np.ones(len(CLASSES)))
    pools = pools or findPools()
    a = np.array([pools[m] for m in MACHINES]) if abundance == "measured" else np.ones(M)
    s = 10 ** rng.uniform(-W / 2, W / 2, size=M)
    primary = np.array([int(np.flatnonzero(cls == k)[0]) for k in range(len(CLASSES))])
    # preference cost for overflow routing: primary cheapest, then list order;
    # small deterministic jitter breaks LP ties
    rank = np.zeros(R)
    for k in range(len(CLASSES)):
        idx = np.flatnonzero(cls == k)
        rank[idx] = np.arange(len(idx))
    cost = 1.0 + rank + 1e-3 * rng.random(R)
    return {"seed": seed, "routes": routes, "U": U, "cls": cls, "d": d,
            "c": a * s, "a": a, "s": s, "primary": primary, "cost": cost}


# ---------------------------------------------------------------- solvers

RETRIES = [0]   # per-process count of LP solves that needed a fallback method


def solveCapLP(obj, U, c, Aeq, beq, n_extra=0):
    """LP with machine-capacity rows U^T x <= c, numerically conditioned:
    each row is divided by its capacity (rhs = 1), and a machine at zero
    capacity is removed by fixing every route that uses it to 0. variables
    are R route fluxes followed by n_extra free-standing nonnegative columns.
    on solver failure the other HiGHS methods are tried and counted."""
    on = c > 0
    dead = (U[:, ~on] > 0).any(axis=1)
    A = np.hstack([(U[:, on] / c[on]).T, np.zeros((int(on.sum()), n_extra))])
    bounds = [(0, 0) if dz else (0, None) for dz in dead] + [(0, None)] * n_extra
    for method in ("highs", "highs-ipm", "highs-ds"):
        res = linprog(obj, A_ub=A if A.shape[0] else None, b_ub=np.ones(A.shape[0]) if A.shape[0] else None,
                      A_eq=Aeq, b_eq=beq, bounds=bounds, method=method)
        if res.status == 0:
            res.on = on
            return res
        RETRIES[0] += 1
    res.on = on
    return res


def solveLamC(inst, c=None, duals=False):
    """lambda*_C: max concurrent flow LP. variables (x_1..x_R, lambda)."""
    U, cls, d = inst["U"], inst["cls"], inst["d"]
    c = inst["c"] if c is None else c
    R, M = U.shape
    K = len(d)
    Aeq = np.zeros((K, R + 1))
    for r in range(R):
        Aeq[cls[r], r] = 1.0
    Aeq[:, R] = -d
    obj = np.zeros(R + 1)
    obj[R] = -1.0
    res = solveCapLP(obj, U, c, Aeq, np.zeros(K), n_extra=1)
    if res.status != 0:
        raise RuntimeError(f"LP failed: {res.message}")
    lam = max(0.0, -res.fun)
    if duals:
        # marginals of the capacity rows: d(lambda*)/d(c_m) >= 0 = bottleneck prices
        y = np.zeros(M)
        y[res.on] = -res.ineqlin.marginals / c[res.on]   # undo the row scaling
        return lam, res.x[:R], y
    return lam


def solveLamA(inst, c=None):
    """pooled bucket: summed constraint, each class uses its cheapest route."""
    c = inst["c"] if c is None else c
    U, cls, d = inst["U"], inst["cls"], inst["d"]
    per = np.array([U[cls == k].sum(axis=1).min() for k in range(len(d))])
    return float(c.sum() / (d @ per))


def solveLamB(inst, c=None):
    """fixed primary routes, static partition of shared machines."""
    c = inst["c"] if c is None else c
    L = inst["d"] @ inst["U"][inst["primary"]]
    on = L > 0
    return float(np.min(c[on] / L[on]))


SOLVERS = {"A": solveLamA, "B": solveLamB, "C": solveLamC}


# ---------------------------------------------------------------- knockdowns

def findKnockdowns(inst, arch, kappa):
    """W_i = lambda*(machine i at kappa*c_i)/lambda*_intact for every i, and
    W_ij for every pair. returns (lam0, W1 (M,), W2 dict (i,j)->W)."""
    f = SOLVERS[arch]
    lam0 = f(inst)
    M = len(inst["c"])
    W1 = np.empty(M)
    for i in range(M):
        c = inst["c"].copy()
        c[i] *= kappa
        W1[i] = f(inst, c) / lam0
    W2 = {}
    for i, j in itertools.combinations(range(M), 2):
        c = inst["c"].copy()
        c[i] *= kappa
        c[j] *= kappa
        W2[(i, j)] = f(inst, c) / lam0
    return lam0, W1, W2


def findPairStats(W1, W2, rho, tol=1e-9):
    """synthetic lethality and interaction statistics at operating point
    rho = lambda_0 / lambda*_intact (feasibility lost iff W < rho)."""
    sl = beyond_min = synergy = 0
    eps = []
    for (i, j), w in W2.items():
        if W1[i] >= rho and W1[j] >= rho and w < rho - tol:
            sl += 1
        if w < min(W1[i], W1[j]) - 1e-7:
            beyond_min += 1          # not explainable by the weakest link
        e = w - W1[i] * W1[j]        # multiplicative epistasis
        eps.append(e)
        if e < -1e-7:
            synergy += 1
    return {"synthetic_lethal_pairs": sl, "pairs_beyond_weakest_link": beyond_min,
            "synergistic_pairs": synergy, "mean_epistasis": float(np.mean(eps)),
            "min_epistasis": float(np.min(eps))}


def findStructuralCuts(inst, kmax=3):
    """minimal machine sets whose knockout (capacity 0) makes every route of
    some class unusable -> lambda*_C = 0 at ANY load. combinatorial, no LP."""
    U, cls = inst["U"], inst["cls"]
    used = [set(np.flatnonzero(U[r] > 0)) for r in range(len(U))]
    routesOf = [[used[r] for r in np.flatnonzero(cls == k)] for k in range(len(inst["d"]))]

    def kills(S):
        return any(all(rs & S for rs in rk) for k, rk in enumerate(routesOf) if inst["d"][k] > 0)

    cuts = []
    for size in range(1, kmax + 1):
        for S in itertools.combinations(range(U.shape[1]), size):
            S = set(S)
            if any(c <= S for c in cuts):
                continue                   # not minimal
            if kills(S):
                cuts.append(S)
    return [tuple(sorted(c)) for c in cuts]


def findStructuralCutsB(inst):
    """architecture B: every machine on any primary route is a single cut."""
    L = inst["d"] @ inst["U"][inst["primary"]]
    return [(i,) for i in np.flatnonzero(L > 0)]


# ---------------------------------------------------------------- saturation order

def solveMinCost(inst, lam, c=None, cap_scale=1.0):
    """overflow routing at load lam: cheapest (preferred) routes first.
    returns x or None if infeasible."""
    U, cls, d, cost = inst["U"], inst["cls"], inst["d"], inst["cost"]
    c = inst["c"] if c is None else c
    R = len(U)
    Aeq = np.zeros((len(d), R))
    for r in range(R):
        Aeq[cls[r], r] = 1.0
    res = solveCapLP(cost, U, c * cap_scale, Aeq, lam * d)
    return res.x if res.status == 0 else None


def findSaturationEvents(inst, n=200, tol=1e-7):
    """sweep lam/lambda*_C over (0, 1); record when each machine first becomes
    saturated (usage >= c (1 - tol)) under preference routing. distinct
    first-saturation levels before lambda* are soft events (rerouting kinks):
    viability is unaffected until lambda*_C."""
    lamC = solveLamC(inst)
    first = {}
    grid = np.arange(1, n) / n             # strictly below lambda*
    for f in grid:
        x = solveMinCost(inst, f * lamC)
        if x is None:                      # cannot happen below lambda* except numerically
            raise RuntimeError(f"min-cost LP failed below lambda* at f={f}")
        usage = inst["U"].T @ x
        for m in np.flatnonzero(usage >= inst["c"] * (1 - tol)):
            first.setdefault(int(m), float(f))
    levels = sorted(set(first.values()))
    order = sorted(first, key=lambda m: first[m])
    return {"n_events_before_failure": len(levels), "n_machines_saturated_before_failure": len(first),
            "first_event_frac": levels[0] if levels else None,
            "order": [(MACHINES[m], first[m]) for m in order]}


# ---------------------------------------------------------------- Pareto (optional objectives)

# ILLUSTRATIVE ATP cost per unit machine use. relative ordering only (ATP-free
# holdases < DnaK cycle < GroEL cycle < AAA+ unfoldases). NOT measured values.
ATP_ILLUSTRATIVE = {"TF": 0, "DnaK": 1, "DnaJ": 0, "GrpE": 0, "GroEL": 7, "GroES": 0,
                    "HtpG": 2, "ClpB": 50, "IbpA": 0, "IbpB": 0, "Lon": 50, "ClpX": 50,
                    "ClpP": 0, "HslU": 50, "HslV": 0}


def findParetoFront(inst, load_frac=0.5, n=25):
    """eps-constraint front at fixed load lam = load_frac * lambda*_C:
    minimize ATP subject to max utilization <= t, sweeping t from its minimum
    (= load_frac, the max-reserve point) to 1. also reports the sequestered
    (toxicity proxy) flux on each front point."""
    lamC = solveLamC(inst)
    atp = inst["U"] @ np.array([ATP_ILLUSTRATIVE[m] for m in MACHINES])
    seq = np.array([float(inst["U"][r, MIDX["IbpA"]] > 0) for r in range(len(inst["U"]))])
    U, cls, d = inst["U"], inst["cls"], inst["d"]
    Aeq = np.zeros((len(d), len(U)))
    for r in range(len(U)):
        Aeq[cls[r], r] = 1.0
    lam = load_frac * lamC
    pts = []
    for t in np.linspace(load_frac * (1 + 1e-9), 1.0, n):
        res = solveCapLP(atp, U, t * inst["c"], Aeq, lam * d)
        if res.status != 0:
            continue
        pts.append({"max_utilization": float(t), "reserve_margin": float(1 - t),
                    "atp": float(res.fun), "sequestered_flux": float(seq @ res.x),
                    "lam": float(lam)})
    return pts
