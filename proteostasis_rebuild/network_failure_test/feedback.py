"""feedback extension: does overload of one module, acting on another, create
genuine cascade thresholds?

base fact (model.py): without feedback every constraint is linear in
(x, lambda), the feasible lambda set is [0, lambda*], and the unserved load
grows continuously past lambda*. a discontinuity (a jump in unserved load at
a load where nothing new has been added) or hysteresis needs a positive
feedback loop. two kinds are implemented, both as a mechanistic, myopic fixed
point (the cell does not re-optimize against its own feedback):

  allocation step, given demand D and effective capacity c_eff:
      min  sum_k w_k z_k + eps * cost.x
      s.t. sum_{r in k} x_r + z_k = D_k,  U^T x <= c_eff,  x, z >= 0
      z_k = unserved load of class k (w_k distinct, fixed: a triage order)

  'titration'  (overload-triggered)
      transfer: unserved folding load becomes aggregate load,
                D_AGG = lambda d_AGG + gamma * sum_{k folding} z_k
      capacity: unprocessed aggregates sequester DnaK,
                c_eff[DnaK] = max(0, c[DnaK] - beta * z_AGG)
  'utilization' (pre-overload, threshold theta < 1)
      transfer as above
      capacity: for each declared coupling (src -> tgt),
                c_eff[tgt] *= 1 - beta * clip((u_src - theta) / (1 - theta), 0, 1)
                u_src = usage_src / c_eff[src]

viability: all z = 0. a 'jump' is a discontinuity of the unserved fraction
Z(lambda) = sum z / sum D on the upward continuation, confirmed by bisection
refinement (the drop must survive a 2^-14 shrink of the interval).
"""
import numpy as np

from model import CLASSES, FOLDING, MACHINES, MIDX, solveCapLP, solveLamC

FOLD_IDX = [CLASSES.index(k) for k in FOLDING]
AGG = CLASSES.index("AGG")
DNAK = MIDX["DnaK"]
TRIAGE = 1.0 + 1e-3 * np.arange(len(CLASSES))   # fixed, distinct class weights


def solveServed(inst, D, c_eff):
    U, cls, cost = inst["U"], inst["cls"], inst["cost"]
    R, K = len(U), len(D)
    Aeq = np.zeros((K, R + K))
    for r in range(R):
        Aeq[cls[r], r] = 1.0
    Aeq[:, R:] = np.eye(K)
    obj = np.concatenate([1e-6 * cost / cost.max(), TRIAGE])
    res = solveCapLP(obj, U, np.maximum(c_eff, 0), Aeq, D, n_extra=K)
    if res.status != 0:
        raise RuntimeError(res.message)
    return res.x[:R], res.x[R:]


def findEffective(inst, lam, state, fb):
    """demand and effective capacity implied by the current state."""
    z, usage = state
    D = lam * inst["d"].copy()
    D[AGG] += fb["gamma"] * z[FOLD_IDX].sum()
    c = inst["c"].copy()
    if fb["kind"] == "titration":
        c[DNAK] = max(0.0, c[DNAK] - fb["beta"] * z[AGG])
    elif fb["kind"] == "utilization":
        base = inst["c"]
        for src, tgt in fb["couplings"]:
            u = usage[src] / base[src] if base[src] > 0 else 0.0
            g = np.clip((u - fb["theta"]) / (1 - fb["theta"]), 0, 1)
            c[tgt] *= 1 - fb["beta"] * g
    return D, np.maximum(c, 0.0)


def solveFixedPoint(inst, lam, fb, state=None, maxit=300, tol=1e-9):
    """iterate state -> allocation(state). undamped first; the step is halved
    whenever the residual stops shrinking (oscillation). convergence is judged
    on the residual of the undamped map, relative to the load scale."""
    K, M = len(inst["d"]), len(inst["c"])
    z, usage = state if state is not None else (np.zeros(K), np.zeros(M))
    eta, prev = 1.0, np.inf
    scale = max(lam * 1.0, 1e-300)
    for it in range(maxit):
        D, c = findEffective(inst, lam, (z, usage), fb)
        x, zn = solveServed(inst, D, c)
        un = inst["U"].T @ x
        res = (np.max(np.abs(zn - z)) + np.max(np.abs(un - usage))) / scale
        if res < tol:
            z, usage = zn, un
            D, _ = findEffective(inst, lam, (z, usage), fb)
            return (z, usage), float(z.sum() / D.sum()), True
        if res > 0.9 * prev:
            eta = max(eta / 2, 1 / 64)
        prev = res
        z = z + eta * (zn - z)
        usage = usage + eta * (un - usage)
    D, _ = findEffective(inst, lam, (z, usage), fb)
    return (z, usage), float(z.sum() / D.sum()), False


def findCascade(inst, fb, fmin=0.2, fmax=3.0, n=141, jump=0.1, refine=14, confirm=0.05):
    """upward then downward continuation in lambda/lambda*_C.
    returns onset of unserved load, confirmed jumps, hysteresis, convergence."""
    lamC = solveLamC(inst)
    grid = np.linspace(fmin, fmax, n)
    up, states, nonconv = [], [], 0
    st = None
    for f in grid:
        st, Z, ok = solveFixedPoint(inst, f * lamC, fb, st)
        up.append(Z)
        states.append(st)
        nonconv += not ok
    down = [None] * n
    st = states[-1]
    for i in range(n - 1, -1, -1):
        st, Z, ok = solveFixedPoint(inst, grid[i] * lamC, fb, st)
        down[i] = Z
        nonconv += not ok
    up, down = np.array(up), np.array(down)
    jumps = []
    for i in range(n - 1):
        if up[i + 1] - up[i] > jump:
            lo, hi, slo = grid[i], grid[i + 1], states[i]
            for _ in range(refine):          # bisect, continuing from the left state
                mid = 0.5 * (lo + hi)
                smid, Zm, _ = solveFixedPoint(inst, mid * lamC, fb, slo)
                if Zm - up[i] > jump / 2:
                    hi = mid
                else:
                    lo, slo = mid, smid
            _, Zlo, _ = solveFixedPoint(inst, lo * lamC, fb, slo)
            _, Zhi, _ = solveFixedPoint(inst, hi * lamC, fb, slo)
            if Zhi - Zlo > confirm:
                jumps.append({"at_frac": float(hi), "dZ": float(Zhi - Zlo)})
    onset = next((float(f) for f, Z in zip(grid, up) if Z > 1e-6), None)
    return {"onset_frac": onset, "n_jumps": len(jumps), "jumps": jumps,
            "hysteresis": float(np.max(np.abs(up - down))),
            "nonconverged_points": int(nonconv), "Z_up": up.tolist(), "Z_down": down.tolist(),
            "grid": grid.tolist()}


def makeFeedback(rng, kind):
    """broad random feedback parameters."""
    fb = {"kind": kind, "gamma": float(10 ** rng.uniform(-1, 1)),
          "beta": float(10 ** rng.uniform(-1, 2))}   # broad: give the cascade every chance
    if kind == "utilization":
        fb["beta"] = float(rng.uniform(0.1, 1.0))       # fractional capacity loss
        fb["theta"] = float(rng.uniform(0.3, 0.95))
        n = int(rng.integers(1, 5))
        pairs = set()
        while len(pairs) < n:
            s, t = rng.choice(len(MACHINES), 2, replace=False)
            pairs.add((int(s), int(t)))
        fb["couplings"] = sorted(pairs)
    return fb
