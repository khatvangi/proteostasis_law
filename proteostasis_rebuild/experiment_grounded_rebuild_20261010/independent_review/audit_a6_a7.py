"""independent audit of EQUATIONS.md sections A6 (sequestration saddle-node) and A7 (maintenance).

all coefficients here are hypothetical and dimensionless. nothing in this file is an E. coli
estimate and nothing here is biological validation. run:

    python audit_a6_a7.py            # prints results, writes results_a6_a7.json
    python -m unittest test_audit_a6_a7 -v

A6 model as written in EQUATIONS.md:
    U     = L / (k_c M_free + k_a + mu)
    dA/dt = k_a U - (mu + delta) A
    M_free = M_tot - beta A                       (written form, can go negative)
    M_free = max(M_tot - beta A, 0)               (saturated form audited here)
notation: a0 = k_a + mu, a = k_c M_tot + a0, b = k_c beta, m = mu + delta, A_max = M_tot / beta.
"""

import json
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp


# ---------------------------------------------------------------------------
# symbolic part (sympy): exact domain restriction of the saddle-node
# ---------------------------------------------------------------------------
def deriveA6Symbolic():
    kc, Mt, be, ka, mu, de, L, A = sp.symbols("k_c M_tot beta k_a mu delta L A", positive=True)
    a0 = ka + mu
    a = kc * Mt + a0
    b = kc * be
    m = mu + de
    A_max = Mt / be

    # saddle-node of the written (unsaturated) field: f = 0 and df/dA = 0
    f = ka * L / (a - b * A) - m * A
    sn = sp.solve([f, sp.diff(f, A)], [A, L], dict=True)
    sn = [d for d in sn if sp.simplify(d[A] - a / (2 * b)) == 0]
    A_c, L_crit = sn[0][A], sn[0][L]

    # load at which the low branch reaches full titration (M_free = 0) on the saturated model
    L_sat = sp.simplify(m * A_max * a0 / ka)

    out = {
        "A_c": sp.simplify(A_c),
        "L_crit": sp.simplify(L_crit),
        "A_max": A_max,
        "L_sat": L_sat,
        # A_max - A_c: positive iff k_c M_tot > k_a + mu
        "A_max_minus_A_c": sp.factor(sp.simplify(A_max - A_c)),
        # distance from titration point to the pole of the written formula
        "pole_minus_A_max": sp.simplify(a / b - A_max),
        # L_crit - L_sat: a perfect square, so L_sat <= L_crit always
        "L_crit_minus_L_sat": sp.factor(sp.simplify(L_crit - L_sat)),
        # high-aggregate state on the saturated branch (M_free = 0)
        "A_hi": sp.solve(ka * L / a0 - m * A, A)[0],
    }
    # explicit checks that the factorizations say what the report says
    expected_gap = (kc * Mt - ka - mu) / (2 * kc * be)
    out["A_gap_ok"] = sp.simplify(out["A_max_minus_A_c"] - expected_gap) == 0
    expected_L_gap = m * (kc * Mt - ka - mu) ** 2 / (4 * kc * be * ka)
    out["L_gap_ok"] = sp.simplify(out["L_crit_minus_L_sat"] - expected_L_gap) == 0
    out["pole_ok"] = sp.simplify(out["pole_minus_A_max"] - a0 / (kc * be)) == 0
    return out


# ---------------------------------------------------------------------------
# numeric part: the saturated A6 model
# ---------------------------------------------------------------------------
def unpackA6(p):
    a0 = p["k_a"] + p["mu"]
    return dict(a0=a0, a=p["k_c"] * p["M_tot"] + a0, b=p["k_c"] * p["beta"],
                m=p["mu"] + p["delta"], A_max=p["M_tot"] / p["beta"])


def criticalLoads(p):
    q = unpackA6(p)
    L_crit = q["m"] * q["a"] ** 2 / (4 * q["b"] * p["k_a"])
    L_sat = q["m"] * q["A_max"] * q["a0"] / p["k_a"]
    return L_crit, L_sat


def fieldA6(A, L, p, saturated=True, delta_of_free=None):
    """dA/dt. delta_of_free(M_free) lets disaggregation depend on free chaperone (sensitivity test)."""
    M_free = p["M_tot"] - p["beta"] * A
    if saturated:
        M_free = max(M_free, 0.0)
    den = p["k_c"] * M_free + p["k_a"] + p["mu"]
    d = p["delta"] if delta_of_free is None else delta_of_free(M_free)
    return p["k_a"] * L / den - (p["mu"] + d) * A


def steadyStatesA6(L, p):
    """all steady states of the saturated model, from the two pieces in closed form.
    returns list of (A, stable) sorted by A."""
    q = unpackA6(p)
    roots = []
    # low piece, 0 <= A <= A_max:  b m A^2 - a m A + k_a L = 0
    disc = (q["a"] * q["m"]) ** 2 - 4 * q["b"] * q["m"] * p["k_a"] * L
    if disc >= 0:
        for s in (-1.0, 1.0):
            A = (q["a"] * q["m"] + s * np.sqrt(disc)) / (2 * q["b"] * q["m"])
            if 0.0 <= A <= q["A_max"]:
                den = q["a"] - q["b"] * A
                slope = p["k_a"] * L * q["b"] / den ** 2 - q["m"]
                roots.append((float(A), bool(slope < 0)))
    # high piece, A > A_max (M_free = 0): linear, slope -m
    A_hi = p["k_a"] * L / (q["a0"] * q["m"])
    if A_hi > q["A_max"]:
        roots.append((float(A_hi), bool(q["m"] > 0)))
    return sorted(roots)


def countSignChanges(L, p, n=200001, A_top=None):
    """brute-force cross-check of the closed-form root count: sign changes of dA/dt on a grid."""
    q = unpackA6(p)
    A_top = A_top or 5.0 * max(q["A_max"], p["k_a"] * L / (q["a0"] * q["m"]))
    grid = np.linspace(0.0, A_top, n)
    g = np.array([fieldA6(x, L, p) for x in grid])
    # compare g > 0 rather than sign(g): an exact zero on the grid would otherwise count twice
    return int(np.sum((g[1:] > 0) != (g[:-1] > 0)))


def integrateA6(L, p, A0=0.0, t_end=None, saturated=True, delta_of_free=None):
    q = unpackA6(p)
    t_end = t_end or 400.0 / q["m"]
    pole = q["a"] / q["b"]

    def rhs(t, y):
        return [fieldA6(y[0], L, p, saturated, delta_of_free)]

    # stop the written (unsaturated) model before its pole, where dA/dt -> +inf
    def nearPole(t, y):
        return y[0] - 0.999 * pole
    nearPole.terminal = not saturated   # the saturated model has no pole

    def titrated(t, y):
        return y[0] - q["A_max"]

    sol = solve_ivp(rhs, (0.0, t_end), [A0], events=[nearPole, titrated],
                    rtol=1e-10, atol=1e-12, max_step=t_end / 200)
    return dict(A_final=float(sol.y[0, -1]), t_final=float(sol.t[-1]),
                hit_pole=(not saturated) and len(sol.t_events[0]) > 0,
                t_titrated=float(sol.t_events[1][0]) if len(sol.t_events[1]) else None,
                M_free_final_written=float(p["M_tot"] - p["beta"] * sol.y[0, -1]))


def hysteresisSweep(p, n=161):
    """quasi-static load sweep up then down; each point integrated from the previous state."""
    L_crit, L_sat = criticalLoads(p)
    Ls = np.linspace(0.5 * L_sat, 1.5 * L_crit, n)
    up, A = [], 0.0
    for L in Ls:
        A = integrateA6(L, p, A0=A)["A_final"]
        up.append(A)
    down = []
    for L in Ls[::-1]:
        A = integrateA6(L, p, A0=A)["A_final"]
        down.append(A)
    down = down[::-1]
    A_max = unpackA6(p)["A_max"]
    # load at which the up-sweep first sits above titration, and down-sweep last sits above it
    L_up_jump = float(Ls[np.argmax(np.array(up) > A_max)])
    above = np.where(np.array(down) > A_max)[0]
    L_down_jump = float(Ls[above[0]]) if len(above) else None
    return dict(L_crit=L_crit, L_sat=L_sat, L_up_jump=L_up_jump, L_down_jump=L_down_jump,
                dL=float(Ls[1] - Ls[0]), A_up_max=float(max(up)))


def randomSweep(n=20000, seed=20261010):
    """random positive parameters; checks the topology claims of the saturated model."""
    rng = np.random.default_rng(seed)
    bad = {"no_root": 0, "more_than_3": 0, "stable_ne_unstable_plus_1": 0,
           "above_Lcrit_not_unique_high": 0, "fold_when_kcM_le_a0": 0}
    n_bistable = 0
    for _ in range(n):
        p = {k: float(10 ** rng.uniform(-2, 2)) for k in ("k_c", "M_tot", "beta", "k_a", "mu", "delta")}
        L_crit, L_sat = criticalLoads(p)
        L = float(L_crit * 10 ** rng.uniform(-2, 1))
        r = steadyStatesA6(L, p)
        ns = sum(1 for _, st in r if st)
        if len(r) == 0:
            bad["no_root"] += 1
        if len(r) > 3:
            bad["more_than_3"] += 1
        if ns != len(r) - ns + 1:
            bad["stable_ne_unstable_plus_1"] += 1
        if L > L_crit and not (len(r) == 1 and r[0][0] > unpackA6(p)["A_max"]):
            bad["above_Lcrit_not_unique_high"] += 1
        if len(r) == 3:
            n_bistable += 1
            if p["k_c"] * p["M_tot"] <= p["k_a"] + p["mu"]:
                bad["fold_when_kcM_le_a0"] += 1
    return dict(n=n, seed=seed, violations=bad, n_bistable=n_bistable)


# ---------------------------------------------------------------------------
# A7: moving quota boundary and native degradation
# ---------------------------------------------------------------------------
def deriveA7Symbolic():
    t = sp.symbols("t")
    ku, kdN, P = sp.symbols("k_u k_dN P", positive=True)
    N = sp.Function("N")(t)
    mu = sp.Function("mu")(t)
    Nreq = sp.Function("N_req")
    # with native degradation added (k_dN = 0 recovers A1.3)
    Ndot = P - (mu + ku + kdN) * N
    phi = N - Nreq(mu)
    phidot = sp.diff(phi, t).subs(sp.Derivative(N, t), Ndot)
    # on the quota face N = N_req(mu)
    phidot_face = sp.simplify(phidot.subs(N, Nreq(mu)))

    # growth-rate-dependent feasibility with proportional synthesis s = s0 mu, P = y s0 mu,
    # N_req = n0 + n1 mu: quota  y s0 mu >= (mu + k_u)(n0 + n1 mu)
    y, s0, n0, n1, m_ = sp.symbols("y s_0 n_0 n_1 mu_", positive=True)
    quad = sp.expand((m_ + ku) * (n0 + n1 * m_) - y * s0 * m_)        # feasible where quad <= 0
    disc = sp.discriminant(quad, m_)
    y_crit = sp.solve(sp.Eq(disc, 0), y)
    y_crit = [r for r in y_crit if sp.simplify(r - (sp.sqrt(n0) + sp.sqrt(n1 * ku)) ** 2 / s0) == 0]
    mu_star = sp.simplify(sp.solve(sp.diff(quad, m_), m_)[0].subs(y, y_crit[0]))
    return dict(phidot_face=phidot_face, quad=quad, y_crit=y_crit[0] if y_crit else None,
                mu_star=sp.simplify(sp.powsimp(sp.sqrtdenest(mu_star), force=True)),
                quota_mu0_with_kdN=sp.limit((m_ + ku + kdN) * (n0 + n1 * m_), m_, 0))


def movingQuotaCounterexample():
    """2-state toy: N with constant productive flux P; mu relaxes upward (nutrient upshift).
    P is set so A7.1 holds with equality at t = 0 (dN/dt = 0 on the face). required native
    N_req = n0 + n1 mu rises, so the state leaves V immediately although A7.1 held."""
    ku, n0, n1, mu0, mu_t, r = 0.1, 1.0, 2.0, 0.5, 1.0, 1.0
    Nreq = lambda mu: n0 + n1 * mu
    P = (mu0 + ku) * Nreq(mu0)                           # A7.1 holds with equality at x0
    rhs = lambda t, x: [P - (x[1] + ku) * x[0], r * (mu_t - x[1])]
    sol = solve_ivp(rhs, (0, 2.0), [Nreq(mu0), mu0], dense_output=True, rtol=1e-10, atol=1e-12)
    ts = np.linspace(0, 2.0, 401)
    N, mu = sol.sol(ts)
    gap = N - Nreq(mu)
    return dict(params=dict(k_u=ku, n0=n0, n1=n1, mu0=mu0, mu_target=mu_t, relax=r, P=P),
                A71_slack_at_t0=P - (mu0 + ku) * Nreq(mu0),
                phidot_at_t0=-n1 * r * (mu_t - mu0),
                min_gap=float(gap.min()), gap_t_0p1=float(sol.sol(0.1)[0] - Nreq(sol.sol(0.1)[1])))


def feasibleGrowthInterval(y, s0, n0, n1, ku, s_basal=0.0):
    """mu values with y (s_basal + s0 mu) >= (mu + k_u)(n0 + n1 mu)."""
    # n1 mu^2 + (n0 + n1 ku - y s0) mu + (ku n0 - y s_basal) <= 0
    c = [n1, n0 + n1 * ku - y * s0, ku * n0 - y * s_basal]
    r = np.roots(c)
    r = np.sort(r[np.abs(r.imag) < 1e-12].real)
    if len(r) < 2:
        return None
    lo, hi = max(r[0], 0.0), r[1]
    return (float(lo), float(hi)) if hi >= lo and hi >= 0 else None


# ---------------------------------------------------------------------------
# parameter sets used in the report (hypothetical)
# checks.py uses a=3, b=1, k_a=1, mu+delta=1 without fixing how a splits into k_c M_tot and
# k_a+mu. the two sets below share those four numbers and differ only in that split.
# ---------------------------------------------------------------------------
P_FOLD = dict(k_c=1.0, M_tot=1.75, beta=1.0, k_a=1.0, mu=0.25, delta=0.75)   # k_c M_tot > k_a+mu
P_NOFOLD = dict(k_c=1.0, M_tot=1.10, beta=1.0, k_a=1.0, mu=0.90, delta=0.10)  # k_c M_tot < k_a+mu


def main():
    res = {}
    s6 = deriveA6Symbolic()
    res["A6_symbolic"] = {k: str(v) for k, v in s6.items()}

    for name, p in (("fold", P_FOLD), ("nofold", P_NOFOLD)):
        L_crit, L_sat = criticalLoads(p)
        block = dict(params=p, L_crit=L_crit, L_sat=L_sat, A_max=unpackA6(p)["A_max"],
                     A_c=unpackA6(p)["a"] / (2 * unpackA6(p)["b"]))
        for f in (0.9, 0.99, 1.1):
            L = f * L_crit
            block[f"steady_states_at_{f}_Lcrit"] = [(round(A, 6), st) for A, st in steadyStatesA6(L, p)]
            block[f"sign_changes_at_{f}_Lcrit"] = countSignChanges(L, p)
        L = 1.1 * L_crit
        block["integrate_saturated_1.1Lcrit"] = integrateA6(L, p)
        block["integrate_written_1.1Lcrit"] = integrateA6(L, p, saturated=False)
        block["A_hi_closed_form_1.1Lcrit"] = p["k_a"] * L / ((p["k_a"] + p["mu"]) * (p["mu"] + p["delta"]))
        res[f"A6_{name}"] = block

    res["A6_hysteresis_fold"] = hysteresisSweep(P_FOLD)
    res["A6_random_sweep"] = randomSweep()

    # zero loss on the titrated branch: the only true runaway
    p0 = dict(P_FOLD, mu=0.0, delta=0.0)
    res["A6_zero_loss"] = integrateA6(2.0 * criticalLoads(dict(P_FOLD))[0], p0, t_end=200.0)

    # disaggregation that needs free chaperone: delta = d_deg + d_dis * M_free / M_tot
    for d_deg in (0.0, 0.05):
        p = dict(P_FOLD)
        dfun = lambda Mf, d_deg=d_deg: d_deg + 0.7 * Mf / p["M_tot"]
        L = 1.1 * criticalLoads(p)[0]
        r = integrateA6(L, p, delta_of_free=dfun)
        r["A_hi_expected"] = p["k_a"] * L / ((p["k_a"] + p["mu"]) * (p["mu"] + d_deg))
        res[f"A6_delta_needs_free_chaperone_d_deg_{d_deg}"] = r

    s7 = deriveA7Symbolic()
    res["A7_symbolic"] = {k: str(v) for k, v in s7.items()}
    res["A7_moving_quota_counterexample"] = movingQuotaCounterexample()
    # feasibility in mu with proportional synthesis (hypothetical numbers)
    s0, n0, n1, ku = 10.0, 1.0, 2.0, 0.1
    y_c = (np.sqrt(n0) + np.sqrt(n1 * ku)) ** 2 / s0
    res["A7_growth_interval"] = dict(
        params=dict(s0=s0, n0=n0, n1=n1, k_u=ku), y_crit=float(y_c),
        mu_star=float(np.sqrt(ku * n0 / n1)),
        interval_at_1p5_ycrit=feasibleGrowthInterval(1.5 * y_c, s0, n0, n1, ku),
        interval_at_1p01_ycrit=feasibleGrowthInterval(1.01 * y_c, s0, n0, n1, ku),
        interval_at_0p99_ycrit=feasibleGrowthInterval(0.99 * y_c, s0, n0, n1, ku),
        interval_basal0p2_at_1p5_ycrit=feasibleGrowthInterval(1.5 * y_c, s0, n0, n1, ku, s_basal=0.2),
    )

    out = Path(__file__).with_name("results_a6_a7.json")
    out.write_text(json.dumps(res, indent=2, default=str))
    print(json.dumps(res, indent=2, default=str))
    return res


if __name__ == "__main__":
    main()
