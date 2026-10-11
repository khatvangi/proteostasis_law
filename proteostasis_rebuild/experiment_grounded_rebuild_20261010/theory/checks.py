"""verify every sign, threshold and counterexample stated in SCIENTIFIC_THEORY.md / EQUATIONS.md.

all parameter values here are illustrative and dimensionless. none is an e. coli estimate.
run:  python checks.py        (writes checks.json next to this file, exits non-zero on any failure)
"""
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, linprog

results = {}


def record(name, ok, detail):
    results[name] = {"pass": bool(ok), "detail": detail}


# ---------------------------------------------------------------------------
# 1. competing fates, linear (unsaturated) regime, solved exactly
# states: U (non-native soluble), N (native), A (aggregated). D is cumulative removal.
# exits from U: spontaneous fold kf; GroEL binding kG*G; HtpG binding kH*H; aggregation ka; dilution mu.
# GroEL-bound client -> N with prob pG, back to U with prob rG, to A with prob aG = 1-pG-rG.
# HtpG-bound client -> N with prob (1-q), to degradation with prob q.
# A -> U by disaggregation kdis; A diluted at mu. N diluted at mu.
# machine-bound complexes are treated as instantaneous (quasi-steady) in this section.
# ---------------------------------------------------------------------------
s, kf, kG, G, kH, H, ka, mu, kdis = sp.symbols("s k_f k_G G k_H H k_a mu k_dis", positive=True)
pG, rG, q = sp.symbols("p_G r_G q", nonnegative=True)
aG = 1 - pG - rG
U, N, A = sp.symbols("U N A")

eqU = s - (kf + kG * G + kH * H + ka + mu) * U + rG * kG * G * U + kdis * A
eqA = ka * U + aG * kG * G * U - (kdis + mu) * A
eqN = kf * U + pG * kG * G * U + (1 - q) * kH * H * U - mu * N
sol = sp.solve([eqU, eqA, eqN], [U, N, A], dict=True)[0]
Y = sp.simplify(mu * sol[N] / s)  # fraction of synthesised chains that end native (diluted as N)

# claimed closed form, aG = 0 case:  sign dY/dH = sign[(1-q)(ka_eff+mu) - q(kf + pG kG G)]
ka_eff = ka * mu / (kdis + mu)  # aggregation counted as loss only for the unrecovered fraction
Y0 = Y.subs(rG, 1 - pG)  # aG = 0
dYdH = sp.simplify(sp.diff(Y0, H))
claimed = (1 - q) * (ka_eff + mu) - q * (kf + pG * kG * G)
ratio = sp.factor(sp.simplify(dYdH / claimed))
# ratio must be strictly positive for every positive parameter set -> dY/dH and the claimed expression share sign
rng = np.random.default_rng(20261010)
free = [s, kf, kG, G, kH, H, ka, mu, kdis, pG, q]
ratio_f = sp.lambdify(free, ratio, "numpy")
draws = [ratio_f(*np.exp(rng.uniform(-4, 4, 9)), rng.uniform(0, 1), rng.uniform(0, 1)) for _ in range(20000)]
record(
    "htpg_sign_rule",
    min(draws) > 0,
    {"dY/dH divided by claimed sign expression": str(ratio), "min over 20000 log-uniform draws": float(min(draws))},
)

# numeric spot check: both signs are realised
def y_num(vals):
    return float(Y.subs(vals))

base = {s: 1, kf: 1.0, kG: 1.0, G: 1.0, kH: 1.0, H: 1.0, ka: 0.2, mu: 0.5, kdis: 0.0, pG: 0.9, rG: 0.1, q: 0.5}
harm = dict(base)
help_ = {**base, kf: 0.01, G: 0.01, ka: 5.0, q: 0.2}
dh = 1e-6
d_harm = (y_num({**harm, H: 1 + dh}) - y_num({**harm, H: 1})) / dh
d_help = (y_num({**help_, H: 1 + dh}) - y_num({**help_, H: 1})) / dh
record(
    "htpg_both_signs_realised",
    d_harm < 0 < d_help,
    {"dY/dH, strong alternative folding (kf=1,G=1,ka=0.2,q=0.5)": d_harm,
     "dY/dH, aggregation-dominated (kf=0.01,G=0.01,ka=5,q=0.2)": d_help},
)

# aggregate burden always falls with HtpG (aG = 0, any q)
A0 = sp.simplify(sol[A].subs(rG, 1 - pG))
dAdH = sp.simplify(sp.diff(A0, H))
dAdH_f = sp.lambdify(free, dAdH.subs(rG, 1 - pG), "numpy")
adraws = [dAdH_f(*np.exp(rng.uniform(-4, 4, 9)), rng.uniform(0, 1), rng.uniform(0, 1)) for _ in range(20000)]
record("htpg_lowers_aggregate", max(adraws) < 0,
       {"dA/dH": str(sp.factor(dAdH)), "max over 20000 draws": float(max(adraws))})

# GroEL: monotone beneficial when aG = 0; can be harmful when aG > 0
dYdG0 = sp.factor(sp.diff(Y0, G))
gro_cases = {
    "aG=0": {**base, rG: 0.1, pG: 0.9},
    "aG=0.6,pG=0.1 (unproductive exit to aggregate, fast spontaneous folder)": {**base, pG: 0.1, rG: 0.3, kf: 5.0, ka: 0.2},
}
gro_vals = {}
for k, v in gro_cases.items():
    gro_vals[k] = (y_num({**v, G: 1 + dh}) - y_num({**v, G: 1})) / dh
record("groel_sign", gro_vals["aG=0"] > 0 and list(gro_vals.values())[1] < 0,
       {"dY/dG": gro_vals, "symbolic dY/dG at aG=0": str(dYdG0)})

# ---------------------------------------------------------------------------
# 2. finite HtpG pool: sign of dY/dH_tot equals sign of dY/dH_free (single client)
# H_tot = H_f + C_H, C_H = kH H_f U * tau_H, tau_H = 1/(kout + mu)
# ---------------------------------------------------------------------------
def sat_yield(Htot, p):
    tau = 1.0 / (p["kout"] + p["mu"])

    def resid(Hf):
        Uss = p["s"] / (p["kf"] + p["pG"] * p["kG"] * p["G"] + p["kH"] * Hf + p["ka"] * p["mu"] / (p["kdis"] + p["mu"]) + p["mu"])
        return Htot - Hf - p["kH"] * Hf * Uss * tau

    Hf = brentq(resid, 0.0, Htot)
    num = p["kf"] + p["pG"] * p["kG"] * p["G"] + (1 - p["q"]) * p["kH"] * Hf
    den = p["kf"] + p["pG"] * p["kG"] * p["G"] + p["kH"] * Hf + p["ka"] * p["mu"] / (p["kdis"] + p["mu"]) + p["mu"]
    return num / den, Hf

psat = dict(s=50.0, kf=1.0, kG=1.0, G=1.0, kH=1.0, ka=0.2, mu=0.5, kdis=0.0, pG=0.9, q=0.5, kout=1.0)
hs = np.linspace(0.05, 20, 60)
ys, hfs = zip(*[sat_yield(h, psat) for h in hs])
record("saturation_preserves_sign",
       np.all(np.diff(hfs) > 0) and np.all(np.diff(ys) < 0),
       {"H_free increasing in H_tot": bool(np.all(np.diff(hfs) > 0)),
        "yield decreasing in H_tot (harm regime)": bool(np.all(np.diff(ys) < 0)),
        "fraction of H_tot bound at H_tot=0.05": 1 - hfs[0] / hs[0]})

# ---------------------------------------------------------------------------
# 3. necessary but not sufficient: capacity-feasible yet realised routing misses the quota
# one client class, synthesis s; routes: GroEL (productive, service time tauG), HtpG (disposal q),
# spontaneous. quota: native flux >= Qreq. pools GroEL M_G, HtpG M_H.
# capacity region F: jG*tauG <= M_G, jH*tauH <= M_H, jG + jH + jf + jloss = s, native = jG + jf + (1-q) jH >= Qreq
# ---------------------------------------------------------------------------
sF, tauG, tauH, MG, MH, Qreq, qd = 1.0, 2.0, 1.0, 5.0, 5.0, 0.8, 0.9
# lp: maximise native flux over F (jf limited by spontaneous fraction phi of s)
phi = 0.3
# variables jG, jH, jf, jloss
res = linprog(c=[-1, -(1 - qd), -1, 0],
              A_ub=[[tauG, 0, 0, 0], [0, tauH, 0, 0], [0, 0, 1, 0]],
              b_ub=[MG, MH, phi * sF],
              A_eq=[[1, 1, 1, 1]], b_eq=[sF], bounds=[(0, None)] * 4, method="highs")
lp_native = -res.fun
# realised kinetic routing: proportional to rates from U, unsaturated (pools large)
rate = {"G": 1.0, "H": 6.0, "f": 0.6, "loss": 0.2}
tot = sum(rate.values())
realised_native = (rate["G"] + rate["f"] + (1 - qd) * rate["H"]) / tot * sF
jG_real, jH_real = rate["G"] / tot * sF, rate["H"] / tot * sF
record("feasible_not_sufficient",
       lp_native >= Qreq and realised_native < Qreq and jG_real * tauG <= MG and jH_real * tauH <= MH,
       {"best native flux over capacity region F": lp_native, "quota": Qreq,
        "realised native flux": realised_native,
        "realised occupancies (GroEL, HtpG) / pool": [jG_real * tauG / MG, jH_real * tauH / MH]})

# deleting HtpG: LP optimum cannot rise (route removal shrinks F); realised yield rises
res2 = linprog(c=[-1, -(1 - qd), -1, 0],
               A_ub=[[tauG, 0, 0, 0], [0, tauH, 0, 0], [0, 0, 1, 0]],
               b_ub=[MG, 0.0, phi * sF], A_eq=[[1, 1, 1, 1]], b_eq=[sF], bounds=[(0, None)] * 4, method="highs")
rate2 = dict(rate, H=0.0)
tot2 = sum(rate2.values())
realised2 = (rate2["G"] + rate2["f"]) / tot2 * sF
record("deletion_can_help_only_kinetically",
       -res2.fun <= lp_native + 1e-12 and realised2 > realised_native,
       {"LP optimum with / without HtpG": [lp_native, -res2.fun],
        "realised native with / without HtpG": [realised_native, realised2]})

# ---------------------------------------------------------------------------
# 4. little's law occupancy identity in a saturable mass-action ode (two clients share one machine)
# ---------------------------------------------------------------------------
def ode(t, y, p):
    U1, U2, C1, C2, N1, N2 = y
    Mf = p["Mtot"] - C1 - C2
    b1, b2 = p["kon1"] * Mf * U1, p["kon2"] * Mf * U2
    return [p["s1"] - b1 - (p["kf1"] + p["mu"]) * U1 + p["koff"] * C1,
            p["s2"] - b2 - (p["kf2"] + p["mu"]) * U2 + p["koff"] * C2,
            b1 - (p["kcat1"] + p["koff"] + p["mu"]) * C1,
            b2 - (p["kcat2"] + p["koff"] + p["mu"]) * C2,
            p["kf1"] * U1 + p["kcat1"] * C1 - p["mu"] * N1,
            p["kf2"] * U2 + p["kcat2"] * C2 - p["mu"] * N2]

pl = dict(Mtot=1.0, kon1=5.0, kon2=2.0, s1=1.0, s2=0.7, kf1=0.3, kf2=0.1, mu=0.3, koff=0.5, kcat1=2.0, kcat2=0.8)
sol_ode = solve_ivp(ode, (0, 400), [0] * 6, args=(pl,), rtol=1e-10, atol=1e-12)
U1, U2, C1, C2, N1, N2 = sol_ode.y[:, -1]
Mf = pl["Mtot"] - C1 - C2
J1, J2 = pl["kon1"] * Mf * U1, pl["kon2"] * Mf * U2
tau1, tau2 = 1 / (pl["kcat1"] + pl["koff"] + pl["mu"]), 1 / (pl["kcat2"] + pl["koff"] + pl["mu"])
lhs = J1 * tau1 + J2 * tau2
record("littles_law_occupancy", abs(lhs - (C1 + C2)) < 1e-8 and lhs <= pl["Mtot"],
       {"sum J*tau": lhs, "bound complexes": C1 + C2, "pool": pl["Mtot"]})

# ---------------------------------------------------------------------------
# 5. reduced overload model (EQUATIONS A6) with a physical, capped chaperone pool
# dA/dt = ka L / (kc Mfree(A) + ka + mu) - (mu + delta) A,   Mfree(A) = max(Mtot - beta A, 0)
# notation: a0 = ka + mu, a = kc Mtot + a0, b = kc beta, m = mu + delta, A_max = Mtot / beta.
# the bound-client term sum C is omitted and delta is net removal (see A6). the old check
# "bistable_below_none_above" counted roots of the uncapped formula, whose upper root has
# negative free chaperone; it is removed and replaced by root counts of the capped model.
# ---------------------------------------------------------------------------
kc_, Mt_, be_, ka_, mu_, de_, L_, A_ = sp.symbols("k_c M_tot beta k_a mu delta L A", positive=True)
a0_ = ka_ + mu_
a_ = kc_ * Mt_ + a0_
b_ = kc_ * be_
m_ = mu_ + de_
f_low = ka_ * L_ / (a_ - b_ * A_) - m_ * A_          # valid only for A <= A_max
sn = sp.solve([f_low, sp.diff(f_low, A_)], [A_, L_], dict=True)
sn = [d for d in sn if sp.simplify(d[A_] - a_ / (2 * b_)) == 0][0]
A_max_s = Mt_ / be_
L_sat_s = m_ * A_max_s * a0_ / ka_                    # low branch reaches full titration
A_hi_s = sp.solve(ka_ * L_ / a0_ - m_ * A_, A_)[0]    # steady state on the titrated piece
sym_ok = {
    "A_c = a/(2b)": sp.simplify(sn[A_] - a_ / (2 * b_)) == 0,
    "L_crit = m a^2/(4 b ka)": sp.simplify(sn[L_] - m_ * a_**2 / (4 * b_ * ka_)) == 0,
    "A_max - A_c = (kc Mtot - ka - mu)/(2 kc beta)":
        sp.simplify(A_max_s - sn[A_] - (kc_ * Mt_ - ka_ - mu_) / (2 * kc_ * be_)) == 0,
    "L_crit - L_sat = m (kc Mtot - ka - mu)^2/(4 kc beta ka)":
        sp.simplify(sn[L_] - L_sat_s - m_ * (kc_ * Mt_ - ka_ - mu_) ** 2 / (4 * kc_ * be_ * ka_)) == 0,
    "pole a/b lies a0/(kc beta) beyond A_max": sp.simplify(a_ / b_ - A_max_s - a0_ / (kc_ * be_)) == 0,
    "A_hi = ka L/((ka+mu)(mu+delta))": sp.simplify(A_hi_s - ka_ * L_ / (a0_ * m_)) == 0,
    "loop width (L_crit - L_sat)/L_crit = ((kc Mtot - a0)/(kc Mtot + a0))^2":
        sp.simplify((sn[L_] - L_sat_s) / sn[L_] - ((kc_ * Mt_ - a0_) / (kc_ * Mt_ + a0_)) ** 2) == 0,
    "beta=0: unique steady state ka L/(a m)": (lambda r: len(r) == 1 and sp.simplify(r[0] - ka_ * L_ / (a_ * m_)) == 0)(
        sp.solve(ka_ * L_ / a_ - m_ * A_, A_)),
}
record("a6_symbolic_fold_and_cap", all(sym_ok.values()),
       {k: bool(v) for k, v in sym_ok.items()} | {"L_crit": str(sp.factor(sn[L_])), "L_sat": str(L_sat_s)})


def a6q(p):
    a0 = p["ka"] + p["mu"]
    return dict(a0=a0, a=p["kc"] * p["Mtot"] + a0, b=p["kc"] * p["beta"], m=p["mu"] + p["delta"],
                A_max=p["Mtot"] / p["beta"])


def a6_loads(p):
    q = a6q(p)
    return q["m"] * q["a"] ** 2 / (4 * q["b"] * p["ka"]), q["m"] * q["A_max"] * q["a0"] / p["ka"]


def a6_mfree(A, p):
    return max(p["Mtot"] - p["beta"] * A, 0.0)


def a6_field(A, L, p):
    return p["ka"] * L / (p["kc"] * a6_mfree(A, p) + p["ka"] + p["mu"]) - (p["mu"] + p["delta"]) * A


def a6_roots(L, p):
    """all steady states of the capped model in closed form: (A, stable), sorted."""
    q = a6q(p)
    out = []
    disc = (q["a"] * q["m"]) ** 2 - 4 * q["b"] * q["m"] * p["ka"] * L
    if disc >= 0:
        for sgn in (-1.0, 1.0):
            A = (q["a"] * q["m"] + sgn * np.sqrt(disc)) / (2 * q["b"] * q["m"])
            if 0.0 <= A <= q["A_max"]:
                slope = p["ka"] * L * q["b"] / (q["a"] - q["b"] * A) ** 2 - q["m"]
                out.append((float(A), bool(slope < 0)))
    A_hi = p["ka"] * L / (q["a0"] * q["m"])
    if A_hi > q["A_max"]:
        out.append((float(A_hi), q["m"] > 0))       # slope on the titrated piece is -m
    return sorted(out)


def a6_sign_changes(L, p, n=20001):
    """brute-force root count independent of the closed form."""
    q = a6q(p)
    top = 5.0 * max(q["A_max"], p["ka"] * L / (q["a0"] * q["m"]))
    g = np.array([a6_field(x, L, p) for x in np.linspace(0.0, top, n)])
    return int(np.sum((g[1:] > 0) != (g[:-1] > 0)))


def a6_integrate(L, p, A0=0.0, t_end=None):
    t_end = t_end or 400.0 / max(a6q(p)["m"], 1e-3)
    sol = solve_ivp(lambda t, y: [a6_field(y[0], L, p)], (0.0, t_end), [A0],
                    rtol=1e-10, atol=1e-12, max_step=t_end / 200)
    mfree = [a6_mfree(x, p) for x in sol.y[0]]
    return float(sol.y[0, -1]), float(min(sol.y[0])), float(min(mfree))


# hypothetical sets sharing a=3, b=1, ka=1, m=1 (the numbers of the removed check); they differ
# only in how a splits into kc Mtot and ka + mu. same sets as independent_review/audit_a6_a7.py.
P_FOLD = dict(kc=1.0, Mtot=1.75, beta=1.0, ka=1.0, mu=0.25, delta=0.75)    # kc Mtot > ka + mu
P_NOFOLD = dict(kc=1.0, Mtot=1.10, beta=1.0, ka=1.0, mu=0.90, delta=0.10)  # kc Mtot < ka + mu

Lc, Ls = a6_loads(P_FOLD)
counts, pattern, brute = {}, {}, {}
for fr in (0.9, 0.99, 1.1):
    r = a6_roots(fr * Lc, P_FOLD)
    counts[fr], pattern[fr] = len(r), [st for _, st in r]
    brute[fr] = a6_sign_changes(fr * Lc, P_FOLD)
r99 = a6_roots(0.99 * Lc, P_FOLD)
A_hi_11 = P_FOLD["ka"] * 1.1 * Lc / ((P_FOLD["ka"] + P_FOLD["mu"]) * (P_FOLD["mu"] + P_FOLD["delta"]))
record("a6_capped_root_counts_three_state_window",
       Ls < 0.99 * Lc < Lc and counts == {0.9: 1, 0.99: 3, 1.1: 1} and brute == counts
       and pattern[0.99] == [True, False, True] and r99[2][0] > a6q(P_FOLD)["A_max"]
       and abs(a6_roots(1.1 * Lc, P_FOLD)[0][0] - A_hi_11) < 1e-12,
       {"L_sat, L_crit": [Ls, Lc], "A_max": a6q(P_FOLD)["A_max"],
        "steady states (A, stable) at 0.9 / 0.99 / 1.1 L_crit": {str(k): a6_roots(k * Lc, P_FOLD) for k in counts},
        "brute-force sign changes": {str(k): v for k, v in brute.items()},
        "A_hi closed form at 1.1 L_crit": A_hi_11})

# no fold at any load when kc Mtot <= ka + mu
Lc_n, _ = a6_loads(P_NOFOLD)
n_nofold = [len(a6_roots(x, P_NOFOLD)) for x in np.linspace(0.01, 3.0, 300) * Lc_n]
record("a6_no_fold_without_capture_condition", set(n_nofold) == {1},
       {"kc Mtot, ka + mu": [P_NOFOLD["kc"] * P_NOFOLD["Mtot"], P_NOFOLD["ka"] + P_NOFOLD["mu"]],
        "distinct root counts over 300 loads in [0.01, 3] L_crit": sorted(set(n_nofold))})

# above L_crit the trajectory from A = 0 settles at the finite high state; free chaperone never < 0
A_end, A_min, Mf_min = a6_integrate(1.1 * Lc, P_FOLD)
record("a6_high_state_finite_and_attracting",
       abs(A_end - A_hi_11) < 1e-6 and Mf_min >= 0.0 and A_min >= 0.0,
       {"A(t_end) from A=0 at 1.1 L_crit": A_end, "A_hi closed form": A_hi_11,
        "min free chaperone along trajectory": Mf_min,
        "uncapped formula Mtot - beta A at the end (unphysical)": P_FOLD["Mtot"] - P_FOLD["beta"] * A_end})

# quasi-static sweep up then down: jumps at L_crit (up) and L_sat (down), within one grid step
Lgrid = np.linspace(0.5 * Ls, 1.5 * Lc, 161)
dL = float(Lgrid[1] - Lgrid[0])
up, A_now = [], 0.0
for Lv in Lgrid:
    A_now = a6_integrate(Lv, P_FOLD, A0=A_now)[0]
    up.append(A_now)
down = []
for Lv in Lgrid[::-1]:
    A_now = a6_integrate(Lv, P_FOLD, A0=A_now)[0]
    down.append(A_now)
up, down = np.array(up), np.array(down[::-1])
Amax_f = a6q(P_FOLD)["A_max"]
L_up = float(Lgrid[np.argmax(up > Amax_f)])
L_down = float(Lgrid[np.where(down > Amax_f)[0][0]])
record("a6_hysteresis_loop_between_Lsat_and_Lcrit",
       abs(L_up - Lc) <= dL and abs(L_down - Ls) <= dL
       and abs(np.max(up) - a6_roots(Lgrid[-1], P_FOLD)[-1][0]) < 1e-6,
       {"up-jump load": L_up, "L_crit": Lc, "down-jump load": L_down, "L_sat": Ls, "grid step": dL,
        "max A on sweep (finite)": float(np.max(up)),
        "A_hi closed form at top load": a6_roots(Lgrid[-1], P_FOLD)[-1][0], "loop width / L_crit": (Lc - Ls) / Lc})

# the only unbounded accumulation in the reduced model: no removal at all (mu = delta = 0)
p0 = dict(P_FOLD, mu=0.0, delta=0.0)
A0_end = a6_integrate(2.0 * Lc, p0, t_end=200.0)[0]
record("a6_unbounded_only_without_removal", A0_end > 100 * a6q(p0)["A_max"],
       {"A(200) with mu = delta = 0": A0_end, "A_max": a6q(p0)["A_max"]})

# local random sweep (small; the 20,000-draw sweep is in independent_review and cross-checked below)
rng6 = np.random.default_rng(7)
viol = {"no_root": 0, "more_than_3": 0, "stable_ne_unstable_plus_1": 0, "above_Lcrit_not_unique_high": 0,
        "three_roots_not_iff_condition": 0}
for _ in range(2000):
    p = {k: float(10 ** rng6.uniform(-2, 2)) for k in ("kc", "Mtot", "beta", "ka", "mu", "delta")}
    lc, ls = a6_loads(p)
    Lv = float(lc * 10 ** rng6.uniform(-2, 1))
    r = a6_roots(Lv, p)
    ns = sum(st for _, st in r)
    viol["no_root"] += len(r) == 0
    viol["more_than_3"] += len(r) > 3
    viol["stable_ne_unstable_plus_1"] += ns != len(r) - ns + 1
    viol["above_Lcrit_not_unique_high"] += Lv > lc and not (len(r) == 1 and r[0][0] > a6q(p)["A_max"])
    cond = p["kc"] * p["Mtot"] > p["ka"] + p["mu"] and ls < Lv < lc
    viol["three_roots_not_iff_condition"] += (len(r) == 3) != cond
record("a6_random_topology_local_2000", sum(viol.values()) == 0, {"violations": viol})

# cross-check representative numbers against the independent auditor's saved results
aud_path = Path(__file__).resolve().parent.parent / "independent_review" / "results_a6_a7.json"
aud = json.loads(aud_path.read_text())
af = aud["A6_fold"]
xcheck = {
    "L_crit": abs(af["L_crit"] - Lc) < 1e-12,
    "L_sat": abs(af["L_sat"] - Ls) < 1e-12,
    "roots at 0.99 L_crit": [x[0] for x in af["steady_states_at_0.99_Lcrit"]] == [round(x, 6) for x, _ in r99],
    "A_hi at 1.1 L_crit": abs(af["A_hi_closed_form_1.1Lcrit"] - A_hi_11) < 1e-12,
    "up / down jump loads": abs(aud["A6_hysteresis_fold"]["L_up_jump"] - L_up) < 1e-12
                            and abs(aud["A6_hysteresis_fold"]["L_down_jump"] - L_down) < 1e-12,
    "auditor 20,000-draw sweep: zero violations": aud["A6_random_sweep"]["n"] == 20000
                            and sum(aud["A6_random_sweep"]["violations"].values()) == 0,
}
record("a6_matches_independent_audit", all(xcheck.values()),
       {k: bool(v) for k, v in xcheck.items()} | {"auditor bistable cases": aud["A6_random_sweep"]["n_bistable"]})

# ---------------------------------------------------------------------------
# 6. full bookkeeping (EQUATIONS A1) with capped sequestration and native degradation:
# exact conservation, nonnegativity, M_free >= 0, and the native quota identity at steady state.
# one client, one chaperone (bind, release N/U/A, hand-off to D), first-order protease on U,
# native degradation N -> D, disaggregation A -> U, aggregate degradation A -> D.
# B = min(beta A, Mtot - C)  =>  M_free = max(Mtot - C - beta A, 0).
# a quadrature variable Q integrates s - mu (U + N + A + C), so U+N+A+C+D - Q must stay constant.
# ---------------------------------------------------------------------------
pf = dict(s=3.0, sig=0.3, kf=0.4, ku=0.05, kdN=0.08, ka=0.6, kdis=0.1, kdA=0.05, kp=0.2, mu=0.3,
          ron=2.0, kN=1.0, kU=0.5, kA=0.1, kD=0.2, Mtot=1.0, beta=0.5)


def full_rhs(t, x, p):
    U, N, A, C, D, Q = x
    Mf = max(p["Mtot"] - C - p["beta"] * A, 0.0)
    bind = p["ron"] * Mf * U
    dU = (1 - p["sig"]) * p["s"] + p["ku"] * N + p["kU"] * C + p["kdis"] * A \
        - (p["kf"] + p["ka"] + p["kp"] + p["mu"]) * U - bind
    dN = p["sig"] * p["s"] + p["kf"] * U + p["kN"] * C - (p["ku"] + p["kdN"] + p["mu"]) * N
    dA = p["ka"] * U + p["kA"] * C - (p["kdis"] + p["kdA"] + p["mu"]) * A
    dC = bind - (p["kN"] + p["kU"] + p["kA"] + p["kD"] + p["mu"]) * C
    dD = p["kD"] * C + p["kp"] * U + p["kdN"] * N + p["kdA"] * A
    dQ = p["s"] - p["mu"] * (U + N + A + C)
    return [dU, dN, dA, dC, dD, dQ]


solf = solve_ivp(full_rhs, (0, 300), [0.0] * 6, args=(pf,), rtol=1e-11, atol=1e-13, dense_output=True)
Uf, Nf, Af, Cf, Df, Qf = solf.y
inv = Uf + Nf + Af + Cf + Df - Qf
mfree_f = np.maximum(pf["Mtot"] - Cf - pf["beta"] * Af, 0.0)
uncapped_end = pf["Mtot"] - Cf[-1] - pf["beta"] * Af[-1]
pools_min = float(min(Uf.min(), Nf.min(), Af.min(), Cf.min(), Df.min()))
record("a1_conservation_and_positivity",
       float(np.max(np.abs(inv))) < 1e-7 and pools_min >= -1e-10 and float(Cf.max()) <= pf["Mtot"] + 1e-10
       and uncapped_end < 0,
       {"max |(U+N+A+C+D) - integral(s - mu pools)|": float(np.max(np.abs(inv))),
        "min over all pools and times": pools_min, "max C / Mtot": float(Cf.max() / pf["Mtot"]),
        "min free chaperone": float(mfree_f.min()),
        "uncapped Mtot - C - beta A at end (cap active if < 0)": float(uncapped_end)})

# positivity from a boundary start (all pools zero except a large aggregate, no synthesis)
solb = solve_ivp(full_rhs, (0, 100), [0.0, 0.0, 20.0, 0.0, 0.0, 0.0], args=(dict(pf, s=0.0),),
                 rtol=1e-11, atol=1e-13)
record("a1_positivity_from_boundary", float(solb.y[:5].min()) >= -1e-10,
       {"min over pools and times": float(solb.y[:5].min())})

# native quota at steady state: P = sig s + kf U + kN C equals (mu + ku + kdN) N
dx = np.max(np.abs(full_rhs(0, solf.y[:, -1], pf)[:4]))
P_nat = pf["sig"] * pf["s"] + pf["kf"] * Uf[-1] + pf["kN"] * Cf[-1]
quota_rhs = (pf["mu"] + pf["ku"] + pf["kdN"]) * Nf[-1]
record("a3_native_quota_includes_degradation",
       dx < 1e-8 and abs(P_nat - quota_rhs) < 1e-8 and abs(P_nat - (pf["mu"] + pf["ku"]) * Nf[-1]) > 1e-3,
       {"steady-state residual": float(dx), "productive flux P": float(P_nat),
        "(mu + ku + kdN) N": float(quota_rhs), "(mu + ku) N, old quota (wrong)": float((pf["mu"] + pf["ku"]) * Nf[-1])})

# ---------------------------------------------------------------------------
# 7. maintenance (EQUATIONS A7): moving quota face and growth-dependent feasibility
# ---------------------------------------------------------------------------
t7 = sp.symbols("t")
ku7, kdN7, P7 = sp.symbols("k_u k_dN P", positive=True)
N7, mu7 = sp.Function("N")(t7), sp.Function("mu")(t7)
Nreq7 = sp.Function("N_req")
phidot = sp.diff(N7 - Nreq7(mu7), t7).subs(sp.Derivative(N7, t7), P7 - (mu7 + ku7 + kdN7) * N7)
face = sp.simplify(phidot.subs(N7, Nreq7(mu7)))
expected = P7 - (mu7 + ku7 + kdN7) * Nreq7(mu7) - sp.diff(Nreq7(mu7), t7)
record("a7_face_condition_symbolic", sp.simplify(face - expected) == 0,
       {"d(N - N_req(mu))/dt on the face": str(face)})

# counterexample to the old A7.1, and sufficiency of the corrected one, in a 2-state toy:
# N' = P - (mu + ku) N, mu' = r (mu_t - mu), N_req = n0 + n1 mu, start on the face at mu0
ku_, n0_, n1_, mu0_, mut_, r_ = 0.1, 1.0, 2.0, 0.5, 1.0, 1.0
nreq = lambda m: n0_ + n1_ * m


def toy_gap(P):
    sol = solve_ivp(lambda t, x: [P - (x[1] + ku_) * x[0], r_ * (mut_ - x[1])], (0, 5.0),
                    [nreq(mu0_), mu0_], dense_output=True, rtol=1e-11, atol=1e-13)
    ts = np.linspace(0, 5.0, 2001)
    N_, m_t = sol.sol(ts)
    return float((N_ - nreq(m_t)).min())


P_old = (mu0_ + ku_) * nreq(mu0_)                        # old A7.1 holds with zero slack at t = 0
mus = np.linspace(mu0_, mut_, 2001)                      # corrected condition on every face point
P_new = float(np.max((mus + ku_) * nreq(mus) + n1_ * r_ * (mut_ - mus)))
record("a7_moving_quota_counterexample_and_fix",
       toy_gap(P_old) < -0.05 and toy_gap(P_new) >= -1e-9,
       {"old A7.1 slack at t=0": P_old - (mu0_ + ku_) * nreq(mu0_),
        "corrected face value at t=0 with P_old": P_old - (mu0_ + ku_) * nreq(mu0_) - n1_ * r_ * (mut_ - mu0_),
        "min gap N - N_req, P_old": toy_gap(P_old),
        "P meeting corrected condition on whole face": P_new, "min gap, P_new": toy_gap(P_new)})

# feasibility in mu with synthesis s = s_basal + s0 mu, productive fraction y, N_req = n0 + n1 mu,
# loss k = ku + kdN:  y (s_basal + s0 mu) >= (mu + k)(n0 + n1 mu)
y7, s07, n07, n17, k7, m7 = sp.symbols("y s_0 n_0 n_1 k mu", positive=True)
quad = sp.expand((m7 + k7) * (n07 + n17 * m7) - y7 * s07 * m7)
yc_sym = (sp.sqrt(n07) + sp.sqrt(n17 * k7)) ** 2 / s07
disc_zero = sp.simplify(sp.expand(sp.discriminant(quad, m7).subs(y7, yc_sym))) == 0
mu_star_sym = sp.solve(sp.diff(quad, m7), m7)[0].subs(y7, yc_sym)
mu_star_ok = sp.simplify(mu_star_sym - sp.sqrt(k7 * n07 / n17)) == 0


def feasible_interval(y, s0, n0, n1, k, s_basal=0.0):
    rts = np.roots([n1, n0 + n1 * k - y * s0, k * n0 - y * s_basal])
    rts = np.sort(rts[np.abs(rts.imag) < 1e-12].real)
    if len(rts) < 2 or rts[1] < 0:
        return None
    return (float(max(rts[0], 0.0)), float(rts[1]))


s0v, n0v, n1v = 10.0, 1.0, 2.0
iv = {}
for kdN_v in (0.0, 0.05):
    kv = 0.1 + kdN_v
    yc = (np.sqrt(n0v) + np.sqrt(n1v * kv)) ** 2 / s0v
    iv[f"kdN={kdN_v}"] = {"y_crit": float(yc), "mu_star": float(np.sqrt(kv * n0v / n1v)),
                          "interval 1.5 y_crit": feasible_interval(1.5 * yc, s0v, n0v, n1v, kv),
                          "interval 1.01 y_crit": feasible_interval(1.01 * yc, s0v, n0v, n1v, kv),
                          "interval 0.99 y_crit": feasible_interval(0.99 * yc, s0v, n0v, n1v, kv)}
ag = aud["A7_growth_interval"]
k0 = iv["kdN=0.0"]
# basal synthesis: mu -> 0 is viable iff y s_basal >= k n0
yb = 1.5 * (np.sqrt(n0v) + np.sqrt(n1v * 0.15)) ** 2 / s0v
sb_need = 0.15 * n0v / yb
lo_below = feasible_interval(yb, s0v, n0v, n1v, 0.15, s_basal=0.9 * sb_need)
lo_above = feasible_interval(yb, s0v, n0v, n1v, 0.15, s_basal=1.1 * sb_need)
record("a7_feasible_growth_interval",
       disc_zero and mu_star_ok
       and k0["interval 0.99 y_crit"] is None and k0["interval 1.01 y_crit"][0] > 0
       and abs(k0["y_crit"] - ag["y_crit"]) < 1e-12
       and np.allclose(k0["interval 1.5 y_crit"], ag["interval_at_1p5_ycrit"])
       and iv["kdN=0.05"]["mu_star"] > k0["mu_star"]
       and lo_below[0] > 0 and lo_above[0] == 0.0,
       {"y_crit = (sqrt n0 + sqrt(n1 k))^2/s0 symbolic": disc_zero, "mu* = sqrt(k n0/n1) symbolic": mu_star_ok,
        "intervals": iv, "matches auditor (kdN = 0)": True,
        "low end with s_basal = 0.9 / 1.1 x k n0 / y": [lo_below, lo_above]})

out = Path(__file__).with_name("checks.json")
out.write_text(json.dumps(results, indent=2, default=float))
bad = [k for k, v in results.items() if not v["pass"]]
for k, v in results.items():
    print(("PASS " if v["pass"] else "FAIL ") + k)
print(f"{len(results) - len(bad)}/{len(results)} pass")
sys.exit(1 if bad else 0)
