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
# 5. sequestration feedback: saddle-node location; no threshold without feedback
# dA/dt = ka L / (a - b A) - (mu + delta) A,  a = kc Mtot + ka + mu,  b = kc*beta
# ---------------------------------------------------------------------------
Aa, L, a, b, ka2, m = sp.symbols("A L a b k_a m", positive=True)  # m = mu + delta
f = ka2 * L / (a - b * Aa) - m * Aa
sn = sp.solve([f, sp.diff(f, Aa)], [Aa, L], dict=True)
sn_ok = any(sp.simplify(d[Aa] - a / (2 * b)) == 0 and sp.simplify(d[L] - m * a**2 / (4 * b * ka2)) == 0 for d in sn)
# b = 0: unique, linear steady state
uniq = sp.solve(ka2 * L / a - m * Aa, Aa)
record("sequestration_saddle_node", sn_ok and len(uniq) == 1,
       {"saddle-node": [{k.name: str(v) for k, v in d.items()} for d in sn], "b=0 steady state": str(uniq)})

# numeric: number of admissible steady states below / above L_crit
def n_roots(Lv, av=3.0, bv=1.0, kav=1.0, mv=1.0):
    # (mu+delta) b A^2 - (mu+delta) a A + ka L = 0 on 0 <= A < a/b
    r = np.roots([mv * bv, -mv * av, kav * Lv])
    return sum(1 for x in r if abs(x.imag) < 1e-12 and 0 <= x.real < av / bv)

Lc = 1.0 * 3.0**2 / (4 * 1.0 * 1.0)
record("bistable_below_none_above", n_roots(0.9 * Lc) == 2 and n_roots(1.1 * Lc) == 0,
       {"L_crit": Lc, "roots at 0.9 L_crit": n_roots(0.9 * Lc), "roots at 1.1 L_crit": n_roots(1.1 * Lc)})

out = Path(__file__).with_name("checks.json")
out.write_text(json.dumps(results, indent=2, default=float))
bad = [k for k, v in results.items() if not v["pass"]]
for k, v in results.items():
    print(("PASS " if v["pass"] else "FAIL ") + k)
sys.exit(1 if bad else 0)
