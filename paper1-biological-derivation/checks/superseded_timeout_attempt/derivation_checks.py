"""
paper 1 derivation build: machine checks for every transition claimed in
files 01-10. run from this directory:

    PYTHONDONTWRITEBYTECODE=1 python derivation_checks.py

writes derivation_checks.json and exits 1 if any check fails.

design
------
one flux table (name -> rate, donor, receivers) defines the most general model
of the build (file 08). every earlier stage is the same table with fluxes
switched off, so the object checked at stage k is literally a sub-table of the
object checked at stage k+1. equations are assembled from the table
(stoichiometry x flux - dilution), never typed by hand.

numerical checks use random parameters drawn log-uniformly over wide ranges.
they test mathematics only. no parameter value here is a claim about E. coli.
the source work orders (proteostasis_rebuild/WO-02, WO-04) are imported
read-only to confirm that the stages reproduce their equations.
"""
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.optimize import brentq, minimize_scalar

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
REBUILD = HERE.parents[1] / "proteostasis_rebuild"
RESULTS = {}


def record(key, ok, **info):
    RESULTS[key] = {"pass": bool(ok), **info}
    print(("PASS " if ok else "FAIL ") + key, info if info else "")


# ---------------------------------------------------------------- symbols
STATE_NAMES = ["N", "U", "B", "A", "C", "CA", "Z"]
X = sp.symbols(STATE_NAMES, nonnegative=True)
N, U, B, A, C, CA, Z = X
PNAMES = ["s_P", "eps", "s_C0", "s_C1", "K_sig", "eps_C", "mu", "k_on", "k_off",
          "k_cat", "phi", "k_d", "k_a", "k_dis", "k_dA", "k_mis", "k_onA",
          "k_offA", "k_dcat", "k_inact", "k_dZ"]
P = sp.symbols(PNAMES, nonnegative=True)
(s_P, eps, s_C0, s_C1, K_sig, eps_C, mu, k_on, k_off, k_cat, phi, k_d, k_a,
 k_dis, k_dA, k_mis, k_onA, k_offA, k_dcat, k_inact, k_dZ) = P

# chaperone synthesis flux. s_C1 = 0 is the constant-synthesis case; s_C1 > 0 is
# the sigma32-titration induction of file 07: free sigma32 = S_T K_sig/(K_sig + C)
s_C = s_C0 + s_C1 * K_sig / (K_sig + C)

# flux name -> (rate, {pool: stoichiometric coefficient}). empty dict = sink.
FLUX = {
    "syn_N":   ((1 - eps) * s_P,          {"N": 1}),
    "syn_U":   (eps * s_P,                {"U": 1}),
    "syn_C":   ((1 - eps_C) * s_C,        {"C": 1}),
    "syn_Cx":  (eps_C * s_C,              {"U": 1}),   # mistranslated chaperone chains
    "unfold":  (k_mis * N,                {"N": -1, "U": 1}),
    "degU":    (k_d * U,                  {"U": -1}),
    "agg":     (k_a * U**2,               {"U": -1, "A": 1}),
    "dis":     (k_dis * A,                {"A": -1, "U": 1}),
    "degA":    (k_dA * A,                 {"A": -1}),
    "bind":    (k_on * C * U,             {"C": -1, "U": -1, "B": 1}),
    "rel":     (k_off * B,                {"B": -1, "C": 1, "U": 1}),
    "cycN":    (phi * k_cat * B,          {"B": -1, "C": 1, "N": 1}),
    "cycU":    ((1 - phi) * k_cat * B,    {"B": -1, "C": 1, "U": 1}),
    "bindA":   (k_onA * C * A,            {"C": -1, "A": -1, "CA": 1}),
    "relA":    (k_offA * CA,              {"CA": -1, "C": 1, "A": 1}),
    "dcat":    (k_dcat * CA,              {"CA": -1, "C": 1, "U": 1}),
    "inact":   (k_inact * CA,             {"CA": -1, "A": 1, "Z": 1}),
    "degZ":    (k_dZ * Z,                 {"Z": -1}),
}

# stage -> fluxes present. dilution -mu*x acts on every state of the stage.
STAGE1 = ["syn_N", "syn_U", "unfold", "degU", "agg", "dis", "degA"]
STAGE2 = STAGE1 + ["syn_C", "bind", "rel", "cycN", "cycU"]
STAGES = {
    "S1": (STAGE1, ["N", "U", "A"]),
    "S2": (STAGE2, ["N", "U", "B", "A", "C"]),
    "FA": (STAGE2 + ["bindA", "relA"], ["N", "U", "B", "A", "C", "CA"]),
    "FB": (STAGE2 + ["bindA", "relA", "dcat"], ["N", "U", "B", "A", "C", "CA"]),
    "FC": (STAGE2 + ["bindA", "relA", "inact", "degZ"], ["N", "U", "B", "A", "C", "CA", "Z"]),
    "FCi": (STAGE2 + ["syn_Cx"], ["N", "U", "B", "A", "C"]),       # influx-mode control
    "F07": (STAGE2, ["N", "U", "B", "A", "C"]),                       # adaptive via s_C1
    "ALL": (list(FLUX), STATE_NAMES),
}
SYM = dict(zip(STATE_NAMES, X))

# parameter restrictions that define each stage (applied after assembly)
RESTRICT = {
    "S1": {s_C0: 0, s_C1: 0, eps_C: 0},
    "S2": {s_C1: 0, eps_C: 0},
    "FA": {s_C1: 0, eps_C: 0},
    "FB": {s_C1: 0, eps_C: 0},
    "FC": {s_C1: 0, eps_C: 0},
    "FCi": {s_C1: 0},
    "F07": {eps_C: 0},
    "ALL": {},
}


def assemble(stage):
    fluxes, states = STAGES[stage]
    rhs = {s: sp.Integer(0) for s in states}
    for f in fluxes:
        rate, st = FLUX[f]
        for s, c in st.items():
            if s in rhs:
                rhs[s] += c * rate
    for s in states:
        rhs[s] += -mu * SYM[s]
    sub = RESTRICT[stage]
    return [sp.expand(rhs[s].subs(sub)) for s in states], [SYM[s] for s in states]


# ------------------------------------------------ K1 donor/receiver audit
def k1_donor_receiver():
    """every internal flux debits exactly the pools it names as donors and
    credits the named receivers; sources have no donor pool inside the model
    (their donor is translation), sinks have no receiver. a flux whose rate is
    positive but debits nothing and is not a named source would be a phantom."""
    sources = {"syn_N", "syn_U", "syn_C", "syn_Cx"}
    sinks = {"degU", "degA", "degZ"}
    bad = []
    for f, (rate, st) in FLUX.items():
        donors = [s for s, c in st.items() if c < 0]
        recv = [s for s, c in st.items() if c > 0]
        if f in sources and donors:
            bad.append(f)
        if f in sinks and recv:
            bad.append(f)
        if f not in sources and not donors:
            bad.append(f + ":no-donor")
        if f not in sinks and not recv:
            bad.append(f + ":no-receiver")
        # the rate of an internal flux must vanish when its donor is empty
        for d in donors:
            if sp.simplify(rate.subs(SYM[d], 0)) != 0:
                bad.append(f + ":rate-nonzero-at-empty-donor-" + d)
    record("K1_donor_receiver", not bad, problems=bad)


# ------------------------------------------------ K2 conservation identities
def k2_conservation():
    ok = {}
    for stage in STAGES:
        rhs, xs = assemble(stage)
        d = dict(zip([str(x) for x in xs], rhs))
        g = lambda k: d.get(k, 0)
        prot = g("N") + g("U") + g("B") + g("A") + g("CA")
        PT = sum(SYM[s] for s in ["N", "U", "B", "A", "CA"] if s in d)
        sc = s_C.subs(RESTRICT[stage])
        ec = eps_C.subs(RESTRICT[stage]) if eps_C in RESTRICT[stage] else eps_C
        claimP = s_P + ec * sc - k_d * U - k_dA * (A if "A" in d else 0) - mu * PT
        okP = sp.simplify(prot - claimP) == 0
        if "C" in d:
            F = sum(SYM[s] for s in ["C", "B", "CA"] if s in d)
            chap = g("C") + g("B") + g("CA")
            claimF = (1 - ec) * sc - mu * F - (k_inact * CA if "Z" in d else 0)
            okF = sp.simplify(chap - claimF) == 0
        else:
            okF = True
        if "Z" in d:
            okZ = sp.simplify(g("Z") - (k_inact * CA - (k_dZ + mu) * Z)) == 0
        else:
            okZ = True
        ok[stage] = bool(okP and okF and okZ)
    record("K2_conservation_all_stages", all(ok.values()), per_stage=ok)


# ------------------------------------------------ K3 trace-back to WO-02 / WO-04
def k3_traceback():
    sys.path.insert(0, str(REBUILD / "WO-01"))
    sys.path.insert(0, str(REBUILD / "WO-02"))
    import model as m2
    rhs, xs = assemble("S2")
    mp = dict(zip(m2.X, xs))
    mp.update({getattr(m2, n): (s_C0 if n == "s_C" else globals()[n]) for n in m2.PARAM_NAMES})
    same2 = all(sp.simplify(rhs[i] - m2.RHS_SYM[i].subs(mp, simultaneous=True)) == 0 for i in range(5))
    sys.path.insert(0, str(REBUILD / "WO-04"))
    import bifurcation as b4
    rhsB, xsB = assemble("FB")
    # WO-04 order N U B A C CA equals this stage's order
    mp4 = dict(zip(b4.XS, xsB))
    mp4.update({getattr(b4, n): (s_C0 if n == "s_C" else globals()[n]) for n in b4.PNAMES})
    same4 = all(sp.simplify(rhsB[i] - b4.RHS[i].subs(mp4, simultaneous=True)) == 0 for i in range(6))
    record("K3_stage2_equals_WO02_and_FB_equals_WO04", same2 and same4,
           WO02_match=bool(same2), WO04_match=bool(same4))


# ------------------------------------------------ manifold coordinates
CT, Cf = sp.symbols("C_T C_f", nonnegative=True)


def manifold_jacobian(stage):
    """jacobian in coordinates where free chaperone is eliminated through the
    conserved functional total: C = C_T - B - CA (C_T a constant on the
    invariant manifold C_T = s_C/mu). returns (J, ys)."""
    rhs, xs = assemble(stage)
    names = [str(x) for x in xs]
    keep = [i for i, s in enumerate(names) if s not in ("C", "Z")]
    sub = {C: CT - B - (CA if "CA" in names else 0)}
    F = sp.Matrix([rhs[i].subs(sub) for i in keep])
    ys = [xs[i] for i in keep]
    return sp.expand(F.jacobian(ys)), ys


def nonneg_poly(expr, pos_sub):
    """true if expr, after substitution, is a polynomial in nonnegative symbols
    with nonnegative coefficients (a sufficient sign test)."""
    e = sp.expand(sp.together(expr.subs(pos_sub)))
    num, den = sp.fraction(e)
    num = sp.expand(num)
    if num == 0:
        return True
    poly = sp.Poly(num, *sorted(num.free_symbols, key=str))
    return all(c >= 0 for c in poly.coeffs()) and sp.Poly(den, *sorted(den.free_symbols, key=str) or [mu]).coeffs()[0] > 0


def k4_metzler_and_column_sums():
    """stage 2 (and influx self-damage FCi): off-diagonals are nonnegative on
    the physical domain C = C_T - B >= 0, and protein-weighted column sums are
    -(mu + degradation). the substitution C_T = B + C_f with C_f >= 0 encodes
    the domain."""
    out = {}
    for stage in ("S1", "S2", "FCi"):
        if stage == "S1":
            rhs, xs = assemble("S1")
            J = sp.Matrix(rhs).jacobian(xs)
            ys = xs
            dom = {}
        else:
            J, ys = manifold_jacobian(stage)
            dom = {CT: B + Cf}
        n = len(ys)
        off = [(str(ys[i]), str(ys[j]), J[i, j]) for i in range(n) for j in range(n) if i != j]
        neg = [(a, b, str(e)) for a, b, e in off if not nonneg_poly(e, dom)]
        colsum = [sp.simplify(sum(J[i, j] for i in range(n))) for j in range(n)]
        out[stage] = {"negative_offdiagonals": neg,
                      "column_sums": {str(ys[j]): str(colsum[j]) for j in range(n)}}
    expect = {"N": -mu, "U": -k_d - mu, "B": -mu, "A": -k_dA - mu}
    okcols = all(sp.simplify(sp.sympify(out[s]["column_sums"][k]) - v) == 0
                 for s in ("S2", "FCi") for k, v in expect.items())
    okmetz = all(not out[s]["negative_offdiagonals"] for s in out)
    record("K4_stage2_metzler_and_column_sums", okcols and okmetz, detail=out)


def k5_loss_of_metzler():
    """feedbacks A, B, C: list the off-diagonal entries that can be negative."""
    out = {}
    for stage in ("FA", "FB", "FC"):
        J, ys = manifold_jacobian(stage)
        dom = {CT: B + CA + Cf}
        n = len(ys)
        neg = [(str(ys[i]), str(ys[j]), str(J[i, j])) for i in range(n) for j in range(n)
               if i != j and not nonneg_poly(J[i, j], dom)]
        cols = {str(ys[j]): str(sp.simplify(sum(J[i, j] for i in range(n)))) for j in range(n)}
        out[stage] = {"negative_offdiagonals": neg, "column_sums": cols}
    expected = {("B", "CA"), ("CA", "B")}
    ok = all({(a, b) for a, b, _ in out[s]["negative_offdiagonals"]} == expected for s in out)
    record("K5_feedbacks_break_metzler_only_via_chaperone_competition", ok, detail=out)


# ------------------------------------------------ K6 stage-1 closed form
def k6_stage1():
    Lam = eps * s_P + k_mis * (1 - eps) * s_P / (k_mis + mu)
    kap = k_a * (k_dA + mu) / (k_dis + k_dA + mu)
    Us = 2 * Lam / ((k_d + mu) + sp.sqrt((k_d + mu)**2 + 4 * kap * Lam))
    rhs, xs = assemble("S1")
    Ns = (1 - eps) * s_P / (k_mis + mu)
    As = k_a * U**2 / (k_dis + k_dA + mu)
    fU = rhs[1].subs({N: Ns, A: As})
    ok_balance = sp.simplify(fU - (Lam - (k_d + mu) * U - kap * U**2)) == 0
    # the closed form solves kap U^2 + (k_d+mu) U - Lam = 0
    ok_root = sp.simplify(sp.expand(kap * Us**2 + (k_d + mu) * Us - Lam).rewrite(sp.sqrt)) == 0 or \
        abs(float((kap * Us**2 + (k_d + mu) * Us - Lam).subs(
            {eps: 0.03, s_P: 0.5, k_mis: 1e-5, mu: 2e-4, k_a: 1e-3, k_dA: 1e-5, k_dis: 4e-4, k_d: 3e-4}))) < 1e-15
    record("K6_stage1_quadratic_closed_form", ok_balance and ok_root)


# ------------------------------------------------ reduction for the general model
def lu(rng, lo, hi):
    return float(np.exp(rng.uniform(np.log(lo), np.log(hi))))


def sample(rng, stage="ALL"):
    p = {"eps": lu(rng, 1e-4, 0.5), "mu": lu(rng, 1e-5, 1e-3), "k_on": lu(rng, 1e-2, 10),
         "k_off": lu(rng, 1e-2, 10), "k_cat": lu(rng, 1e-3, 1), "phi": rng.uniform(0.05, 1),
         "k_d": lu(rng, 1e-5, 1e-2), "k_a": lu(rng, 1e-6, 1e-1), "k_dis": lu(rng, 1e-5, 1e-2),
         "k_dA": lu(rng, 1e-7, 1e-3), "k_mis": lu(rng, 1e-7, 1e-4), "k_onA": lu(rng, 1e-3, 10),
         "k_offA": lu(rng, 1e-3, 10), "k_dcat": lu(rng, 1e-4, 1), "k_inact": lu(rng, 1e-5, 1e-2),
         "k_dZ": lu(rng, 1e-5, 1e-3), "eps_C": lu(rng, 1e-4, 0.1), "K_sig": lu(rng, 0.1, 100)}
    P_T, C_T = lu(rng, 300, 5000), lu(rng, 1, 100)
    p["s_P"] = p["mu"] * P_T
    p["s_C0"] = p["mu"] * C_T
    p["s_C1"] = p["mu"] * lu(rng, 0.1, 100)
    for k, v in RESTRICT[stage].items():
        p[str(k)] = float(v)
    if stage == "FA":
        p["k_dcat"] = p["k_inact"] = 0.0
    if stage == "FB":
        p["k_inact"] = 0.0
    if stage == "FC":
        p["k_dcat"] = 0.0
    if stage in ("S2", "FCi", "F07"):
        p["k_onA"] = 0.0
    return p


def pv(p):
    return np.array([p[k] for k in PNAMES], dtype=float)


RHS_ALL, _ = assemble("ALL")
_f = sp.lambdify((X, P), RHS_ALL, "numpy")
_J = sp.lambdify((X, P), sp.Matrix(RHS_ALL).jacobian(X), "numpy")
_fp = {k: sp.lambdify((X, P), [sp.diff(r, sp.Symbol(k, nonnegative=True)) for r in RHS_ALL], "numpy")
       for k in ("eps", "s_C0")}


def sC(Cv, p):
    return p["s_C0"] + p["s_C1"] * p["K_sig"] / (p["K_sig"] + Cv)


def reduce_at(Uv, p):
    """steady-state values of every pool except U, given U (file 08 reduction).
    C solves a scalar equation whose left side increases and right side
    decreases in C, so the root is unique."""
    KM = (p["k_off"] + p["k_cat"] + p["mu"]) / p["k_on"]
    lamA = p["k_dis"] + p["k_dA"] + p["mu"]
    if p["k_onA"] > 0:
        KAe = (p["k_offA"] + p["k_dcat"] + p["k_inact"] + p["mu"]) / p["k_onA"]
    else:
        KAe = np.inf

    def A_of(Cv):
        return p["k_a"] * Uv**2 / (lamA + (p["k_dcat"] + p["mu"]) * Cv / KAe)

    def res(Cv):
        Av = A_of(Cv)
        return (p["mu"] * Cv * (1 + Uv / KM + Av / KAe) + p["k_inact"] * Cv * Av / KAe
                - (1 - p["eps_C"]) * sC(Cv, p))
    hi = 1.0
    while res(hi) < 0:
        hi *= 2
    Cv = brentq(res, 0.0, hi, xtol=1e-300, rtol=1e-15, maxiter=500)
    Av = A_of(Cv)
    Bv = Cv * Uv / KM
    CAv = Cv * Av / KAe
    Nv = ((1 - p["eps"]) * p["s_P"] + p["phi"] * p["k_cat"] * Bv) / (p["k_mis"] + p["mu"])
    Zv = p["k_inact"] * CAv / (p["k_dZ"] + p["mu"])
    return np.array([Nv, Uv, Bv, Av, Cv, CAv, Zv])


def Gfun(Uv, p):
    return float(_f(reduce_at(Uv, p), pv(p))[1])


def Lam(Uv, p):
    x = reduce_at(Uv, p)
    return (p["eps"] * p["s_P"] + p["k_mis"] * (1 - p["eps"]) * p["s_P"] / (p["k_mis"] + p["mu"])
            + p["eps_C"] * sC(x[4], p))


def Rfun(Uv, p):
    x = reduce_at(Uv, p)
    rho = p["phi"] * p["k_cat"] * p["mu"] / (p["k_mis"] + p["mu"]) + p["mu"]
    return rho * x[2] + (p["k_d"] + p["mu"]) * Uv + (p["k_dA"] + p["mu"]) * x[3] + p["mu"] * x[5]


def dfun(fn, Uv, p, h=1e-6):
    return (fn(Uv * (1 + h), p) - fn(Uv * (1 - h), p)) / (2 * h * Uv)


# ------------------------------------------------ K7 master identity G = Lambda - R
def k7_master_identity(n=300):
    rng = np.random.default_rng(7)
    worst = 0.0
    for _ in range(n):
        p = sample(rng, "ALL")
        Uv = lu(rng, 1e-3, 1e3)
        x = reduce_at(Uv, p)
        f = np.asarray(_f(x, pv(p)), dtype=float)
        scale = max(Lam(Uv, p), 1e-300)
        others = np.delete(f, 1)
        worst = max(worst, float(np.max(np.abs(others))) / scale,
                    abs(f[1] - (Lam(Uv, p) - Rfun(Uv, p))) / scale)
    record("K7_master_identity_G_equals_Lambda_minus_R", worst < 1e-9, max_rel_residual=worst, n=n)


# ------------------------------------------------ K8 Schur identity and sign of det J_yy
def k8_schur(n=300):
    rng = np.random.default_rng(8)
    worst, min_detyy_sign, stable_violations = 0.0, 1.0, 0
    for _ in range(n):
        p = sample(rng, "ALL")
        Uv = lu(rng, 1e-2, 1e2)
        x = reduce_at(Uv, p)
        J = np.asarray(_J(x, pv(p)), dtype=float)
        Jyy = np.delete(np.delete(J, 1, 0), 1, 1)
        dG = dfun(Gfun, Uv, p)
        lhs, rhs = np.linalg.det(J), np.linalg.det(Jyy) * dG
        worst = max(worst, abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1e-300))
        min_detyy_sign = min(min_detyy_sign, np.sign(np.linalg.det(Jyy)))
        # a point with dG/dU > 0 must have an eigenvalue with positive real part
        if dG > 0 and np.max(np.linalg.eigvals(J).real) <= 0:
            stable_violations += 1
    ok = worst < 1e-4 and min_detyy_sign > 0 and stable_violations == 0
    record("K8_detJ_equals_detJyy_times_dGdU", ok, max_rel_err=worst, n=n,
           det_Jyy_always_positive=bool(min_detyy_sign > 0),
           positive_slope_but_stable=stable_violations,
           note="finite-difference dG/dU (h=1e-6), so tolerance 1e-4")


# ------------------------------------------------ K9 contraction consequence (stage 2)
def k9_stage2_eigs(n=400):
    rng = np.random.default_rng(9)
    worst = -np.inf
    for _ in range(n):
        p = sample(rng, "S2")
        # unique root of G on (0, s_P/mu)
        Umax = p["s_P"] / p["mu"]
        Ur = brentq(lambda u: Gfun(u, p), 1e-14 * Umax, Umax, xtol=1e-300, rtol=1e-14)
        x = reduce_at(Ur, p)
        ev = np.linalg.eigvals(np.asarray(_J(x, pv(p)), dtype=float))
        worst = max(worst, float(np.max(ev.real) / p["mu"]))
    # contraction in the weighted l1 norm bounds every eigenvalue by -mu
    record("K9_stage2_all_eigenvalues_le_minus_mu", worst <= -1 + 1e-6,
           max_real_eig_over_mu=worst, n=n)


# ------------------------------------------------ K10 sequestration peak and the Gamma bound
def k10_sequestration_limit():
    u_, m_, al, KM_, KA_, lam_, ka_ = sp.symbols("u m alpha K_M K_A lambda_A k_a", positive=True)
    Uu = sp.Symbol("U", positive=True)
    # limit A = k_a U^2/lambda_A, B = C_T U/(K_M + U + K_M A/K_A)
    Bexpr = CT * Uu / (KM_ + Uu + KM_ * ka_ * Uu**2 / (lam_ * KA_))
    crit = sp.solve(sp.diff(Bexpr, Uu), Uu)
    Upeak = [c for c in crit if c.is_positive is not False]
    A_at_peak = sp.simplify(ka_ * Upeak[0]**2 / lam_)
    ok_peak = sp.simplify(A_at_peak - KA_) == 0
    # h(u; m) = m (1 - u^2)/(m + u + m u^2)^2 ; sup over m of -h at fixed u
    g = (u_**2 - 1) / (u_ * (1 + u_**2))
    ustar = sp.sqrt(2 + sp.sqrt(5))
    gmax = sp.nsimplify(g.subs(u_, ustar))
    Gamma_min = float(4 / g.subs(u_, ustar))
    m_star = float((ustar / (1 + ustar**2)).evalf())
    # numerical confirmation: minimise -h over (u, m) directly
    h = lambda uu, mm: mm * (1 - uu**2) / (mm + uu + mm * uu**2)**2
    grid_u = np.geomspace(1.0001, 100, 4000)
    best = max(float(np.max(-h(grid_u, mm))) for mm in np.geomspace(1e-3, 1e2, 2001))
    ok_bound = abs(1 / best - Gamma_min) / Gamma_min < 1e-3

    # Gamma_crit(m) at beta = 0, and the exact finite-beta test in the limit model
    def gamma_crit(mm):
        r = minimize_scalar(lambda uu: h(uu, mm), bounds=(1.0, 1e3), method="bounded",
                            options={"xatol": 1e-12})
        return 1.0 / (-r.fun)
    gc = {f"{mm:g}": gamma_crit(mm) for mm in (0.05, 0.1, 0.2, m_star, 0.5, 1.0, 2.0, 5.0)}
    record("K10_sequestration_peak_at_A_eq_KA_and_Gamma_bound", ok_peak and ok_bound,
           A_at_peak=str(A_at_peak), u_star=float(ustar.evalf()), g_max=float(gmax.evalf()),
           Gamma_min=Gamma_min, m_star=m_star, Gamma_crit_of_m=gc)
    return Gamma_min, m_star, gamma_crit


def k11_limit_vs_exact(Gamma_min, m_star, gamma_crit):
    """build exact feedback-A parameter sets deep in the limit
    mu C_T << K_A lambda_A, with beta ~ 0, and Gamma just above / below the
    critical value; count folds of the EXACT R(U) (sign changes of R')."""
    out = {}
    for tag, factor in (("above", 1.05), ("below", 0.95)):
        p = {k: 0.0 for k in PNAMES}
        p.update({"mu": 1e-6, "k_d": 1e-3, "k_dA": 0.0, "k_mis": 0.0, "phi": 1.0,
                  "k_on": 1.0, "k_off": 0.0, "k_onA": 1.0, "k_offA": 1.0, "k_dis": 1e-2,
                  "k_a": 1e-2, "eps": 0.01, "eps_C": 0.0})
        p["k_cat"] = 1.0
        KM = (p["k_off"] + p["k_cat"] + p["mu"]) / p["k_on"]
        KA = (p["k_offA"] + p["mu"]) / p["k_onA"]
        lamA = p["k_dis"] + p["k_dA"] + p["mu"]
        # choose k_a so that m = K_M/U_p = m_star
        Up = KM / m_star
        p["k_a"] = KA * lamA / Up**2
        rho = p["phi"] * p["k_cat"] * p["mu"] / (p["k_mis"] + p["mu"]) + p["mu"]
        CTv = factor * Gamma_min * (p["k_d"] + p["mu"]) * Up / rho
        p["s_C0"] = p["mu"] * CTv
        p["s_P"] = 1.0
        beta = p["k_a"] * (p["k_dA"] + p["mu"]) / lamA * Up / (p["k_d"] + p["mu"])
        Ug = np.geomspace(1e-3 * Up, 1e3 * Up, 6000)
        dR = np.array([dfun(Rfun, uu, p) for uu in Ug])
        out[tag] = {"Gamma_over_Gamma_min": factor, "beta": beta,
                    "limit_ratio_muCT_over_KA_lamA": p["mu"] * CTv / (KA * lamA),
                    "n_sign_changes_dR": int(np.sum(np.sign(dR[:-1]) != np.sign(dR[1:]))),
                    "min_dR_over_kd": float(np.min(dR) / (p["k_d"] + p["mu"]))}
    ok = out["above"]["n_sign_changes_dR"] == 2 and out["below"]["n_sign_changes_dR"] == 0
    record("K11_Gamma_threshold_reproduced_in_exact_model_within_limit", ok, detail=out)


def k12_wo04_crosswalk(n=400):
    """WO-04 V1 sampler (read-only import): among sets with 3 steady states,
    how many satisfy the limit-derived necessary condition Gamma > Gamma_min?
    the exact model is NOT in the limit everywhere, so violations are possible
    and are reported, not hidden."""
    sys.path.insert(0, str(REBUILD / "WO-04"))
    import bifurcation as b4
    rng = np.random.default_rng(104)
    Gmin = RESULTS["K10_sequestration_peak_at_A_eq_KA_and_Gamma_bound"]["Gamma_min"]
    multi, satisfied, limit_ratio = 0, 0, []
    gam_multi, gam_single = [], []
    for _ in range(n):
        q = b4.sample(rng, "V1")
        rs, _ = b4.roots(q, n=3000)
        KA = (q["k_offA"] + q["mu"]) / q["k_onA"]
        lamA = q["k_dis"] + q["k_dA"] + q["mu"]
        Up = np.sqrt(KA * lamA / q["k_a"])
        rho = q["phi"] * q["k_cat"] * q["mu"] / (q["k_mis"] + q["mu"]) + q["mu"]
        Gam = rho * (q["s_C"] / q["mu"]) / ((q["k_d"] + q["mu"]) * Up)
        if len(rs) == 3:
            multi += 1
            gam_multi.append(Gam)
            satisfied += Gam > Gmin
            limit_ratio.append(q["s_C"] / (KA * lamA))
        else:
            gam_single.append(Gam)
    record("K12_WO04_V1_multistable_sets_meet_Gamma_condition", satisfied == multi,
           n=n, n_multistable=multi, n_satisfying=int(satisfied),
           min_Gamma_multistable=float(min(gam_multi)) if gam_multi else None,
           frac_single_with_Gamma_above_min=float(np.mean(np.array(gam_single) > Gmin)),
           max_limit_ratio_in_multistable=float(max(limit_ratio)) if limit_ratio else None,
           note="necessary, not sufficient: many single-root sets also exceed Gamma_min")


# ------------------------------------------------ K13 self-damage: exact K_seq substitution
def k13_selfdamage_Kseq(n=200):
    rng = np.random.default_rng(13)
    worst = 0.0
    for _ in range(n):
        p = sample(rng, "FC")
        Uv = lu(rng, 1e-2, 1e2)
        x = reduce_at(Uv, p)
        KM = (p["k_off"] + p["k_cat"] + p["mu"]) / p["k_on"]
        KAe = (p["k_offA"] + p["k_inact"] + p["mu"]) / p["k_onA"]
        Kseq = KAe / (1 + p["k_inact"] / p["mu"])
        C_formula = (p["s_C0"] / p["mu"]) / (1 + Uv / KM + x[3] / Kseq)
        worst = max(worst, abs(C_formula - x[4]) / x[4])
    record("K13_selfdamage_equals_sequestration_with_Kseq", worst < 1e-10, max_rel_err=worst)


# ------------------------------------------------ K14 feedback B: quadratic for A(U)
def k14_feedbackB_quadratic(n=200):
    rng = np.random.default_rng(14)
    worst = 0.0
    for _ in range(n):
        p = sample(rng, "FB")
        Uv = lu(rng, 1e-2, 1e2)
        x = reduce_at(Uv, p)
        KM = (p["k_off"] + p["k_cat"] + p["mu"]) / p["k_on"]
        KAe = (p["k_offA"] + p["k_dcat"] + p["mu"]) / p["k_onA"]
        lamA = p["k_dis"] + p["k_dA"] + p["mu"]
        CT0 = p["s_C0"] / p["mu"]
        pi_ = (p["k_dcat"] + p["mu"]) / (p["k_offA"] + p["k_dcat"] + p["mu"])
        a2 = lamA / KAe
        a1 = lamA * (1 + Uv / KM) + p["k_onA"] * pi_ * CT0 - p["k_a"] * Uv**2 / KAe
        a0 = -p["k_a"] * Uv**2 * (1 + Uv / KM)
        Aq = (-a1 + np.sqrt(a1 * a1 - 4 * a2 * a0)) / (2 * a2)
        worst = max(worst, abs(Aq - x[3]) / x[3])
    record("K14_feedbackB_aggregate_quadratic", worst < 1e-8, max_rel_err=worst)


# ------------------------------------------------ K15 adaptive: C_T(U) increasing, R' > 0
def k15_adaptive_monotone(n=200):
    rng = np.random.default_rng(15)
    minslope, min_dCT = np.inf, np.inf
    for _ in range(n):
        p = sample(rng, "F07")
        Ug = np.geomspace(1e-3, 1e3, 200)
        R = np.array([Rfun(u, p) for u in Ug])
        CTs = np.array([reduce_at(u, p)[2] + reduce_at(u, p)[4] for u in Ug])
        minslope = min(minslope, float(np.min(np.diff(R))))
        min_dCT = min(min_dCT, float(np.min(np.diff(CTs))))
    record("K15_adaptive_alone_R_increasing_and_CT_increasing", minslope > 0 and min_dCT > 0,
           min_dR=minslope, min_dCT=min_dCT)


# ------------------------------------------------ K16 linear response = scalar formula
def k16_linear_response(n=200):
    rng = np.random.default_rng(16)
    worst, used = 0.0, 0
    for _ in range(n):
        p = sample(rng, "ALL")
        Umax = (p["s_P"] + p["eps_C"] * sC(0, p)) / p["mu"]
        grid = np.geomspace(1e-10 * Umax, Umax, 400)
        g = np.array([Gfun(u, p) for u in grid])
        idx = np.where(np.sign(g[:-1]) != np.sign(g[1:]))[0]
        if len(idx) == 0:
            continue
        Ur = brentq(lambda u: Gfun(u, p), grid[idx[0]], grid[idx[0] + 1], xtol=1e-300, rtol=1e-14)
        x = reduce_at(Ur, p)
        J = np.asarray(_J(x, pv(p)), dtype=float)
        slope = dfun(Rfun, Ur, p) - dfun(Lam, Ur, p)
        for c in ("eps", "s_C0"):
            fc = np.asarray(_fp[c](x, pv(p)), dtype=float)
            full = -np.linalg.solve(J, fc)[1]
            # scalar: (dLam/dc - dR/dc at fixed U) / (R' - Lam')
            h = 1e-6 * p[c]
            q1, q0 = {**p, c: p[c] + h}, {**p, c: p[c] - h}
            num = ((Lam(Ur, q1) - Rfun(Ur, q1)) - (Lam(Ur, q0) - Rfun(Ur, q0))) / (2 * h)
            scal = num / slope
            worst = max(worst, abs(full - scal) / max(abs(full), 1e-300))
        used += 1
    record("K16_minus_Jinv_fc_equals_scalar_response", worst < 1e-4, max_rel_err=worst, n_used=used)


# ------------------------------------------------ K17 flux-level additivity
def k17_additivity():
    fA, xA = assemble("FA")
    rhs_all, xs = assemble("ALL")
    # increments relative to FA, computed as the flux terms themselves
    dB = {s: 0 for s in STATE_NAMES}
    dC = {s: 0 for s in STATE_NAMES}
    for s, c in FLUX["dcat"][1].items():
        dB[s] += c * FLUX["dcat"][0]
    for f in ("inact", "degZ"):
        for s, c in FLUX[f][1].items():
            dC[s] += c * FLUX[f][0]
    fAd = dict(zip([str(x) for x in xA], fA))
    ok_ABC = all(sp.simplify(
        (rhs_all[i] - (fAd.get(s, -mu * SYM[s]) + dB[s] + dC[s])).subs({s_C1: 0, eps_C: 0})) == 0
        for i, s in enumerate(STATE_NAMES))
    # adaptive x influx: the chaperone-synthesis terms carry a product eps_C*s_C1
    dCeq = rhs_all[STATE_NAMES.index("C")]
    inter = sp.simplify(sp.diff(dCeq, eps_C, s_C1))
    ok_inter = inter != 0
    record("K17_ABC_flux_additive_but_adaptive_x_influx_not", ok_ABC and ok_inter,
           cross_term_in_dC=str(inter))


# ------------------------------------------------ K18 stage 2 knee: susceptibility ratio
def k18_knee():
    u_, G0, b0 = sp.symbols("u Gamma_0 beta_0", positive=True)
    r = G0 * u_ / (1 + u_) + u_ + b0 * u_**2
    slope0 = sp.diff(r, u_).subs(u_, 0)
    ok = sp.simplify(slope0 - (G0 + 1)) == 0
    # rescue share of total exit, and its limit at large u
    share = sp.limit((G0 * u_ / (1 + u_)) / r, u_, sp.oo)
    record("K18_stage2_dimensionless_exit_curve", ok and share == 0,
           slope_at_zero=str(slope0), rescue_share_large_u=str(share))


if __name__ == "__main__":
    k1_donor_receiver()
    k2_conservation()
    k3_traceback()
    k4_metzler_and_column_sums()
    k5_loss_of_metzler()
    k6_stage1()
    k7_master_identity()
    k8_schur()
    k9_stage2_eigs()
    Gmin, mstar, gc = k10_sequestration_limit()
    k11_limit_vs_exact(Gmin, mstar, gc)
    k12_wo04_crosswalk()
    k13_selfdamage_Kseq()
    k14_feedbackB_quadratic()
    k15_adaptive_monotone()
    k16_linear_response()
    k17_additivity()
    k18_knee()
    (HERE / "derivation_checks.json").write_text(json.dumps(RESULTS, indent=1, default=str))
    n_fail = sum(not v["pass"] for v in RESULTS.values())
    print(f"\n{len(RESULTS) - n_fail}/{len(RESULTS)} checks pass")
    sys.exit(1 if n_fail else 0)
