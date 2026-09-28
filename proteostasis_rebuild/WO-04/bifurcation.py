"""
WO-04 dynamics, stability, bifurcation.

three model variants, all conservative (WO-02 extended with one species, the
chaperone-aggregate complex CA, which holds one chaperone and one monomer-
equivalent of aggregate):

  V0  conservative core (WO-02). CA never forms (k_onA = 0).
  V1  sequestration: chaperone binds aggregate (C + A <-> CA) and is lost to
      folding; disaggregation stays a constant k_dis A.
  V2  chaperone-dependent disaggregation: CA -> C + U at k_dcat (DnaK/ClpB-like),
      and no chaperone-independent disaggregation (k_dis = 0).

full ODE (uM, s), states N U B A C CA:
  dN  = (1-eps)s_P + phi k_cat B - k_mis N - mu N
  dU  = eps s_P + k_mis N + ((1-phi)k_cat + k_off)B - k_on C U - k_d U - k_a U^2
        + k_dis A + k_dcat CA - mu U
  dB  = k_on C U - (k_off + k_cat + mu) B
  dA  = k_a U^2 - (k_dis + k_dA + mu) A - k_onA C A + k_offA CA
  dC  = s_C + (k_off + k_cat) B + (k_offA + k_dcat) CA - k_on C U - k_onA C A - mu C
  dCA = k_onA C A - (k_offA + k_dcat + mu) CA
totals: P_T = N+U+B+A+CA, C_T = C+B+CA.

steady-state reduction (exact at steady state, not a QSS approximation):
  C_T = s_C/mu, B = C U/K_M, CA = C A/K_Ae, C = C_T/(1 + U/K_M + A/K_Ae)
  K_M = (k_off+k_cat+mu)/k_on,  K_Ae = (k_offA+k_dcat+mu)/k_onA
  A(U) solves k_a U^2 = (k_dis+k_dA+mu) A + (k_dcat+mu) CA   (unique: RHS increasing in A)
  N = ((1-eps)s_P + phi k_cat B)/(k_mis+mu)
  G(U) = eps s_P + k_mis N - (phi k_cat + mu) B - (k_d+mu) U - k_a U^2 + k_dis A + k_dcat CA
steady states <-> roots of G on U >= 0.
"""
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.optimize import brentq, fsolve

sys.dont_write_bytecode = True

NAMES = ["N", "U", "B", "A", "C", "CA"]
PNAMES = ["s_P", "eps", "s_C", "mu", "k_on", "k_off", "k_cat", "phi", "k_d",
          "k_a", "k_dis", "k_dA", "k_mis", "k_onA", "k_offA", "k_dcat"]
XS = sp.symbols(NAMES, nonnegative=True)
PS = sp.symbols(PNAMES, nonnegative=True)
N, U, B, A, C, CA = XS
(s_P, eps, s_C, mu, k_on, k_off, k_cat, phi, k_d, k_a, k_dis, k_dA, k_mis,
 k_onA, k_offA, k_dcat) = PS

RHS = [
    (1 - eps) * s_P + phi * k_cat * B - k_mis * N - mu * N,
    eps * s_P + k_mis * N + ((1 - phi) * k_cat + k_off) * B - k_on * C * U - k_d * U
    - k_a * U**2 + k_dis * A + k_dcat * CA - mu * U,
    k_on * C * U - (k_off + k_cat + mu) * B,
    k_a * U**2 - (k_dis + k_dA + mu) * A - k_onA * C * A + k_offA * CA,
    s_C + (k_off + k_cat) * B + (k_offA + k_dcat) * CA - k_on * C * U - k_onA * C * A - mu * C,
    k_onA * C * A - (k_offA + k_dcat + mu) * CA,
]
_f = sp.lambdify((XS, PS), RHS, "numpy")
_J = sp.lambdify((XS, PS), sp.Matrix(RHS).jacobian(XS), "numpy")


def conservation_ok():
    dPT = RHS[0] + RHS[1] + RHS[2] + RHS[3] + RHS[5]
    dCT = RHS[4] + RHS[2] + RHS[5]
    okP = sp.simplify(dPT - (s_P - k_d * U - k_dA * A - mu * (N + U + B + A + CA))) == 0
    okC = sp.simplify(dCT - (s_C - mu * (C + B + CA))) == 0
    return okP, okC


def reduces_to_wo02():
    """with k_onA = k_dcat = 0 and CA = 0 the first five equations are WO-02's."""
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "WO-02"))
    import model as m2
    sub = {k_onA: 0, k_dcat: 0, CA: 0}
    mp = dict(zip(m2.X, XS[:5]))
    mp.update({getattr(m2, n): globals()[n] for n in m2.PARAM_NAMES})
    return all(sp.simplify(RHS[i].subs(sub) - m2.RHS_SYM[i].subs(mp)) == 0 for i in range(5))


def pv(d):
    return np.array([float(d[k]) for k in PNAMES])


def f(x, p):
    return np.asarray(_f(x, pv(p)), dtype=float)


def jac(x, p):
    return np.asarray(_J(x, pv(p)), dtype=float)


# ------------------------------------------------------------ reduction
def _K(p):
    KM = (p["k_off"] + p["k_cat"] + p["mu"]) / p["k_on"]
    KAe = (p["k_offA"] + p["k_dcat"] + p["mu"]) / p["k_onA"] if p["k_onA"] > 0 else np.inf
    return KM, KAe


def solve_A(Ug, p, iters=110):
    """vectorised bisection for A(U); monotone so the bracket always holds."""
    Ug = np.asarray(Ug, dtype=float)
    CT = p["s_C"] / p["mu"]
    KM, KAe = _K(p)
    lin = p["k_dis"] + p["k_dA"] + p["mu"]
    target = p["k_a"] * Ug**2
    lo = np.zeros_like(Ug)
    hi = target / lin
    if not np.isfinite(KAe):
        return hi
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        Cm = CT / (1 + Ug / KM + mid / KAe)
        val = lin * mid + (p["k_dcat"] + p["mu"]) * Cm * mid / KAe - target
        hi = np.where(val > 0, mid, hi)
        lo = np.where(val > 0, lo, mid)
    return 0.5 * (lo + hi)


def state_from_U(Ug, p):
    Ug = np.asarray(Ug, dtype=float)
    CT = p["s_C"] / p["mu"]
    KM, KAe = _K(p)
    Ag = solve_A(Ug, p)
    Cg = CT / (1 + Ug / KM + (Ag / KAe if np.isfinite(KAe) else 0.0))
    Bg = Cg * Ug / KM
    CAg = Cg * Ag / KAe if np.isfinite(KAe) else np.zeros_like(Ug)
    Ng = ((1 - p["eps"]) * p["s_P"] + p["phi"] * p["k_cat"] * Bg) / (p["k_mis"] + p["mu"])
    return Ng, Bg, Ag, Cg, CAg


def G(Ug, p):
    Ng, Bg, Ag, Cg, CAg = state_from_U(Ug, p)
    Ug = np.asarray(Ug, dtype=float)
    return (p["eps"] * p["s_P"] + p["k_mis"] * Ng - (p["phi"] * p["k_cat"] + p["mu"]) * Bg
            - (p["k_d"] + p["mu"]) * Ug - p["k_a"] * Ug**2 + p["k_dis"] * Ag + p["k_dcat"] * CAg)


def roots(p, n=6000):
    """all roots of G on (0, U_max], U_max = s_P/mu (no pool can exceed it)."""
    Umax = p["s_P"] / p["mu"]
    grid = np.geomspace(1e-12 * Umax, Umax, n)
    g = G(grid, p)
    idx = np.where(np.sign(g[:-1]) * np.sign(g[1:]) < 0)[0]
    out = []
    for i in idx:
        r = brentq(lambda u: float(G(np.array([u]), p)[0]), grid[i], grid[i + 1],
                   xtol=1e-15 * grid[i + 1], rtol=1e-14, maxiter=500)
        out.append(r)
    return out, (g[0] > 0)


def full_state(Ur, p):
    Ng, Bg, Ag, Cg, CAg = state_from_U(np.array([Ur]), p)
    return np.array([Ng[0], Ur, Bg[0], Ag[0], Cg[0], CAg[0]])


def classify(x, p):
    ev = np.linalg.eigvals(jac(x, p))
    return {"max_real_eig": float(np.max(ev.real)),
            "n_unstable": int(np.sum(ev.real > 0)),
            "stable": bool(np.all(ev.real < 0)),
            "det": float(np.linalg.det(jac(x, p))),
            "min_abs_real_eig": float(np.min(np.abs(ev.real)))}


# ------------------------------------------------------------ sampling domain
def sample(rng, variant):
    lu = lambda lo, hi: float(np.exp(rng.uniform(np.log(lo), np.log(hi))))
    p = {"eps": lu(1e-4, 0.5), "mu": lu(1e-5, 1e-3), "P_T": lu(300, 5000),
         "C_T": lu(1, 100), "k_on": lu(1e-2, 10), "k_off": lu(1e-2, 10),
         "k_cat": lu(1e-3, 1), "phi": rng.uniform(0.05, 1), "k_d": lu(1e-5, 1e-2),
         "k_a": lu(1e-6, 1e-1), "k_dis": lu(1e-5, 1e-2), "k_dA": lu(1e-7, 1e-3),
         "k_mis": lu(1e-7, 1e-4), "k_onA": 0.0, "k_offA": 1.0, "k_dcat": 0.0}
    if variant in ("V1", "V2"):
        p["k_onA"] = lu(1e-3, 10)
        p["k_offA"] = lu(1e-3, 10)
    if variant == "V2":
        p["k_dcat"] = lu(1e-4, 1)
        p["k_dis"] = 0.0
    p["s_P"] = p["mu"] * p["P_T"]
    p["s_C"] = p["mu"] * p["C_T"]
    return p


def fold_in_eps(p, eps_lo=1e-5, eps_hi=0.999, n=250):
    """scan eps; return eps values where the root count changes."""
    es = np.geomspace(eps_lo, eps_hi, n)
    counts = [len(roots({**p, "eps": e}, n=3000)[0]) for e in es]
    changes = [(es[i], es[i + 1], counts[i], counts[i + 1])
               for i in range(n - 1) if counts[i] != counts[i + 1]]
    return es, counts, changes


def locate_fold(p, e_lo, e_hi):
    """solve G = 0, dG/dU = 0 for (U, eps) starting from a bracketing interval
    in eps where the count changes. returns (U*, eps*)."""
    # bisection on eps using root count
    c_lo = len(roots({**p, "eps": e_lo}, n=3000)[0])
    for _ in range(60):
        e_mid = np.sqrt(e_lo * e_hi)
        c_mid = len(roots({**p, "eps": e_mid}, n=3000)[0])
        if c_mid == c_lo:
            e_lo = e_mid
        else:
            e_hi = e_mid
    # the two merging roots at the multi-root side
    side = e_lo if c_lo > len(roots({**p, "eps": e_hi}, n=3000)[0]) else e_hi
    rs = sorted(roots({**p, "eps": side}, n=3000)[0])
    gaps = np.diff(rs)
    j = int(np.argmin(gaps))
    U0 = 0.5 * (rs[j] + rs[j + 1])

    def eqs(z):
        u, le = z
        q = {**p, "eps": float(np.exp(le))}
        h = 1e-7 * u
        g = float(G(np.array([u]), q)[0])
        dg = float((G(np.array([u + h]), q)[0] - G(np.array([u - h]), q)[0]) / (2 * h))
        return [g / (p["s_P"]), dg / (p["mu"] + p["k_d"] + p["k_a"] * u + 1e-300)]
    sol = fsolve(eqs, [U0, np.log(side)], xtol=1e-13)
    return float(sol[0]), float(np.exp(sol[1])), (rs[j], rs[j + 1], side)


# ------------------------------------------------------------ G4.3 analytic (V0)
def analytic_v0():
    """analytical steady-state count for the conservative core V0 (k_onA = 0).

    step 1 (necessity + uniqueness of the reduction): each of dN, dB, dA, dCA
    and dC_T is linear in one unknown with a strictly positive coefficient
    (mu > 0), so sympy's solve returns exactly one expression for N, B, A, CA
    and C in terms of U. every steady state therefore lies on this curve.
    step 2 (sufficiency): on that curve dN = dB = dA = dC = dCA = 0
    identically and dU equals G(U).
    step 3 (count): dG/dU equals an explicit sum of terms that are each <= 0
    for nonnegative parameters, one of which, -(k_d + mu), is < 0 for mu > 0.
    G(0) = eps s_P + k_mis (1-eps) s_P/(k_mis+mu) > 0 whenever eps > 0, and
    G -> -inf. hence exactly one root on U >= 0: exactly one steady state."""
    CT, KM = sp.symbols("C_T K_M", positive=True)
    sub0 = {k_onA: 0}
    # step 1: solve each balance for its own unknown
    Bs = sp.solve(sp.Eq(RHS[2], 0), B)
    CAs = sp.solve(sp.Eq(RHS[5].subs(sub0), 0), CA)
    As = sp.solve(sp.Eq(RHS[3].subs(sub0).subs(CA, 0), 0), A)
    Ns = sp.solve(sp.Eq(RHS[0], 0), N)
    Ctot = sp.solve(sp.Eq(s_C - mu * CT, 0), CT)         # dC_T/dt = s_C - mu C_T
    Cs = sp.solve(sp.Eq(C + Bs[0] - CT, 0), C)          # C_T = C + B (CA = 0)
    single = all(len(s) == 1 for s in (Bs, CAs, As, Ns, Ctot, Cs))
    ca_zero = sp.simplify(CAs[0]) == 0
    km = (k_off + k_cat + mu) / k_on
    Cexpr = sp.simplify(Cs[0].subs(CT, s_C / mu))
    Bexpr = sp.simplify(Bs[0].subs(C, Cexpr))
    Aexpr = As[0]
    Nexpr = sp.simplify(Ns[0].subs(B, Bexpr))
    curve = {N: Nexpr, B: Bexpr, A: Aexpr, C: Cexpr, CA: 0}
    # step 2
    zero_others = all(sp.simplify(RHS[i].subs(sub0).subs(curve)) == 0 for i in (0, 2, 3, 4, 5))
    Gsym = sp.simplify(RHS[1].subs(sub0).subs(curve))
    # step 3: the claimed sign decomposition of dG/dU
    CT_ = s_C / mu
    claimed = (-(phi * k_cat * mu / (k_mis + mu) + mu) * CT_ * km / (km + U) ** 2
               - (k_d + mu)
               - 2 * k_a * U * (k_dA + mu) / (k_dis + k_dA + mu))
    dG_ok = sp.simplify(sp.diff(Gsym, U) - claimed) == 0
    G0 = sp.simplify(Gsym.subs(U, 0))
    G0_ok = sp.simplify(G0 - (eps * s_P + k_mis * (1 - eps) * s_P / (k_mis + mu))) == 0
    return {"reduction_unique": bool(single and ca_zero),
            "curve_satisfies_other_balances": bool(zero_others),
            "dG_dU_decomposition_holds": bool(dG_ok),
            "dG_dU": str(claimed), "G0": str(G0), "G0_formula_holds": bool(G0_ok),
            "conclusion": "exactly one steady state for all nonnegative parameters "
                          "with mu > 0 and eps > 0 (or k_mis > 0)"}


# ------------------------------------------------------------ G4.5 continuation
def eps_coeff(p):
    """G(U, eps) = eps * eps_coeff(p) + H(U); the U-dependent pools do not depend on eps."""
    return p["s_P"] * p["mu"] / (p["k_mis"] + p["mu"])


def eps_of_U(Ug, p):
    """exact equilibrium curve: every steady state is (U, eps(U))."""
    H = G(Ug, {**p, "eps": 0.0})
    return -H / eps_coeff(p)


def continue_branch(p, n=4000):
    """continuation of the full equilibrium curve parametrised by U (never
    turns in U). returns the physical part 0 < eps < 1, the sign of det(J) of
    the full 6-D system at each point, and the folds (extrema of eps(U))
    refined by brentq on d eps/dU."""
    Umax = p["s_P"] / p["mu"]
    Ug = np.geomspace(1e-9 * Umax, Umax, n)
    e = eps_of_U(Ug, p)
    keep = (e > 0) & (e < 1)
    Ug, e = Ug[keep], e[keep]
    dets = np.array([np.sign(np.linalg.det(jac(full_state(u, {**p, "eps": ee}), {**p, "eps": ee})))
                     for u, ee in zip(Ug, e)])

    def de(u):
        h = 1e-6 * u
        return float((eps_of_U(np.array([u + h]), p)[0] - eps_of_U(np.array([u - h]), p)[0]) / (2 * h))

    d = np.array([de(u) for u in Ug])
    folds = []
    for i in np.where(np.sign(d[:-1]) * np.sign(d[1:]) < 0)[0]:
        uf = brentq(de, Ug[i], Ug[i + 1], xtol=1e-14 * Ug[i + 1], rtol=1e-13)
        ef = float(eps_of_U(np.array([uf]), p)[0])
        q = {**p, "eps": ef}
        xf = full_state(uf, q)
        ev = np.linalg.eigvals(jac(xf, q))
        # reference scale: slowest eigenvalue at the regular curve points
        # either side of the fold (grid neighbours, not the fold itself)
        ref = min(classify(full_state(Ug[j], {**p, "eps": e[j]}), {**p, "eps": e[j]})["min_abs_real_eig"]
                  for j in (max(i - 5, 0), min(i + 6, len(Ug) - 1)))
        folds.append({"U_fold": float(uf), "eps_fold": ef,
                      "G_rel_at_fold": float(G(np.array([uf]), q)[0] / p["s_P"]),
                      "min_abs_eig_at_fold": float(np.min(np.abs(ev))),
                      "ref_min_abs_eig_nearby": float(ref),
                      "det_sign_left": float(dets[max(i - 5, 0)]),
                      "det_sign_right": float(dets[min(i + 6, len(Ug) - 1)]),
                      "kind": "max" if d[i] > 0 else "min"})
    det_changes = [(float(Ug[i]), float(Ug[i + 1])) for i in np.where(dets[:-1] * dets[1:] < 0)[0]]
    return {"n_points": int(len(Ug)), "folds": folds, "det_sign_changes": det_changes,
            "U": Ug, "eps": e, "det_sign": dets}
