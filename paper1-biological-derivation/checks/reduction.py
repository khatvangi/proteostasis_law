"""
numeric steady-state reduction for stages 2 and 4 (derivation in
02_FINITE_POOLS_AND_CAPACITY.md, 03_BASELINE_ANALYSIS.md, 04_FEEDBACK_A_...md).

at steady state every pool except U is fixed by its own balance:
  machine totals      C_T = s_C/mu, D_T = s_D/mu, E_T = s_E/mu   (mu > 0)
  chaperone split     C = C_T/w, B = C u, CA = C a,  w = 1 + u + a,
                      u = U/K_M, a = A/K_Ae  (a = 0 before stage 4)
  protease complex    DU = D_T U/(K_MD + U)
  disaggregase        EA = E_T A/(K_ME + A)
  aggregate           F_A(U, A) = k_a U^2 - lin A - (k_catE+mu) EA - mu CA = 0
  native              N = ((1-eps) s_P + phi k_cat B)/(k_mis + mu)
and the U balance on that manifold is the organising function

  G(U) = S_0 - kappa_B B - mu U - (k_dA+mu) A - mu CA
             - [(k_catD+mu) DU + mu EA]        (finite protease/disaggregase)
             - [k_d U]                          (first-order limit, 'C' stages)

  S_0 = s_P (k_mis + eps mu)/(k_mis + mu),  kappa_B = phi k_cat mu/(k_mis+mu) + mu

the slope G_U is returned term by term (T_B, T_DU, T_U, T_A, T_EA, T_CA) so a
fold can be attributed to the term that changes sign.

`finite=True` uses the protease and disaggregase cycles (stages S2, S4);
`finite=False` uses first-order k_d U and k_dis A (S2C = WO-02, S4C = WO-04 V1).
`seq=True` switches on feedback A (C + A <-> CA).
"""
import numpy as np


def consts(p):
    mu = p["mu"]
    K = {"K_M": (p["k_off"] + p["k_cat"] + mu) / p["k_on"]}
    K["K_Ae"] = (p["k_offA"] + mu) / p["k_onA"] if p.get("k_onA", 0) > 0 else np.inf
    if p.get("k_onD", 0) > 0:
        K["K_MD"] = (p["k_offD"] + p["k_catD"] + mu) / p["k_onD"]
    if p.get("k_onE", 0) > 0:
        K["K_ME"] = (p["k_offE"] + p["k_catE"] + mu) / p["k_onE"]
    k_mis = p["k_mis"]
    K["S_0"] = p["s_P"] * (k_mis + p["eps"] * mu) / (k_mis + mu)
    K["kappa_B"] = p["phi"] * p["k_cat"] * mu / (k_mis + mu) + mu
    K["G_eps"] = p["s_P"] * mu / (k_mis + mu)
    return K


def _pools(Ug, Ag, p, K, finite, seq):
    CT = p["C_T"]
    u = Ug / K["K_M"]
    a = Ag / K["K_Ae"] if seq else np.zeros_like(Ag)
    w = 1.0 + u + a
    out = {"C": CT / w, "B": CT * u / w, "CA": CT * a / w, "w": w, "u": u, "a": a}
    if finite:
        out["DU"] = p["D_T"] * Ug / (K["K_MD"] + Ug)
        out["EA"] = p["E_T"] * Ag / (K["K_ME"] + Ag)
    else:
        out["DU"] = np.zeros_like(Ug)
        out["EA"] = np.zeros_like(Ag)
    return out


def _lin(p, finite):
    return p["k_dA"] + p["mu"] + (0.0 if finite else p["k_dis"])


def F_A(Ug, Ag, p, K, finite, seq):
    q = _pools(Ug, Ag, p, K, finite, seq)
    return (p["k_a"] * Ug**2 - _lin(p, finite) * Ag
            - ((p["k_catE"] + p["mu"]) * q["EA"] if finite else 0.0) - p["mu"] * q["CA"])


def solve_A(Ug, p, K, finite, seq, iters=200):
    """A(U) from F_A = 0. F_A is strictly decreasing in A (proved symbolically in
    check_stage4_fold.py), so bisection on [0, k_a U^2 / lin] always brackets."""
    Ug = np.atleast_1d(np.asarray(Ug, dtype=float))
    if not finite and seq:
        return _A_quadratic(Ug, p, K)
    lo = np.zeros_like(Ug)
    hi = p["k_a"] * Ug**2 / _lin(p, finite)
    if not finite and not seq:
        return hi                                   # linear: exact
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        pos = F_A(Ug, mid, p, K, finite, seq) > 0
        lo = np.where(pos, mid, lo)
        hi = np.where(pos, hi, mid)
    return 0.5 * (lo + hi)


def _A_quadratic(Ug, p, K):
    """S4C (= WO-04 V1): multiplying F_A = 0 by (w0 + A/K_Ae), w0 = 1 + U/K_M,
    gives  (L/K) A^2 + (L w0 + mu C_T/K - k_a U^2/K) A - k_a U^2 w0 = 0
    with L = k_dis + k_dA + mu, K = K_Ae. the constant term is negative, so
    exactly one positive root. written cancellation-free."""
    L, Kae, CT = _lin(p, False), K["K_Ae"], p["C_T"]
    w0 = 1.0 + Ug / K["K_M"]
    qa = L / Kae
    qb = L * w0 + p["mu"] * CT / Kae - p["k_a"] * Ug**2 / Kae
    qc = -p["k_a"] * Ug**2 * w0
    disc = np.sqrt(qb * qb - 4 * qa * qc)
    with np.errstate(divide="ignore", invalid="ignore"):   # np.where evaluates both branches
        return np.where(qb > 0, -2 * qc / (qb + disc), (-qb + disc) / (2 * qa))


def G(Ug, p, finite, seq, eps=None):
    K = consts(p if eps is None else {**p, "eps": eps})
    Ug = np.atleast_1d(np.asarray(Ug, dtype=float))
    Ag = solve_A(Ug, p, K, finite, seq)
    q = _pools(Ug, Ag, p, K, finite, seq)
    mu = p["mu"]
    g = K["S_0"] - K["kappa_B"] * q["B"] - mu * Ug - (p["k_dA"] + mu) * Ag - mu * q["CA"]
    if finite:
        g = g - (p["k_catD"] + mu) * q["DU"] - mu * q["EA"]
    else:
        g = g - p["k_d"] * Ug
    return g


def slope_terms(Ug, p, finite, seq):
    """G_U = T_B + T_DU + T_U + T_A + T_EA + T_CA, each from the implicit
    function theorem on F_A (A' = -F_U/F_AA)."""
    K = consts(p)
    Ug = np.atleast_1d(np.asarray(Ug, dtype=float))
    Ag = solve_A(Ug, p, K, finite, seq)
    q = _pools(Ug, Ag, p, K, finite, seq)
    mu, CT, KM, w = p["mu"], p["C_T"], K["K_M"], q["w"]
    Kae = K["K_Ae"] if seq else np.inf
    B_U = CT * (1 + q["a"]) / (KM * w**2)
    B_A = -CT * q["u"] / (Kae * w**2) if seq else 0.0 * Ug
    CA_U = -CT * q["a"] / (KM * w**2) if seq else 0.0 * Ug
    CA_A = CT * (1 + q["u"]) / (Kae * w**2) if seq else 0.0 * Ug
    if finite:
        EA_A = p["E_T"] * K["K_ME"] / (K["K_ME"] + Ag)**2
        DU_U = p["D_T"] * K["K_MD"] / (K["K_MD"] + Ug)**2
    else:
        EA_A = 0.0 * Ug
        DU_U = 0.0 * Ug
    F_U = 2 * p["k_a"] * Ug - mu * CA_U
    F_AA = -_lin(p, finite) - ((p["k_catE"] + mu) * EA_A if finite else 0.0) - mu * CA_A
    Ap = -F_U / F_AA
    T = {"T_B": -K["kappa_B"] * (B_U + B_A * Ap),
         "T_DU": -(p["k_catD"] + mu) * DU_U if finite else 0.0 * Ug,
         "T_U": -(mu + (0.0 if finite else p["k_d"])) + 0.0 * Ug,
         "T_A": -(p["k_dA"] + mu) * Ap,
         "T_EA": -mu * EA_A * Ap,
         "T_CA": -mu * (CA_U + CA_A * Ap)}
    T["G_U"] = sum(T[k] for k in ("T_B", "T_DU", "T_U", "T_A", "T_EA", "T_CA"))
    T.update({"A": Ag, "Aprime": Ap, "dB": B_U + B_A * Ap, "B": q["B"], "C": q["C"],
              "CA": q["CA"], "a": q["a"], "u": q["u"]})
    return T


def full_state(Ur, p, stage, eps=None):
    """reconstruct every pool of a stage from a root U* (order model.stage_states)."""
    finite = stage in ("S2", "S4")
    seq = stage in ("S4", "S4C")
    pe = p if eps is None else {**p, "eps": eps}
    K = consts(pe)
    Ug = np.array([float(Ur)])
    Ag = solve_A(Ug, pe, K, finite, seq)
    q = _pools(Ug, Ag, pe, K, finite, seq)
    Nn = ((1 - pe["eps"]) * pe["s_P"] + pe["phi"] * pe["k_cat"] * q["B"]) / (pe["k_mis"] + pe["mu"])
    v = {"N": Nn[0], "U": Ur, "B": q["B"][0], "C": q["C"][0], "A": Ag[0], "CA": q["CA"][0]}
    if finite:
        v.update({"DU": q["DU"][0], "D": pe["D_T"] - q["DU"][0],
                  "EA": q["EA"][0], "E": pe["E_T"] - q["EA"][0]})
    from model import stage_states
    return np.array([v[n] for n in stage_states(stage)])


def roots(p, finite, seq, n=20000, eps=None):
    """sign changes of G on (0, U_max], U_max = s_P/mu (no pool exceeds the
    balanced client total), refined by brentq."""
    from scipy.optimize import brentq
    Umax = p["s_P"] / p["mu"]
    # lower end fixed in absolute terms as well, so that small mu (huge Umax)
    # cannot push the grid above a root at small U
    grid = np.geomspace(min(1e-14 * Umax, 1e-12), Umax, n)
    g = G(grid, p, finite, seq, eps)
    idx = np.where(np.sign(g[:-1]) * np.sign(g[1:]) < 0)[0]
    return [brentq(lambda x: float(G(x, p, finite, seq, eps)[0]), grid[i], grid[i + 1],
                   xtol=1e-15 * grid[i + 1], rtol=1e-15, maxiter=500) for i in idx]
