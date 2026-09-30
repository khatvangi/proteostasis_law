"""
permanent checks: stages 1-3 (baseline steady-state reduction).

  B1  stage 1 closed form: on the steady-state manifold dU/dt = S_0 - alpha U
      - beta U^2 exactly; one positive root; dU*/deps formula
  B2  stage 2 organising identity: (dU/dt + dA/dt) on the manifold where
      N, B, C, DU, D, EA, E are at their own steady states equals
        G = S_0 - kappa_B B - (k_catD+mu) DU - mu U - (k_dA+mu) A - mu EA
      and equals the client total balance dP_T/dt on that manifold
  B3  slope: dG/dU = sum of terms each <= 0, with -mu < 0 (symbolic, via the
      implicit-function derivative A' = -F_U/F_AA and sign-definite numerators)
  B4  nesting stage 2 -> stage 1: exact in the limit k_cat, k_catD, k_catE -> inf
      (machine-bound client fraction -> 0); the low-load (U << K_M) limit leaves
      a residual bound-client dilution term, which is reported, not hidden
  B5  reproduction of WO-04 V0: this code's G equals WO-04's G on sampled sets,
      the WO-04 V0 slope formula holds, and the seed-100 scan gives 1 root in
      1000/1000 sets
  B6  stage 2 numerics: one root per sampled set, full 9-pool residual ~ 0,
      all Jacobian eigenvalues in the left half-plane, analytic slope equals
      finite differences, every slope term <= 0 on a grid; long integrations
      land on the reduced root
  B7  mu = 0 against mu > 0: at mu = 0 every chaperone parameter drops out of
      the U balance; with k_dA = 0 a finite ceiling s_crit exists, the unique
      root escapes to infinity as s_P -> s_crit with G_U < 0 throughout (no
      saddle-node); with mu > 0 a unique root exists for every s_P
"""
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import model as m  # noqa: E402
import reduction as rd  # noqa: E402
from check_conservation import load_source  # noqa: E402

CT, DT, ET = sp.symbols("C_T D_T E_T", positive=True)
Kc = (m.k_off + m.k_cat + m.mu) / m.k_on
KD = (m.k_offD + m.k_catD + m.mu) / m.k_onD
KE = (m.k_offE + m.k_catE + m.mu) / m.k_onE
S0 = m.s_P * (m.k_mis + m.eps * m.mu) / (m.k_mis + m.mu)
kB = m.phi * m.k_cat * m.mu / (m.k_mis + m.mu) + m.mu


def all_nonneg_coeffs(expr):
    """true if expr is a ratio of polynomials in nonnegative symbols whose
    numerator and denominator have only nonnegative coefficients."""
    e = sp.together(sp.expand(expr))
    num, den = sp.fraction(e)
    for part in (sp.expand(num), sp.expand(den)):
        if part == 0:
            continue
        if part.is_number:
            if part < 0:
                return False
            continue
        if any(c < 0 for c in sp.Poly(part, *part.free_symbols).coeffs()):
            return False
    return True


def manifold_S2():
    """pools of stage 2 at their own steady states, as functions of (U, A)."""
    B_ = CT * m.U / (Kc + m.U)
    DU_ = DT * m.U / (KD + m.U)
    EA_ = ET * m.A / (KE + m.A)
    N_ = ((1 - m.eps) * m.s_P + m.phi * m.k_cat * B_) / (m.k_mis + m.mu)
    return {m.B: B_, m.C: CT - B_, m.DU: DU_, m.D: DT - DU_, m.EA: EA_, m.E: ET - EA_,
            m.N: N_, m.s_C: m.mu * CT, m.s_D: m.mu * DT, m.s_E: m.mu * ET}


def b1_stage1():
    rhs = m.assemble("S1")
    Ns = sp.solve(rhs["N"], m.N)
    As = sp.solve(rhs["A"], m.A)
    single = len(Ns) == 1 and len(As) == 1
    G1 = sp.simplify(rhs["U"].subs({m.N: Ns[0], m.A: As[0]}))
    alpha = m.k_r * m.mu / (m.k_mis + m.mu) + m.k_d + m.mu
    beta = m.k_a * (m.k_dA + m.mu) / (m.k_dis + m.k_dA + m.mu)
    form_ok = sp.simplify(G1 - (S0 - alpha * m.U - beta * m.U**2)) == 0
    # product of the two roots of beta U^2 + alpha U - S_0 is -S_0/beta < 0:
    # exactly one positive root whenever S_0 > 0 and beta > 0. verified on
    # abstract symbols (alpha, beta, S_0 do not depend on U), which avoids
    # simplifying a nested square root of the full parameter expression
    a_, b_, s_ = sp.symbols("alpha beta S_0", positive=True)
    Ustar = 2 * s_ / (a_ + sp.sqrt(a_**2 + 4 * b_ * s_))
    root_ok = sp.simplify(sp.radsimp(b_ * Ustar**2 + a_ * Ustar - s_)) == 0
    qa, qb, qc = sp.Poly(b_ * m.U**2 + a_ * m.U - s_, m.U).all_coeffs()
    vieta_ok = sp.simplify(qc / qa + s_ / b_) == 0      # root product = -S_0/beta < 0
    # implicit differentiation: only S_0 depends on eps, dS_0/deps = s_P mu/(k_mis+mu)
    dS0 = sp.simplify(sp.diff(S0, m.eps) - m.s_P * m.mu / (m.k_mis + m.mu)) == 0
    eps_free = not ({m.eps} & (alpha.free_symbols | beta.free_symbols))
    # total client on the manifold: s_P - k_d U - k_dA A - mu P_T equals G1, so
    # at the steady state P_T = (s_P - degradation)/mu, not s_P/mu
    PT = Ns[0] + m.U + As[0]
    PT_ok = sp.simplify(m.s_P - m.k_d * m.U - m.k_dA * As[0] - m.mu * PT - G1) == 0
    return {"pass": bool(single and form_ok and root_ok and vieta_ok and dS0 and eps_free and PT_ok),
            "reduction_unique": single, "G1_equals_S0_minus_alphaU_minus_betaU2": form_ok,
            "closed_form_root_verified": root_ok, "vieta_product_negative": vieta_ok,
            "dUstar_deps = s_P mu/((k_mis+mu)(alpha+2 beta U*))": bool(dS0 and eps_free),
            "client_balance_on_manifold_equals_G1": PT_ok,
            "alpha": str(alpha), "beta": str(beta)}


def b2_identity_S2():
    rhs = m.assemble("S2")
    M = manifold_S2()
    zero_ok = all(sp.simplify(rhs[n].subs(M)) == 0 for n in ("N", "B", "C", "DU", "D", "EA", "E"))
    GU = rhs["U"].subs(M) + rhs["A"].subs(M)
    claim = S0 - kB * M[m.B] - (m.k_catD + m.mu) * M[m.DU] - m.mu * m.U - (m.k_dA + m.mu) * m.A - m.mu * M[m.EA]
    identity = sp.simplify(GU - claim) == 0
    PTdot = (m.expected_total_rate("client", "S2")).subs(M)
    massbal = sp.simplify(GU - PTdot) == 0
    # the U balance alone differs from G by exactly the aggregate balance F_A
    FA = rhs["A"].subs(M)
    return {"pass": bool(zero_ok and identity and massbal), "other_balances_zero_on_manifold": zero_ok,
            "G_identity_residual_zero": identity, "G_equals_client_balance_on_manifold": massbal,
            "F_A": str(sp.simplify(FA)), "G": str(claim)}


def b3_slope_S2():
    M = manifold_S2()
    rhs = m.assemble("S2")
    FA = sp.simplify(rhs["A"].subs(M))
    GA = S0 - kB * M[m.B] - (m.k_catD + m.mu) * M[m.DU] - m.mu * m.U - (m.k_dA + m.mu) * m.A - m.mu * M[m.EA]
    FU, FAA = sp.diff(FA, m.U), sp.diff(FA, m.A)
    Ap = -FU / FAA
    total = sp.diff(GA, m.U) + sp.diff(GA, m.A) * Ap
    terms = {"T_B": -kB * sp.diff(M[m.B], m.U), "T_DU": -(m.k_catD + m.mu) * sp.diff(M[m.DU], m.U),
             "T_U": -m.mu, "T_A": -(m.k_dA + m.mu) * Ap, "T_EA": -m.mu * sp.diff(M[m.EA], m.A) * Ap}
    decomposition = sp.simplify(total - sum(terms.values())) == 0
    Ap_pos = all_nonneg_coeffs(FU) and all_nonneg_coeffs(-FAA)
    # T_A and T_EA contain A' = F_U/(-F_AA); their sign follows from the factor
    # signs (a ratio of two sign-definite factors), tested factor by factor
    signs = {"T_B": all_nonneg_coeffs(-terms["T_B"]), "T_DU": all_nonneg_coeffs(-terms["T_DU"]),
             "T_U": True,
             "T_A": all_nonneg_coeffs(m.k_dA + m.mu) and Ap_pos,
             "T_EA": all_nonneg_coeffs(sp.diff(M[m.EA], m.A)) and Ap_pos}
    # A(U) is unique: F_A strictly decreasing in A
    return {"pass": bool(decomposition and all(signs.values()) and Ap_pos),
            "decomposition_holds": decomposition, "each_term_nonpositive": signs,
            "F_U>=0_and_F_AA<0": Ap_pos,
            "consequence": "G_U <= -mu < 0: exactly one steady state for mu > 0; dU*/deps <= s_P/(k_mis+mu)"}


def b4_limits():
    rhs = m.assemble("S2")
    M = manifold_S2()
    G2 = S0 - kB * M[m.B] - (m.k_catD + m.mu) * M[m.DU] - m.mu * m.U - (m.k_dA + m.mu) * m.A - m.mu * M[m.EA]
    FA2 = rhs["A"].subs(M)
    alpha_r = {m.k_r: m.phi * m.k_on * CT, m.k_d: m.k_onD * DT, m.k_dis: m.k_onE * ET}
    G1 = (S0 - (m.k_r * m.mu / (m.k_mis + m.mu) + m.k_d + m.mu) * m.U - (m.k_dA + m.mu) * m.A).subs(alpha_r)
    FA1 = (m.k_a * m.U**2 - (m.k_dis + m.k_dA + m.mu) * m.A).subs(alpha_r)
    lamG, lamF = G2, FA2
    for kc in (m.k_cat, m.k_catD, m.k_catE):
        lamG = sp.limit(lamG, kc, sp.oo)
        lamF = sp.limit(lamF, kc, sp.oo)
    exact_G = sp.simplify(lamG - G1) == 0
    exact_F = sp.simplify(lamF - FA1) == 0
    # low-load linearisation (architecture wording): first order in U and A
    KMs, KMDs, KMEs = sp.symbols("K_M K_MD K_ME", positive=True)
    G2k = G2.subs({Kc: KMs, KD: KMDs, KE: KMEs})
    G2k = S0 - kB * CT * m.U / (KMs + m.U) - (m.k_catD + m.mu) * DT * m.U / (KMDs + m.U) \
        - m.mu * m.U - (m.k_dA + m.mu) * m.A - m.mu * ET * m.A / (KMEs + m.A)
    t = sp.Symbol("t", positive=True)
    lin2 = sp.series(G2k.subs({m.U: t * m.U, m.A: t * m.A}), t, 0, 2).removeO().subs(t, 1)
    G1lin = (S0 - (m.k_r * m.mu / (m.k_mis + m.mu) + m.k_d + m.mu) * m.U - (m.k_dA + m.mu) * m.A).subs(
        {m.k_r: m.phi * m.k_cat * CT / KMs, m.k_d: m.k_catD * DT / KMDs})
    resid = sp.simplify(sp.expand(lin2 - G1lin))
    expected_resid = -m.mu * (CT / KMs + DT / KMDs) * m.U - m.mu * ET / KMEs * m.A
    resid_ok = sp.simplify(resid - expected_resid) == 0
    return {"pass": bool(exact_G and exact_F and resid_ok),
            "exact_limit_kcat_to_inf_G": exact_G, "exact_limit_kcat_to_inf_F_A": exact_F,
            "limit_rates": "k_r = phi k_on C_T, k_d = k_onD D_T, k_dis = k_onE E_T (capture-limited)",
            "low_load_residual": str(resid),
            "low_load_residual_is_bound_client_dilution": resid_ok}


def b5_wo04_v0(n=1000, seed=100):
    bf = load_source("WO-04", "bifurcation", "wo04_bifurcation")
    rng = np.random.default_rng(seed)
    counts, maxdiff, maxslope = [], 0.0, 0.0
    for k in range(n):
        p = bf.sample(rng, "V0")
        rs = rd.roots(p, finite=False, seq=False)
        counts.append(len(rs))
        if k < 100:
            Ug = np.geomspace(1e-6, p["s_P"] / p["mu"], 50)
            g_ours = rd.G(Ug, p, False, False)
            g_wo4 = bf.G(Ug, p)
            maxdiff = max(maxdiff, float(np.max(np.abs(g_ours - g_wo4)) / p["s_P"]))
            KM = (p["k_off"] + p["k_cat"] + p["mu"]) / p["k_on"]
            wo4 = (-(p["phi"] * p["k_cat"] * p["mu"] / (p["k_mis"] + p["mu"]) + p["mu"]) * p["C_T"] * KM / (KM + Ug)**2
                   - (p["k_d"] + p["mu"]) - 2 * p["k_a"] * Ug * (p["k_dA"] + p["mu"]) / (p["k_dis"] + p["k_dA"] + p["mu"]))
            ours = rd.slope_terms(Ug, p, False, False)["G_U"]
            maxslope = max(maxslope, float(np.max(np.abs(ours - wo4) / np.abs(wo4))))
    hist = {int(a): int(b) for a, b in zip(*np.unique(counts, return_counts=True))}
    ok = hist == {1: n} and maxdiff < 1e-9 and maxslope < 1e-9
    return {"pass": ok, "count_hist": hist, "wo04_reported": {"1": 1000},
            "max_rel_G_difference_vs_WO04": maxdiff, "max_rel_slope_difference_vs_WO04_formula": maxslope}


def b6_stage2_numerics(n=300, seed=21, n_int=8):
    rng = np.random.default_rng(seed)
    names, f, J = m.lambdas("S2")
    hist, max_res, n_unstable, max_fd, n_pos_term, n_fd_points = {}, 0.0, 0, 0.0, 0, 0
    for _ in range(n):
        p = m.sample_stage2(rng)
        rs = rd.roots(p, finite=True, seq=False)
        hist[len(rs)] = hist.get(len(rs), 0) + 1
        for r in rs:
            x = rd.full_state(r, p, "S2")
            res = float(np.max(np.abs(np.asarray(f(x, m.pvec(p)), float))) / p["s_P"])
            max_res = max(max_res, res)
            ev = np.linalg.eigvals(np.asarray(J(x, m.pvec(p)), float))
            n_unstable += int(np.any(ev.real >= 0))
        Ug = np.geomspace(1e-8, 0.5 * p["s_P"] / p["mu"], 40)
        T = rd.slope_terms(Ug, p, True, False)
        # Richardson-extrapolated central difference (truncation O(h^4)); compared
        # only where the roundoff estimate is small. G is a difference of terms of
        # size ~S_0, so eps_mach * S_0 / h is its roundoff scale
        h = 1e-3 * Ug
        d1 = (rd.G(Ug + h, p, True, False) - rd.G(Ug - h, p, True, False)) / (2 * h)
        d2 = (rd.G(Ug + 2 * h, p, True, False) - rd.G(Ug - 2 * h, p, True, False)) / (4 * h)
        fd = (4 * d1 - d2) / 3
        S0v = rd.consts(p)["S_0"]
        ok_pts = 1e-15 * S0v / (h * np.abs(T["G_U"])) < 1e-8
        n_fd_points += int(ok_pts.sum())
        if ok_pts.any():
            max_fd = max(max_fd, float(np.max(np.abs(fd - T["G_U"])[ok_pts] / np.abs(T["G_U"])[ok_pts])))
        n_pos_term += sum(int(np.any(T[k] > 0)) for k in ("T_B", "T_DU", "T_U", "T_A", "T_EA"))
    # long integrations from a far initial state must land on the reduced root
    rng2 = np.random.default_rng(seed + 1)
    max_int = 0.0
    for _ in range(n_int):
        p = m.sample_stage2(rng2)
        r = rd.roots(p, True, False)[0]
        xs = rd.full_state(r, p, "S2")
        idx = {nm: i for i, nm in enumerate(names)}
        x0 = np.zeros(len(names))
        x0[idx["N"]] = 0.3 * p["P_T"]
        x0[idx["A"]] = 0.3 * p["P_T"]
        for mname, tot in (("C", "C_T"), ("D", "D_T"), ("E", "E_T")):
            x0[idx[mname]] = p[tot]
        pv = m.pvec(p)
        sol = solve_ivp(lambda t, y: np.asarray(f(y, pv), float), (0, 80 / p["mu"]), x0,
                        method="LSODA", jac=lambda t, y: np.asarray(J(y, pv), float),
                        rtol=1e-10, atol=1e-12 * p["P_T"])
        max_int = max(max_int, float(np.max(np.abs(sol.y[:, -1] - xs)) / p["P_T"]))
    ok = (hist == {1: n} and max_res < 1e-8 and n_unstable == 0 and max_fd < 1e-6 and n_fd_points > 0.5 * 40 * n
          and n_pos_term == 0 and max_int < 1e-6)
    return {"pass": ok, "count_hist": hist, "max_rel_residual_full_9pool": max_res,
            "n_roots_with_eigenvalue_Re>=0": n_unstable, "max_rel_slope_vs_finite_difference": max_fd,
            "n_slope_points_compared": n_fd_points,
            "n_sets_with_any_positive_slope_term": n_pos_term,
            "n_integrations": n_int, "max_integration_vs_reduction_rel_PT": max_int}


def b7_ceiling_vs_fold():
    out = {}
    # (a) symbolic: at mu = 0 (machine totals now measured constants) the U
    # balance on the manifold contains no chaperone parameter
    M = manifold_S2()
    G2 = S0 - kB * M[m.B] - (m.k_catD + m.mu) * M[m.DU] - m.mu * m.U - (m.k_dA + m.mu) * m.A - m.mu * M[m.EA]
    FA2 = m.assemble("S2")["A"].subs(M)
    G0 = sp.simplify(G2.subs(m.mu, 0))
    F0 = sp.simplify(FA2.subs(m.mu, 0))
    chap = {CT, m.k_on, m.k_off, m.k_cat, m.phi}
    out["mu0_G_free_of_chaperone"] = not (G0.free_symbols & chap)
    out["mu0_FA_free_of_chaperone"] = not (F0.free_symbols & chap)
    out["mu0_G"] = str(G0)
    # with k_mis = 0 as well, the native pool has no exit: dN/dt > 0 forever
    dN0 = m.assemble("S2")["N"].subs({m.mu: 0, m.k_mis: 0})
    ra, rb = sp.symbols("r_a r_b", positive=True)
    out["mu0_kmis0_native_has_no_steady_state"] = bool(
        all_nonneg_coeffs(dN0.subs(m.eps, ra / (ra + rb))) and sp.simplify(dN0.subs(m.B, 0)) != 0)
    # (b) closed-form ceiling at mu = 0, k_dA = 0
    VD, VE, KD0, KE0, ka = 2.0, 0.5, 5.0, 3.0, 1e-3
    Uc = np.sqrt(VE / ka)
    s_crit = VD * Uc / (KD0 + Uc)
    out["mu0_s_crit"] = s_crit
    base = {"eps": 0.05, "k_mis": 1e-5, "k_a": ka, "k_dA": 0.0, "phi": 0.8,
            "k_on": 1.0, "k_off": 1.0, "k_cat": 0.1, "C_T": 20.0,
            "k_onD": 1.0, "k_offD": KD0 - 1.0, "k_catD": 1.0, "D_T": VD,
            "k_onE": 1.0, "k_offE": KE0 - 1.0, "k_catE": 1.0, "E_T": VE}
    # mu = 0 branch: U* = K_MD s_P/(V_D - s_P), slope -V_D K_MD/(K_MD+U)^2
    sweep = []
    for r in (0.5, 0.9, 0.99, 0.999, 0.9999):
        sP = r * s_crit
        U0 = KD0 * sP / (VD - sP)
        slope0 = -VD * KD0 / (KD0 + U0)**2
        A0 = KE0 * ka * U0**2 / (VE - ka * U0**2)
        sweep.append({"s_P/s_crit": r, "U*": U0, "A*": A0, "G_U": slope0,
                      "aggregate_mode_rate": -VE * KE0 / (KE0 + A0)**2})
    out["mu0_branch"] = sweep
    out["mu0_branch_unique_and_GU_negative"] = all(s["G_U"] < 0 for s in sweep)
    out["mu0_A_diverges_at_s_crit"] = sweep[-1]["A*"] > 1e3 * sweep[0]["A*"]
    # (c) mu > 0: a unique root for every s_P, including above the mu = 0 ceiling
    mu_rows = []
    for mu_ in (1e-4, 1e-6, 1e-8):
        for r in (0.5, 0.99, 1.5, 3.0):
            p = {**base, "mu": mu_, "s_P": r * s_crit}
            rs = rd.roots(p, finite=True, seq=False, n=40000)
            Ug = np.geomspace(1e-10, p["s_P"] / mu_, 400)
            maxGU = float(np.max(rd.slope_terms(Ug, p, True, False)["G_U"]))
            mu_rows.append({"mu": mu_, "s_P/s_crit": r, "n_roots": len(rs),
                            "U*": rs[0] if rs else None, "max_G_U_on_grid": maxGU})
    out["mu_pos_rows"] = mu_rows
    out["mu_pos_always_one_root_no_tangency"] = all(r["n_roots"] == 1 and r["max_G_U_on_grid"] < 0 for r in mu_rows)
    ok = (out["mu0_G_free_of_chaperone"] and out["mu0_FA_free_of_chaperone"]
          and out["mu0_kmis0_native_has_no_steady_state"] and out["mu0_branch_unique_and_GU_negative"]
          and out["mu0_A_diverges_at_s_crit"] and out["mu_pos_always_one_root_no_tangency"])
    return {"pass": bool(ok), **out}


def run():
    return {"B1_stage1_closed_form": b1_stage1(), "B2_identity_S2": b2_identity_S2(),
            "B3_slope_S2": b3_slope_S2(), "B4_limits": b4_limits(), "B5_wo04_v0": b5_wo04_v0(),
            "B6_stage2_numerics": b6_stage2_numerics(), "B7_ceiling_vs_fold": b7_ceiling_vs_fold()}


if __name__ == "__main__":
    import json
    r = run()
    print(json.dumps({k: v["pass"] for k, v in r.items()}, indent=1))
    sys.exit(0 if all(v["pass"] for v in r.values()) else 1)
