"""
permanent checks: stage 4, feedback A (aggregate-mediated chaperone sequestration).

symbolic
  F1  on the manifold where N, B, C, CA, DU, D, EA, E are at their own steady
      states, dU/dt + dA/dt = G4 = G2 - mu CA, and equals the client balance
  F2  allocation: C_T = C + B + CA with C = C_T/(1 + U/K_M + A/K_Ae) and
      K_Ae = (k_offA + mu)/k_onA  (dilution enters; it is not a K_d)
  F3  A(U) is still unique and increasing (F_U > 0, F_AA < 0)
  F4  dB/dU = C_T (1 + a - U da/dU)/(K_M w^2): B falls with load iff the
      elasticity of A in U exceeds 1 + C/CA
  F5  eps enters G only through S_0: every manifold pool is eps-free, so the
      fold condition G_U = 0 is one equation in U alone
  F6  at mu = 0 the U balance and the aggregate balance contain no chaperone
      parameter (k_onA, k_offA, C_T, K_M, ...): mechanism A cannot move U* or
      make a fold at mu = 0
numeric
  F7  independent reproduction of WO-04 V1 (seed 101, 1000 sets): this code's
      G (organising form, closed-form quadratic A) against WO-04's count and G
  F8  every fold of every multistable set: located by G_U = 0, eps* from
      G = 0; slope decomposed term by term; full 6-pool Jacobian has a zero
      eigenvalue there and det J changes sign across it
  F9  mu sweep on stored examples (s_P, C_T held): folds disappear as mu -> 0
  F10 nesting of the fold into the finite-protease/disaggregase stage S4:
      the fold location converges to the S4C (= V1) fold as k_cat -> inf
"""
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.optimize import brentq

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import model as m  # noqa: E402
import reduction as rd  # noqa: E402
from check_baseline import CT, DT, ET, Kc, KD, KE, S0, kB, all_nonneg_coeffs  # noqa: E402
from check_conservation import load_source  # noqa: E402

KA = (m.k_offA + m.mu) / m.k_onA


def manifold_S4():
    w = 1 + m.U / Kc + m.A / KA
    C_ = CT / w
    B_ = C_ * m.U / Kc
    CA_ = C_ * m.A / KA
    DU_ = DT * m.U / (KD + m.U)
    EA_ = ET * m.A / (KE + m.A)
    N_ = ((1 - m.eps) * m.s_P + m.phi * m.k_cat * B_) / (m.k_mis + m.mu)
    return {m.C: C_, m.B: B_, m.CA: CA_, m.DU: DU_, m.D: DT - DU_, m.EA: EA_, m.E: ET - EA_,
            m.N: N_, m.s_C: m.mu * CT, m.s_D: m.mu * DT, m.s_E: m.mu * ET}


def G4_expr(M):
    return (S0 - kB * M[m.B] - (m.k_catD + m.mu) * M[m.DU] - m.mu * m.U
            - (m.k_dA + m.mu) * m.A - m.mu * M[m.EA] - m.mu * M[m.CA])


def f1_f2_identity():
    rhs = m.assemble("S4")
    M = manifold_S4()
    zero = {n: sp.simplify(rhs[n].subs(M)) == 0 for n in ("N", "B", "C", "CA", "DU", "D", "EA", "E")}
    GU = rhs["U"].subs(M) + rhs["A"].subs(M)
    ident = sp.simplify(GU - G4_expr(M)) == 0
    massbal = sp.simplify(GU - m.expected_total_rate("client", "S4").subs(M)) == 0
    alloc = sp.simplify(M[m.C] + M[m.B] + M[m.CA] - CT) == 0
    # the complex balance fixes K_Ae with dilution included: k_onA C A = (k_offA + mu) CA
    kae = sp.simplify(rhs["CA"].subs(M)) == 0
    return {"pass": bool(all(zero.values()) and ident and massbal and alloc and kae),
            "other_balances_zero": zero, "G4_identity": ident, "G4_equals_client_balance": massbal,
            "allocation_C+B+CA=C_T": alloc, "K_Ae=(k_offA+mu)/k_onA": kae,
            "G4": "S_0 - kappa_B B - (k_catD+mu) DU - mu U - (k_dA+mu) A - mu EA - mu CA"}


def f3_f4_f5_f6():
    M = manifold_S4()
    FA = m.assemble("S4")["A"].subs(M)
    FU, FAA = sp.diff(FA, m.U), sp.diff(FA, m.A)
    # F_U = 2 k_a U - mu dCA/dU with dCA/dU < 0; F_AA = -(k_dA+mu) - (k_catE+mu) dEA/dA - mu dCA/dA
    dCA_U, dCA_A = sp.diff(M[m.CA], m.U), sp.diff(M[m.CA], m.A)
    FU_ok = sp.simplify(FU - (2 * m.k_a * m.U - m.mu * dCA_U)) == 0 and all_nonneg_coeffs(-dCA_U)
    FAA_ok = (sp.simplify(FAA + (m.k_dA + m.mu) + (m.k_catE + m.mu) * sp.diff(M[m.EA], m.A)
                          + m.mu * dCA_A) == 0 and all_nonneg_coeffs(dCA_A)
              and all_nonneg_coeffs(sp.diff(M[m.EA], m.A)))
    # F4: B criterion with a = A/K_Ae as an arbitrary function of U
    Kms, Kas, Cts = sp.symbols("K_M K_Ae C_T", positive=True)
    af = sp.Function("a")(m.U)
    Bf = Cts * (m.U / Kms) / (1 + m.U / Kms + af)
    w = 1 + m.U / Kms + af
    crit = sp.simplify(sp.diff(Bf, m.U) - Cts * (1 + af - m.U * sp.diff(af, m.U)) / (Kms * w**2)) == 0
    # F5: eps only in S_0
    pools_eps_free = all(m.eps not in M[k].free_symbols for k in (m.B, m.C, m.CA, m.DU, m.EA))
    FA_eps_free = m.eps not in sp.simplify(FA).free_symbols
    G = G4_expr(M)
    Geps = sp.simplify(sp.diff(G, m.eps) - m.s_P * m.mu / (m.k_mis + m.mu)) == 0
    # F6: mu = 0 (machine totals are measured constants)
    G0 = sp.simplify(G.subs(m.mu, 0))
    F0 = sp.simplify(FA.subs(m.mu, 0))
    chap = {CT, m.k_on, m.k_off, m.k_cat, m.phi, m.k_onA, m.k_offA}
    mu0 = not ((G0.free_symbols | F0.free_symbols) & chap)
    return {"pass": bool(FU_ok and FAA_ok and crit and pools_eps_free and FA_eps_free and Geps and mu0),
            "F3_F_U>0": FU_ok, "F3_F_AA<0": FAA_ok, "F4_B_slope_formula": crit,
            "F4_criterion": "dB/dU < 0  <=>  dlnA/dlnU > 1 + K_Ae/A = 1 + C/CA",
            "F5_pools_eps_free": pools_eps_free and FA_eps_free, "F5_G_eps_constant": Geps,
            "F6_mu0_free_of_chaperone_and_binding": mu0, "F6_mu0_G": str(G0), "F6_mu0_F_A": str(F0)}


def eps_curve(Ug, p, finite=False):
    """exact equilibrium curve: G = eps G_eps + H(U), so eps(U) = -H(U)/G_eps."""
    H = rd.G(Ug, p, finite, True, eps=0.0)
    return -H / rd.consts(p)["G_eps"]


def folds_of(p, finite=False, n=20000):
    Umax = p["s_P"] / p["mu"]
    Ug = np.geomspace(min(1e-12 * Umax, 1e-12), Umax, n)
    gu = rd.slope_terms(Ug, p, finite, True)["G_U"]
    out = []
    for i in np.where(np.sign(gu[:-1]) * np.sign(gu[1:]) < 0)[0]:
        uf = brentq(lambda x: float(rd.slope_terms(np.array([x]), p, finite, True)["G_U"][0]),
                    Ug[i], Ug[i + 1], xtol=1e-14 * Ug[i + 1], rtol=1e-15)
        out.append(uf)
    return out


def f7_wo04_v1(n=1000, seed=101):
    bf = load_source("WO-04", "bifurcation", "wo04_bifurcation")
    rng = np.random.default_rng(seed)
    counts, multi, maxdiff = [], [], 0.0
    for k in range(n):
        p = bf.sample(rng, "V1")
        rs = rd.roots(p, finite=False, seq=True)
        counts.append(len(rs))
        if len(rs) >= 2:
            multi.append((k, p, rs))
        if k < 100:
            Ug = np.geomspace(1e-6, p["s_P"] / p["mu"], 50)
            d = np.abs(rd.G(Ug, p, False, True) - bf.G(Ug, p)) / p["s_P"]
            maxdiff = max(maxdiff, float(np.max(d)))
    hist = {int(a): int(b) for a, b in zip(*np.unique(counts, return_counts=True))}
    return hist, multi, maxdiff


def classify_fold(uf, p):
    names, f, J = m.lambdas("S4C")
    ef = float(eps_curve(np.array([uf]), p)[0])
    q = {**p, "eps": ef}
    T = rd.slope_terms(np.array([uf]), q, False, True)
    x = rd.full_state(uf, q, "S4C")
    ev = np.linalg.eigvals(np.asarray(J(x, m.pvec(q)), float))
    dets = []
    for s in (1 - 1e-3, 1 + 1e-3):
        u2 = uf * s
        e2 = float(eps_curve(np.array([u2]), p)[0])
        x2 = rd.full_state(u2, {**p, "eps": e2}, "S4C")
        dets.append(np.sign(np.linalg.det(np.asarray(J(x2, m.pvec({**p, "eps": e2})), float))))
    elas = float(T["Aprime"][0] * uf / T["A"][0]) if T["A"][0] > 0 else np.nan
    CoverCA = float(T["C"][0] / T["CA"][0]) if T["CA"][0] > 0 else np.inf
    terms = {k: float(T[k][0]) for k in ("T_B", "T_U", "T_A", "T_CA")}
    scale = max(abs(v) for v in terms.values())
    return {"U_fold": uf, "eps_fold": ef, "physical": bool(0 < ef < 1),
            "terms": terms, "G_U_rel": float(T["G_U"][0]) / scale,
            "positive_terms": [k for k, v in terms.items() if v > 1e-12 * scale],
            "dB_dU_negative": bool(T["dB"][0] < 0),
            "elasticity_A": elas, "one_plus_C_over_CA": 1 + CoverCA,
            "criterion_consistent": bool((T["dB"][0] < 0) == (elas > 1 + CoverCA)),
            "CA_over_CT": float(T["CA"][0] / p["C_T"]), "CA_over_C": 1 / CoverCA,
            "min_abs_eig_over_mu": float(np.min(np.abs(ev)) / p["mu"]),
            "det_sign_change": bool(dets[0] != dets[1])}


def f8_folds(multi):
    rows = []
    for k, p, rs in multi:
        for uf in folds_of(p):
            r = classify_fold(uf, p)
            r["idx"] = k
            rows.append(r)
    phys = [r for r in rows if r["physical"]]
    per_set = {}
    for r in phys:
        per_set[r["idx"]] = per_set.get(r["idx"], 0) + 1
    summary = {
        "n_multistable_sets": len(multi), "n_folds_total": len(rows), "n_folds_physical": len(phys),
        "physical_folds_per_set_hist": {int(a): int(b) for a, b in zip(*np.unique(list(per_set.values()), return_counts=True))} if per_set else {},
        "n_sets_with_no_physical_fold": len(multi) - len(per_set),
        "all_folds_G_U_zero": all(abs(r["G_U_rel"]) < 1e-8 for r in rows),
        "all_folds_zero_eigenvalue": all(r["min_abs_eig_over_mu"] < 1e-5 for r in rows),
        "all_folds_det_sign_change": all(r["det_sign_change"] for r in rows),
        "n_folds_T_B_positive": sum("T_B" in r["positive_terms"] for r in rows),
        "n_folds_T_CA_positive": sum("T_CA" in r["positive_terms"] for r in rows),
        "n_folds_only_T_CA_positive": sum(r["positive_terms"] == ["T_CA"] for r in rows),
        "n_folds_T_U_or_T_A_positive": sum(("T_U" in r["positive_terms"]) or ("T_A" in r["positive_terms"]) for r in rows),
        "all_criterion_consistent": all(r["criterion_consistent"] for r in rows),
        "min_CA_over_CT_at_fold": min(r["CA_over_CT"] for r in rows) if rows else None,
        "min_CA_over_C_at_fold": min(r["CA_over_C"] for r in rows) if rows else None,
        "median_CA_over_CT_at_fold": float(np.median([r["CA_over_CT"] for r in rows])) if rows else None,
    }
    # cross-check the per-example physical fold count against WO-04's stored
    # continuation results (wo04_results.json, read-only)
    import json
    stored = json.load(open(HERE.parents[1] / "proteostasis_rebuild" / "WO-04" / "wo04_results.json"))
    wo4 = {ex["idx"]: sum(f["physical"] for f in ex["folds"]) for ex in stored["G4.5_continuation"]["V1"]}
    ours = {i: sum(1 for r in phys if r["idx"] == i) for i in wo4}
    summary["wo04_stored_physical_folds_by_idx"] = wo4
    summary["this_code_physical_folds_by_idx"] = ours
    summary["agrees_with_wo04_stored"] = ours == wo4
    return summary, rows


def f9_mu_sweep(multi, n_examples=5):
    """lower mu with s_P and C_T held (a mathematical sweep, not a growth
    condition). the fold needs kappa_B |B'| to beat k_d + mu; kappa_B ~
    phi k_cat mu/k_mis once mu << k_mis, so the folds must vanish, but only far
    below k_mis. reported: the largest mu at which each example has no fold."""
    out = []
    for k, p, _ in multi[:n_examples]:
        row = {"idx": k, "mu_native": p["mu"], "k_mis": p["k_mis"], "k_d": p["k_d"],
               "phi_kcat_over_kmis": p["phi"] * p["k_cat"] / p["k_mis"], "n_folds_by_mu": {}}
        first_zero = None
        for e in range(0, 26, 2):
            q = {**p, "mu": p["mu"] * 10.0**(-e)}        # s_P and C_T held fixed
            nf = len(folds_of(q, n=8000))
            row["n_folds_by_mu"][f"{q['mu']:.3e}"] = nf
            if nf == 0 and first_zero is None:
                first_zero = q["mu"]
        row["largest_mu_without_fold"] = first_zero
        out.append(row)
    vanish = all(list(r["n_folds_by_mu"].values())[-1] == 0 for r in out)
    return {"rows": out, "folds_absent_at_smallest_mu": vanish}


def f10_nesting(multi, n_examples=3):
    out = []
    for k, p, _ in multi[:n_examples]:
        ref = sorted(folds_of(p))
        row = {"idx": k, "U_folds_S4C": ref, "rel_shift": {}}
        for lam in (1e1, 1e3, 1e5):
            q = {**p, "D_T": 1.0, "k_onD": p["k_d"], "k_offD": 0.0, "k_catD": lam,
                 "E_T": 1.0, "k_onE": max(p["k_dis"], 1e-300), "k_offE": 0.0, "k_catE": lam}
            fs = sorted(folds_of(q, finite=True))
            if len(fs) == len(ref) and ref:
                row["rel_shift"][f"{lam:.0e}"] = float(max(abs(a - b) / b for a, b in zip(fs, ref)))
            else:
                row["rel_shift"][f"{lam:.0e}"] = f"fold count {len(fs)} vs {len(ref)}"
        out.append(row)
    conv = all(isinstance(r["rel_shift"]["1e+05"], float) and r["rel_shift"]["1e+05"] < 1e-3 for r in out)
    return {"rows": out, "converges_to_S4C_fold": conv}


def run():
    res = {"F1_F2_identity_allocation": f1_f2_identity(), "F3_to_F6_symbolic": f3_f4_f5_f6()}
    hist, multi, maxdiff = f7_wo04_v1()
    res["F7_wo04_v1"] = {"pass": hist == {1: 946, 3: 54} and maxdiff < 1e-9, "count_hist": hist,
                         "wo04_reported": {"1": 946, "3": 54}, "max_rel_G_difference_vs_WO04": maxdiff}
    summ, rows = f8_folds(multi)
    res["F8_folds"] = {"pass": bool(summ["all_folds_G_U_zero"] and summ["all_folds_zero_eigenvalue"]
                                    and summ["all_folds_det_sign_change"] and summ["all_criterion_consistent"]
                                    and summ["n_folds_T_B_positive"] + summ["n_folds_only_T_CA_positive"] == summ["n_folds_total"]
                                    and summ["n_folds_T_U_or_T_A_positive"] == 0
                                    and summ["agrees_with_wo04_stored"]),
                       **summ, "rows": rows}
    sw = f9_mu_sweep(multi)
    res["F9_mu_sweep"] = {"pass": sw["folds_absent_at_smallest_mu"], **sw}
    ns = f10_nesting(multi)
    res["F10_nesting_S4_to_S4C"] = {"pass": ns["converges_to_S4C_fold"], **ns}
    return res


if __name__ == "__main__":
    import json
    r = run()
    print(json.dumps({k: v["pass"] for k, v in r.items()}, indent=1))
    sys.exit(0 if all(v["pass"] for v in r.values()) else 1)
