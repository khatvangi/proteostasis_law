#!/usr/bin/env python3
"""WO-03 analysis: gates G3.1, G3.2, G3.3, G3.5. writes wo03_results.json."""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import chaperone as ch  # noqa: E402
m = ch.m


def g31_equilibrium():
    C_T, K_d = 50.0, 1.0
    M_T = np.linspace(0, 300, 3001)
    res = ch.mass_balance_residuals(C_T, M_T[1:], K_d)
    Cf, Mf, b = ch.free_exact(C_T, 50.0, K_d)
    approx = ch.free_legacy(C_T, 50.0, K_d)
    fold = lambda c: c / (c + K_d)   # legacy v_fold shape with k_obs_max = 1
    return {"max_residual_chaperone": float(res["chaperone"].max()),
            "max_residual_client": float(res["client"].max()),
            "max_rel_residual_mass_action": float(res["mass_action"].max()),
            "at_M_T_50": {"C_f_exact": Cf, "C_b_exact": b, "C_f_legacy": approx,
                          "fold_ratio_exact_over_legacy": fold(Cf) / fold(approx)},
            "audit_C_f_reference": 6.58872, "audit_fold_ratio_reference": 1.75382}


def g32_cycle():
    """closed driven cycle: steady-state bound complex vs quadratic with K_M
    (predicted) and with K_d (equilibrium benchmark)."""
    C_T, M_T, k_on, k_off = 50.0, 50.0, 1.0, 1.0
    rows = []
    for k_cat in (1e-3, 1e-2, 1e-1, 1.0, 10.0):
        sol = solve_ivp(ch.cycle_closed_rhs, (0, 1e4), [C_T, M_T, 0.0],
                        args=(k_on, k_off, k_cat), method="LSODA",
                        rtol=1e-12, atol=1e-14)
        B_ode = float(sol.y[2, -1])
        KM = ch.K_M(k_on, k_off, k_cat)
        B_KM = float(ch.bound_exact(C_T, M_T, KM))
        B_Kd = float(ch.bound_exact(C_T, M_T, k_off / k_on))
        rows.append({"k_cat_over_k_off": k_cat / k_off, "K_M": KM, "K_d": k_off / k_on,
                     "B_ode": B_ode, "B_pred_K_M": B_KM, "B_pred_K_d": B_Kd,
                     "rel_err_K_M": abs(B_ode - B_KM) / B_ode,
                     "rel_err_K_d": abs(B_ode - B_Kd) / B_ode})

    # four-state DnaK-like cycle. illustrative rates only: ATP state fast and
    # weak (K_dT = 10 uM), ADP state slow and tight (K_dD = 1 uM).
    base = dict(k_onT=10.0, k_offT=100.0, k_onD=0.01, k_offD=0.01,
                k_h=10.0, k_h0=0.01, k_ex=0.1)
    K_dT = base["k_offT"] / base["k_onT"]
    K_dD = base["k_offD"] / base["k_onD"]
    U = 1.0
    occ_drv, _ = ch.four_state_occupancy(U, **base)
    # undriven control: choose k_h so the cycle obeys detailed balance
    # (Kolmogorov: k_onT k_h k_offD k_ex = k_h0 k_onD k_offT k_ex)
    k_h_db = base["k_h0"] * base["k_onD"] * base["k_offT"] / (base["k_onT"] * base["k_offD"])
    occ_eq, _ = ch.four_state_occupancy(U, **{**base, "k_h": k_h_db})
    four = {"K_dT": K_dT, "K_dD": K_dD, "U": U,
            "K_eff_driven": ch.K_eff(U, occ_drv), "occupancy_driven": occ_drv,
            "k_h_detailed_balance": k_h_db,
            "K_eff_detailed_balance": ch.K_eff(U, occ_eq), "occupancy_detailed_balance": occ_eq}
    four["driven_tighter_than_both_states"] = four["K_eff_driven"] < min(K_dT, K_dD)
    four["undriven_within_state_range"] = (min(K_dT, K_dD) - 1e-9
                                           <= four["K_eff_detailed_balance"]
                                           <= max(K_dT, K_dD) + 1e-9)
    return {"simple_cycle": rows, "four_state": four}


def g33_competition():
    okP, okC = ch.ext_conservation()
    sym_same, sym_stays = ch.ext_reduces_symbolically()
    p = m.scenario_params()
    ext = {**p, "nu_c": 0.0, "k_onX": 1.0, "k_offX": 1.0, "k_catX": 0.05,
           "k_fX": 0.01, "k_xu": 1e-3}
    # reduction test: nu_c = 0 and X = BX = 0 must reproduce WO-02 exactly
    x0 = np.array([2900.0, 1.0, 0.5, 0.1, 49.5])
    T = 5.0 / p["mu"]
    te = np.linspace(0, T, 200)
    s_base = solve_ivp(m.rhs, (0, T), x0, args=(m.pvec(p),), method="LSODA",
                       jac=m.jac, rtol=1e-11, atol=1e-12, t_eval=te)
    s_ext = solve_ivp(ch.ext_rhs, (0, T), np.append(x0, [0.0, 0.0]),
                      args=(ch.ext_pvec(ext),), method="LSODA", jac=ch.ext_jac,
                      rtol=1e-11, atol=1e-12, t_eval=te)
    red_err = float(np.max(np.abs(s_ext.y[:5] - s_base.y) / np.maximum(np.abs(s_base.y), 1e-6)))
    ext_nascent_max = float(np.max(np.abs(s_ext.y[5:])))

    # theta as an output: nascent client load swept via nu_c
    rows = []
    for nu in (0.0, 0.05, 0.1, 0.2, 0.4):
        q = ch.ext_pvec({**ext, "nu_c": nu})
        xs, resid = ch.steady_state(ch.ext_rhs, ch.ext_jac, q,
                                    np.append(x0, [1.0, 1.0]))
        N, U, B, A, C, X, BX = xs
        CT = C + B + BX
        rows.append({"nu_c": nu, "theta_BX_over_C_T": BX / CT, "C_free": C,
                     "U": U, "A": A, "B": B, "A_over_P_T": A / (N + U + B + A + X + BX),
                     "max_abs_rhs": resid})
    return {"P_T_identity": okP, "C_T_identity": okC,
            "symbolic_reduction_equations_identical": sym_same,
            "symbolic_reduction_nascent_stays_zero": sym_stays,
            "reduction_max_rel_err": red_err, "reduction_nascent_max": ext_nascent_max,
            "theta_sweep": rows}


def g35_qss():
    """finite-pool tQSSA of the WO-02 model (B eliminated via K_M) against the
    full ODE, from the same initial condition with B at its QSS value."""
    p = m.scenario_params(eps=0.04)
    pv = m.pvec(p)
    C_T = p["C_T"]
    KM = ch.K_M(p["k_on"], p["k_off"], p["k_cat"], p["mu"])

    def reduced(t, y):
        N, W, A = y
        B = ch.bound_exact(C_T, W, KM)
        U = W - B
        dN = (1 - p["eps"]) * p["s_P"] + p["phi"] * p["k_cat"] * B - p["k_mis"] * N - p["mu"] * N
        dW = (p["eps"] * p["s_P"] + p["k_mis"] * N - p["phi"] * p["k_cat"] * B
              - p["k_d"] * U - p["k_a"] * U * U + p["k_dis"] * A - p["mu"] * W)
        dA = p["k_a"] * U * U - (p["k_dis"] + p["k_dA"]) * A - p["mu"] * A
        return [dN, dW, dA]

    W0, N0, A0 = 5.0, 2990.0, 0.0
    B0 = float(ch.bound_exact(C_T, W0, KM))
    x0 = np.array([N0, W0 - B0, B0, A0, C_T - B0])
    T = 10.0 / p["mu"]
    te = np.linspace(0, T, 400)
    full = solve_ivp(m.rhs, (0, T), x0, args=(pv,), method="LSODA", jac=m.jac,
                     rtol=1e-11, atol=1e-12, t_eval=te)
    red = solve_ivp(reduced, (0, T), [N0, W0, A0], method="LSODA",
                    rtol=1e-11, atol=1e-12, t_eval=te)
    U_full, A_full = full.y[1], full.y[3]
    B_red = ch.bound_exact(C_T, red.y[1], KM)
    U_red, A_red = red.y[1] - B_red, red.y[2]
    sl = te > 0.1 / p["mu"]    # compare after the first 10% of a doubling
    return {"K_M": KM, "K_d": p["k_off"] / p["k_on"],
            "max_rel_err_U": float(np.max(np.abs(U_red - U_full)[sl] / U_full[sl])),
            "max_rel_err_A": float(np.max(np.abs(A_red - A_full)[sl] / np.maximum(A_full[sl], 1e-12))),
            "steady_U_full": float(U_full[-1]), "steady_A_over_P_T": float(A_full[-1] / full.y[:4, -1].sum())}


def run():
    r = {"G3.1": g31_equilibrium(), "G3.2": g32_cycle(),
         "G3.3": g33_competition(), "G3.5": g35_qss()}
    a = r["G3.1"]["at_M_T_50"]
    sc = r["G3.2"]["simple_cycle"]
    r["pass"] = {
        "G3.1": (r["G3.1"]["max_residual_chaperone"] < 1e-10
                 and r["G3.1"]["max_residual_client"] < 1e-10
                 and r["G3.1"]["max_rel_residual_mass_action"] < 1e-10
                 and abs(a["C_f_exact"] - 6.58872) < 5e-6
                 and abs(a["fold_ratio_exact_over_legacy"] - 1.75382) < 5e-6),
        "G3.2": (all(x["rel_err_K_M"] < 1e-6 for x in sc)
                 and sc[-1]["rel_err_K_d"] > 0.01),
        "G3.3": (r["G3.3"]["P_T_identity"] and r["G3.3"]["C_T_identity"]
                 and r["G3.3"]["symbolic_reduction_equations_identical"]
                 and r["G3.3"]["symbolic_reduction_nascent_stays_zero"]
                 and r["G3.3"]["reduction_max_rel_err"] < 1e-6
                 and r["G3.3"]["reduction_nascent_max"] == 0.0),
        "G3.5": r["G3.5"]["max_rel_err_U"] < 1e-2 and r["G3.5"]["max_rel_err_A"] < 1e-2,
    }
    r["pass"] = {k: bool(v) for k, v in r["pass"].items()}
    return r


if __name__ == "__main__":
    r = run()
    print(json.dumps(r, indent=1, default=float))
    (HERE / "wo03_results.json").write_text(json.dumps(r, indent=2, default=float))
    sys.exit(0 if all(r["pass"].values()) else 1)
