#!/usr/bin/env python3
"""
WO-02 gates G2.1-G2.5. writes proofs.json.

G2.1 symbolic conservation of P_T and C_T
G2.2 dimensional check of every RHS term (WO-01 checker)
G2.3 forward invariance of the nonnegative orthant, symbolically
G2.4 numerical integration: totals conserved, no negative states
G2.5 legacy Phi inflow: donor test and phantom mass rate at its operating point
"""
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import model as m  # noqa: E402


def g21_conservation():
    dPT = sum(m.PROTEIN_W[s] * r for s, r in zip(m.STATE_NAMES, m.RHS_SYM))
    dCT = sum(m.CHAP_W[s] * r for s, r in zip(m.STATE_NAMES, m.RHS_SYM))
    want_P = m.s_P - m.k_d * m.U - m.k_dA * m.A - m.mu * m.P_T
    want_C = m.s_C - m.mu * m.C_T
    return {"P_T_identity": sp.simplify(dPT - want_P) == 0,
            "C_T_identity": sp.simplify(dCT - want_C) == 0}


def g22_dimensions():
    units = dict(m.var.UNITS)
    bad = []
    for s, r in zip(m.STATE_NAMES, m.RHS_SYM):
        for term in sp.Add.make_args(sp.expand(r)):
            try:
                got = m.un.unit_of(term, units)
                if not got.same(m.un.UM_PER_S):
                    bad.append(f"d{s}/dt term {term}: {got}")
            except m.un.DimensionError as e:
                bad.append(f"d{s}/dt term {term}: {e}")
    return {"n_bad_terms": len(bad), "bad": bad}


def g23_invariance():
    """on face x_i = 0, dx_i/dt must be >= 0 for all nonnegative states and
    parameters. eps and phi live in [0,1], so write eps = a/(a+b), phi =
    c/(c+d) with a,b,c,d >= 0; then (a+b)(c+d) * face expression must be a
    polynomial with only nonnegative coefficients in nonnegative symbols."""
    a, b, c, d = sp.symbols("a b c d", nonnegative=True)
    out = {}
    for s, x, r in zip(m.STATE_NAMES, m.X, m.RHS_SYM):
        face = r.subs(x, 0).subs({m.eps: a / (a + b), m.phi: c / (c + d)})
        num = sp.expand(sp.cancel(face * (a + b) * (c + d)))
        if num == 0:
            out[s] = True
            continue
        coeffs = sp.Poly(num, *num.free_symbols).coeffs()
        out[s] = all(co >= 0 for co in coeffs)
    return out


def g24_numeric(n=200, seed=20260927):
    rng = np.random.default_rng(seed)
    worst_P, worst_C, worst_neg, n_fail = 0.0, 0.0, 0.0, 0

    def logu(lo, hi):
        return float(np.exp(rng.uniform(np.log(lo), np.log(hi))))

    for _ in range(n):
        p = {"s_P": logu(1e-2, 1e1), "eps": rng.uniform(0, 1),
             "s_C": logu(1e-3, 1e-1), "mu": logu(1e-5, 1e-3),
             "k_on": logu(1e-2, 1e1), "k_off": logu(1e-2, 1e1),
             "k_cat": logu(1e-3, 1e0), "phi": rng.uniform(0, 1),
             "k_d": logu(1e-5, 1e-2), "k_a": logu(1e-5, 1e-1),
             "k_dis": logu(1e-5, 1e-2), "k_dA": logu(1e-6, 1e-3),
             "k_mis": logu(1e-7, 1e-4)}
        pv = m.pvec(p)
        x0 = rng.uniform(0, 100, size=5)
        # augmented state: 5 species + Q, where dQ/dt is the claimed total
        # protein balance written independently of the species equations
        def aug(t, y):
            dx = m.rhs(t, y[:5], pv)
            dQ = p["s_P"] - p["k_d"] * y[1] - p["k_dA"] * y[3] - p["mu"] * y[5]
            return np.append(dx, dQ)
        T = 3.0 / p["mu"]
        y0 = np.append(x0, x0[:4].sum())
        sol = solve_ivp(aug, (0, T), y0, method="LSODA", rtol=1e-10, atol=1e-12,
                        t_eval=np.linspace(0, T, 50))
        if not sol.success:
            n_fail += 1
            continue
        Y = sol.y
        PT = Y[:4].sum(0)
        CT_num = Y[2] + Y[4]
        CT_exact = p["s_C"] / p["mu"] + (x0[2] + x0[4] - p["s_C"] / p["mu"]) * np.exp(-p["mu"] * sol.t)
        worst_P = max(worst_P, float(np.max(np.abs(PT - Y[5]) / np.maximum(np.abs(Y[5]), 1e-9))))
        worst_C = max(worst_C, float(np.max(np.abs(CT_num - CT_exact) / np.maximum(CT_exact, 1e-9))))
        worst_neg = min(worst_neg, float(Y[:5].min()))
    return {"n": n, "n_integration_failures": n_fail,
            "max_rel_err_P_T": worst_P, "max_rel_err_C_T": worst_C,
            "most_negative_state": worst_neg}


def g25_phi_phantom():
    """legacy inflow J_bare*Phi(P): the (Phi-1) part has no donor. quantify it
    at the legacy's own operating point (usage-weighted mu, baseline params)."""
    sys.path.insert(0, "/storage/kiran-stuff/proteostasis_law/envelope-paper/scripts/vendor")
    import two_pool_ode as L
    p = L.Params()
    f = 6.334247974475959e-4
    J = f * p.N_prot * (1 - p.S_avg) * p.p_baseline / p.T_gen_s
    P_star, A_star = L.steady_state(J, p)
    P_dag, J_crit, mech, _ = L.saddle_node_operational(L.J_curve_two, L.A_qs, p)
    phi_star = float(L.phi(P_star, p))
    phi_dag = float(L.phi(P_dag, p))
    # donor test, computed: sum the legacy pool equations (WO-01 transcription
    # of two_pool_ode.py:19-20). anything left after removing the named
    # external source J_bare and the named sinks (degradation, folding to the
    # untracked native pool, aggregate clearance) is flux with no donor.
    sys.path.insert(0, str(HERE.parent / "WO-01"))
    import legacy_units_audit as LA
    total = LA.dPdt + LA.dAdt
    named = LA.J_bare - LA.R - LA.k_clear * LA.A
    residual = sp.simplify(total - named)
    expected_phantom = sp.simplify(LA.J_bare * (LA.phi - 1))
    donor_found = residual == 0
    return {"P_star": P_star, "Phi_at_P_star": phi_star,
            "phantom_fraction_of_inflow_at_P_star": (phi_star - 1) / phi_star,
            "P_dagger": P_dag, "Phi_at_P_dagger": phi_dag,
            "phantom_fraction_of_inflow_at_P_dagger": (phi_dag - 1) / phi_dag,
            "phantom_rate_at_P_dagger_per_s": J_crit * (phi_dag - 1),
            "unaccounted_flux_symbolic": str(residual),
            "unaccounted_equals_J_bare_times_Phi_minus_1":
                sp.simplify(residual - expected_phantom) == 0,
            "donor_pool_found": donor_found}


def run():
    res = {"G2.1": g21_conservation(), "G2.2": g22_dimensions(),
           "G2.3": g23_invariance(), "G2.4": g24_numeric(), "G2.5": g25_phi_phantom()}
    res["pass"] = {
        "G2.1": all(res["G2.1"].values()),
        "G2.2": res["G2.2"]["n_bad_terms"] == 0,
        "G2.3": all(res["G2.3"].values()),
        "G2.4": (res["G2.4"]["n_integration_failures"] == 0
                 and res["G2.4"]["max_rel_err_P_T"] < 1e-6
                 and res["G2.4"]["max_rel_err_C_T"] < 1e-6
                 and res["G2.4"]["most_negative_state"] > -1e-9),
        "G2.5": (not res["G2.5"]["donor_pool_found"]
                 and res["G2.5"]["unaccounted_equals_J_bare_times_Phi_minus_1"]),
    }
    return res


if __name__ == "__main__":
    r = run()
    print(json.dumps(r, indent=2, default=str))
    (HERE / "proofs.json").write_text(json.dumps(r, indent=2, default=str))
    sys.exit(0 if all(r["pass"].values()) else 1)
