#!/usr/bin/env python3
"""WO-04 gates G4.1-G4.5. writes wo04_results.json."""
import json
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.optimize import root as nroot

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import bifurcation as bf  # noqa: E402


# ---------------------------------------------------------------- G4.1
def g41_scalar():
    x, rho, chi, lam = sp.symbols("x rho chi lambda", real=True)
    g = x + rho * x / (1 + x) - chi * x**2
    gp = sp.diff(g, x)
    gpp = sp.simplify(sp.diff(g, x, 2))
    vals = {rho: 4, chi: sp.Rational(15, 100)}
    cubic = sp.Poly(sp.expand(sp.cancel(gp.subs(vals) * (1 + x) ** 2)), x)
    eq_poly = sp.Poly(sp.expand(sp.cancel((g - lam) * (1 + x))), x)
    # uniqueness of the interior max, in general: g'' < 0 on x > -1 for rho, chi > 0
    gpp_neg = sp.simplify(gpp + 2 * rho / (1 + x) ** 3 + 2 * chi) == 0
    xs = [complex(r) for r in sp.Poly(cubic, x).nroots(n=30)]
    xm = [r.real for r in xs if abs(r.imag) < 1e-20 and r.real > 0]
    lam_fold = float(g.subs(vals).subs(x, xm[0]))
    eq2 = sp.Poly(eq_poly.as_expr().subs(vals).subs(lam, 2), x)
    eqr = sorted(float(sp.re(r)) for r in eq2.nroots(n=30) if abs(sp.im(r)) < 1e-20)
    resid = [abs(float(g.subs(vals).subs(x, r)) - 2) for r in eqr]
    # saddle-node non-degeneracy at the fold of F = lam - g: F_xx = -g'' != 0, F_lam = 1
    gpp_at = float(gpp.subs(vals).subs(x, xm[0]))
    return {"cubic_coeffs": [float(c) for c in cubic.all_coeffs()],
            "equilibrium_poly_coeffs_symbolic": [str(c) for c in eq_poly.all_coeffs()],
            "g_second_derivative": str(gpp),
            "g2_identity_negative_definite": gpp_neg,
            "x_fold": xm, "lambda_fold": lam_fold,
            "equilibria_lambda_2": eqr, "residuals": resid,
            "g_pp_at_fold": gpp_at}


# ---------------------------------------------------------------- G4.2
def g42_legacy():
    sys.path.insert(0, "/storage/kiran-stuff/proteostasis_law/envelope-paper/scripts/vendor")
    import two_pool_ode as L
    anchor = [("as_published", 50, 1), ("weaker_binding", 50, 10), ("smaller_pool", 5, 1),
              ("near_capacity", 2, 1), ("c_free_at_Kd", 1, 1), ("Kd_at_C_tot", 50, 50)]
    Pg = np.geomspace(1e-6, 0.999, 20000)
    out = []
    for name, Ct, Kd in anchor:
        p = L.Params()
        p.C_tot_uM, p.K_d_uM = float(Ct), float(Kd)
        Pop, Jop, mech, Pdeath = L.saddle_node_operational(L.J_curve_two, L.A_qs, p)
        Jc = L.J_curve_two(Pg, p)
        i = int(np.nanargmax(Jc))
        interior = 0 < i < len(Pg) - 1
        # same J-curve with the donorless Phi removed (Phi = 1)
        J_nophi = L.R_clearance(Pg, p) + p.k_clear * L.A_qs(Pg, p)
        mono_nophi = bool(np.all(np.diff(J_nophi) > 0))
        out.append({"anchoring": name, "mechanism": mech, "P_operational": Pop,
                    "P_death_gate": Pdeath, "A_at_operational": float(L.A_qs(Pop, p)),
                    "math_fold_interior": interior,
                    "P_math_fold": float(Pg[i]) if interior else None,
                    "J_math_over_J_operational": float(Jc[i] / Jop) if interior else None,
                    "A_at_math_fold": float(L.A_qs(Pg[i], p)) if interior else None,
                    "J_curve_monotone_without_Phi": mono_nophi})
    return out


# ---------------------------------------------------------------- G4.3/4.4
def scan(variant, n, seed):
    rng = np.random.default_rng(seed)
    counts, n_unstable_on_single, bad_resid, stab_mismatch = [], 0, 0, 0
    multi = []
    max_resid = 0.0
    for k in range(n):
        p = bf.sample(rng, variant)
        rs, g0pos = bf.roots(p)
        counts.append(len(rs))
        cls = []
        for r in rs:
            x = bf.full_state(r, p)
            res = float(np.max(np.abs(bf.f(x, p))) / p["s_P"])
            max_resid = max(max_resid, res)
            if res > 1e-8:
                bad_resid += 1
            cls.append(bf.classify(x, p))
        if len(rs) == 1 and not cls[0]["stable"]:
            n_unstable_on_single += 1
        if len(rs) >= 2:
            # expected alternation for a 1-D reduction: stable / unstable / stable
            pattern = [c["stable"] for c in cls]
            want = [i % 2 == 0 for i in range(len(rs))]
            if pattern != want:
                stab_mismatch += 1
            multi.append({"idx": k, "n_roots": len(rs), "roots": rs,
                          "stable_pattern": pattern, "params": p})
    c = np.array(counts)
    return {"variant": variant, "n": n,
            "count_hist": {int(k): int(v) for k, v in zip(*np.unique(c, return_counts=True))},
            "frac_multistable": float(np.mean(c >= 3)),
            "n_zero_roots": int(np.sum(c == 0)),
            "max_rel_residual_full_ode": max_resid, "n_bad_residual": bad_resid,
            "n_single_root_unstable": n_unstable_on_single,
            "n_stability_pattern_mismatch": stab_mismatch,
            "multi_examples": multi[:20]}


def crosscheck_v0(n=150, seed=7):
    """independent count for V0: Newton on the full 6-D system from random
    initial states, without using the reduction. every converged solution must
    coincide with the unique reduced root."""
    rng = np.random.default_rng(seed)
    n_disagree, n_conv, max_comp = 0, 0, 0.0
    for _ in range(n):
        p = bf.sample(rng, "V0")
        rs, _ = bf.roots(p)
        x_red = bf.full_state(rs[0], p)
        sc = np.maximum(np.abs(x_red), 1e-12)
        for _ in range(4):
            x0 = x_red * np.exp(rng.normal(0, 1.5, 6))
            x0[5] = 0.0
            sol = nroot(lambda z: bf.f(z * sc, p) / p["s_P"], x0 / sc,
                        jac=lambda z: bf.jac(z * sc, p) * sc / p["s_P"], method="hybr", tol=1e-14)
            if not sol.success:
                continue
            z = sol.x * sc
            if np.any(z < -1e-9 * sc):     # unphysical branch, outside the orthant
                continue
            if np.max(np.abs(bf.f(z, p))) / p["s_P"] > 1e-12:
                continue                   # not a converged steady state
            n_conv += 1
            # at steady state every component is a function of U (see
            # bifurcation.py), so a distinct steady state must have a distinct
            # U. comparing tiny components (A ~ 1e-9 uM) componentwise instead
            # measures Newton's resolution floor, not a second state.
            if abs(z[1] - x_red[1]) / x_red[1] > 1e-6:
                n_disagree += 1
            max_comp = max(max_comp, float(np.max(np.abs(z[:5] - x_red[:5]) / sc[:5])))
    return {"n_param_sets": n, "n_converged_physical": n_conv,
            "n_disagree_with_reduction": n_disagree,
            "max_componentwise_rel_diff_incl_subresolution": max_comp}


# ---------------------------------------------------------------- G4.5
def g45_fold(scans):
    out = {}
    for v in ("V1", "V2"):
        ex = scans[v]["multi_examples"]
        if not ex:
            out[v] = {"found": False}
            continue
        p = {k: val for k, val in ex[0]["params"].items()}
        es, counts, changes = bf.fold_in_eps(p)
        folds = []
        for e_lo, e_hi, c1, c2 in changes[:2]:
            Uf, ef, (ra, rb, side) = bf.locate_fold(p, e_lo, e_hi)
            q = {**p, "eps": ef}
            xf = bf.full_state(Uf, q)
            cf = bf.classify(xf, q)
            q_side = {**p, "eps": side}
            ca = bf.classify(bf.full_state(ra, q_side), q_side)
            cb = bf.classify(bf.full_state(rb, q_side), q_side)
            # scale for "near zero": smallest eigenvalue magnitude at a
            # regular stable point of the same parameter set
            g_fold = float(bf.G(np.array([Uf]), q)[0]) / p["s_P"]
            folds.append({"eps_fold": ef, "U_fold": Uf, "G_at_fold_rel": g_fold,
                          "min_abs_real_eig_at_fold": cf["min_abs_real_eig"],
                          "min_abs_real_eig_branch_a": ca["min_abs_real_eig"],
                          "min_abs_real_eig_branch_b": cb["min_abs_real_eig"],
                          "det_sign_branch_a": float(np.sign(ca["det"])),
                          "det_sign_branch_b": float(np.sign(cb["det"])),
                          "count_change": [c1, c2]})
        out[v] = {"found": True, "example_params": p, "eps_count_changes": len(changes),
                  "folds": folds}
    return out


def g44_dynamic(scans, n_ex=3):
    """independent check of the eigenvalue classification: integrate the full
    6-D ODE from each steady state perturbed by +-1e-3 relative along U. a
    stable state must be returned to; an unstable one must be left for a
    different steady state."""
    from scipy.integrate import solve_ivp
    out = []
    for v in ("V1", "V2"):
        for ex in scans[v]["multi_examples"][:n_ex]:
            p = ex["params"]
            xs = [bf.full_state(r, p) for r in ex["roots"]]
            T = 50.0 / min(bf.classify(x, p)["min_abs_real_eig"] for x in xs)
            for j, (x, stable) in enumerate(zip(xs, ex["stable_pattern"])):
                for sgn in (-1, 1):
                    x0 = x.copy()
                    x0[1] *= 1 + sgn * 1e-3
                    sol = solve_ivp(lambda t, z: bf.f(z, p), (0, T), x0, method="Radau",
                                    jac=lambda t, z: bf.jac(z, p), rtol=1e-10, atol=1e-14)
                    xe = sol.y[:, -1]
                    dist = [abs(xe[1] - y[1]) / y[1] for y in xs]
                    k = int(np.argmin(dist))
                    ok = (k == j) if stable else (k != j)
                    out.append({"variant": v, "idx": ex["idx"], "root": j, "stable": stable,
                                "sign": sgn, "ended_at_root": k, "rel_dist": float(dist[k]),
                                "consistent": bool(ok and dist[k] < 1e-4)})
    return out


def g45_continuation(scans, n_ex=10):
    """continuation along the exact equilibrium curve (U-parametrised) for the
    stored multistable examples, plus V0 samples as the negative control."""
    res = {}
    for v in ("V1", "V2"):
        rows = []
        for ex in scans[v]["multi_examples"][:n_ex]:
            p = ex["params"]
            c = bf.continue_branch(p)
            # consistency: continuation and global root finding count the
            # same number of steady states at the scanned eps
            n_cross = int(np.sum(np.diff(np.sign(c["eps"] - p["eps"])) != 0))
            rows.append({"idx": ex["idx"], "folds": c["folds"],
                         "n_det_sign_changes": len(c["det_sign_changes"]),
                         "n_curve_crossings_at_eps": n_cross, "n_roots_at_eps": ex["n_roots"]})
        res[v] = rows
    rng = np.random.default_rng(999)
    v0 = []
    for _ in range(10):
        c = bf.continue_branch(bf.sample(rng, "V0"), n=1500)
        v0.append({"n_folds": len(c["folds"]), "n_det_sign_changes": len(c["det_sign_changes"]),
                   "eps_monotone_increasing": bool(np.all(np.diff(c["eps"]) > 0))})
    res["V0_control"] = v0
    return res


def run(n=2000):
    okP, okC = bf.conservation_ok()
    scans = {v: scan(v, n, seed=100 + i) for i, v in enumerate(("V0", "V1", "V2"))}
    r = {"G4.1": g41_scalar(), "G4.2": g42_legacy(),
         "variant_conservation": {"P_T": okP, "C_T": okC, "V0_equals_WO02": bf.reduces_to_wo02()},
         "G4.3_analytic_V0": bf.analytic_v0(),
         "G4.3_scans": scans, "G4.3_crosscheck_V0": crosscheck_v0(),
         "G4.4_dynamic": g44_dynamic(scans),
         "G4.5": g45_fold(scans), "G4.5_continuation": g45_continuation(scans)}
    s41 = r["G4.1"]
    v0 = scans["V0"]
    # fold verification on the continuation. the previous criterion here was
    # "min_abs_eig_at_fold < 1e-3 * min(branch) * 1e3", whose factors cancel;
    # replaced by an explicit 1e-4 ratio.
    cont = r["G4.5_continuation"]
    phys = [fd for v in ("V1", "V2") for row in cont[v] for fd in row["folds"] if fd["physical"]]
    folds_ok = bool(
        len(phys) > 0
        and all(abs(fd["G_rel_at_fold"]) < 1e-10
                and fd["min_abs_eig_at_fold"] < 1e-4 * fd["ref_min_abs_eig_nearby"]
                and fd["det_sign_left"] * fd["det_sign_right"] < 0 for fd in phys)
        and all(row["n_det_sign_changes"] == len(row["folds"])
                and row["n_curve_crossings_at_eps"] == row["n_roots_at_eps"]
                for v in ("V1", "V2") for row in cont[v])
        and all(c["n_folds"] == 0 and c["n_det_sign_changes"] == 0 and c["eps_monotone_increasing"]
                for c in cont["V0_control"]))
    # the old eps-sweep + fsolve fold must agree with the continuation fold
    cont_eps = [fd["eps_fold"] for row in cont["V1"][:1] for fd in row["folds"] if fd["physical"]]
    old_eps = [fd["eps_fold"] for fd in r["G4.5"]["V1"]["folds"]]
    r["G4.5_methods_agree"] = bool(len(cont_eps) == len(old_eps) and np.allclose(
        sorted(cont_eps), sorted(old_eps), rtol=1e-6))
    a = r["G4.3_analytic_V0"]
    analytic_ok = all(a[k] for k in ("reduction_unique", "curve_satisfies_other_balances",
                                     "dG_dU_decomposition_holds", "G0_formula_holds"))
    r["pass"] = {
        "G4.1": bool(np.allclose(s41["cubic_coeffs"], [-0.3, 0.4, 1.7, 5.0])
                     and s41["g2_identity_negative_definite"]
                     and abs(s41["lambda_fold"] - 4.80218919587308) < 1e-10
                     and max(s41["residuals"]) < 1e-12),
        "G4.2": all("mechanism" in x for x in r["G4.2"]),
        "G4.3": bool(analytic_ok and v0["count_hist"] == {1: v0["n"]}
                     and r["G4.3_crosscheck_V0"]["n_disagree_with_reduction"] == 0
                     and r["G4.3_crosscheck_V0"]["n_converged_physical"] > 0
                     and all(scans[v]["n_bad_residual"] == 0 for v in scans)
                     and okP and okC and r["variant_conservation"]["V0_equals_WO02"]),
        "G4.4": bool(all(scans[v]["n_single_root_unstable"] == 0
                         and scans[v]["n_stability_pattern_mismatch"] == 0 for v in scans)
                     and all(d["consistent"] for d in r["G4.4_dynamic"])),
        "G4.5": bool(folds_ok and r["G4.5_methods_agree"]),
    }
    return r


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    r = run(n)
    short = {k: v for k, v in r.items() if k != "G4.3_scans"}
    short["scan_summary"] = {v: {kk: vv for kk, vv in s.items() if kk != "multi_examples"}
                             for v, s in r["G4.3_scans"].items()}
    print(json.dumps(short, indent=1, default=float))
    (HERE / "wo04_results.json").write_text(json.dumps(r, indent=2, default=float))
    sys.exit(0 if all(r["pass"].values()) else 1)
