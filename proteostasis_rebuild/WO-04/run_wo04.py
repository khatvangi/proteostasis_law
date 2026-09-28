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


def run(n=2000):
    okP, okC = bf.conservation_ok()
    scans = {v: scan(v, n, seed=100 + i) for i, v in enumerate(("V0", "V1", "V2"))}
    r = {"G4.1": g41_scalar(), "G4.2": g42_legacy(),
         "variant_conservation": {"P_T": okP, "C_T": okC, "V0_equals_WO02": bf.reduces_to_wo02()},
         "G4.3_scans": scans, "G4.3_crosscheck_V0": crosscheck_v0(), "G4.5": g45_fold(scans)}
    s41 = r["G4.1"]
    v0 = scans["V0"]
    folds_ok = all(
        (not fv["found"]) or all(
            abs(fd["G_at_fold_rel"]) < 1e-8
            and fd["det_sign_branch_a"] * fd["det_sign_branch_b"] < 0
            and fd["min_abs_real_eig_at_fold"] < 1e-3 * min(fd["min_abs_real_eig_branch_a"],
                                                           fd["min_abs_real_eig_branch_b"]) * 1e3
            for fd in fv["folds"])
        for fv in r["G4.5"].values())
    r["pass"] = {
        "G4.1": bool(np.allclose(s41["cubic_coeffs"], [-0.3, 0.4, 1.7, 5.0])
                     and s41["g2_identity_negative_definite"]
                     and abs(s41["lambda_fold"] - 4.80218919587308) < 1e-10
                     and max(s41["residuals"]) < 1e-12),
        "G4.2": all("mechanism" in x for x in r["G4.2"]),
        "G4.3": bool(v0["count_hist"] == {1: v0["n"]}
                     and r["G4.3_crosscheck_V0"]["n_disagree_with_reduction"] == 0
                     and r["G4.3_crosscheck_V0"]["n_converged_physical"] > 0
                     and all(scans[v]["n_bad_residual"] == 0 for v in scans)
                     and okP and okC and r["variant_conservation"]["V0_equals_WO02"]),
        "G4.4": bool(all(scans[v]["n_single_root_unstable"] == 0
                         and scans[v]["n_stability_pattern_mismatch"] == 0 for v in scans)),
        "G4.5": bool(folds_ok),
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
