"""
permanent checks: bookkeeping and conservation for stages 0, 1, 2 and 4.

  C0  every flux is internal (column sum 0 on every total) or crosses the
      boundary exactly as the BOUNDARY table declares (stage 0 rule)
  C1  every rhs term has units uM/s (dimension exponents, stage 1, 2, 4)
  C2  totals: dX_T/dt equals source - boundary sinks - mu X_T symbolically,
      for the client total and each machine total, in every stage
  C3  nonnegative orthant: on each face x_i = 0 the rhs of x_i has only
      nonnegative terms (eps, phi written as ratios of nonnegatives)
  C4  term-for-term nesting against the source work orders (read-only import):
      S2C == WO-02 model.RHS_SYM, S4C == WO-04 V1 RHS, S4C(k_onA=0) == V0
      and S4 -> S2 at k_onA = 0
  C5  term ledger: every flux and dilution term has one row, and each row's
      donor and receiver are exactly the pools the stoichiometry debits and credits
  C6  numerical: random integrations of S2 and S4 conserve all totals
      (rel. < 1e-6 against closed forms) and no pool goes negative
"""
import csv
import sys
from pathlib import Path

import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import model as m  # noqa: E402

REBUILD = HERE.parents[1] / "proteostasis_rebuild"
LEDGER = HERE.parent / "EQUATION_TERM_LEDGER.tsv"


def load_source(wo, stem, alias):
    """import a source work-order module read-only under a unique name (WO-02's
    file is also called model.py). bytecode writing is off, so nothing is written."""
    import importlib.util
    if alias in sys.modules:
        return sys.modules[alias]
    sys.path.insert(0, str(REBUILD / "WO-01"))
    spec = importlib.util.spec_from_file_location(alias, REBUILD / wo / f"{stem}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[alias] = mod
    spec.loader.exec_module(mod)
    return mod


def c0_flux_bookkeeping():
    bad = []
    for f, (_, st, _) in m.FLUX.items():
        for kind, w in m.WEIGHTS.items():
            net = sum(c * w.get(p, 0) for p, c in st.items())
            if net != m.BOUNDARY[kind].get(f, 0):
                bad.append((f, kind, net))
    return {"pass": not bad, "violations": bad, "n_fluxes": len(m.FLUX)}


def c1_units():
    bad = []
    Lc, Tt = sp.symbols("Lc Tt", positive=True)
    scale = {m.XS[n]: m.XS[n] * Lc for n in m.STATE_NAMES}
    for n in m.PARAM_NAMES:
        ce, te = m.DIMS[n]
        scale[m.PS[n]] = m.PS[n] * Lc**ce * Tt**te
    n_terms = 0
    for st in m.STAGES:
        for pool, ex in m.assemble(st).items():
            for term in sp.Add.make_args(sp.expand(ex)):
                n_terms += 1
                r = sp.powsimp(sp.expand(term.subs(scale, simultaneous=True) / (Lc / Tt)))
                if r.has(Lc) or r.has(Tt):
                    bad.append((st, pool, str(term)))
    return {"pass": not bad, "violations": bad, "n_terms_checked": n_terms}


def c2_totals():
    out = {}
    for st in m.STAGES:
        rhs = m.assemble(st)
        for kind, w in m.WEIGHTS.items():
            pools = [p for p in w if p in rhs]
            if not pools:
                continue
            lhs = sum(w[p] * rhs[p] for p in pools)
            out[f"{st}:{kind}"] = sp.simplify(lhs - m.expected_total_rate(kind, st)) == 0
    return {"pass": all(out.values()), "identities": out}


def c3_faces():
    a, b, c, d = sp.symbols("a b c d", positive=True)
    ratio = {m.eps: a / (a + b), m.phi: c / (c + d)}
    bad = []
    for st in m.STAGES:
        for pool, ex in m.assemble(st).items():
            face = sp.together(ex.subs(m.XS[pool], 0).subs(ratio))
            num = sp.expand(sp.numer(face))
            if num == 0:
                continue
            coeffs = sp.Poly(num, *num.free_symbols).coeffs()
            if any(cf < 0 for cf in coeffs):
                bad.append((st, pool))
    return {"pass": not bad, "violations": bad}


def c4_nesting():
    """compare assembled stages with the source work orders' own sympy rhs."""
    res = {}
    wo2 = load_source("WO-02", "model", "wo02_model")
    bf = load_source("WO-04", "bifurcation", "wo04_bifurcation")

    def sub_names(expr, mod_syms):
        return expr.subs({s: (m.XS.get(s.name) or m.PS[s.name]) for s in expr.free_symbols
                          if s.name in m.XS or s.name in m.PS}, simultaneous=True)

    s2c = m.assemble("S2C")
    order2 = wo2.STATE_NAMES
    res["S2C_equals_WO02"] = all(
        sp.simplify(s2c[n] - sub_names(wo2.RHS_SYM[i], None)) == 0 for i, n in enumerate(order2))
    s4c = m.assemble("S4C")
    order4 = bf.NAMES
    v1 = {bf.k_dcat: 0}
    res["S4C_equals_WO04_V1"] = all(
        sp.simplify(s4c[n] - sub_names(bf.RHS[i].subs(v1), None)) == 0
        for i, n in enumerate(order4))
    v0 = {bf.k_dcat: 0, bf.k_onA: 0}
    s4c0 = {k: v.subs({m.k_onA: 0}) for k, v in s4c.items()}
    res["S4C_kon0_equals_WO04_V0"] = all(
        sp.simplify(s4c0[n] - sub_names(bf.RHS[i].subs(v0), None)) == 0
        for i, n in enumerate(order4))
    # S4 -> S2 exactly when the binding step is switched off and CA = 0
    s4 = m.assemble("S4")
    s2 = m.assemble("S2")
    res["S4_kon0_CA0_equals_S2"] = all(
        sp.simplify(s4[n].subs({m.k_onA: 0, m.CA: 0}) - s2[n]) == 0 for n in s2) and \
        sp.simplify(s4["CA"].subs({m.k_onA: 0, m.CA: 0})) == 0
    return {"pass": all(res.values()), **res}


def donors_receivers(f):
    rate, st, _ = m.FLUX[f]
    don = sorted(p for p, c in st.items() if c < 0)
    rec = sorted(p for p, c in st.items() if c > 0)
    if m.BOUNDARY["client"].get(f) == 1 or f.startswith("syn_"):
        don = ["boundary: synthesis"]
    if m.BOUNDARY["client"].get(f) == -1:
        rec = sorted(rec + ["boundary: degradation"])
    return don, rec


def c5_ledger():
    rows = list(csv.DictReader(open(LEDGER), delimiter="\t"))
    by_id = {r["term_id"]: r for r in rows}
    missing, wrong = [], []
    for f in m.FLUX:
        r = by_id.get(f)
        if r is None:
            missing.append(f)
            continue
        don, rec = donors_receivers(f)
        got_d = sorted(x.strip() for x in r["donor_pool"].split("+"))
        got_r = sorted(x.strip() for x in r["receiver_pool"].split("+"))
        if got_d != don or got_r != rec:
            wrong.append((f, got_d, don, got_r, rec))
        if int(r["stage_introduced"]) != m.FLUX[f][2]:
            wrong.append((f, "stage", r["stage_introduced"], m.FLUX[f][2]))
    for n in m.STATE_NAMES:
        r = by_id.get("dil_" + n)
        if r is None:
            missing.append("dil_" + n)
        elif "dilution" not in r["receiver_pool"]:
            wrong.append(("dil_" + n, r["receiver_pool"]))
    blank = [r["term_id"] for r in rows if r["kind"] in ("flux", "dilution")
             and (not r["donor_pool"].strip() or not r["receiver_pool"].strip()
                  or not r["rate_law_justification"].strip()
                  or not r["parameter_status_WO07"].strip())]
    banned = [r["term_id"] for r in rows if r["kind"] in ("flux", "dilution")
              and r["rate_law_justification"].strip().lower() in ("phenomenological", "effective", "standard")]
    ok = not (missing or wrong or blank or banned)
    return {"pass": ok, "n_rows": len(rows), "missing": missing, "wrong": wrong,
            "blank_required_field": blank, "banned_justification": banned}


def c6_integrations(n=24, seed=11):
    rng = np.random.default_rng(seed)
    worst = {"client": 0.0, "machine": 0.0, "min_state_rel": 0.0}
    fails = 0
    for stage in ("S2", "S4"):
        names, f, J = m.lambdas(stage)
        idx = {nm: i for i, nm in enumerate(names)}
        for _ in range(n):
            p = m.sample_stage2(rng)
            if stage == "S4":
                p["k_onA"], p["k_offA"] = m.lu(rng, 1e-3, 10), m.lu(rng, 1e-3, 10)
            pv = m.pvec(p)
            x0 = np.zeros(len(names))
            # start away from steady state: all client native at half balance,
            # machines free at twice their balanced totals
            x0[idx["N"]] = 0.5 * p["P_T"]
            for mname, tot in (("C", "C_T"), ("D", "D_T"), ("E", "E_T")):
                x0[idx[mname]] = 2 * p[tot]
            # augmented variable Q integrates the CLAIMED client balance
            wc = np.array([m.WEIGHTS["client"].get(nm, 0) for nm in names], float)
            kdA, kcD, mu_ = p["k_dA"], p["k_catD"], p["mu"]

            def rhs(t, y):
                x = y[:-1]
                dq = p["s_P"] - kcD * x[idx["DU"]] - kdA * x[idx["A"]] - mu_ * y[-1]
                return np.append(np.asarray(f(x, pv), float), dq)

            def jac(t, y):
                Jx = np.asarray(J(y[:-1], pv), float)
                Jq = np.zeros(len(names))
                Jq[idx["DU"]], Jq[idx["A"]] = -kcD, -kdA
                top = np.hstack([Jx, np.zeros((len(names), 1))])
                return np.vstack([top, np.append(Jq, -mu_)])

            y0 = np.append(x0, wc @ x0)
            T = 3.0 / mu_
            sol = solve_ivp(rhs, (0, T), y0, method="LSODA", jac=jac, rtol=1e-10,
                            atol=1e-12 * p["P_T"], dense_output=False)
            if not sol.success:
                fails += 1
                continue
            X = sol.y[:-1]
            Pt = wc @ X
            ec = float(np.max(np.abs(Pt - sol.y[-1]) / np.maximum(np.abs(sol.y[-1]), 1e-300)))
            em = 0.0
            for kind, tot in (("chaperone", "C_T"), ("protease", "D_T"), ("disaggregase", "E_T")):
                w = np.array([m.WEIGHTS[kind].get(nm, 0) for nm in names], float)
                exact = p[tot] + (2 * p[tot] - p[tot]) * np.exp(-mu_ * sol.t)
                em = max(em, float(np.max(np.abs(w @ X - exact) / exact)))
            worst["client"] = max(worst["client"], ec)
            worst["machine"] = max(worst["machine"], em)
            worst["min_state_rel"] = min(worst["min_state_rel"], float(X.min() / p["P_T"]))
    ok = fails == 0 and worst["client"] < 1e-6 and worst["machine"] < 1e-6 and worst["min_state_rel"] > -1e-9
    return {"pass": ok, "n_per_stage": n, "solver_failures": fails, **worst,
            "solver": "LSODA rtol 1e-10, atol 1e-12 P_T, horizon 3/mu"}


def c7_negative_controls():
    """each check must FAIL on a planted banned construction (architecture 0.1).
    the flux table is patched temporarily and restored in finally."""
    import copy
    saved = (copy.copy(m.FLUX), copy.deepcopy(m.STAGES))
    out = {}
    try:
        # (a) legacy donorless amplification: extra inflow into U with no donor
        m.FLUX["phi_amp"] = (sp.Symbol("J_bare") * (sp.Symbol("Phi") - 1), {"U": 1}, 1)
        m.STAGES["S2"] = m.STAGES["S2"] + ["phi_amp"]
        out["donorless_source_caught_C0"] = not c0_flux_bookkeeping()["pass"]
        out["donorless_source_caught_C2"] = not c2_totals()["pass"]
        m.FLUX.pop("phi_amp")
        m.STAGES["S2"] = saved[1]["S2"]
        # (b) +chi U^2 SOURCE of non-native protein (aggregation written as a source)
        m.FLUX["chiU2"] = (sp.Symbol("chi") * m.U**2, {"U": 1}, 1)
        m.STAGES["S1"] = m.STAGES["S1"] + ["chiU2"]
        out["chiU2_source_caught_C0"] = not c0_flux_bookkeeping()["pass"]
        m.FLUX.pop("chiU2")
        m.STAGES["S1"] = saved[1]["S1"]
        # (c) unit error: aggregation written k_a U^3
        m.FLUX["agg"] = (m.k_a * m.U**3, {"U": -1, "A": 1}, 1)
        out["unit_error_caught_C1"] = not c1_units()["pass"]
        m.FLUX["agg"] = saved[0]["agg"]
        # (d) a term that differs from WO-02 must break the nesting identity
        m.FLUX["rel"] = (2 * m.k_off * m.B, {"B": -1, "C": 1, "U": 1}, 2)
        out["perturbed_term_caught_C4"] = not c4_nesting()["S2C_equals_WO02"]
        m.FLUX["rel"] = saved[0]["rel"]
        # (e) a sink that does not vanish on its own face (degradation of U at U = 0)
        m.FLUX["degU"] = (m.k_d * (m.U + 1), {"U": -1}, 1)
        out["face_violation_caught_C3"] = not c3_faces()["pass"]
    finally:
        m.FLUX.clear()
        m.FLUX.update(saved[0])
        m.STAGES.clear()
        m.STAGES.update(saved[1])
    out["restored"] = c0_flux_bookkeeping()["pass"] and c4_nesting()["pass"]
    return {"pass": all(out.values()), **out}


def run():
    return {"C0_flux_bookkeeping": c0_flux_bookkeeping(), "C1_units": c1_units(),
            "C2_totals": c2_totals(), "C3_faces": c3_faces(), "C4_nesting": c4_nesting(),
            "C5_ledger": c5_ledger(), "C6_integrations": c6_integrations(),
            "C7_negative_controls": c7_negative_controls()}


if __name__ == "__main__":
    import json
    r = run()
    print(json.dumps({k: v["pass"] for k, v in r.items()}, indent=1))
    sys.exit(0 if all(v["pass"] for v in r.values()) else 1)
