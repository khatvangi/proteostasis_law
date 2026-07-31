#!/usr/bin/env python
"""experiment A -- baseline sanity and invariant-domain tests.

this is the falsification harness pointed at the model's own internals. it does
not test the theory against biology; it tests whether the numerical object we
are about to sweep in experiments B-D is the object the theory specifies.

checks, each with an explicit tolerance and a pass/fail verdict:

  A1  every rate law and both rhs terms carry concentration time^-1, and a
      deliberately mis-dimensioned rate constant is REJECTED (negative control)
  A2  the nondimensionalization is an exact rescaling of the dimensional model
  A3  the rapid-equilibrium free-pool solution is unique (knaster-tarski
      bracket) at every sampled state
  A4  the algebraic closure reproduces the totals it was built from
  A5  d(u+a)/dt equals influx minus removal exactly (transfer stoichiometry)
  A6  analytic jacobian matches central differences
  A7  the nonnegative orthant is forward invariant under integration
  A8  j above the analytic removal ceiling admits no bounded state
  A9  blind multi-start and warm-started continuation find the same equilibria
  A10 a repeated deterministic slice reproduces bit-for-bit

usage:  python scripts/run_experiment_a.py --config configs/phase1/experiment_a.json
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from proteostasis import (Params, ModelError, jacobian, numericalJacobian, fluxes,
                          solveFreePoolsCertified, poolBalanceResiduals,
                          massBalanceResidual, removalCeiling, findEquilibria,
                          solveEquilibriumFrom, traceBranch, simulate,
                          defaultInitialConditions, allocationParams)
from proteostasis import units as U
from proteostasis.provenance import (loadConfig, writeTable, writeProvenance,
                                     hashObject, canonicalJson)
from proteostasis.sweeps import latinHypercube, paramsFromSample, rangesFromConfig

TOL = dict(units=0.0, rescale=1e-6, unique=1e-10, closure=1e-10,
           massbal=1e-11, jacobian=1e-5, positivity=0.0, ceiling=0.0,
           crossvalidate=1e-4, determinism=0.0)


def _verdict(name: str, metric: float, tol: float, detail: str = "") -> dict:
    ok = bool(np.isfinite(metric) and metric <= tol)
    return dict(check=name, metric=float(metric), tolerance=float(tol),
                passed=ok, detail=detail)


# --------------------------------------------------------------------------- A1
def checkDimensions() -> list:
    """positive control: the real model is dimensionally consistent.
    negative control: a wrong-dimension rate constant must RAISE."""
    dp = U.DimensionalParams()
    q = lambda v: U.Quantity(v, U.CONC)
    U.dimensionalRhs(q(0.1), q(0.05), q(0.3), q(0.3), dp)
    U.dimensionalConservationCheck(q(0.1), q(0.05), q(0.3), q(0.3), dp)
    rows = [_verdict("A1_dimensional_consistency", 0.0, TOL["units"],
                     "all rate laws and both rhs terms are amount volume^-1 time^-1")]

    # negative control: declare the nucleation coefficient as a plain rate
    caught = False
    try:
        bad = U.Quantity(1e-3, U.RATE) * (q(0.1) ** dp.m)     # wrong dimension
        (U.Quantity(dp.J, U.FLUX) - bad)
    except U.DimensionError:
        caught = True
    rows.append(_verdict("A1_negative_control_rejects_bad_units",
                         0.0 if caught else 1.0, TOL["units"],
                         "mis-dimensioned k_n is rejected by term addition"))
    return rows


# --------------------------------------------------------------------------- A2
def checkNondimensionalization(t_end_dimensional: float = 2.0e4, n: int = 60) -> list:
    """integrate the dimensional and nondimensional models and compare.

    if any scaling exponent in `toNondimensional` were wrong, the two
    trajectories would diverge; a units-only inspection could not catch that.

    both systems are evaluated at MATCHED times (t and t/tau) rather than
    interpolated onto a common grid -- interpolation error would otherwise
    dominate the comparison and be mistaken for a scaling error. the two sides
    also use different free-pool solvers (bounded least squares on the full
    4-unknown system vs 2-d newton on the eliminated system), so agreement is
    evidence about the algebra rather than about one solver reproducing itself.
    """
    from scipy.integrate import solve_ivp
    from proteostasis.model import rhsVector

    dp = U.DimensionalParams()
    p, phi, tau = U.toNondimensional(dp)
    U0, A0 = 0.05 * phi, 0.01 * phi
    t_dim = np.linspace(0.0, t_end_dimensional, n)

    solD = solve_ivp(lambda t, y: U.dimensionalRhsNumeric(max(y[0], 0.0), max(y[1], 0.0), dp),
                     (0.0, t_end_dimensional), [U0, A0], t_eval=t_dim,
                     method="LSODA", rtol=1e-11, atol=1e-16)
    solN = solve_ivp(lambda t, y: rhsVector((max(y[0], 0.0), max(y[1], 0.0)), p),
                     (0.0, t_end_dimensional / tau), [U0 / phi, A0 / phi],
                     t_eval=t_dim / tau, method="LSODA", rtol=1e-11, atol=1e-16)

    scale = max(float(np.max(np.abs(solD.y))), 1e-12)
    err = float(np.max(np.abs(solN.y * phi - solD.y)) / scale)
    return [_verdict("A2_nondimensionalization_is_exact_rescaling", err, TOL["rescale"],
                     f"phi={phi:.6g}, tau={tau:.6g}, {n} matched time points")]


# --------------------------------------------------------------------------- A3-A6
def checkStateGrid(param_sets: list, states: np.ndarray) -> tuple:
    """free-pool uniqueness, closure, mass balance and jacobian over a grid."""
    rows = []
    worst = dict(unique_gap=0.0, closure=0.0, massbal=0.0, jac=0.0)
    for k, p in enumerate(param_sets):
        for (u, a) in states:
            _, _, _, _, cert = solveFreePoolsCertified(float(u), float(a), p)
            worst["unique_gap"] = max(worst["unique_gap"], cert["gap"])
            res = poolBalanceResiduals(float(u), float(a), p)
            worst["closure"] = max(worst["closure"], max(abs(v) for v in res.values()))
            mb = abs(massBalanceResidual(float(u), float(a), p))   # already relative
            worst["massbal"] = max(worst["massbal"], mb)
            J, Jn = jacobian(float(u), float(a), p), numericalJacobian(float(u), float(a), p)
            denom = max(float(np.max(np.abs(Jn))), 1e-12)
            worst["jac"] = max(worst["jac"], float(np.max(np.abs(J - Jn))) / denom)
            rows.append(dict(param_set=k, u=float(u), a=float(a),
                             unique_gap=cert["gap"], unique=cert["unique"],
                             closure_max=max(abs(v) for v in res.values()),
                             massbal_rel=mb,
                             jac_rel_err=float(np.max(np.abs(J - Jn))) / denom))
    verdicts = [
        _verdict("A3_free_pool_solution_unique", worst["unique_gap"], TOL["unique"],
                 "least and greatest fixed points of the monotone map coincide"),
        _verdict("A4_pool_closure_mass_balance", worst["closure"], TOL["closure"],
                 "reconstructed totals equal c_tot, d_tot, u, a"),
        _verdict("A5_ode_mass_balance", worst["massbal"], TOL["massbal"],
                 "d(u+a)/dt == influx - refold - degrade_u - degrade_a"),
        _verdict("A6_analytic_jacobian", worst["jac"], TOL["jacobian"],
                 "implicit-function jacobian vs central differences"),
    ]
    return verdicts, pd.DataFrame(rows)


# --------------------------------------------------------------------------- A7
def checkPositivity(param_sets: list, seed: int, t_end: float) -> tuple:
    rows, worst_leak, n_neg = [], 0.0, 0
    for k, p in enumerate(param_sets):
        for (u0, a0) in defaultInitialConditions(scale=1.0, n=12, seed=seed + k):
            tr = simulate(p, u0, a0, t_end=t_end, n_out=200, blowup=1.0e4)
            leak = max(-min(tr.min_u, tr.min_a), 0.0)
            worst_leak = max(worst_leak, leak)
            n_neg += int(tr.status == "negative")
            row = dict(param_set=k, u0=u0, a0=a0, leak=leak)
            row.update(tr.summary())
            rows.append(row)
    return ([_verdict("A7_nonnegative_orthant_invariant", worst_leak, TOL["positivity"],
                      f"{len(rows)} trajectories, {n_neg} negative-event terminations")],
            pd.DataFrame(rows))


# --------------------------------------------------------------------------- A8
def checkRemovalCeiling(param_sets: list, t_end: float) -> tuple:
    """j above the analytic ceiling must produce unbounded growth."""
    rows, violations = [], 0
    for k, p in enumerate(param_sets):
        ceiling = removalCeiling(p)
        pv = p.with_(j=1.5 * ceiling)
        tr = simulate(pv, 0.0, 0.0, t_end=t_end, n_out=200, blowup=1.0e4)
        eqs = findEquilibria(pv, n_grid=7)
        bad = (len(eqs) > 0) or (tr.status not in ("blowup",))
        violations += int(bad)
        rows.append(dict(param_set=k, ceiling=ceiling, j_tested=pv.j,
                         n_equilibria=len(eqs), status=tr.status,
                         final_burden=tr.final_u + tr.final_a, violation=bad))
    return ([_verdict("A8_no_bounded_state_above_removal_ceiling", float(violations),
                      TOL["ceiling"], "j > c_tot + (rho_U+rho_A) d_tot implies escape")],
            pd.DataFrame(rows))


# --------------------------------------------------------------------------- A9
def checkEquilibriumCrossValidation(base: Params, j_values: np.ndarray) -> tuple:
    """the blind grid finder and warm-started continuation must agree."""
    rows, worst = [], 0.0
    branch = traceBranch(base, "j", j_values)
    for bp in branch:
        p = base.with_(j=bp.value)
        blind = [e for e in findEquilibria(p, n_grid=9) if e.stable]
        cont = bp.eq
        if not blind and cont is None:
            rows.append(dict(j=bp.value, blind_n=0, agree=True, rel_diff=0.0))
            continue
        if bool(blind) != bool(cont is not None):
            rows.append(dict(j=bp.value, blind_n=len(blind), agree=False, rel_diff=np.inf))
            worst = np.inf
            continue
        b = blind[0]
        rel = max(abs(b.u - cont.u) / max(b.u, 1e-30),
                  abs(b.a - cont.a) / max(b.a, 1e-30))
        worst = max(worst, rel)
        rows.append(dict(j=bp.value, blind_n=len(blind), agree=rel <= TOL["crossvalidate"],
                         rel_diff=rel, u_blind=b.u, u_cont=cont.u,
                         a_blind=b.a, a_cont=cont.a))
    return ([_verdict("A9_blind_and_continuation_agree", worst, TOL["crossvalidate"],
                      f"{len(rows)} influx values compared")], pd.DataFrame(rows))


# --------------------------------------------------------------------------- A10
def checkDeterminism(base: Params, seed: int, ranges) -> list:
    """the same seed must reproduce the same sample and the same results."""
    def slice_():
        samples = latinHypercube(ranges, 12, seed=seed)
        out = []
        for s in samples:
            p = paramsFromSample(s, base)
            eq = solveEquilibriumFrom(1e-3, 1e-5, p.with_(j=1e-3))
            out.append(dict(sample=s, eq=None if eq is None else eq.asDict()))
        return out
    h1, h2 = hashObject(slice_()), hashObject(slice_())
    return [_verdict("A10_deterministic_rerun", 0.0 if h1 == h2 else 1.0,
                     TOL["determinism"], f"sha256={h1[:16]}")]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", default=None)
    args = ap.parse_args()

    cfg = loadConfig(args.config)
    outdir = Path(args.outdir or cfg["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    base = Params(**cfg.get("base_params", {})).validate()
    ranges = rangesFromConfig(cfg.get("param_ranges"))
    param_sets = [base] + [paramsFromSample(s, base)
                           for s in latinHypercube(ranges, cfg["n_param_sets"], cfg["seed"])]
    states = np.array(cfg["state_grid"], dtype=float)
    print(f"[A] {len(param_sets)} parameter sets x {len(states)} states", flush=True)

    verdicts, tables = [], {}
    verdicts += checkDimensions()
    print("[A] A1 done", flush=True)
    verdicts += checkNondimensionalization()
    print("[A] A2 done", flush=True)
    v, tables["state_grid"] = checkStateGrid(param_sets, states)
    verdicts += v
    print("[A] A3-A6 done", flush=True)
    v, tables["positivity"] = checkPositivity(param_sets, cfg["seed"], cfg["t_end"])
    verdicts += v
    print("[A] A7 done", flush=True)
    v, tables["removal_ceiling"] = checkRemovalCeiling(param_sets, cfg["t_end"])
    verdicts += v
    print("[A] A8 done", flush=True)
    v, tables["equilibrium_crossvalidation"] = checkEquilibriumCrossValidation(
        base, np.array(cfg["crossvalidation_j"], dtype=float))
    verdicts += v
    print("[A] A9 done", flush=True)
    verdicts += checkDeterminism(base, cfg["seed"], ranges)
    print("[A] A10 done", flush=True)

    outputs = {}
    for name, df in tables.items():
        outputs[f"{name}.tsv"] = writeTable(df, outdir / f"{name}.tsv")
    vdf = pd.DataFrame(verdicts)
    outputs["checks.tsv"] = writeTable(vdf, outdir / "checks.tsv")

    summary = dict(experiment="A", n_checks=len(verdicts),
                   n_failed=int((~vdf["passed"]).sum()),
                   all_passed=bool(vdf["passed"].all()),
                   checks=verdicts, runtime_s=time.time() - t0)
    with open(outdir / "checks.json", "w") as fh:
        fh.write(canonicalJson({k: v for k, v in summary.items() if k != "runtime_s"}))
        fh.write("\n")
    from proteostasis.provenance import hashFile
    outputs["checks.json"] = hashFile(outdir / "checks.json")
    writeProvenance(outdir, cfg, outputs, extra=dict(runtime_s=time.time() - t0))

    for v in verdicts:
        print(f"  {'PASS' if v['passed'] else 'FAIL'}  {v['check']:<45} "
              f"metric={v['metric']:.3e} tol={v['tolerance']:.1e}")
    print(f"[A] {summary['n_checks'] - summary['n_failed']}/{summary['n_checks']} passed "
          f"in {summary['runtime_s']:.1f}s -> {outdir}")
    return 0 if summary["all_passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
