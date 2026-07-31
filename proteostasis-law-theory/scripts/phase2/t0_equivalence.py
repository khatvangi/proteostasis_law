#!/usr/bin/env python3
"""T0 -- the epsilon -> 0 right-hand-side and jacobian identity test.

T0 is the mandatory first test of phase 2A.  it is the cheapest possible
falsification of the boron <-> nitrogen mapping: no sweep, no root finder, no
integrator.  if T0 fails, every downstream label comparison is void and the
mapping must be re-derived before a single percentage is compared.

what is asserted
----------------
1. at epsilon = 1e-6 the two vector fields agree to relative discrepancy < 1e-5
   on a fixed 12-state grid (the grid is `state_grid` from
   `configs/phase1/experiment_a.json`, reused verbatim);
2. the same at epsilon = 1e-6 for the two ANALYTIC jacobians;
3. each model's analytic jacobian matches its own central-difference jacobian,
   so a passing comparison cannot be two identical mistakes;
4. the discrepancy scales as O(epsilon) across {1e-6, 1e-3, 1e-2, 1e-1}.

why O(epsilon) and not O(epsilon**2) or O(1)
--------------------------------------------
on the ladder kappa_ca = kappa_cu and kappa_da = kappa_du, so boron's two
sequestration factors collapse to one scalar

    su = sa = 1 + epsilon * sigma,      sigma = cf/c_tot + df/d_tot in (0, 2]

and uf = u/(1 + epsilon*sigma) = u*(1 - epsilon*sigma) + O(epsilon**2).  the
free-limit model is the same expression with su = sa = 1.  therefore

    F_free - F_boron = epsilon*sigma*[u dF/duf + a dF/daf] + O(epsilon**2)

which is exactly first order, with a coefficient that vanishes nowhere generic.
an observed exponent of 2 would mean the leading term cancels, i.e. the mapping
has aligned the wrong pair of models; an exponent of 0 would mean a term is
mismatched independently of epsilon.  either falsifies section 3 of
`theory/MATCHED_IMPLEMENTATION_PROTOCOL.md`.

tolerances are derived a priori in section 5 of that document, not tuned to the
observed numbers.
"""

from __future__ import annotations

import argparse
import json
import math
import platform
import sys
from pathlib import Path
from typing import Dict, List, Sequence

import numpy as np

_HERE = Path(__file__).resolve().parent
if str(_HERE.parent) not in sys.path:                       # scripts/ on the path
    sys.path.insert(0, str(_HERE.parent))

from phase2.mapping import T0_EPSILONS, boronToNitrogen, bindingConstants   # noqa: E402
from phase2 import nitrogen_limit                                           # noqa: E402

REPO_ROOT = _HERE.parents[1]
EXPERIMENT_A = REPO_ROOT / "configs" / "phase1" / "experiment_a.json"

# --- pass criteria (see MATCHED_IMPLEMENTATION_PROTOCOL.md section 5) -------
ANCHOR_EPSILON = 1e-6
RHS_ANCHOR_TOL = 1e-5          # task requirement
JAC_ANCHOR_TOL = 1e-5
#: analytic vs central-difference jacobian of the SAME model, matrix-scaled.
#: the observed discrepancy is bounded by the difference stencil's own
#: truncation error, ~1e-9 relative; 1e-6 leaves ~2 orders of headroom while
#: still failing on any sign or factor error.
SELF_JACOBIAN_TOL = 1e-6
SLOPE_TARGET = 1.0
SLOPE_TOL_FULL = 0.05          # least squares over all four epsilons
SLOPE_TOL_ASYMPTOTIC = 0.01    # two smallest epsilons only
R2_MIN = 0.999


def loadBaseline() -> Dict[str, float]:
    """boron's phase-1 baseline parameters, read from the shipped config.

    T0 deliberately anchors on `experiment_a.json` rather than on a table typed
    into this file: if the config moves, T0 moves with it and cannot silently
    test a stale parameter set.
    """
    cfg = json.loads(EXPERIMENT_A.read_text())
    return dict(cfg["base_params"]), [tuple(s) for s in cfg["state_grid"]]


def matchedPair(base: Dict[str, float], epsilon: float):
    """(boron Params, nitrogen params) at one rung of the ladder."""
    from proteostasis.model import Params

    kb = bindingConstants(epsilon, base["c_tot"], base["d_tot"])
    kw = dict(base)
    kw.update({k: kb[k] for k in ("kappa_cu", "kappa_ca", "kappa_du", "kappa_da")})
    p = Params(**kw).validate()
    q = boronToNitrogen(p, epsilon)
    return p, q


#: relative floor for componentwise jacobian ratios, as a fraction of the
#: largest entry of the reference matrix.  an entry below this is dynamically
#: negligible and, more importantly, can be exactly zero -- see `_relativeMatrix`.
JAC_ENTRY_FLOOR_FRACTION = 1e-6


def _relative(delta: np.ndarray, reference: np.ndarray, floor: float) -> float:
    """componentwise |delta| / max(|reference|, floor), reduced by max.

    componentwise rather than norm-wise because a norm can hide a single badly
    mismatched term behind a large sibling.  the floor is the influx scale j,
    which is what stops a near-zero component from manufacturing a divergent
    ratio at states where the field genuinely vanishes.
    """
    denom = np.maximum(np.abs(reference), floor)
    return float(np.max(np.abs(delta) / denom))


def _relativeEntrywise(delta: np.ndarray, reference: np.ndarray) -> float:
    """componentwise ratio with the floor set to a fraction of the matrix scale.

    used for the CROSS-MODEL jacobian comparison, where the componentwise form
    is wanted (a mismatched small entry must not hide behind a large one) but an
    absolute-zero entry must not manufacture an infinite ratio.
    """
    scale = max(float(np.max(np.abs(reference))), 1e-300)
    denom = np.maximum(np.abs(reference), JAC_ENTRY_FLOOR_FRACTION * scale)
    return float(np.max(np.abs(delta) / denom))


def _relativeMatrix(delta: np.ndarray, reference: np.ndarray) -> float:
    """max |delta| / max |reference|, i.e. a matrix-scaled relative discrepancy.

    used for the SELF checks (analytic jacobian vs central difference of the
    same model).  a componentwise relative error is not merely inconvenient
    there, it is undefined: at the state (0, 0) the analytic entry
    d(da/dt)/du is EXACTLY zero (the nucleation term alpha_n*u**m has zero slope
    at the origin for m > 1) while the difference stencil -- clamped one-sided
    at the boundary by `model.numericalJacobian` -- returns alpha_n*h.  no
    relative statement about an exact zero exists.

    the cost of the matrix scaling is that an error smaller than
    (tolerance x largest entry) in a small entry is not detected.  at the
    tolerance used here that bound is ~1e-6 of the leading jacobian entry, which
    is far below any error that could change an eigenvalue sign.
    """
    scale = max(float(np.max(np.abs(reference))), 1e-300)
    return float(np.max(np.abs(delta)) / scale)


def compareAtEpsilon(base: Dict[str, float], grid: Sequence, epsilon: float) -> Dict:
    """all four comparisons at one epsilon, plus the per-state worst offender."""
    from proteostasis import model as boron_model

    p, q = matchedPair(base, epsilon)
    floor = float(p.j)

    rhs_rel: List[float] = []
    jac_rel: List[float] = []
    self_b: List[float] = []
    self_n: List[float] = []
    per_state = []

    for (u, a) in grid:
        fb = boron_model.rhsVector((u, a), p)
        fn = nitrogen_limit.rhsVector((u, a), q)
        jb = boron_model.jacobian(float(u), float(a), p)
        jn = nitrogen_limit.jacobian(float(u), float(a), q)
        nb = boron_model.numericalJacobian(float(u), float(a), p)
        nn = nitrogen_limit.numericalJacobian(float(u), float(a), q)

        r_rhs = _relative(fn - fb, fb, floor)
        r_jac = _relativeEntrywise(jn - jb, jb)
        r_sb = _relativeMatrix(nb - jb, jb)
        r_sn = _relativeMatrix(nn - jn, jn)

        rhs_rel.append(r_rhs)
        jac_rel.append(r_jac)
        self_b.append(r_sb)
        self_n.append(r_sn)
        per_state.append({
            "u": float(u), "a": float(a),
            "rhs_rel": r_rhs, "jac_rel": r_jac,
            "boron_self_jac_rel": r_sb, "free_self_jac_rel": r_sn,
            "F_boron": [float(v) for v in fb], "F_free": [float(v) for v in fn],
            "dF": [float(v) for v in (fn - fb)],
        })

    worst = int(np.argmax(rhs_rel))
    return {
        "epsilon": float(epsilon),
        "boron_params": p.asDict(),
        "nitrogen_params": q.asDict(),
        "influx_floor": floor,
        "rhs_rel_max": float(np.max(rhs_rel)),
        "jac_rel_max": float(np.max(jac_rel)),
        "boron_self_jac_rel_max": float(np.max(self_b)),
        "free_self_jac_rel_max": float(np.max(self_n)),
        "worst_state": {"u": per_state[worst]["u"], "a": per_state[worst]["a"]},
        "per_state": per_state,
    }


def logLogSlope(eps: Sequence[float], disc: Sequence[float]):
    """least-squares slope and R^2 of log(discrepancy) against log(epsilon)."""
    x = np.log(np.asarray(eps, dtype=float))
    y = np.log(np.asarray(disc, dtype=float))
    slope, intercept = np.polyfit(x, y, 1)
    pred = slope * x + intercept
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return float(slope), float(intercept), float(r2)


def runT0(epsilons: Sequence[float] = T0_EPSILONS) -> Dict:
    base, grid = loadBaseline()
    cells = [compareAtEpsilon(base, grid, e) for e in epsilons]

    eps = [c["epsilon"] for c in cells]
    rhs_d = [c["rhs_rel_max"] for c in cells]
    jac_d = [c["jac_rel_max"] for c in cells]

    order = np.argsort(eps)
    eps_s = [eps[i] for i in order]
    rhs_s = [rhs_d[i] for i in order]
    jac_s = [jac_d[i] for i in order]

    slope_rhs, _, r2_rhs = logLogSlope(eps_s, rhs_s)
    slope_jac, _, r2_jac = logLogSlope(eps_s, jac_s)
    asym_rhs = (math.log(rhs_s[1] / rhs_s[0]) / math.log(eps_s[1] / eps_s[0]))
    asym_jac = (math.log(jac_s[1] / jac_s[0]) / math.log(eps_s[1] / eps_s[0]))

    anchor = next(c for c in cells if abs(c["epsilon"] - ANCHOR_EPSILON) < 1e-30)

    checks = {
        "anchor_rhs_below_tol": anchor["rhs_rel_max"] < RHS_ANCHOR_TOL,
        "anchor_jac_below_tol": anchor["jac_rel_max"] < JAC_ANCHOR_TOL,
        "boron_analytic_matches_numeric": all(
            c["boron_self_jac_rel_max"] < SELF_JACOBIAN_TOL for c in cells),
        "free_analytic_matches_numeric": all(
            c["free_self_jac_rel_max"] < SELF_JACOBIAN_TOL for c in cells),
        "rhs_slope_full_is_linear": abs(slope_rhs - SLOPE_TARGET) <= SLOPE_TOL_FULL,
        "rhs_slope_asymptotic_is_linear": abs(asym_rhs - SLOPE_TARGET) <= SLOPE_TOL_ASYMPTOTIC,
        "rhs_loglog_linear": r2_rhs >= R2_MIN,
        "jac_slope_full_is_linear": abs(slope_jac - SLOPE_TARGET) <= SLOPE_TOL_FULL,
        "jac_slope_asymptotic_is_linear": abs(asym_jac - SLOPE_TARGET) <= SLOPE_TOL_ASYMPTOTIC,
        "jac_loglog_linear": r2_jac >= R2_MIN,
    }

    return {
        "test": "T0_epsilon_zero_identity",
        "passed": all(checks.values()),
        "checks": checks,
        "criteria": {
            "anchor_epsilon": ANCHOR_EPSILON,
            "rhs_anchor_tol": RHS_ANCHOR_TOL,
            "jac_anchor_tol": JAC_ANCHOR_TOL,
            "self_jacobian_tol": SELF_JACOBIAN_TOL,
            "slope_target": SLOPE_TARGET,
            "slope_tol_full": SLOPE_TOL_FULL,
            "slope_tol_asymptotic": SLOPE_TOL_ASYMPTOTIC,
            "r2_min": R2_MIN,
        },
        "scaling": {
            "epsilons": eps_s,
            "rhs_rel_max": rhs_s,
            "jac_rel_max": jac_s,
            "rhs_slope_loglog": slope_rhs, "rhs_r2": r2_rhs,
            "rhs_slope_asymptotic": asym_rhs,
            "jac_slope_loglog": slope_jac, "jac_r2": r2_jac,
            "jac_slope_asymptotic": asym_jac,
        },
        "cells": cells,
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "node": platform.node(),
        },
        "state_grid_source": str(EXPERIMENT_A.relative_to(REPO_ROOT)),
    }


def failureReport(report: Dict) -> str:
    """exact mismatched terms, for the stop-work report demanded on failure."""
    lines = ["T0 FAILED -- downstream classification is void.", ""]
    for name, ok in report["checks"].items():
        lines.append(f"  {'PASS' if ok else 'FAIL'}  {name}")
    lines.append("")
    for c in report["cells"]:
        if c["rhs_rel_max"] < RHS_ANCHOR_TOL and c["jac_rel_max"] < JAC_ANCHOR_TOL:
            continue
        lines.append(f"epsilon = {c['epsilon']:g}: rhs_rel_max = {c['rhs_rel_max']:.6e}, "
                     f"jac_rel_max = {c['jac_rel_max']:.6e}")
        worst = max(c["per_state"], key=lambda s: s["rhs_rel"])
        lines.append(f"  worst state (u, a) = ({worst['u']:g}, {worst['a']:g})")
        lines.append(f"    F_boron = {worst['F_boron']}")
        lines.append(f"    F_free  = {worst['F_free']}")
        lines.append(f"    dF      = {worst['dF']}")
    return "\n".join(lines)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--json", type=Path, default=None,
                    help="write the full report here")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    report = runT0()
    if args.json:
        args.json.parent.mkdir(parents=True, exist_ok=True)
        args.json.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    if not args.quiet:
        s = report["scaling"]
        print(f"T0 {'PASSED' if report['passed'] else 'FAILED'}")
        print(f"{'epsilon':>10}  {'rhs_rel_max':>14}  {'jac_rel_max':>14}")
        for e, r, j in zip(s["epsilons"], s["rhs_rel_max"], s["jac_rel_max"]):
            print(f"{e:>10.0e}  {r:>14.6e}  {j:>14.6e}")
        print(f"rhs slope (log-log, 4 pt) = {s['rhs_slope_loglog']:.6f}  "
              f"R^2 = {s['rhs_r2']:.8f}  asymptotic = {s['rhs_slope_asymptotic']:.6f}")
        print(f"jac slope (log-log, 4 pt) = {s['jac_slope_loglog']:.6f}  "
              f"R^2 = {s['jac_r2']:.8f}  asymptotic = {s['jac_slope_asymptotic']:.6f}")
        for name, ok in report["checks"].items():
            print(f"  {'PASS' if ok else 'FAIL'}  {name}")
        if not report["passed"]:
            print()
            print(failureReport(report))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
