"""independently transcribed FREE-SUBSTRATE (epsilon = 0) right-hand side.

this is a local transcription of nitrogen's
`/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory-nitrogen-check/
src/proteostasis_model.py` (sha256 recorded in
`scripts/phase2/NITROGEN_SOURCE.sha256`), kept here deliberately rather than
imported across hosts so that:

  * the T0 test is runnable on either machine with no ssh dependency;
  * the benchmark's free-limit arm and its boron arm run in one process, under
    one numpy/scipy build, so a version difference cannot masquerade as a model
    difference.

fidelity of the transcription is not asserted, it is TESTED:
`tests/phase2/test_nitrogen_transcription.py` compares this module against a
reference dump produced by executing the real nitrogen module on nitrogen
(`data/phase2/nitrogen_reference.json`), on the same states and parameters.

state convention: u, a are FREE binding-competent concentrations.  there is no
substrate sequestration factor; su = sa = 1 identically.  that is the epsilon = 0
face of the boron family, not an approximation to it.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

from .mapping import NitrogenParams


def resourceFractions(u, a, q: NitrogenParams):
    """free chaperone / degradation fractions and every bound pool.

    transcribed from proteostasis_model.resource_fractions.  because u and a
    are free, this closed form is exact -- no implicit solve is required, which
    is the entire structural difference from boron.
    """
    cq = 1.0 + q.n_load + q.o_load + q.c_u * u + q.c_a * a
    dq = 1.0 + q.d_u * u + q.d_a * a
    cf = q.c_tot / cq
    df = q.d_tot / dq
    pools = {
        "C_f": cf, "C_N": cf * q.n_load, "C_U": cf * q.c_u * u,
        "C_A": cf * q.c_a * a, "C_O": cf * q.o_load,
        "D_f": df, "D_U": df * q.d_u * u, "D_A": df * q.d_a * a,
    }
    return cf, df, pools


def fluxes(u, a, q: NitrogenParams) -> Dict[str, float]:
    """individual mass fluxes, named to line up with `model.fluxes`."""
    cf, df, _ = resourceFractions(u, a, q)
    return {
        "influx": q.j,
        "refold": q.ref_capacity * cf * u / (q.ref_k + u),
        "degrade_u": q.deg_u_capacity * df * u / (q.deg_u_k + u),
        "nucleate": q.nucleation * u ** q.m,
        "grow": q.growth * u * a,
        "disaggregate": q.disaggregation * cf * a / (q.disaggregation_k + a),
        "degrade_a": q.deg_a_capacity * df * a / (q.deg_a_k + a),
        "uf": u, "af": a, "cf": cf, "df": df,
    }


def rhs(u, a, q: NitrogenParams) -> Tuple[float, float]:
    """du/dt, da/dt.  signature matches `model.rhs` for interchangeability."""
    f = fluxes(u, a, q)
    du = (f["influx"] - f["refold"] - f["degrade_u"] - f["nucleate"]
          - f["grow"] + f["disaggregate"])
    da = f["nucleate"] + f["grow"] - f["disaggregate"] - f["degrade_a"]
    return du, da


def rhsVector(x, q: NitrogenParams) -> np.ndarray:
    du, da = rhs(x[0], x[1], q)
    return np.array([du, da], dtype=float)


def jacobian(u, a, q: NitrogenParams) -> np.ndarray:
    """analytic jacobian.  transcribed from proteostasis_model.jacobian.

    dtype follows the state so a complex-step check works unchanged; the
    closed-form resource derivatives are dcf/du = -cf**2 * c_u / c_tot, which is
    just the quotient rule on cf = c_tot/cq.
    """
    cf, df, _ = resourceFractions(u, a, q)
    cfu, cfa = -cf * cf * q.c_u / q.c_tot, -cf * cf * q.c_a / q.c_tot
    dfu, dfa = -df * df * q.d_u / q.d_tot, -df * df * q.d_a / q.d_tot
    ref_u = q.ref_capacity * (cfu * u / (q.ref_k + u)
                              + cf * q.ref_k / (q.ref_k + u) ** 2)
    ref_a = q.ref_capacity * cfa * u / (q.ref_k + u)
    du_u = q.deg_u_capacity * (dfu * u / (q.deg_u_k + u)
                               + df * q.deg_u_k / (q.deg_u_k + u) ** 2)
    du_a = q.deg_u_capacity * dfa * u / (q.deg_u_k + u)
    nuc_u = q.nucleation * q.m * u ** (q.m - 1)
    dis_u = q.disaggregation * cfu * a / (q.disaggregation_k + a)
    dis_a = q.disaggregation * (cfa * a / (q.disaggregation_k + a)
                                + cf * q.disaggregation_k / (q.disaggregation_k + a) ** 2)
    da_u = q.deg_a_capacity * dfu * a / (q.deg_a_k + a)
    da_a = q.deg_a_capacity * (dfa * a / (q.deg_a_k + a)
                               + df * q.deg_a_k / (q.deg_a_k + a) ** 2)
    grow_u, grow_a = q.growth * a, q.growth * u
    return np.array([[-ref_u - du_u - nuc_u - grow_u + dis_u,
                      -ref_a - du_a - grow_a + dis_a],
                     [nuc_u + grow_u - dis_u - da_u,
                      grow_a - dis_a - da_a]],
                    dtype=np.result_type(np.asarray([u, a]), float))


def numericalJacobian(u: float, a: float, q: NitrogenParams,
                      h_rel: float = 1e-6) -> np.ndarray:
    """central-difference jacobian.

    the step rule mirrors `model.numericalJacobian` exactly (including the
    clamp at zero, which makes the boundary stencil one-sided) so that a
    boron-vs-nitrogen comparison of NUMERICAL jacobians is not contaminated by
    a differencing-scheme difference.
    """
    x = np.array([u, a], dtype=float)
    out = np.zeros((2, 2))
    for k in range(2):
        h = h_rel * max(abs(x[k]), 1e-3)
        xp, xm = x.copy(), x.copy()
        xp[k] += h
        xm[k] = max(xm[k] - h, 0.0)
        step = xp[k] - xm[k]
        out[:, k] = (rhsVector(xp, q) - rhsVector(xm, q)) / step
    return out


def poolBalanceResiduals(u: float, a: float, q: NitrogenParams) -> Dict[str, float]:
    """resource-conservation residuals; each entry should be ~0."""
    cf, df, pools = resourceFractions(u, a, q)
    scale = max(1.0, q.c_tot, q.d_tot)
    c_rebuilt = pools["C_f"] + pools["C_N"] + pools["C_U"] + pools["C_A"] + pools["C_O"]
    d_rebuilt = pools["D_f"] + pools["D_U"] + pools["D_A"]
    return {"c_total": (c_rebuilt - q.c_tot) / scale,
            "d_total": (d_rebuilt - q.d_tot) / scale}


# ---------------------------------------------------------------------------
# equilibrium / stability / integration, transcribed for protocol parity
# ---------------------------------------------------------------------------


def classifyStability(x, q: NitrogenParams, tol: float = 1e-8):
    eig = np.linalg.eigvals(jacobian(x[0], x[1], q))
    if np.all(np.real(eig) < -tol):
        label = "stable"
    elif np.any(np.real(eig) > tol):
        label = "unstable"
    else:
        label = "nonhyperbolic"
    return label, eig


def solveEquilibrium(q: NitrogenParams, guess=(0.05, 0.02), tol: float = 1e-10):
    """transcribed from proteostasis_model.solve_equilibrium (linear coords)."""
    sol = root(lambda z: rhsVector(z, q), np.asarray(guess, dtype=float),
               jac=lambda z: jacobian(z[0], z[1], q), tol=tol)
    residual = float(np.linalg.norm(rhsVector(sol.x, q), ord=np.inf))
    return {"x": sol.x, "residual": residual,
            "success": bool(sol.success and np.isfinite(residual)),
            "message": str(sol.message), "nfev": int(sol.nfev)}


def integrate(x0, q: NitrogenParams, t_span=(0.0, 100.0), *, rtol: float = 1e-9,
              atol: float = 1e-11, max_step=np.inf, dense_output: bool = False):
    """transcribed from proteostasis_model.integrate (DOP853, no events)."""
    return solve_ivp(lambda t, x: rhsVector(x, q), t_span,
                     np.asarray(x0, dtype=float), method="DOP853", rtol=rtol,
                     atol=atol, max_step=max_step, dense_output=dense_output)


def inAdmissibleRegion(x, q: NitrogenParams,
                       thresholds: Optional[Tuple[float, float]] = None) -> bool:
    """componentwise D = {u < H_u, a < H_a, u,a >= 0}."""
    hu, ha = thresholds or (q.threshold_u, q.threshold_a)
    return bool(np.all(np.asarray(x) >= 0) and x[0] < hu and x[1] < ha)


def removalCeiling(q: NitrogenParams) -> float:
    return q.removalCeiling()


def residualScale(q: NitrogenParams) -> float:
    return q.residualScale()


def freeCoordinates(u: float, a: float, q: NitrogenParams) -> Tuple[float, float]:
    """free substrate at a state.  identically the state itself in this
    convention -- provided so the benchmark can call the same accessor on both
    model forms and record both coordinate systems for every equilibrium."""
    return float(u), float(a)
