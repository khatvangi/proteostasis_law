"""minimal nondimensional two-state conserved-resource proteostasis model.

state (both are TOTAL pools, i.e. free plus machinery-bound):
    u   soluble damaged / misfolded monomer
    a   aggregate burden, carried in monomer-equivalent units so that transfer
        between u and a is mass-conserving 1:1

equations (nondimensional; concentration scale phi = C_tot + D_tot, time scale
tau = 1/k_ref):

    du/dt = j
            - cf * uf/(kappa_ref + uf)          refolding to native
            - rho_U * df * uf/(kappa_u + uf)    soluble degradation
            - alpha_n * uf**m                   nucleation, m > 1
            - alpha_g * uf * af                 aggregate growth
            + alpha_d * cf * af/(kappa_dis+af)  disaggregation returns to u

    da/dt = alpha_n * uf**m
            + alpha_g * uf * af
            - alpha_d * cf * af/(kappa_dis+af)
            - rho_A * df * af/(kappa_a + af)    aggregate clearance

rapid-equilibrium conserved pools, written in FREE concentrations only
(DECISIONS.md D004 -- total substrate is never inserted into a free-resource
binding formula):

    uf = u / (1 + cf/kappa_cu + df/kappa_du)
    af = a / (1 + cf/kappa_ca + df/kappa_da)
    cf = c_tot / (1 + nu + uf/kappa_cu + af/kappa_ca)
    df = d_tot / (1 + uf/kappa_du + af/kappa_da)

`nu = N_free/K_N` is the ordinary nascent-chain chaperone occupancy load. it is
a separate parameter and contributes NO damage influx; it acts purely by
consuming chaperone capacity, which is the mechanism the theory requires (see
theory/DYNAMICAL_SYSTEM.md).

structural facts used and tested elsewhere in this package:

  * the nonnegative orthant is forward invariant: at u=0 we have uf=0 so
    du/dt = j + disaggregation >= 0; at a=0 we have af=0 so da/dt =
    alpha_n*uf**m >= 0.
  * mass balance: d(u+a)/dt = j - v_ref - v_degU - v_degA exactly. the internal
    u<->a transfer cancels.
  * removal is bounded above by c_tot + (rho_U + rho_A)*d_tot, because
    cf <= c_tot, df <= d_tot and every saturating factor is <= 1. therefore
    j > c_tot + (rho_U + rho_A)*d_tot admits NO bounded state: burden grows
    without bound. this is an analytic necessary condition on feasibility.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict, replace
from typing import Dict, Optional, Tuple

import numpy as np

# free-pool solver tolerances
_BIND_TOL = 1e-13
_BIND_MAX_NEWTON = 60
_BIND_MAX_FIXEDPOINT = 200000
_TINY = 1e-300


class ModelError(ValueError):
    """raised for invalid parameters or a failed free-pool solve."""


@dataclass(frozen=True)
class Params:
    """nondimensional parameters. every field is strictly positive except that
    `nu` may be zero (no ordinary nascent-chain load)."""

    # --- loads -------------------------------------------------------------
    j: float = 0.05        # site-resolved damage influx, scaled
    nu: float = 0.5        # ordinary nascent-chain chaperone occupancy load

    # --- conserved rescue pools (scaled so c_tot + d_tot = 1 at baseline) ---
    c_tot: float = 0.6     # total chaperone
    d_tot: float = 0.4     # total degradation machinery

    # --- kinetics ----------------------------------------------------------
    rho_U: float = 1.0     # soluble degradation turnover / refolding turnover
    rho_A: float = 0.5     # aggregate clearance turnover / refolding turnover
    alpha_n: float = 0.5   # nucleation coefficient
    alpha_g: float = 1.0   # aggregate growth coefficient
    alpha_d: float = 0.3   # disaggregation coefficient
    m: float = 2.0         # nucleation order, must exceed 1

    # --- michaelis / dissociation constants (scaled) ------------------------
    kappa_ref: float = 0.5
    kappa_u: float = 0.5
    kappa_a: float = 0.5
    kappa_dis: float = 0.5
    kappa_cu: float = 0.3
    kappa_ca: float = 0.3
    kappa_du: float = 0.5
    kappa_da: float = 0.5

    def validate(self) -> "Params":
        for name in ("c_tot", "d_tot", "rho_U", "rho_A", "alpha_n", "alpha_g",
                     "alpha_d", "kappa_ref", "kappa_u", "kappa_a", "kappa_dis",
                     "kappa_cu", "kappa_ca", "kappa_du", "kappa_da"):
            if not (getattr(self, name) > 0.0 and np.isfinite(getattr(self, name))):
                raise ModelError(f"parameter '{name}' must be finite and positive")
        if self.j < 0.0 or not np.isfinite(self.j):
            raise ModelError("influx j must be finite and nonnegative")
        if self.nu < 0.0 or not np.isfinite(self.nu):
            raise ModelError("nascent occupancy nu must be finite and nonnegative")
        if not self.m > 1.0:
            raise ModelError("nucleation order m must exceed 1 (m<=1 breaks the fold structure)")
        return self

    def asDict(self) -> Dict[str, float]:
        return asdict(self)

    def with_(self, **kw) -> "Params":
        return replace(self, **kw)

    @property
    def rescueTotal(self) -> float:
        return self.c_tot + self.d_tot

    @property
    def chaperoneShare(self) -> float:
        """rescue allocation chi: fraction of rescue machinery that is chaperone."""
        return self.c_tot / (self.c_tot + self.d_tot)


def paramsFromDict(d: Dict[str, float]) -> Params:
    known = {f: d[f] for f in Params.__dataclass_fields__ if f in d}
    return Params(**known).validate()


def allocationParams(base: Params, chi: float, rescue_total: Optional[float] = None) -> Params:
    """re-express rescue as (allocation chi, total rescue R).

    chi is the chaperone share; 1-chi goes to degradation. this is the knob
    used by experiment B ('rescue allocation') and by the capacity-knockdown
    arm of experiment D.
    """
    if not (0.0 < chi < 1.0):
        raise ModelError("rescue allocation chi must lie strictly in (0,1)")
    R = base.rescueTotal if rescue_total is None else rescue_total
    if not R > 0.0:
        raise ModelError("total rescue must be positive")
    return base.with_(c_tot=chi * R, d_tot=(1.0 - chi) * R)


# ---------------------------------------------------------------------------
# conserved-pool algebra
# ---------------------------------------------------------------------------


def _bindingResidual(cf: float, df: float, u: float, a: float, p: Params):
    """residual of the rapid-equilibrium pool balances plus its 2x2 jacobian.

    unknowns are the free resource concentrations (cf, df). free substrate is
    eliminated analytically:  uf = u/su,  af = a/sa.
    """
    su = 1.0 + cf / p.kappa_cu + df / p.kappa_du
    sa = 1.0 + cf / p.kappa_ca + df / p.kappa_da
    uf = u / su
    af = a / sa

    # partials of free substrate w.r.t. free resources
    duf_dcf = -uf / (su * p.kappa_cu)
    duf_ddf = -uf / (su * p.kappa_du)
    daf_dcf = -af / (sa * p.kappa_ca)
    daf_ddf = -af / (sa * p.kappa_da)

    qc = 1.0 + p.nu + uf / p.kappa_cu + af / p.kappa_ca
    qd = 1.0 + uf / p.kappa_du + af / p.kappa_da

    r1 = cf * qc - p.c_tot
    r2 = df * qd - p.d_tot

    j11 = qc + cf * (duf_dcf / p.kappa_cu + daf_dcf / p.kappa_ca)
    j12 = cf * (duf_ddf / p.kappa_cu + daf_ddf / p.kappa_ca)
    j21 = df * (duf_dcf / p.kappa_du + daf_dcf / p.kappa_da)
    j22 = qd + df * (duf_ddf / p.kappa_du + daf_ddf / p.kappa_da)

    aux = dict(su=su, sa=sa, uf=uf, af=af, duf_dcf=duf_dcf, duf_ddf=duf_ddf,
               daf_dcf=daf_dcf, daf_ddf=daf_ddf)
    return (r1, r2), (j11, j12, j21, j22), aux


def _fixedPointMap(cf: float, df: float, u: float, a: float, p: Params):
    """one application of the rapid-equilibrium map T(cf, df).

    T is monotone INCREASING in both arguments (raising a free resource lowers
    free substrate, which lowers occupancy, which raises the free resource) and
    maps [0,c_tot] x [0,d_tot] into itself. knaster-tarski therefore guarantees
    a least and a greatest fixed point, obtainable by iterating from (0,0) and
    from (c_tot,d_tot). this is the uniqueness certificate used by
    `solveFreePoolsCertified`.
    """
    uf = u / (1.0 + cf / p.kappa_cu + df / p.kappa_du)
    af = a / (1.0 + cf / p.kappa_ca + df / p.kappa_da)
    return (p.c_tot / (1.0 + p.nu + uf / p.kappa_cu + af / p.kappa_ca),
            p.d_tot / (1.0 + uf / p.kappa_du + af / p.kappa_da))


def _iterateMonotone(u: float, a: float, p: Params, fromBelow: bool,
                     tol: float = _BIND_TOL, maxiter: int = _BIND_MAX_FIXEDPOINT):
    """monotone iteration to the least (fromBelow) or greatest fixed point."""
    cf, df = (0.0, 0.0) if fromBelow else (p.c_tot, p.d_tot)
    for k in range(maxiter):
        cf_new, df_new = _fixedPointMap(cf, df, u, a, p)
        delta = max(abs(cf_new - cf), abs(df_new - df))
        cf, df = cf_new, df_new
        if delta <= tol * max(1.0, p.c_tot, p.d_tot):
            return cf, df, k + 1
    return cf, df, maxiter


def solveFreePools(u: float, a: float, p: Params,
                   guess: Optional[Tuple[float, float]] = None):
    """solve jointly for (uf, af, cf, df) at total burdens (u, a).

    safeguarded newton on the 2-d residual with a monotone fixed-point fallback.
    returns (uf, af, cf, df).
    """
    if u < 0.0 or a < 0.0:
        raise ModelError(f"negative total burden passed to free-pool solve: u={u}, a={a}")

    lo_c, hi_c = 1e-16 * p.c_tot, p.c_tot
    lo_d, hi_d = 1e-16 * p.d_tot, p.d_tot
    if guess is None:
        cf = p.c_tot / (1.0 + p.nu)
        df = p.d_tot
    else:
        cf = min(max(guess[0], lo_c), hi_c)
        df = min(max(guess[1], lo_d), hi_d)

    (r1, r2), jac, aux = _bindingResidual(cf, df, u, a, p)
    norm = max(abs(r1), abs(r2))
    for _ in range(_BIND_MAX_NEWTON):
        if norm <= _BIND_TOL * max(1.0, p.c_tot, p.d_tot):
            break
        j11, j12, j21, j22 = jac
        det = j11 * j22 - j12 * j21
        if not np.isfinite(det) or abs(det) < 1e-14:
            break
        dc = -(j22 * r1 - j12 * r2) / det
        dd = -(-j21 * r1 + j11 * r2) / det

        # damped, box-projected step
        step = 1.0
        accepted = False
        for _ in range(40):
            cf_try = min(max(cf + step * dc, lo_c), hi_c)
            df_try = min(max(df + step * dd, lo_d), hi_d)
            (t1, t2), jac_try, aux_try = _bindingResidual(cf_try, df_try, u, a, p)
            n_try = max(abs(t1), abs(t2))
            if np.isfinite(n_try) and n_try < norm:
                cf, df, r1, r2, jac, aux, norm = cf_try, df_try, t1, t2, jac_try, aux_try, n_try
                accepted = True
                break
            step *= 0.5
        if not accepted:
            break

    if norm > 1e-9 * max(1.0, p.c_tot, p.d_tot):
        # newton stalled; fall back to the guaranteed-convergent monotone iteration
        cf, df, _ = _iterateMonotone(u, a, p, fromBelow=True)
        (r1, r2), jac, aux = _bindingResidual(cf, df, u, a, p)
        if max(abs(r1), abs(r2)) > 1e-6 * max(1.0, p.c_tot, p.d_tot):
            raise ModelError(f"free-pool solve failed at u={u:.6g}, a={a:.6g}")

    return aux["uf"], aux["af"], cf, df


def solveFreePoolsCertified(u: float, a: float, p: Params, tol: float = 1e-10):
    """solve the free pools AND certify uniqueness.

    iterates the monotone map from below and from above. by knaster-tarski the
    two limits bracket every fixed point, so if they agree to `tol` the
    rapid-equilibrium solution is unique at this state.

    returns (uf, af, cf, df, certificate_dict).
    """
    lo_c, lo_d, n_lo = _iterateMonotone(u, a, p, fromBelow=True)
    hi_c, hi_d, n_hi = _iterateMonotone(u, a, p, fromBelow=False)
    gap = max(abs(hi_c - lo_c), abs(hi_d - lo_d))
    unique = bool(gap <= tol * max(1.0, p.c_tot, p.d_tot))

    cf, df = 0.5 * (lo_c + hi_c), 0.5 * (lo_d + hi_d)
    (r1, r2), _, aux = _bindingResidual(cf, df, u, a, p)
    cert = dict(least=(lo_c, lo_d), greatest=(hi_c, hi_d), gap=gap, unique=unique,
                residual=max(abs(r1), abs(r2)), iters_below=n_lo, iters_above=n_hi)
    return aux["uf"], aux["af"], cf, df, cert


def poolBalanceResiduals(u: float, a: float, p: Params) -> Dict[str, float]:
    """mass-balance residuals of the algebraic closure, for invariant tests.

    each entry should be ~0: the reconstructed totals must equal the totals we
    started from, and the reconstructed pools must equal c_tot / d_tot.
    """
    uf, af, cf, df = solveFreePools(u, a, p)
    u_rebuilt = uf * (1.0 + cf / p.kappa_cu + df / p.kappa_du)
    a_rebuilt = af * (1.0 + cf / p.kappa_ca + df / p.kappa_da)
    c_rebuilt = cf * (1.0 + p.nu + uf / p.kappa_cu + af / p.kappa_ca)
    d_rebuilt = df * (1.0 + uf / p.kappa_du + af / p.kappa_da)
    scale = max(1.0, u, a, p.c_tot, p.d_tot)
    return {
        "u_total": (u_rebuilt - u) / scale,
        "a_total": (a_rebuilt - a) / scale,
        "c_total": (c_rebuilt - p.c_tot) / scale,
        "d_total": (d_rebuilt - p.d_tot) / scale,
    }


# ---------------------------------------------------------------------------
# right-hand side, fluxes, jacobian
# ---------------------------------------------------------------------------


def fluxes(u: float, a: float, p: Params,
           guess: Optional[Tuple[float, float]] = None) -> Dict[str, float]:
    """all individual mass fluxes at a state; used for mass-balance checks."""
    uf, af, cf, df = solveFreePools(u, a, p, guess)
    return {
        "influx": p.j,
        "refold": cf * uf / (p.kappa_ref + uf),
        "degrade_u": p.rho_U * df * uf / (p.kappa_u + uf),
        "nucleate": p.alpha_n * uf ** p.m,
        "grow": p.alpha_g * uf * af,
        "disaggregate": p.alpha_d * cf * af / (p.kappa_dis + af),
        "degrade_a": p.rho_A * df * af / (p.kappa_a + af),
        "uf": uf, "af": af, "cf": cf, "df": df,
    }


def rhs(u: float, a: float, p: Params,
        guess: Optional[Tuple[float, float]] = None) -> Tuple[float, float]:
    """du/dt, da/dt."""
    f = fluxes(u, a, p, guess)
    du = (f["influx"] - f["refold"] - f["degrade_u"] - f["nucleate"]
          - f["grow"] + f["disaggregate"])
    da = f["nucleate"] + f["grow"] - f["disaggregate"] - f["degrade_a"]
    return du, da


def rhsVector(x, p: Params) -> np.ndarray:
    du, da = rhs(float(x[0]), float(x[1]), p)
    return np.array([du, da], dtype=float)


def massBalanceResidual(u: float, a: float, p: Params) -> float:
    """relative residual of d(u+a)/dt == influx - refold - degrade_u - degrade_a.

    the internal u->a transfer (nucleation + growth - disaggregation) must
    cancel exactly, so a nonzero residual means broken transfer stoichiometry.

    the residual is normalised by the LARGEST flux entering the sum, not by the
    surviving removal fluxes. at high burden the growth term can exceed the net
    rate by orders of magnitude, and a cancellation-limited difference can only
    be judged against the magnitude of the terms that cancelled.
    """
    f = fluxes(u, a, p)
    du = (f["influx"] - f["refold"] - f["degrade_u"] - f["nucleate"]
          - f["grow"] + f["disaggregate"])
    da = f["nucleate"] + f["grow"] - f["disaggregate"] - f["degrade_a"]
    expected = f["influx"] - f["refold"] - f["degrade_u"] - f["degrade_a"]
    scale = max(abs(f[k]) for k in ("influx", "refold", "degrade_u", "degrade_a",
                                    "nucleate", "grow", "disaggregate"))
    return ((du + da) - expected) / max(scale, 1e-300)


def jacobian(u: float, a: float, p: Params,
             guess: Optional[Tuple[float, float]] = None) -> np.ndarray:
    """analytic jacobian d(du/dt, da/dt)/d(u, a).

    the rhs depends on (u,a) only through the free concentrations, which are
    defined implicitly by the pool balances. the implicit function theorem gives

        d(cf,df)/d(u,a) = -J_R^{-1} * dR/d(u,a)

    and free substrate then follows by the chain rule.
    """
    uf, af, cf, df = solveFreePools(u, a, p, guess)
    (_, _), (j11, j12, j21, j22), aux = _bindingResidual(cf, df, u, a, p)
    su, sa = aux["su"], aux["sa"]

    # explicit dependence of the pool residuals on the total burdens
    r1_u = cf / (su * p.kappa_cu)
    r1_a = cf / (sa * p.kappa_ca)
    r2_u = df / (su * p.kappa_du)
    r2_a = df / (sa * p.kappa_da)

    det = j11 * j22 - j12 * j21
    if not np.isfinite(det) or abs(det) < 1e-300:
        raise ModelError("singular binding jacobian; cannot differentiate free pools")
    inv = np.array([[j22, -j12], [-j21, j11]], dtype=float) / det
    dres = np.array([[r1_u, r1_a], [r2_u, r2_a]], dtype=float)
    d_resources = -inv @ dres                       # rows: cf, df; cols: u, a
    dcf_du, dcf_da = d_resources[0]
    ddf_du, ddf_da = d_resources[1]

    duf_dcf, duf_ddf = aux["duf_dcf"], aux["duf_ddf"]
    daf_dcf, daf_ddf = aux["daf_dcf"], aux["daf_ddf"]

    duf_du = 1.0 / su + duf_dcf * dcf_du + duf_ddf * ddf_du
    duf_da = duf_dcf * dcf_da + duf_ddf * ddf_da
    daf_du = daf_dcf * dcf_du + daf_ddf * ddf_du
    daf_da = 1.0 / sa + daf_dcf * dcf_da + daf_ddf * ddf_da

    # partials of the rhs w.r.t. the free concentrations
    ref_den = (p.kappa_ref + uf) ** 2
    degu_den = (p.kappa_u + uf) ** 2
    dis_den = (p.kappa_dis + af) ** 2
    dega_den = (p.kappa_a + af) ** 2
    nuc_prime = p.alpha_n * p.m * (uf ** (p.m - 1.0)) if uf > 0.0 else 0.0

    f1_uf = (-cf * p.kappa_ref / ref_den - p.rho_U * df * p.kappa_u / degu_den
             - nuc_prime - p.alpha_g * af)
    f1_af = -p.alpha_g * uf + p.alpha_d * cf * p.kappa_dis / dis_den
    f1_cf = -uf / (p.kappa_ref + uf) + p.alpha_d * af / (p.kappa_dis + af)
    f1_df = -p.rho_U * uf / (p.kappa_u + uf)

    f2_uf = nuc_prime + p.alpha_g * af
    f2_af = (p.alpha_g * uf - p.alpha_d * cf * p.kappa_dis / dis_den
             - p.rho_A * df * p.kappa_a / dega_den)
    f2_cf = -p.alpha_d * af / (p.kappa_dis + af)
    f2_df = -p.rho_A * af / (p.kappa_a + af)

    j_uu = f1_uf * duf_du + f1_af * daf_du + f1_cf * dcf_du + f1_df * ddf_du
    j_ua = f1_uf * duf_da + f1_af * daf_da + f1_cf * dcf_da + f1_df * ddf_da
    j_au = f2_uf * duf_du + f2_af * daf_du + f2_cf * dcf_du + f2_df * ddf_du
    j_aa = f2_uf * duf_da + f2_af * daf_da + f2_cf * dcf_da + f2_df * ddf_da
    return np.array([[j_uu, j_ua], [j_au, j_aa]], dtype=float)


def numericalJacobian(u: float, a: float, p: Params, h_rel: float = 1e-6) -> np.ndarray:
    """central-difference jacobian, for validating the analytic one."""
    x = np.array([u, a], dtype=float)
    out = np.zeros((2, 2))
    for k in range(2):
        h = h_rel * max(abs(x[k]), 1e-3)
        xp, xm = x.copy(), x.copy()
        xp[k] += h
        xm[k] = max(xm[k] - h, 0.0)
        step = xp[k] - xm[k]
        out[:, k] = (rhsVector(xp, p) - rhsVector(xm, p)) / step
    return out


# ---------------------------------------------------------------------------
# analytic bounds
# ---------------------------------------------------------------------------


def removalCeiling(p: Params) -> float:
    """supremum of total removal flux over the whole nonnegative orthant.

    cf <= c_tot, df <= d_tot and every saturating factor is < 1, so
        v_ref + v_degU + v_degA < c_tot + (rho_U + rho_A) * d_tot.
    if j exceeds this, d(u+a)/dt > 0 everywhere and NO bounded state exists.
    """
    return p.c_tot + (p.rho_U + p.rho_A) * p.d_tot


def influxAdmitsNoBoundedState(p: Params) -> bool:
    """analytic sufficient condition for infeasibility."""
    return p.j > removalCeiling(p)
