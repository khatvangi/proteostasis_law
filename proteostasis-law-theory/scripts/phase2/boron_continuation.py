"""analytic continuation of the boron field to a neighbourhood of the orthant.

why this module exists
----------------------
nitrogen's root protocol runs `scipy.optimize.root(method='hybr')` in LINEAR
coordinates.  hybr's line search freely probes states with a < 0.  nitrogen's
model is a rational expression and evaluates there without complaint; boron's
`solveFreePools` raises `ModelError`, because the rapid-equilibrium conservation
solve is stated for nonnegative TOTAL burdens.

running the nitrogen protocol against the boron model therefore fails on ~100%
of samples, which would delete cell 1 -- the one cell that makes the whole
factorial a code-equivalence test.

three fixes were considered.

  (a) let the error stand.  destroys cell 1.  rejected.
  (b) clamp every evaluation at max(x, 0), as `simulate.py` does for its own
      integrator probes.  MEASURED and rejected: applying the clamp symmetrically
      to both arms changes 247-265 of 400 free-limit results at every epsilon
      tested, i.e. the clamp is a large uncontrolled protocol change, not a
      no-op.  it would make cell 2 stop reproducing nitrogen.
  (c) continue the boron field analytically past the boundary.  adopted.

(c) is legitimate because the guard in `solveFreePools` is a DOMAIN ASSERTION,
not a singularity.  the binding residual

    r1 = cf*(1 + nu + uf/kappa_cu + af/kappa_ca) - c_tot
    r2 = df*(1 +      uf/kappa_du + af/kappa_da) - d_tot,   uf = u/su, af = a/sa

is a smooth function of (u, a) in a neighbourhood of the nonnegative orthant, and
its jacobian is nonsingular there, so the implicit function theorem gives a
unique smooth continuation.  what is lost past the boundary is the
knaster-tarski UNIQUENESS certificate (the fixed-point map stops being monotone
for negative substrate), not existence -- so the continuation is solved by
newton from the interior solution rather than by monotone iteration, and it is
only ever used for transient probes.

the discipline that keeps cell 1 honest
---------------------------------------
`solveFreePoolsExtended` DELEGATES to the shipped `model.solveFreePools`
whenever u >= 0 and a >= 0.  every initial guess of the nitrogen protocol is
positive and every converged root it accepts is nonnegative (it rejects
x < -1e-8), so the root, its residual, its jacobian and its eigenvalues all come
from the shipped boron code path.  the continuation is reached only by
intermediate line-search iterates.  `tests/phase2/test_continuation.py` asserts
the delegation is exact on nonnegative states and that the continuation matches
central differences on negative ones.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple

import numpy as np

from proteostasis import model as _bm

_TOL = 1e-13
_MAX_NEWTON = 80


def _isNonnegative(u: float, a: float) -> bool:
    return u >= 0.0 and a >= 0.0


def solveFreePoolsExtended(u: float, a: float, p) -> Tuple[float, float, float, float]:
    """(uf, af, cf, df), continued to negative burdens.

    on the nonnegative orthant this is `model.solveFreePools` verbatim.
    """
    if _isNonnegative(u, a):
        return _bm.solveFreePools(u, a, p)

    # newton on the same 2-d residual, started from the burden-free solution,
    # with no box projection (the box [0, c_tot] is a consequence of nonnegative
    # substrate and does not hold on the continuation).
    cf = p.c_tot / (1.0 + p.nu)
    df = p.d_tot
    (r1, r2), jac, aux = _bm._bindingResidual(cf, df, u, a, p)
    norm = max(abs(r1), abs(r2))
    for _ in range(_MAX_NEWTON):
        if norm <= _TOL * max(1.0, p.c_tot, p.d_tot):
            break
        j11, j12, j21, j22 = jac
        det = j11 * j22 - j12 * j21
        if not np.isfinite(det) or abs(det) < 1e-14:
            raise _bm.ModelError(
                f"singular continued binding jacobian at u={u:.6g}, a={a:.6g}")
        dc = -(j22 * r1 - j12 * r2) / det
        dd = -(-j21 * r1 + j11 * r2) / det
        step, accepted = 1.0, False
        for _ in range(40):
            (t1, t2), jac_try, aux_try = _bm._bindingResidual(
                cf + step * dc, df + step * dd, u, a, p)
            n_try = max(abs(t1), abs(t2))
            if np.isfinite(n_try) and n_try < norm:
                cf, df = cf + step * dc, df + step * dd
                r1, r2, jac, aux, norm = t1, t2, jac_try, aux_try, n_try
                accepted = True
                break
            step *= 0.5
        if not accepted:
            break
    if not np.isfinite(norm) or norm > 1e-8 * max(1.0, p.c_tot, p.d_tot):
        raise _bm.ModelError(
            f"continued free-pool solve failed at u={u:.6g}, a={a:.6g}")
    return aux["uf"], aux["af"], cf, df


def fluxesExtended(u: float, a: float, p) -> Dict[str, float]:
    """the same seven fluxes as `model.fluxes`, on the continued pools.

    the expressions are transcribed from `model.fluxes` unchanged; only the
    free-pool solve differs.  `nucleate` uses uf**m, which for negative uf is
    real only when m is an integer -- the benchmark pins m = 2.0, and anything
    else is refused rather than silently returning nan.
    """
    if _isNonnegative(u, a):
        return _bm.fluxes(u, a, p)
    uf, af, cf, df = solveFreePoolsExtended(u, a, p)
    if uf < 0.0 and float(p.m) != int(p.m):
        raise _bm.ModelError(
            f"continuation to uf<0 needs integer nucleation order; m={p.m}")
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


def rhsVectorExtended(x, p) -> np.ndarray:
    u, a = float(x[0]), float(x[1])
    if _isNonnegative(u, a):
        return _bm.rhsVector((u, a), p)
    f = fluxesExtended(u, a, p)
    du = (f["influx"] - f["refold"] - f["degrade_u"] - f["nucleate"]
          - f["grow"] + f["disaggregate"])
    da = f["nucleate"] + f["grow"] - f["disaggregate"] - f["degrade_a"]
    return np.array([du, da], dtype=float)


def jacobianExtended(u: float, a: float, p) -> np.ndarray:
    """`model.jacobian` transcribed onto the continued free-pool solve.

    every line after the pool solve is identical to `model.jacobian`; the
    derivation (implicit function theorem on the pool balances, then the chain
    rule through the free concentrations) is unchanged because the residual is
    the same smooth function on either side of the boundary.
    """
    if _isNonnegative(u, a):
        return _bm.jacobian(u, a, p)

    uf, af, cf, df = solveFreePoolsExtended(u, a, p)
    if uf < 0.0 and float(p.m) != int(p.m):
        raise _bm.ModelError(
            f"continuation to uf<0 needs integer nucleation order; m={p.m}")
    (_, _), (j11, j12, j21, j22), aux = _bm._bindingResidual(cf, df, u, a, p)
    su, sa = aux["su"], aux["sa"]

    r1_u = cf / (su * p.kappa_cu)
    r1_a = cf / (sa * p.kappa_ca)
    r2_u = df / (su * p.kappa_du)
    r2_a = df / (sa * p.kappa_da)

    det = j11 * j22 - j12 * j21
    if not np.isfinite(det) or abs(det) < 1e-300:
        raise _bm.ModelError("singular continued binding jacobian")
    inv = np.array([[j22, -j12], [-j21, j11]], dtype=float) / det
    dres = np.array([[r1_u, r1_a], [r2_u, r2_a]], dtype=float)
    d_resources = -inv @ dres
    dcf_du, dcf_da = d_resources[0]
    ddf_du, ddf_da = d_resources[1]

    duf_dcf, duf_ddf = aux["duf_dcf"], aux["duf_ddf"]
    daf_dcf, daf_ddf = aux["daf_dcf"], aux["daf_ddf"]

    duf_du = 1.0 / su + duf_dcf * dcf_du + duf_ddf * ddf_du
    duf_da = duf_dcf * dcf_da + duf_ddf * ddf_da
    daf_du = daf_dcf * dcf_du + daf_ddf * ddf_du
    daf_da = 1.0 / sa + daf_dcf * dcf_da + daf_ddf * ddf_da

    ref_den = (p.kappa_ref + uf) ** 2
    degu_den = (p.kappa_u + uf) ** 2
    dis_den = (p.kappa_dis + af) ** 2
    dega_den = (p.kappa_a + af) ** 2
    # model.jacobian writes `if uf > 0 else 0`; on the continuation the honest
    # derivative of alpha_n*uf**m at uf < 0 with integer m is the same formula.
    nuc_prime = p.alpha_n * p.m * (uf ** (p.m - 1.0)) if uf != 0.0 else 0.0

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


def continuationUsed(u: float, a: float) -> bool:
    """True when an evaluation at this state left the shipped code path."""
    return not _isNonnegative(float(u), float(a))
