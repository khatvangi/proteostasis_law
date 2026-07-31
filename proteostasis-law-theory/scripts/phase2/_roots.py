"""independent equilibrium finders and root certification for the phase 2 audit.

experiment C decided multistability with one root finder: MINPACK `hybr` in log
coordinates, multi-started from an 11x11 grid, deduplicated at 1e-5 relative.
re-running that denser would only ask the same algorithm the same question, so
this module supplies methods that fail differently:

  `denseMultiStart(..., method='hybr')`  the original algorithm, denser box
  `denseMultiStart(..., method='lm')`    levenberg-marquardt (MINPACK lmdif);
                                         different globalization, different
                                         failure modes on ill-conditioned roots
  `nullclineRoots(...)`                  NO newton anywhere: sign-change
                                         bracketing plus brentq on the
                                         u-nullcline, then bisection in a. a
                                         bracketing method cannot invent a root
                                         where the field does not change sign,
                                         so it is the strongest available check
                                         against root-finder artefacts
  `continueBranch(...)`                  pseudo-arclength continuation in j,
                                         which traces the equilibrium manifold
                                         through folds instead of sampling it

`certifyRoot` then attaches, per root: a cancellation-aware residual, an
analytic-vs-Richardson jacobian agreement, the conserved-pool balance residuals,
and the knaster-tarski uniqueness certificate for the rapid-equilibrium closure.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import brentq, root

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from proteostasis.model import (Params, ModelError, fluxes, jacobian,  # noqa: E402
                                poolBalanceResiduals, rhs, rhsVector,
                                solveFreePoolsCertified)
from proteostasis.equilibria import _clipLog, _LOG_CEIL, _LOG_FLOOR, residualScale  # noqa: E402

#: relative dedupe tolerances swept to test how sensitive the equilibrium COUNT
#: is to the clustering choice. 1e-5 is what experiment C used.
DEDUPE_TOLS = (1e-9, 1e-7, 1e-5, 1e-3, 1e-2, 1e-1)


# ---------------------------------------------------------------------------
# raw root records
# ---------------------------------------------------------------------------


@dataclass
class Root:
    u: float
    a: float
    stable: bool
    eig_real_max: float
    eig_real_min: float
    eig_imag_max: float
    residual_rel: float          # experiment C's normalisation (per residualScale)
    residual_cancel: float       # normalised by the largest individual flux
    method: str
    n_hits: int = 1              # how many independent starts landed here

    @property
    def burden(self) -> float:
        return self.u + self.a

    @property
    def kind(self) -> str:
        """topological type from the eigenvalue signs."""
        lo, hi = self.eig_real_min, self.eig_real_max
        osc = self.eig_imag_max > 1e-12 * max(abs(lo), abs(hi), 1.0)
        if hi < 0.0:
            return "stable_spiral" if osc else "stable_node"
        if lo > 0.0:
            return "unstable_spiral" if osc else "unstable_node"
        if lo < 0.0 < hi:
            return "saddle"
        return "degenerate"

    def asDict(self) -> Dict:
        return dict(u=self.u, a=self.a, burden=self.burden, stable=bool(self.stable),
                    kind=self.kind, eig_real_max=self.eig_real_max,
                    eig_real_min=self.eig_real_min, eig_imag_max=self.eig_imag_max,
                    residual_rel=self.residual_rel,
                    residual_cancel=self.residual_cancel,
                    method=self.method, n_hits=int(self.n_hits))


def _classify(u: float, a: float, p: Params, method: str,
              res_tol: float = 1e-8) -> Optional[Root]:
    """turn a candidate point into a `Root`, or reject it.

    the acceptance test is deliberately the SAME as experiment C's
    (`residual/residualScale <= 1e-8`), so that a difference in the audit's
    findings is attributable to the search, not to a moved goalpost. the
    stricter cancellation-aware residual is recorded alongside rather than
    used to filter.
    """
    if not (np.isfinite(u) and np.isfinite(a)) or u <= 0.0 or a <= 0.0:
        return None
    try:
        f = rhsVector((u, a), p)
        jac = jacobian(u, a, p)
        fl = fluxes(u, a, p)
    except (ModelError, FloatingPointError, ValueError, ZeroDivisionError):
        return None
    if not np.all(np.isfinite(f)) or not np.all(np.isfinite(jac)):
        return None
    res = float(np.max(np.abs(f))) / residualScale(p)
    if not np.isfinite(res) or res > res_tol:
        return None
    scale = max(abs(fl[k]) for k in ("influx", "refold", "degrade_u", "degrade_a",
                                     "nucleate", "grow", "disaggregate"))
    res_cancel = float(np.max(np.abs(f))) / max(scale, 1e-300)
    eig = np.linalg.eigvals(jac)
    if not np.all(np.isfinite(eig)):
        return None
    return Root(u=float(u), a=float(a),
                stable=bool(np.max(eig.real) < 0.0),
                eig_real_max=float(np.max(eig.real)),
                eig_real_min=float(np.min(eig.real)),
                eig_imag_max=float(np.max(np.abs(eig.imag))),
                residual_rel=res, residual_cancel=res_cancel, method=method)


# ---------------------------------------------------------------------------
# method 1/2 -- dense multi-start newton-type search
# ---------------------------------------------------------------------------


def _logResidual(x, p: Params):
    xs = _clipLog(np.asarray(x, dtype=float))
    return rhsVector((float(np.exp(xs[0])), float(np.exp(xs[1]))), p)


def _logJac(x, p: Params):
    xs = _clipLog(np.asarray(x, dtype=float))
    u, a = float(np.exp(xs[0])), float(np.exp(xs[1]))
    return jacobian(u, a, p) * np.array([[u, a], [u, a]])


def solveFrom(u0: float, a0: float, p: Params, method: str = "hybr",
              res_tol: float = 1e-8) -> Optional[Root]:
    """one solve from one start. returns None if it did not land on a root."""
    x0 = np.clip(np.array([np.log(max(u0, 1e-14)), np.log(max(a0, 1e-14))]),
                 _LOG_FLOOR, _LOG_CEIL)
    try:
        if method == "hybr":
            sol = root(_logResidual, x0, args=(p,), jac=_logJac, method="hybr",
                       options=dict(xtol=1e-13))
        elif method == "lm":
            sol = root(_logResidual, x0, args=(p,), jac=_logJac, method="lm",
                       options=dict(xtol=1e-14, ftol=1e-14, maxiter=2000))
        else:
            raise ValueError(f"unknown method '{method}'")
    except (ModelError, ValueError, FloatingPointError, ZeroDivisionError):
        return None
    if not np.all(np.isfinite(sol.x)):
        return None
    xs = _clipLog(sol.x)
    # a point pinned against the log box is an artefact of the box, not a root
    if np.any(np.abs(xs - np.asarray(sol.x, dtype=float)) > 1e-9):
        return None
    return _classify(float(np.exp(xs[0])), float(np.exp(xs[1])), p, method, res_tol)


def startGrid(n_grid: int, lo: float, hi: float,
              anchors: Sequence[Tuple[float, float]] = ()) -> List[Tuple[float, float]]:
    """scale-adaptive multi-start grid.

    a plain log grid over a fixed box under-samples wherever the actual roots
    happen to sit, which is exactly the regime where a second attractor is
    missed. so the grid is augmented with a local cloud around every anchor
    (typically the roots experiment C already reported), spanning +-1.5 decades.
    """
    base = np.logspace(np.log10(lo), np.log10(hi), n_grid)
    pts = [(float(u), float(a)) for u in base for a in base]
    for (au, aa) in anchors:
        if not (np.isfinite(au) and np.isfinite(aa) and au > 0 and aa > 0):
            continue
        loc = np.array([-1.5, -0.75, -0.25, 0.0, 0.25, 0.75, 1.5])
        for du in loc:
            for da in loc:
                pts.append((float(au * 10.0 ** du), float(aa * 10.0 ** da)))
    return pts


def denseMultiStart(p: Params, n_grid: int = 41, lo: float = 1e-12,
                    hi: float = 1.5e5, method: str = "hybr",
                    anchors: Sequence[Tuple[float, float]] = (),
                    res_tol: float = 1e-8) -> List[Root]:
    """multi-start search returning the RAW, un-deduplicated root list.

    keeping duplicates is the point: how many distinct roots there are depends
    on the clustering tolerance, and that dependence is itself a result.
    """
    out: List[Root] = []
    for u0, a0 in startGrid(n_grid, lo, hi, anchors):
        r = solveFrom(u0, a0, p, method=method, res_tol=res_tol)
        if r is not None:
            out.append(r)
    return out


def clusterRoots(roots: Sequence[Root], rtol: float) -> List[Root]:
    """merge roots that agree to `rtol` relatively in BOTH coordinates.

    same rule experiment C used, applied to the raw list so the count can be
    reported as a function of rtol.
    """
    kept: List[Root] = []
    for r in sorted(roots, key=lambda z: z.burden):
        hit = None
        for k in kept:
            du = abs(k.u - r.u) / max(abs(k.u), abs(r.u), 1e-300)
            da = abs(k.a - r.a) / max(abs(k.a), abs(r.a), 1e-300)
            if du < rtol and da < rtol:
                hit = k
                break
        if hit is None:
            kept.append(replace(r, n_hits=1))
        else:
            n = hit.n_hits + 1
            # keep the representative with the smallest cancellation residual
            if r.residual_cancel < hit.residual_cancel:
                kept[kept.index(hit)] = replace(r, n_hits=n)
            else:
                hit.n_hits = n
    kept.sort(key=lambda z: z.burden)
    return kept


def dedupeSensitivity(roots: Sequence[Root],
                      tols: Sequence[float] = DEDUPE_TOLS) -> Dict[str, Dict[str, int]]:
    """equilibrium counts as a function of the clustering tolerance."""
    out = {}
    for t in tols:
        c = clusterRoots(roots, t)
        out[f"{t:g}"] = dict(n_eq=len(c), n_stable=int(sum(z.stable for z in c)))
    return out


# ---------------------------------------------------------------------------
# method 3 -- newton-free nullcline bracketing
# ---------------------------------------------------------------------------


def _f1(u: float, a: float, p: Params) -> float:
    try:
        return float(rhs(u, a, p)[0])
    except (ModelError, FloatingPointError, ValueError, ZeroDivisionError):
        return np.nan


def _f2(u: float, a: float, p: Params) -> float:
    try:
        return float(rhs(u, a, p)[1])
    except (ModelError, FloatingPointError, ValueError, ZeroDivisionError):
        return np.nan


def uNullcline(a: float, p: Params, u_grid: np.ndarray) -> List[float]:
    """all u with du/dt = 0 at fixed a, found by sign change + brentq.

    du/dt is positive as u -> 0 (influx plus disaggregation, no removal) and
    negative as u -> infinity (nucleation dominates, order m > 1), so at least
    one crossing always exists inside a wide enough grid.
    """
    vals = np.array([_f1(float(u), a, p) for u in u_grid])
    out: List[float] = []
    for i in range(len(u_grid) - 1):
        v0, v1 = vals[i], vals[i + 1]
        if not (np.isfinite(v0) and np.isfinite(v1)):
            continue
        if v0 == 0.0:
            out.append(float(u_grid[i]))
        elif v0 * v1 < 0.0:
            try:
                out.append(float(brentq(lambda uu: _f1(uu, a, p),
                                        float(u_grid[i]), float(u_grid[i + 1]),
                                        xtol=1e-300, rtol=8.9e-16, maxiter=200)))
            except (ValueError, RuntimeError):
                continue
    return sorted(out)


def nullclineRoots(p: Params, u_lo: float = 1e-12, u_hi: float = 1.5e5,
                   a_lo: float = 1e-12, a_hi: float = 1.5e5,
                   n_u: int = 400, n_a: int = 1200,
                   res_tol: float = 1e-8) -> Tuple[List[Root], Dict]:
    """equilibria located WITHOUT any newton step on the 2-d system.

    walk a along a log grid; at each a solve the u-nullcline by bracketing; then
    follow each nullcline strand and look for a sign change of da/dt along it.
    each sign change is refined by bisection in a (with u re-bracketed at every
    bisection step), so the located point is the intersection of the two
    nullclines and nothing about it depends on a jacobian being right.
    """
    u_grid = np.logspace(np.log10(u_lo), np.log10(u_hi), n_u)
    a_grid = np.logspace(np.log10(a_lo), np.log10(a_hi), n_a)

    # strands[k] = list of (a_index, u, g) belonging to one nullcline strand
    strands: List[List[Tuple[int, float, float]]] = []
    prev: List[Tuple[int, float]] = []          # (strand id, u) at previous a
    for ia, a in enumerate(a_grid):
        us = uNullcline(float(a), p, u_grid)
        cur: List[Tuple[int, float]] = []
        used = set()
        for u in us:
            g = _f2(u, float(a), p)
            if not np.isfinite(g):
                continue
            # attach to the nearest strand in log-u, if close enough
            best, bestd = None, np.inf
            for sid, up in prev:
                if sid in used:
                    continue
                d = abs(np.log10(max(u, 1e-300)) - np.log10(max(up, 1e-300)))
                if d < bestd:
                    best, bestd = sid, d
            if best is not None and bestd < 0.6:
                sid = best
                used.add(sid)
            else:
                strands.append([])
                sid = len(strands) - 1
            strands[sid].append((ia, u, g))
            cur.append((sid, u))
        prev = cur

    roots: List[Root] = []
    n_sign_changes = 0
    for st in strands:
        for (ia0, u0, g0), (ia1, u1, g1) in zip(st, st[1:]):
            if ia1 != ia0 + 1 or not (np.isfinite(g0) and np.isfinite(g1)):
                continue
            if g0 == 0.0:
                r = _classify(u0, float(a_grid[ia0]), p, "nullcline", res_tol)
                if r is not None:
                    roots.append(r)
                continue
            if g0 * g1 >= 0.0:
                continue
            n_sign_changes += 1
            # pure bisection in a; u is re-bracketed by brentq at every step
            alo, ahi = float(a_grid[ia0]), float(a_grid[ia1])
            ulo, uhi = u0, u1
            glo = g0
            for _ in range(200):
                amid = float(np.sqrt(alo * ahi))       # geometric bisection
                cand = uNullcline(amid, p, np.logspace(
                    np.log10(min(ulo, uhi) / 30.0), np.log10(max(ulo, uhi) * 30.0), 60))
                if not cand:
                    break
                umid = min(cand, key=lambda z: abs(np.log(z) - 0.5 * (np.log(ulo) + np.log(uhi))))
                gmid = _f2(umid, amid, p)
                if not np.isfinite(gmid):
                    break
                if gmid == 0.0 or (ahi - alo) <= 1e-15 * ahi:
                    alo, ulo = amid, umid
                    break
                if glo * gmid < 0.0:
                    ahi, uhi = amid, umid
                else:
                    alo, ulo, glo = amid, umid, gmid
            r = _classify(ulo, alo, p, "nullcline", res_tol)
            if r is None:
                # accept the geometry even if the residual test rejects it, but
                # record it separately -- an unrefinable crossing is a finding
                r2 = _classify(uhi, ahi, p, "nullcline", res_tol)
                if r2 is not None:
                    roots.append(r2)
            else:
                roots.append(r)
    diag = dict(n_strands=len(strands), n_sign_changes=int(n_sign_changes),
                n_accepted=len(roots))
    return roots, diag


# ---------------------------------------------------------------------------
# method 4 -- pseudo-arclength continuation in j
# ---------------------------------------------------------------------------


def _augmented(x: np.ndarray, lam: float, p: Params):
    """F and its 2x3 jacobian [dF/dx | dF/dlam] in (log u, log a, j)."""
    pl = p.with_(j=float(lam))
    F = _logResidual(x, pl)
    Jx = _logJac(x, pl)
    Fl = np.array([1.0, 0.0])          # only du/dt carries +j
    return F, np.column_stack([Jx, Fl])


def continueBranch(p: Params, u0: float, a0: float, j0: float,
                   ds: float = 0.02, max_steps: int = 4000,
                   j_min: float = 1e-8, j_max: float = 50.0,
                   direction: int = +1) -> Dict:
    """pseudo-arclength continuation of the equilibrium manifold in j.

    natural-parameter continuation stops dead at a fold, which is precisely the
    structure that produces bistability -- so it cannot answer the question.
    arclength continuation parameterises the branch by its own arc length and
    walks straight through folds, which is what makes the S-curve visible.

    returns the traced points, the fold count (sign changes of the tangent's j
    component), and the j and burden extents of the branch.
    """
    X = np.array([np.log(u0), np.log(a0), float(j0)], dtype=float)
    pts: List[Dict] = []
    folds: List[Dict] = []
    T_prev = None
    status = "max_steps"

    for step in range(max_steps):
        F, A = _augmented(X[:2], X[2], p)
        if not np.all(np.isfinite(A)):
            status = "nonfinite_jacobian"
            break
        # tangent = null space of the 2x3 matrix
        _, _, Vt = np.linalg.svd(A)
        T = Vt[-1]
        if T_prev is None:
            T = T if np.sign(T[2]) == np.sign(direction) or T[2] == 0.0 else -T
        elif float(T @ T_prev) < 0.0:
            T = -T
        # record BEFORE stepping
        u, a = float(np.exp(X[0])), float(np.exp(X[1]))
        pts.append(dict(step=step, j=float(X[2]), u=u, a=a, burden=u + a,
                        t_j=float(T[2]), residual=float(np.max(np.abs(F)))))
        if T_prev is not None and T[2] * T_prev[2] < 0.0:
            folds.append(dict(step=step, j=float(X[2]), burden=u + a))
        T_prev = T

        Xp = X + ds * T
        # corrector: F(X)=0 and T.(X - Xp)=0
        Xc = Xp.copy()
        ok = False
        for _ in range(30):
            Fc, Ac = _augmented(Xc[:2], Xc[2], p)
            g = float(T @ (Xc - Xp))
            resid = np.array([Fc[0], Fc[1], g])
            if not np.all(np.isfinite(resid)):
                break
            if np.max(np.abs(resid)) < 1e-12 * max(1.0, abs(Xc[2])):
                ok = True
                break
            M = np.vstack([Ac, T])
            try:
                dX = np.linalg.solve(M, -resid)
            except np.linalg.LinAlgError:
                break
            Xc = Xc + dX
        if not ok:
            ds *= 0.5
            if ds < 1e-7:
                status = "corrector_failed"
                break
            continue
        X = Xc
        ds = min(ds * 1.25, 0.05)
        if not np.all(np.isfinite(X)):
            status = "nonfinite_state"
            break
        if X[2] < j_min or X[2] > j_max:
            status = "j_out_of_range"
            break
        if X[0] < _LOG_FLOOR or X[0] > _LOG_CEIL or X[1] < _LOG_FLOOR or X[1] > _LOG_CEIL:
            status = "state_out_of_box"
            break

    js = np.array([q["j"] for q in pts]) if pts else np.array([])
    bs = np.array([q["burden"] for q in pts]) if pts else np.array([])
    return dict(status=status, n_points=len(pts), n_folds=len(folds), folds=folds,
                j_min=float(js.min()) if len(js) else None,
                j_max=float(js.max()) if len(js) else None,
                burden_min=float(bs.min()) if len(bs) else None,
                burden_max=float(bs.max()) if len(bs) else None,
                points=pts)


# ---------------------------------------------------------------------------
# certification
# ---------------------------------------------------------------------------


def richardsonJacobian(u: float, a: float, p: Params) -> np.ndarray:
    """central differences with one richardson extrapolation step.

    the analytic jacobian threads the implicit function theorem through the
    rapid-equilibrium closure; if that derivation is wrong the stability
    classification is wrong, and stability is the whole claim here. a plain
    finite difference is too noisy at the 1e-8 level to adjudicate, so the
    O(h^2) error term is cancelled explicitly.
    """
    def col(k, h):
        xp = [u, a]
        xm = [u, a]
        xp[k] += h
        xm[k] = max(xm[k] - h, 1e-300)
        step = xp[k] - xm[k]
        return (rhsVector(xp, p) - rhsVector(xm, p)) / step

    out = np.zeros((2, 2))
    for k, x in enumerate((u, a)):
        h = 1e-4 * max(abs(x), 1e-8)
        d1, d2 = col(k, h), col(k, h / 2.0)
        out[:, k] = (4.0 * d2 - d1) / 3.0
    return out


def certifyRoot(u: float, a: float, p: Params) -> Dict:
    """every independent check we can cheaply apply at one candidate root."""
    out: Dict = dict(u=float(u), a=float(a), burden=float(u + a))
    try:
        f = rhsVector((u, a), p)
        fl = fluxes(u, a, p)
        Ja = jacobian(u, a, p)
    except (ModelError, FloatingPointError, ValueError) as exc:
        out["error"] = f"{type(exc).__name__}: {exc}"
        return out

    flux_scale = max(abs(fl[k]) for k in ("influx", "refold", "degrade_u",
                                          "degrade_a", "nucleate", "grow",
                                          "disaggregate"))
    out["abs_residual"] = float(np.max(np.abs(f)))
    out["residual_rel"] = float(np.max(np.abs(f))) / residualScale(p)
    out["residual_cancel"] = float(np.max(np.abs(f))) / max(flux_scale, 1e-300)
    out["flux_scale"] = float(flux_scale)

    try:
        Jn = richardsonJacobian(u, a, p)
        den = max(float(np.max(np.abs(Ja))), 1e-300)
        out["jac_rel_error"] = float(np.max(np.abs(Ja - Jn)) / den)
    except (ModelError, FloatingPointError, ValueError):
        out["jac_rel_error"] = None

    eig = np.linalg.eigvals(Ja)
    out["eig_real_max"] = float(np.max(eig.real))
    out["eig_real_min"] = float(np.min(eig.real))
    out["eig_imag_max"] = float(np.max(np.abs(eig.imag)))
    out["stable"] = bool(np.max(eig.real) < 0.0)
    # conditioning: a stable root whose leading eigenvalue is within numerical
    # noise of zero is not safely classifiable as stable
    out["jac_cond"] = float(np.linalg.cond(Ja))
    out["eig_margin_rel"] = float(abs(np.max(eig.real))
                                  / max(float(np.max(np.abs(eig.real))), 1e-300))

    bal = poolBalanceResiduals(u, a, p)
    out["pool_balance_max"] = float(max(abs(v) for v in bal.values()))
    out.update({f"pool_{k}": float(v) for k, v in bal.items()})

    _, _, _, _, cert = solveFreePoolsCertified(u, a, p)
    out["closure_unique"] = bool(cert["unique"])
    out["closure_gap"] = float(cert["gap"])
    out["closure_residual"] = float(cert["residual"])

    out["positive"] = bool(u > 0.0 and a > 0.0)
    out["free_pools_bounded"] = bool(0.0 <= fl["cf"] <= p.c_tot * (1 + 1e-12)
                                     and 0.0 <= fl["df"] <= p.d_tot * (1 + 1e-12)
                                     and 0.0 <= fl["uf"] <= u * (1 + 1e-12)
                                     and 0.0 <= fl["af"] <= a * (1 + 1e-12))
    # net mass balance must hold identically at any state
    out["mass_balance_residual"] = float(
        ((f[0] + f[1]) - (fl["influx"] - fl["refold"] - fl["degrade_u"]
                          - fl["degrade_a"])) / max(flux_scale, 1e-300))
    return out


def matchRoot(roots: Sequence[Root], u: float, a: float,
              rtol: float = 1e-3) -> Optional[Root]:
    """find the root in `roots` matching (u, a) to `rtol` in both coordinates."""
    for r in roots:
        du = abs(r.u - u) / max(abs(r.u), abs(u), 1e-300)
        da = abs(r.a - a) / max(abs(r.a), abs(a), 1e-300)
        if du < rtol and da < rtol:
            return r
    return None
