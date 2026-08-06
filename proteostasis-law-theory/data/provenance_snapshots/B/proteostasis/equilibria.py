"""equilibria, local stability, branch continuation and fold detection.

equilibria are solved in log coordinates (log u, log a). that is legitimate
here because the nonnegative boundary carries no equilibrium when j > 0:

  * at u = 0 we have uf = 0, so du/dt = j + disaggregation > 0;
  * at a = 0 we have af = 0, so da/dt = alpha_n * uf**m, which vanishes only
    when uf = 0, i.e. only when u = 0.

log coordinates therefore lose nothing and enforce positivity for free.

two independent equilibrium finders are provided and cross-checked in
experiment A:
  `findEquilibria`  -- multi-start root finding from a log grid (robust, blind)
  `traceBranch`     -- natural-parameter continuation with warm starts (fast,
                       follows one branch and sees where it ends)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.optimize import root

from .model import Params, ModelError, jacobian, rhsVector, removalCeiling

_LOG_FLOOR = -30.0     # exp(-30) ~ 9e-14
_LOG_CEIL = 12.0       # exp(12)  ~ 1.6e5


@dataclass(frozen=True)
class Equilibrium:
    u: float
    a: float
    eig_real_max: float
    eig_real_min: float
    eig_imag_max: float
    residual: float
    stable: bool

    @property
    def burden(self) -> float:
        return self.u + self.a

    def asDict(self) -> Dict[str, float]:
        return dict(u=self.u, a=self.a, burden=self.burden,
                    eig_real_max=self.eig_real_max, eig_real_min=self.eig_real_min,
                    eig_imag_max=self.eig_imag_max, residual=self.residual,
                    stable=bool(self.stable))


def _clipLog(x: np.ndarray) -> np.ndarray:
    """keep the solver inside a numerically representable burden box.

    hybr explores freely in log space and will otherwise evaluate uf**m at
    exp(700). clipping makes the residual flat outside the box, so the solve
    simply fails there and the candidate is discarded -- which is the intent.
    """
    return np.clip(np.nan_to_num(x, nan=0.0, posinf=_LOG_CEIL, neginf=_LOG_FLOOR),
                   _LOG_FLOOR, _LOG_CEIL)


def _logResidual(x: np.ndarray, p: Params):
    xs = _clipLog(np.asarray(x, dtype=float))
    u, a = float(np.exp(xs[0])), float(np.exp(xs[1]))
    return rhsVector((u, a), p)


def _logJacobian(x: np.ndarray, p: Params):
    xs = _clipLog(np.asarray(x, dtype=float))
    u, a = float(np.exp(xs[0])), float(np.exp(xs[1]))
    return jacobian(u, a, p) * np.array([[u, a], [u, a]])


def residualScale(p: Params) -> float:
    """characteristic flux magnitude, used to make residual tests scale free."""
    return max(p.j, 1e-30, 1e-6 * removalCeiling(p))


def classifyState(u: float, a: float, p: Params, res_tol: float = 1e-8) -> Optional[Equilibrium]:
    """evaluate residual and local stability at a candidate equilibrium."""
    try:
        f = rhsVector((u, a), p)
        jac = jacobian(u, a, p)
    except (ModelError, FloatingPointError, ValueError):
        return None
    res = float(np.max(np.abs(f))) / residualScale(p)
    if not np.isfinite(res) or res > res_tol:
        return None
    eig = np.linalg.eigvals(jac)
    if not np.all(np.isfinite(eig)):
        return None
    return Equilibrium(u=u, a=a,
                       eig_real_max=float(np.max(eig.real)),
                       eig_real_min=float(np.min(eig.real)),
                       eig_imag_max=float(np.max(np.abs(eig.imag))),
                       residual=res,
                       stable=bool(np.max(eig.real) < 0.0))


def solveEquilibriumFrom(u0: float, a0: float, p: Params,
                         res_tol: float = 1e-8) -> Optional[Equilibrium]:
    """single newton/hybrid solve in log coordinates from one starting point."""
    x0 = np.array([np.log(max(u0, 1e-13)), np.log(max(a0, 1e-13))])
    x0 = np.clip(x0, _LOG_FLOOR, _LOG_CEIL)
    try:
        sol = root(_logResidual, x0, args=(p,), jac=_logJacobian, method="hybr",
                   options=dict(xtol=1e-13))
    except (ModelError, ValueError, FloatingPointError):
        return None
    if not np.all(np.isfinite(sol.x)):
        return None
    xs = _clipLog(sol.x)
    # a solution pinned to the clipping box is an artefact of the box, not a root
    if np.any(np.abs(xs - np.asarray(sol.x, dtype=float)) > 1e-9):
        return None
    return classifyState(float(np.exp(xs[0])), float(np.exp(xs[1])), p, res_tol=res_tol)


def findEquilibria(p: Params, n_grid: int = 9, lo: float = 1e-7, hi: float = 1e4,
                   res_tol: float = 1e-8, dedupe_rtol: float = 1e-5) -> List[Equilibrium]:
    """blind multi-start search over a logarithmic grid of initial guesses.

    deliberately independent of continuation so the two can be cross-validated.
    """
    grid = np.logspace(np.log10(lo), np.log10(hi), n_grid)
    found: List[Equilibrium] = []
    for u0 in grid:
        for a0 in grid:
            eq = solveEquilibriumFrom(u0, a0, p, res_tol=res_tol)
            if eq is None:
                continue
            dup = False
            for e in found:
                du = abs(e.u - eq.u) / max(abs(e.u), abs(eq.u), 1e-30)
                da = abs(e.a - eq.a) / max(abs(e.a), abs(eq.a), 1e-30)
                if du < dedupe_rtol and da < dedupe_rtol:
                    dup = True
                    break
            if not dup:
                found.append(eq)
    found.sort(key=lambda e: e.burden)
    return found


def lowestStableEquilibrium(p: Params, **kw) -> Optional[Equilibrium]:
    eqs = [e for e in findEquilibria(p, **kw) if e.stable]
    return eqs[0] if eqs else None


# ---------------------------------------------------------------------------
# continuation
# ---------------------------------------------------------------------------


@dataclass
class BranchPoint:
    value: float
    eq: Optional[Equilibrium]


def traceBranch(base: Params, param: str, values: Sequence[float],
                start: Optional[Tuple[float, float]] = None,
                res_tol: float = 1e-8, confirm_loss: bool = True) -> List[BranchPoint]:
    """natural-parameter continuation of the low stable branch along `param`.

    each step uses adaptive substepping (`_advance`) rather than a single
    warm-started solve. that matters: a naive one-shot warm start fails
    routinely when the branch moves by an order of magnitude between grid
    points, and reports a branch loss that experiment A's blind cross-check
    then contradicts.

    when `confirm_loss` is set, an apparent loss is re-tested with the
    independent blind finder before being recorded as absence of an
    equilibrium.
    """
    out: List[BranchPoint] = []
    prev_v: Optional[float] = None
    prev_eq: Optional[Equilibrium] = None
    if start is not None:
        prev_eq = classifyState(start[0], start[1], base.with_(
            **{param: float(values[0])}), res_tol=res_tol)
        prev_v = float(values[0])

    for v in values:
        v = float(v)
        eq = None
        if prev_eq is not None and prev_v is not None:
            eq = _advance(base, param, prev_v, prev_eq, v, res_tol)
        if eq is None:
            try:
                p = base.with_(**{param: v}).validate()
            except ModelError:
                p = None
            if p is not None and (prev_eq is None or confirm_loss):
                eq = lowestStableEquilibrium(p)
        out.append(BranchPoint(value=v, eq=eq))
        if eq is not None:
            prev_v, prev_eq = v, eq
    return out


@dataclass
class FoldResult:
    param: str
    found: bool
    fold_value: Optional[float]
    bracket: Optional[Tuple[float, float]]
    eq_at_fold: Optional[Equilibrium]
    reason: str
    n_steps: int

    def asDict(self) -> Dict:
        d = dict(param=self.param, found=bool(self.found), fold_value=self.fold_value,
                 bracket_lo=None if self.bracket is None else self.bracket[0],
                 bracket_hi=None if self.bracket is None else self.bracket[1],
                 reason=self.reason, n_steps=self.n_steps)
        if self.eq_at_fold is not None:
            d.update({f"fold_{k}": v for k, v in self.eq_at_fold.asDict().items()})
        return d


#: multiplicative offsets applied to a failed warm start before giving up.
#: cheap insurance against hybr stalling on a perfectly healthy branch.
_WARM_LADDER = ((1.0, 1.0), (2.0, 1.0), (0.5, 1.0), (1.0, 2.0), (1.0, 0.5),
                (5.0, 5.0), (0.2, 0.2))


def stableEquilibriumAt(base: Params, param: str, value: float,
                        warm: Optional[Tuple[float, float]],
                        res_tol: float = 1e-8, blind: bool = False,
                        blind_grid: int = 9) -> Optional[Equilibrium]:
    """find the low stable equilibrium at one parameter value.

    escalating effort: warm start, then a ladder of perturbed warm starts, then
    (only if `blind`) an independent multi-start grid search. the blind stage is
    what prevents a solver stall from being mistaken for a fold, but it costs
    ~n_grid**2 root solves, so it is reserved for the moment a branch is about
    to be declared lost.
    """
    try:
        p = base.with_(**{param: float(value)}).validate()
    except ModelError:
        return None
    if warm is not None:
        for su, sa in _WARM_LADDER:
            eq = solveEquilibriumFrom(warm[0] * su, warm[1] * sa, p, res_tol=res_tol)
            if eq is not None and eq.stable:
                return eq
    if blind:
        return lowestStableEquilibrium(p, n_grid=blind_grid)
    return None


def _advance(base: Params, param: str, v_from: float, eq_from: Equilibrium,
             v_to: float, res_tol: float, depth: int = 0,
             max_depth: int = 6) -> Optional[Equilibrium]:
    """continue from v_from to v_to, subdividing the step on failure.

    subdivision is what lets a coarse outer march coexist with a branch that
    moves by orders of magnitude in a single nominal step.
    """
    eq = stableEquilibriumAt(base, param, v_to, (eq_from.u, eq_from.a), res_tol)
    if eq is not None or depth >= max_depth:
        return eq
    mid = 0.5 * (v_from + v_to)
    eq_mid = _advance(base, param, v_from, eq_from, mid, res_tol, depth + 1, max_depth)
    if eq_mid is None:
        return None
    return _advance(base, param, mid, eq_mid, v_to, res_tol, depth + 1, max_depth)


def findFold(base: Params, param: str, lo: float, hi: float,
             n_march: int = 24, n_bisect: int = 40, res_tol: float = 1e-8,
             scale: str = "log", tol_rel: float = 1e-6,
             blind_gap_rel: float = 1e-2, blind_grid: int = 7) -> FoldResult:
    """locate the loss of the low-burden stable branch along `param`.

    (1) confirm a stable low-burden equilibrium at `lo` with the blind
        multi-start finder;
    (2) march upward with warm starts and adaptive substepping;
    (3) before declaring the branch lost, require an INDEPENDENT blind
        multi-start search at that value to also find nothing -- a solver stall
        must not be reported as a fold;
    (4) bisect the resulting bracket.

    the returned value is the largest parameter value at which a stable
    low-burden equilibrium could still be located. that is a saddle-node only
    if the branch ends by colliding with the upper unstable equilibrium;
    experiments B and C additionally record the equilibrium count either side,
    so the interpretation is checked rather than assumed.
    """
    if not (hi > lo > 0.0) and scale == "log":
        scale = "linear"
    try:
        p_lo = base.with_(**{param: float(lo)}).validate()
    except ModelError as exc:
        return FoldResult(param, False, None, None, None, f"invalid params: {exc}", 0)

    eq0 = lowestStableEquilibrium(p_lo)
    if eq0 is None:
        return FoldResult(param, False, None, None, None, "no stable low branch at lo", 0)

    if scale == "log":
        grid = np.logspace(np.log10(lo), np.log10(hi), n_march + 1)
    else:
        grid = np.linspace(lo, hi, n_march + 1)

    steps = 0
    last_ok_v, last_ok_eq = float(lo), eq0
    first_bad_v = None
    for v in grid[1:]:
        eq = _advance(base, param, last_ok_v, last_ok_eq, float(v), res_tol)
        steps += 1
        if eq is None:
            # independent confirmation before calling this a fold
            eq = stableEquilibriumAt(base, param, float(v), None, res_tol,
                                     blind=True, blind_grid=blind_grid)
            if eq is None:
                first_bad_v = float(v)
                break
        last_ok_v, last_ok_eq = float(v), eq

    if first_bad_v is None:
        return FoldResult(param, False, None, (last_ok_v, float(hi)), last_ok_eq,
                          "branch survives to hi", steps)

    a_, b_ = last_ok_v, first_bad_v
    eq_ok = last_ok_eq
    for _ in range(n_bisect):
        if (b_ - a_) <= tol_rel * max(1e-12, abs(b_)):
            break
        mid = 0.5 * (a_ + b_)
        eq = _advance(base, param, a_, eq_ok, mid, res_tol, max_depth=5)
        if eq is None and (b_ - a_) > blind_gap_rel * max(1e-12, abs(b_)):
            # while the bracket is still wide, a stalled warm start is more
            # likely than a genuine branch end; pay for the blind check
            eq = stableEquilibriumAt(base, param, mid, None, res_tol,
                                     blind=True, blind_grid=blind_grid)
        steps += 1
        if eq is not None:
            a_, eq_ok = mid, eq
        else:
            b_ = mid

    return FoldResult(param, True, a_, (a_, b_), eq_ok, "branch terminates", steps)


def eigenvalueAlongBranch(base: Params, param: str, values: Sequence[float],
                          start: Optional[Tuple[float, float]] = None) -> List[Dict]:
    """leading eigenvalue along a continued branch.

    used to test the critical-slowing-down prediction (PREDICTIONS.md #5): the
    leading eigenvalue should approach zero as the branch approaches its fold.
    """
    rows = []
    for bp in traceBranch(base, param, values, start=start):
        if bp.eq is None:
            rows.append(dict(value=bp.value, exists=False, u=np.nan, a=np.nan,
                             burden=np.nan, eig_real_max=np.nan, stable=False))
        else:
            rows.append(dict(value=bp.value, exists=True, u=bp.eq.u, a=bp.eq.a,
                             burden=bp.eq.burden, eig_real_max=bp.eq.eig_real_max,
                             stable=bp.eq.stable))
    return rows
