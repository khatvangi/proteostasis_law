"""trajectory integration with positivity monitoring and blow-up detection.

the integrator is stiff-capable (Radau, analytic jacobian) because the model
mixes fast resource re-equilibration with slow burden accumulation. two events
terminate a run:

  * `blowup`  -- total burden exceeds a ceiling; the trajectory has escaped and
                 no bounded state is reachable from this initial condition.
  * `negative`-- either state falls below a small negative tolerance. this must
                 never fire: the nonnegative orthant is forward invariant, so a
                 firing means either an rhs sign error or an integrator
                 tolerance failure, and the tests treat it as a hard failure.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from scipy.integrate import solve_ivp

from .model import Params, ModelError, jacobian, rhsVector

NEG_TOL = -1e-9        # states below this are treated as leakage, not round-off


@dataclass
class Trajectory:
    t: np.ndarray
    u: np.ndarray
    a: np.ndarray
    status: str            # "converged" | "blowup" | "negative" | "timeout" | "error"
    min_u: float
    min_a: float
    final_u: float
    final_a: float
    final_rate: float      # max |d/dt| at the final point, relative to influx scale
    message: str = ""

    @property
    def burden(self) -> np.ndarray:
        return self.u + self.a

    def summary(self) -> Dict:
        return dict(status=self.status, min_u=self.min_u, min_a=self.min_a,
                    final_u=self.final_u, final_a=self.final_a,
                    final_burden=self.final_u + self.final_a,
                    final_rate=self.final_rate, t_end=float(self.t[-1]),
                    message=self.message)


def simulate(p: Params, u0: float, a0: float, t_end: float = 5.0e4,
             n_out: int = 400, blowup: float = 1.0e6,
             rtol: float = 1e-9, atol: float = 1e-12) -> Trajectory:
    """integrate from one initial condition."""
    p.validate()
    if u0 < 0.0 or a0 < 0.0:
        raise ModelError("initial condition must be nonnegative")

    def f(_t, y):
        # clip only for the solver's internal probing; leakage is still recorded
        return rhsVector((max(y[0], 0.0), max(y[1], 0.0)), p)

    def jac(_t, y):
        return jacobian(max(y[0], 0.0), max(y[1], 0.0), p)

    def evBlow(_t, y):
        return blowup - (y[0] + y[1])
    evBlow.terminal = True
    evBlow.direction = -1.0

    def evNeg(_t, y):
        return min(y[0], y[1]) - NEG_TOL
    evNeg.terminal = True
    evNeg.direction = -1.0

    t_eval = np.linspace(0.0, t_end, n_out)
    try:
        sol = solve_ivp(f, (0.0, t_end), [float(u0), float(a0)], method="Radau",
                        jac=jac, t_eval=t_eval, events=(evBlow, evNeg),
                        rtol=rtol, atol=atol, dense_output=False)
    except (ModelError, ValueError, FloatingPointError) as exc:
        z = np.array([0.0])
        return Trajectory(z, np.array([u0]), np.array([a0]), "error",
                          u0, a0, u0, a0, np.nan, str(exc))

    t = np.concatenate([sol.t, np.concatenate(sol.t_events)]) if any(
        len(e) for e in sol.t_events) else sol.t
    ys = sol.y
    if any(len(e) for e in sol.y_events):
        extra = np.concatenate([e for e in sol.y_events if len(e)], axis=0).T
        ys = np.concatenate([ys, extra], axis=1)
    order = np.argsort(t)
    t, ys = t[order], ys[:, order]
    u, a = ys[0], ys[1]

    if len(sol.t_events[1]):
        status = "negative"
    elif len(sol.t_events[0]):
        status = "blowup"
    elif not sol.success:
        status = "error"
    else:
        status = "timeout"

    final_rate = np.nan
    if status in ("timeout", "converged"):
        try:
            fv = rhsVector((max(u[-1], 0.0), max(a[-1], 0.0)), p)
            scale = max(p.j, 1e-30)
            final_rate = float(np.max(np.abs(fv)) / scale)
            if final_rate < 1e-6:
                status = "converged"
        except ModelError:
            pass

    return Trajectory(t=t, u=u, a=a, status=status,
                      min_u=float(np.min(u)), min_a=float(np.min(a)),
                      final_u=float(u[-1]), final_a=float(a[-1]),
                      final_rate=final_rate, message=sol.message or "")


def defaultInitialConditions(scale: float = 1.0, n: int = 12,
                             seed: int = 0) -> List[Tuple[float, float]]:
    """a spread of biologically relevant initial conditions.

    includes both boundary corners (0,0), (x,0), (0,x) -- the invariant-domain
    test needs those, since that is where positivity is most fragile.
    """
    rng = np.random.default_rng(seed)
    fixed = [(0.0, 0.0), (scale, 0.0), (0.0, scale), (scale, scale),
             (1e-8, 1e-8), (10.0 * scale, 10.0 * scale)]
    rand = [(float(x), float(y)) for x, y in
            10.0 ** rng.uniform(-6.0, 1.0, size=(max(n - len(fixed), 0), 2)) * scale]
    return fixed + rand


def basinScan(p: Params, ics: Sequence[Tuple[float, float]], **kw) -> List[Dict]:
    """run many initial conditions and report the reached outcome for each."""
    rows = []
    for u0, a0 in ics:
        tr = simulate(p, u0, a0, **kw)
        row = dict(u0=float(u0), a0=float(a0))
        row.update(tr.summary())
        rows.append(row)
    return rows


def recoveryTime(p: Params, eq_u: float, eq_a: float, kick: float = 0.05,
                 t_end: float = 5.0e4) -> Optional[float]:
    """time to return within 1/e of a small displacement from an equilibrium.

    a direct dynamic readout for the critical-slowing-down prediction, computed
    from the trajectory rather than from the eigenvalue, so the two can be
    compared instead of assumed equal.
    """
    d0 = kick * max(eq_u + eq_a, 1e-12)
    tr = simulate(p, eq_u + d0, eq_a, t_end=t_end, n_out=2000)
    if tr.status not in ("converged", "timeout"):
        return None
    dist = np.hypot(tr.u - eq_u, tr.a - eq_a)
    if dist[0] <= 0.0:
        return None
    target = dist[0] / np.e
    below = np.flatnonzero(dist <= target)
    return float(tr.t[below[0]]) if len(below) else None
