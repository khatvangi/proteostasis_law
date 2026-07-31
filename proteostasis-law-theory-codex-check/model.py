"""Independent nondimensional two-state finite-resource U/A model."""
from __future__ import annotations

from dataclasses import dataclass, asdict
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares


@dataclass(frozen=True)
class Params:
    J: float = 0.35
    N: float = 0.5
    Ctot: float = 1.4
    Dtot: float = 1.1
    KN: float = 1.0
    KCU: float = 0.8
    KCA: float = 1.2
    KDU: float = 0.9
    KDA: float = 1.1
    kref: float = 1.5
    Kref: float = 0.7
    kU: float = 1.0
    KU: float = 0.8
    kn: float = 0.04
    m: float = 2.0
    kg: float = 0.03
    kdis: float = 0.65
    Kdis: float = 0.9
    kA: float = 0.75
    KA: float = 1.0

    def dict(self):
        return asdict(self)


def resources(y, p: Params):
    U, A = np.asarray(y, dtype=float)
    cden = 1.0 + p.N / p.KN + U / p.KCU + A / p.KCA
    dden = 1.0 + U / p.KDU + A / p.KDA
    return p.Ctot / cden, p.Dtot / dden


def rhs(_t, y, p: Params):
    U, A = np.asarray(y, dtype=float)
    Cf, Df = resources(y, p)
    ref = p.kref * Cf * U / (p.Kref + U)
    deg_u = p.kU * Df * U / (p.KU + U)
    nuc = p.kn * U**p.m
    growth = p.kg * U * A
    dis = p.kdis * Cf * A / (p.Kdis + A)
    deg_a = p.kA * Df * A / (p.KA + A)
    return np.array([p.J-ref-deg_u-nuc-growth+dis,
                     nuc+growth-dis-deg_a])


def jacobian(y, p: Params):
    U, A = np.asarray(y, dtype=float)
    cden = 1.0 + p.N/p.KN + U/p.KCU + A/p.KCA
    dden = 1.0 + U/p.KDU + A/p.KDA
    Cf, Df = p.Ctot/cden, p.Dtot/dden
    dC = np.array([-p.Ctot/(cden*cden*p.KCU), -p.Ctot/(cden*cden*p.KCA)])
    dD = np.array([-p.Dtot/(dden*dden*p.KDU), -p.Dtot/(dden*dden*p.KDA)])
    sr, su, sd, sa = U/(p.Kref+U), U/(p.KU+U), A/(p.Kdis+A), A/(p.KA+A)
    dr = p.kref*(dC*sr + Cf*np.array([p.Kref/(p.Kref+U)**2, 0.0]))
    du = p.kU*(dD*su + Df*np.array([p.KU/(p.KU+U)**2, 0.0]))
    dn = np.array([p.kn*p.m*U**(p.m-1), 0.0])
    dg = np.array([p.kg*A, p.kg*U])
    dd = p.kdis*(dC*sd + Cf*np.array([0.0, p.Kdis/(p.Kdis+A)**2]))
    da = p.kA*(dD*sa + Df*np.array([0.0, p.KA/(p.KA+A)**2]))
    return np.vstack((-dr-du-dn-dg+dd, dn+dg-dd-da))


def finite_difference_jacobian(y, p: Params, rel_step=1e-6):
    y = np.asarray(y, float)
    out = np.empty((2, 2))
    for j in range(2):
        h = rel_step * max(1.0, abs(y[j]))
        yp, ym = y.copy(), y.copy()
        yp[j] += h
        ym[j] -= h
        out[:, j] = (rhs(0, yp, p)-rhs(0, ym, p))/(2*h)
    return out


def equilibria(p: Params, search_max=100.0, grid=9, residual_tol=1e-8):
    """Find positive equilibria via deterministic log-coordinate multistart."""
    roots = []
    starts = np.linspace(-8.0, np.log(search_max), grid)
    for lu in starts:
        for la in starts:
            fit = least_squares(lambda z: rhs(0, np.exp(z), p), [lu, la],
                                bounds=(-20.0, np.log(search_max)),
                                xtol=1e-11, ftol=1e-11, gtol=1e-11,
                                max_nfev=1200)
            y = np.exp(fit.x)
            if np.linalg.norm(rhs(0, y, p), ord=np.inf) <= residual_tol:
                if not any(np.linalg.norm(np.log(y/r)) < 1e-5 for r in roots):
                    roots.append(y)
    return sorted(roots, key=lambda z: (z[0]+z[1], z[0]))


def classify_equilibrium(y, p: Params, eig_tol=1e-7):
    eig = np.linalg.eigvals(jacobian(y, p))
    mx = float(np.max(eig.real))
    label = "stable" if mx < -eig_tol else ("unstable" if mx > eig_tol else "nonhyperbolic")
    return label, eig


def integrate(p: Params, y0, t_final=500.0, escape=200.0):
    def escape_event(_t, y):
        return escape - np.max(y)
    escape_event.terminal = True
    escape_event.direction = -1
    sol = solve_ivp(lambda t, y: rhs(t, y, p), (0, t_final), y0, method="LSODA",
                    rtol=2e-8, atol=1e-10, events=escape_event)
    ok = bool(sol.success and np.all(np.isfinite(sol.y)) and np.min(sol.y) >= -1e-8)
    escaped = bool(sol.t_events and len(sol.t_events[0]))
    return sol, ok, escaped

