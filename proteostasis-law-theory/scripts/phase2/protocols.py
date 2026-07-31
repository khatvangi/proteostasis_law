"""the two root-search + attraction protocols, blind to which model they solve.

the whole point of the 2x2 factorial is to separate the MODEL-FORM effect from
the SOLVER effect.  that only works if the same protocol code runs against both
model forms, so both protocols here take a `models.ModelAdapter` and touch the
model only through its interface.

protocol P_BORON  -- transcribed from `scripts/proteostasis/equilibria.py`
                     (`findEquilibria`, 9x9 log grid over [1e-7, 1e4], hybr in
                     log coordinates, RELATIVE residual tolerance) plus the
                     attractor confirmation of `run_experiment_c._confirmsAttractor`
                     (+/-10% kick, 5% relative return, Radau with analytic
                     jacobian and blow-up / negativity events).

protocol P_NITROGEN -- transcribed from
                     nitrogen:.../src/lhs_sweep.py::classify_one (four fixed
                     linear guesses, ABSOLUTE residual tolerance 1e-7,
                     +/-0.1% kick, ABSOLUTE 1e-3 return, DOP853 to t = 150).

`tests/phase2/test_protocol_transcription.py` asserts that P_BORON's root stage
reproduces `equilibria.findEquilibria` exactly on the boron model, and that
P_NITROGEN reproduces nitrogen's `classify_one` labels on a reference dump taken
from the real nitrogen module.  neither transcription is trusted on inspection.

labels
------
`stable_attractor`, `no_bounded_attractor_operational`, `numerical_failure`.

deliberately NOT `stable_subthreshold` / `stable_overthreshold`: the threshold
split depends on which admissibility criterion is applied, and there are four
defensible ones here (burden vs componentwise) x (total vs free coordinates).
mixing them silently is exactly the confound this benchmark exists to remove, so
the protocols report the attractor and its coordinates, and the criteria are
applied downstream and reported side by side.

`no_bounded_attractor_operational` is an OPERATIONAL label.  it means no stable
attracting state was found within this finite numerical protocol.  it is not a
mathematical nonexistence result.
"""

from __future__ import annotations

from typing import Dict, Optional, Sequence, Tuple

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root

from .models import ModelAdapter

LABEL_STABLE = "stable_attractor"
LABEL_NONE = "no_bounded_attractor_operational"
LABEL_FAIL = "numerical_failure"

# --- boron protocol constants (equilibria.py:29-30, 122-123; simulate.py:25) --
_LOG_FLOOR = -30.0
_LOG_CEIL = 12.0
BORON_GRID = 9
BORON_LO, BORON_HI = 1e-7, 1e4
BORON_RES_TOL = 1e-8            # RELATIVE, divided by adapter.residualScale()
BORON_DEDUPE_RTOL = 1e-5
BORON_KICK = 0.1                # configs/phase1/experiment_c.json:165
BORON_REL_TOL = 0.05            # configs/phase1/experiment_c.json:166
BORON_T_END = 5.0e4             # configs/phase1/experiment_c.json:167
BORON_BLOWUP = 1.0e6
BORON_NEG_TOL = -1e-9

# --- nitrogen protocol constants (lhs_sweep.py:51-69) ------------------------
NITROGEN_GUESSES = ((0.01, 0.01), (0.1, 0.05), (0.5, 0.2), (1.0, 1.0))
NITROGEN_RES_TOL = 1e-7         # ABSOLUTE
NITROGEN_DEDUPE = 1e-5
NITROGEN_EIG_TOL = 1e-8
NITROGEN_KICK = 1e-3
NITROGEN_RETURN_TOL = 1e-3      # ABSOLUTE
NITROGEN_T_END = 150.0


# ---------------------------------------------------------------------------
# shared helpers
# ---------------------------------------------------------------------------


def _eigen(adapter: ModelAdapter, x) -> np.ndarray:
    return np.linalg.eigvals(adapter.jacobian(x))


def _record(adapter: ModelAdapter, x, residual: float, eig: np.ndarray) -> Dict:
    uf, af = adapter.freeCoordinates(x)
    return {
        "u": float(x[0]), "a": float(x[1]),
        "u_free": uf, "a_free": af,
        "residual": float(residual),
        "eig_real_max": float(np.max(eig.real)),
        "eig_real_min": float(np.min(eig.real)),
        "eig_imag_absmax": float(np.max(np.abs(eig.imag))),
    }


# ---------------------------------------------------------------------------
# protocol P_BORON
# ---------------------------------------------------------------------------


def _clipLog(x: np.ndarray) -> np.ndarray:
    return np.clip(np.nan_to_num(x, nan=0.0, posinf=_LOG_CEIL, neginf=_LOG_FLOOR),
                   _LOG_FLOOR, _LOG_CEIL)


def _logResidual(x, adapter: ModelAdapter) -> np.ndarray:
    xs = _clipLog(np.asarray(x, dtype=float))
    return adapter.rhsVector((float(np.exp(xs[0])), float(np.exp(xs[1]))))


def _logJacobian(x, adapter: ModelAdapter) -> np.ndarray:
    xs = _clipLog(np.asarray(x, dtype=float))
    u, a = float(np.exp(xs[0])), float(np.exp(xs[1]))
    return adapter.jacobian((u, a)) * np.array([[u, a], [u, a]])


def _classifyState(adapter: ModelAdapter, u: float, a: float,
                   res_tol: float = BORON_RES_TOL) -> Optional[Dict]:
    """relative-residual acceptance, matching `equilibria.classifyState`."""
    try:
        f = adapter.rhsVector((u, a))
        eig = _eigen(adapter, (u, a))
    except (ValueError, FloatingPointError, np.linalg.LinAlgError):
        return None
    res = float(np.max(np.abs(f))) / adapter.residualScale()
    if not np.isfinite(res) or res > res_tol:
        return None
    if not np.all(np.isfinite(eig)):
        return None
    rec = _record(adapter, (u, a), res, eig)
    rec["stable"] = bool(np.max(eig.real) < 0.0)
    return rec


def _solveFrom(adapter: ModelAdapter, u0: float, a0: float,
               res_tol: float = BORON_RES_TOL) -> Optional[Dict]:
    x0 = np.clip(np.array([np.log(max(u0, 1e-13)), np.log(max(a0, 1e-13))]),
                 _LOG_FLOOR, _LOG_CEIL)
    try:
        sol = root(_logResidual, x0, args=(adapter,), jac=_logJacobian,
                   method="hybr", options=dict(xtol=1e-13))
    except (ValueError, FloatingPointError, np.linalg.LinAlgError):
        return None
    if not np.all(np.isfinite(sol.x)):
        return None
    xs = _clipLog(sol.x)
    # a solution pinned to the clipping box is an artefact of the box, not a root
    if np.any(np.abs(xs - np.asarray(sol.x, dtype=float)) > 1e-9):
        return None
    return _classifyState(adapter, float(np.exp(xs[0])), float(np.exp(xs[1])), res_tol)


def findEquilibriaBoron(adapter: ModelAdapter, n_grid: int = BORON_GRID,
                        lo: float = BORON_LO, hi: float = BORON_HI,
                        res_tol: float = BORON_RES_TOL,
                        dedupe_rtol: float = BORON_DEDUPE_RTOL) -> list:
    """blind multi-start over a logarithmic grid; sorted by total burden."""
    grid = np.logspace(np.log10(lo), np.log10(hi), n_grid)
    found: list = []
    for u0 in grid:
        for a0 in grid:
            eq = _solveFrom(adapter, float(u0), float(a0), res_tol)
            if eq is None:
                continue
            dup = False
            for e in found:
                du = abs(e["u"] - eq["u"]) / max(abs(e["u"]), abs(eq["u"]), 1e-30)
                da = abs(e["a"] - eq["a"]) / max(abs(e["a"]), abs(eq["a"]), 1e-30)
                if du < dedupe_rtol and da < dedupe_rtol:
                    dup = True
                    break
            if not dup:
                found.append(eq)
    found.sort(key=lambda e: e["u"] + e["a"])
    return found


def _simulateBoron(adapter: ModelAdapter, u0: float, a0: float,
                   t_end: float = BORON_T_END, blowup: float = BORON_BLOWUP):
    """Radau with analytic jacobian and the two terminal events of simulate.py."""
    def f(_t, y):
        return adapter.rhsVector((max(y[0], 0.0), max(y[1], 0.0)))

    def jac(_t, y):
        return adapter.jacobian((max(y[0], 0.0), max(y[1], 0.0)))

    def evBlow(_t, y):
        return blowup - (y[0] + y[1])
    evBlow.terminal = True
    evBlow.direction = -1.0

    def evNeg(_t, y):
        return min(y[0], y[1]) - BORON_NEG_TOL
    evNeg.terminal = True
    evNeg.direction = -1.0

    sol = solve_ivp(f, (0.0, t_end), [float(u0), float(a0)], method="Radau",
                    jac=jac, events=(evBlow, evNeg), rtol=1e-9, atol=1e-12)
    if len(sol.t_events[1]):
        status = "negative"
    elif len(sol.t_events[0]):
        status = "blowup"
    elif not sol.success:
        status = "error"
    else:
        status = "timeout"
    return sol, status


def _confirmsAttractorBoron(adapter: ModelAdapter, eq: Dict,
                            t_end: float = BORON_T_END) -> Tuple[bool, str]:
    """displace both ways and require BOTH trajectories to come back.

    transcribed from `run_experiment_c._confirmsAttractor`: one-sided return is
    not enough, because a trajectory that falls off a fold onto the low branch
    also 'returns' in one direction.
    """
    for s in (1.0 + BORON_KICK, 1.0 - BORON_KICK):
        try:
            sol, status = _simulateBoron(adapter, eq["u"] * s, eq["a"] * s, t_end)
        except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
            return False, f"integrator error: {type(exc).__name__}"
        if status not in ("converged", "timeout"):
            return False, status
        fu, fa = float(sol.y[0, -1]), float(sol.y[1, -1])
        back = (abs(fu - eq["u"]) <= BORON_REL_TOL * max(eq["u"], 1e-12)
                and abs(fa - eq["a"]) <= BORON_REL_TOL * max(eq["a"], 1e-12))
        if not back:
            return False, "did not return"
    return True, "confirmed"


def classifyBoron(adapter: ModelAdapter, t_end: float = BORON_T_END) -> Dict:
    """P_BORON: dense log multistart, then dynamic confirmation of the lowest
    stable equilibrium."""
    try:
        roots = findEquilibriaBoron(adapter)
    except Exception as exc:                                    # noqa: BLE001
        return {"label": LABEL_FAIL, "status": f"{type(exc).__name__}: {exc}",
                "root": None, "n_roots": 0, "n_stable": 0}
    stable = [r for r in roots if r["stable"] and r["u"] >= 0.0 and r["a"] >= 0.0]
    if not stable:
        return {"label": LABEL_NONE, "status": "no stable root",
                "root": None, "n_roots": len(roots), "n_stable": 0}
    eq = stable[0]                                              # lowest burden
    ok, status = _confirmsAttractorBoron(adapter, eq, t_end)
    return {"label": LABEL_STABLE if ok else LABEL_NONE, "status": status,
            "root": eq, "n_roots": len(roots), "n_stable": len(stable)}


# ---------------------------------------------------------------------------
# protocol P_NITROGEN
# ---------------------------------------------------------------------------


def _solveEquilibriumNitrogen(adapter: ModelAdapter, guess) -> Dict:
    sol = root(lambda z: adapter.rhsVector(z), np.asarray(guess, dtype=float),
               jac=lambda z: adapter.jacobian(z), tol=1e-10)
    residual = float(np.linalg.norm(adapter.rhsVector(sol.x), ord=np.inf))
    return {"x": np.asarray(sol.x, dtype=float), "residual": residual,
            "success": bool(sol.success and np.isfinite(residual)),
            "message": str(sol.message)}


def classifyNitrogen(adapter: ModelAdapter,
                     guesses: Sequence = NITROGEN_GUESSES) -> Dict:
    """P_NITROGEN: four fixed linear guesses, then the +/-0.1% return check.

    the return check is retained exactly as nitrogen wrote it, including the
    known near-vacuity for sub-threshold roots (an ABSOLUTE 1e-3 tolerance
    against a `1e-3 * ||x*||` displacement).  weakening or strengthening it here
    would change the solver factor and destroy the factorial.
    """
    roots: list = []
    try:
        for guess in guesses:
            r = _solveEquilibriumNitrogen(adapter, guess)
            if not (r["success"] and np.all(np.isfinite(r["x"]))
                    and np.all(r["x"] >= -1e-8) and r["residual"] < NITROGEN_RES_TOL):
                continue
            if any(np.linalg.norm(r["x"] - old["x"]) < NITROGEN_DEDUPE for old in roots):
                continue
            eig = _eigen(adapter, r["x"])
            if np.all(np.real(eig) < -NITROGEN_EIG_TOL):
                r["stable"] = "stable"
            elif np.any(np.real(eig) > NITROGEN_EIG_TOL):
                r["stable"] = "unstable"
            else:
                r["stable"] = "nonhyperbolic"
            r["eig"] = eig
            roots.append(r)
    except Exception as exc:                                    # noqa: BLE001
        return {"label": LABEL_FAIL, "status": f"{type(exc).__name__}: {exc}",
                "root": None, "n_roots": 0, "n_stable": 0}

    stable = [r for r in roots if r["stable"] == "stable" and np.all(r["x"] >= 0)]
    if not stable:
        return {"label": LABEL_NONE, "status": "no stable root",
                "root": None, "n_roots": len(roots), "n_stable": 0}

    r = min(stable, key=lambda q: np.linalg.norm(q["x"]))
    xstar = r["x"]
    rec = _record(adapter, xstar, r["residual"], r["eig"])
    perturbed = np.maximum(xstar * (1.0 + np.array([NITROGEN_KICK, -NITROGEN_KICK])), 1e-10)
    try:
        sol = solve_ivp(lambda t, x: adapter.rhsVector(x), (0.0, NITROGEN_T_END),
                        perturbed, method="DOP853", rtol=2e-8, atol=1e-10)
    except Exception as exc:                                    # noqa: BLE001
        return {"label": LABEL_FAIL, "status": f"{type(exc).__name__}: {exc}",
                "root": rec, "n_roots": len(roots), "n_stable": len(stable)}
    endpoint = sol.y[:, -1]
    bad = (not sol.success or not np.all(np.isfinite(endpoint))
           or np.any(endpoint < -1e-7) or np.any(endpoint > 1e6)
           or np.linalg.norm(endpoint - xstar) > NITROGEN_RETURN_TOL)
    rec["u_end"], rec["a_end"] = float(endpoint[0]), float(endpoint[1])
    return {"label": LABEL_NONE if bad else LABEL_STABLE,
            "status": str(sol.message), "root": rec,
            "n_roots": len(roots), "n_stable": len(stable)}


PROTOCOLS = {"boron": classifyBoron, "nitrogen": classifyNitrogen}
