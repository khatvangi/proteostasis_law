"""phase 3 follow-ups: uniqueness, dilution laws, thresholds, and a denominator.

this module closes the four items left open by `theory/FOLD_THEOREM.md`.

1. UNIQUENESS / ORDERING.  the theorem identifies constrained critical points of
   total removal on the aggregate nullcline; it does not say how many there are.
   `criticalPointsOnNullcline` enumerates the WHOLE nullcline -- both roots in
   `a` at each `u`, not just the first -- and counts sign changes of det J along
   it, so "the" fold can be checked for being unique rather than assumed.

2. THE DILUTION LAW.  `dilution.py` used one functional guess. proteome
   partitioning implies a different shape: diverting a fixed protein fraction to
   damage costs growth LINEARLY, not hyperbolically. `LINEAR` here is that second
   form. this is a second physiologically-motivated FORM, not a calibration --
   no growth-burden curve has been measured or fitted anywhere in this project.

3. THE THRESHOLD ACROSS THE BOX.  the dilution rate at which the boundary is
   destroyed was located at one parameter point. `thresholdSweep` locates it
   across draws by continuation in mu, which is far cheaper than rescanning.

4. A DENOMINATOR THAT SURVIVES DIVISION.  the old margin divided by the enzymatic
   capacity bound, which stops being an upper bound once cells divide. splitting
   critical influx into its enzymatic and dilution parts gives the exact identity

        j_crit = C_enz . phi_enz / (1 - delta)

   with  phi_enz = R_enzymatic(u*,a*) / C_enz   the fraction of enzymatic
   capacity in use at collapse, and  delta = R_dilution / j_crit  the share of
   disposal done by division. both are dimensionless and lie in [0,1). the
   divergence of j_crit as dilution rises is exactly `delta -> 1`, so the escape
   observed in `dilution.py` is quantitative rather than merely observed.

CLAIM LABELS
  Mathematical  : the j_crit identity in (4); it is algebra.
  Computational : every count, threshold and distribution.
  Empirical     : nothing. no organism data is used.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import brentq, root

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import dilution as D  # noqa: E402
from fold_theorem import KINETIC_FIELDS, phase1RunDir  # noqa: E402


# ---------------------------------------------------------------------------
# 2. a second dilution law, motivated by proteome partitioning
# ---------------------------------------------------------------------------


class LinearGrowth(D.Growth):
    """mu = mu0 . max(0, 1 - (u+a)/k_mu).

    bacterial growth laws make growth rate proportional to the ribosomal
    proteome fraction, so diverting a fixed fraction of the proteome to damage
    and its handling costs growth LINEARLY and reaches zero at a finite burden.
    the hyperbolic form in `dilution.Growth` instead approaches zero only
    asymptotically, so it never fully arrests.

    this is a second functional FORM, not a fit. no growth-burden relation has
    been measured in this project, and `k_mu` remains a free parameter.
    """

    def rate(self, u: float, a: float) -> float:
        if self.mu0 == 0.0:
            return 0.0
        if not np.isfinite(self.k_mu):
            return self.mu0
        return self.mu0 * max(0.0, 1.0 - (u + a) / self.k_mu)

    def rateGradient(self, u: float, a: float) -> Tuple[float, float]:
        if self.mu0 == 0.0 or not np.isfinite(self.k_mu):
            return 0.0, 0.0
        if (u + a) >= self.k_mu:
            return 0.0, 0.0          # clamped at zero growth; one-sided
        d = -self.mu0 / self.k_mu
        return d, d


# ---------------------------------------------------------------------------
# 1. enumerate the whole nullcline and every critical point on it
# ---------------------------------------------------------------------------


def allNullclineRootsAt(u: float, p: M.Params, g: D.Growth,
                        a_hi: float = 1e4, n: int = 400) -> List[float]:
    """every root in `a` of the aggregate equation at fixed u, not just the first.

    the curve {G = 0} is a graph over u with (generically) two branches that
    merge at a turning point, so collecting all roots at each u traces the whole
    curve rather than only the branch `lowerNullclineADil` happens to reach.
    """
    grid = np.concatenate(([0.0], np.geomspace(1e-12, a_hi, n)))
    vals = []
    for a in grid:
        try:
            vals.append(D.aggregateGDil(u, float(a), p, g))
        except M.ModelError:
            vals.append(np.nan)
    roots = []
    for i in range(len(grid) - 1):
        f0, f1 = vals[i], vals[i + 1]
        if not (np.isfinite(f0) and np.isfinite(f1)):
            continue
        if f0 == 0.0:
            roots.append(float(grid[i]))
        elif f0 * f1 < 0.0:
            try:
                roots.append(float(brentq(
                    lambda x: D.aggregateGDil(u, x, p, g),
                    grid[i], grid[i + 1], xtol=1e-15, rtol=8.9e-16)))
            except (ValueError, M.ModelError):
                pass
    return roots


def criticalPointsOnNullcline(p: M.Params, g: Optional[D.Growth] = None,
                              u_lo: float = 1e-5, u_hi: float = 50.0,
                              n_u: int = 220) -> Dict[str, object]:
    """trace the whole nullcline and count constrained critical points on it.

    returns the ordered curve, the det J values along it, and the located roots.
    a single sign change means the fold the theorem picks out is the only one.
    """
    g = g or D.Growth(0.0)
    lower, upper = [], []
    for u in np.geomspace(u_lo, u_hi, n_u):
        rs = allNullclineRootsAt(float(u), p, g)
        if not rs:
            continue
        lower.append((float(u), min(rs)))
        if len(rs) > 1:
            upper.append((float(u), max(rs)))

    # order along the curve: up the lower branch, back down the upper branch
    curve = lower + upper[::-1]
    dets = []
    for u, a in curve:
        try:
            dets.append(float(np.linalg.det(D.jacobianDil(u, a, p, g))))
        except (M.ModelError, np.linalg.LinAlgError):
            dets.append(np.nan)

    crossings = []
    for i in range(len(dets) - 1):
        d0, d1 = dets[i], dets[i + 1]
        if np.isfinite(d0) and np.isfinite(d1) and d0 * d1 < 0.0:
            crossings.append((curve[i], curve[i + 1]))

    return {"n_lower": len(lower), "n_upper": len(upper),
            "n_points": len(curve), "n_sign_changes": len(crossings),
            "crossings": crossings, "curve": curve, "dets": dets}


# ---------------------------------------------------------------------------
# 4. the decomposition that survives division
# ---------------------------------------------------------------------------


def boundaryDecomposition(p: M.Params, g: D.Growth) -> Optional[Dict[str, float]]:
    """split the critical influx into enzymatic and dilution shares.

    identity:  j_crit = C_enz . phi_enz / (1 - delta)
    """
    out = D.foldSolveDil(p, g)
    if out is None:
        return None
    j_crit, u_s, a_s = out
    f = M.fluxes(u_s, a_s, p)
    enzymatic = f["refold"] + f["degrade_u"] + f["degrade_a"]
    dil = g.rate(u_s, a_s) * (u_s + a_s)
    C_enz = M.removalCeiling(p)
    phi_enz = enzymatic / C_enz
    delta = dil / j_crit if j_crit > 0 else np.nan
    return {"j_crit": j_crit, "u_star": u_s, "a_star": a_s,
            "C_enz": C_enz, "phi_enz": phi_enz, "delta": delta,
            "identity_rhs": C_enz * phi_enz / (1.0 - delta),
            "identity_rel_err": abs(C_enz * phi_enz / (1.0 - delta) - j_crit)
                                / max(j_crit, 1e-300)}


# ---------------------------------------------------------------------------
# 3. where the boundary is destroyed, across the parameter box
# ---------------------------------------------------------------------------


def _solveSeeded(p: M.Params, g: D.Growth, seed: Tuple[float, float]):
    def residual(x):
        u, a = float(np.exp(x[0])), float(np.exp(x[1]))
        try:
            return [D.aggregateGDil(u, a, p, g),
                    float(np.linalg.det(D.jacobianDil(u, a, p, g)))]
        except (M.ModelError, np.linalg.LinAlgError):
            return [1e6, 1e6]

    s = root(residual, [np.log(seed[0]), np.log(max(seed[1], 1e-12))],
             method="hybr", options={"xtol": 1e-13})
    if not s.success:
        return None
    u, a = float(np.exp(s.x[0])), float(np.exp(s.x[1]))
    if not (np.isfinite(u) and np.isfinite(a) and u > 0.0):
        return None
    if abs(D.aggregateGDil(u, a, p, g)) > 1e-9:
        return None
    if abs(float(np.linalg.det(D.jacobianDil(u, a, p, g)))) > 1e-7:
        return None
    return u, a


def dilutionThreshold(p: M.Params, mu_hi: float = 4.0, n_step: int = 26,
                      n_bisect: int = 22) -> Optional[Dict[str, float]]:
    """largest CONSTANT dilution rate at which a boundary still exists.

    continuation in mu, seeded from the previous solution, then bisection on the
    last bracket. returns None if no boundary exists even at mu = 0.
    """
    base = D.foldSolveDil(p, D.Growth(0.0))
    if base is None:
        return None
    seed = (base[1], base[2])
    lo, hi = 0.0, None
    last = None
    for mu in np.geomspace(1e-4, mu_hi, n_step):
        got = _solveSeeded(p, D.Growth(mu0=float(mu)), seed)
        if got is None:
            hi = float(mu)
            break
        lo, seed = float(mu), got
        last = got
    if hi is None:
        return {"threshold": float("inf"), "bracketed": False,
                "j_at_lo": None, "delta_at_lo": None}

    seed_lo = seed
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        got = _solveSeeded(p, D.Growth(mu0=mid), seed_lo)
        if got is None:
            hi = mid
        else:
            lo, seed_lo = mid, got

    g_lo = D.Growth(mu0=lo)
    dec = boundaryDecomposition(p, g_lo)
    return {"threshold": lo, "bracketed": True,
            "j_at_lo": None if dec is None else dec["j_crit"],
            "delta_at_lo": None if dec is None else dec["delta"],
            "a_star_at_lo": None if last is None else seed_lo[1]}


def thresholdSweep(k: int = 30, seed: int = 5) -> pd.DataFrame:
    """dilution thresholds across draws from the phase 1 parameter box."""
    run = phase1RunDir()
    c = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    draws = c.sample(n=k, random_state=seed)

    rows = []
    for idx, (_, r) in enumerate(draws.iterrows()):
        kin = {f: float(r["p_" + f]) for f in KINETIC_FIELDS}
        chi = float(r["p_chi"])
        try:
            p = M.Params(**kin).with_(nu=float(r["p_nu"]), c_tot=chi,
                                      d_tot=1.0 - chi).validate()
        except M.ModelError:
            continue
        try:
            t = dilutionThreshold(p)
        except (M.ModelError, ValueError, np.linalg.LinAlgError):
            continue
        if t is None:
            continue
        rows.append({"draw": idx, "threshold": t["threshold"],
                     "bracketed": t["bracketed"],
                     "delta_at_threshold": t["delta_at_lo"],
                     "j_at_threshold": t["j_at_lo"],
                     "ceiling": M.removalCeiling(p)})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 1b. what the multiple critical points MEAN: dilution makes the system bistable
# ---------------------------------------------------------------------------


def settle(j: float, x0: Tuple[float, float], p: M.Params, g: D.Growth,
           t_end: float = 5e4) -> Tuple[float, float]:
    """integrate to steady state from a given initial state."""
    from scipy.integrate import solve_ivp
    pj = p.with_(j=float(j)).validate()

    def field(_t, x):
        return list(D.rhsDil(max(x[0], 0.0), max(x[1], 0.0), pj, g))

    s = solve_ivp(field, [0.0, t_end], list(x0), method="Radau",
                  rtol=1e-10, atol=1e-13)
    return float(s.y[0, -1]), float(s.y[1, -1])


def hysteresisSweep(p: M.Params, g: D.Growth, js) -> Dict[str, object]:
    """sweep j up from zero burden, then back down from the attained state.

    when dilution admits a BOUNDED high-burden attractor, loss of the low-burden
    branch is a jump to that attractor rather than divergence, and the up and
    down sweeps disagree over an interval. that interval is the bistable window
    and should be bracketed by the two constrained critical points.
    """
    js = list(js)
    x, up = (0.0, 0.0), {}
    for j in js:
        x = settle(j, x, p, g)
        up[j] = x
    x, down = up[js[-1]], {}
    for j in reversed(js):
        x = settle(j, x, p, g)
        down[j] = x
    hyst = [j for j in js
            if abs(up[j][1] - down[j][1]) / max(down[j][1], 1e-12) > 0.01]
    return {"up": up, "down": down, "hysteretic_j": hyst,
            "window": (min(hyst), max(hyst)) if hyst else None}


def main() -> int:
    p = M.Params().validate()

    print("1. UNIQUENESS -- critical points on the whole nullcline")
    for label, g in (("no dilution", D.Growth(0.0)),
                     ("constant mu=0.04", D.Growth(0.04)),
                     ("hyperbolic mu0=0.3", D.Growth(0.3, 0.5)),
                     ("linear mu0=0.3", LinearGrowth(0.3, 2.0))):
        r = criticalPointsOnNullcline(p, g)
        print("   %-20s points=%-5d lower=%-4d upper=%-4d det J sign changes = %d"
              % (label, r["n_points"], r["n_lower"], r["n_upper"],
                 r["n_sign_changes"]))
    print()

    print("2. DILUTION LAW -- hyperbolic vs linear (proteome-partitioning shape)")
    print("   %-10s %-16s %-16s %-10s" % ("mu0", "hyperbolic j_crit", "linear j_crit", "both?"))
    for mu0 in (0.0, 0.05, 0.1, 0.3, 0.6, 1.0):
        h = D.foldSolveDil(p, D.Growth(mu0, 0.5))
        l_ = D.foldSolveDil(p, LinearGrowth(mu0, 2.0))
        print("   %-10.3g %-16s %-16s %-10s" % (
            mu0,
            "-" if h is None else "%.6f" % h[0],
            "-" if l_ is None else "%.6f" % l_[0],
            "yes" if (h is not None and l_ is not None) else
            ("neither" if (h is None and l_ is None) else "DIFFER")))
    print()

    print("3./4. DECOMPOSITION  j_crit = C_enz . phi_enz / (1 - delta)")
    print("   %-10s %-11s %-11s %-11s %-11s" % ("mu", "j_crit", "phi_enz", "delta", "identity err"))
    for mu in (0.0, 0.02, 0.04, 0.06, 0.08):
        d = boundaryDecomposition(p, D.Growth(mu0=mu))
        if d is None:
            print("   %-10.3g  none" % mu)
            continue
        print("   %-10.3g %-11.6f %-11.6f %-11.6f %-11.2e"
              % (mu, d["j_crit"], d["phi_enz"], d["delta"], d["identity_rel_err"]))
    print()

    print("1b. WHAT THE EXTRA CRITICAL POINT MEANS -- bistability under dilution")
    g = D.Growth(0.04)
    pts = []
    for seed in ((0.45, 0.35), (0.16, 1.95)):
        got = _solveSeeded(p, g, seed)
        if got is not None:
            pts.append((D.removalTotal(got[0], got[1], p, g), got[0], got[1]))
    pts.sort()
    for j, u, a in pts:
        print("   critical point: j=%.6f  u=%.6f  a=%.6f" % (j, u, a))
    h = hysteresisSweep(p, g, [0.10, 0.14, 0.155, 0.158, 0.16, 0.17, 0.18,
                               0.19, 0.194, 0.196, 0.21])
    print("   bistable window from the sweeps : %s" % (h["window"],))
    if len(pts) == 2 and h["window"]:
        print("   window lies inside [j_lower, j_upper] = [%.6f, %.6f] : %s"
              % (pts[0][0], pts[1][0],
                 pts[0][0] <= h["window"][0] and h["window"][1] <= pts[1][0]))
    print("   loss of the low branch is a JUMP to a bounded high-aggregate state,")
    print("   not divergence: at j=0.21 the attained state is u=%.4f a=%.4f"
          % h["up"][0.21])
    print()

    print("3. THRESHOLD ACROSS THE PARAMETER BOX")
    S = thresholdSweep()
    if len(S) == 0:
        print("   SKIP: phase 1 run root absent")
        return 0
    fin = S[np.isfinite(S["threshold"])]
    print("   draws with a boundary at mu = 0        : %d" % len(S))
    print("   of which the boundary is destroyed     : %d" % int(S["bracketed"].sum()))
    print("   threshold mu, quantiles p10/p50/p90    : %.4g / %.4g / %.4g"
          % tuple(fin["threshold"].quantile([.1, .5, .9])))
    print("   threshold spans                        : %.2f decades"
          % np.log10(fin["threshold"].max() / fin["threshold"].min()))
    d = S["delta_at_threshold"].dropna()
    print("   dilution share delta AT the threshold  : median %.4f  (p10 %.4f, p90 %.4f)"
          % (d.median(), d.quantile(.1), d.quantile(.9)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
