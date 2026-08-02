"""phase 3: growth dilution, and what it does to the fold.

WHY THIS EXISTS
---------------
the two-state model in `scripts/proteostasis/model.py` has no cell division. in a
growing cell every species is diluted by volume expansion, and for most E. coli
proteins dilution outpaces proteolysis. a "disposal capacity" that omits it is
therefore the wrong object, and `theory/FOLD_THEOREM.md` records this as the most
serious limit on the phase 3 result.

this module adds dilution WITHOUT touching the frozen upstream model: the
rapid-equilibrium binding algebra is unaffected by dilution, so the diluted field
is just

        du/dt = (model du/dt) - mu(u,a).u
        da/dt = (model da/dt) - mu(u,a).a

THE THEOREM SURVIVES
--------------------
both structural facts still hold:

  (i)  `j` still enters du/dt and nowhere else, so the aggregate nullcline is
       still a fixed curve in burden space;
  (ii) the internal u <-> a transfer still cancels, giving

           du/dt + da/dt = j - R_tot ,     R_tot = R + mu(u,a).(u + a)

so the same row operation gives `det J = -( grad R_tot x grad G_dil )` and the
fold is still the constrained maximum of TOTAL removal -- now including dilution
-- on the aggregate nullcline. this holds for ANY dilution law mu(u,a); nothing
below depends on the particular form.

WHAT CHANGES
------------
the removal ceiling does not survive. `R_tot` contains `mu.(u+a)`, which is
unbounded in burden, so the analytic bound A8 (`j > c_tot + (rho_U+rho_A).d_tot`
admits no bounded state) is FALSE once cells divide. whether a fold survives at
all is then a question about the dilution law, and `foldSurvival` below answers
it numerically rather than by assumption.

CLAIM LABELS
------------
  Mathematical  : the identity, and the mu -> 0 reduction to the frozen model.
  Computational : every fold-survival number; a property of this model and these
                  parameter ranges.
  Empirical     : nothing. no organism data is used.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
from scipy.optimize import brentq, root

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "scripts"))

from proteostasis import model as M  # noqa: E402


@dataclass(frozen=True)
class Growth:
    """dilution law.

    `mu0`   dilution rate at zero burden, in the model's time units (1/k_ref).
    `k_mu`  burden at which the growth rate has halved. `inf` means CONSTANT
            dilution, i.e. no feedback from burden onto growth -- the case that
            turns out to abolish the fold.
    """

    mu0: float = 0.0
    k_mu: float = float("inf")

    def rate(self, u: float, a: float) -> float:
        if self.mu0 == 0.0:
            return 0.0
        if not np.isfinite(self.k_mu):
            return self.mu0
        return self.mu0 / (1.0 + (u + a) / self.k_mu)

    def rateGradient(self, u: float, a: float) -> Tuple[float, float]:
        """d(mu)/du and d(mu)/da. equal, since mu depends on burden only."""
        if self.mu0 == 0.0 or not np.isfinite(self.k_mu):
            return 0.0, 0.0
        m = self.rate(u, a)
        d = -(m * m) / (self.mu0 * self.k_mu)
        return d, d


# ---------------------------------------------------------------------------
# diluted vector field
# ---------------------------------------------------------------------------


def rhsDil(u: float, a: float, p: M.Params, g: Growth) -> Tuple[float, float]:
    du, da = M.rhs(u, a, p)
    mu = g.rate(u, a)
    return du - mu * u, da - mu * a


def aggregateGDil(u: float, a: float, p: M.Params, g: Growth) -> float:
    return rhsDil(u, a, p, g)[1]


def removalTotal(u: float, a: float, p: M.Params, g: Growth) -> float:
    """R_tot = enzymatic removal + dilution.

    this is the quantity that equals j at every equilibrium once cells divide.
    """
    f = M.fluxes(u, a, p)
    enzymatic = f["refold"] + f["degrade_u"] + f["degrade_a"]
    return enzymatic + g.rate(u, a) * (u + a)


def jacobianDil(u: float, a: float, p: M.Params, g: Growth) -> np.ndarray:
    """analytic jacobian of the diluted field.

    the dilution vector field is (-mu.u, -mu.a), whose jacobian is

        -[[ mu + u.mu_u ,      u.mu_a      ],
          [     a.mu_u  ,  mu + a.mu_a     ]]

    for CONSTANT mu this collapses to -mu.I, so the diluted jacobian is
    `J - mu.I` and the saddle-node condition det(J - mu.I) = 0 says exactly that
    mu is an eigenvalue of the undiluted jacobian.
    """
    J = M.jacobian(u, a, p)
    mu = g.rate(u, a)
    mu_u, mu_a = g.rateGradient(u, a)
    D = -np.array([[mu + u * mu_u, u * mu_a],
                   [a * mu_u, mu + a * mu_a]], dtype=float)
    return J + D


def numericalJacobianDil(u: float, a: float, p: M.Params, g: Growth,
                         h_rel: float = 1e-6) -> np.ndarray:
    out = np.zeros((2, 2))
    for k, (du_, da_) in enumerate(((h_rel * u, 0.0), (0.0, h_rel * a))):
        plus = rhsDil(u + du_, a + da_, p, g)
        minus = rhsDil(u - du_, a - da_, p, g)
        step = 2.0 * (du_ + da_)
        out[:, k] = [(plus[0] - minus[0]) / step, (plus[1] - minus[1]) / step]
    return out


def determinantIdentityDil(u: float, a: float, p: M.Params, g: Growth,
                           h_rel: float = 1e-6) -> Dict[str, float]:
    """check det J_dil == -(grad R_tot x grad G_dil)."""
    detJ = float(np.linalg.det(jacobianDil(u, a, p, g)))

    def grad(fn):
        hu, ha = h_rel * u, h_rel * a
        return ((fn(u + hu, a) - fn(u - hu, a)) / (2 * hu),
                (fn(u, a + ha) - fn(u, a - ha)) / (2 * ha))

    Ru, Ra = grad(lambda uu, aa: removalTotal(uu, aa, p, g))
    Gu, Ga = grad(lambda uu, aa: aggregateGDil(uu, aa, p, g))
    cross = Ru * Ga - Ra * Gu
    scale = max(abs(detJ), abs(cross), 1e-300)
    return {"det_J": detJ, "minus_cross": -cross,
            "rel_err": abs(detJ + cross) / scale}


# ---------------------------------------------------------------------------
# folds under dilution
# ---------------------------------------------------------------------------


def lowerNullclineADil(u: float, p: M.Params, g: Growth,
                       a_hi: float = 1e6) -> Optional[float]:
    """first root in `a` of the diluted aggregate equation at fixed u."""
    if aggregateGDil(u, 0.0, p, g) <= 0.0:
        return 0.0
    a_prev, a = 0.0, 1e-10
    while a < a_hi:
        try:
            if aggregateGDil(u, a, p, g) <= 0.0:
                return float(brentq(lambda x: aggregateGDil(u, x, p, g),
                                    a_prev, a, xtol=1e-15, rtol=8.9e-16))
        except (M.ModelError, ValueError):
            return None
        a_prev, a = a, a * 1.5
    return None


def foldSolveDil(p: M.Params, g: Growth, u_lo: float = 1e-7, u_hi: float = 200.0,
                 n_scan: int = 160) -> Optional[Tuple[float, float, float]]:
    """solve {G_dil = 0, det J_dil = 0}; return (j_crit, u*, a*) or None."""
    us = np.geomspace(u_lo, u_hi, n_scan)
    prev, guess = None, None
    for u in us:
        a = lowerNullclineADil(u, p, g)
        if a is None:
            break
        try:
            d = float(np.linalg.det(jacobianDil(u, a, p, g)))
        except (M.ModelError, np.linalg.LinAlgError):
            break
        if prev is not None and np.sign(d) != np.sign(prev[1]):
            guess = (float(np.sqrt(prev[0] * u)), a)
            break
        prev = (u, d)
    if guess is None:
        # det J need not change sign along the scan: on this model the lower
        # nullcline branch TERMINATES (experiment B records reason "branch
        # terminates") rather than crossing, so the last valid point is the
        # bracket. this mirrors `fold_theorem.foldSolve`; dropping it reports
        # "no fold" even at mu = 0, where a fold demonstrably exists.
        if prev is None:
            return None
        a0 = lowerNullclineADil(prev[0], p, g)
        guess = (prev[0], a0 if a0 else 1e-6)

    def residual(x):
        uu, aa = float(np.exp(x[0])), float(np.exp(x[1]))
        try:
            return [aggregateGDil(uu, aa, p, g),
                    float(np.linalg.det(jacobianDil(uu, aa, p, g)))]
        except (M.ModelError, np.linalg.LinAlgError):
            return [1e6, 1e6]

    sol = root(residual, [np.log(guess[0]), np.log(max(guess[1], 1e-12))],
               method="hybr", options={"xtol": 1e-13})
    if not sol.success:
        return None
    u_s, a_s = float(np.exp(sol.x[0])), float(np.exp(sol.x[1]))
    if not (np.isfinite(u_s) and np.isfinite(a_s) and u_s > 0.0):
        return None
    return removalTotal(u_s, a_s, p, g), u_s, a_s


def foldSurvival(p: Optional[M.Params] = None,
                 mus=(0.0, 1e-4, 1e-3, 1e-2, 3e-2, 0.1, 0.3),
                 k_mu: float = float("inf")) -> Dict[str, list]:
    """does a fold survive as dilution is turned up?

    `k_mu = inf` is constant dilution (no growth feedback); a finite `k_mu` makes
    growth slow as burden accumulates, which is what real cells do.
    """
    p = (p or M.Params()).validate()
    ceiling = M.removalCeiling(p)
    rows = []
    for mu in mus:
        g = Growth(mu0=float(mu), k_mu=k_mu)
        out = foldSolveDil(p, g)
        rows.append({
            "mu": float(mu),
            "k_mu": k_mu,
            "fold_exists": out is not None,
            "j_crit": None if out is None else out[0],
            "u_star": None if out is None else out[1],
            "a_star": None if out is None else out[2],
            "j_over_enzymatic_ceiling": None if out is None else out[0] / ceiling,
        })
    return {"ceiling": ceiling, "rows": rows}


def main() -> int:
    p = M.Params().validate()
    print("base parameters, enzymatic removal ceiling = %.4f" % M.removalCeiling(p))
    print()

    print("REDUCTION  mu = 0 must reproduce the frozen model exactly")
    g0 = Growth(mu0=0.0)
    worst = 0.0
    for u, a in ((0.05, 0.01), (0.3, 0.1), (1.0, 0.5)):
        d0, d1 = M.rhs(u, a, p), rhsDil(u, a, p, g0)
        worst = max(worst, abs(d0[0] - d1[0]), abs(d0[1] - d1[1]))
    print("  max |rhs - rhs_dil| at mu=0 : %.3e" % worst)
    print()

    print("IDENTITY  det J_dil == -(grad R_tot x grad G_dil)")
    for label, g in (("constant mu=0.02", Growth(0.02)),
                     ("feedback mu0=0.2, k_mu=0.5", Growth(0.2, 0.5))):
        errs = [determinantIdentityDil(u, a, p, g)["rel_err"]
                for u, a in ((0.05, 0.01), (0.3, 0.1), (1.0, 0.5))]
        print("  %-28s median rel err %.3e" % (label, float(np.median(errs))))
    print()

    print("ANALYTIC vs NUMERICAL jacobian of the diluted field")
    for label, g in (("constant mu=0.02", Growth(0.02)),
                     ("feedback mu0=0.2, k_mu=0.5", Growth(0.2, 0.5))):
        e = max(float(np.max(np.abs(jacobianDil(u, a, p, g)
                                    - numericalJacobianDil(u, a, p, g))))
                for u, a in ((0.05, 0.01), (0.3, 0.1)))
        print("  %-28s max abs diff %.3e" % (label, e))
    print()

    print("FOLD SURVIVAL under CONSTANT dilution (no growth feedback)")
    r = foldSurvival(p)
    print("   %-10s %-12s %-14s %s" % ("mu", "fold?", "j_crit", "j_crit/ceiling"))
    for row in r["rows"]:
        print("   %-10.4g %-12s %-14s %s" % (
            row["mu"], row["fold_exists"],
            "-" if row["j_crit"] is None else "%.6f" % row["j_crit"],
            "-" if row["j_crit"] is None else "%.4f" % row["j_over_enzymatic_ceiling"]))
    print()

    print("FOLD SURVIVAL with GROWTH FEEDBACK (k_mu = 0.5)")
    r2 = foldSurvival(p, k_mu=0.5)
    print("   %-10s %-12s %-14s %s" % ("mu0", "fold?", "j_crit", "j_crit/ceiling"))
    for row in r2["rows"]:
        print("   %-10.4g %-12s %-14s %s" % (
            row["mu"], row["fold_exists"],
            "-" if row["j_crit"] is None else "%.6f" % row["j_crit"],
            "-" if row["j_crit"] is None else "%.4f" % row["j_over_enzymatic_ceiling"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
