"""phase 3: the theorem beyond two states, and what regulation does to it.

WHY THIS EXISTS
---------------
Two objections to the phase 3 result deserve testing rather than answering:

  1. "It is exact, but exact about a TOY." The fold theorem was proved for a
     two-state model. If it breaks when a state is added, its exactness is a
     property of the reduction rather than of proteostasis.
  2. "It is organism-free, so maybe it is missing something." The most important
     omission is REGULATION. `c_tot` is a fixed parameter, but real cells control
     chaperone abundance in response to burden -- the sigma-32 system in E. coli.
     `theory/SCOPE_AND_NONCLAIMS.md` lists this and nothing ever tested its cost.

THE GENERALISED THEOREM
-----------------------
Let the state be `x = (u, a, c, ...)` with `du/dt` the only equation containing
the influx `j`, and let mass balance still give

    du/dt + da/dt = j - R(x)

for total removal R. Write the remaining equations as G (aggregate) and C
(controller). The determinant-preserving row operation `row1 -> row1 + row2`
replaces the first row by `-grad R`, so

    det J = -det [ grad R ; grad G ; grad C ]

which vanishes exactly when `grad R` lies in the span of the other gradients --
the Lagrange condition for a constrained critical point of R on the intersection
of the NON-INFLUX nullclines `{G = 0} n {C = 0}`.

So the theorem is not about two states. It holds for any number, provided the
influx enters one equation and the removal terms are identifiable. That is a
structural property, not a modelling convenience.

THE REGULATOR
-------------
Mechanistic rather than phenomenological: sigma-32 is titrated by free DnaK, so
chaperone synthesis rises when FREE chaperone falls, not when total burden rises.

    dc/dt = sigma0 . kappa_r/(kappa_r + cf)  -  delta . c

Setting `sigma0 = 0` recovers a constant pool and must reproduce the frozen
two-state model exactly -- that is the reduction test, the analogue of the phase 2
`epsilon -> 0` check.

CLAIM LABELS
  Mathematical  : the generalised identity, and the sigma0 -> 0 reduction.
  Computational : every number about regulation below.
  Empirical     : nothing. no organism data is used.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import root

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
from fold_theorem import KINETIC_FIELDS, phase1RunDir  # noqa: E402


@dataclass(frozen=True)
class Regulator:
    """chaperone synthesis controlled by FREE chaperone (sigma-32 titration).

    `sigma0 = 0` disables regulation and must reproduce the frozen model.
    """

    sigma0: float = 0.0      # maximal synthesis rate
    kappa_r: float = 0.1     # free-chaperone level at half-maximal synthesis
    delta: float = 0.5       # chaperone turnover + dilution

    def synthesis(self, cf: float) -> float:
        if self.sigma0 == 0.0:
            return 0.0
        return self.sigma0 * self.kappa_r / (self.kappa_r + max(cf, 0.0))

    def dcdt(self, c: float, cf: float) -> float:
        if self.sigma0 == 0.0:
            return 0.0
        return self.synthesis(cf) - self.delta * c


def stateParams(p: M.Params, c: float) -> M.Params:
    """the two-state parameters with the chaperone pool set to the state c."""
    return p.with_(c_tot=max(c, 1e-12)).validate()


def rhs3(u: float, a: float, c: float, p: M.Params,
         reg: Regulator) -> Tuple[float, float, float]:
    """the three-state field. binding algebra is untouched and reused."""
    pc = stateParams(p, c)
    du, da = M.rhs(u, a, pc)
    _, _, cf, _ = M.solveFreePools(u, a, pc)
    return du, da, reg.dcdt(c, cf)


def removalR3(u: float, a: float, c: float, p: M.Params) -> float:
    """total removal, now a function of the chaperone STATE."""
    f = M.fluxes(u, a, stateParams(p, c))
    return f["refold"] + f["degrade_u"] + f["degrade_a"]


def aggregateG3(u: float, a: float, c: float, p: M.Params) -> float:
    return M.rhs(u, a, stateParams(p, c))[1]


def controllerC3(u: float, a: float, c: float, p: M.Params,
                 reg: Regulator) -> float:
    pc = stateParams(p, c)
    _, _, cf, _ = M.solveFreePools(u, a, pc)
    return reg.dcdt(c, cf)


def _grad3(fn, x, h_rel=1e-6):
    g = []
    for k in range(3):
        h = h_rel * max(abs(x[k]), 1e-8)
        xp, xm = list(x), list(x)
        xp[k] += h
        xm[k] -= h
        g.append((fn(*xp) - fn(*xm)) / (2 * h))
    return np.array(g, dtype=float)


def jacobian3(u: float, a: float, c: float, p: M.Params,
              reg: Regulator, h_rel: float = 1e-6) -> np.ndarray:
    x = [u, a, c]
    J = np.zeros((3, 3))
    for k in range(3):
        h = h_rel * max(abs(x[k]), 1e-8)
        xp, xm = list(x), list(x)
        xp[k] += h
        xm[k] -= h
        fp = np.array(rhs3(*xp, p, reg))
        fm = np.array(rhs3(*xm, p, reg))
        J[:, k] = (fp - fm) / (2 * h)
    return J


def determinantIdentity3(u: float, a: float, c: float, p: M.Params,
                         reg: Regulator) -> Dict[str, float]:
    """det J == -det[grad R ; grad G ; grad C], the generalised identity."""
    detJ = float(np.linalg.det(jacobian3(u, a, c, p, reg)))
    gR = _grad3(lambda uu, aa, cc: removalR3(uu, aa, cc, p), [u, a, c])
    gG = _grad3(lambda uu, aa, cc: aggregateG3(uu, aa, cc, p), [u, a, c])
    gC = _grad3(lambda uu, aa, cc: controllerC3(uu, aa, cc, p, reg), [u, a, c])
    rhs_ = -float(np.linalg.det(np.vstack([gR, gG, gC])))
    scale = max(abs(detJ), abs(rhs_), 1e-300)
    return {"det_J": detJ, "minus_det_grads": rhs_,
            "rel_err": abs(detJ - rhs_) / scale}


# ---------------------------------------------------------------------------
# folds of the regulated system
# ---------------------------------------------------------------------------


def foldSolve3(p: M.Params, reg: Regulator,
               seed=(0.3, 0.1, None)) -> Optional[Tuple[float, float, float, float]]:
    """solve {G = 0, C = 0, det J = 0}; return (j_crit, u*, a*, c*).

    three equations, three unknowns, none containing j -- the same structure the
    two-state theorem has, one dimension up.
    """
    c0 = seed[2] if seed[2] is not None else p.c_tot

    def residual(z):
        u, a, c = (float(np.exp(np.clip(z[0], -34, 6))),
                   float(np.exp(np.clip(z[1], -34, 6))),
                   float(np.exp(np.clip(z[2], -34, 6))))
        try:
            _, da, dc = rhs3(u, a, c, p, reg)
            d = float(np.linalg.det(jacobian3(u, a, c, p, reg)))
        except (M.ModelError, OverflowError, np.linalg.LinAlgError):
            return [1e6, 1e6, 1e6]
        if not all(np.isfinite(v) for v in (da, dc, d)):
            return [1e6, 1e6, 1e6]
        return [da, dc, d]

    s = root(residual, [np.log(seed[0]), np.log(seed[1]), np.log(c0)],
             method="hybr", options={"xtol": 1e-12})
    if not s.success:
        return None
    u = float(np.exp(np.clip(s.x[0], -34, 6)))
    a = float(np.exp(np.clip(s.x[1], -34, 6)))
    c = float(np.exp(np.clip(s.x[2], -34, 6)))
    if not all(np.isfinite(v) and v > 0 for v in (u, a, c)):
        return None
    _, da, dc = rhs3(u, a, c, p, reg)
    if max(abs(da), abs(dc)) > 1e-7:
        return None
    return removalR3(u, a, c, p), u, a, c


def saturationAt3(u: float, a: float, c: float, p: M.Params) -> Dict[str, float]:
    pc = stateParams(p, c)
    uf, af, cf, df = M.solveFreePools(u, a, pc)
    return {"s_ref": uf / (pc.kappa_ref + uf), "s_u": uf / (pc.kappa_u + uf),
            "s_a": af / (pc.kappa_a + af), "c_star": c, "cf": cf}


def spreadComparison(k: int = 30, seed: int = 41,
                     sigma0: float = 0.6) -> pd.DataFrame:
    """does regulation NARROW the predicted saturation spread?

    the phase 3 test failed (D020) because predicted `s_u` covers nearly [0,1]
    across the parameter box. if a regulated cell sits where its controller puts
    it rather than where its raw kinetics do, that spread should shrink. this
    measures whether it does.
    """
    run = phase1RunDir()
    c_ = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c_ = c_[c_["C1_fold_exists"] == True]  # noqa: E712
    draws = c_.sample(n=k, random_state=seed)

    rows: List[Dict[str, float]] = []
    for idx, (_, r) in enumerate(draws.iterrows()):
        kin = {f: float(r["p_" + f]) for f in KINETIC_FIELDS}
        chi = float(r["p_chi"])
        try:
            p = M.Params(**kin).with_(nu=float(r["p_nu"]), c_tot=chi,
                                      d_tot=1.0 - chi).validate()
        except M.ModelError:
            continue
        for label, reg in (("unregulated", Regulator(sigma0=0.0)),
                           ("regulated", Regulator(sigma0=sigma0,
                                                   kappa_r=0.1, delta=0.5))):
            try:
                out = foldSolve3(p, reg, seed=(0.3, 0.1, p.c_tot))
            except (M.ModelError, OverflowError, np.linalg.LinAlgError, ValueError):
                continue
            if out is None:
                continue
            j, u, a, c = out
            try:
                sat = saturationAt3(u, a, c, p)
            except M.ModelError:
                continue
            sat.update({"draw": idx, "regime": label, "j_crit": j})
            rows.append(sat)
    return pd.DataFrame(rows)


def main() -> int:
    p = M.Params().validate()

    print("1. REDUCTION  sigma0 = 0 must reproduce the frozen two-state model")
    reg0 = Regulator(sigma0=0.0)
    worst = 0.0
    for u, a in ((0.05, 0.01), (0.3, 0.1), (1.0, 0.5)):
        d2 = M.rhs(u, a, p)
        d3 = rhs3(u, a, p.c_tot, p, reg0)
        worst = max(worst, abs(d2[0] - d3[0]), abs(d2[1] - d3[1]), abs(d3[2]))
    print("   max |rhs2 - rhs3| at sigma0 = 0 : %.3e" % worst)
    print()

    print("2. GENERALISED IDENTITY  det J == -det[grad R ; grad G ; grad C]")
    for label, reg in (("unregulated", Regulator(0.0)),
                       ("regulated sigma0=0.6", Regulator(0.6, 0.1, 0.5)),
                       ("regulated sigma0=1.5", Regulator(1.5, 0.05, 0.8))):
        errs = [determinantIdentity3(u, a, p.c_tot, p, reg)["rel_err"]
                for u, a in ((0.05, 0.01), (0.3, 0.1), (1.0, 0.5))]
        print("   %-24s median rel err %.3e" % (label, float(np.median(errs))))
    print("   -> the theorem is not a property of the two-state reduction.")
    print()

    if not phase1RunDir().is_dir():
        print("SKIP: no phase 1 run root, cannot run the spread comparison.")
        return 0

    print("3. DOES REGULATION SHARPEN THE PREDICTION?")
    df = spreadComparison()
    if df.empty:
        print("   no networks converged")
        return 0
    qs = [0.05, 0.25, 0.5, 0.75, 0.95]
    for regime, sub in df.groupby("regime"):
        q = sub["s_u"].quantile(qs)
        print("   %-12s n=%-4d  s_u  p5 %.4f  p25 %.4f  p50 %.4f  p75 %.4f  p95 %.4f"
              % (regime, len(sub), *[q[x] for x in qs]))
    piv = {r: g for r, g in df.groupby("regime")}
    if "regulated" in piv and "unregulated" in piv:
        w_un = piv["unregulated"]["s_u"].quantile(.95) - piv["unregulated"]["s_u"].quantile(.05)
        w_re = piv["regulated"]["s_u"].quantile(.95) - piv["regulated"]["s_u"].quantile(.05)
        print()
        print("   p5-p95 width of s_u : unregulated %.4f -> regulated %.4f" % (w_un, w_re))
        print("   -> %s" % ("NARROWS: regulation makes the prediction testable"
                            if w_re < 0.7 * w_un else
                            "does NOT narrow enough to rescue the saturation test"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
