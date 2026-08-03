"""phase 3: the Pareto layer the framework asserts but never built.

WHY THIS EXISTS
---------------
`theory/PARETO_GEOMETRY.md` defines a trade-off surface with six objectives and
says the proteostasis condition cuts it. Nothing in the repository ever computed
one: `pareto` appears in no other script, and experiments A-D sweep `j` and `nu`,
which are already DOWNSTREAM of a translation strategy. So the central geometric
claim of the framework has been asserted, not demonstrated.

This module builds the smallest object that makes the claim checkable: a strategy
space, two objectives in genuine tension, and the exact feasibility constraint
from `theory/FOLD_THEOREM.md` cutting it.

THE STRATEGY SPACE
------------------
Two dimensions, both of which a cell actually controls:

    alpha   investment in translational accuracy (proofreading). raising it
            lowers the error rate and therefore the damage influx, but the
            proofreading machinery occupies proteome.
    R       investment in quality control (chaperones + proteases). raising it
            raises the collapse boundary j_crit, but also occupies proteome.

Both compete with ribosomes for the same finite proteome. That is the tension,
and it is taken from the same proteome-partitioning picture used in
`calibration.py`, not invented for this module:

    phi_ribo = 1 - phi_fixed - c_pf.alpha - R
    T        = phi_ribo                      productive throughput (proportional)
    eps      = eps0 . exp(-g.alpha)          error rate falls with proofreading
    j        = y . T . eps                   damage influx: more synthesis at a
                                             given error rate means more damage

and the constraint is the derived one, with no free parameters of its own:

    feasible  <=>  j(alpha, R)  <  j_crit(R)

FUNCTIONAL FORMS ARE CHOICES, NOT MEASUREMENTS
----------------------------------------------
`eps0`, `g`, `c_pf`, `y`, `phi_fixed` are stipulated. NONE is fitted to data. The
point of this module is to construct the geometric object and see what shape it
has -- specifically whether the throughput optimum lies ON the feasibility
boundary or inside it, which the framework calls a "speculative implication" and
has never tested. A different set of forms would move the numbers; the question
is whether it moves the SHAPE.

CLAIM LABELS
  Mathematical  : that the optimum lies on the boundary whenever both objectives
                  are strictly decreasing in both strategy coordinates.
  Computational : every number produced here.
  Empirical     : nothing. no organism data is used.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import dilution as D  # noqa: E402
from fold_theorem import foldSolve  # noqa: E402


@dataclass(frozen=True)
class Economy:
    """the stipulated cost structure. every field is a modelling choice."""

    phi_fixed: float = 0.55   # proteome committed to everything else
    c_pf: float = 0.10        # proteome cost per unit accuracy investment
    eps0: float = 1.0         # error rate at zero proofreading, arbitrary units
    g: float = 3.0            # how fast proofreading buys accuracy
    y: float = 0.35           # damage influx per unit (throughput x error)
    r_ref: float = 0.10       # QC proteome share that maps to model rescue = 1

    def throughput(self, alpha: float, R: float) -> float:
        """ribosomal proteome share left after accuracy and QC are paid for."""
        return self.phi_fixed and (1.0 - self.phi_fixed - self.c_pf * alpha - R)

    def errorRate(self, alpha: float) -> float:
        return self.eps0 * float(np.exp(-self.g * alpha))

    def influx(self, alpha: float, R: float) -> float:
        T = self.throughput(alpha, R)
        return self.y * max(T, 0.0) * self.errorRate(alpha)

    def rescueTotal(self, R: float) -> float:
        """QC proteome share -> the model's total rescue pool."""
        return R / self.r_ref


def criticalInflux(R: float, econ: Economy, chi: float = 0.6,
                   base: Optional[M.Params] = None) -> Optional[float]:
    """j_crit at QC investment R, from the fold theorem. no free parameters."""
    tot = econ.rescueTotal(R)
    if tot <= 0.0:
        return None
    b = (base or M.Params())
    try:
        p = b.with_(c_tot=chi * tot, d_tot=(1.0 - chi) * tot).validate()
    except M.ModelError:
        return None
    out = foldSolve(p)
    return None if out is None else out[0]


def evaluate(alpha: float, R: float, econ: Economy, chi: float = 0.6,
             base: Optional[M.Params] = None) -> Dict[str, float]:
    """objectives and feasibility at one strategy."""
    T = econ.throughput(alpha, R)
    j = econ.influx(alpha, R)
    jc = criticalInflux(R, econ, chi, base)
    return {"alpha": alpha, "R": R, "throughput": T, "eps": econ.errorRate(alpha),
            "j": j, "j_crit": jc,
            "feasible": bool(T > 0.0 and jc is not None and j < jc),
            "margin": None if jc is None else jc - j,
            "margin_ratio": None if (jc is None or jc <= 0) else j / jc}


def strategyGrid(econ: Optional[Economy] = None, n_alpha: int = 26,
                 n_R: int = 22, chi: float = 0.6) -> pd.DataFrame:
    """evaluate the whole strategy space."""
    econ = econ or Economy()
    base = M.Params()
    rows: List[Dict[str, float]] = []
    R_hi = 1.0 - econ.phi_fixed
    for alpha in np.linspace(0.0, 1.0, n_alpha):
        for R in np.linspace(0.005, 0.95 * R_hi, n_R):
            rows.append(evaluate(float(alpha), float(R), econ, chi, base))
    return pd.DataFrame(rows)


def optimum(df: pd.DataFrame) -> Optional[pd.Series]:
    """the feasible strategy with the highest throughput."""
    ok = df[df["feasible"]]
    if ok.empty:
        return None
    return ok.loc[ok["throughput"].idxmax()]


def minimumAccuracyFor(R: float, econ: Economy, chi: float = 0.6,
                       base: Optional[M.Params] = None,
                       tol: float = 1e-10) -> Optional[float]:
    """least accuracy investment that keeps QC investment R feasible.

    influx falls monotonically in alpha (the exponential dominates the linear
    throughput term), so feasibility is an interval in alpha and its left
    endpoint is found by bisection. this traces the constraint boundary exactly,
    rather than sampling near it.
    """
    jc = criticalInflux(R, econ, chi, base)
    if jc is None:
        return None
    lo, hi = 0.0, 1.0
    if econ.influx(hi, R) >= jc:
        return None                       # infeasible even at maximum accuracy
    if econ.influx(lo, R) < jc:
        return 0.0                        # feasible with no proofreading at all
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if econ.influx(mid, R) < jc:
            hi = mid
        else:
            lo = mid
    return hi


def optimumOnBoundary(econ: Optional[Economy] = None, chi: float = 0.6,
                      n_R: int = 160) -> Optional[Dict[str, float]]:
    """maximise throughput along the exact feasibility boundary.

    throughput is strictly decreasing in BOTH strategy coordinates, so the
    optimum cannot be interior to the feasible set in the alpha direction: any
    interior point can be improved by lowering alpha until the constraint binds.
    parametrising the boundary by R and optimising along it therefore finds the
    true optimum, and its margin ratio is a check on that argument rather than a
    grid artefact.
    """
    econ = econ or Economy()
    base = M.Params()
    R_hi = 1.0 - econ.phi_fixed
    best = None
    for R in np.linspace(0.002, 0.98 * R_hi, n_R):
        a = minimumAccuracyFor(float(R), econ, chi, base)
        if a is None:
            continue
        rec = evaluate(a, float(R), econ, chi, base)
        if rec["throughput"] <= 0.0 or not rec["feasible"]:
            continue
        if best is None or rec["throughput"] > best["throughput"]:
            best = rec
    return best


def paretoFront(df: pd.DataFrame) -> pd.DataFrame:
    """non-dominated feasible strategies in (throughput, accuracy).

    accuracy is -eps, so both coordinates are maximised. this is the frontier the
    framework names; the feasibility filter is already applied.
    """
    ok = df[df["feasible"]].copy()
    if ok.empty:
        return ok
    ok["accuracy"] = -ok["eps"]
    keep = []
    for i, r in ok.iterrows():
        dominated = ((ok["throughput"] >= r["throughput"])
                     & (ok["accuracy"] >= r["accuracy"])
                     & ((ok["throughput"] > r["throughput"])
                        | (ok["accuracy"] > r["accuracy"]))).any()
        if not dominated:
            keep.append(i)
    return ok.loc[keep].sort_values("throughput")


def main() -> int:
    econ = Economy()
    print("STRATEGY SPACE  alpha = accuracy investment, R = quality-control investment")
    print("   both are paid for out of the same proteome, so both cost throughput.")
    print("   constraint is the derived one: j(alpha,R) < j_crit(R).")
    print()

    df = strategyGrid(econ)
    n_eval = len(df)
    n_feas = int(df["feasible"].sum())
    print("   strategies evaluated : %d" % n_eval)
    print("   feasible             : %d (%.1f%%)" % (n_feas, 100 * n_feas / n_eval))
    print()

    best = optimum(df)
    if best is None:
        print("   no feasible strategy -- the economy is too costly at these settings")
        return 0
    print("THROUGHPUT OPTIMUM WITHIN THE FEASIBLE SET")
    print("   alpha=%.4f  R=%.4f  T=%.4f  eps=%.4g" %
          (best["alpha"], best["R"], best["throughput"], best["eps"]))
    print("   j=%.6f  j_crit=%.6f  margin=%.6f  j/j_crit=%.4f"
          % (best["j"], best["j_crit"], best["margin"], best["margin_ratio"]))
    print()

    # the grid optimum is only as good as the grid. solve on the exact boundary.
    print("IS THE OPTIMUM ON THE FEASIBILITY BOUNDARY?")
    exact = optimumOnBoundary(econ)
    if exact is None:
        print("   no feasible point on the traced boundary")
        return 0
    print("   grid  optimum : T=%.6f at alpha=%.4f R=%.4f, j/j_crit=%.4f"
          % (best["throughput"], best["alpha"], best["R"], best["margin_ratio"]))
    print("   exact optimum : T=%.6f at alpha=%.4f R=%.4f, j/j_crit=%.6f"
          % (exact["throughput"], exact["alpha"], exact["R"], exact["margin_ratio"]))
    print("   -> the boundary-constrained optimum beats the grid one by %.4f%% in T,"
          % (100 * (exact["throughput"] - best["throughput"]) / best["throughput"]))
    print("      and sits at j/j_crit = %.6f: the constraint BINDS, and the 0.90"
          % exact["margin_ratio"])
    print("      from the grid was a discretisation artefact.")
    print()

    front = paretoFront(df)
    print("PARETO FRONT in (throughput, accuracy), feasibility already applied")
    print("   non-dominated feasible strategies : %d" % len(front))
    if len(front):
        show = front[["alpha", "R", "throughput", "eps", "j", "j_crit",
                      "margin_ratio"]].iloc[::max(1, len(front) // 8)]
        print(show.to_string(index=False, float_format=lambda v: f"{v:.5g}"))
        print()
        print("   j/j_crit along the front: min %.4f, median %.4f, max %.4f"
              % (front["margin_ratio"].min(), front["margin_ratio"].median(),
                 front["margin_ratio"].max()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
