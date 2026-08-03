"""phase 3: spatial sequestration, the mechanism named by a FAILED post-diction.

PROVENANCE, which matters here more than usual
----------------------------------------------
D026 tried to explain Lindner et al. 2008 (aging and rejuvenation by asymmetric
segregation of aggregates) and failed. The failure named this mechanism: the
observation is not merely that a high-burden state exists, it is that aggregates
are SPATIALLY SEQUESTERED at a pole. The two-state model is well-mixed.

This is the first mechanism in the project that arrived from an observation
rather than from inspecting the model, and the prediction was written down before
this file was run -- see D028, committed in the preceding commit, and
`notes/POSTDICTION_PROTOCOL.md` for why.

THE MECHANISM
-------------
Aggregate splits into REACTIVE `a_r` and SEQUESTERED `a_s`:

    du/dt   = j - v_ref - v_degU - v_nuc - v_grow + v_dis          ( - mu.u )
    da_r/dt = v_nuc + v_grow - v_dis - v_degA
              - k_seq.a_r^q + k_rel.a_s                            ( - mu.a_r )
    da_s/dt = k_seq.a_r^q - k_rel.a_s                              ( - mu.a_s )

every flux on the first two lines is evaluated at `(u, a_r)`, never at
`a_r + a_s`. That is the entire mechanism: `a_s` does not bind chaperone or
protease, does not nucleate, does not grow, and is not cleared -- it can only
leave by releasing back to `a_r`. Sequestration changes KINETICS, not
bookkeeping, so mass balance is untouched:

    du/dt + da_r/dt + da_s/dt = j - R            ( - mu.(u + a_r + a_s) )

because the `k_seq`/`k_rel` pair cancels between the two aggregate equations.
`k_seq = 0` reduces to the two-state model exactly.

ONE MODELLING CHOICE, DECLARED
------------------------------
Growth rate is taken to respond to TOTAL burden `u + a_r + a_s`, i.e. sequestered
aggregate still costs growth. It is protein diverted from the functional
proteome whether or not it is soluble, and Lindner's old-pole cells lose
reproductive ability while visibly carrying inclusion bodies. The alternative --
sequestration as a free lunch, `mu = f(u + a_r)` -- is a different mechanism and
is run too, and reported, per protocol rule 3.

D024 APPLIES, BUT IT IS VERIFIED RATHER THAN ASSUMED
----------------------------------------------------
Three states with influx in one equation and mass balance intact, so

    det J = -det[ grad R ; grad G ; grad C ]

with `G = da_r/dt` and `C = da_s/dt`. `identity3s` checks it BEFORE any number
from this model is interpreted, per D028.

CLAIM LABELS
  Mathematical  : the identity, inherited from D024.
  Computational : every number below.
  Empirical     : the comparison against Lindner et al., which is a post-diction
                  against a published measurement, not a new observation.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import dilution as D  # noqa: E402
import boundary_structure as B  # noqa: E402
import fold_theorem as FT  # noqa: E402
from postdiction_aging import LINDNER2008  # noqa: E402

# the preregistered match band (D028). BOTH edges are declared in advance; the
# upper one is the weaker of the two and D028 says so.
BAND_LO, BAND_HI = 0.30, 0.60

# regimes disqualified in advance (D028, protocol rule 4)
DISQUALIFIED = ("constant",)


@dataclass(frozen=True)
class Sequestration:
    """`k_seq = 0` is the two-state model, exactly."""

    k_seq: float = 0.0
    k_rel: float = 0.0
    q: float = 1.0          # sequestration order; 1 linear, >1 superlinear

    def validate(self) -> "Sequestration":
        for name in ("k_seq", "k_rel"):
            v = getattr(self, name)
            if not (v >= 0.0 and np.isfinite(v)):
                raise M.ModelError(f"'{name}' must be finite and nonnegative")
        if not (self.q >= 1.0 and np.isfinite(self.q)):
            raise M.ModelError("sequestration order q must be finite and >= 1")
        return self

    def fluxes(self, a_r: float, a_s: float) -> Tuple[float, float]:
        """(sequestration flux, release flux). both nonnegative."""
        if self.k_seq == 0.0 and self.k_rel == 0.0:
            return 0.0, 0.0
        ar = max(a_r, 0.0)
        seq = self.k_seq * (ar ** self.q if self.q != 1.0 else ar)
        return float(seq), float(self.k_rel * max(a_s, 0.0))


# ---------------------------------------------------------------------------
# the three-state field
# ---------------------------------------------------------------------------


def rhs3s(u: float, a_r: float, a_s: float, p: M.Params, seq: Sequestration,
          g: Optional[D.Growth] = None,
          s_costs_growth: bool = True) -> Tuple[float, float, float]:
    """(du/dt, da_r/dt, da_s/dt). all reactive fluxes evaluated at (u, a_r)."""
    du, da = M.rhs(u, a_r, p)
    fs, fr = seq.fluxes(a_r, a_s)
    du_, dar, das = du, da - fs + fr, fs - fr
    if g is not None and g.mu0 != 0.0:
        burden_a = a_r + a_s if s_costs_growth else a_r
        mu = g.rate(u, burden_a)
        du_ -= mu * u
        dar -= mu * a_r
        das -= mu * a_s
    return float(du_), float(dar), float(das)


def removalR3s(u: float, a_r: float, a_s: float, p: M.Params,
               g: Optional[D.Growth] = None,
               s_costs_growth: bool = True) -> float:
    """total removal. dilution of ALL THREE states counts as disposal (D010)."""
    R = FT.removalR(u, a_r, p)
    if g is not None and g.mu0 != 0.0:
        burden_a = a_r + a_s if s_costs_growth else a_r
        R += g.rate(u, burden_a) * (u + a_r + a_s)
    return float(R)


def aggregateG3s(u: float, a_r: float, a_s: float, p: M.Params,
                 seq: Sequestration, g=None, s_costs_growth: bool = True) -> float:
    return rhs3s(u, a_r, a_s, p, seq, g, s_costs_growth)[1]


def controllerC3s(u: float, a_r: float, a_s: float, p: M.Params,
                  seq: Sequestration, g=None, s_costs_growth: bool = True) -> float:
    return rhs3s(u, a_r, a_s, p, seq, g, s_costs_growth)[2]


def _grad3(fn, x, h_rel: float = 1e-4) -> np.ndarray:
    out = []
    for i in range(3):
        h = h_rel * max(abs(x[i]), 1e-8)
        xp, xm = list(x), list(x)
        xp[i] += h
        xm[i] -= h
        out.append((fn(*xp) - fn(*xm)) / (2.0 * h))
    return np.array(out)


def jacobian3s(u: float, a_r: float, a_s: float, p: M.Params,
               seq: Sequestration, g=None, s_costs_growth: bool = True,
               h_rel: float = 1e-4) -> np.ndarray:
    x = [u, a_r, a_s]
    cols = []
    for i in range(3):
        h = h_rel * max(abs(x[i]), 1e-8)
        xp, xm = list(x), list(x)
        xp[i] += h
        xm[i] -= h
        fp = np.array(rhs3s(*xp, p, seq, g, s_costs_growth))
        fm = np.array(rhs3s(*xm, p, seq, g, s_costs_growth))
        cols.append((fp - fm) / (2.0 * h))
    return np.column_stack(cols)


def identity3s(u: float, a_r: float, a_s: float, p: M.Params,
               seq: Sequestration, g=None,
               s_costs_growth: bool = True) -> Dict[str, float]:
    """D024 on the extended system: det J == -det[grad R; grad G; grad C].

    verified BEFORE anything else here is interpreted, per D028.
    """
    detJ = float(np.linalg.det(jacobian3s(u, a_r, a_s, p, seq, g, s_costs_growth)))
    gR = _grad3(lambda a, b, c: removalR3s(a, b, c, p, g, s_costs_growth),
                [u, a_r, a_s])
    gG = _grad3(lambda a, b, c: aggregateG3s(a, b, c, p, seq, g, s_costs_growth),
                [u, a_r, a_s])
    gC = _grad3(lambda a, b, c: controllerC3s(a, b, c, p, seq, g, s_costs_growth),
                [u, a_r, a_s])
    rhs_ = -float(np.linalg.det(np.vstack([gR, gG, gC])))
    scale = max(abs(detJ), abs(rhs_), 1e-300)
    return {"det_J": detJ, "minus_det_grads": rhs_,
            "rel_err": abs(detJ - rhs_) / scale}


# ---------------------------------------------------------------------------
# attractors and the measured quantity
# ---------------------------------------------------------------------------


def settle3(j: float, x0, p: M.Params, seq: Sequestration, g: D.Growth,
            s_costs_growth: bool = True, t_end: float = 5e4):
    from scipy.integrate import solve_ivp
    pj = p.with_(j=float(j)).validate()

    def field(_t, x):
        return list(rhs3s(max(x[0], 0.0), max(x[1], 0.0), max(x[2], 0.0),
                          pj, seq, g, s_costs_growth))

    s = solve_ivp(field, [0.0, t_end], list(x0), method="Radau",
                  rtol=1e-10, atol=1e-13)
    return tuple(float(s.y[i, -1]) for i in range(3))


def _isEquilibrium(x, j, p, seq, g, s_costs_growth) -> bool:
    """D016: a large state is not an attractor. check against the vector field."""
    if not all(np.isfinite(v) for v in x):
        return False
    pj = p.with_(j=float(j)).validate()
    try:
        f = rhs3s(x[0], x[1], x[2], pj, seq, g, s_costs_growth)
    except (M.ModelError, OverflowError):
        return False
    return max(abs(v) for v in f) <= 1e-7 * max(1.0, sum(abs(v) for v in x))


def attractorPair(j: float, p: M.Params, seq: Sequestration, g: D.Growth,
                  s_costs_growth: bool = True) -> Optional[Dict[str, object]]:
    """locate the low- and high-burden attractors and the reproductive loss.

    the quantity returned as `reproductive_loss` is `1 - mu_high/mu_low`, the
    OLD-POLE lineage against the NEW-POLE lineage, which is the form Lindner et
    al. report ("relative to the new-pole progeny").
    """
    low = settle3(j, (0.0, 0.0, 0.0), p, seq, g, s_costs_growth)
    big = max(50.0 * max(low[1], 1e-6), 5.0)
    high = settle3(j, (low[0], big, big), p, seq, g, s_costs_growth)
    if not (_isEquilibrium(low, j, p, seq, g, s_costs_growth)
            and _isEquilibrium(high, j, p, seq, g, s_costs_growth)):
        return None

    tot_lo = low[1] + low[2]
    tot_hi = high[1] + high[2]
    bistable = tot_hi > 1.5 * max(tot_lo, 1e-12)
    ba_lo = tot_lo if s_costs_growth else low[1]
    ba_hi = tot_hi if s_costs_growth else high[1]
    m_lo, m_hi = g.rate(low[0], ba_lo), g.rate(high[0], ba_hi)
    if m_lo <= 0.0:
        return None
    return {"j": j, "bistable": bool(bistable),
            "a_total_low": tot_lo, "a_total_high": tot_hi,
            "a_seq_high": high[2],
            "seq_share_high": high[2] / max(tot_hi, 1e-30),
            "aggregate_ratio": tot_hi / max(tot_lo, 1e-12),
            "reproductive_loss": float(1.0 - m_hi / m_lo) if bistable else 0.0}


def inBand(loss: float) -> bool:
    """the preregistered verdict, D028. not to be widened after the fact."""
    return bool(BAND_LO <= loss <= BAND_HI)


# ---------------------------------------------------------------------------
# the run
# ---------------------------------------------------------------------------

# every regime tried is listed here and every one is reported, per protocol
# rule 3. "constant" is disqualified in advance and is run only to show what it
# would have claimed.
GROWTH_LAWS = {
    "constant": lambda mu0, k: D.Growth(mu0, float("inf")),
    "hyperbolic": lambda mu0, k: D.Growth(mu0, k),
    "linear_arrest": lambda mu0, k: B.LinearGrowth(mu0, k),
}

SEQ_GRID = (
    Sequestration(0.0, 0.0, 1.0),        # control: the two-state model
    Sequestration(0.1, 0.01, 1.0),
    Sequestration(1.0, 0.01, 1.0),
    Sequestration(1.0, 0.1, 1.0),
    Sequestration(10.0, 0.1, 1.0),
    Sequestration(0.1, 0.01, 2.0),       # superlinear
    Sequestration(1.0, 0.01, 2.0),
    Sequestration(1.0, 0.1, 2.0),
)

GROWTH_GRID = ((0.05, 0.5), (0.1, 0.5), (0.2, 0.5), (0.1, 1.0), (0.3, 1.0),
               (0.2, 2.0))


def scan(p: Optional[M.Params] = None, s_costs_growth: bool = True,
         laws=tuple(GROWTH_LAWS), seq_grid=SEQ_GRID,
         growth_grid=GROWTH_GRID, fracs=(0.90, 0.95, 0.98, 0.995)) -> pd.DataFrame:
    """every (growth law, sequestration setting, influx) cell. nothing dropped."""
    p = (p or M.Params()).validate()
    rows: List[Dict[str, object]] = []
    for law in laws:
        for mu0, k_mu in growth_grid:
            g = GROWTH_LAWS[law](mu0, k_mu)
            out = D.foldSolveDil(p, g)
            if out is None:
                rows.append({"law": law, "mu0": mu0, "k_mu": k_mu,
                             "k_seq": np.nan, "k_rel": np.nan, "q": np.nan,
                             "j_over_jcrit": np.nan, "fold": False,
                             "settled": False, "bistable": False,
                             "reproductive_loss": np.nan, "in_band": False})
                continue
            j_crit = out[0]
            for seq in seq_grid:
                for frac in fracs:
                    j = frac * j_crit
                    try:
                        r = attractorPair(j, p, seq, g, s_costs_growth)
                    except (M.ModelError, ValueError, OverflowError):
                        r = None
                    base = {"law": law, "mu0": mu0, "k_mu": k_mu,
                            "k_seq": seq.k_seq, "k_rel": seq.k_rel, "q": seq.q,
                            "j_over_jcrit": frac, "fold": True}
                    if r is None:
                        rows.append({**base, "settled": False, "bistable": False,
                                     "reproductive_loss": np.nan, "in_band": False})
                        continue
                    rows.append({**base, "settled": True,
                                 "bistable": r["bistable"],
                                 "a_total_low": r["a_total_low"],
                                 "a_total_high": r["a_total_high"],
                                 "seq_share_high": r["seq_share_high"],
                                 "aggregate_ratio": r["aggregate_ratio"],
                                 "reproductive_loss": r["reproductive_loss"],
                                 "in_band": bool(r["bistable"]
                                                 and inBand(r["reproductive_loss"]))})
    return pd.DataFrame(rows)


def verdict(df: pd.DataFrame) -> Dict[str, object]:
    """the preregistered verdict. a pass requires a QUALIFIED growth law."""
    ok = df[~df["law"].isin(DISQUALIFIED)]
    dis = df[df["law"].isin(DISQUALIFIED)]
    return {
        "n_cells": int(len(df)),
        "n_qualified": int(len(ok)),
        "qualified_bistable": int(ok["bistable"].sum()),
        "qualified_in_band": int(ok["in_band"].sum()),
        "disqualified_in_band": int(dis["in_band"].sum()),
        "passes": bool(ok["in_band"].any()),
        "band": (BAND_LO, BAND_HI),
        "observed": LINDNER2008["reproductive_loss_old_pole"],
    }


def main() -> int:
    p = M.Params().validate()
    print("POST-DICTION: spatial sequestration vs Lindner et al. 2008")
    print("  doi:%s (PMID %s), via PubMed" % (LINDNER2008["doi"], LINDNER2008["pmid"]))
    print("  preregistered band (D028): [%.2f, %.2f] on 1 - mu_high/mu_low"
          % (BAND_LO, BAND_HI))
    print("  disqualified in advance:", ", ".join(DISQUALIFIED))

    print("\n--- D024 identity on the extended system (verified before anything else)")
    errs = []
    for seq in SEQ_GRID:
        for g in (None, D.Growth(0.1, 0.5), B.LinearGrowth(0.1, 0.5)):
            errs.append(identity3s(0.12, 0.08, 0.05, p, seq, g)["rel_err"])
    print("  n = %d   median rel_err = %.3e   max = %.3e"
          % (len(errs), float(np.median(errs)), float(np.max(errs))))

    for costs in (True, False):
        tag = ("a_s costs growth" if costs else
               "a_s is FREE (alternative mechanism, reported per rule 3)")
        print("\n=== regime: %s" % tag)
        df = scan(p, s_costs_growth=costs)
        v = verdict(df)
        print("  cells %d (qualified %d)  bistable %d  in band %d  "
              "| disqualified-law in band %d"
              % (v["n_cells"], v["n_qualified"], v["qualified_bistable"],
                 v["qualified_in_band"], v["disqualified_in_band"]))
        s = df[df["settled"] & df["bistable"]]
        if not s.empty:
            print(s.groupby(["law", "q"])["reproductive_loss"]
                   .agg(["count", "min", "median", "max"]).to_string())
        print("  VERDICT:", "PASS" if v["passes"] else "FAIL")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
