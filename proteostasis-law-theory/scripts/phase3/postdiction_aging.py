"""phase 3 post-diction: asymmetric division, aging, and rejuvenation.

THE OBSERVATION (not ours, and it long predates this model)
-----------------------------------------------------------
Lindner AB, Madden R, Demarez A, Stewart EJ, Taddei F (2008) "Asymmetric
segregation of protein aggregates is associated with cellular aging and
rejuvenation." PNAS 105(8):3076-3081, doi:10.1073/pnas.0708931105, PMID 18287048,
via PubMed. E. coli growing under NON-STRESSED conditions spontaneously
accumulate protein aggregates; these partition to old-pole cells, which lose
">30% of reproductive ability"; the new-pole progeny, "devoid of parental
inclusion bodies, exhibit rejuvenation."

WHY IT IS A TEST OF THIS MODEL RATHER THAN A DECORATION
-------------------------------------------------------
"Rejuvenation" is only a meaningful category if the system is BISTABLE. In a
monostable system, halving a daughter's aggregate load merely places it lower on
the same continuous branch, and it relaxes straight back to the one attractor --
there is nothing to be rejuvenated INTO. Two attractors separated by a separatrix
are what make inheritance of a low-burden state possible at all.

The model produces exactly that structure (D013), and -- this is the point --
only under ONE of the two candidate growth-arrest laws:

  hyperbolic arrest (growth approaches zero asymptotically) -> dilution never
      switches off -> the high-burden state is BOUNDED -> bistable.
  linear arrest (growth reaches exactly zero at finite burden) -> dilution
      switches off -> the high-burden state RUNS AWAY -> not bistable.

D015 identified that shape as the single most valuable unmeasured quantity and
noted the measurement cannot settle it, since the data constrain the slope three
decades below the arrest burden. The aging literature bears on it from a
completely different direction: an observed, heritable, viable high-aggregate
state is evidence for boundedness, hence for the hyperbolic branch.

WHAT THIS MODULE COMPUTES
-------------------------
The quantitative version. From the high-burden attractor, what fraction of
aggregate must a daughter shed to land in the low-burden basin? That is the
model's rejuvenation threshold, and it is directly comparable to a measured
segregation asymmetry.

CLAIM LABELS
  Computational : every number here.
  Empirical     : the Lindner observation, which is quoted, not produced here.
                  no organism data is analysed.
"""

from __future__ import annotations

import sys
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
import boundary_structure as B  # noqa: E402

LINDNER2008 = {
    "doi": "10.1073/pnas.0708931105",
    "pmid": "18287048",
    "reproductive_loss_old_pole": 0.30,   # ">30%" as reported
    "condition": "non-stressed growth",
}


def basinOf(u0: float, a0: float, j: float, p: M.Params, g: D.Growth,
            low_ref, high_ref, t_end: float = 5e4) -> Optional[str]:
    """which attractor the state (u0, a0) relaxes to."""
    x = B.settle(j, (u0, a0), p.with_(j=float(j)).validate() and p, g, t_end)
    if not all(np.isfinite(v) for v in x):
        return None
    d_low = abs(x[1] - low_ref[1]) / max(low_ref[1], 1e-12)
    d_high = abs(x[1] - high_ref[1]) / max(high_ref[1], 1e-12)
    return "low" if d_low < d_high else "high"


def rejuvenationThreshold(p: M.Params, g: D.Growth, j: float,
                          n_bisect: int = 26) -> Optional[Dict[str, float]]:
    """fraction of aggregate a daughter must shed to escape the high state.

    both attractors are located first by settling from zero burden (low) and from
    a large aggregate load (high). the daughter's state is the high attractor
    with aggregate scaled by (1 - f); f is bisected for the escape point.
    """
    pj = p.with_(j=float(j)).validate()
    low = B.settle(j, (0.0, 0.0), p, g)
    high = B.settle(j, (low[0], max(50.0 * max(low[1], 1e-6), 5.0)), p, g)
    if not all(np.isfinite(v) for v in low + high):
        return None
    if high[1] <= 1.5 * max(low[1], 1e-12):
        return None                      # monostable at this influx

    # verify both are genuine equilibria
    for st in (low, high):
        du, da = D.rhsDil(st[0], st[1], pj, g)
        if max(abs(du), abs(da)) > 1e-7 * max(1.0, st[0] + st[1]):
            return None

    lo_f, hi_f = 0.0, 1.0                # f = 0 stays high, f = 1 sheds all
    for _ in range(n_bisect):
        f = 0.5 * (lo_f + hi_f)
        landed = B.settle(j, (high[0], high[1] * (1.0 - f)), p, g)
        d_low = abs(landed[1] - low[1]) / max(low[1], 1e-12)
        d_high = abs(landed[1] - high[1]) / max(high[1], 1e-12)
        if d_low < d_high:
            hi_f = f                     # escaped
        else:
            lo_f = f                     # still high
    return {"j": j, "a_low": low[1], "a_high": high[1],
            "u_low": low[0], "u_high": high[0],
            "shed_fraction_needed": hi_f,
            "aggregate_ratio": high[1] / max(low[1], 1e-12)}


def scanInfluxes(p: Optional[M.Params] = None, g: Optional[D.Growth] = None,
                 n: int = 7) -> pd.DataFrame:
    """the rejuvenation threshold across the bistable window."""
    p = (p or M.Params()).validate()
    g = g or D.Growth(0.04)
    out = D.foldSolveDil(p, g)
    if out is None:
        return pd.DataFrame()
    j_crit = out[0]
    rows: List[Dict[str, float]] = []
    for j in np.linspace(0.83 * j_crit, 0.995 * j_crit, n):
        try:
            r = rejuvenationThreshold(p, g, float(j))
        except (M.ModelError, ValueError, OverflowError):
            continue
        if r is not None:
            rows.append(r)
    return pd.DataFrame(rows)


def main() -> int:
    p = M.Params().validate()
    print("POST-DICTION: aging and rejuvenation by asymmetric segregation")
    print("   observation: Lindner et al. 2008 PNAS, doi:%s (PMID %s), via PubMed"
          % (LINDNER2008["doi"], LINDNER2008["pmid"]))
    print("   aggregates accumulate in old-pole cells under %s, costing >%d%% of"
          % (LINDNER2008["condition"],
             100 * LINDNER2008["reproductive_loss_old_pole"]))
    print("   reproductive ability; new-pole progeny 'exhibit rejuvenation'.")
    print()
    print("   rejuvenation is only a meaningful category if the system is BISTABLE:")
    print("   in a monostable system a daughter with less aggregate simply relaxes")
    print("   back to the single attractor. so the observation bears on which")
    print("   growth-arrest law holds -- the question D015 left open.")
    print()

    print("A. HYPERBOLIC arrest (dilution never switches off) -> bounded high state")
    df = scanInfluxes(p, D.Growth(0.04))
    if df.empty:
        print("   no bistable window found")
    else:
        show = df[["j", "a_low", "a_high", "aggregate_ratio",
                   "shed_fraction_needed"]]
        print(show.to_string(index=False, float_format=lambda v: f"{v:.5g}"))
        print()
        print("   aggregate ratio high/low : %.1f - %.1f fold"
              % (df["aggregate_ratio"].min(), df["aggregate_ratio"].max()))
        print("   shed fraction to rejuvenate : %.3f - %.3f (median %.3f)"
              % (df["shed_fraction_needed"].min(),
                 df["shed_fraction_needed"].max(),
                 df["shed_fraction_needed"].median()))
    print()

    print("B. LINEAR arrest (dilution switches off) -> no bounded high state")
    import boundary_structure as BS
    dfl = scanInfluxes(p, BS.LinearGrowth(0.3, 1.5625))
    print("   bistable influxes found : %d" % len(dfl))
    print("   -> %s" % ("no bistability, so no rejuvenation is possible under this law"
                        if dfl.empty else "unexpected: see table above"))
    print()
    print("   CONCLUSION. an observed heritable low-burden state after asymmetric")
    print("   segregation is evidence for a BOUNDED high-burden attractor, hence")
    print("   for asymptotic rather than complete growth arrest. that is evidence")
    print("   about the unmeasured quantity of D015, arriving from the aging")
    print("   literature rather than from a proteostasis experiment.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
