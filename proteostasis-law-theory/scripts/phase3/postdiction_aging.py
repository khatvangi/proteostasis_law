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

The model produces that structure under SOME dilution laws (D013), so the
expectation going in was that the aging literature would bear on D015's open
question -- whether growth arrest is asymptotic or complete -- from a direction
no proteostasis experiment could reach.

It does bear on it, but not the way expected. See the findings below: the
expectation was wrong, and the regime that supplies bistability is the one that
gets the measured quantity wrong.

WHAT THIS MODULE FOUND -- A NEGATIVE, AND A POINTER
---------------------------------------------------
The logical point above stands on its own: rejuvenation is only a coherent
category in a bistable system. But the model does NOT robustly supply that
bistability once growth responds to burden.

  * under CONSTANT dilution (`k_mu = inf`) the system is bistable and a daughter
    must shed 35-82% of its aggregate to escape the high state. but constant
    dilution means growth does not respond to burden at all, so it predicts ZERO
    reproductive loss -- flatly contradicting the >30% that is the observation.
  * under HYPERBOLIC FEEDBACK, which is the physiologically appropriate law, four
    of six settings tested are MONOSTABLE. where bistability does appear it
    carries a predicted reproductive loss of 48-95% (median 51%), more severe
    than the reported figure.
  * under LINEAR arrest there is no bounded high state at all, so no bistability
    and no rejuvenation.

So the model does not post-dict aging and rejuvenation. Reporting it as a success
would have required quoting the constant-dilution regime, which is exactly the
regime that predicts the wrong answer for the quantity actually measured.

WHAT THE FAILURE POINTS AT
--------------------------
The Lindner observation is not merely that a high-burden state exists; it is that
aggregates are SPATIALLY SEQUESTERED at a pole and segregated asymmetrically at
division. This model is well-mixed. Sequestration into an inclusion body removes
aggregate from the reactive pool -- it changes the kinetics, not just the
bookkeeping -- and that is a mechanism capable of producing a stable high-burden
state without requiring the growth law to do the work.

The honest reading is therefore that the aging literature identifies a MISSING
MECHANISM rather than confirming the present one, and that spatial sequestration
is the specific candidate.

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


def feedbackScan(p: Optional[M.Params] = None) -> pd.DataFrame:
    """is the system bistable once growth actually responds to burden?

    constant dilution is not a physiological law: it says growth is unaffected by
    damage, which contradicts the very observation being explained. this scans
    hyperbolic FEEDBACK laws instead and records the predicted reproductive loss
    of the high state, which is the quantity Lindner et al. measured.
    """
    p = (p or M.Params()).validate()
    rows: List[Dict[str, float]] = []
    for mu0, k_mu in ((0.05, 0.5), (0.1, 0.5), (0.2, 0.5),
                      (0.1, 1.0), (0.3, 1.0), (0.2, 2.0)):
        g = D.Growth(mu0, k_mu)
        out = D.foldSolveDil(p, g)
        if out is None:
            continue
        j_crit = out[0]
        for frac in (0.90, 0.95, 0.98, 0.995):
            j = frac * j_crit
            low = B.settle(j, (0.0, 0.0), p, g)
            high = B.settle(j, (low[0], max(50 * max(low[1], 1e-6), 5.0)), p, g)
            if not all(np.isfinite(v) for v in low + high):
                continue
            pj = p.with_(j=j).validate()
            if any(max(map(abs, D.rhsDil(s[0], s[1], pj, g)))
                   > 1e-7 * max(1.0, s[0] + s[1]) for s in (low, high)):
                continue
            bistable = high[1] > 1.5 * max(low[1], 1e-12)
            m_lo, m_hi = g.rate(*low) / mu0, g.rate(*high) / mu0
            rows.append({"mu0": mu0, "k_mu": k_mu, "j_over_jcrit": frac,
                         "bistable": bistable,
                         "a_low": low[1], "a_high": high[1],
                         "reproductive_loss": 1.0 - m_hi / m_lo if bistable else 0.0})
    return pd.DataFrame(rows)


def main() -> int:
    p = M.Params().validate()
    print("POST-DICTION ATTEMPT: aging and rejuvenation by asymmetric segregation")
    print("   observation: Lindner et al. 2008 PNAS, doi:%s (PMID %s), via PubMed"
          % (LINDNER2008["doi"], LINDNER2008["pmid"]))
    print("   aggregates accumulate in old-pole cells under %s, costing >%d%% of"
          % (LINDNER2008["condition"],
             100 * LINDNER2008["reproductive_loss_old_pole"]))
    print("   reproductive ability; new-pole progeny 'exhibit rejuvenation'.")
    print()
    print("   THE LOGICAL POINT, which stands independently: rejuvenation is only a")
    print("   coherent category in a BISTABLE system. in a monostable one a daughter")
    print("   with less aggregate simply relaxes back to the single attractor.")
    print()

    print("A. CONSTANT dilution -- bistable, but growth cannot respond to burden")
    df = scanInfluxes(p, D.Growth(0.04))
    if not df.empty:
        print("   aggregate ratio high/low    : %.1f - %.1f fold"
              % (df["aggregate_ratio"].min(), df["aggregate_ratio"].max()))
        print("   shed fraction to rejuvenate : %.3f - %.3f (median %.3f)"
              % (df["shed_fraction_needed"].min(), df["shed_fraction_needed"].max(),
                 df["shed_fraction_needed"].median()))
    print("   predicted reproductive loss : 0.000  <- k_mu = inf means growth is")
    print("                                 unaffected by burden BY CONSTRUCTION,")
    print("                                 contradicting the >30%% that was measured.")
    print()

    print("B. HYPERBOLIC FEEDBACK -- the physiologically appropriate law")
    fb = feedbackScan(p)
    if fb.empty:
        print("   no settings converged")
        return 0
    n_bi = int(fb["bistable"].sum())
    laws = fb.groupby(["mu0", "k_mu"])["bistable"].any()
    print("   settings tested          : %d  (bistable at some influx in %d)"
          % (len(laws), int(laws.sum())))
    print("   influx points evaluated  : %d  (bistable in %d)" % (len(fb), n_bi))
    if n_bi:
        bi = fb[fb["bistable"]]
        print("   predicted reproductive loss of the high state : %.3f - %.3f (median %.3f)"
              % (bi["reproductive_loss"].min(), bi["reproductive_loss"].max(),
                 bi["reproductive_loss"].median()))
        print("   observed (Lindner)                            : > %.2f"
              % LINDNER2008["reproductive_loss_old_pole"])
    print()

    print("C. LINEAR arrest -- no bounded high state at all")
    import boundary_structure as BS
    dfl = scanInfluxes(p, BS.LinearGrowth(0.3, 1.5625))
    print("   bistable influxes found : %d" % len(dfl))
    print()

    print("VERDICT: the model does NOT post-dict aging and rejuvenation.")
    print("   bistability is generic only under constant dilution, which predicts")
    print("   zero reproductive loss; under feedback it is rare and the loss is")
    print("   more severe than observed; under linear arrest it does not occur.")
    print()
    print("   WHAT THE FAILURE POINTS AT: the observation is about aggregates being")
    print("   SPATIALLY SEQUESTERED at a pole and segregated asymmetrically. this")
    print("   model is well-mixed. sequestration removes aggregate from the reactive")
    print("   pool, changing the kinetics rather than the bookkeeping, and is a")
    print("   candidate mechanism for a stable high-burden state that does not make")
    print("   the growth law do the work. that is a missing mechanism, identified")
    print("   by a failed post-diction rather than by a confirmed one.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
