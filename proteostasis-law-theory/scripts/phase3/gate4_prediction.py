"""phase 3: the quantitative prediction Gate 4 is preregistered against.

Gate 4 (`empirical/GATE4_PROPOSAL.md`) tests H1' -- that at the viability
boundary the ClpXP arm is far from saturation -- against H0, capacity
exhaustion. To preregister it, the THRESHOLD must come from the model rather
than be chosen, so this module computes the predicted distribution of `s_u`, the
soluble-degradation saturation fraction at the collapse point.

Two populations are computed and they answer different questions:

  NON-DIVIDING   `s_u` at the fold across the phase 1 parameter box. this is the
                 distribution already reported as median 0.155 in Consequence 2.
  DIVIDING       `s_u` at the fold under the CALIBRATED growth-burden law, which
                 is the regime an actual experiment sits in. growth dilution is
                 a disposal channel (D010), so it changes where the boundary is
                 and therefore what saturation state is reached there.

The dividing population is the one Gate 4 must be judged against. Reporting the
non-dividing figure for a dividing-cell experiment would repeat the arm-quoting
error D019 warns about, in a different coordinate.

CLAIM LABELS
  Computational : every number here. it is a model prediction, not a measurement.
  Empirical     : nothing. no organism data is used.
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
import calibration as C  # noqa: E402
from fold_theorem import (KINETIC_FIELDS, foldStateFromSampleRow,  # noqa: E402
                          paramsFromSampleRow, phase1RunDir)

QUANTILES = [0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99]


def saturationAt(u: float, a: float, p: M.Params) -> Dict[str, float]:
    """the three saturation fractions at a state."""
    uf, af, cf, df = M.solveFreePools(u, a, p)
    return {"s_ref": uf / (p.kappa_ref + uf),
            "s_u": uf / (p.kappa_u + uf),
            "s_a": af / (p.kappa_a + af)}


def nonDividingPopulation(run: Optional[Path] = None) -> pd.DataFrame:
    """`s_u` at the recorded phase 1 folds, no dilution."""
    run = run or phase1RunDir()
    c = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    c = c[pd.to_numeric(c["fold_burden"], errors="coerce").notna()]
    rows: List[Dict[str, float]] = []
    for _, r in c.iterrows():
        try:
            p = paramsFromSampleRow(r)
            u, a = foldStateFromSampleRow(r)
            rows.append(saturationAt(u, a, p))
        except (M.ModelError, ValueError, KeyError):
            continue
    return pd.DataFrame(rows)


def dividingPopulation(run: Optional[Path] = None, k: int = 40,
                       p_qcs=(0.01, 0.02, 0.05), mu0s=(0.01, 0.05, 0.1),
                       seed: int = 17) -> pd.DataFrame:
    """`s_u` at the fold under the calibrated growth-burden law.

    kinetics are drawn from the phase 1 box; the growth law is the measured-slope
    linear form from `calibration.py`, swept over the two conversions that were
    never obtained.
    """
    run = run or phase1RunDir()
    c = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    draws = c.sample(n=k, random_state=seed)

    rows: List[Dict[str, float]] = []
    for idx, (_, r) in enumerate(draws.iterrows()):
        kin = {f: float(r["p_" + f]) for f in KINETIC_FIELDS}
        chi = float(r["p_chi"])
        for p_qc in p_qcs:
            for mu0 in mu0s:
                try:
                    p = M.Params(**kin).with_(nu=float(r["p_nu"]), c_tot=chi,
                                              d_tot=1.0 - chi).validate()
                    g = C.calibratedGrowth(mu0, p_qc)
                    out = D.foldSolveDil(p, g)
                except (M.ModelError, ValueError):
                    continue
                if out is None:
                    continue
                _, u_s, a_s = out
                try:
                    rec = saturationAt(u_s, a_s, p)
                except M.ModelError:
                    continue
                rec.update({"draw": idx, "p_qc": p_qc, "mu0": mu0})
                rows.append(rec)
    return pd.DataFrame(rows)


def thresholdFrom(df: pd.DataFrame, col: str = "s_u",
                  q: float = 0.99) -> float:
    """the preregistered rejection threshold for H1'.

    set at the model's own upper quantile: if the measured saturation fraction at
    the boundary exceeds the value that only 1 % of modelled networks reach, the
    far-from-saturation claim is rejected. this is derived from the prediction
    rather than chosen for convenience, and it is deliberately GENEROUS to H1 --
    a threshold at the median would reject on half the model's own draws.
    """
    return float(df[col].quantile(q))


def summarise(df: pd.DataFrame, name: str) -> Dict[str, float]:
    out = {"population": name, "n": len(df)}
    for c_ in ("s_ref", "s_u", "s_a"):
        if c_ in df:
            for q in QUANTILES:
                out["%s_p%g" % (c_, 100 * q)] = float(df[c_].quantile(q))
    return out


def main() -> int:
    run = phase1RunDir()
    if not run.is_dir():
        print("SKIP: no phase 1 run root (results/ is gitignored).")
        return 0

    nd = nonDividingPopulation(run)
    dv = dividingPopulation(run)

    print("PREDICTED SATURATION AT THE COLLAPSE BOUNDARY")
    print()
    for name, df in (("non-dividing (phase 1 box)", nd),
                     ("DIVIDING (calibrated growth law)", dv)):
        print("%s   n = %d" % (name, len(df)))
        tab = df[["s_ref", "s_u", "s_a"]].quantile(QUANTILES).T
        tab.columns = ["p1", "p5", "p25", "p50", "p75", "p95", "p99"]
        print(tab.to_string(float_format=lambda v: f"{v:.4f}"))
        print()

    t_nd = thresholdFrom(nd)
    t_dv = thresholdFrom(dv)
    print("K1 THRESHOLD for H1', derived as the p99 of the predicted s_u")
    print("   non-dividing p99 : %.4f" % t_nd)
    print("   DIVIDING     p99 : %.4f   <- the one Gate 4 uses" % t_dv)
    print()
    print("   H0 (capacity exhaustion) predicts s_u -> 1.")
    print("   separation between the dividing p99 and H0 : %.4f" % (1.0 - t_dv))
    print()
    frac = float((dv["s_u"] > 0.5).mean())
    print("   fraction of dividing draws with s_u > 0.5 : %.4f" % frac)
    print("   fraction of dividing draws with s_u > %.2f : %.4f"
          % (t_dv, float((dv["s_u"] > t_dv).mean())))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
