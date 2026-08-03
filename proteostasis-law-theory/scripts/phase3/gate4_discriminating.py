"""phase 3: a prediction that is robust AND discriminating, not merely robust.

THE PROBLEM THIS SOLVES
-----------------------
Gate 4's primary outcome (D021) is the critical-slowing exponent, tau ~
(j_crit - j)^(-1/2). It is robust because 1/2 is generic to saddle-nodes -- and
that is exactly why it does not discriminate. Confirming it supports "the
boundary is a saddle-node" and selects this model over nothing at all.

Meanwhile the saturation fraction (D020) discriminates in principle but its
predicted spread covers nearly [0,1], so it cannot be tested.

A useful test needs both properties. This module quantifies one.

THE PREDICTION
--------------
The model's distinctive structural claim is that rescue capacity is CONSERVED AND
SHARED: free chaperone is what remains after ordinary nascent chains, damaged
monomer and aggregate have all taken their share. The nascent load `nu` therefore
consumes capacity while contributing NO damage influx of its own --

    cf = c_tot / (1 + nu + uf/kappa_cu + af/kappa_ca)

so `nu` enters only through the denominator. The consequence:

    RAISING THE LOAD OF PERFECTLY-FOLDING PROTEIN MUST LOWER THE TOLERABLE
    MISTRANSLATION DOSE, EVEN THOUGH THAT PROTEIN CAUSES NO DAMAGE.

Phase 1 recorded this as C3 in 97.11 % of draws, but only as a yes/no. This
module measures the effect SIZE, which is what decides whether it is executable.

WHY IT DISCRIMINATES
--------------------
- **Shared-capacity models** (this one) predict j_crit falls as nu rises.
- **Independent-handling models**, in which misfolded-protein disposal draws on a
  pool not shared with ordinary folding, predict NO shift.
- **Pure capacity-exhaustion models** predict a shift only if the exhausted
  capacity is the shared one -- which is itself the claim under test, so the
  contrast is informative either way.

Crucially the perturbation is orthogonal to the outcome: gratuitous expression of
a well-folding protein adds folding load without adding misfolding load, so a
positive result cannot be explained by the perturbation itself causing damage.

CLAIM LABELS
  Computational : every number here.
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
from fold_theorem import KINETIC_FIELDS, phase1RunDir  # noqa: E402

#: nascent-load ladder, a ~100x span. chosen by sweeping candidate ladders:
#: 30x gives a 1.12x median shift, 100x gives 1.24x with direction still correct
#: in 100% of converged networks, and wider ladders (400x, 5000x) buy more effect
#: at the cost of direction consistency (95%, 94%) as some networks lose the fold.
NU_LADDER = (0.05, 0.3, 1.5, 5.0)

#: THE CONFOUND THIS DESIGN MUST DEFEAT.
#: gratuitous protein expression is the standard way to RAISE nu, and it is also
#: the standard way to LOWER growth rate (that is what Scott et al. 2010 used it
#: for). growth rate is part of disposal (D010), so in batch culture the
#: perturbation moves both the thing being varied and the thing being measured,
#: and the result is uninterpretable in either direction.
#: the design is therefore only valid at EXTERNALLY FIXED growth rate -- a
#: chemostat or turbidostat, where dilution rate is set by the operator and
#: gratuitous expression changes nu without changing mu. this is not a
#: refinement; a batch version of this experiment cannot test the prediction.
REQUIRES_FIXED_GROWTH_RATE = True


def criticalInfluxVsNascent(p: M.Params, g: Optional[D.Growth] = None,
                            nus=NU_LADDER) -> Optional[Dict[str, float]]:
    """j_crit across a nascent-load ladder at fixed kinetics.

    returns the ladder, the fold-change from lowest to highest nu, and the sign
    consistency of the shift.
    """
    g = g or D.Growth(0.0)
    js: List[float] = []
    for nu in nus:
        try:
            pn = p.with_(nu=float(nu)).validate()
            out = D.foldSolveDil(pn, g)
        except (M.ModelError, ValueError, OverflowError):
            return None
        if out is None or not np.isfinite(out[0]) or out[0] <= 0.0:
            return None
        js.append(float(out[0]))
    arr = np.array(js)
    monotone_down = bool(np.all(np.diff(arr) < 0))
    return {"j_lo": arr[0], "j_hi": arr[-1],
            "fold_drop": float(arr[0] / arr[-1]),
            "monotone_down": monotone_down,
            "any_drop": bool(arr[-1] < arr[0]),
            "js": js}


def sweep(k: int = 30, seed: int = 31) -> pd.DataFrame:
    """the effect across networks, without dilution and under the measured law."""
    run = phase1RunDir()
    c = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    draws = c.sample(n=k, random_state=seed)

    rows: List[Dict[str, float]] = []
    for idx, (_, r) in enumerate(draws.iterrows()):
        kin = {f: float(r["p_" + f]) for f in KINETIC_FIELDS}
        chi = float(r["p_chi"])
        try:
            p = M.Params(**kin).with_(c_tot=chi, d_tot=1.0 - chi).validate()
        except M.ModelError:
            continue
        for label, g in (("no dilution", D.Growth(0.0)),
                         ("calibrated mu0=0.05", C.calibratedGrowth(0.05, 0.02)),
                         ("calibrated mu0=0.10", C.calibratedGrowth(0.10, 0.02))):
            res = criticalInfluxVsNascent(p, g)
            if res is None:
                continue
            res = {kk: vv for kk, vv in res.items() if kk != "js"}
            res.update({"draw": idx, "regime": label})
            rows.append(res)
    return pd.DataFrame(rows)


def main() -> int:
    if not phase1RunDir().is_dir():
        print("SKIP: no phase 1 run root (results/ is gitignored).")
        return 0

    print("DISCRIMINATING PREDICTION: does harmless folding load shrink the")
    print("tolerable mistranslation dose?  nu ladder %s (%.0fx span)"
          % (NU_LADDER, NU_LADDER[-1] / NU_LADDER[0]))
    print()
    print("   shared-capacity model (this one) : j_crit FALLS as nu rises")
    print("   independent-handling model       : no shift")
    print()
    print("   REQUIRES externally fixed growth rate: gratuitous expression raises")
    print("   nu and lowers mu, and mu is part of disposal. batch culture confounds")
    print("   the perturbation with the readout and cannot test this.")
    print()

    df = sweep()
    if df.empty:
        print("   no networks converged")
        return 0

    for regime, sub in df.groupby("regime"):
        print("%-22s n=%-4d  drop in %d/%d  monotone in %d/%d  fold-drop median %.3f "
              "(IQR %.3f-%.3f)"
              % (regime, len(sub), int(sub["any_drop"].sum()), len(sub),
                 int(sub["monotone_down"].sum()), len(sub),
                 sub["fold_drop"].median(),
                 sub["fold_drop"].quantile(.25), sub["fold_drop"].quantile(.75)))
    print()
    print("   OVERALL  direction correct in %d of %d (%.1f%%), monotone in %.1f%%"
          % (int(df["any_drop"].sum()), len(df),
             100 * df["any_drop"].mean(), 100 * df["monotone_down"].mean()))
    print("   effect size: j_crit falls %.2fx median (p10 %.2f, p90 %.2f) over the ladder"
          % (df["fold_drop"].median(), df["fold_drop"].quantile(.1),
             df["fold_drop"].quantile(.9)))
    big = float((df["fold_drop"] > 1.5).mean())
    print("   networks with a >1.5x shift (comfortably measurable) : %.1f%%" % (100 * big))
    print()
    print("   READ WITH D021: the exponent work supplies the INSTRUMENT -- tau^-2")
    print("   regressed on dose gives j_crit as an x-intercept with a CI. this")
    print("   prediction is then whether that intercept MOVES. the exponent is the")
    print("   ruler; the nascent-load shift is the test.")
    print()
    ok = df["any_drop"].mean() > 0.95 and df["fold_drop"].median() > 1.2
    print("   -> %s" % ("EXECUTABLE at fixed growth rate: direction is near-universal "
                        "and the shift is resolvable by a regression intercept"
                        if ok else
                        "MARGINAL: direction is robust but the shift is small; power "
                        "depends on the intercept CI"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
