"""phase 3: a MEASURED growth-burden relation, and what it does to the results.

PROVENANCE -- the one measured anchor
-------------------------------------
Geiler-Samerotte KA, Dion MF, Budnik BA, Wang SM, Hartl DL, Drummond DA (2011)
"Misfolded proteins impose a dosage-dependent fitness cost and trigger a
cytosolic unfolded protein response in yeast." PNAS 108(2):680-685.
PMID 21187411, PMC3021021, doi:10.1073/pnas.1017570108. Retrieved via PubMed.

The reported quantity, verbatim from the abstract: "a fitness cost that increases
with misfolded protein abundance, up to as much as a 3.2% growth rate reduction
when misfolded YFP represents less than 0.1% of total cellular protein."

That is a real, dosage-resolved measurement of growth rate against MISFOLDED
protein burden -- the exact relation `dilution.Growth` previously guessed. Two
things must be said about it immediately:

  * it is YEAST, not E. coli. no equivalent dosage-resolved bacterial measurement
    was located.
  * "less than 0.1%" is an UPPER bound on the abundance, so the slope derived
    from it is a LOWER bound on steepness and the arrest burden an UPPER bound.

WHAT IS STILL NOT MEASURED
--------------------------
Using the anchor inside the model needs two conversions, and NEITHER was obtained
as a number in this work. They are therefore SWEPT, not assumed:

  P_qc   the chaperone + protease share of total cellular protein. the model's
         concentration scale is phi = C_tot + D_tot, so a model burden of 1
         equals a proteome fraction of P_qc. the canonical source for this is
         Schmidt et al. (2016) Nat Biotechnol 34:104-110, PMID 26641532,
         doi:10.1038/nbt.3418 -- the value was NOT extracted from it here.
  mu0    the dilution rate in model time units, i.e. growth rate divided by the
         refolding turnover k_ref. no k_ref was obtained.

So this module does NOT deliver "the calibrated model". It delivers the honest
version: the measured slope, propagated through an explicitly swept conversion,
answering ONE question -- do the phase 3 conclusions survive calibration, or do
they depend on where in that uncertainty the truth lies?

THE DERIVED LAW
---------------
Relative growth reduction is linear in misfolded proteome fraction with slope

    S = 0.032 / 0.001 = 32 per unit proteome fraction   (a lower bound)

so growth would arrest at a misfolded fraction of 1/S = 0.03125. Converting,

    mu(u,a) = mu0 . max(0, 1 - (u+a)/k_mu),      k_mu = 1 / (S . P_qc)

which is the LINEAR form already implemented as `boundary_structure.LinearGrowth`
-- the shape is what proteome partitioning predicts and what this measurement
supports. Only `k_mu` is new, and it is now derived rather than chosen.

CLAIM LABELS
  Mathematical  : the conversion algebra.
  Computational : every simulation result below.
  Empirical     : the slope S, which is measured -- IN YEAST, at one dose range.
                  no claim is made that it transfers to E. coli.
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

# --- the measured anchor, as reported -------------------------------------
GS2011 = {
    "doi": "10.1073/pnas.1017570108",
    "pmid": "21187411",
    "growth_rate_reduction": 0.032,        # 3.2 %
    "misfolded_proteome_fraction": 0.001,  # "less than 0.1 %" -> upper bound
    "organism": "Saccharomyces cerevisiae",
}

#: slope of relative growth reduction per unit misfolded proteome fraction.
#: a LOWER bound, because the abundance figure it divides by is an upper bound.
SLOPE_PER_PROTEOME_FRACTION = (GS2011["growth_rate_reduction"]
                               / GS2011["misfolded_proteome_fraction"])

#: misfolded proteome fraction at which the linear law reaches zero growth.
ARREST_PROTEOME_FRACTION = 1.0 / SLOPE_PER_PROTEOME_FRACTION


def kMuFromProteomeShare(p_qc: float) -> float:
    """burden (in model units) at which the measured law arrests growth.

    a model burden of 1 equals a proteome fraction of `p_qc`, since the model's
    concentration scale is the total rescue pool.
    """
    if not (0.0 < p_qc < 1.0):
        raise ValueError("p_qc must be a proteome fraction in (0,1)")
    return ARREST_PROTEOME_FRACTION / p_qc


def calibratedGrowth(mu0: float, p_qc: float) -> B.LinearGrowth:
    """the measured law, converted into model units."""
    return B.LinearGrowth(mu0=float(mu0), k_mu=kMuFromProteomeShare(p_qc))


def sweepCalibration(p_qcs=(0.005, 0.01, 0.02, 0.03, 0.05, 0.10),
                     mu0s=(0.001, 0.01, 0.05, 0.1, 0.3),
                     p: Optional[M.Params] = None) -> pd.DataFrame:
    """does a collapse boundary survive under the MEASURED law?

    swept over the two unmeasured conversions. the question is whether the
    qualitative phase 3 conclusions depend on where in that range the truth sits.
    """
    p = (p or M.Params()).validate()
    rows: List[Dict[str, object]] = []
    for p_qc in p_qcs:
        k_mu = kMuFromProteomeShare(p_qc)
        for mu0 in mu0s:
            g = calibratedGrowth(mu0, p_qc)
            dec = B.boundaryDecomposition(p, g)
            rows.append({
                "p_qc": p_qc, "k_mu": k_mu, "mu0": mu0,
                "boundary": dec is not None,
                "j_crit": None if dec is None else dec["j_crit"],
                "phi_enz": None if dec is None else dec["phi_enz"],
                "delta": None if dec is None else dec["delta"],
            })
    return pd.DataFrame(rows)


def compareToGuess(p: Optional[M.Params] = None,
                   p_qcs=(0.005, 0.01, 0.02, 0.05, 0.10)) -> pd.DataFrame:
    """calibrated k_mu against the values guessed before any data was consulted.

    the guesses were k_mu = 0.5 (hyperbolic, `dilution.py`) and 2.0 (linear,
    `boundary_structure.py`). this says whether they were in the right regime.
    """
    rows = []
    for p_qc in p_qcs:
        rows.append({"p_qc": p_qc, "k_mu_calibrated": kMuFromProteomeShare(p_qc),
                     "k_mu_guess_hyperbolic": 0.5, "k_mu_guess_linear": 2.0})
    return pd.DataFrame(rows)


def bistabilityUnderCalibration(p_qc: float = 0.02, mu0: float = 0.05,
                                p: Optional[M.Params] = None,
                                n_j: int = 8) -> Dict[str, object]:
    """is the system still bistable when the growth law is the measured one?

    the answer turns on the HIGH-burden shape of the growth law, which the
    measurement does not constrain. the linear form arrests growth completely
    beyond `k_mu`, which switches dilution off and leaves the high-burden state
    unbounded; the hyperbolic form only approaches zero, so dilution keeps
    bounding it. hence the settled-state check is the point of this function --
    an unsettled upper branch is a runaway, not an attractor.
    """
    p = (p or M.Params()).validate()
    g = calibratedGrowth(mu0, p_qc)
    dec = B.boundaryDecomposition(p, g)
    if dec is None:
        return {"boundary": False}
    j = dec["j_crit"]
    js = [round(x, 6) for x in np.geomspace(0.1 * j, 1.25 * j, n_j)]
    h = B.hysteresisSweep(p, g, js)
    return {"boundary": True, "j_crit": j, "k_mu": g.k_mu, "mu0": mu0,
            "all_settled": h["all_settled"],
            "down_settled": sum(1 for x in js if h["settled_down"][x]),
            "n_j": len(js),
            "window": h["window"], "hysteretic": bool(h["hysteretic_j"]),
            "top_state": h["up"][js[-1]]}


def main() -> int:
    print("MEASURED ANCHOR  Geiler-Samerotte et al. 2011 PNAS,")
    print("                 doi:%s (PMID %s), via PubMed" % (GS2011["doi"], GS2011["pmid"]))
    print("   %.1f%% growth-rate reduction at <%.1f%% of total cellular protein misfolded (%s)"
          % (100 * GS2011["growth_rate_reduction"],
             100 * GS2011["misfolded_proteome_fraction"], GS2011["organism"]))
    print("   -> slope %.1f per unit proteome fraction (a LOWER bound)"
          % SLOPE_PER_PROTEOME_FRACTION)
    print("   -> linear arrest at a misfolded proteome fraction of %.4f (an UPPER bound)"
          % ARREST_PROTEOME_FRACTION)
    print()

    print("CALIBRATED k_mu vs the values guessed before any data was consulted")
    print(compareToGuess().to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print()

    print("DOES A BOUNDARY SURVIVE UNDER THE MEASURED LAW?")
    S = sweepCalibration()
    piv = S.pivot_table(index="p_qc", columns="mu0", values="boundary",
                        aggfunc="first")
    print(piv.to_string())
    n_ok = int(S["boundary"].sum())
    print("   boundary present in %d of %d calibrated cells" % (n_ok, len(S)))
    ok = S[S["boundary"]]
    print("   j_crit  range %.4f - %.4f" % (ok["j_crit"].min(), ok["j_crit"].max()))
    print("   phi_enz range %.4f - %.4f  (the enzymatic condition)"
          % (ok["phi_enz"].min(), ok["phi_enz"].max()))
    print("   delta   range %.4f - %.4f  (dilution share of disposal)"
          % (ok["delta"].min(), ok["delta"].max()))
    print()

    print("IS IT STILL BISTABLE UNDER THE MEASURED LAW?")
    print("   (an unsettled DOWN branch means a runaway, not a second attractor)")
    for p_qc, mu0 in ((0.02, 0.05), (0.02, 0.1), (0.05, 0.1)):
        r = bistabilityUnderCalibration(p_qc, mu0)
        if not r["boundary"]:
            print("   p_qc=%.3f mu0=%.3g : no boundary" % (p_qc, mu0))
            continue
        print("   p_qc=%.3f mu0=%.3g k_mu=%.3f : j_crit=%.5f  down-branch settled %d/%d  window=%s"
              % (p_qc, mu0, r["k_mu"], r["j_crit"],
                 r["down_settled"], r["n_j"], r["window"]))
    print()
    print("   CONTRAST -- the hyperbolic form at the same dilution rate:")
    ph = M.Params().validate()
    h = B.hysteresisSweep(ph, D.Growth(0.04), [0.10, 0.155, 0.17, 0.19, 0.196, 0.21])
    print("      all endpoints settled %s, window %s, top state u=%.4f a=%.4f"
          % (h["all_settled"], h["window"], *h["up"][0.21]))
    print("   so bistability is a property of the growth law's HIGH-burden shape,")
    print("   which the measurement does not constrain: it only fixes the slope")
    print("   at <0.1%% misfolded, three decades below the arrest burden.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
