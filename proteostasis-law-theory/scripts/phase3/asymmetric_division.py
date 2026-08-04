"""phase 3, the framing test: was bistability ever required?

THE CLAIM UNDER TEST
--------------------
Three post-dictions failed in the same direction, all because the model's bistable
high state is 7.5x to 254x more aggregate-laden than the cell Lindner measured
(D029). Before adding a fourth mechanism, test whether the requirement was real.

D026's surviving claim was: "rejuvenation is only a coherent category in a
BISTABLE system, since in a monostable one a daughter with less aggregate relaxes
straight back to the single attractor." That argument assumes the old-pole cell is
SITTING in a second basin. It is not. It inherits a physical inclusion body at
every division. That is a CONTINUOUSLY RENEWED PERTURBATION, not an attractor, and
a monostable system driven by a renewed perturbation has a stationary offset with
no separatrix anywhere.

THE MECHANISM
-------------
No sequestration, no second attractor, no new state variable. The two-state
diluted model under the calibrated hyperbolic law, plus one thing: division
partitions aggregate ASYMMETRICALLY.

Between divisions the state follows `dilution.rhsDil`, whose `-mu.x` terms are
dilution by volume growth. Division itself halves the volume. For a well-mixed
species split evenly, concentration is unchanged, which is exactly what continuous
dilution already encodes -- so a SYMMETRIC split adds nothing, and `f = 0.5` is the
exact control. For an asymmetric split the daughter receiving fraction `f` of the
aggregate in half the volume starts at concentration `2f . a`:

    old-pole daughter :  a -> 2f.a          f > 0.5
    new-pole daughter :  a -> 2(1-f).a
    soluble monomer   :  u -> u             well mixed, splits evenly

`f = 0.5` gives `a -> a` for both, i.e. the plain diluted model.

WHAT IS MEASURED, AND WHY IT IS THE SAME QUANTITY
-------------------------------------------------
Division is defined by volume doubling: the generation ends when the accumulated
`integral of mu dt` reaches `ln 2`. The interdivision time `T` is therefore an
OUTPUT, and the growth rate `ln2/T` is the observable Lindner reports rather than
a parameter. The old-pole lineage is iterated to its stationary cycle; at the
division out of that cycle the two daughters are followed for one generation each
and scored as

    Delta(GR_old - GR_new) / GR_mean

against the data-derived band from `postdiction_aging.aggregateAttributableLoss()`.
The band is FIXED and is not re-derived here.

CLAIM LABELS
  Computational : every number below.
  Empirical     : the Lindner comparison, a post-diction against a published
                  measurement.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import dilution as D  # noqa: E402
import calibration as C  # noqa: E402
from postdiction_aging import aggregateAttributableLoss  # noqa: E402

LN2 = float(np.log(2.0))

# the band is fixed (D028) and imported, never recomputed here
BAND_LO, BAND_HI = aggregateAttributableLoss()

# f = 0.5 is the mechanism-OFF control and must not carry a pass (protocol rule 6)
CONTROL_F = 0.5


def calibratedHyperbolic(mu0: float, p_qc: float) -> D.Growth:
    """the D015 conversion, applied to the HYPERBOLIC form rather than the linear.

    `calibration.calibratedGrowth` returns the linear-arrest law. The framing test
    needs the hyperbolic one, which never sets mu to exactly zero, so a
    reproductive loss of exactly 1.000 cannot be produced by the law itself.
    """
    return D.Growth(mu0=float(mu0), k_mu=C.kMuFromProteomeShare(p_qc))


def oneGeneration(x0: Tuple[float, float], p: M.Params, g: D.Growth,
                  t_max_mult: float = 200.0) -> Optional[Dict[str, float]]:
    """integrate until the cell has doubled in volume; return the end state and T.

    the third variable accumulates `integral of mu dt`; the generation ends when it
    reaches ln 2. so the interdivision time is an OUTPUT of the burden, which is
    what makes `ln2/T` comparable to a measured growth rate.
    """
    if g.mu0 <= 0.0:
        return None
    t_max = t_max_mult * LN2 / g.mu0

    def field(_t, y):
        u, a = max(y[0], 0.0), max(y[1], 0.0)
        try:
            du, da = D.rhsDil(u, a, p, g)
        except (M.ModelError, OverflowError):
            return [0.0, 0.0, 0.0]
        return [du, da, g.rate(u, a)]

    def doubled(_t, y):
        return y[2] - LN2
    doubled.terminal = True
    doubled.direction = 1.0

    s = solve_ivp(field, [0.0, t_max], [x0[0], x0[1], 0.0], method="Radau",
                  events=doubled, rtol=1e-10, atol=1e-13)
    if not s.t_events[0].size:
        return None                       # never doubled: growth arrested
    T = float(s.t_events[0][0])
    if not (np.isfinite(T) and T > 0.0):
        return None
    return {"u": float(s.y[0, -1]), "a": float(s.y[1, -1]), "T": T,
            "growth_rate": LN2 / T}


def oldPoleCycle(p: M.Params, g: D.Growth, f: float, n_gen: int = 400,
                 tol: float = 1e-10) -> Optional[Dict[str, float]]:
    """iterate the old-pole lineage to its stationary division cycle.

    the map is: integrate one generation, then keep fraction `f` of the aggregate
    in half the volume. a fixed point of that map is the stationary old-pole cell.
    """
    x = (1e-4, 1e-6)
    prev = None
    for k in range(n_gen):
        out = oneGeneration(x, p, g)
        if out is None:
            return None
        x_end = (out["u"], out["a"])
        x = (x_end[0], 2.0 * f * x_end[1])          # old-pole daughter
        if prev is not None:
            scale = max(1.0, abs(x[0]) + abs(x[1]))
            if max(abs(x[0] - prev[0]), abs(x[1] - prev[1])) / scale < tol:
                return {"u_start": x[0], "a_start": x[1],
                        "u_end": x_end[0], "a_end": x_end[1],
                        "T": out["T"], "generations": k + 1, "converged": True}
        prev = x
    return {"u_start": x[0], "a_start": x[1], "u_end": np.nan, "a_end": np.nan,
            "T": np.nan, "generations": n_gen, "converged": False}


def lineageContrast(p: M.Params, g: D.Growth, f: float) -> Optional[Dict[str, float]]:
    """the measured quantity: Delta(GR_old - GR_new)/GR_mean at one division.

    both daughters come from the SAME mother, the stationary old-pole cell, which
    is how Lindner et al. pair them.
    """
    cyc = oldPoleCycle(p, g, f)
    if cyc is None or not cyc["converged"]:
        return None
    ue, ae = cyc["u_end"], cyc["a_end"]
    old = oneGeneration((ue, 2.0 * f * ae), p, g)
    new = oneGeneration((ue, 2.0 * (1.0 - f) * ae), p, g)
    if old is None or new is None:
        return None
    gr_old, gr_new = old["growth_rate"], new["growth_rate"]
    gr_mean = 0.5 * (gr_old + gr_new)
    if gr_mean <= 0.0:
        return None
    return {"f": f, "a_mother_end": ae, "a_old_start": 2.0 * f * ae,
            "a_new_start": 2.0 * (1.0 - f) * ae,
            "gr_old": gr_old, "gr_new": gr_new, "gr_mean": gr_mean,
            "aging_effect": (gr_new - gr_old) / gr_mean,
            "generations_to_cycle": cyc["generations"]}


def inBand(x: float) -> bool:
    return bool(BAND_LO <= x <= BAND_HI)


# ---------------------------------------------------------------------------
# the sweep
# ---------------------------------------------------------------------------

F_LADDER = (0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.99)
PQC_LADDER = (0.005, 0.01, 0.02, 0.03, 0.05, 0.10)
MU0_LADDER = (0.01, 0.05, 0.1, 0.3)
J_FRACS = (0.2, 0.5, 0.8)


def sweep(p: Optional[M.Params] = None, f_ladder=F_LADDER,
          pqc_ladder=PQC_LADDER, mu0_ladder=MU0_LADDER,
          j_fracs=J_FRACS) -> pd.DataFrame:
    """every (f, p_qc, mu0, j) cell. nothing dropped, per protocol rule 3."""
    p = (p or M.Params()).validate()
    rows: List[Dict[str, object]] = []
    for p_qc in pqc_ladder:
        for mu0 in mu0_ladder:
            g = calibratedHyperbolic(mu0, p_qc)
            fold = D.foldSolveDil(p, g)
            if fold is None:
                rows.append({"p_qc": p_qc, "mu0": mu0, "fold": False,
                             "j_frac": np.nan, "f": np.nan,
                             "aging_effect": np.nan, "in_band": False})
                continue
            for jf in j_fracs:
                pj = p.with_(j=jf * fold[0]).validate()
                for f in f_ladder:
                    try:
                        r = lineageContrast(pj, g, f)
                    except (M.ModelError, ValueError, OverflowError):
                        r = None
                    base = {"p_qc": p_qc, "mu0": mu0, "fold": True,
                            "j_frac": jf, "f": f}
                    if r is None:
                        rows.append({**base, "aging_effect": np.nan,
                                     "in_band": False, "settled": False})
                        continue
                    rows.append({**base, "settled": True,
                                 "aging_effect": r["aging_effect"],
                                 "gr_old": r["gr_old"], "gr_new": r["gr_new"],
                                 "a_old_start": r["a_old_start"],
                                 "a_new_start": r["a_new_start"],
                                 "in_band": inBand(r["aging_effect"])})
    return pd.DataFrame(rows)


def verdict(df: pd.DataFrame) -> Dict[str, object]:
    """rule 6: the control f = 0.5 must not carry the pass."""
    if "aging_effect" in df:
        df = df.assign(in_band=[bool(np.isfinite(x)) and inBand(float(x))
                                for x in df["aging_effect"].fillna(np.nan)])
    on = df[df["f"] > CONTROL_F]
    ctrl = df[df["f"] == CONTROL_F]
    hit = on[on["in_band"]]
    return {
        "n_cells": int(len(df)),
        "settled": int(df["settled"].sum()) if "settled" in df else 0,
        "in_band_mechanism_on": int(on["in_band"].sum()),
        "in_band_control": int(ctrl["in_band"].sum()),
        "mechanism_passes": bool(on["in_band"].any()),
        "f_in_band": (float(hit["f"].min()), float(hit["f"].max()))
                     if len(hit) else None,
        "band": (BAND_LO, BAND_HI),
    }


def main() -> int:
    pd.set_option("display.width", 150)
    print("FRAMING TEST: does the observation require bistability at all?")
    print("  band (fixed, D028): [%.5f, %.5f]" % (BAND_LO, BAND_HI))
    print("  mechanism: asymmetric partitioning at division, f > 0.5")
    print("  control  : f = 0.5, which is the plain diluted model exactly")
    df = sweep()
    v = verdict(df)
    for k in ("n_cells", "settled", "in_band_mechanism_on", "in_band_control",
              "mechanism_passes", "f_in_band"):
        print("  %-22s %s" % (k, v[k]))

    s = df[df.get("settled", False) == True]  # noqa: E712
    if not s.empty:
        print("\n  aging effect by f (median over all settings):")
        print(s.groupby("f")["aging_effect"]
               .agg(["count", "min", "median", "max"]).to_string())
        print("\n  control f = 0.5 must be ~0:")
        print(s[s["f"] == CONTROL_F]["aging_effect"].describe().to_string())
    print("\n  VERDICT:", "IN BAND" if v["mechanism_passes"] else "NOT IN BAND")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


# ---------------------------------------------------------------------------
# the parameter-free form: what the model actually requires
# ---------------------------------------------------------------------------

# Tomoyasu T, Mogk A, Langen H, Goloubinoff P, Bukau B (2001) Mol Microbiol
# 40(2):397-413, doi:10.1046/j.1365-2958.2001.02383.x, PMID 11309122, via PubMed.
# "In DeltarpoH mutants ... 5-10% and 20-30% of total protein aggregated at 30
#  degrees C and 42 degrees C respectively"; "In rpoH+ cells, DnaK depletion did
#  not lead to protein aggregation at 30 degrees C".
# So the WILD-TYPE aggregate fraction is an upper bound (undetected), not a value.
TOMOYASU2001 = {
    "doi": "10.1046/j.1365-2958.2001.02383.x",
    "pmid": "11309122",
    "rpoH_null_30C": (0.05, 0.10),
    "wild_type_30C": None,          # undetected: a bound, not a measurement
}

# (!) SUPERSEDED BY D033. `MEASURED_F = 1.0` treats the visible inclusion body as
# the WHOLE of `a`. It is not: `a` is TOTAL aggregate, and the IbpA-marked focus is
# one object inside that pool. Diffuse and small oligomeric species are not in the
# focus and partition roughly evenly. Kept only so the D031 numbers remain
# reproducible; use `fEffFromBeta` instead.
MEASURED_F = 1.0

# beta = (aggregate mass in the visible inclusion body) / (total aggregate a).
# D031 assumed beta = 1 without saying so.
BETA_LADDER = (1.0, 0.75, 0.5, 0.25)


def partitionMultipliers(beta: float) -> Tuple[float, float]:
    """(old-pole, new-pole) concentration multipliers at division.

    A fraction `beta` of the aggregate sits in the indivisible focus and goes
    entirely to the old-pole daughter; the remaining `1 - beta` is well mixed and
    splits evenly. Each daughter has half the volume, so

        old : 2.[beta.a + 0.5(1-beta).a] = (1 + beta).a
        new : 2.[        0.5(1-beta).a ] = (1 - beta).a
    """
    b = float(beta)
    return 1.0 + b, 1.0 - b


def fEffFromBeta(beta: float) -> float:
    """the two-compartment rule IS the scalar rule, reparameterised.

    `2.f_eff = 1 + beta`, so `f_eff = 0.5 + 0.5.beta`. This is not a new model and
    adds no dynamics -- which is the point. `f` was never a free knob to be pinned
    at 1; it is DETERMINED by beta, and beta is a physical quantity that has to be
    measured. beta = 1 recovers D031's f = 1; beta = 0 recovers the f = 0.5
    control, at which the aging effect is exactly zero.
    """
    return 0.5 + 0.5 * float(beta)

# D015's measured slope: growth-rate loss per unit misfolded proteome fraction.
SLOPE_PER_PROTEOME_FRACTION = 32.0


# damping measured by `dampingForBeta` at each of these beta, recorded in D033.
# it is stored rather than recomputed at call time because each evaluation solves a
# ladder of folds and takes minutes; `tests/phase3/test_asymmetric_division.py`
# recomputes it and asserts these values, so it is a cache and not an orphan.
DAMPING_BY_BETA = ((0.05, 0.3544), (0.25, 0.3546), (0.50, 0.3549),
                   (0.75, 0.3552), (1.00, 0.3462))


def dampingAtBeta(beta: float) -> float:
    """the measured damping at `beta`, interpolated between measured points.

    Flat at 0.355 over the whole range except the approach to beta = 1, which is
    the only place the measurement moves. Using one number for all beta -- as an
    earlier version of the figure did -- puts the figure 1.5% off the section 8.4
    table at every row, which is a caption-against-text divergence with no cause
    other than convenience.
    """
    xs = [b for b, _ in DAMPING_BY_BETA]
    ys = [d for _, d in DAMPING_BY_BETA]
    return float(np.interp(float(beta), xs, ys))


def requiredAggregateFractionBeta(beta: float, damping: float,
                                  band=(BAND_LO, BAND_HI)):
    """the inverted requirement at focus share `beta`.

    `B_old - B_new = (1+beta).a - (1-beta).a = 2.beta.a`, so in the linear regime
    `aging = 64.beta.a_fraction`, damped by the measured exact/linear ratio. The
    requirement is therefore proportional to `1/beta` and rises without bound as
    the focus share falls -- LOWER beta needs MORE aggregate to make the same
    lineage difference.
    """
    lo, hi = band
    lead = 2.0 * float(beta) * SLOPE_PER_PROTEOME_FRACTION * float(damping)
    if lead <= 0.0:
        return (float("inf"), float("inf"))
    return (lo / lead, hi / lead)


def requiredAggregateFraction(band=(BAND_LO, BAND_HI), f: float = MEASURED_F,
                              damping: Optional[float] = None):
    """invert the band into the one quantity the prediction depends on.

    WHY THIS IS THE PARAMETER-FREE FORM. At small burden the hyperbolic law gives

        aging = (B_old - B_new)/k_mu

    and `k_mu = ARREST_PROTEOME_FRACTION/p_qc` while a model burden of 1 equals a
    proteome fraction of `p_qc`, so

        (B_old - B_new)/k_mu  ==  32 . (B_old - B_new)_as_proteome_fraction

    identically. **`p_qc` cancels.** `mu0` never appears. Both enter only through
    what the stationary aggregate load IS, so taking that load from measurement
    removes them, and the prediction collapses onto a single measurable number:
    the aggregate as a fraction of the proteome.

    At `f = 1` the daughters start at `2.a_end` and `0`, so the difference is
    `2.a_end` and `aging = 64 . a_end_fraction` in the linear regime. `damping`
    carries the exact-versus-linear correction measured by the machinery (the
    daughter's load relaxes during its generation, so the time-averaged
    difference is smaller than the initial one).
    """
    lo, hi = band
    lead = 2.0 * (2.0 * f - 1.0) * SLOPE_PER_PROTEOME_FRACTION
    d = 1.0 if damping is None else float(damping)
    return (lo / (lead * d), hi / (lead * d))


def dampingFactor(df: pd.DataFrame) -> float:
    """measured ratio of the exact aging effect to its linear-regime estimate."""
    import calibration as C
    s = df[(df.get("settled", False) == True) & (df["f"] > CONTROL_F)]  # noqa: E712
    if s.empty:
        return float("nan")
    k = np.array([C.kMuFromProteomeShare(p) for p in s["p_qc"]])
    lin = (s["a_old_start"].to_numpy() - s["a_new_start"].to_numpy()) / k
    return float(np.median(s["aging_effect"].to_numpy() / lin))


def dampingForBeta(beta: float, p: Optional[M.Params] = None,
                   pqc_ladder=PQC_LADDER, mu0_ladder=(0.05, 0.1),
                   j_fracs=(0.5, 0.8)) -> float:
    """exact/linear ratio, measured AT this beta rather than carried over.

    D032 quoted a single damping of 0.4386 measured under `f = 1`. Under
    two-compartment partitioning the new-pole daughter is no longer
    aggregate-free, so its relaxation over its own generation differs and the
    damping is beta-dependent. It is recomputed here per beta.
    """
    p = (p or M.Params()).validate()
    f = fEffFromBeta(beta)
    ratios = []
    for p_qc in pqc_ladder:
        g = calibratedHyperbolic(mu0_ladder[0], p_qc)
        for mu0 in mu0_ladder:
            g = calibratedHyperbolic(mu0, p_qc)
            fold = D.foldSolveDil(p, g)
            if fold is None:
                continue
            for jf in j_fracs:
                pj = p.with_(j=jf * fold[0]).validate()
                try:
                    r = lineageContrast(pj, g, f)
                except (M.ModelError, ValueError, OverflowError):
                    continue
                if r is None:
                    continue
                lin = (r["a_old_start"] - r["a_new_start"]) / g.k_mu
                if lin > 0.0:
                    ratios.append(r["aging_effect"] / lin)
    return float(np.median(ratios)) if ratios else float("nan")


def betaTable(betas=BETA_LADDER, p: Optional[M.Params] = None) -> pd.DataFrame:
    """the beta-indexed requirement, with damping recomputed at each beta."""
    rows = []
    for b in betas:
        dmp = dampingForBeta(b, p)
        lo, hi = requiredAggregateFractionBeta(b, dmp)
        rpo_lo, rpo_hi = TOMOYASU2001["rpoH_null_30C"]
        rows.append({"beta": b, "f_eff": fEffFromBeta(b), "damping": dmp,
                     "required_lo": lo, "required_hi": hi,
                     "required_lo_pct": 100 * lo, "required_hi_pct": 100 * hi,
                     "ratio_to_rpoH_null_lo": rpo_lo / hi if hi > 0 else np.nan,
                     "ratio_to_rpoH_null_hi": rpo_hi / lo if lo > 0 else np.nan})
    return pd.DataFrame(rows)
