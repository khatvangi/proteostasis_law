"""phase 3, antecedent check A1: does the machinery damage itself, and does the
determinant identity survive it?

THE WORRY
---------
The whole derivation runs on two structural facts (`fold_theorem.py`):

  (i)  the damage influx `j` enters `du/dt` and nowhere else;
  (ii) mass balance `du/dt + da/dt = j - R` is exact.

Fact (i) is stated as though influx and clearance capacity were independent. In a
cell they are not. Chaperones and proteases are themselves translated, at the
same per-codon error rate that produces the damage, so raising `j` degrades the
machinery that clears it. If that coupling breaks

        det J = -( grad R x grad G )

then the theorem holds for a model class that EXCLUDES cells, and that belongs at
the top of any manuscript rather than in a referee report.

THE TEST
--------
Enzymatic capacity is made error-dependent with a form that cannot go negative:

        C_enz(load) = C_0 / (1 + eps . load),        eps = 1 / j_ref

applied to BOTH conserved pools, since both are translated. `eps = 0` must
recover the frozen model exactly (asserted at 0.0, as the dilution work does).
`j_ref` is NOT tuned to any biological value -- this is a structural check.

Two modes, and the distinction is the whole point:

  INFLUX mode   load = j.  This is the stated worry: capacity falls with the
                error rate. `j` is a PARAMETER, so this changes the equations but
                not their dependence on the state.
  BURDEN mode   load = u + a.  Capacity falls with the damage currently carried
                (occupancy, sequestration of machinery into inclusions). Here the
                capacity enters `grad R` and `grad G` themselves. This is the
                mode that can actually break the identity, and it is run for that
                reason even though the stated worry is the influx one.

WHAT THE COUPLING DOES BREAK, EVEN WHERE THE IDENTITY SURVIVES
--------------------------------------------------------------
Fact (i) had a second half: the aggregate nullcline `{G = 0}` is a FIXED curve
that does not move when the load changes. Under influx-mode self-damage it moves,
because the pools that appear in `G` now depend on `j`. So `j_crit = R(u*,a*)`
stops being an evaluation and becomes a self-consistency condition, and the fold
solve grows from two equations in `(u,a)` to three in `(u,a,j)`:

        G(u,a; j) = 0,    det J(u,a; j) = 0,    R(u,a; j) = j .

That is a real loss -- the shortcut, not the theorem.

CLAIM LABELS (see theory/SCOPE_AND_NONCLAIMS.md)
  Mathematical  : the determinant identity and why it is indifferent to how the
                  parameters depend on `j`.
  Computational : every number below -- a property of this model over the
                  phase 1 parameter box.
  Empirical     : nothing. `eps` is not fitted to any measurement.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import root

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import fold_theorem as FT  # noqa: E402

INFLUX = "influx"
BURDEN = "burden"

# ladder spanning four decades of eps. at j ~ 0.2 this runs from a 0.2% capacity
# loss to a 95% loss, so the far end is not a perturbation.
EPS_LADDER = (1e-2, 1e-1, 1e0, 1e1, 1e2)


@dataclass(frozen=True)
class SelfDamage:
    """capacity self-damage. `eps = 0` is the frozen model, exactly."""

    eps: float = 0.0
    mode: str = INFLUX

    def validate(self) -> "SelfDamage":
        if not (self.eps >= 0.0 and np.isfinite(self.eps)):
            raise M.ModelError("eps must be finite and nonnegative")
        if self.mode not in (INFLUX, BURDEN):
            raise M.ModelError(f"unknown self-damage mode '{self.mode}'")
        return self

    def factor(self, j: float, u: float, a: float) -> float:
        """C_enz/C_0 in [0,1]. bounded below by 0 and never negative."""
        if self.eps == 0.0:
            return 1.0
        load = j if self.mode == INFLUX else (u + a)
        return 1.0 / (1.0 + self.eps * max(load, 0.0))


def effectiveParams(p: M.Params, u: float, a: float, sd: SelfDamage) -> M.Params:
    """the params actually in force at this state, with both pools scaled.

    at eps = 0 this returns `p` unchanged, which is what makes the ladder a
    strict generalisation rather than a reparameterisation.
    """
    f = sd.factor(p.j, u, a)
    if f == 1.0:
        return p
    return p.with_(c_tot=p.c_tot * f, d_tot=p.d_tot * f)


# ---------------------------------------------------------------------------
# the two scalar fields, and the state derivatives, under self-damage
# ---------------------------------------------------------------------------


def removalRsd(u: float, a: float, p: M.Params, sd: SelfDamage) -> float:
    return FT.removalR(u, a, effectiveParams(p, u, a, sd))


def aggregateGsd(u: float, a: float, p: M.Params, sd: SelfDamage) -> float:
    return FT.aggregateG(u, a, effectiveParams(p, u, a, sd))


def rhsSd(u: float, a: float, p: M.Params, sd: SelfDamage) -> Tuple[float, float]:
    du, da = M.rhs(u, a, effectiveParams(p, u, a, sd))
    return float(du), float(da)


def _grad(fn, u: float, a: float, p: M.Params, sd: SelfDamage, h_rel: float = 1e-6):
    hu, ha = h_rel * u, h_rel * a
    gu = (fn(u + hu, a, p, sd) - fn(u - hu, a, p, sd)) / (2.0 * hu)
    ga = (fn(u, a + ha, p, sd) - fn(u, a - ha, p, sd)) / (2.0 * ha)
    return gu, ga


def numericalJacobianSd(u: float, a: float, p: M.Params, sd: SelfDamage,
                        h_rel: float = 1e-6) -> np.ndarray:
    """central-difference jacobian of the state equations at FIXED j.

    an analytic jacobian is unavailable in BURDEN mode because the model package
    assumes the conserved pools are constants. differencing both sides of the
    identity keeps the comparison apples-to-apples: the eps = 0 row of every
    ladder is the differencing floor, and is what the eps > 0 rows are read
    against.
    """
    hu, ha = h_rel * u, h_rel * a
    du_p, da_p = rhsSd(u + hu, a, p, sd)
    du_m, da_m = rhsSd(u - hu, a, p, sd)
    col_u = ((du_p - du_m) / (2.0 * hu), (da_p - da_m) / (2.0 * hu))
    du_p, da_p = rhsSd(u, a + ha, p, sd)
    du_m, da_m = rhsSd(u, a - ha, p, sd)
    col_a = ((du_p - du_m) / (2.0 * ha), (da_p - da_m) / (2.0 * ha))
    return np.array([[col_u[0], col_a[0]], [col_u[1], col_a[1]]])


def identityAt(u: float, a: float, p: M.Params, sd: SelfDamage,
               h_rel: float = 1e-6) -> Dict[str, float]:
    """check det J == -(grad R x grad G) under self-damage, both sides numerical.

    TWO error measures, and the second is the one to read.

    `rel_err` divides by max(|det J|, |cross|), which VANISHES at a saddle-node.
    Holding the evaluation state fixed while raising eps moves that state away
    from criticality, so |det J| grows and `rel_err` falls -- the identity would
    appear to IMPROVE under stronger coupling, which is an artefact of the
    denominator, not a property of the identity.

    `rel_err_grad` divides by |grad R| . |grad G|, which stays O(1) through the
    fold (only the SINE of the angle between them vanishes there). That is the
    scale-stable measure and it is what the ladder slope is fitted on.
    """
    detJ = float(np.linalg.det(numericalJacobianSd(u, a, p, sd, h_rel)))
    Ru, Ra = _grad(removalRsd, u, a, p, sd, h_rel)
    Gu, Ga = _grad(aggregateGsd, u, a, p, sd, h_rel)
    cross = Ru * Ga - Ra * Gu
    err = abs(detJ + cross)
    nR, nG = float(np.hypot(Ru, Ra)), float(np.hypot(Gu, Ga))
    return {"det_J": detJ, "minus_cross": -cross,
            "rel_err": err / max(abs(detJ), abs(cross), 1e-300),
            "rel_err_grad": err / max(nR * nG, 1e-300),
            "sin_angle": abs(cross) / max(nR * nG, 1e-300),
            "capacity_factor": sd.factor(p.j, u, a)}


# ---------------------------------------------------------------------------
# the fold, now a self-consistency problem in (u, a, j)
# ---------------------------------------------------------------------------


def foldSolveSd(p: M.Params, sd: SelfDamage,
                seed: Optional[Tuple[float, float, float]] = None
                ) -> Optional[Tuple[float, float, float]]:
    """solve {G = 0, det J = 0, R = j} for (j_crit, u*, a*) under self-damage.

    the frozen model solves only the first two, because neither contains `j`,
    and then READS OFF j_crit = R(u*,a*). that shortcut dies here: capacity
    depends on `j`, so the nullcline moves with the load and the third equation
    has to be imposed rather than evaluated.
    """
    if seed is None:
        base = FT.foldSolve(p)
        if base is None:
            return None
        seed = base
    j0, u0, a0 = seed

    def residual(x):
        u, a, j = (float(np.exp(np.clip(x[0], -40.0, 8.0))),
                   float(np.exp(np.clip(x[1], -40.0, 8.0))),
                   float(np.exp(np.clip(x[2], -40.0, 8.0))))
        pj = p.with_(j=j)
        try:
            pe = effectiveParams(pj, u, a, sd)
            return [FT.aggregateG(u, a, pe),
                    float(np.linalg.det(numericalJacobianSd(u, a, pj, sd))),
                    FT.removalR(u, a, pe) - j]
        except (M.ModelError, OverflowError, np.linalg.LinAlgError,
                FloatingPointError):
            return [1e6, 1e6, 1e6]

    s = root(residual, [np.log(max(u0, 1e-14)), np.log(max(a0, 1e-14)),
                        np.log(max(j0, 1e-14))],
             method="hybr", options={"xtol": 1e-12})
    if not s.success:
        return None
    u, a, j = (float(np.exp(np.clip(s.x[0], -40.0, 8.0))),
               float(np.exp(np.clip(s.x[1], -40.0, 8.0))),
               float(np.exp(np.clip(s.x[2], -40.0, 8.0))))
    if not all(np.isfinite(v) and v > 0.0 for v in (u, a, j)):
        return None
    if max(abs(v) for v in residual([np.log(u), np.log(a), np.log(j)])) > 1e-7:
        return None
    return j, u, a


# ---------------------------------------------------------------------------
# critical slowing under self-damage
# ---------------------------------------------------------------------------


def influxCeilingSelfDamaged(p: M.Params, sd: SelfDamage) -> float:
    """largest influx admitting ANY bounded state, in closed form (INFLUX mode).

    The frozen model has the analytic necessary condition `j <= C_0`, with
    `C_0 = c_tot + (rho_U + rho_A).d_tot` (model.py). Under influx self-damage
    every removal flux carries the factor `1/(1 + eps.j)`, so the condition
    becomes `j <= C_0/(1 + eps.j)`, i.e.

        eps.j^2 + j - C_0 <= 0   ->   j <= ( sqrt(1 + 4.eps.C_0) - 1 ) / (2.eps)

    which tends to `sqrt(C_0/eps)` for large eps. So self-damage converts a LINEAR
    capacity ceiling into a SQUARE-ROOT one -- doubling the machinery buys only
    sqrt(2) in tolerable error rate once the machinery is itself error-prone.

    This is exact and it is new, but it is LOOSE: measured against the solved
    folds it is never violated and never binding (see D027). It is reported as a
    necessary condition, not as the boundary.

    In BURDEN mode the factor depends on the state, so no closed form follows;
    the frozen ceiling still bounds, since the factor never exceeds 1.
    """
    C0 = M.removalCeiling(p)
    if sd.mode != INFLUX or sd.eps == 0.0:
        return float(C0)
    return float((np.sqrt(1.0 + 4.0 * sd.eps * C0) - 1.0) / (2.0 * sd.eps))


def stableEquilibriumSd(j: float, p: M.Params, sd: SelfDamage,
                        guess=(1e-3, 1e-6)) -> Optional[Tuple[float, float, float]]:
    """low-burden equilibrium at influx j, with capacity evaluated at that j."""
    pj = p.with_(j=float(j))

    def res(x):
        u = float(np.exp(np.clip(x[0], -34.0, 6.0)))
        a = float(np.exp(np.clip(x[1], -34.0, 6.0)))
        try:
            du, da = rhsSd(u, a, pj, sd)
        except (M.ModelError, OverflowError, FloatingPointError):
            return [1e6, 1e6]
        if not (np.isfinite(du) and np.isfinite(da)):
            return [1e6, 1e6]
        return [du, da]

    s = root(res, [np.log(max(guess[0], 1e-12)), np.log(max(guess[1], 1e-14))],
             method="hybr", options={"xtol": 1e-13})
    if not s.success:
        return None
    u = float(np.exp(np.clip(s.x[0], -34.0, 6.0)))
    a = float(np.exp(np.clip(s.x[1], -34.0, 6.0)))
    try:
        du, da = rhsSd(u, a, pj, sd)
        if max(abs(du), abs(da)) > 1e-9 * max(1.0, u + a):
            return None
        ev = np.linalg.eigvals(numericalJacobianSd(u, a, pj, sd))
    except (M.ModelError, OverflowError, np.linalg.LinAlgError):
        return None
    lead = float(np.max(ev.real))
    if lead >= 0.0:
        return None
    return u, a, lead


def slowingExponentSd(p: M.Params, sd: SelfDamage, n: int = 14,
                      lo: float = 1e-4, hi: float = 3e-2) -> Optional[Dict[str, float]]:
    """fit log|lambda| against log(j_crit - j). saddle-node normal form gives 1/2.

    a departure from 1/2 would be a NEW prediction, not a repair: it would say
    self-damage makes collapse steeper than a generic saddle-node.
    """
    out = foldSolveSd(p, sd)
    if out is None:
        return None
    j_crit, u_f, a_f = out
    ds, lams = [], []
    guess = (u_f, max(a_f, 1e-12))
    for rel in np.geomspace(lo, hi, n):
        j = j_crit * (1.0 - rel)
        eq = stableEquilibriumSd(j, p, sd, guess)
        if eq is None:
            continue
        u, a, lead = eq
        guess = (u, max(a, 1e-14))
        ds.append(j_crit - j)
        lams.append(abs(lead))
    if len(ds) < 6:
        return None
    x, y = np.log10(np.array(ds)), np.log10(np.array(lams))
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    r2 = 1.0 - float(np.sum(resid ** 2) / max(np.sum((y - y.mean()) ** 2), 1e-300))
    return {"j_crit": j_crit, "n_points": len(ds), "slope": float(slope),
            "r2": r2, "tau_exponent": -float(slope)}


# ---------------------------------------------------------------------------
# the ladder, run over the phase 1 parameter box
# ---------------------------------------------------------------------------


def _networks(k: int, seed: int) -> list:
    """k sampled networks from the phase 1 experiment C table, else the base point."""
    run = FT.phase1RunDir()
    path = run / "C" / "samples.tsv"
    if not path.is_file():
        return [M.Params().validate()]
    c = pd.read_csv(path, sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    c = c[pd.to_numeric(c["fold_burden"], errors="coerce").notna()]
    out = []
    for _, r in c.sample(n=min(k, len(c)), random_state=seed).iterrows():
        try:
            out.append((FT.paramsFromSampleRow(r), FT.foldStateFromSampleRow(r)))
        except (M.ModelError, ValueError, KeyError):
            continue
    return out


# differencing step. 1e-6 is the repo's habit elsewhere; here a LARGE step is
# correct, and the reason is structural rather than a tolerance tweak.
#
# `du/dt = j - R - G` holds POINTWISE, so the central-difference operator -- being
# linear -- reproduces the row operation exactly at every step size: the
# differenced d(du/dt)/du equals -(differenced R_u) - (differenced G_u) with no
# remainder. The identity therefore carries NO truncation term, only roundoff,
# which falls as h grows. `refinementCheck` measures slope -0.97 in h across four
# decades with no V-shaped minimum, which is that prediction confirmed.
#
# The consequence is worth stating plainly: this numerical check cannot fail
# unless mass balance itself fails, so it certifies the IMPLEMENTATION, and the
# analytic argument carries the theorem. See D027.
IDENTITY_H_REL = 1e-2


def identityLadder(k: int = 20, seed: int = 11, mode: str = INFLUX,
                   ladder=EPS_LADDER, h_rel: float = IDENTITY_H_REL) -> pd.DataFrame:
    """median relative error of the identity as a function of eps.

    the eps = 0 row is the control: it is the differencing floor, and the
    eps > 0 rows are only interpretable against it.
    """
    nets = _networks(k, seed)
    rows = []
    for eps in (0.0,) + tuple(ladder):
        sd = SelfDamage(eps=eps, mode=mode).validate()
        errs = []
        for item in nets:
            if isinstance(item, tuple):
                p, (u, a) = item
            else:
                p, (u, a) = item, (0.05, 0.05)
            try:
                errs.append(identityAt(u, a, p, sd, h_rel=h_rel))
            except (M.ModelError, OverflowError, np.linalg.LinAlgError):
                continue
        if not errs:
            continue
        rows.append({"mode": mode, "eps": eps, "n": len(errs),
                     "rel_err_median": float(np.median([e["rel_err"] for e in errs])),
                     "rel_err_grad_median":
                         float(np.median([e["rel_err_grad"] for e in errs])),
                     "rel_err_grad_max":
                         float(np.max([e["rel_err_grad"] for e in errs])),
                     "capacity_factor":
                         float(np.median([e["capacity_factor"] for e in errs]))})
    return pd.DataFrame(rows)


def ladderSlope(df: pd.DataFrame, col: str = "rel_err_grad_median") -> Dict[str, float]:
    """log-log slope of median identity error against eps, over eps > 0.

    slope ~ 0  -> the error does not grow with the coupling: identity survives
    slope ~ 1  -> first-order valid only; the coupling degrades it linearly
    """
    d = df[df["eps"] > 0.0]
    if len(d) < 3:
        return {"slope": float("nan"), "r2": float("nan"), "floor": float("nan")}
    x = np.log10(d["eps"].to_numpy())
    y = np.log10(np.maximum(d[col].to_numpy(), 1e-300))
    slope, b = np.polyfit(x, y, 1)
    resid = y - (slope * x + b)
    denom = max(float(np.sum((y - y.mean()) ** 2)), 1e-300)
    floor = df.loc[df["eps"] == 0.0, col]
    return {"slope": float(slope),
            "r2": 1.0 - float(np.sum(resid ** 2) / denom),
            "floor": float(floor.iloc[0]) if len(floor) else float("nan")}


def refinementCheck(k: int = 20, seed: int = 11, mode: str = BURDEN,
                    eps: float = 100.0,
                    steps=(1e-4, 1e-5, 1e-6, 1e-7)) -> pd.DataFrame:
    """is a residual trend in the ladder identity breakdown, or differencing error?

    THE DISCRIMINATOR. Central differencing has truncation error O(h^2). If the
    residual shrinks like h^2 as the step is refined, it is the numerics; if it
    is flat in h, it is the identity. Strong state-dependence of capacity (large
    eps in BURDEN mode) inflates the third derivatives that set the truncation
    constant, which is exactly the regime where the two explanations have to be
    told apart rather than assumed.
    """
    nets = _networks(k, seed)
    sd = SelfDamage(eps=eps, mode=mode).validate()
    rows = []
    for h in steps:
        errs = []
        for item in nets:
            p, (u, a) = item if isinstance(item, tuple) else (item, (0.05, 0.05))
            try:
                errs.append(identityAt(u, a, p, sd, h_rel=h)["rel_err_grad"])
            except (M.ModelError, OverflowError, np.linalg.LinAlgError):
                continue
        if errs:
            rows.append({"mode": mode, "eps": eps, "h_rel": h, "n": len(errs),
                         "rel_err_grad_median": float(np.median(errs))})
    df = pd.DataFrame(rows)
    if len(df) >= 3:
        x = np.log10(df["h_rel"].to_numpy())
        y = np.log10(np.maximum(df["rel_err_grad_median"].to_numpy(), 1e-300))
        df.attrs["slope_in_h"] = float(np.polyfit(x, y, 1)[0])
    return df


def directionLadder(k: int = 8, seed: int = 11, mode: str = INFLUX,
                    ladder=EPS_LADDER, with_exponent: bool = True,
                    n_sub: int = 6) -> pd.DataFrame:
    """does self-damage lower j_crit, and does the approach steepen?

    reported independently of the identity: even where the identity survives
    exactly, the boundary can move and the exponent can change.

    CONTINUATION, not re-seeding. each rung is seeded from the previous rung's
    solution, with `n_sub` intermediate steps, so a reported "no fold" means the
    branch was followed until it ended rather than that a solve launched from the
    frozen fold failed to reach a distant root. seeding every rung from eps = 0
    loses folds the continuation finds.
    """
    nets = _networks(k, seed)
    rows = []
    for item in nets:
        p = item[0] if isinstance(item, tuple) else item
        base = foldSolveSd(p, SelfDamage(0.0, mode).validate())
        if base is None:
            continue
        e0 = (slowingExponentSd(p, SelfDamage(0.0, mode).validate())
              if with_exponent else None)
        cur, prev_eps = base, 0.0
        for eps in ladder:
            sd = SelfDamage(eps=eps, mode=mode).validate()
            out = None
            if cur is not None:
                sub = np.linspace(prev_eps, eps, n_sub + 1)[1:]
                walk = cur
                for e in sub:
                    step = foldSolveSd(p, SelfDamage(float(e), mode).validate(),
                                       seed=walk)
                    if step is None:
                        walk = None
                        break
                    walk = step
                out = walk
            if out is None:                       # fall back to a direct solve
                out = foldSolveSd(p, sd, seed=base)
            cur, prev_eps = out, eps
            if out is None:
                rows.append({"mode": mode, "eps": eps, "fold_exists": False,
                             "j_crit": np.nan, "j_crit_ratio": np.nan,
                             "tau_exponent": np.nan,
                             "tau_exponent_frozen": (e0 or {}).get("tau_exponent",
                                                                   np.nan)})
                continue
            ex = slowingExponentSd(p, sd) if with_exponent else None
            # the identity AT the self-consistent fold: sin_angle must vanish
            # there, which is what makes the solved point a saddle-node at all
            idf = identityAt(out[1], out[2], p.with_(j=out[0]), sd)
            rows.append({"mode": mode, "eps": eps, "fold_exists": True,
                         "j_crit": out[0], "j_crit_ratio": out[0] / base[0],
                         "fold_sin_angle": idf["sin_angle"],
                         "fold_rel_err_grad": idf["rel_err_grad"],
                         "tau_exponent": (ex or {}).get("tau_exponent", np.nan),
                         "tau_exponent_frozen": (e0 or {}).get("tau_exponent",
                                                               np.nan)})
    return pd.DataFrame(rows)


def main():
    pd.set_option("display.width", 130)
    for mode in (INFLUX, BURDEN):
        print(f"\n=== identity ladder, mode = {mode} "
              f"({'j is a parameter' if mode == INFLUX else 'capacity is state-dependent'})")
        df = identityLadder(mode=mode)
        print(df.to_string(index=False))
        print("  log-log fit over eps > 0:", ladderSlope(df))

    print("\n=== direction: does self-damage lower j_crit and steepen the approach?")
    for mode in (INFLUX, BURDEN):
        d = directionLadder(mode=mode)
        if d.empty:
            print(f"  {mode}: no folds solved")
            continue
        g = d.groupby("eps").agg(
            n=("fold_exists", "size"),
            folds=("fold_exists", "sum"),
            j_crit_ratio_median=("j_crit_ratio", "median"),
            fold_sin_angle_max=("fold_sin_angle", "max"),
            tau_median=("tau_exponent", "median"))
        frozen = d.groupby("eps")["tau_exponent_frozen"].median().median()
        print(f"\n  mode = {mode} (frozen exponent median {frozen:.4f})")
        print(g.to_string())


if __name__ == "__main__":
    main()
