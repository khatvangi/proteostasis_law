"""phase 3: the fold theorem, and phi decomposed into its two physical deficits.

THE THEOREM
-----------
write `G(u,a) = da/dt` and let `R(u,a)` be the TOTAL removal flux

    R = v_ref + v_degU + v_degA .

two structural facts about the model in `scripts/proteostasis/model.py`:

  (i)  the influx `j` enters `du/dt` and NOWHERE else.  the aggregate nullcline
       `{G = 0}` is therefore a curve in burden space that does not move when
       the load changes.
  (ii) mass balance (tested as A5) gives `du/dt + da/dt = j - R` exactly, because
       the internal u <-> a transfer cancels.

apply the determinant-preserving row operation `row1 -> row1 + row2` to the
equilibrium jacobian:

    det J = det [ d(du/dt)/du  d(du/dt)/da ] = det [ -R_u  -R_a ]
                [ d(da/dt)/du  d(da/dt)/da ]       [  G_u   G_a ]

            = -( R_u G_a - R_a G_u )  =  -( grad R x grad G ) .

three consequences follow immediately:

  * equilibria are exactly the points of the fixed curve `{G = 0}` at which
    `R = j`;
  * a saddle-node (`det J = 0`) occurs exactly where `grad R` is parallel to
    `grad G` -- i.e. at a CONSTRAINED CRITICAL POINT of R restricted to `{G=0}`;
  * hence  `j_crit = R(u*, a*)` at the first such point reached from the
    low-burden end, and

        THE COLLAPSE BOUNDARY IS WHERE TOTAL REMOVAL STOPS RESPONDING
        TO BURDEN ALONG THE AGGREGATE NULLCLINE.

the proven statement is CRITICAL POINT, not maximum.  an earlier version of this
docstring claimed "constrained maximum"; that is withdrawn.  R rises
monotonically along the first-root branch, and the solved fold has dG/da > 0, so
the critical point lies past the nullcline's turning point.  see the precision
note in `theory/FOLD_THEOREM.md`.

this is theorem-level: it is an identity, not a property of a sample.  it also
turns fold location into a 2x2 root solve in `(u, a)` with no continuation
sweep in `j`, which is what `foldSolve` below does.

WHAT phi IS
-----------
`phi = j_crit / removalCeiling`.  from the identity,

    R       = cf.s_ref + rho_U.df.s_u + rho_A.df.s_a
    ceiling = c_tot    + rho_U.d_tot  + rho_A.d_tot

so phi falls below 1 for exactly two reasons, and `phiDecomposition` separates
them by counterfactual:

  AVAILABILITY  cf < c_tot and df < d_tot -- rescue machinery is sequestered on
                nascent chains and on substrate at the collapse point;
  SATURATION    s_ref, s_u, s_a < 1      -- the michaelis factors are far from
                their asymptote when collapse happens.

CLAIM LABELS (see theory/SCOPE_AND_NONCLAIMS.md)
------------
  Mathematical  : the determinant identity and the constrained-critical-point
                  characterisation.  independent of sampling.
  Computational : every number produced by `verifyAgainstRun` and
                  `nestedInvariance` -- properties of this model over the
                  phase 1 parameter box.
  Empirical     : nothing here.  no organism data is used.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import brentq, root

REPO_ROOT = Path(__file__).resolve().parents[2]
if str(REPO_ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "scripts"))

from proteostasis import model as M  # noqa: E402

# the 14 kinetic parameters, i.e. everything except the loads j, nu and the
# rescue pools c_tot, d_tot
KINETIC_FIELDS = ("alpha_d", "alpha_g", "alpha_n", "kappa_a", "kappa_ca",
                  "kappa_cu", "kappa_da", "kappa_dis", "kappa_du", "kappa_ref",
                  "kappa_u", "m", "rho_A", "rho_U")


# ---------------------------------------------------------------------------
# the two scalar fields of the theorem
# ---------------------------------------------------------------------------


def removalR(u: float, a: float, p: M.Params) -> float:
    """total removal flux R = refolding + soluble degradation + aggregate clearance.

    this is the quantity that equals j at every equilibrium, by mass balance.
    """
    f = M.fluxes(u, a, p)
    return f["refold"] + f["degrade_u"] + f["degrade_a"]


def aggregateG(u: float, a: float, p: M.Params) -> float:
    """G = da/dt. the aggregate nullcline {G = 0} does not depend on j."""
    return M.rhs(u, a, p)[1]


def _centralGradient(fn, u: float, a: float, p: M.Params, h_rel: float = 1e-6):
    """central-difference gradient of a scalar field, relative step in each coord."""
    hu, ha = h_rel * u, h_rel * a
    gu = (fn(u + hu, a, p) - fn(u - hu, a, p)) / (2.0 * hu)
    ga = (fn(u, a + ha, p) - fn(u, a - ha, p)) / (2.0 * ha)
    return gu, ga


def determinantIdentity(u: float, a: float, p: M.Params) -> Dict[str, float]:
    """check det J == -(grad R x grad G) at a state.

    the analytic jacobian comes from the model package; the cross product is
    built from central differences, so the residual floor is the differencing
    error (~1e-7 relative), not the identity.
    """
    detJ = float(np.linalg.det(M.jacobian(u, a, p)))
    Ru, Ra = _centralGradient(removalR, u, a, p)
    Gu, Ga = _centralGradient(aggregateG, u, a, p)
    cross = Ru * Ga - Ra * Gu
    scale = max(abs(detJ), abs(cross), 1e-300)
    nR, nG = float(np.hypot(Ru, Ra)), float(np.hypot(Gu, Ga))
    return {
        "det_J": detJ,
        "minus_cross": -cross,
        "rel_err": abs(detJ + cross) / scale,
        # sin of the angle between the gradients: 0 exactly at the saddle-node
        "sin_angle": abs(cross) / max(nR * nG, 1e-300),
        "G": aggregateG(u, a, p),
    }


# ---------------------------------------------------------------------------
# solving for the fold directly from the theorem
# ---------------------------------------------------------------------------


def lowerNullclineA(u: float, p: M.Params, a_hi: float = 1e5) -> Optional[float]:
    """first root in `a` of G(u, .) = 0 at fixed u -- the branch stable in a.

    G > 0 at a = 0 (nucleation is the only term) and G > 0 for large a (growth
    outruns saturating removal), so a root exists only where the middle dips
    negative. the FIRST crossing has dG/da < 0.
    returns 0.0 if G <= 0 already at a = 0, or None if no crossing is found.

    THIS IS A BRACKETING HEURISTIC ONLY. it does NOT identify the branch the fold
    lives on: the solved fold state has dG/da > 0 and lies past the turning point
    of the nullcline, so it is not a first-root point. `foldSolve` uses this only
    to seed the root solve, and the solve itself is on {G = 0, det J = 0}.
    """
    if aggregateG(u, 0.0, p) <= 0.0:
        return 0.0
    a_prev, a = 0.0, 1e-10
    while a < a_hi:
        try:
            if aggregateG(u, a, p) <= 0.0:
                return float(brentq(lambda x: aggregateG(u, x, p), a_prev, a,
                                    xtol=1e-15, rtol=8.9e-16))
        except (M.ModelError, ValueError):
            return None
        a_prev, a = a, a * 1.5
    return None


def foldSolve(p: M.Params, u_lo: float = 1e-7, u_hi: float = 100.0,
              n_scan: int = 120) -> Optional[Tuple[float, float, float]]:
    """solve {G = 0, det J = 0} for the fold state; return (j_crit, u*, a*).

    this is the theorem used as an algorithm: no sweep in j is required because
    neither equation contains j. the coarse scan only supplies a bracket.
    """
    us = np.geomspace(u_lo, u_hi, n_scan)
    prev, guess = None, None
    for u in us:
        a = lowerNullclineA(u, p)
        if a is None:
            break
        try:
            d = float(np.linalg.det(M.jacobian(u, a, p)))
        except (M.ModelError, np.linalg.LinAlgError):
            break
        if prev is not None and np.sign(d) != np.sign(prev[1]):
            guess = (float(np.sqrt(prev[0] * u)), a)
            break
        prev = (u, d)
    if guess is None:
        if prev is None:
            return None
        a0 = lowerNullclineA(prev[0], p)
        guess = (prev[0], a0 if a0 else 1e-6)

    def residual(x):
        u, a = float(np.exp(x[0])), float(np.exp(x[1]))
        try:
            return [aggregateG(u, a, p),
                    float(np.linalg.det(M.jacobian(u, a, p)))]
        except (M.ModelError, np.linalg.LinAlgError):
            return [1e6, 1e6]

    sol = root(residual, [np.log(guess[0]), np.log(max(guess[1], 1e-12))],
               method="hybr", options={"xtol": 1e-13})
    if not sol.success:
        return None
    u_s, a_s = float(np.exp(sol.x[0])), float(np.exp(sol.x[1]))
    if not (np.isfinite(u_s) and np.isfinite(a_s) and u_s > 0.0 and a_s >= 0.0):
        return None
    return removalR(u_s, a_s, p), u_s, a_s


# ---------------------------------------------------------------------------
# decomposing phi
# ---------------------------------------------------------------------------


def phiDecomposition(u: float, a: float, p: M.Params) -> Dict[str, float]:
    """split phi at a fold state into the availability and saturation deficits.

    `phi_if_saturated`  forces every michaelis factor to 1 but keeps the real,
                        sequestered free pools -- so only availability is left.
    `phi_if_full_pools` keeps the real saturation but restores cf -> c_tot and
                        df -> d_tot -- so only saturation is left.
    """
    uf, af, cf, df = M.solveFreePools(u, a, p)
    s_ref = uf / (p.kappa_ref + uf)
    s_u = uf / (p.kappa_u + uf)
    s_a = af / (p.kappa_a + af)

    ceiling = M.removalCeiling(p)
    R_true = cf * s_ref + p.rho_U * df * s_u + p.rho_A * df * s_a
    R_noSat = cf + p.rho_U * df + p.rho_A * df
    R_noSeq = p.c_tot * s_ref + p.rho_U * p.d_tot * s_u + p.rho_A * p.d_tot * s_a
    return {
        "phi": R_true / ceiling,
        "phi_if_saturated": R_noSat / ceiling,
        "phi_if_full_pools": R_noSeq / ceiling,
        "cf_frac": cf / p.c_tot, "df_frac": df / p.d_tot,
        "s_ref": s_ref, "s_u": s_u, "s_a": s_a,
    }


# ---------------------------------------------------------------------------
# parameter reconstruction from the phase 1 experiment C table
# ---------------------------------------------------------------------------


def paramsFromSampleRow(row) -> M.Params:
    """rebuild a Params from one experiment C row.

    experiment C samples the rescue ALLOCATION `p_chi` at a fixed total rescue of
    1.0, which is the Params default c_tot + d_tot. the rebuild check in
    `verifyAgainstRun` is what certifies this reading: if the split were wrong
    the reconstructed phi would not reproduce `C2_fold_to_ceiling_ratio`.
    """
    kin = {f: float(row["p_" + f]) for f in KINETIC_FIELDS}
    chi = float(row["p_chi"])
    return M.Params(**kin).with_(nu=float(row["p_nu"]),
                                 c_tot=chi, d_tot=1.0 - chi).validate()


def foldStateFromSampleRow(row) -> Tuple[float, float]:
    """recover (u*, a*) from the recorded burden and aggregate fraction."""
    burden = float(row["fold_burden"])
    frac = float(row["fold_aggregate_fraction"])
    return burden * (1.0 - frac), burden * frac


# ---------------------------------------------------------------------------
# artefact-dependent verification
# ---------------------------------------------------------------------------


def phase1RunDir() -> Path:
    """resolve the phase 1 run root; results/ is gitignored so it may be absent."""
    env = os.environ.get("PHASE1_RUN_DIR")
    if env and Path(env).is_dir():
        return Path(env)
    root = REPO_ROOT / "results" / "phase1"
    cands = sorted(root.glob("run_*")) if root.is_dir() else []
    return cands[-1] if cands else root / "__missing__"


def verifyAgainstRun(run: Optional[Path] = None, n_identity: int = 20) -> Dict[str, float]:
    """recompute every phase 3 headline from the phase 1 artefacts."""
    run = run or phase1RunDir()
    b = pd.read_csv(run / "B" / "fold_boundary.tsv", sep="\t")
    b = b[b["found"] == True]  # noqa: E712
    base = M.Params()

    # (1) the determinant identity at recorded fold states
    ident, para = [], []
    for _, r in b.sample(n=n_identity, random_state=1).iterrows():
        p = M.allocationParams(base.with_(nu=float(r["nu"])), float(r["chi"])).validate()
        d = determinantIdentity(float(r["fold_u"]), float(r["fold_a"]), p)
        ident.append(d["rel_err"])
        # parallelism should scale with how close the bracketed state is to the
        # true saddle-node, measured by the recorded leading eigenvalue
        para.append((d["sin_angle"], abs(float(r["fold_eig_real_max"]))))

    # (2) the solver reproduces the continuation-derived fold
    soln = []
    for _, r in b.sample(n=n_identity, random_state=3).iterrows():
        p = M.allocationParams(base.with_(nu=float(r["nu"])), float(r["chi"])).validate()
        out = foldSolve(p)
        if out is not None:
            soln.append(abs(out[0] - float(r["fold_value"])) / float(r["fold_value"]))

    # (3) phi rebuilds from first principles across every experiment C fold
    c = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    c = c[pd.to_numeric(c["fold_burden"], errors="coerce").notna()]
    rebuild, dec = [], []
    for _, r in c.iterrows():
        try:
            p = paramsFromSampleRow(r)
            u, a = foldStateFromSampleRow(r)
            d = phiDecomposition(u, a, p)
        except (M.ModelError, ValueError, KeyError):
            continue
        rebuild.append(abs(d["phi"] - float(r["C2_fold_to_ceiling_ratio"])))
        dec.append(d)
    D = pd.DataFrame(dec)
    short = 1.0 - D["phi"]

    return {
        "n_identity": len(ident),
        "identity_rel_err_median": float(np.median(ident)),
        "parallelism_tracks_eigenvalue": float(np.corrcoef(
            np.log10([x[0] for x in para]), np.log10([x[1] for x in para]))[0, 1]),
        "n_solver": len(soln),
        "solver_rel_err_max": float(np.max(soln)),
        "n_phi": len(rebuild),
        "phi_rebuild_err_median": float(np.median(rebuild)),
        "phi_rebuild_err_max": float(np.max(rebuild)),
        "phi_median": float(D["phi"].median()),
        "s_ref_median": float(D["s_ref"].median()),
        "s_u_median": float(D["s_u"].median()),
        "s_a_median": float(D["s_a"].median()),
        "shortfall_share_saturation": float(
            ((D["phi_if_saturated"] - D["phi"]) / short).median()),
        "shortfall_share_sequestration": float(
            ((D["phi_if_full_pools"] - D["phi"]) / short).median()),
    }


def nestedInvariance(run: Optional[Path] = None, k_theta: int = 12,
                     n_nu: int = 7, n_chi: int = 7, seed: int = 11) -> pd.DataFrame:
    """the properly nested design: k kinetic draws x a (nu, chi) load grid.

    the phase 1 latin hypercube varied loads and kinetics together, so it cannot
    say whether phi is a property of the network or of how the network is
    loaded. this crosses them explicitly. affordable only because `foldSolve`
    replaced the continuation sweep.
    """
    run = run or phase1RunDir()
    c = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    draws = c.sample(n=k_theta, random_state=seed)

    nus = np.geomspace(0.01, 20.0, n_nu)
    chis = np.linspace(0.15, 0.85, n_chi)
    rec = []
    for idx, (_, r) in enumerate(draws.iterrows()):
        kin = {f: float(r["p_" + f]) for f in KINETIC_FIELDS}
        for nu in nus:
            for chi in chis:
                try:
                    p = M.Params(**kin).with_(nu=float(nu), c_tot=float(chi),
                                              d_tot=1.0 - float(chi)).validate()
                    out = foldSolve(p)
                except (M.ModelError, ValueError):
                    continue
                if out is None or not np.isfinite(out[0]) or out[0] <= 0.0:
                    continue
                rec.append(dict(theta=idx, nu=float(nu), chi=float(chi),
                                phi=out[0] / M.removalCeiling(p)))
    D = pd.DataFrame(rec)
    D["lphi"] = np.log10(D["phi"])
    return D


def invarianceSummary(D: pd.DataFrame, n_cells: int = 49) -> Dict[str, float]:
    """within-theta vs between-theta spread, and the split by load coordinate."""
    full = D.groupby("theta").filter(lambda g: len(g) >= 0.6 * n_cells)
    g = full.groupby("theta")["lphi"]
    within, between = g.std(ddof=1).median(), g.mean().std(ddof=1)
    chi_at_nu = full.groupby(["theta", "nu"])["lphi"].std(ddof=1).groupby("theta").median()
    nu_at_chi = full.groupby(["theta", "chi"])["lphi"].std(ddof=1).groupby("theta").median()
    return {
        "n_theta": int(full["theta"].nunique()),
        "n_folds": int(len(full)),
        "within_sd": float(within),
        "between_sd": float(between),
        "variance_ratio": float(between ** 2 / within ** 2),
        "spread_across_chi_at_fixed_nu": float(10 ** (2 * chi_at_nu.median())),
        "spread_across_nu_at_fixed_chi": float(10 ** (2 * nu_at_chi.median())),
        "spread_between_theta": float(10 ** (2 * between)),
    }


def main() -> int:
    run = phase1RunDir()
    if not run.is_dir():
        print("SKIP: no phase 1 run root found (results/ is gitignored). "
              "set PHASE1_RUN_DIR to verify.")
        return 0
    print("phase 1 run: %s" % run)
    print()
    v = verifyAgainstRun(run)
    print("THEOREM  det J == -(grad R x grad G)")
    print("  states checked                    : %d" % v["n_identity"])
    print("  median relative error             : %.3e  (central-difference floor)"
          % v["identity_rel_err_median"])
    print("  corr(log sin-angle, log |eig|)    : %+.4f  (-> residual is bracket "
          "tolerance, not the identity)" % v["parallelism_tracks_eigenvalue"])
    print()
    print("SOLVER  {G = 0, det J = 0} vs the continuation sweep")
    print("  folds compared                    : %d" % v["n_solver"])
    print("  max relative error in j_crit      : %.3e" % v["solver_rel_err_max"])
    print()
    print("phi  rebuilt from first principles")
    print("  folds rebuilt                     : %d" % v["n_phi"])
    print("  median | max absolute error       : %.3e | %.3e"
          % (v["phi_rebuild_err_median"], v["phi_rebuild_err_max"]))
    print("  median phi                        : %.4f" % v["phi_median"])
    print("  saturation at collapse s_ref/s_u/s_a: %.4f / %.4f / %.4f"
          % (v["s_ref_median"], v["s_u_median"], v["s_a_median"]))
    print("  shortfall share, saturation       : %.1f%%"
          % (100 * v["shortfall_share_saturation"]))
    print("  shortfall share, sequestration    : %.1f%%"
          % (100 * v["shortfall_share_sequestration"]))
    print()
    D = nestedInvariance(run)
    s = invarianceSummary(D)
    print("NESTED phi-INVARIANCE  %d networks, %d folds solved" % (s["n_theta"], s["n_folds"]))
    print("  spread across chi at fixed nu     : %.2fx" % s["spread_across_chi_at_fixed_nu"])
    print("  spread across nu at fixed chi     : %.2fx" % s["spread_across_nu_at_fixed_chi"])
    print("  spread between networks           : %.2fx" % s["spread_between_theta"])
    print("  variance ratio between/within     : %.1f" % s["variance_ratio"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
