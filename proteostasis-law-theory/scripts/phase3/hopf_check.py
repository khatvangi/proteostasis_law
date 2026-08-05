"""Block 1b verification: are the tr J >= 0 branches real, or a tracing artefact?

`genericity.py:branchToFold` reported tr J >= 0 somewhere on the stable branch
for 101 of 2765 kinetic-box networks, with a maximum of +6.97. Three of the four
Block 1 numbers before it were method artefacts, and this one has two candidate
mechanisms that would produce exactly the same output:

  (A) the segment ordering. `branchToFold` decides which end of the contour is
      low burden by `u + a`, a sum of two quantities on different scales, and
      then keeps everything from that end to the fold. If it picks the wrong
      end it walks the unstable arc.
  (B) contour resolution. The vertices come from a 150x150 field contour and are
      linearly interpolated between grid nodes, so they do not lie on {G = 0}.
      tr J is evaluated wherever they happen to land.

THIS MODULE REMOVES BOTH, then applies a check that shares no failure mode with
the tracer at all:

  1. ordering, removed.   The accessible branch is the maximal contiguous run of
     det J > 0 adjacent to the solved fold. No heuristic decides which end is
     which; connectivity to the fold decides it.
  2. resolution, removed. Every vertex is projected onto {G = 0} along grad G
     before its Jacobian is evaluated, so each point IS an equilibrium: with
     G = 0 and j := R, both du/dt = -G = 0 and da/dt = G = 0 exactly.
  3. dynamics.  At the point of largest tr J, integrate the full nonlinear
     system from a small perturbation and measure the growth exponent. If the
     equilibrium is unstable the perturbation grows, at a rate the Jacobian
     predicts. This test never touches the contour, the ordering, or the field
     grid, which is what the last four incidents all lacked.

A tr J = 0 crossing with det J > 0 has pure imaginary eigenvalues, so det at the
crossing is reported beside it.

THAT IS NECESSARY AND NOT SUFFICIENT, and this docstring used to say "a Hopf
point by definition, not an inference", which is false. Pure imaginary
eigenvalues at a point say nothing about whether the pair CROSSES the axis
there (transversality) or what the cubic term does (the first Lyapunov
coefficient). Both are computed in `hopf_lyapunov.py`, per crossing rather than
per network, and only a crossing that passes both is called a Hopf bifurcation.

The escape statistic below is subject to the same caution and it is sharper than
it looks. `d_escape` is 1e-2 of the state scale and the perturbation is 1e-6 of
it, so "escaped" means the perturbation grew ten-thousandfold -- which a
SUPERCRITICAL Hopf produces just as reliably as a subcritical one, by growth
onto a limit cycle of finite amplitude. The escape rate therefore measures
departure from a small neighbourhood and CANNOT discriminate criticality. It was
read that way once; 120 of 145 crossings turn out to be supercritical despite
102 of 104 networks escaping.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from proteostasis import model as M  # noqa: E402
import fold_theorem as FT  # noqa: E402
import genericity as GEN  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "computed"


# ---------------------------------------------------------------------------
# the branch, without the ordering heuristic and with every point on {G = 0}
# ---------------------------------------------------------------------------


def _contourPaths(p: M.Params, u_star: float, a_star: float, n: int):
    UU, AA, GG, _, _, _ = GEN._fields(p, max(4.0 * u_star, 1e-2),
                                      max(6.0 * a_star, 1e-3), n)
    if not np.isfinite(GG).any():
        return None, (UU, AA)
    fig = plt.figure()
    try:
        cs = fig.gca().contour(UU, AA, GG, levels=[0.0])
        paths = [pp.vertices for cc in cs.collections for pp in cc.get_paths()] \
            if hasattr(cs, "collections") else [pp.vertices for pp in cs.get_paths()]
    except Exception:
        plt.close(fig)
        return None, (UU, AA)
    plt.close(fig)
    return (paths or None), (UU, AA)


def branchProfile(p: M.Params, u_star: float, a_star: float,
                  n: int = 150) -> Optional[Dict]:
    """the accessible branch as a profile, ordered along the curve.

    Returns the projected points with j, tr J and det J, the fold index, and the
    contiguous det J > 0 run that touches the fold. `None` on a trace failure,
    which is counted rather than dropped.
    """
    paths, (UU, AA) = _contourPaths(p, u_star, a_star, n)
    if paths is None:
        return None

    best, best_d = None, np.inf
    for v in paths:
        d = float(np.min(np.hypot(v[:, 0] - u_star, v[:, 1] - a_star)))
        if d < best_d:
            best, best_d = v, d
    if best is None or len(best) < 8:
        return None

    rows, n_proj_fail = [], 0
    for uu, aa in best:
        q = GEN._projectOntoNullcline(float(uu), float(aa), p)
        if q is None:
            n_proj_fail += 1
            rows.append(None)
            continue
        u, a = q
        try:
            J = M.jacobian(u, a, p)
            rows.append({
                "u": u, "a": a,
                "j": FT.removalR(u, a, p),
                "tr_J": float(J[0, 0] + J[1, 1]),
                "det_J": float(np.linalg.det(J)),
                "G": float(FT.aggregateG(u, a, p)),
            })
        except (M.ModelError, np.linalg.LinAlgError, OverflowError):
            n_proj_fail += 1
            rows.append(None)

    keep = [i for i, r in enumerate(rows) if r is not None]
    if len(keep) < 8:
        return None
    D = pd.DataFrame([rows[i] for i in keep]).reset_index(drop=True)

    k_fold = int(np.argmin(np.hypot(D["u"] - u_star, D["a"] - a_star)))

    # the accessible branch: the maximal contiguous det J > 0 run adjacent to the
    # fold. connectivity decides it; nothing guesses which end is low burden.
    pos = (D["det_J"] > 0.0).to_numpy()
    seeds = [k for k in (k_fold, k_fold - 1, k_fold + 1)
             if 0 <= k < len(pos) and pos[k]]
    if not seeds:
        return None
    s = seeds[0]
    lo = s
    while lo - 1 >= 0 and pos[lo - 1]:
        lo -= 1
    hi = s
    while hi + 1 < len(pos) and pos[hi + 1]:
        hi += 1
    B = D.iloc[lo:hi + 1].reset_index(drop=True)
    if len(B) < 5:
        return None

    return {
        "profile": D, "branch": B, "k_fold": k_fold,
        "lo": lo, "hi": hi,
        "n_proj_fail": n_proj_fail, "path_count": len(paths),
        "u_min": float(UU.min()), "a_min": float(AA.min()),
        "u_max": float(UU.max()), "a_max": float(AA.max()),
    }


def branchSummary(p: M.Params, u_star: float, a_star: float,
                  n: int = 150) -> Optional[Dict[str, float]]:
    """one row per network: does the refined branch still cross tr J = 0?"""
    out = branchProfile(p, u_star, a_star, n=n)
    if out is None:
        return None
    B, D = out["branch"], out["profile"]

    k_lo = int(B["j"].idxmin())          # the low-burden end of the branch
    k_hi = int(B["tr_J"].idxmax())       # where tr J is largest
    j_crit = float(FT.removalR(u_star, a_star, p))

    # the first tr J >= 0 point, walking from the low-burden end toward the fold
    order = B["j"].to_numpy().argsort()
    tr_sorted = B["tr_J"].to_numpy()[order]
    det_sorted = B["det_J"].to_numpy()[order]
    j_sorted = B["j"].to_numpy()[order]
    hit = np.nonzero(tr_sorted >= 0.0)[0]
    k_x = int(hit[0]) if hit.size else -1

    return {
        "tr_max": float(B["tr_J"].max()),
        "tr_at_low_j": float(B["tr_J"].iloc[k_lo]),
        "det_at_low_j": float(B["det_J"].iloc[k_lo]),
        "j_lo": float(B["j"].min()),
        "j_hi": float(B["j"].max()),
        "j_crit": j_crit,
        "n_points": int(len(B)),
        "n_cross": int((B["tr_J"] >= 0.0).sum()),
        "u_at_trmax": float(B["u"].iloc[k_hi]),
        "a_at_trmax": float(B["a"].iloc[k_hi]),
        "det_at_trmax": float(B["det_J"].iloc[k_hi]),
        "j_at_trmax": float(B["j"].iloc[k_hi]),
        "G_absmax": float(B["G"].abs().max()),
        # where the first crossing sits, as a fraction of the way to j_crit
        "j_at_first_cross": float(j_sorted[k_x]) if k_x >= 0 else np.nan,
        "det_at_first_cross": float(det_sorted[k_x]) if k_x >= 0 else np.nan,
        # is the branch pinned against the window edge?
        "at_u_edge": int(bool((B["u"] <= 1.001 * out["u_min"]).any())),
        "at_a_edge": int(bool((B["a"] <= 1.001 * out["a_min"]).any())),
        # the fold must be the j-maximum of the accessible branch. where it is
        # not, the contiguous det J > 0 run adjacent to the fold runs the OTHER
        # way -- a second singular point terminates it -- and the branch this
        # walk returns is not the one the influx reaches first. those networks
        # are ambiguous by multiplicity and are reported as such, not counted.
        "fold_is_j_max": int(bool(j_crit > 0 and B["j"].max() <= 1.001 * j_crit
                                  and B["j"].min() < 0.999 * j_crit)),
        # where the low-burden end of the branch stops, and the two cases are not
        # the same fact. hitting the window FLOOR (u = 1e-4, a = 1e-6) is the
        # branch reaching genuinely near-zero burden -- benign, and it is what
        # 325 of 325 load-grid branches do. hitting the CEILING means the curve
        # left through the high-burden side and the low-j part was never seen, so
        # "tr J at the low-burden end" is about the edge of the window rather
        # than about low burden. lumping the two reported 325 of 325 truncated.
        "lowj_at_floor": int(bool(
            B["u"].iloc[k_lo] <= 1.001 * out["u_min"]
            or B["a"].iloc[k_lo] <= 1.001 * out["a_min"])),
        "lowj_at_ceiling": int(bool(
            B["u"].iloc[k_lo] >= 0.999 * out["u_max"]
            or B["a"].iloc[k_lo] >= 0.999 * out["a_max"])),
        "frac_of_path": float(len(B) / max(len(D), 1)),
        "n_proj_fail": int(out["n_proj_fail"]),
        "path_count": int(out["path_count"]),
    }


# ---------------------------------------------------------------------------
# the check that shares no failure mode with the tracer
# ---------------------------------------------------------------------------


T_CAP = 2.0e4     # hard horizon; see integrationTest for why it is not optional


def integrationTest(p: M.Params, u0: float, a0: float,
                    rel: float = 1e-6, n_efold: float = 25.0,
                    n_t: int = 8000) -> Optional[Dict]:
    """perturb the equilibrium and integrate the full nonlinear system.

    The state is made an exact equilibrium by setting j := R(u0, a0): with G = 0
    there, du/dt = -G = 0 and da/dt = G = 0. If the linearisation says unstable
    the perturbation must leave, at a rate -- and, for a complex pair, with a
    period -- the Jacobian predicts. That comparison is what makes this a test
    rather than a restatement.

    THREE THINGS THE FIRST VERSION GOT WRONG, all of them measurement rather
    than physics, and all of them visible only because the whole population was
    run instead of the five worst cases:

      * a FIXED horizon. With `lambda_max` spanning 5.8e-05 to 3.6, t = 50 is
        170 e-folds for one network and 0.005 for another. The horizon has to
        scale with the growth rate it is measuring.
      * growth read at the ENDPOINT. A spiral oscillates while it grows, so a
        slowly growing one can end nearer the equilibrium than it started. Three
        networks were scored "did not grow" on that alone. Take the maximum over
        the trajectory.
      * a slope fitted over a PARTIAL period. For a complex pair log|delta|
        wobbles around the linear trend, so a window shorter than a few periods
        returns a slope that is mostly phase. Fit only when the window covers at
        least three periods, and report it as not evaluable otherwise rather
        than reporting a wrong number.

    The period is measured too, from the spacing of the peaks of |delta(t)|, and
    compared against 2*pi/omega. It is a second quantity the linearisation
    predicts and the nonlinear integration can refuse.
    """
    from scipy.integrate import solve_ivp

    try:
        j_eq = FT.removalR(u0, a0, p)
        q = p.with_(j=float(j_eq)).validate()
        res = float(np.hypot(*M.rhsVector([u0, a0], q)))
        J = M.jacobian(u0, a0, q)
    except (M.ModelError, np.linalg.LinAlgError, OverflowError, ValueError):
        return None
    scale = max(abs(u0), abs(a0), 1e-300)
    if res / scale > 1e-8:
        return None                      # not an equilibrium; report, do not test

    ev = np.linalg.eigvals(J)
    lam = float(np.max(ev.real))
    omega = float(np.max(np.abs(ev.imag)))
    is_complex = bool(omega > 1e-12)

    # the horizon carries a fixed number of e-folds of whatever rate is present,
    # and at least a few oscillation periods when there is a period. it is also
    # CAPPED: `25/lambda` is 4.3e+05 for the slowest network here, and a horizon
    # that long is not merely expensive, it is spent entirely past the escape,
    # where `solveFreePools` can run to its 2e5-iteration fixed-point cap on
    # every right-hand side call. where the cap binds, the network is reported
    # as too slow to observe rather than given a verdict.
    t_want = n_efold / max(abs(lam), 1e-12)
    if is_complex:
        t_want = max(t_want, 8.0 * 2.0 * np.pi / omega)
    t_end = float(min(t_want, T_CAP))
    horizon_short = bool(t_want > T_CAP)

    def f(t, x):
        if x[0] <= 0.0 or x[1] <= 0.0 or not np.all(np.isfinite(x)):
            return np.zeros(2)
        try:
            return M.rhsVector(x, q)
        except (M.ModelError, OverflowError):
            return np.zeros(2)

    # leaving the linear neighbourhood ANSWERS the question, so stop there.
    # integrating past it buys nothing and costs everything.
    d_escape = 1e-2 * scale

    def leave(t, x):
        return float(np.hypot(x[0] - u0, x[1] - a0)) - d_escape
    leave.terminal = True
    leave.direction = 1.0

    x0 = np.array([u0 * (1.0 + rel), a0 * (1.0 + rel)])
    ts = np.linspace(0.0, t_end, n_t)
    try:
        sol = solve_ivp(f, (0.0, t_end), x0, t_eval=ts, events=leave,
                        method="LSODA", rtol=1e-11, atol=1e-15 * scale)
    except Exception:
        return None
    # a fast escape leaves few sampled points, and dropping those would discard
    # exactly the most unstable networks -- the selection this whole check exists
    # to avoid. keep them; the exponent fit reports itself unevaluable instead.
    if not sol.success or sol.y.shape[1] < 5:
        return None

    d = np.hypot(sol.y[0] - u0, sol.y[1] - a0)
    t = sol.t
    d0 = float(d[0])
    if d0 <= 0.0 or not np.all(np.isfinite(d)):
        return None

    # the linear regime: everything up to escape, which is where the run stopped
    idx = np.nonzero((d < d_escape) & (d > 0.0))[0]
    slope, n_per, per_meas, estimator = np.nan, np.nan, np.nan, "none"
    if idx.size >= 20 and t[idx[-1]] > t[idx[0]]:
        tw, dw = t[idx], d[idx]
        span = float(tw[-1] - tw[0])
        n_per = span / (2.0 * np.pi / omega) if is_complex else np.inf

        # the peaks of |delta|. NOTE the factor of two: writing
        # delta(t) = e^(sigma t)(v cos wt + w sin wt) gives
        # |delta|^2 = e^(2 sigma t)(A + B cos 2wt + C sin 2wt), so the DISTANCE
        # oscillates at 2*omega and peaks every pi/omega, not every 2*pi/omega.
        # the first version predicted 2*pi/omega and the measurement came back at
        # exactly half it, twice; the measurement was right.
        pk = np.nonzero((dw[1:-1] > dw[:-2]) & (dw[1:-1] > dw[2:]))[0] + 1
        if is_complex and pk.size >= 3:
            per_meas = float(np.median(np.diff(tw[pk])))

        # TWO estimators, because the two eigenvalue types fail differently and a
        # single whole-window least squares is wrong for both. On a complex pair
        # log|delta| wobbles about the linear trend, so a window spanning a
        # non-integer number of periods is biased by the wobble; successive PEAKS
        # sit at identical phase, so fitting those is unbiased. On a real pair
        # there is no wobble but the subdominant mode contaminates the early
        # window; the late half is where the leading mode dominates. Measured on
        # 25 networks: whole-window least squares gave median relative error
        # 1.1e-01 on real pairs (1 of 13 within 5%), where the late half gives
        # 2.0e-02 (10 of 13) and peaks give 7.2e-04 on complex pairs (12 of 12).
        if is_complex and pk.size >= 3:
            slope = float(np.polyfit(tw[pk], np.log(dw[pk]), 1)[0])
            estimator = "peak"
        elif not is_complex:
            h = idx.size // 2
            if idx.size - h >= 8:
                slope = float(np.polyfit(tw[h:], np.log(dw[h:]), 1)[0])
                estimator = "late_half"

    return {
        "lambda_max": lam,
        "omega": omega,
        "complex_pair": is_complex,
        "slope": slope,
        "n_periods": float(n_per),
        "estimator": estimator,
        "period_measured": per_meas,
        "period_predicted": float(np.pi / omega) if is_complex else np.nan,
        # the run STOPS at escape, which caps d.max()/d0 near d_escape/d0 ~ 1e2,
        # so a 100x threshold could never fire. 10x is inside the cap and still
        # far outside anything a stable equilibrium does (controls sit at 1.0).
        "grew": bool(d.max() > 10.0 * d0),
        "ratio_max": float(d.max() / d0),
        "ratio_end": float(d[-1] / d0),
        "escaped": bool(sol.t_events is not None and len(sol.t_events[0]) > 0),
        "horizon_short": horizon_short,
        "t_end": t_end,
        "t_last": float(sol.t[-1]),
        "eq_residual": float(res / scale),
        "j_eq": float(j_eq),
        "n_fit": int(idx.size),
    }


# ---------------------------------------------------------------------------
# drivers
# ---------------------------------------------------------------------------


def parameterCorner(run: Path, S: pd.DataFrame, top: int = 6) -> pd.DataFrame:
    """where in the kinetic box the crossers sit, by per-parameter median ratio.

    A numerical artefact has no reason to concentrate on particular kinetic
    parameters. This is the check that the crossers are a REGION rather than a
    scatter, and it is the transferable half of the result -- the incidence rate
    is a property of the sampling box, the location is a property of the model.
    """
    c = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    c = c[pd.to_numeric(c["fold_burden"], errors="coerce").notna()]
    c = c.assign(name=[f"draw{i}" for i in c.index])
    T = S[S["traced"] == True]  # noqa: E712
    m = c.merge(T[["name", "tr_max", "fold_is_j_max"]], on="name")
    m["cross"] = (m["tr_max"] >= 0.0) & (m["fold_is_j_max"] == 1)
    cols = [x for x in c.columns if x.startswith("p_")]
    g = m.groupby("cross")[cols].median().T
    g.columns = ["no_cross", "cross"]
    g["ratio"] = g["cross"] / g["no_cross"]
    g["abs_log2"] = np.abs(np.log2(g["ratio"]))
    g = g.sort_values("abs_log2", ascending=False)
    g.attrs["n_cross"] = int(m["cross"].sum())
    g.attrs["n_total"] = int(len(m))
    return g.head(top)


def _summaryWorker(item):
    name, p, u, a, n = item
    try:
        s = branchSummary(p, u, a, n=n)
    except Exception:
        s = None
    if s is None:
        return {"name": name, "traced": False}
    s["name"] = name
    s["traced"] = True
    return s


def _integrationWorker(item):
    name, p, u, a = item
    try:
        r = integrationTest(p, u, a)
    except Exception:
        r = None
    if r is None:
        return {"name": name, "tested": False}
    r["name"] = name
    r["tested"] = True
    return r


def _pool(items, fn, workers: Optional[int] = None):
    import multiprocessing as mp
    if workers is None:
        workers = max(1, min(56, (os.cpu_count() or 2) - 8))
    if workers == 1 or len(items) < 2:
        return [fn(i) for i in items]
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=workers) as pool:
        return pool.map(fn, items, chunksize=4)


def main() -> int:
    run = FT.phase1RunDir()
    if not (run / "C" / "samples.tsv").exists():
        print("SKIP: no phase 1 run root")
        return 0
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    pops = {
        "load_grid": GEN.loadGridStates(run),
        "kinetic_box": GEN.kineticBoxStates(run),
    }

    for pop, states in pops.items():
        print(f"\n\nREFINED BRANCH  {pop}: {len(states)} networks")
        rows = _pool([(nm, p, u, a, 150) for nm, p, u, a in states],
                     _summaryWorker)
        S = pd.DataFrame(rows)
        S.to_csv(OUT_DIR / f"hopf_refined_{pop}.tsv", sep="\t", index=False)
        T = S[S["traced"] == True]  # noqa: E712
        print(f"  traced {len(T)}, trace failures {len(S) - len(T)}")
        if T.empty:
            continue
        cross = T[T["tr_max"] >= 0.0]
        print(f"  max tr J on the branch: median {T['tr_max'].median():+.3e}  "
              f"p99 {T['tr_max'].quantile(0.99):+.3e}  max {T['tr_max'].max():+.3e}")
        print(f"  networks with tr J >= 0 anywhere: {len(cross)} of {len(T)}")
        print(f"  tr J at the LOW-BURDEN end: median "
              f"{T['tr_at_low_j'].median():+.3e}  max {T['tr_at_low_j'].max():+.3e}"
              f"  >= 0 in {int((T['tr_at_low_j'] >= 0).sum())}")
        print(f"  det J at the LOW-BURDEN end: min {T['det_at_low_j'].min():.3e}"
              f"  <= 0 in {int((T['det_at_low_j'] <= 0).sum())}")
        print(f"  |G| on the projected branch: max {T['G_absmax'].max():.3e}")
        print(f"  branch touches the window's u floor in "
              f"{int(T['at_u_edge'].sum())}, a floor in {int(T['at_a_edge'].sum())}")
        print(f"  projection failures per network: max {T['n_proj_fail'].max():.0f}"
              f"  total {T['n_proj_fail'].sum():.0f}")

        print(f"  fold is NOT the j-maximum of the branch in "
              f"{int((T['fold_is_j_max'] == 0).sum())} of {len(T)} "
              "(ambiguous by multiplicity)")
        print(f"  low-burden end of the branch reaches the window FLOOR "
              f"(near-zero burden) in {int(T['lowj_at_floor'].sum())} of {len(T)}; "
              f"leaves through the CEILING (low-j part unseen) in "
              f"{int(T['lowj_at_ceiling'].sum())}")

        if cross.empty:
            continue
        clean = cross[cross["fold_is_j_max"] == 1]
        print(f"\n  the {len(cross)} crossers, {len(clean)} of them on branches "
              "whose j-maximum is the fold:")
        frac = clean["j_at_first_cross"] / clean["j_crit"]
        print(f"    j_first_cross / j_crit: min {frac.min():.3f}  "
              f"median {frac.median():.3f}  max {frac.max():.3f}")
        print(f"    det J at the first crossing: min "
              f"{clean['det_at_first_cross'].min():.3e}  "
              f"median {clean['det_at_first_cross'].median():.3e}  "
              f"<= 0 in {int((clean['det_at_first_cross'] <= 0).sum())}")
        print(f"    of these, branch left through the CEILING so the low-j part "
              f"is unseen: {int(clean['lowj_at_ceiling'].sum())} -- the crossing "
              "still precedes the fold, but where it begins is not determined")

        # --- the independent check, at the point of largest tr J --------------
        # run on the crossers AND on a control block of non-crossers. a test that
        # reports growth at a stable equilibrium reports nothing at all, and the
        # only way to know is to give it equilibria that must not grow.
        byname = {nm: (p, u, a) for nm, p, u, a in states}
        items = [("+" + r["name"], byname[r["name"]][0],
                  float(r["u_at_trmax"]), float(r["a_at_trmax"]))
                 for _, r in clean.iterrows()]
        ctrl = T[T["tr_max"] < 0.0]
        ctrl = ctrl.iloc[::max(1, len(ctrl) // 200)]
        items += [("-" + r["name"], byname[r["name"]][0],
                   float(r["u_at_trmax"]), float(r["a_at_trmax"]))
                  for _, r in ctrl.iterrows()]
        I = pd.DataFrame(_pool(items, _integrationWorker))
        I.to_csv(OUT_DIR / f"hopf_integration_{pop}.tsv", sep="\t", index=False)
        I["is_cross"] = I["name"].str.startswith("+")
        ok = I[I["tested"] == True]  # noqa: E712
        print(f"\n  INTEGRATION at the tr J maximum: {len(ok)} of {len(items)} "
              f"were exact equilibria and testable "
              f"({int(ok['is_cross'].sum())} crossers, "
              f"{int((~ok['is_cross']).sum())} stable controls)")
        for lab, sub in (("crossers", ok[ok["is_cross"]]),
                         ("controls", ok[~ok["is_cross"]])):
            if sub.empty:
                continue
            print(f"    {lab}: left the linear neighbourhood in "
                  f"{int(sub['escaped'].sum())} of {len(sub)}; grew >10x in "
                  f"{int(sub['grew'].sum())}; max |delta| ratio median "
                  f"{sub['ratio_max'].median():.4g}  max {sub['ratio_max'].max():.4g}")
            print(f"      horizon cap bound (too slow to observe in t={T_CAP:g}) "
                  f"in {int(sub['horizon_short'].sum())}")
        cr = ok[ok["is_cross"]]
        fit = cr[cr["slope"].notna()]
        if not fit.empty:
            rel = (fit["slope"] - fit["lambda_max"]).abs() \
                / fit["lambda_max"].abs().clip(lower=1e-12)
            print(f"    exponent vs max Re(lambda), {len(fit)} of {len(cr)} "
                  f"evaluable: median rel. diff {rel.median():.3e}  "
                  f"p90 {rel.quantile(0.90):.3e}  "
                  f"within 5% in {int((rel < 0.05).sum())}")
            for est, sub in fit.groupby("estimator"):
                e = (sub["slope"] - sub["lambda_max"]).abs() \
                    / sub["lambda_max"].abs().clip(lower=1e-12)
                print(f"      via {est:10s} n={len(e):3d}  median "
                      f"{e.median():.3e}  within 5% in {int((e < 0.05).sum())}")
        per = cr[cr["period_measured"].notna()]
        if not per.empty:
            pr = (per["period_measured"] - per["period_predicted"]).abs() \
                / per["period_predicted"]
            print(f"    oscillation period vs pi/omega, {len(per)} measurable: "
                  f"median rel. diff {pr.median():.3e}  "
                  f"within 5% in {int((pr < 0.05).sum())}")

        if pop == "kinetic_box":
            g = parameterCorner(run, S)
            g.to_csv(OUT_DIR / "hopf_parameter_corner.tsv", sep="\t")
            print(f"\n  WHERE the crossers sit: {g.attrs['n_cross']} of "
                  f"{g.attrs['n_total']} = "
                  f"{100 * g.attrs['n_cross'] / g.attrs['n_total']:.2f}%")
            for k, r in g.iterrows():
                d = "higher" if r["ratio"] > 1 else "lower"
                f = r["ratio"] if r["ratio"] > 1 else 1.0 / r["ratio"]
                print(f"    {k:12s} median {f:6.2f}x {d}  "
                      f"({r['no_cross']:.4g} -> {r['cross']:.4g})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
