"""C5a: does ONE cycle span the window, or is l1 being interpolated across it?

Section 7 currently says a stable limit cycle "is born as the influx enters the
window and is destroyed as it leaves". That is `l1 < 0` computed at the two
ENDPOINTS and an interpolation across the middle. It does not follow. Between the
two Hopf points the cycle could be destroyed and recreated, could collide with
another cycle at a fold of cycles, or could lose stability without the
equilibrium regaining it. Nothing computed so far looks inside the band.

THE TEST. For each resolved window network, take influx values spanning
`[j_H1, j_H2]`, locate the branch equilibrium at each, perturb it, and integrate
the full nonlinear system until the envelope SATURATES. Classify by that
envelope:

    cycle      the envelope stops changing between blocks
    fixed      collapses to a point (same equilibrium, or a different one)
    divergent  the envelope grows without bound, or the state leaves the orthant
    not_converged  neither, inside the block budget -- reported, not classified

A window spanned by one cycle gives `cycle` at every interior influx. Anything
else is a finding, and the point of running it is that both outcomes are
publishable and only one of them is what Section 7 currently asserts.

WHY THE ENVELOPE AND NOT THE ESCAPE FLAG. `hopf_check.integrationTest` stops at
`d_escape`, a ten-thousandfold amplification, which a stable cycle reaches just
as a divergence does (decision D058). This integration deliberately has NO
terminal event: it runs past that point, because the question is where the
trajectory ends up, not whether it left.

WHY STAGED, AND NOT ONE LONG SOLVE. The first version integrated once to
`min(60/sigma, 5e4)`. The cap bound 42 of 266 points, INCLUDING every point it
then called divergent, and the worst of them saw 3.1 e-folds of a requested 60.
`growth = e_late/e_mid > 3` cannot separate a slow spiral still on its way out
from a genuine divergence when the run stopped at a fifth of the horizon it asked
for -- and the truncation was worst exactly where the physics is slowest, at the
smallest `max(tr J)`, which is the closest approach to Hopf-pair annihilation and
the most interesting network in the set. Blocks fix both halves: the run
continues from the previous final state until the envelope stops moving, and each
block carries its own `t_eval` so 40 samples per period is held whatever the
total elapsed time. Non-convergence inside the budget is now its own outcome
rather than being absorbed into `divergent`.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import root

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import fold_theorem as FT  # noqa: E402
import genericity as GEN  # noqa: E402

COMPUTED = REPO_ROOT / "data" / "computed"
OUT = COMPUTED / "hopf_window_fate.json"
OUT_TSV = COMPUTED / "hopf_window_fate.tsv"

PERIODS_PER_BLOCK = 100.0   # samples/period is held at 4000/100 = 40
EFOLDS_PER_BLOCK_MAX = 200.0  # a runaway guard only -- see below
EFOLDS_NO_OSC = 20.0        # block length when there is no oscillation to count
EFOLDS_TARGET = 120.0       # twice the 60 the single-solve version asked for
N_BLOCKS_MIN, N_BLOCKS_MAX = 8, 240
SETTLE_TOL = 0.02           # envelope change between blocks, twice running
AMP_FLOOR = 1e-8            # below this the envelope is a fixed point
AMP_CEIL = 1e6              # above this it is a divergence

# WHY THE E-FOLD GUARD IS LOOSE. At 40 e-folds it bound the block length
# wherever `sigma/omega` was large, cutting the block to a median of 14
# oscillations -- and an envelope measured over a non-integer number of periods
# moves between blocks on phase alone. 110 of 266 points then failed a 1% settle
# test while carrying a SMALLER period dispersion than the points that passed
# (median cv 1.5e-3 against 8.1e-3), which is a settled orbit being called
# unsettled. The guard is now far enough out to be inactive in practice;
# divergence is caught inside the block by the orthant and finiteness checks,
# which is where it should have been caught in the first place.


def equilibriumAtInflux(p: M.Params, j: float, u0: float, a0: float):
    """solve {G = 0, R = j} near (u0, a0): the branch equilibrium at this load."""
    def residual(z):
        u, a = float(np.exp(z[0])), float(np.exp(z[1]))
        try:
            return [FT.aggregateG(u, a, p), FT.removalR(u, a, p) - j]
        except (M.ModelError, OverflowError):
            return [1e6, 1e6]

    s = root(residual, [np.log(max(u0, 1e-300)), np.log(max(a0, 1e-300))],
             method="hybr", options={"xtol": 1e-13})
    if not s.success:
        return None
    u, a = float(np.exp(s.x[0])), float(np.exp(s.x[1]))
    if not (np.isfinite(u) and np.isfinite(a) and u > 0 and a > 0):
        return None
    q = p.with_(j=j).validate()
    if float(np.hypot(*M.rhsVector([u, a], q))) / max(u, a, 1e-300) > 1e-8:
        return None
    return u, a


def _envelope(U, A):
    """peak-to-peak of each coordinate relative to its own mean, whichever larger."""
    du = (U.max() - U.min()) / max(abs(U.mean()), 1e-300)
    da = (A.max() - A.min()) / max(abs(A.mean()), 1e-300)
    return float(max(du, da))


def _peaks(A):
    return np.nonzero((A[1:-1] > A[:-2]) & (A[1:-1] >= A[2:]))[0] + 1


def _envelopeSnapped(U, A):
    """the envelope of the late trajectory over WHOLE oscillations.

    Anchoring both ends on a maximum of `a` makes the window an integer number
    of periods, so the block-to-block comparison measures amplitude and not the
    phase the block happened to end on. Falls back to the plain second half when
    there are too few peaks to snap to.
    """
    h0 = len(A) // 2
    pk = _peaks(A)
    pk = pk[pk >= h0]
    if pk.size >= 3:
        i, jx = int(pk[0]), int(pk[-1])
        return _envelope(U[i:jx + 1], A[i:jx + 1]), i, jx
    return _envelope(U[h0:], A[h0:]), h0, len(A) - 1


def fateAt(p: M.Params, j: float, u_eq: float, a_eq: float,
           rel: float = 1e-6) -> dict | None:
    """integrate in blocks until the envelope saturates, then classify it."""
    q = p.with_(j=float(j)).validate()
    try:
        J = M.jacobian(u_eq, a_eq, q)
    except (M.ModelError, np.linalg.LinAlgError):
        return None
    ev = np.linalg.eigvals(J)
    sigma = float(np.max(ev.real))
    omega = float(np.max(np.abs(ev.imag)))
    scale = max(abs(u_eq), abs(a_eq), 1e-300)
    s_abs = max(abs(sigma), 1e-12)

    # one block: 100 oscillations, so the envelope always spans whole periods
    if omega > 1e-12:
        t_blk = float(min(PERIODS_PER_BLOCK * 2.0 * np.pi / omega,
                          EFOLDS_PER_BLOCK_MAX / s_abs))
    else:
        t_blk = float(EFOLDS_NO_OSC / s_abs)
    n_blocks = int(np.ceil(EFOLDS_TARGET / s_abs / t_blk))
    n_blocks = int(np.clip(n_blocks, N_BLOCKS_MIN, N_BLOCKS_MAX))

    def f(t, x):
        if x[0] <= 0.0 or x[1] <= 0.0 or not np.all(np.isfinite(x)):
            return np.zeros(2)
        try:
            return M.rhsVector(x, q)
        except (M.ModelError, OverflowError, ValueError):
            return np.zeros(2)

    x = np.array([u_eq * (1.0 + rel), a_eq * (1.0 + rel)], float)
    envs, t_total, fate, last = [], 0.0, None, None
    n_settled = 0

    for _ in range(n_blocks):
        try:
            sol = solve_ivp(f, (0.0, t_blk), x, method="LSODA", rtol=1e-10,
                            atol=1e-14 * scale, dense_output=False,
                            t_eval=np.linspace(0.0, t_blk, 4000))
        except Exception:
            return None
        if not sol.success or sol.y.shape[1] < 100:
            return None
        U, A, T = sol.y[0], sol.y[1], sol.t
        t_total += t_blk

        if not (np.all(np.isfinite(U)) and np.all(np.isfinite(A))):
            fate, last = "divergent", (U, A, T)
            break
        if U.min() <= 0.0 or A.min() <= 0.0:
            fate, last = "divergent", (U, A, T)
            break

        # the envelope of the SECOND half of the block, snapped to whole
        # oscillations: the first half still carries the hand-off transient from
        # the previous block's final state
        e, i0, i1 = _envelopeSnapped(U, A)
        h = slice(i0, i1 + 1)
        envs.append(e)
        x = np.array([U[-1], A[-1]], float)
        last = (U, A, T)

        if e > AMP_CEIL:
            fate = "divergent"
            break
        if e < AMP_FLOOR:
            d0 = float(np.hypot(U[-1] - u_eq, A[-1] - a_eq)) / scale
            fate = "fixed_same" if d0 < 1e-4 else "fixed_other"
            break
        if len(envs) >= 2:
            prev = envs[-2]
            if abs(e - prev) / max(e, prev, 1e-300) < SETTLE_TOL:
                n_settled += 1
                if n_settled >= 2:      # two consecutive quiet blocks
                    fate = "cycle"
                    break
            else:
                n_settled = 0

    if fate is None:
        fate = "not_converged"
    if last is None:
        return None

    U, A, T = last
    e_last, i0, i1 = _envelopeSnapped(U, A)
    h = slice(i0, i1 + 1)
    e_late = envs[-1] if envs else e_last
    e_mid = envs[-2] if len(envs) >= 2 else float("nan")
    growth = float(e_late / e_mid) if len(envs) >= 2 else float("nan")
    first = slice(0, len(T) // 2)
    drift = (abs(U[h].mean() - U[first].mean())
             / max(abs(U[first].mean()), 1e-300))

    # periodicity: successive maxima of a(t) in the final block. a settled cycle
    # has both a stable period and a stable peak height; a slow transient has
    # neither, and the escape flag cannot tell them apart (D058).
    aa, tt = A[h], T[h]
    pk = _peaks(aa)
    period = period_cv = peak_cv = np.nan
    n_peaks = int(pk.size)
    if n_peaks >= 4:
        gaps = np.diff(tt[pk])
        period = float(np.median(gaps))
        period_cv = float(np.std(gaps) / max(abs(np.mean(gaps)), 1e-300))
        hts = aa[pk]
        peak_cv = float(np.std(hts) / max(abs(np.mean(hts)), 1e-300))
    period_pred = float(2.0 * np.pi / omega) if omega > 1e-12 else np.nan

    return {
        "fate": fate, "j": float(j), "sigma": sigma, "omega": omega,
        "amp_late": float(e_late), "amp_mid": float(e_mid),
        "growth": growth, "drift": float(drift),
        "u_eq": u_eq, "a_eq": a_eq,
        # the horizon is now reported, because the version that hid it inside a
        # constant classified 42 truncated points as if they had converged
        "t_end": float(t_total), "n_blocks_run": int(len(envs)),
        "n_blocks_budget": int(n_blocks), "t_block": float(t_blk),
        "efolds": float(t_total * s_abs),
        "converged": bool(fate != "not_converged"),
        "u_min": float(U[h].min()), "u_max": float(U[h].max()),
        # section 9's observable, computed rather than predicted: does the mean
        # burden over the cycle stay near the equilibrium while the amplitude
        # grows, or does the cycle sit at a higher mean burden?
        "a_amp_rel": float((A[h].max() - A[h].min()) / max(a_eq, 1e-300)),
        "a_mean_over_eq": float(A[h].mean() / max(a_eq, 1e-300)),
        "u_mean_over_eq": float(U[h].mean() / max(u_eq, 1e-300)),
        "n_peaks": n_peaks, "period": period, "period_cv": period_cv,
        "peak_cv": peak_cv, "period_predicted": period_pred,
    }


def _worker(item):
    """ONE (network, influx) point. Flattened so the pool sees 266 units of work
    rather than 38; with 7 sequential integrations per unit the wall time was
    the slowest chain, not the slowest point."""
    name, p, u_seed, a_seed, j, fr, j_crit = item
    eq = equilibriumAtInflux(p, j, u_seed, a_seed)
    if eq is None:
        return {"name": name, "frac": float(fr), "j": float(j),
                "fate": "no_equilibrium"}
    r = fateAt(p, j, eq[0], eq[1])
    if r is None:
        return {"name": name, "frac": float(fr), "j": float(j),
                "fate": "not_evaluable"}
    r.update({"name": name, "frac": float(fr),
              "j_over_jcrit": float(j / j_crit)})
    return r


def _seedWorker(item):
    """branch points nearest each target influx, as seeds for the equilibrium solve.

    Seeding every influx from the FOLD state failed at 145 of 266 points: the
    window can sit far below the fold and hybr was started in the wrong place.
    The branch trace already knows where the equilibrium is at each load.
    """
    import hopf_check as HC
    name, p, u_f, a_f, jH1, jH2, j_crit = item
    try:
        out = HC.branchProfile(p, u_f, a_f, n=150)
    except Exception:
        out = None
    fracs = np.linspace(0.08, 0.92, 7)
    items = []
    if out is None:
        for fr in fracs:
            items.append((name, p, u_f, a_f, jH1 + fr * (jH2 - jH1), fr, j_crit))
        return items
    B = out["branch"]
    jb = B["j"].to_numpy(float)
    for fr in fracs:
        j = jH1 + fr * (jH2 - jH1)
        k = int(np.argmin(np.abs(jb - j)))
        items.append((name, p, float(B["u"].iloc[k]), float(B["a"].iloc[k]),
                      float(j), float(fr), j_crit))
    return items


def run(workers: int | None = None, max_networks: int | None = None) -> dict:
    Z = pd.read_csv(COMPUTED / "hopf_zero_counts.tsv", sep="\t")
    W = Z[(Z["group"] == "window") & (Z["parity_ok"] == True)]  # noqa: E712
    W = W[W["j_H2"].notna()]
    if max_networks:
        W = W.head(max_networks)

    run_dir = FT.phase1RunDir()
    c = pd.read_csv(run_dir / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    c = c[pd.to_numeric(c["fold_burden"], errors="coerce").notna()]
    byname = {}
    for i, r in c.iterrows():
        byname[f"draw{i}"] = r

    items = []
    for _, w in W.iterrows():
        r = byname.get(w["name"])
        if r is None:
            continue
        try:
            p = FT.paramsFromSampleRow(r)
            u0, a0 = FT.foldStateFromSampleRow(r)
        except Exception:
            continue
        s = GEN.polishFold(p, u0, a0) or FT.foldSolve(p)
        if s is None:
            continue
        items.append((w["name"], p, float(s[1]), float(s[2]),
                      float(w["j_H1"]), float(w["j_H2"]), float(w["j_crit"])))

    if workers is None:
        workers = max(1, min(56, (os.cpu_count() or 2) - 8))
    with mp.get_context("fork").Pool(processes=workers) as pool:
        seeded = pool.map(_seedWorker, items, chunksize=1)
        flat = [x for sub in seeded for x in sub]
        rows = pool.map(_worker, flat, chunksize=1)
    D = pd.DataFrame(rows)
    D.to_csv(OUT_TSV, sep="\t", index=False)

    counts = D["fate"].value_counts().to_dict()
    _EVAL = ("cycle", "divergent", "fixed_same", "fixed_other")

    def verdict(s):
        """classify a network by its EVALUABLE points only.

        The first run called a network "no_cycle" when every point failed to
        yield an equilibrium, which reported a method failure as a result about
        the model. Not-evaluable is now its own verdict and is counted, and so
        is a point whose envelope had not settled inside the block budget.
        """
        e = s[s.isin(_EVAL)]
        if e.empty:
            return "none_evaluable"
        if (e == "cycle").all():
            return "all_evaluable_cycle"
        return "no_cycle" if (e == "cycle").sum() == 0 else "mixed"

    per_net = D.groupby("name")["fate"].apply(verdict)

    out = {
        "n_networks": int(len(items)),
        "n_points": int(len(D)),
        "fate_counts": {k: int(v) for k, v in counts.items()},
        "networks_all_evaluable_points_cycle": int(
            (per_net == "all_evaluable_cycle").sum()),
        "networks_mixed": int((per_net == "mixed").sum()),
        "networks_no_cycle": int((per_net == "no_cycle").sum()),
        "networks_none_evaluable": int((per_net == "none_evaluable").sum()),
        "points_per_network_evaluable_median": float(
            D[D["fate"].isin(_EVAL)].groupby("name").size().median()),
    }
    if "converged" in D.columns:
        out["n_not_converged"] = int((D["fate"] == "not_converged").sum())
        fin = D[D["efolds"].notna()]
        out["efolds_min"] = float(fin["efolds"].min())
        out["efolds_median"] = float(fin["efolds"].median())
        out["n_at_block_budget"] = int(
            (fin["n_blocks_run"] >= fin["n_blocks_budget"]).sum())
    ev = D[D["fate"].isin(_EVAL)]
    if not ev.empty:
        out["frac_cycle_over_evaluable"] = float((ev["fate"] == "cycle").mean())
        out["n_evaluable"] = int(len(ev))
    return out


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-networks", type=int, default=None)
    args = ap.parse_args()
    o = run(max_networks=args.max_networks)
    OUT.write_text(json.dumps(o, indent=2) + "\n")
    print(json.dumps(o, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
