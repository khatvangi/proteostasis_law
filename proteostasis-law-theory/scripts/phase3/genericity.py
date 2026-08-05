"""Block 1: the genericity conditions for the converse, and Hopf exclusion.

WHAT THEOREM 1 ACTUALLY PROVES
Saddle-node => grad R parallel to grad G. The converse needs hypotheses, because
`det J = 0` is also consistent with a transcritical or pitchfork bifurcation, a
cusp, a double-zero eigenvalue, or a curve of equilibria. The word "exactly" in
part (3) claims an iff that was not established. This module verifies the standard
genericity conditions at every recorded fold state, so the converse can be stated
with its hypotheses rather than asserted.

  (G1) constraint regularity      |grad G| != 0
  (G2) simple zero eigenvalue     tr J != 0   (2D: eigenvalues are {0, tr J})
  (G3) nondegeneracy              d2R/ds2 != 0 along {G = 0} at the critical point
  (G4) parameter transversality   w . dF/dj != 0, w the left null vector of J
  (G5) n-state constraint rank    the non-influx constraint gradients are independent

WHY (G4) IS THE ONE THAT COULD HAVE FAILED
The influx enters `du/dt` only, so dF/dj = (1, 0). Transversality is therefore
exactly `w_1 != 0`, and nothing in the model forces that.

HOPF EXCLUSION
Locating a fold does not show it is the FIRST loss of stability. A Hopf
bifurcation (tr J = 0 while det J > 0) could precede it on the stable low-burden
branch. This module traces that branch from low burden to the fold and reports
max(tr J) over it. The trace follows the nullcline as a CONTOUR, because
root-finding in `a` at fixed `u` loses the curve at its turning point and the fold
lies past that turn (D034).

MULTIPLICITY
Solving {G = 0, det J = 0} gives candidate singular points. If a network has more
than one, the equations alone do not say which lies on the accessible branch. The
count is reported per network so Corollary 1's "no continuation sweep" can be
scoped to what the ensemble actually shows.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.optimize import brentq

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from proteostasis import model as M  # noqa: E402
import fold_theorem as FT  # noqa: E402


# ---------------------------------------------------------------------------
# 1a. genericity at a recorded fold state
# ---------------------------------------------------------------------------


def _onNullcline(u: float, a_seed: float, p: M.Params,
                 tol: float = 1e-14, n_iter: int = 40) -> Optional[float]:
    """solve G(u, .) = 0 near `a_seed` by Newton, staying on the SAME branch.

    A bracket does not work here. Near the fold G(u, .) has two roots within a
    few percent of `a_seed` -- it runs positive, negative, then positive again --
    so any bracket wide enough to be robust contains both and brentq either fails
    or lands on the wrong branch. Newton from the seed converges to the local one.
    """
    a = float(a_seed)
    for _ in range(n_iter):
        try:
            g = FT.aggregateG(u, a, p)
            h = 1e-7 * max(a, 1e-12)
            ga = (FT.aggregateG(u, a + h, p) - FT.aggregateG(u, a - h, p)) / (2 * h)
        except (M.ModelError, OverflowError):
            return None
        if not (np.isfinite(g) and np.isfinite(ga)) or ga == 0.0:
            return None
        step = g / ga
        a_new = a - step
        if a_new <= 0.0:
            a_new = 0.5 * a
        if abs(a_new - a) <= tol * max(a, 1e-12):
            return float(a_new)
        a = a_new
    return float(a) if abs(FT.aggregateG(u, a, p)) < 1e-10 else None


def _projectOntoNullcline(u: float, a: float, p: M.Params,
                          tol: float = 1e-14, n_iter: int = 60):
    """Newton back onto {G = 0} along grad G, from a point just off the curve.

    Correcting along the GRADIENT rather than along `a` is what makes this work
    through a vertical tangent. Fixing `u` and solving in `a` failed one-sided at
    12 of the 325 load-grid folds, which is D034's failure in a new guise.
    """
    for _ in range(n_iter):
        try:
            g = FT.aggregateG(u, a, p)
            gu, ga = FT._centralGradient(FT.aggregateG, u, a, p)
        except (M.ModelError, OverflowError):
            return None
        n2 = gu * gu + ga * ga
        if not np.isfinite(g) or n2 <= 0.0:
            return None
        du, da = -g * gu / n2, -g * ga / n2
        u, a = u + du, a + da
        if a <= 0.0 or u <= 0.0:
            return None
        if abs(du) + abs(da) <= tol * (abs(u) + abs(a)):
            return float(u), float(a)
    try:
        return (float(u), float(a)) if abs(FT.aggregateG(u, a, p)) < 1e-10 else None
    except (M.ModelError, OverflowError):
        return None


def secondDerivativeAlongNullcline(u: float, a: float, p: M.Params,
                                   h_rel: float = 1e-4) -> Optional[float]:
    """d2R/ds2 along {G = 0} at the critical point, by ARCLENGTH.

    The fold is where dR/ds = 0 on the constraint, so this second derivative is
    what separates a genuine turning point from a degenerate one. Stepping along
    the tangent and correcting along the gradient parametrises by arclength,
    which stays valid where the nullcline runs vertical.
    """
    try:
        gu, ga = FT._centralGradient(FT.aggregateG, u, a, p)
    except (M.ModelError, OverflowError):
        return None
    n = float(np.hypot(gu, ga))
    if n <= 0.0:
        return None
    tu, ta = -ga / n, gu / n
    h = h_rel * max(float(np.hypot(u, a)), 1e-12)

    pts = []
    for k in (-1, 0, 1):
        if k == 0:
            pts.append((u, a))
            continue
        q = _projectOntoNullcline(u + k * h * tu, a + k * h * ta, p)
        if q is None:
            return None
        pts.append(q)

    try:
        vals = [FT.removalR(x, y, p) for x, y in pts]
    except (M.ModelError, OverflowError):
        return None
    s = [0.0]
    for k in (1, 2):
        s.append(s[-1] + float(np.hypot(pts[k][0] - pts[k - 1][0],
                                        pts[k][1] - pts[k - 1][1])))
    (s0, s1, s2), (r0, r1, r2) = s, vals
    if min(s1 - s0, s2 - s1) <= 0.0:
        return None
    d2 = 2.0 * (r0 * (s2 - s1) - r1 * (s2 - s0) + r2 * (s1 - s0)) \
        / ((s1 - s0) * (s2 - s1) * (s2 - s0))
    return float(d2)


def genericityAt(u: float, a: float, p: M.Params) -> Optional[Dict[str, float]]:
    """the five conditions, evaluated at one recorded fold state."""
    try:
        J = M.jacobian(u, a, p)
    except (M.ModelError, np.linalg.LinAlgError):
        return None
    Gu, Ga = FT._centralGradient(FT.aggregateG, u, a, p)
    Ru, Ra = FT._centralGradient(FT.removalR, u, a, p)
    grad_G = float(np.hypot(Gu, Ga))
    grad_R = float(np.hypot(Ru, Ra))
    tr = float(J[0, 0] + J[1, 1])
    det = float(np.linalg.det(J))

    # left null vector of J: w^T J = 0, i.e. the null vector of J^T
    _, sv, vt = np.linalg.svd(J.T)
    w = vt[-1, :]
    # dF/dj = (1, 0) exactly, because j enters du/dt and nothing else
    transversality = float(abs(w[0]) / max(float(np.linalg.norm(w)), 1e-300))

    d2 = secondDerivativeAlongNullcline(u, a, p)
    return {
        "u": u, "a": a,
        "grad_G": grad_G, "grad_R": grad_R,
        "tr_J": tr, "det_J": det,
        "sv_ratio": float(sv[-1] / max(sv[0], 1e-300)),
        "transversality": transversality,
        "d2R_ds2": np.nan if d2 is None else abs(d2),
        "d2R_signed": np.nan if d2 is None else d2,
    }


# ---------------------------------------------------------------------------
# 1b. the stable low-burden branch, traced through the nullcline's turn
# ---------------------------------------------------------------------------


def _fields(p: M.Params, u_hi: float, a_hi: float, n: int):
    us = np.geomspace(1e-4, u_hi, n)
    as_ = np.geomspace(1e-6, a_hi, n)
    UU, AA = np.meshgrid(us, as_)
    GG = np.full_like(UU, np.nan)
    RR = np.full_like(UU, np.nan)
    TR = np.full_like(UU, np.nan)
    DT = np.full_like(UU, np.nan)
    for i in range(UU.shape[0]):
        for k in range(UU.shape[1]):
            u, a = float(UU[i, k]), float(AA[i, k])
            try:
                GG[i, k] = FT.aggregateG(u, a, p)
                RR[i, k] = FT.removalR(u, a, p)
                J = M.jacobian(u, a, p)
                TR[i, k] = J[0, 0] + J[1, 1]
                DT[i, k] = np.linalg.det(J)
            except (M.ModelError, np.linalg.LinAlgError, OverflowError):
                pass
    return UU, AA, GG, RR, TR, DT


def branchToFold(p: M.Params, u_star: float, a_star: float,
                 n: int = 150) -> Optional[Dict[str, float]]:
    """walk {G = 0} from low burden to the fold; report tr J over that segment.

    Returns None when the contour cannot be traced, which is reported as a trace
    failure rather than quietly dropped.
    """
    UU, AA, GG, RR, TR, DT = _fields(p, max(4.0 * u_star, 1e-2),
                                     max(6.0 * a_star, 1e-3), n)
    if not np.isfinite(GG).any():
        return None
    fig = plt.figure()
    try:
        cs = fig.gca().contour(UU, AA, GG, levels=[0.0])
        paths = [pp.vertices for cc in cs.collections for pp in cc.get_paths()] \
            if hasattr(cs, "collections") else [pp.vertices for pp in cs.get_paths()]
    except Exception:
        plt.close(fig)
        return None
    plt.close(fig)
    if not paths:
        return None

    # the segment that passes closest to the solved fold is the accessible branch
    best, best_d = None, np.inf
    for v in paths:
        d = np.min(np.hypot(v[:, 0] - u_star, v[:, 1] - a_star))
        if d < best_d:
            best, best_d = v, d
    if best is None or len(best) < 8:
        return None

    # order the polyline from low burden, then cut at the solved fold
    v = best if best[0, 0] + best[0, 1] < best[-1, 0] + best[-1, 1] else best[::-1]
    k_fold = int(np.argmin(np.hypot(v[:, 0] - u_star, v[:, 1] - a_star)))
    seg = v[:k_fold + 1]
    if len(seg) < 5:
        return None

    rows = []
    for u, a in seg:
        try:
            J = M.jacobian(float(u), float(a), p)
            rows.append((FT.removalR(float(u), float(a), p),
                         float(J[0, 0] + J[1, 1]), float(np.linalg.det(J))))
        except (M.ModelError, np.linalg.LinAlgError, OverflowError):
            continue
    if len(rows) < 5:
        return None
    D = pd.DataFrame(rows, columns=["j", "tr_J", "det_J"])
    # the accessible branch is the stable side: keep points up to the first
    # det J <= 0, which is the fold itself
    stable = D[D["det_J"] > 0.0]
    if stable.empty:
        return None
    # every det J sign change on the WHOLE nullcline, from the same field pass
    n_singular = 0
    for vv in paths:
        dd = []
        for u, a in vv:
            try:
                dd.append(float(np.linalg.det(M.jacobian(float(u), float(a), p))))
            except (M.ModelError, np.linalg.LinAlgError, OverflowError):
                dd.append(np.nan)
        arr = np.asarray(dd)
        good = np.isfinite(arr)
        if good.sum() >= 3:
            sg = np.sign(arr[good])
            n_singular += int((sg[:-1] * sg[1:] < 0).sum())

    return {
        "n_singular": n_singular,
        "n_points": int(len(stable)),
        "tr_max": float(stable["tr_J"].max()),
        "tr_at_fold": float(D["tr_J"].iloc[-1]),
        "tr_crossings": int((stable["tr_J"] >= 0.0).sum()),
        "j_max_on_branch": float(stable["j"].max()),
        "det_min_on_branch": float(stable["det_J"].min()),
        "path_count": len(paths),
    }


def singularPointCount(p: M.Params, u_star: float, a_star: float,
                       n: int = 150) -> Optional[int]:
    """how many points of {G = 0} in the plotted window have det J = 0."""
    UU, AA, GG, RR, TR, DT = _fields(p, max(4.0 * u_star, 1e-2),
                                     max(6.0 * a_star, 1e-3), n)
    fig = plt.figure()
    try:
        cs = fig.gca().contour(UU, AA, GG, levels=[0.0])
        paths = [pp.vertices for cc in cs.collections for pp in cc.get_paths()] \
            if hasattr(cs, "collections") else [pp.vertices for pp in cs.get_paths()]
    except Exception:
        plt.close(fig)
        return None
    plt.close(fig)
    total = 0
    for v in paths:
        dets = []
        for u, a in v:
            try:
                dets.append(float(np.linalg.det(M.jacobian(float(u), float(a), p))))
            except (M.ModelError, np.linalg.LinAlgError, OverflowError):
                dets.append(np.nan)
        d = np.asarray(dets)
        ok = np.isfinite(d)
        if ok.sum() < 3:
            continue
        s = np.sign(d[ok])
        total += int((s[:-1] * s[1:] < 0).sum())
    return total


# ---------------------------------------------------------------------------
# ensemble drivers
# ---------------------------------------------------------------------------

OUT_DIR = REPO_ROOT / "data" / "computed"

# a condition is "violated" when it is not bounded away from zero by this much.
# stated rather than tuned: the identity residual floor is ~1e-9, so anything
# above 1e-6 is unambiguously nonzero and anything below is reported explicitly.
TOL = 1e-6


def loadGridStates(run: Path) -> List[Tuple[str, M.Params, float, float]]:
    b = pd.read_csv(run / "B" / "fold_boundary.tsv", sep="\t")
    b = b[b["found"] == True]  # noqa: E712
    base = M.Params()
    out = []
    for _, r in b.iterrows():
        try:
            p = M.allocationParams(base.with_(nu=float(r["nu"])),
                                   float(r["chi"])).validate()
        except M.ModelError:
            continue
        out.append((f"nu={r['nu']:.4g},chi={r['chi']:.3g}", p,
                    float(r["fold_u"]), float(r["fold_a"])))
    return out


def polishFold(p: M.Params, u0: float, a0: float,
               tol: float = 1e-12, n_iter: int = 60):
    """Newton onto {G = 0, det J = 0} from the RECORDED state as the seed.

    `foldSolve` seeds from its own coarse scan over `lowerNullclineA`, which is a
    first-root heuristic that does not identify the branch the fold lives on. It
    failed to converge for 718 of the 2884 kinetic-box networks -- a 25% loss, and
    one that could correlate with the very geometry a genericity test is about, so
    reporting "no violations over the survivors" would have been a selected
    result. Seeding from the recorded state instead converges from the right
    branch, because experiment C's continuation already put it there.
    """
    x = np.array([np.log(max(u0, 1e-300)), np.log(max(a0, 1e-300))], float)

    def F(z):
        u, a = float(np.exp(z[0])), float(np.exp(z[1]))
        return np.array([FT.aggregateG(u, a, p),
                         float(np.linalg.det(M.jacobian(u, a, p)))], float)

    for _ in range(n_iter):
        try:
            f = F(x)
        except (M.ModelError, np.linalg.LinAlgError, OverflowError):
            return None
        if not np.all(np.isfinite(f)):
            return None
        Jm = np.empty((2, 2))
        for k in range(2):
            h = 1e-7 * max(abs(x[k]), 1.0)
            xp, xm = x.copy(), x.copy()
            xp[k] += h
            xm[k] -= h
            try:
                Jm[:, k] = (F(xp) - F(xm)) / (2 * h)
            except (M.ModelError, np.linalg.LinAlgError, OverflowError):
                return None
        try:
            step = np.linalg.solve(Jm, f)
        except np.linalg.LinAlgError:
            return None
        if not np.all(np.isfinite(step)):
            return None
        # damped, so a bad first step cannot throw the iterate out of the domain
        step = np.clip(step, -0.5, 0.5)
        x = x - step
        if np.max(np.abs(step)) < tol:
            break
    u, a = float(np.exp(x[0])), float(np.exp(x[1]))
    if not (np.isfinite(u) and np.isfinite(a) and u > 0.0 and a > 0.0):
        return None
    try:
        if abs(FT.aggregateG(u, a, p)) > 1e-8:
            return None
        if abs(float(np.linalg.det(M.jacobian(u, a, p)))) > 1e-8:
            return None
    except (M.ModelError, np.linalg.LinAlgError, OverflowError):
        return None
    return FT.removalR(u, a, p), u, a


def kineticBoxStates(run: Path, resolve: bool = True
                     ) -> List[Tuple[str, M.Params, float, float]]:
    """the kinetic box, with the fold RE-SOLVED rather than read off the record.

    The recorded experiment-C fold states are not tight enough to test genericity
    on. |det J| at them has median 2.8e-05 but p99 2.4 and max 52, against a
    load-grid maximum of 8.9e-05 -- the top percentile are not folds at all, and
    "tr J = 0 at a fold" evaluated there measures nothing. `foldSolve` returns the
    state that satisfies {G = 0, det J = 0} for the same parameters, which is what
    the conditions are about. Networks whose solve does not converge are counted
    and reported, never dropped silently.
    """
    c = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    c = c[pd.to_numeric(c["fold_burden"], errors="coerce").notna()]
    out, failed = [], 0
    for i, r in c.iterrows():
        try:
            p = FT.paramsFromSampleRow(r)
            u, a = FT.foldStateFromSampleRow(r)
        except (M.ModelError, ValueError, KeyError):
            failed += 1
            continue
        if resolve:
            s = polishFold(p, u, a) or FT.foldSolve(p)
            if s is None:
                failed += 1
                continue
            _, u, a = s
        out.append((f"draw{i}", p, float(u), float(a)))
    if failed:
        print(f"  kinetic box: {failed} of {len(c)} networks have no solvable "
              "fold state and carry no genericity result; reported, not screened")
    return out


def runGenericity(states) -> pd.DataFrame:
    rows = []
    for name, p, u, a in states:
        g = genericityAt(u, a, p)
        if g is None:
            rows.append({"name": name, "ok": False})
            continue
        g["name"] = name
        g["ok"] = True
        rows.append(g)
    return pd.DataFrame(rows)


def _branchWorker(item):
    """one network, for the process pool. Top level so it pickles."""
    name, p, u, a, n = item
    try:
        b = branchToFold(p, u, a, n=n)
    except Exception:
        b = None
    if b is None:
        return {"name": name, "traced": False}
    b["name"] = name
    b["traced"] = True
    return b


def runBranches(states, n: int = 150, workers: Optional[int] = None
                ) -> pd.DataFrame:
    """one field pass per network, fanned out across cores.

    Each network is an independent 150x150 field evaluation plus a contour walk,
    so this is embarrassingly parallel and was needlessly serial: 2767 networks at
    ~2.2 s each is 100 minutes on one core and about two on this machine. `map`
    preserves input order, so the output is identical to the serial version --
    checked, not assumed, by `--check-serial`.
    """
    import multiprocessing as mp

    if workers is None:
        workers = max(1, min(56, (os.cpu_count() or 2) - 8))
    items = [(name, p, u, a, n) for name, p, u, a in states]
    if workers == 1:
        return pd.DataFrame([_branchWorker(i) for i in items])
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=workers) as pool:
        rows = pool.map(_branchWorker, items, chunksize=4)
    return pd.DataFrame(rows)


def _report(D: pd.DataFrame, col: str, label: str, *, want: str) -> str:
    v = D[col].dropna()
    if v.empty:
        return f"  {label:34s} no finite values"
    if want == "away_from_zero":
        bad = int((v.abs() <= TOL).sum())
        return (f"  {label:34s} n={len(v):5d}  median {v.abs().median():.3e}  "
                f"min |.| {v.abs().min():.3e}  violations {bad}")
    if want == "negative":
        bad = int((v >= 0.0).sum())
        return (f"  {label:34s} n={len(v):5d}  median {v.median():+.3e}  "
                f"max {v.max():+.3e}  crossings {bad}")
    return f"  {label:34s} n={len(v)}"


def main() -> int:
    run = FT.phase1RunDir()
    if not (run / "B" / "fold_boundary.tsv").exists():
        print("SKIP: no phase 1 run root")
        return 0
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    pops = {"load_grid": loadGridStates(run), "kinetic_box": kineticBoxStates(run)}

    print("BLOCK 1a  genericity conditions at the recorded fold states")
    for pop, states in pops.items():
        D = runGenericity(states)
        D.to_csv(OUT_DIR / f"genericity_{pop}.tsv", sep="\t", index=False)
        ok = D[D["ok"] == True]  # noqa: E712
        print(f"\n  {pop}: {len(states)} states, {len(ok)} evaluated, "
              f"{len(D) - len(ok)} not evaluable")
        print(_report(ok, "grad_G", "(G1) |grad G|", want="away_from_zero"))
        print(_report(ok, "tr_J", "(G2) tr J", want="away_from_zero"))
        print(_report(ok, "d2R_signed", "(G3) d2R/ds2 along {G=0}",
                      want="away_from_zero"))
        print(_report(ok, "transversality", "(G4) |w.dF/dj| / |w|",
                      want="away_from_zero"))
        n_nan = int(ok["d2R_signed"].isna().sum())
        if n_nan:
            print(f"  (G3) not evaluable at {n_nan} states "
                  "-- reported, not screened")

    print("\n\nBLOCK 1b  tr J on the stable low-burden branch, and multiplicity")
    for pop, states in pops.items():
        D = runBranches(states)
        D.to_csv(OUT_DIR / f"branch_{pop}.tsv", sep="\t", index=False)
        tr = D[D["traced"] == True]  # noqa: E712
        print(f"\n  {pop}: {len(states)} networks, {len(tr)} traced, "
              f"{len(D) - len(tr)} trace failures")
        if tr.empty:
            continue
        print(_report(tr, "tr_max", "max tr J on the branch", want="negative"))
        print(f"  {'networks with tr J >= 0 anywhere':34s} "
              f"{int((tr['tr_max'] >= 0.0).sum())} of {len(tr)}")
        ns = tr["n_singular"].dropna()
        if not ns.empty:
            print(f"  {'singular points per network':34s} "
                  f"median {ns.median():.0f}  max {ns.max():.0f}  "
                  f">1 in {int((ns > 1).sum())} of {len(ns)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
