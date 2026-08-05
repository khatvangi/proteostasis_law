"""C2 and C3: transversality and the first Lyapunov coefficient, per CROSSING.

`tr J = 0` with `det J > 0` is NECESSARY for a Hopf bifurcation, not sufficient.
`hopf_check.py` called it "a Hopf point by definition"; that is false and is
corrected there in the same commit as this file. A crossing is a Hopf point when,
in addition, the pair crosses the axis with nonzero speed (C2) and the cubic
normal-form coefficient is nonzero (C3). The sign of that coefficient decides
whether the cycle born there is stable.

PER CROSSING, NOT PER NETWORK. Task B7 showed the crossers split in two:

    group F   tr J > 0 at the fold -- one crossing, then the fold
    group W   tr J < 0 at the fold -- two crossings, stability recovered between

A group W network has an instability WINDOW interior to its stable branch, and
its two crossings have opposite transversality sign by construction. Reporting
one number per network would average them.

THE FORMULA. Guckenheimer and Holmes (1983) eq. 3.4.11, for a planar system put
in the form

    x' = -omega*y + f(x,y)
    y' =  omega*x + g(x,y)

    a = (1/16) [ f_xxx + f_xyy + g_xxy + g_yyy ]
      + (1/(16*omega)) [ f_xy (f_xx + f_yy) - g_xy (g_xx + g_yy)
                         - f_xx g_xx + f_yy g_yy ]

`a < 0` is supercritical: a stable limit cycle exists on the side where the
equilibrium is unstable. `a > 0` is subcritical: an unstable cycle shrinks onto
the equilibrium as the parameter rises, and the loss of stability is hard.

THE TRANSFORMATION. At the crossing, `J` has eigenvalues `+/- i*omega` with
`omega = sqrt(det J)`. If `v = p + i q` is an eigenvector for `+i omega`, then
`J p = -omega q` and `J q = omega p`, so with `T = [p | -q]` as columns,
`T^{-1} J T = [[0, -omega], [omega, 0]]`, which is the form the formula wants.
Scaling both columns by one constant preserves that form and rescales `a` by a
positive factor, so the SIGN -- the only thing claimed -- is invariant.

WHY THE SELF-TEST IS NOT OPTIONAL. Every published sign convention for this
formula differs by an author's choice of which rotation is positive, and a sign
error here reverses the paper's conclusion about what happens at the crossing.
`selfTest()` runs the differentiator against the two textbook normal forms whose
coefficients are known in closed form (`a = -1` and `a = +1`) and is asserted
before any model result is computed.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import fold_theorem as FT  # noqa: E402
import genericity as GEN  # noqa: E402

COMPUTED = REPO_ROOT / "data" / "computed"
OUT = COMPUTED / "hopf_lyapunov.json"
OUT_TSV = COMPUTED / "hopf_lyapunov.tsv"


# ---------------------------------------------------------------------------
# the coefficient
# ---------------------------------------------------------------------------


def _derivs(F, h: float):
    """all partials of a planar map at 0, orders 2 and 3, by central differences.

    13-point stencil. F must satisfy F(0,0) = 0 and be evaluable on the stencil;
    the caller supplies a step small enough to stay in the positive orthant.
    """
    def e(x, y):
        return np.asarray(F(x, y), float)

    f0 = e(0.0, 0.0)
    px, mx = e(h, 0.0), e(-h, 0.0)
    py, my = e(0.0, h), e(0.0, -h)
    p2x, m2x = e(2 * h, 0.0), e(-2 * h, 0.0)
    p2y, m2y = e(0.0, 2 * h), e(0.0, -2 * h)
    pp, pm = e(h, h), e(h, -h)
    mp_, mm = e(-h, h), e(-h, -h)

    d = {}
    d["xx"] = (px - 2 * f0 + mx) / h ** 2
    d["yy"] = (py - 2 * f0 + my) / h ** 2
    d["xy"] = (pp - pm - mp_ + mm) / (4 * h ** 2)
    d["xxx"] = (p2x - 2 * px + 2 * mx - m2x) / (2 * h ** 3)
    d["yyy"] = (p2y - 2 * py + 2 * my - m2y) / (2 * h ** 3)
    # d^3/dx^2 dy  =  d/dy of the second x-difference
    d["xxy"] = ((pp - 2 * py + mp_) - (pm - 2 * my + mm)) / (2 * h ** 3)
    d["xyy"] = ((pp - 2 * px + pm) - (mp_ - 2 * mx + mm)) / (2 * h ** 3)
    return d


def lyapunovFromField(F, omega: float, h: float) -> float:
    """G&H 3.4.11 for a field already in canonical form, nonlinear part F."""
    d = _derivs(F, h)
    f, g = 0, 1                       # component indices
    term1 = (d["xxx"][f] + d["xyy"][f] + d["xxy"][g] + d["yyy"][g])
    term2 = (d["xy"][f] * (d["xx"][f] + d["yy"][f])
             - d["xy"][g] * (d["xx"][g] + d["yy"][g])
             - d["xx"][f] * d["xx"][g]
             + d["yy"][f] * d["yy"][g])
    return float(term1 / 16.0 + term2 / (16.0 * omega))


def selfTest() -> dict:
    """the two textbook normal forms, whose coefficients are known exactly.

    supercritical:  x' = -y + x(-(x^2+y^2)),  y' = x + y(-(x^2+y^2))   -> a = -1
    subcritical:    the same with the cubic sign flipped                -> a = +1

    A convention error in the transformation or in the formula reverses these,
    which is exactly the failure this guards.
    """
    out = {}
    for label, s in (("supercritical", -1.0), ("subcritical", +1.0)):
        def F(x, y, s=s):
            r2 = x * x + y * y
            return (s * x * r2, s * y * r2)      # nonlinear part only
        vals = [lyapunovFromField(F, 1.0, hh) for hh in (1e-2, 3e-3, 1e-3)]
        out[label] = {"expected": s, "computed": vals,
                      "sign_correct": all(np.sign(v) == np.sign(s) for v in vals),
                      "max_abs_err": float(max(abs(v - s) for v in vals))}
    return out


def canonicalTransform(J: np.ndarray):
    """T with T^{-1} J T = [[0,-omega],[omega,0]], and omega."""
    w, V = np.linalg.eig(J)
    k = int(np.argmax(w.imag))            # the eigenvalue with +i omega
    omega = float(w[k].imag)
    if not omega > 0.0:
        return None, None
    v = V[:, k]
    p, q = np.real(v), np.imag(v)
    T = np.column_stack([p, -q])
    if abs(np.linalg.det(T)) < 1e-300:
        return None, None
    # one common scale factor: preserves the canonical form, rescales `a` by a
    # positive constant, leaves the sign alone.
    T = T / max(np.linalg.norm(p), np.linalg.norm(q))
    return T, omega


def lyapunovAt(u: float, a: float, p: M.Params, h_rel: float = 3e-3):
    """first Lyapunov coefficient at a crossing state, with a step-size check."""
    try:
        j_eq = FT.removalR(u, a, p)
        q = p.with_(j=float(j_eq)).validate()
        res = float(np.hypot(*M.rhsVector([u, a], q)))
        J = M.jacobian(u, a, q)
    except (M.ModelError, np.linalg.LinAlgError, OverflowError, ValueError):
        return None
    scale = max(abs(u), abs(a), 1e-300)
    if res / scale > 1e-8:
        return None                        # not an equilibrium; do not test

    T, omega = canonicalTransform(np.asarray(J, float))
    if T is None:
        return None
    Tinv = np.linalg.inv(T)

    def F(x, y):
        z = np.array([u, a]) + T @ np.array([x, y])
        if z[0] <= 0.0 or z[1] <= 0.0:
            raise M.ModelError("stencil left the positive orthant")
        return Tinv @ M.rhsVector(z, q)

    vals = []
    for hh in (h_rel * scale, 0.4 * h_rel * scale, 2.5 * h_rel * scale):
        try:
            vals.append(lyapunovFromField(F, omega, hh))
        except (M.ModelError, OverflowError, ValueError):
            vals.append(np.nan)
    v = np.array(vals, float)
    good = v[np.isfinite(v)]
    if good.size == 0:
        return None
    signs_agree = bool(np.all(np.sign(good) == np.sign(good[0])))
    spread = float(np.max(np.abs(good)) / max(np.min(np.abs(good)), 1e-300))
    return {
        "l1": float(good[0]), "l1_all": [float(x) for x in v],
        "omega": omega, "det_J": float(np.linalg.det(J)),
        "tr_J": float(J[0, 0] + J[1, 1]),
        "sign_stable": signs_agree, "magnitude_spread": spread,
        "j_eq": float(j_eq),
    }


# ---------------------------------------------------------------------------
# locating the crossings, and C2
# ---------------------------------------------------------------------------


def _traceCrossings(p: M.Params, u_star: float, a_star: float, n: int = 150):
    """every tr J sign change on the projected branch, refined and characterised."""
    import hopf_check as HC
    out = HC.branchProfile(p, u_star, a_star, n=n)
    if out is None:
        return None
    B = out["branch"].sort_values("j").reset_index(drop=True)
    tr = B["tr_J"].to_numpy(float)
    idx = np.nonzero(np.sign(tr[:-1]) * np.sign(tr[1:]) < 0)[0]

    crossings = []
    for k in idx:
        x0 = np.array([B["u"].iloc[k], B["a"].iloc[k]])
        x1 = np.array([B["u"].iloc[k + 1], B["a"].iloc[k + 1]])
        t0, t1 = tr[k], tr[k + 1]

        # bisection on the segment, PROJECTING onto {G = 0} at each trial so the
        # refined point is on the nullcline and not on the chord.
        lo, hi, flo = 0.0, 1.0, t0
        pt = None
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            z = GEN._projectOntoNullcline(*(x0 + mid * (x1 - x0)), p)
            if z is None:
                break
            try:
                Jm = M.jacobian(z[0], z[1], p)
            except (M.ModelError, np.linalg.LinAlgError):
                break
            tm = float(Jm[0, 0] + Jm[1, 1])
            pt = (z[0], z[1], tm)
            if abs(tm) < 1e-13:
                break
            if np.sign(tm) == np.sign(flo):
                lo, flo = mid, tm
            else:
                hi = mid
        if pt is None:
            continue
        uu, aa, tm = pt
        crossings.append({
            "u": uu, "a": aa, "tr_at_cross": tm,
            "direction": float(np.sign(t1 - t0)),
            "index": len(crossings),
        })
    return crossings, B, out


def _transversality(uu: float, aa: float, p: M.Params, B: pd.DataFrame):
    """d(tr J)/dj at a crossing, from the two nearest branch points.

    Re(lambda) = tr J / 2 along the branch, so the Hopf transversality condition
    is d(tr J)/dj != 0. The branch is already traced; this is a difference on it.
    """
    d = np.hypot(B["u"] - uu, B["a"] - aa).to_numpy()
    order = np.argsort(d)[:4]
    sub = B.iloc[np.sort(order)]
    if len(sub) < 2:
        return np.nan
    jj = sub["j"].to_numpy(float)
    tt = sub["tr_J"].to_numpy(float)
    if np.ptp(jj) <= 0.0:
        return np.nan
    return float(np.polyfit(jj, tt, 1)[0])


def _worker(item):
    name, p, u, a, group = item
    try:
        got = _traceCrossings(p, u, a)
    except Exception:
        return []
    if got is None:
        return []
    crossings, B, _ = got
    rows = []
    for c in crossings:
        r = {"name": name, "group": group, "crossing_index": c["index"],
             "u": c["u"], "a": c["a"], "tr_at_cross": c["tr_at_cross"],
             "direction": c["direction"]}
        r["dtr_dj"] = _transversality(c["u"], c["a"], p, B)
        try:
            j_crit = float(FT.removalR(u, a, p))
        except Exception:
            j_crit = np.nan
        r["j_crit"] = j_crit
        L = lyapunovAt(c["u"], c["a"], p)
        if L is None:
            r.update({"l1": np.nan, "sign_stable": False, "omega": np.nan,
                      "det_J": np.nan, "j_eq": np.nan, "magnitude_spread": np.nan})
        else:
            r.update({k: L[k] for k in ("l1", "sign_stable", "omega", "det_J",
                                        "j_eq", "magnitude_spread")})
        rows.append(r)
    return rows


def run(workers: int | None = None) -> dict:
    st = selfTest()
    if not all(v["sign_correct"] for v in st.values()):
        raise SystemExit(f"Lyapunov self-test FAILED: {json.dumps(st, indent=2)}")

    Z = pd.read_csv(COMPUTED / "hopf_zero_counts.tsv", sep="\t")
    groups = dict(zip(Z["name"], Z["group"]))

    run_dir = FT.phase1RunDir()
    c = pd.read_csv(run_dir / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    c = c[pd.to_numeric(c["fold_burden"], errors="coerce").notna()]

    items = []
    for i, r in c.iterrows():
        nm = f"draw{i}"
        if nm not in groups:
            continue
        try:
            p = FT.paramsFromSampleRow(r)
            u0, a0 = FT.foldStateFromSampleRow(r)
        except Exception:
            continue
        s = GEN.polishFold(p, u0, a0) or FT.foldSolve(p)
        if s is None:
            continue
        items.append((nm, p, float(s[1]), float(s[2]), groups[nm]))

    if workers is None:
        workers = max(1, min(56, (os.cpu_count() or 2) - 8))
    with mp.get_context("fork").Pool(processes=workers) as pool:
        nested = pool.map(_worker, items, chunksize=2)
    D = pd.DataFrame([r for rows in nested for r in rows])
    D.to_csv(OUT_TSV, sep="\t", index=False)

    def summarise(sub: pd.DataFrame) -> dict:
        if sub.empty:
            return {"n": 0}
        ev = sub[sub["l1"].notna() & sub["sign_stable"]]
        tv = sub[sub["dtr_dj"].notna()]
        return {
            "n_crossings": int(len(sub)),
            "n_networks": int(sub["name"].nunique()),
            # C2
            "n_transversality_evaluable": int(len(tv)),
            "abs_dtr_dj_min": float(tv["dtr_dj"].abs().min()) if len(tv) else None,
            "n_dtr_dj_positive": int((tv["dtr_dj"] > 0).sum()),
            "n_dtr_dj_negative": int((tv["dtr_dj"] < 0).sum()),
            # C3
            "n_l1_evaluable": int(len(ev)),
            "n_l1_sign_unstable": int((sub["l1"].notna()
                                       & ~sub["sign_stable"]).sum()),
            "n_l1_not_computed": int(sub["l1"].isna().sum()),
            "n_supercritical": int((ev["l1"] < 0).sum()),
            "n_subcritical": int((ev["l1"] > 0).sum()),
            "abs_l1_min": float(ev["l1"].abs().min()) if len(ev) else None,
            "abs_l1_median": float(ev["l1"].abs().median()) if len(ev) else None,
        }

    out = {"self_test": st,
           "n_networks_attempted": len(items),
           "all": summarise(D)}
    for g in ("terminal", "window"):
        out[g] = summarise(D[D["group"] == g])
    # group W's two crossings, taken separately
    W = D[D["group"] == "window"]
    for k in (0, 1):
        out[f"window_crossing_{k}"] = summarise(W[W["crossing_index"] == k])
    return out


def main() -> int:
    o = run()
    OUT.write_text(json.dumps(o, indent=2) + "\n")
    print(json.dumps(o, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
