"""C5a: does ONE cycle span the window, or is l1 being interpolated across it?

Section 7 currently says a stable limit cycle "is born as the influx enters the
window and is destroyed as it leaves". That is `l1 < 0` computed at the two
ENDPOINTS and an interpolation across the middle. It does not follow. Between the
two Hopf points the cycle could be destroyed and recreated, could collide with
another cycle at a fold of cycles, or could lose stability without the
equilibrium regaining it. Nothing computed so far looks inside the band.

THE TEST. For each resolved window network, take influx values spanning
`[j_H1, j_H2]`, locate the branch equilibrium at each, perturb it, and integrate
the full nonlinear system long enough to see what the trajectory does. Classify
by the envelope of the late trajectory:

    cycle      bounded oscillation with a non-growing envelope
    fixed      settles to a point (which may or may not be where it started)
    divergent  envelope still growing at the horizon, or leaves the orthant

A window spanned by one cycle gives `cycle` at every interior influx. Anything
else is a finding, and the point of running it is that both outcomes are
publishable and only one of them is what Section 7 currently asserts.

WHY THE ENVELOPE AND NOT THE ESCAPE FLAG. `hopf_check.integrationTest` stops at
`d_escape`, a ten-thousandfold amplification, which a stable cycle reaches just
as a divergence does (decision D058). This integration deliberately has NO
terminal event: it runs past that point, because the question is where the
trajectory ends up, not whether it left.
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

T_CAP = 5.0e4


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


def fateAt(p: M.Params, j: float, u_eq: float, a_eq: float,
           rel: float = 1e-6) -> dict | None:
    """integrate past escape and classify the trajectory by its late envelope."""
    q = p.with_(j=float(j)).validate()
    try:
        J = M.jacobian(u_eq, a_eq, q)
    except (M.ModelError, np.linalg.LinAlgError):
        return None
    ev = np.linalg.eigvals(J)
    sigma = float(np.max(ev.real))
    omega = float(np.max(np.abs(ev.imag)))
    scale = max(abs(u_eq), abs(a_eq), 1e-300)

    # long enough to grow AND to settle: many e-folds and many periods
    t_want = 60.0 / max(abs(sigma), 1e-12)
    if omega > 1e-12:
        t_want = max(t_want, 200.0 * 2.0 * np.pi / omega)
    t_end = float(min(t_want, T_CAP))

    def f(t, x):
        if x[0] <= 0.0 or x[1] <= 0.0 or not np.all(np.isfinite(x)):
            return np.zeros(2)
        try:
            return M.rhsVector(x, q)
        except (M.ModelError, OverflowError, ValueError):
            return np.zeros(2)

    x0 = [u_eq * (1.0 + rel), a_eq * (1.0 + rel)]
    try:
        sol = solve_ivp(f, (0.0, t_end), x0, method="LSODA", rtol=1e-10,
                        atol=1e-14 * scale, dense_output=False,
                        t_eval=np.linspace(0.0, t_end, 4000))
    except Exception:
        return None
    if not sol.success or sol.y.shape[1] < 100:
        return None

    U, A, T = sol.y[0], sol.y[1], sol.t
    if not np.all(np.isfinite(U)) or not np.all(np.isfinite(A)):
        return {"fate": "divergent", "reason": "nonfinite", "sigma": sigma,
                "omega": omega, "j": j}

    n = len(T)
    late = slice(int(0.60 * n), n)
    mid = slice(int(0.30 * n), int(0.60 * n))

    def envelope(sl):
        uu, aa = U[sl], A[sl]
        du = (uu.max() - uu.min()) / max(abs(uu.mean()), 1e-300)
        da = (aa.max() - aa.min()) / max(abs(aa.mean()), 1e-300)
        return max(du, da)

    e_late, e_mid = envelope(late), envelope(mid)
    growth = e_late / max(e_mid, 1e-300)
    drift = abs(U[late].mean() - U[mid].mean()) / max(abs(U[mid].mean()), 1e-300)

    if e_late < 1e-8:
        d0 = float(np.hypot(U[-1] - u_eq, A[-1] - a_eq)) / scale
        fate = "fixed_same" if d0 < 1e-4 else "fixed_other"
    elif growth > 3.0 or (growth > 1.2 and drift > 0.1):
        fate = "divergent"
    else:
        fate = "cycle"

    # periodicity: successive maxima of a(t) in the late window. a settled cycle
    # has both a stable period and a stable peak height; a slow transient has
    # neither, and the escape flag cannot tell them apart (D058).
    aa, tt = A[late], T[late]
    pk = np.nonzero((aa[1:-1] > aa[:-2]) & (aa[1:-1] >= aa[2:]))[0] + 1
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
        "growth": float(growth), "drift": float(drift),
        "u_eq": u_eq, "a_eq": a_eq, "t_end": t_end,
        "u_min": float(U[late].min()), "u_max": float(U[late].max()),
        # section 9's observable, computed rather than predicted: does the mean
        # burden over the cycle stay near the equilibrium while the amplitude
        # grows, or does the cycle sit at a higher mean burden?
        "a_amp_rel": float((A[late].max() - A[late].min()) / max(a_eq, 1e-300)),
        "a_mean_over_eq": float(A[late].mean() / max(a_eq, 1e-300)),
        "u_mean_over_eq": float(U[late].mean() / max(u_eq, 1e-300)),
        "n_peaks": n_peaks, "period": period, "period_cv": period_cv,
        "peak_cv": peak_cv, "period_predicted": period_pred,
    }


def _worker(item):
    name, p, u_f, a_f, jH1, jH2, j_crit = item
    rows = []
    fracs = np.linspace(0.08, 0.92, 7)          # interior of the window only
    for fr in fracs:
        j = jH1 + fr * (jH2 - jH1)
        eq = equilibriumAtInflux(p, j, u_f, a_f)
        if eq is None:
            rows.append({"name": name, "frac": float(fr), "j": float(j),
                         "fate": "no_equilibrium"})
            continue
        r = fateAt(p, j, eq[0], eq[1])
        if r is None:
            rows.append({"name": name, "frac": float(fr), "j": float(j),
                         "fate": "not_evaluable"})
            continue
        r.update({"name": name, "frac": float(fr),
                  "j_over_jcrit": float(j / j_crit)})
        rows.append(r)
    return rows


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
        nested = pool.map(_worker, items, chunksize=1)
    D = pd.DataFrame([r for rows in nested for r in rows])
    D.to_csv(OUT_TSV, sep="\t", index=False)

    counts = D["fate"].value_counts().to_dict()
    per_net = D.groupby("name")["fate"].apply(
        lambda s: "all_cycle" if (s == "cycle").all()
        else ("no_cycle" if (s == "cycle").sum() == 0 else "mixed"))

    out = {
        "n_networks": int(len(items)),
        "n_points": int(len(D)),
        "fate_counts": {k: int(v) for k, v in counts.items()},
        "networks_all_interior_points_cycle": int((per_net == "all_cycle").sum()),
        "networks_mixed": int((per_net == "mixed").sum()),
        "networks_no_cycle": int((per_net == "no_cycle").sum()),
    }
    ev = D[D["fate"].isin(["cycle", "divergent", "fixed_same", "fixed_other"])]
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
