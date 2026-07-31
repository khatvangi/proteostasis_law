#!/usr/bin/env python
"""phase 2 step 3 -- basin scans for every case the dense search left standing.

a root finder reporting Re(lambda) < 0 at two separated points is NOT evidence
of bistability. it is evidence of two locally stable roots, which is compatible
with (a) genuine bistability, (b) one attractor plus a numerical twin, and
(c) an attractor whose basin is too thin for any trajectory to reach. this step
separates those.

three probes per evaluation point:

  grid          a deterministic log-spaced grid of initial conditions across the
                whole state box. reproducible, blind, and it measures basin
                VOLUME rather than basin existence.
  eigendirection displacements along each stable root's own eigenvectors at
                several magnitudes. this is the local test: an attractor that
                does not recapture a small displacement along its own slow mode
                is not an attractor.
  saddle        displacements along the intervening saddle's UNSTABLE
                eigendirection, in both signs. this is the decisive test. the
                unstable manifold of a separating saddle must leave along the
                two sides of the separatrix, so if the two sides land on two
                DIFFERENT stable roots, those roots cannot be numerical
                duplicates of each other and both have nonempty basins.

every trajectory is additionally checked for positivity (the nonnegative orthant
is forward invariant, so a negative excursion is an integrator failure) and for
conserved-resource balance at sampled points along the path.

usage:  PHASE2_AUDIT_DIR=... python scripts/phase2/c03_basins.py --workers 20
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _audit_common import (auditRoot, loadProvenance, loadSamples, paramsFromRow,  # noqa: E402
                           parallelMap, runRoot, writeJson, writeTable)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from proteostasis.model import ModelError, jacobian, poolBalanceResiduals  # noqa: E402
from proteostasis.simulate import simulate  # noqa: E402

#: a trajectory endpoint counts as "reached root r" when it matches r to this
#: relative tolerance in both coordinates. deliberately loose (5%) -- the
#: question is which attractor was reached, not where it is.
REACH_RTOL = 0.05

_CTX: dict = {}


def _classifyEndpoint(tr, roots, blowup: float) -> str:
    if tr.status == "negative":
        return "positivity_violation"
    if tr.status == "blowup":
        return "blowup"
    if tr.status == "error":
        return "integration_error"
    for k, r in enumerate(roots):
        du = abs(tr.final_u - r["u"]) / max(abs(r["u"]), 1e-300)
        da = abs(tr.final_a - r["a"]) / max(abs(r["a"]), 1e-300)
        if du < REACH_RTOL and da < REACH_RTOL:
            return f"root{k}"
    if tr.final_u + tr.final_a > 0.5 * blowup:
        return "escaping"
    return "unresolved"


def _runOne(p, u0, a0, roots, cfg) -> dict:
    try:
        tr = simulate(p, float(u0), float(a0), t_end=cfg["t_end"], n_out=cfg["n_out"],
                      blowup=cfg["blowup"], rtol=1e-9, atol=1e-12)
    except (ModelError, ValueError, FloatingPointError) as exc:
        return dict(u0=float(u0), a0=float(a0), outcome="setup_error",
                    message=f"{type(exc).__name__}: {exc}")
    rec = dict(u0=float(u0), a0=float(a0), outcome=_classifyEndpoint(tr, roots, cfg["blowup"]),
               status=tr.status, final_u=tr.final_u, final_a=tr.final_a,
               final_rate=tr.final_rate, min_u=tr.min_u, min_a=tr.min_a)
    # conserved-resource balance sampled along the path, not only at the end
    worst = 0.0
    idx = np.unique(np.linspace(0, len(tr.u) - 1, min(12, len(tr.u))).astype(int))
    for i in idx:
        u, a = max(float(tr.u[i]), 0.0), max(float(tr.a[i]), 0.0)
        try:
            bal = poolBalanceResiduals(u, a, p)
        except (ModelError, ValueError, FloatingPointError):
            worst = np.inf
            break
        worst = max(worst, max(abs(v) for v in bal.values()))
    rec["max_pool_balance_on_path"] = float(worst)
    rec["positivity_ok"] = bool(tr.min_u >= -1e-9 and tr.min_a >= -1e-9)
    return rec


def _eigenDisplacements(p, r, mags) -> list:
    """(u0, a0, label) displacements along a root's own eigenvectors."""
    out = []
    try:
        w, V = np.linalg.eig(jacobian(r["u"], r["a"], p))
    except (ModelError, ValueError, FloatingPointError, np.linalg.LinAlgError):
        return out
    order = np.argsort(w.real)                    # slow (largest real) last
    for rank, k in enumerate(order):
        v = np.real(V[:, k])
        n = np.linalg.norm(v)
        if not np.isfinite(n) or n == 0.0:
            continue
        v = v / n
        for s in (+1.0, -1.0):
            for m in mags:
                d = m * max(r["u"] + r["a"], 1e-12) * s * v
                u0, a0 = r["u"] + d[0], r["a"] + d[1]
                if u0 <= 0.0 or a0 <= 0.0:
                    continue
                out.append((float(u0), float(a0),
                            f"eig{rank}_{'p' if s > 0 else 'm'}_{m:g}"))
    return out


def _task(item):
    idx, frac = item
    cfg, df, bcfg = _CTX["cfg"], _CTX["df"], _CTX["bcfg"]
    pts = _CTX["points"][(idx, frac)]
    row = df.loc[df["sample_index"] == idx].iloc[0]
    rec = dict(sample_index=int(idx), frac=float(frac), error="")
    try:
        p = paramsFromRow(row, cfg["base_params"], cfg.get("rescue_total", 1.0))
        pj = p.with_(j=float(pts["j"]))
        roots = pts["roots"]                       # unanimous roots, burden-ordered
        stable_idx = [k for k, r in enumerate(roots) if r["stable"]]
        saddle_idx = [k for k, r in enumerate(roots) if r["kind"] == "saddle"]
        rec.update(j=float(pts["j"]), n_roots=len(roots),
                   n_stable=len(stable_idx), n_saddle=len(saddle_idx),
                   root_burdens=[r["burden"] for r in roots],
                   root_kinds=[r["kind"] for r in roots])

        # --- probe 1: deterministic log grid -------------------------------
        g = np.logspace(np.log10(bcfg["grid_lo"]), np.log10(bcfg["grid_hi"]),
                        bcfg["n_grid"])
        grid_rows = [_runOne(pj, u0, a0, roots, bcfg) for u0 in g for a0 in g]
        for r in grid_rows:
            r["probe"] = "grid"

        # --- probe 2: eigendirection displacements at each stable root ------
        eig_rows = []
        for k in stable_idx:
            for u0, a0, lab in _eigenDisplacements(pj, roots[k], bcfg["eig_mags"]):
                r = _runOne(pj, u0, a0, roots, bcfg)
                r.update(probe="eigen", from_root=k, label=lab)
                eig_rows.append(r)

        # --- probe 3: the saddle's unstable manifold ------------------------
        sad_rows = []
        for k in saddle_idx:
            for u0, a0, lab in _eigenDisplacements(pj, roots[k], bcfg["saddle_mags"]):
                r = _runOne(pj, u0, a0, roots, bcfg)
                r.update(probe="saddle", from_root=k, label=lab)
                sad_rows.append(r)

        allrows = grid_rows + eig_rows + sad_rows
        rec["n_trajectories"] = len(allrows)
        rec["outcome_counts"] = {k: int(v) for k, v in
                                 pd.Series([r["outcome"] for r in allrows])
                                 .value_counts().items()}
        rec["grid_outcome_counts"] = {k: int(v) for k, v in
                                      pd.Series([r["outcome"] for r in grid_rows])
                                      .value_counts().items()}
        # basin evidence per stable root
        basins = {}
        for k in stable_idx:
            tag = f"root{k}"
            basins[tag] = dict(
                burden=roots[k]["burden"],
                grid_hits=int(sum(r["outcome"] == tag for r in grid_rows)),
                grid_fraction=float(np.mean([r["outcome"] == tag for r in grid_rows])),
                eigen_self_return=int(sum(r["outcome"] == tag and r.get("from_root") == k
                                          for r in eig_rows)),
                eigen_self_total=int(sum(r.get("from_root") == k for r in eig_rows)),
                saddle_hits=int(sum(r["outcome"] == tag for r in sad_rows)),
            )
        rec["basins"] = basins
        # decisive: does one saddle's unstable manifold reach two DIFFERENT
        # stable roots?
        reached = set()
        for k in saddle_idx:
            hits = {r["outcome"] for r in sad_rows
                    if r.get("from_root") == k and r["outcome"].startswith("root")}
            if len({h for h in hits if int(h[4:]) in stable_idx}) >= 2:
                rec["saddle_separates"] = True
            reached |= hits
        rec.setdefault("saddle_separates", False)
        rec["distinct_stable_reached"] = int(len(
            {h for h in reached if h.startswith("root") and int(h[4:]) in stable_idx}))
        rec["distinct_stable_reached_grid"] = int(len(
            {r["outcome"] for r in grid_rows
             if r["outcome"].startswith("root") and int(r["outcome"][4:]) in stable_idx}))
        rec["n_stable_with_nonempty_grid_basin"] = int(
            sum(basins[f"root{k}"]["grid_hits"] > 0 for k in stable_idx))
        rec["n_stable_with_local_recapture"] = int(
            sum(basins[f"root{k}"]["eigen_self_return"] > 0
                and basins[f"root{k}"]["eigen_self_return"]
                == basins[f"root{k}"]["eigen_self_total"] for k in stable_idx))
        rec["any_positivity_violation"] = bool(
            any(not r.get("positivity_ok", True) for r in allrows))
        rec["max_pool_balance_on_paths"] = float(
            max((r.get("max_pool_balance_on_path", 0.0) for r in allrows), default=0.0))
        rec["trajectories"] = allrows
    except Exception as exc:                                      # noqa: BLE001
        import traceback
        rec["error"] = f"{type(exc).__name__}: {exc}\n{traceback.format_exc()[:600]}"
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=None)
    ap.add_argument("--audit-dir", default=None)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--n-grid", type=int, default=21)
    ap.add_argument("--t-end", type=float, default=5.0e4)
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    run = runRoot(args.run)
    out = auditRoot(args.audit_dir)
    df = loadSamples(run)
    cfg = loadProvenance(run)["config"]

    with open(out / "c02_dense_multistability.json") as fh:
        c02 = json.load(fh)

    # basin work is only meaningful where at least two stable roots survived all
    # three finders. anything less has nothing to have a basin.
    points = {}
    for rec in c02["records"]:
        if rec.get("error"):
            continue
        for pt in rec.get("points", []):
            unan = [r for r in pt["roots"] if r["n_methods"] == 3]
            if sum(bool(r["stable"]) for r in unan) >= 2:
                points[(int(rec["sample_index"]), float(pt["frac"]))] = dict(
                    j=pt["j"], roots=sorted(unan, key=lambda r: r["burden"]),
                    group=rec["group"])
    items = sorted(points)
    if args.limit:
        items = items[:args.limit]

    bcfg = dict(n_grid=args.n_grid, grid_lo=1e-8, grid_hi=1.0e5, t_end=args.t_end,
                n_out=200, blowup=1.0e6, eig_mags=(1e-3, 1e-2, 1e-1, 0.5),
                saddle_mags=(1e-6, 1e-4, 1e-2), reach_rtol=REACH_RTOL)
    _CTX.update(cfg=cfg, df=df, bcfg=bcfg, points=points)

    print(f"[c03] {len(items)} (draw, evaluation point) pairs with >=2 unanimous "
          f"stable roots, {args.workers} workers", flush=True)
    t0 = time.time()
    recs = parallelMap(_task, items, args.workers, chunksize=1)
    print(f"[c03] done in {time.time()-t0:.1f}s", flush=True)

    writeJson(dict(basin_config=bcfg, records=recs), out / "c03_basins.json")
    flat = pd.DataFrame([{k: (json.dumps(v) if isinstance(v, (dict, list)) else v)
                          for k, v in r.items() if k != "trajectories"} for r in recs])
    writeTable(flat, out / "c03_basins.tsv")

    ok = [r for r in recs if not r.get("error")]
    print(json.dumps(dict(
        n_pairs=len(recs), n_errors=len(recs) - len(ok),
        n_with_two_nonempty_grid_basins=int(
            sum(r["n_stable_with_nonempty_grid_basin"] >= 2 for r in ok)),
        n_with_saddle_separatrix=int(sum(r["saddle_separates"] for r in ok)),
        n_with_positivity_violation=int(sum(r["any_positivity_violation"] for r in ok)),
    ), indent=2))
    print(f"-> {out/'c03_basins.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
