#!/usr/bin/env python
"""experiment B -- stability maps over damage influx, nascent occupancy and
rescue allocation.

the 3-d map is the product grid (j, nu, chi):

    j    site-resolved damage influx, scaled
    nu   ordinary nascent-chain chaperone occupancy (N_free/K_N) -- a pure
         capacity load carrying NO damage influx
    chi  rescue allocation, the chaperone share of a fixed total rescue pool

for each (nu, chi) the low-burden branch is continued across j; the fold in j
is bracketed; the equilibrium count either side of the fold is recorded with
the independent blind finder so that "fold" is a checked interpretation rather
than an assumed one; and recovery time is measured at several distances from
the boundary to compare against the leading eigenvalue.

three tables are written:
    stability_map.tsv   one row per (j, nu, chi)
    fold_boundary.tsv   one row per (nu, chi)
    slowing_down.tsv    eigenvalue and recovery time vs distance to fold

usage:  python scripts/run_experiment_b.py --config configs/phase1/experiment_b.json
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from proteostasis import (Params, allocationParams, findEquilibria, findFold,
                          traceBranch, recoveryTime)
from proteostasis.provenance import loadConfig, writeTable, writeProvenance, canonicalJson
from proteostasis.sweeps import parallelMap, resolveWorkers

_CTX: dict = {}


def _rowTask(task: tuple):
    """all work for one (nu, chi) cell. module level so it is picklable."""
    nu, chi = task
    cfg = _CTX["cfg"]
    base = _CTX["base"]
    j_grid = _CTX["j_grid"]
    H = cfg["burden_threshold"]

    p0 = allocationParams(base.with_(nu=float(nu)), float(chi),
                          cfg.get("rescue_total", 1.0))

    fold = findFold(p0, "j", cfg["j_lo"], cfg["j_hi"],
                    n_march=cfg["fold_n_march"], tol_rel=cfg["fold_tol_rel"])
    j_fold = fold.fold_value

    # --- map ------------------------------------------------------------
    map_rows = []
    for bp in traceBranch(p0, "j", j_grid):
        r = dict(nu=float(nu), chi=float(chi), j=bp.value,
                 j_fold=np.nan if j_fold is None else j_fold,
                 exists=bp.eq is not None)
        if bp.eq is None:
            r.update(u=np.nan, a=np.nan, burden=np.nan, eig_real_max=np.nan,
                     stable=False, below_threshold=False, aggregate_fraction=np.nan)
        else:
            b = bp.eq.burden
            r.update(u=bp.eq.u, a=bp.eq.a, burden=b, eig_real_max=bp.eq.eig_real_max,
                     stable=bp.eq.stable, below_threshold=bool(b < H),
                     aggregate_fraction=bp.eq.a / max(b, 1e-300))
        map_rows.append(r)

    # --- equilibrium count either side of the bracketed fold --------------
    count_rows = []
    if j_fold is not None:
        for label, jj in (("below", j_fold * 0.95), ("at", j_fold),
                          ("above", min(j_fold * 1.05, cfg["j_hi"]))):
            eqs = findEquilibria(p0.with_(j=jj), n_grid=cfg["count_grid"], hi=1e5)
            count_rows.append(dict(nu=float(nu), chi=float(chi), side=label, j=jj,
                                   n_equilibria=len(eqs),
                                   n_stable=sum(e.stable for e in eqs),
                                   n_unstable=sum(not e.stable for e in eqs),
                                   max_burden_found=max((e.burden for e in eqs),
                                                        default=np.nan)))

    # --- critical slowing down -------------------------------------------
    slow_rows = []
    if j_fold is not None:
        for frac in cfg["slowdown_fractions"]:
            jj = float(frac) * j_fold
            pts = [bp for bp in traceBranch(p0, "j", [cfg["j_lo"], jj]) if bp.eq is not None]
            if not pts or pts[-1].value != jj:
                continue
            eq = pts[-1].eq
            tr = recoveryTime(p0.with_(j=jj), eq.u, eq.a,
                              kick=cfg["recovery_kick"], t_end=cfg["t_end"])
            slow_rows.append(dict(nu=float(nu), chi=float(chi), j=jj,
                                  j_fold=j_fold, distance=1.0 - float(frac),
                                  burden=eq.burden, eig_real_max=eq.eig_real_max,
                                  eig_timescale=-1.0 / eq.eig_real_max if eq.eig_real_max < 0 else np.nan,
                                  recovery_time=np.nan if tr is None else tr))

    fold_row = dict(nu=float(nu), chi=float(chi))
    fold_row.update(fold.asDict())
    fold_row["removal_ceiling"] = p0.rescueTotal and (p0.c_tot + (p0.rho_U + p0.rho_A) * p0.d_tot)
    if j_fold is not None:
        fold_row["fold_to_ceiling_ratio"] = j_fold / fold_row["removal_ceiling"]
    return map_rows, fold_row, count_rows, slow_rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--workers", type=int, default=None)
    args = ap.parse_args()

    cfg = loadConfig(args.config)
    outdir = Path(args.outdir or cfg["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    base = Params(**cfg.get("base_params", {})).validate()
    nu_grid = np.logspace(np.log10(cfg["nu_lo"]), np.log10(cfg["nu_hi"]), cfg["n_nu"])
    chi_grid = np.linspace(cfg["chi_lo"], cfg["chi_hi"], cfg["n_chi"])
    j_grid = np.logspace(np.log10(cfg["j_lo"]), np.log10(cfg["j_hi"]), cfg["n_j"])

    _CTX.update(cfg=cfg, base=base, j_grid=j_grid)
    tasks = [(float(n), float(c)) for n in nu_grid for c in chi_grid]
    workers = resolveWorkers(args.workers or cfg.get("workers"))
    print(f"[B] {len(tasks)} (nu,chi) cells x {len(j_grid)} influx values, "
          f"{workers} workers", flush=True)

    results = parallelMap(_rowTask, tasks, workers, progress_every=25)

    map_rows, fold_rows, count_rows, slow_rows = [], [], [], []
    for m, f, c, s in results:
        map_rows += m
        fold_rows.append(f)
        count_rows += c
        slow_rows += s

    dfm = pd.DataFrame(map_rows)
    dff = pd.DataFrame(fold_rows)
    dfc = pd.DataFrame(count_rows)
    dfs = pd.DataFrame(slow_rows)

    outputs = {
        "stability_map.tsv": writeTable(dfm, outdir / "stability_map.tsv"),
        "fold_boundary.tsv": writeTable(dff, outdir / "fold_boundary.tsv"),
        "equilibrium_counts.tsv": writeTable(dfc, outdir / "equilibrium_counts.tsv"),
        "slowing_down.tsv": writeTable(dfs, outdir / "slowing_down.tsv"),
    }

    # headline model consequences, stated as fractions rather than assertions
    have = dff["fold_value"].notna()
    summary = dict(
        experiment="B",
        n_cells=int(len(dff)),
        n_cells_with_fold=int(have.sum()),
        frac_cells_with_fold=float(have.mean()),
        fold_value_median=float(dff.loc[have, "fold_value"].median()) if have.any() else None,
        fold_to_ceiling_ratio_median=float(dff.loc[have, "fold_to_ceiling_ratio"].median())
        if have.any() and "fold_to_ceiling_ratio" in dff else None,
        nascent_monotone_frac=_monotoneFraction(dff, "nu"),
        allocation_interior_optimum_frac=_interiorOptimumFraction(dff),
        equilibrium_count_pattern=_countPattern(dfc),
        slowdown_spearman=_slowdownStat(dfs),
        n_map_rows=int(len(dfm)),
        frac_map_rows_with_state=float(dfm["exists"].mean()) if len(dfm) else None,
    )
    with open(outdir / "summary.json", "w") as fh:
        fh.write(canonicalJson(summary))
        fh.write("\n")
    from proteostasis.provenance import hashFile
    outputs["summary.json"] = hashFile(outdir / "summary.json")
    writeProvenance(outdir, cfg, outputs, extra=dict(runtime_s=time.time() - t0))

    print(canonicalJson(summary))
    print(f"[B] done in {time.time()-t0:.1f}s -> {outdir}")
    return 0


def _monotoneFraction(dff: pd.DataFrame, by: str) -> float | None:
    """fraction of chi columns in which j_fold decreases monotonically with nu.

    this is the sharp form of the load-coupling claim: ordinary nascent-chain
    occupancy carries no damage at all, so any narrowing of the feasible influx
    window is caused purely by capacity competition.
    """
    ok, tot = 0, 0
    for chi, g in dff.groupby("chi"):
        g = g.sort_values(by).dropna(subset=["fold_value"])
        if len(g) < 3:
            continue
        tot += 1
        ok += int(np.all(np.diff(g["fold_value"].to_numpy()) <= 1e-12))
    return None if tot == 0 else ok / tot


def _interiorOptimumFraction(dff: pd.DataFrame) -> float | None:
    """fraction of nu rows whose feasibility-maximising allocation is interior."""
    ok, tot = 0, 0
    for nu, g in dff.groupby("nu"):
        g = g.sort_values("chi").dropna(subset=["fold_value"])
        if len(g) < 3:
            continue
        tot += 1
        k = int(np.argmax(g["fold_value"].to_numpy()))
        ok += int(0 < k < len(g) - 1)
    return None if tot == 0 else ok / tot


def _countPattern(dfc: pd.DataFrame) -> dict | None:
    """does the equilibrium count go 2 -> 0 across the bracketed fold?"""
    if not len(dfc):
        return None
    piv = dfc.pivot_table(index=["nu", "chi"], columns="side", values="n_equilibria")
    piv = piv.dropna()
    if not len(piv):
        return None
    saddle_node = (piv.get("below", 0) == 2) & (piv.get("above", 0) == 0)
    return dict(n=int(len(piv)),
                frac_two_to_zero=float(saddle_node.mean()),
                frac_any_third_state=float((piv.max(axis=1) > 2).mean()))


def _slowdownStat(dfs: pd.DataFrame) -> dict | None:
    """rank correlation between distance-to-fold and relaxation rate."""
    if not len(dfs):
        return None
    from scipy.stats import spearmanr
    g = dfs.dropna(subset=["distance", "eig_real_max"])
    out = {"n": int(len(g))}
    if len(g) > 5:
        r, p = spearmanr(g["distance"], -g["eig_real_max"])
        out["eig_rate_vs_distance_rho"] = float(r)
        out["eig_rate_vs_distance_p"] = float(p)
    g2 = dfs.dropna(subset=["distance", "recovery_time"])
    if len(g2) > 5:
        r, p = spearmanr(g2["distance"], g2["recovery_time"])
        out["recovery_time_vs_distance_rho"] = float(r)
        out["recovery_time_vs_distance_p"] = float(p)
        out["n_recovery"] = int(len(g2))
    return out


if __name__ == "__main__":
    raise SystemExit(main())
