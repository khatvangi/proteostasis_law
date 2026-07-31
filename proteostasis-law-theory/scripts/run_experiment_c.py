#!/usr/bin/env python
"""experiment C -- global parameter sweep: which qualitative predictions survive?

a latin hypercube over 16 kinetic, affinity and allocation parameters spanning
2-4 orders of magnitude each. for every draw the script asks a fixed list of
yes/no questions and reports the FRACTION of parameter space in which each
holds. the deliverable is a robustness table, not a confirmation: a prediction
that survives in 55% of draws is reported as surviving in 55% of draws.

questions asked per draw:

  C1  does a finite critical influx exist (the low branch terminates below the
      swept ceiling)?
  C2  how far below the analytic removal ceiling does collapse occur? the
      ceiling c_tot + (rho_U+rho_A)*d_tot is what naive capacity accounting
      would predict; j_fold/ceiling measures what the nonlinearity costs.
  C3  does raising ordinary nascent-chain occupancy (which carries NO damage)
      shrink the feasible influx window?
  C4  is there a second STABLE high-burden attractor? theory/SCOPE_AND_NONCLAIMS
      nonclaim 6-7 denies one for the minimal model; this searches for it over a
      wide box instead of assuming its absence.
  C5  does the leading eigenvalue approach zero as the fold is approached
      (critical slowing down)?
  C6  is the feasibility-maximising rescue allocation interior, or is it always
      at a boundary of the allocation range?

usage:  python scripts/run_experiment_c.py --config configs/phase1/experiment_c.json
"""

from __future__ import annotations

import argparse
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from proteostasis import (Params, ModelError, allocationParams, findEquilibria,
                          findFold, lowestStableEquilibrium, removalCeiling)
from proteostasis.equilibria import stableEquilibriumAt
from proteostasis.provenance import loadConfig, writeTable, writeProvenance, canonicalJson
from proteostasis.sweeps import (latinHypercube, paramsFromSample, parallelMap,
                                 rangesFromConfig, resolveWorkers)

_CTX: dict = {}


def _fold(p: Params, cfg: dict, cheap: bool = False):
    kw = cfg["fold_cheap"] if cheap else cfg["fold_full"]
    return findFold(p, "j", cfg["j_lo"], cfg["j_hi"], **kw)


def _sampleTask(item: tuple):
    """evaluate one LHS draw. module level so it is picklable."""
    idx, sample = item
    cfg = _CTX["cfg"]
    base = _CTX["base"]
    row = dict(sample_index=int(idx), error="")
    row.update({f"p_{k}": float(v) for k, v in sample.items()})

    try:
        p = paramsFromSample(sample, base, cfg.get("rescue_total", 1.0))
        row["removal_ceiling"] = removalCeiling(p)

        # is there any viable state at all at the smallest swept influx?
        if lowestStableEquilibrium(p.with_(j=cfg["j_lo"]), n_grid=9) is None:
            row.update(viable_at_j_lo=False, C1_fold_exists=False)
            return row
        row["viable_at_j_lo"] = True

        # --- C1 / C2 ------------------------------------------------------
        fold = _fold(p, cfg)
        jf = fold.fold_value
        row["C1_fold_exists"] = bool(fold.found)
        row["fold_reason"] = fold.reason
        row["j_fold"] = np.nan if jf is None else jf
        if jf is None:
            return row
        row["C2_fold_to_ceiling_ratio"] = jf / row["removal_ceiling"]
        row["C2_collapse_below_ceiling"] = bool(jf < row["removal_ceiling"])
        row["fold_burden"] = fold.eq_at_fold.burden if fold.eq_at_fold else np.nan
        row["fold_aggregate_fraction"] = (
            fold.eq_at_fold.a / max(fold.eq_at_fold.burden, 1e-300)
            if fold.eq_at_fold else np.nan)

        # --- C3 nascent-occupancy load coupling ---------------------------
        nu_hi = p.nu * cfg["nascent_multiplier"]
        fold_nu = _fold(p.with_(nu=nu_hi), cfg, cheap=True)
        row["j_fold_high_nu"] = np.nan if fold_nu.fold_value is None else fold_nu.fold_value
        if fold_nu.fold_value is not None:
            row["C3_nascent_shrinks_window"] = bool(fold_nu.fold_value < jf)
            row["C3_shrink_ratio"] = fold_nu.fold_value / jf

        # --- C4 second stable attractor -----------------------------------
        # a root finder reporting Re(lambda) < 0 is weaker evidence than a
        # trajectory actually converging, and a spurious root would look
        # identical. any candidate is therefore re-tested dynamically.
        n_stable_max, n_eq_max, confirmed = 0, 0, False
        for frac in cfg["attractor_fractions"]:
            eqs = findEquilibria(p.with_(j=frac * jf), n_grid=cfg["attractor_grid"],
                                 lo=1e-8, hi=cfg["attractor_hi"])
            n_stable_max = max(n_stable_max, sum(e.stable for e in eqs))
            n_eq_max = max(n_eq_max, len(eqs))
            stable = sorted((e for e in eqs if e.stable), key=lambda e: e.burden)
            if len(stable) >= 2 and _confirmsAttractor(p.with_(j=frac * jf),
                                                       stable[-1], cfg):
                confirmed = True
        row["max_stable_equilibria"] = int(n_stable_max)
        row["max_equilibria"] = int(n_eq_max)
        row["C4_second_stable_attractor"] = bool(n_stable_max >= 2)
        row["C4_second_attractor_confirmed"] = bool(confirmed)

        # --- C5 critical slowing down -------------------------------------
        eig = {}
        warm = None
        for frac in sorted(cfg["slowdown_fractions"]):
            eq = stableEquilibriumAt(p, "j", frac * jf, warm, blind=(warm is None))
            if eq is None:
                continue
            warm = (eq.u, eq.a)
            eig[frac] = eq.eig_real_max
        if len(eig) >= 2:
            f_far, f_near = min(eig), max(eig)
            row["eig_far"] = eig[f_far]
            row["eig_near"] = eig[f_near]
            row["C5_slowing_down"] = bool(abs(eig[f_near]) < abs(eig[f_far]))
            row["C5_rate_ratio"] = abs(eig[f_near]) / max(abs(eig[f_far]), 1e-300)

        # --- C6 rescue allocation ------------------------------------------
        chis = np.array(cfg["chi_scan"], dtype=float)
        jf_chi = []
        for chi in chis:
            fr = _fold(allocationParams(p, float(chi), cfg.get("rescue_total", 1.0)),
                       cfg, cheap=True)
            jf_chi.append(np.nan if fr.fold_value is None else fr.fold_value)
        jf_chi = np.array(jf_chi, dtype=float)
        ok = np.isfinite(jf_chi)
        for chi, v in zip(chis, jf_chi):
            row[f"j_fold_chi_{chi:g}"] = v
        if ok.sum() >= 3:
            k = int(np.nanargmax(jf_chi))
            row["C6_best_chi"] = float(chis[k])
            row["C6_interior_optimum"] = bool(0 < k < len(chis) - 1)
            row["C6_allocation_span"] = float(np.nanmax(jf_chi) / max(np.nanmin(jf_chi), 1e-300))
    except (ModelError, ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
        row["error"] = f"{type(exc).__name__}: {exc}"
    except Exception as exc:                                  # noqa: BLE001
        row["error"] = f"UNEXPECTED {type(exc).__name__}: {exc}\n{traceback.format_exc()[:400]}"
    return row


def _confirmsAttractor(p: Params, eq, cfg: dict) -> bool:
    """dynamic confirmation that a candidate high-burden equilibrium attracts.

    displace off it in both directions and integrate. it counts only if BOTH
    trajectories return to it -- not to the low branch, and not to escape.
    """
    from proteostasis import simulate
    kick = cfg.get("attractor_kick", 0.1)
    tol = cfg.get("attractor_rel_tol", 0.05)
    for s in (1.0 + kick, 1.0 - kick):
        tr = simulate(p, eq.u * s, eq.a * s, t_end=cfg.get("attractor_t_end", 5.0e4),
                      n_out=200, blowup=cfg["attractor_hi"] * 10.0)
        if tr.status not in ("converged", "timeout"):
            return False
        back = (abs(tr.final_u - eq.u) <= tol * max(eq.u, 1e-12)
                and abs(tr.final_a - eq.a) <= tol * max(eq.a, 1e-12))
        if not back:
            return False
    return True


def _fraction(df: pd.DataFrame, col: str) -> dict | None:
    """report a robustness fraction together with its denominator.

    the denominator matters: 'holds in 92%' of 400 evaluable draws out of 1000
    is a different statement from 'holds in 92%' of 1000.
    """
    if col not in df:
        return None
    s = df[col].dropna()
    if not len(s):
        return None
    return dict(n_evaluable=int(len(s)), n_true=int(s.astype(bool).sum()),
                fraction=float(s.astype(bool).mean()))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--n-samples", type=int, default=None)
    args = ap.parse_args()

    cfg = loadConfig(args.config)
    if args.n_samples:
        cfg["n_samples"] = args.n_samples
    outdir = Path(args.outdir or cfg["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    base = Params(**cfg.get("base_params", {})).validate()
    ranges = rangesFromConfig(cfg.get("param_ranges"))
    samples = latinHypercube(ranges, cfg["n_samples"], cfg["seed"])
    _CTX.update(cfg=cfg, base=base)

    workers = resolveWorkers(args.workers or cfg.get("workers"))
    print(f"[C] {len(samples)} LHS draws over {len(ranges)} parameters, "
          f"{workers} workers", flush=True)

    rows = parallelMap(_sampleTask, list(enumerate(samples)), workers, progress_every=100)
    df = pd.DataFrame(rows)
    outputs = {"samples.tsv": writeTable(df, outdir / "samples.tsv")}

    ok = df["error"] == ""
    summary = dict(
        experiment="C",
        n_samples=int(len(df)),
        n_errors=int((~ok).sum()),
        error_examples=sorted(set(df.loc[~ok, "error"]))[:5],
        frac_viable_at_j_lo=float(df.get("viable_at_j_lo", pd.Series(dtype=float))
                                  .fillna(False).astype(bool).mean()),
        robustness={
            "C1_fold_exists": _fraction(df, "C1_fold_exists"),
            "C2_collapse_below_ceiling": _fraction(df, "C2_collapse_below_ceiling"),
            "C3_nascent_shrinks_window": _fraction(df, "C3_nascent_shrinks_window"),
            "C4_second_stable_attractor": _fraction(df, "C4_second_stable_attractor"),
            "C4_second_attractor_confirmed": _fraction(df, "C4_second_attractor_confirmed"),
            "C5_slowing_down": _fraction(df, "C5_slowing_down"),
            "C6_interior_optimum": _fraction(df, "C6_interior_optimum"),
        },
        quantiles={
            k: (None if k not in df or not df[k].notna().any() else
                {q: float(df[k].quantile(q)) for q in (0.05, 0.25, 0.5, 0.75, 0.95)})
            for k in ("j_fold", "C2_fold_to_ceiling_ratio", "C3_shrink_ratio",
                      "C5_rate_ratio", "C6_allocation_span", "fold_aggregate_fraction")
        },
    )
    with open(outdir / "summary.json", "w") as fh:
        fh.write(canonicalJson(summary))
        fh.write("\n")
    from proteostasis.provenance import hashFile
    outputs["summary.json"] = hashFile(outdir / "summary.json")
    writeProvenance(outdir, cfg, outputs, extra=dict(runtime_s=time.time() - t0))

    print(canonicalJson(summary))
    print(f"[C] done in {time.time()-t0:.1f}s -> {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
