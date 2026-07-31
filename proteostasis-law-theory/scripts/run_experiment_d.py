#!/usr/bin/env python
"""experiment D -- perturbation interaction study.

tests PREDICTIONS.md #1 (burden-by-capacity synthetic interaction) INSIDE the
model, i.e. asks whether the model actually produces the interaction the theory
advertises, and against which null.

design. each background is a parameter draw whose baseline influx is set to a
fixed fraction of its own fold, so every background starts with comparable
headroom rather than comparable absolute influx. from the unperturbed steady
state -- the biologically relevant initial condition for a perturbation
experiment -- the perturbed system is integrated and scored.

three perturbation pairs, each a full factorial:

  1. influx burden  x  total rescue capacity
  2. influx burden  x  chaperone-only knockdown   (degradation untouched)
  3. nascent load   x  total rescue capacity      (a burden that carries NO
                                                   damage influx at all)

three nulls, because the answer depends on the scale the null is stated on and
the theory does not privilege one:

  additive        db_12 == db_1 + db_2                     (burden scale)
  multiplicative  b_12/b_0 == (b_1/b_0)(b_2/b_0)           (log-burden scale)
  bliss           s_12 == s_1 * s_2                        (survival scale,
                  s = clip(1 - burden/H, 0, 1) normalised to the control)

synthetic collapse is recorded separately from any null: both single
perturbations survivable, the combination not.

usage:  python scripts/run_experiment_d.py --config configs/phase1/experiment_d.json
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

from proteostasis import (Params, ModelError, findFold, lowestStableEquilibrium,
                          simulate)
from proteostasis.equilibria import stableEquilibriumAt
from proteostasis.provenance import loadConfig, writeTable, writeProvenance, canonicalJson
from proteostasis.sweeps import (latinHypercube, paramsFromSample, parallelMap,
                                 rangesFromConfig, resolveWorkers)

_CTX: dict = {}

#: (label, burden axis, capacity axis). the capacity axis multiplies a
#: conserved pool, so a knockdown is a genuine reduction of total machinery,
#: not a rate-constant change.
PAIRS = (
    ("influx_x_total_capacity", "j", "rescue_total"),
    ("influx_x_chaperone_only", "j", "c_tot"),
    ("nascent_x_total_capacity", "nu", "rescue_total"),
)


def _perturb(p: Params, axis: str, factor: float) -> Params:
    if axis == "j":
        return p.with_(j=p.j * factor)
    if axis == "nu":
        return p.with_(nu=p.nu * factor)
    if axis == "c_tot":
        return p.with_(c_tot=p.c_tot * factor)
    if axis == "d_tot":
        return p.with_(d_tot=p.d_tot * factor)
    if axis == "rescue_total":
        return p.with_(c_tot=p.c_tot * factor, d_tot=p.d_tot * factor)
    raise ModelError(f"unknown perturbation axis '{axis}'")


def _score(p: Params, u0: float, a0: float, cfg: dict) -> dict:
    """integrate the perturbed system from the pre-perturbation steady state."""
    H = cfg["burden_threshold"]
    b_max = cfg["burden_censor"]
    tr = simulate(p, u0, a0, t_end=cfg["t_end"], n_out=cfg["n_out"],
                  blowup=b_max, rtol=cfg.get("rtol", 1e-8), atol=cfg.get("atol", 1e-11))
    burden = tr.final_u + tr.final_a
    escaped = tr.status in ("blowup", "error")
    # censor rather than record inf, so the null arithmetic stays finite; the
    # censoring is reported so a reader can see how many cells hit the ceiling
    burden_c = float(min(burden, b_max)) if np.isfinite(burden) else b_max
    if escaped:
        burden_c = b_max
    viable = bool((not escaped) and burden_c < H)
    survival = float(np.clip(1.0 - burden_c / H, 0.0, 1.0))
    return dict(status=tr.status, burden=burden_c, censored=bool(burden >= b_max or escaped),
                viable=viable, survival=survival, peak_burden=float(np.max(tr.u + tr.a)))


def _backgroundTask(item: tuple):
    """everything for one parameter background. module level so it is picklable."""
    idx, sample = item
    cfg = _CTX["cfg"]
    base = _CTX["base"]
    rows, meta = [], dict(background=int(idx), error="", usable=False)
    meta.update({f"p_{k}": float(v) for k, v in sample.items()})

    try:
        p_raw = paramsFromSample(sample, base, cfg.get("rescue_total", 1.0))
        if lowestStableEquilibrium(p_raw.with_(j=cfg["j_lo"]), n_grid=9) is None:
            meta["reason"] = "no viable state at j_lo"
            return rows, meta
        fold = findFold(p_raw, "j", cfg["j_lo"], cfg["j_hi"], **cfg["fold_full"])
        if fold.fold_value is None:
            meta["reason"] = f"no fold: {fold.reason}"
            return rows, meta

        j0 = cfg["baseline_fraction"] * fold.fold_value
        p0 = p_raw.with_(j=j0)
        eq0 = stableEquilibriumAt(p0, "j", j0, None, blind=True)
        if eq0 is None:
            meta["reason"] = "baseline equilibrium not found"
            return rows, meta

        base_score = _score(p0, eq0.u, eq0.a, cfg)
        if not base_score["viable"]:
            meta["reason"] = "baseline not viable below threshold"
            return rows, meta

        meta.update(usable=True, j_fold=fold.fold_value, j_baseline=j0,
                    baseline_burden=base_score["burden"],
                    baseline_survival=base_score["survival"],
                    baseline_u=eq0.u, baseline_a=eq0.a)

        burden_levels = np.array(cfg["burden_levels"], dtype=float)
        capacity_levels = np.array(cfg["capacity_levels"], dtype=float)

        for label, b_axis, c_axis in PAIRS:
            # singles first, so the nulls can be assembled without re-simulating
            singles_b = {float(fb): _score(_perturb(p0, b_axis, float(fb)),
                                           eq0.u, eq0.a, cfg) for fb in burden_levels}
            singles_c = {float(fc): _score(_perturb(p0, c_axis, float(fc)),
                                           eq0.u, eq0.a, cfg) for fc in capacity_levels}
            for fb in burden_levels:
                for fc in capacity_levels:
                    pp = _perturb(_perturb(p0, b_axis, float(fb)), c_axis, float(fc))
                    s12 = _score(pp, eq0.u, eq0.a, cfg)
                    s1, s2 = singles_b[float(fb)], singles_c[float(fc)]
                    b0 = base_score["burden"]
                    r = dict(background=int(idx), pair=label,
                             burden_axis=b_axis, capacity_axis=c_axis,
                             burden_factor=float(fb), capacity_factor=float(fc),
                             burden_0=b0, burden_1=s1["burden"], burden_2=s2["burden"],
                             burden_12=s12["burden"],
                             viable_0=base_score["viable"], viable_1=s1["viable"],
                             viable_2=s2["viable"], viable_12=s12["viable"],
                             survival_0=base_score["survival"], survival_1=s1["survival"],
                             survival_2=s2["survival"], survival_12=s12["survival"],
                             censored_12=s12["censored"], status_12=s12["status"])

                    # --- nulls ------------------------------------------
                    add_pred = s1["burden"] + s2["burden"] - b0
                    r["null_additive_pred"] = add_pred
                    r["excess_additive"] = s12["burden"] - add_pred

                    mult_pred = b0 * (s1["burden"] / b0) * (s2["burden"] / b0)
                    r["null_multiplicative_pred"] = mult_pred
                    r["excess_multiplicative"] = s12["burden"] - mult_pred
                    r["log_excess_multiplicative"] = (
                        np.log(max(s12["burden"], 1e-300)) - np.log(max(mult_pred, 1e-300)))

                    s0 = max(base_score["survival"], 1e-12)
                    bliss_pred = s0 * (s1["survival"] / s0) * (s2["survival"] / s0)
                    r["null_bliss_pred"] = bliss_pred
                    r["excess_bliss"] = s12["survival"] - bliss_pred   # negative = worse

                    r["synthetic_collapse"] = bool(
                        s1["viable"] and s2["viable"] and not s12["viable"])
                    r["single_perturbation"] = bool(fb == 1.0 or fc == 1.0)
                    rows.append(r)
    except (ModelError, ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
        meta["error"] = f"{type(exc).__name__}: {exc}"
    except Exception as exc:                                  # noqa: BLE001
        meta["error"] = f"UNEXPECTED {type(exc).__name__}: {exc}\n{traceback.format_exc()[:400]}"
    return rows, meta


def _pairSummary(df: pd.DataFrame) -> dict:
    """interaction statistics restricted to genuinely double perturbations."""
    d = df[~df["single_perturbation"]]
    if not len(d):
        return {}
    out = dict(n_double_cells=int(len(d)),
               n_backgrounds=int(d["background"].nunique()),
               frac_censored=float(d["censored_12"].mean()),
               frac_synthetic_collapse=float(d["synthetic_collapse"].mean()),
               frac_backgrounds_with_collapse=float(
                   d.groupby("background")["synthetic_collapse"].any().mean()))
    for null, col, worse in (("additive", "excess_additive", "greater"),
                             ("multiplicative", "log_excess_multiplicative", "greater"),
                             ("bliss", "excess_bliss", "less")):
        s = d[col].replace([np.inf, -np.inf], np.nan).dropna()
        if not len(s):
            continue
        frac_worse = float((s > 0).mean()) if worse == "greater" else float((s < 0).mean())
        out[f"{null}_median_excess"] = float(s.median())
        out[f"{null}_frac_worse_than_null"] = frac_worse
        out[f"{null}_n"] = int(len(s))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--workers", type=int, default=None)
    ap.add_argument("--n-backgrounds", type=int, default=None)
    args = ap.parse_args()

    cfg = loadConfig(args.config)
    if args.n_backgrounds:
        cfg["n_backgrounds"] = args.n_backgrounds
    outdir = Path(args.outdir or cfg["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    base = Params(**cfg.get("base_params", {})).validate()
    ranges = rangesFromConfig(cfg.get("param_ranges"))
    samples = latinHypercube(ranges, cfg["n_backgrounds"], cfg["seed"])
    _CTX.update(cfg=cfg, base=base)

    workers = resolveWorkers(args.workers or cfg.get("workers"))
    print(f"[D] {len(samples)} backgrounds x {len(PAIRS)} pairs x "
          f"{len(cfg['burden_levels'])}x{len(cfg['capacity_levels'])} cells, "
          f"{workers} workers", flush=True)

    results = parallelMap(_backgroundTask, list(enumerate(samples)), workers,
                          progress_every=10)
    rows, metas = [], []
    for r, m in results:
        rows += r
        metas.append(m)

    df = pd.DataFrame(rows)
    dfm = pd.DataFrame(metas)
    outputs = {
        "interactions.tsv": writeTable(df, outdir / "interactions.tsv"),
        "backgrounds.tsv": writeTable(dfm, outdir / "backgrounds.tsv"),
    }

    summary = dict(
        experiment="D",
        n_backgrounds=int(len(dfm)),
        n_usable_backgrounds=int(dfm["usable"].sum()) if "usable" in dfm else 0,
        n_errors=int((dfm["error"] != "").sum()) if "error" in dfm else 0,
        unusable_reasons=(dfm.loc[~dfm["usable"], "reason"].value_counts().to_dict()
                          if "reason" in dfm and (~dfm["usable"]).any() else {}),
        n_cells=int(len(df)),
        burden_threshold=cfg["burden_threshold"],
        by_pair={label: _pairSummary(df[df["pair"] == label]) for label, _, _ in PAIRS}
        if len(df) else {},
        overall=_pairSummary(df) if len(df) else {},
    )
    with open(outdir / "summary.json", "w") as fh:
        fh.write(canonicalJson(summary))
        fh.write("\n")
    from proteostasis.provenance import hashFile
    outputs["summary.json"] = hashFile(outdir / "summary.json")
    writeProvenance(outdir, cfg, outputs, extra=dict(runtime_s=time.time() - t0))

    print(canonicalJson(summary))
    print(f"[D] done in {time.time()-t0:.1f}s -> {outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
