#!/usr/bin/env python
"""phase 2 step 1 -- verify that experiment C completed cleanly, and reproduce
every headline count in `summary.json` directly from `samples.tsv`.

what this checks, in order of how badly it would matter if it failed:

  1. the two output files still hash to what `provenance.json` recorded, so the
     table being audited is the table the run produced;
  2. the recorded config hashes to its own recorded hash, and to the config
     still on disk at `configs/phase1/experiment_c.json`;
  3. the run reported no per-draw errors, and the table agrees;
  4. every robustness fraction and denominator in `summary.json` is reproduced
     from the raw table by an independently written counter;
  5. the parameter reconstruction used by the rest of phase 2 is correct, proven
     against the `removal_ceiling` column the run itself stored.

usage:  PHASE2_AUDIT_DIR=results/phase2/audit_... python scripts/phase2/c01_verify_run.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _audit_common import (auditRoot, checkReconstruction, loadProvenance,  # noqa: E402
                           loadSamples, runRoot, sha256File, writeJson)

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from proteostasis.provenance import hashObject, loadConfig  # noqa: E402

#: the boolean columns whose fractions summary.json reports.
ROBUSTNESS_COLS = ("C1_fold_exists", "C2_collapse_below_ceiling",
                   "C3_nascent_shrinks_window", "C4_second_stable_attractor",
                   "C4_second_attractor_confirmed", "C5_slowing_down",
                   "C6_interior_optimum")


def recountFraction(df: pd.DataFrame, col: str):
    """independent reimplementation of run_experiment_c._fraction.

    written from the definition ('of the draws where this question was
    evaluable, in how many was the answer yes') rather than by copying the
    original, so a bug shared with the original cannot cancel out.
    """
    if col not in df.columns:
        return None
    s = df[col].dropna()
    if len(s) == 0:
        return None
    truth = s.map(lambda v: bool(v) if not isinstance(v, str)
                  else v.strip().lower() in ("true", "1", "yes"))
    return dict(n_evaluable=int(len(s)), n_true=int(truth.sum()),
                fraction=float(truth.mean()))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=None)
    ap.add_argument("--audit-dir", default=None)
    args = ap.parse_args()

    run = runRoot(args.run)
    out = auditRoot(args.audit_dir)
    report: dict = {"run": str(run), "step": "c01_verify_run"}

    # --- 1. completion and checksums --------------------------------------
    cdir = run / "C"
    prov = loadProvenance(run)
    files = {}
    for name, recorded in prov["output_sha256"].items():
        path = cdir / name
        files[name] = dict(exists=path.exists(),
                           recorded_sha256=recorded,
                           actual_sha256=sha256File(path) if path.exists() else None,
                           size_bytes=path.stat().st_size if path.exists() else None)
        files[name]["matches"] = files[name]["actual_sha256"] == recorded
    report["output_files"] = files
    report["all_checksums_match"] = all(f["matches"] for f in files.values())

    # --- 2. config provenance ---------------------------------------------
    disk_cfg = loadConfig(str(Path(__file__).resolve().parents[2]
                              / "configs" / "phase1" / "experiment_c.json"))
    report["config"] = dict(
        recorded_hash=prov["config_hash"],
        recorded_config_rehashes_to_recorded_hash=(
            hashObject(prov["config"]) == prov["config_hash"]),
        disk_config_hash=hashObject(disk_cfg),
        disk_config_matches_run=(hashObject(disk_cfg) == prov["config_hash"]),
    )
    report["git"] = prov["git"]
    report["environment"] = prov["environment"]
    report["runtime_s"] = prov["extra"]["runtime_s"]

    # --- 3. completion markers --------------------------------------------
    stdout = (run / "ops" / "experiment_C.stdout.log")
    stderr = (run / "ops" / "experiment_C.stderr.log")
    stdout_text = stdout.read_text() if stdout.exists() else ""
    report["process"] = dict(
        stdout_bytes=stdout.stat().st_size if stdout.exists() else None,
        stderr_bytes=stderr.stat().st_size if stderr.exists() else None,
        stderr_empty=(stderr.exists() and stderr.stat().st_size == 0),
        stdout_has_done_line=("[C] done in" in stdout_text),
    )

    # --- 4. reproduce the summary from the raw table ----------------------
    df = loadSamples(run)
    with open(cdir / "summary.json") as fh:
        summary = json.load(fh)

    errs = df["error"].astype(str).str.strip()
    n_err = int((errs != "").sum())
    recomputed = dict(
        n_samples=int(len(df)),
        n_errors=n_err,
        error_examples=sorted(set(errs[errs != ""]))[:5],
        frac_viable_at_j_lo=float(df["viable_at_j_lo"].fillna(False)
                                  .astype(bool).mean()),
        robustness={c: recountFraction(df, c) for c in ROBUSTNESS_COLS},
    )
    mismatches = []
    for key in ("n_samples", "n_errors", "frac_viable_at_j_lo"):
        if summary.get(key) != recomputed[key]:
            mismatches.append(dict(field=key, summary=summary.get(key),
                                   recomputed=recomputed[key]))
    for c in ROBUSTNESS_COLS:
        s, r = summary["robustness"].get(c), recomputed["robustness"][c]
        if s is None or r is None:
            if s is not r:
                mismatches.append(dict(field=c, summary=s, recomputed=r))
            continue
        for k in ("n_evaluable", "n_true"):
            if int(s[k]) != int(r[k]):
                mismatches.append(dict(field=f"{c}.{k}", summary=s[k], recomputed=r[k]))
        if abs(float(s["fraction"]) - r["fraction"]) > 1e-12:
            mismatches.append(dict(field=f"{c}.fraction", summary=s["fraction"],
                                   recomputed=r["fraction"]))
    report["summary_recomputed"] = recomputed
    report["summary_mismatches"] = mismatches
    report["summary_reproduces"] = (len(mismatches) == 0)

    # --- structural sanity of the table -----------------------------------
    report["table"] = dict(
        n_rows=int(len(df)),
        n_cols=int(df.shape[1]),
        sample_index_unique=bool(df["sample_index"].is_unique),
        sample_index_contiguous=bool(
            df["sample_index"].sort_values().tolist() == list(range(len(df)))),
        sample_index_ordered=bool(df["sample_index"].is_monotonic_increasing),
        n_viable_at_j_lo=int(df["viable_at_j_lo"].fillna(False).astype(bool).sum()),
        n_fold_found=int(df["C1_fold_exists"].fillna(False).astype(bool).sum()),
        n_j_fold_finite=int(np.isfinite(df["j_fold"]).sum()),
        max_stable_equilibria_counts={
            str(k): int(v) for k, v in
            df["max_stable_equilibria"].value_counts(dropna=False).items()},
        max_equilibria_counts={
            str(k): int(v) for k, v in
            df["max_equilibria"].value_counts(dropna=False).items()},
        fold_reason_counts={str(k): int(v) for k, v in
                            df["fold_reason"].value_counts(dropna=False).items()},
    )

    # the C4 denominator must be exactly the fold-evaluable subset
    c4_eval = df["C4_second_stable_attractor"].notna()
    report["table"]["c4_denominator_equals_fold_subset"] = bool(
        (c4_eval == df["C1_fold_exists"].fillna(False).astype(bool)).all())

    # --- 5. parameter reconstruction --------------------------------------
    report["reconstruction"] = checkReconstruction(
        df, prov["config"]["base_params"], prov["config"].get("rescue_total", 1.0))

    # --- candidate list for the next steps --------------------------------
    cand = df[df["max_stable_equilibria"].fillna(0) >= 2]
    report["candidates"] = dict(
        n_with_two_or_more_stable=int(len(cand)),
        n_confirmed_by_phase1=int(df["C4_second_attractor_confirmed"]
                                  .fillna(False).astype(bool).sum()),
        sample_indices=[int(i) for i in cand["sample_index"].tolist()],
    )

    report["verdict"] = dict(
        completed_cleanly=bool(report["process"]["stderr_empty"]
                               and report["process"]["stdout_has_done_line"]
                               and n_err == 0),
        provenance_intact=bool(report["all_checksums_match"]
                               and report["config"]["disk_config_matches_run"]),
        summary_reproduces=report["summary_reproduces"],
    )

    writeJson(report, out / "c01_verify_run.json")
    print(json.dumps(report["verdict"], indent=2))
    print(f"candidates (>=2 stable): {report['candidates']['n_with_two_or_more_stable']}")
    print(f"phase1-confirmed:        {report['candidates']['n_confirmed_by_phase1']}")
    print(f"mismatches:              {len(mismatches)}")
    print(f"-> {out/'c01_verify_run.json'}")
    return 0 if all(report["verdict"].values()) else 1


if __name__ == "__main__":
    raise SystemExit(main())
