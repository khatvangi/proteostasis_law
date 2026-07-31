#!/usr/bin/env python
"""phase 2 step 7 -- record the exact verified state of experiment D, and audit
it independently if and only if it has finished.

this script is strictly read-only with respect to experiment D. it inspects the
systemd transient unit, the live process tree, the output directory and the
logs. it never signals, restarts or duplicates the run, and it never waits.

one structural fact worth knowing before reading its output: `run_experiment_d`
collects its results with `parallelMap`, which is `Pool.map` -- every background
must finish before ANY row is written. so an empty output directory is the
expected state at every moment before completion and carries no information
about progress. the only observable progress signals are accumulated CPU time
and the worker count, and this script records both rather than guessing.

usage:  PHASE2_AUDIT_DIR=... python scripts/phase2/d01_status.py
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _audit_common import auditRoot, runRoot, sha256File, writeJson  # noqa: E402

UNIT = "proteostasis-phase1-20260731T162946-0500-D.service"

#: the three perturbation pairs, restated here rather than imported, so the
#: audit does not inherit a relabelling from the code it is auditing.
PAIRS = ("influx_x_total_capacity", "influx_x_chaperone_only",
         "nascent_x_total_capacity")


def _sh(*args) -> str:
    try:
        return subprocess.run(args, capture_output=True, text=True,
                              timeout=30).stdout.strip()
    except Exception as exc:                                      # noqa: BLE001
        return f"<failed: {type(exc).__name__}: {exc}>"


def systemdState(unit: str) -> dict:
    props = ("Id", "Description", "ActiveState", "SubState", "Result", "MainPID",
             "ExecMainPID", "ExecMainStatus", "ExecMainStartTimestamp",
             "ExecMainExitTimestamp", "CPUUsageNSec", "MemoryCurrent",
             "MemoryPeak", "TasksCurrent", "NRestarts")
    txt = _sh("systemctl", "--user", "show", unit, "--no-pager",
              *[f"--property={k}" for k in props])
    out = {}
    for line in txt.splitlines():
        if "=" in line:
            k, v = line.split("=", 1)
            out[k] = v
    return out


def processTree(pid: int) -> dict:
    if pid <= 0:
        return dict(alive=False)
    p = Path(f"/proc/{pid}")
    if not p.exists():
        return dict(alive=False)
    try:
        status = dict(
            line.split(":", 1) for line in
            (p / "status").read_text().splitlines() if ":" in line)
        status = {k: v.strip() for k, v in status.items()}
    except OSError:
        return dict(alive=False)
    kids = _sh("pgrep", "-P", str(pid)).split()
    child = []
    for c in kids:
        try:
            child.append(dict(
                pid=int(c),
                utime_s=float(Path(f"/proc/{c}/stat").read_text().split()[13]) / 100.0,
                rss_kb=int(Path(f"/proc/{c}/status").read_text()
                           .split("VmRSS:")[1].split()[0])
                if "VmRSS:" in Path(f"/proc/{c}/status").read_text() else None))
        except (OSError, IndexError, ValueError):
            child.append(dict(pid=int(c), utime_s=None, rss_kb=None))
    return dict(
        alive=True, pid=pid, name=status.get("Name"), state=status.get("State"),
        ppid=status.get("PPid"), threads=status.get("Threads"),
        vm_rss=status.get("VmRSS"),
        cmdline=(p / "cmdline").read_text().replace("\x00", " ").strip(),
        n_children=len(kids), children=child,
        elapsed_s=float(_sh("ps", "-o", "etimes=", "-p", str(pid)) or 0),
    )


# ---------------------------------------------------------------------------
# independent recomputation, used only once D has finished
# ---------------------------------------------------------------------------


def pairSummary(d: pd.DataFrame) -> dict:
    """reimplementation of run_experiment_d._pairSummary from its definitions.

    written from the docstring's statement of each null, not by copying the
    original, so a sign error shared with the original cannot cancel out.
    """
    d = d[~d["single_perturbation"].astype(bool)]
    if not len(d):
        return {}
    out = dict(n_double_cells=int(len(d)),
               n_backgrounds=int(d["background"].nunique()),
               frac_censored=float(d["censored_12"].astype(bool).mean()),
               frac_synthetic_collapse=float(d["synthetic_collapse"].astype(bool).mean()),
               frac_backgrounds_with_collapse=float(
                   d.groupby("background")["synthetic_collapse"]
                   .apply(lambda s: bool(s.astype(bool).any())).mean()))
    for null, col, worse_is in (("additive", "excess_additive", "greater"),
                                ("multiplicative", "log_excess_multiplicative", "greater"),
                                ("bliss", "excess_bliss", "less")):
        s = d[col].replace([np.inf, -np.inf], np.nan).dropna()
        if not len(s):
            continue
        out[f"{null}_median_excess"] = float(s.median())
        out[f"{null}_frac_worse_than_null"] = float(
            (s > 0).mean() if worse_is == "greater" else (s < 0).mean())
        out[f"{null}_n"] = int(len(s))
    return out


def recomputeNulls(d: pd.DataFrame) -> dict:
    """rebuild the excess columns themselves from the raw burden/survival columns.

    checking the summary against the stored excess columns would only prove the
    aggregation is right. the nulls are where a sign or scale error would live,
    so they are recomputed from burden_0/1/2/12 and survival_0/1/2/12.
    """
    b0, b1, b2, b12 = (d[f"burden_{k}"].astype(float) for k in ("0", "1", "2", "12"))
    s0, s1, s2, s12 = (d[f"survival_{k}"].astype(float) for k in ("0", "1", "2", "12"))
    add = b1 + b2 - b0
    mult = b0 * (b1 / b0) * (b2 / b0)
    bliss = s0.clip(lower=1e-12) * (s1 / s0.clip(lower=1e-12)) * (s2 / s0.clip(lower=1e-12))
    def worst(a, b):
        return float(np.nanmax(np.abs(np.asarray(a, float) - np.asarray(b, float))))
    return dict(
        max_abs_diff_excess_additive=worst(d["excess_additive"], b12 - add),
        max_abs_diff_excess_multiplicative=worst(d["excess_multiplicative"], b12 - mult),
        max_abs_diff_log_excess_multiplicative=worst(
            d["log_excess_multiplicative"],
            np.log(b12.clip(lower=1e-300)) - np.log(mult.clip(lower=1e-300))),
        max_abs_diff_excess_bliss=worst(d["excess_bliss"], s12 - bliss),
        synthetic_collapse_matches=bool(
            (d["synthetic_collapse"].astype(bool)
             == (d["viable_1"].astype(bool) & d["viable_2"].astype(bool)
                 & ~d["viable_12"].astype(bool))).all()),
        single_perturbation_matches=bool(
            (d["single_perturbation"].astype(bool)
             == ((d["burden_factor"] == 1.0) | (d["capacity_factor"] == 1.0))).all()),
    )


def auditCompleted(ddir: Path) -> dict:
    with open(ddir / "provenance.json") as fh:
        prov = json.load(fh)
    files = {}
    for name, recorded in prov["output_sha256"].items():
        path = ddir / name
        files[name] = dict(exists=path.exists(), recorded_sha256=recorded,
                           actual_sha256=sha256File(path) if path.exists() else None)
        files[name]["matches"] = files[name]["actual_sha256"] == recorded

    df = pd.read_csv(ddir / "interactions.tsv", sep="\t")
    dfm = pd.read_csv(ddir / "backgrounds.tsv", sep="\t")
    with open(ddir / "summary.json") as fh:
        summary = json.load(fh)

    recomputed = dict(
        experiment="D",
        n_backgrounds=int(len(dfm)),
        n_usable_backgrounds=int(dfm["usable"].astype(bool).sum()),
        n_errors=int((dfm["error"].fillna("").astype(str).str.strip() != "").sum()),
        n_cells=int(len(df)),
        by_pair={label: pairSummary(df[df["pair"] == label]) for label in PAIRS},
        overall=pairSummary(df),
    )
    mism = []
    for k in ("n_backgrounds", "n_usable_backgrounds", "n_errors", "n_cells"):
        if summary.get(k) != recomputed[k]:
            mism.append(dict(field=k, summary=summary.get(k), recomputed=recomputed[k]))
    for scope in ("overall",):
        for k, v in recomputed[scope].items():
            sv = summary.get(scope, {}).get(k)
            if sv is None or (isinstance(v, float) and not np.isclose(sv, v, rtol=1e-9,
                                                                     atol=1e-12)) \
               or (isinstance(v, int) and sv != v):
                mism.append(dict(field=f"{scope}.{k}", summary=sv, recomputed=v))
    for label in PAIRS:
        for k, v in recomputed["by_pair"][label].items():
            sv = summary.get("by_pair", {}).get(label, {}).get(k)
            if sv is None or (isinstance(v, float) and not np.isclose(sv, v, rtol=1e-9,
                                                                     atol=1e-12)) \
               or (isinstance(v, int) and sv != v):
                mism.append(dict(field=f"by_pair.{label}.{k}", summary=sv, recomputed=v))

    return dict(output_files=files,
                all_checksums_match=all(f["matches"] for f in files.values()),
                config_hash=prov["config_hash"], git=prov["git"],
                environment=prov["environment"], runtime_s=prov["extra"]["runtime_s"],
                summary_recomputed=recomputed, summary_mismatches=mism,
                summary_reproduces=(len(mism) == 0),
                null_recomputation=recomputeNulls(df),
                censoring=dict(
                    frac_cells_censored=float(df["censored_12"].astype(bool).mean()),
                    status_counts={str(k): int(v) for k, v in
                                   df["status_12"].value_counts().items()}),
                unusable_reasons={str(k): int(v) for k, v in
                                  dfm.loc[~dfm["usable"].astype(bool), "reason"]
                                  .value_counts().items()}
                if "reason" in dfm else {})


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=None)
    ap.add_argument("--audit-dir", default=None)
    ap.add_argument("--unit", default=UNIT)
    ap.add_argument("--out", default="d01_status.json")
    args = ap.parse_args()

    run = runRoot(args.run)
    out = auditRoot(args.audit_dir)
    ddir = run / "D"

    sd = systemdState(args.unit)
    pid = int(sd.get("MainPID", "0") or 0)
    state = dict(
        observed_at=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        unit=args.unit, systemd=sd, process=processTree(pid),
        outdir=str(ddir),
        outdir_contents=sorted(
            dict(name=q.name, size_bytes=q.stat().st_size,
                 mtime=time.strftime("%Y-%m-%dT%H:%M:%S%z",
                                     time.localtime(q.stat().st_mtime)))
            for q in ddir.iterdir()) if ddir.exists() else None,
    )
    for tag in ("stdout", "stderr"):
        f = run / "ops" / f"experiment_D.{tag}.log"
        state[tag] = dict(exists=f.exists(),
                          size_bytes=f.stat().st_size if f.exists() else None,
                          text=f.read_text()[-4000:] if f.exists() else None)

    expected = ("interactions.tsv", "backgrounds.tsv", "summary.json",
                "provenance.json")
    present = {q.name for q in ddir.iterdir()} if ddir.exists() else set()
    complete = (all(e in present for e in expected)
                and sd.get("ActiveState") in ("inactive", "failed")
                and "[D] done in" in (state["stdout"]["text"] or ""))
    state["expected_outputs"] = list(expected)
    state["present_outputs"] = sorted(present)
    state["complete"] = bool(complete)
    state["still_running"] = bool(sd.get("SubState") == "running"
                                  and state["process"].get("alive"))
    if complete:
        state["audit"] = auditCompleted(ddir)
    else:
        state["audit"] = None
        state["note"] = (
            "experiment D has not produced outputs. run_experiment_d collects "
            "results with Pool.map, so nothing is written until every background "
            "finishes; an empty output directory is therefore the expected state "
            "throughout the run and is NOT evidence of a problem. no result, "
            "partial or otherwise, exists to audit.")

    writeJson(state, out / args.out)
    print(json.dumps({k: state[k] for k in
                      ("observed_at", "complete", "still_running",
                       "present_outputs")}, indent=2))
    print(json.dumps({k: sd.get(k) for k in
                      ("ActiveState", "SubState", "MainPID", "CPUUsageNSec",
                       "TasksCurrent", "ExecMainStartTimestamp")}, indent=2))
    if state["audit"]:
        print(json.dumps({k: state["audit"][k] for k in
                          ("all_checksums_match", "summary_reproduces")}, indent=2))
    print(f"-> {out/args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
