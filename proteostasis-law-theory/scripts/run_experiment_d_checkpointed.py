#!/usr/bin/env python
"""experiment D -- checkpointed, timeout-bounded recovery runner.

WHY THIS EXISTS. `run_experiment_d.py` evaluates all backgrounds through a single
`Pool.map`, which returns only when every task has returned. one pathological
background -- a stiff parameter draw whose integration never reaches `t_end` in
reasonable time -- therefore holds the entire experiment hostage and leaves zero
partial output on disk. this runner changes ONLY that: how work is scheduled,
checkpointed and bounded in time.

WHAT IS NOT CHANGED. the science is not reimplemented here. the equations, the
LHS sample matrix, the perturbation axes and levels, the burden threshold and
censor, the integrator tolerances, the three nulls and the summary statistics
are all reached by *importing* `run_experiment_d` and calling its own
`_backgroundTask` and `_pairSummary`. there is no second copy of the model to
drift out of step with the first.

TIMEOUT SEMANTICS (important, and deliberately not a failure mode). a background
that exceeds the wall limit is recorded as `unresolved_timeout`. it is NOT a
numerical failure, NOT an error, and NOT counted in `n_errors`. it contributes no
rows to `interactions.tsv` and is excluded from every interaction summary. the
honest reading of an unresolved background is "this draw was not evaluated within
the budget", which is a statement about the budget, not about the model.

usage
  # controller (default): run every background, checkpointing as it goes
  python scripts/run_experiment_d_checkpointed.py \
      --config configs/phase1/experiment_d.json --outdir OUT \
      --concurrency 16 --background-timeout 3600

  # one isolated background, in this process (what the controller spawns)
  python scripts/run_experiment_d_checkpointed.py \
      --config configs/phase1/experiment_d.json --outdir OUT --background-index 7

  # re-merge existing checkpoints without computing anything
  python scripts/run_experiment_d_checkpointed.py \
      --config configs/phase1/experiment_d.json --outdir OUT --merge-only
"""

from __future__ import annotations

import argparse
import json
import os
import signal
import subprocess
import sys
import time
from collections import deque
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from proteostasis.provenance import (REPO_ROOT, canonicalJson, hashFile, hashObject,
                                     loadConfig, writeProvenance, writeTable)
from proteostasis.sweeps import latinHypercube, rangesFromConfig
from proteostasis.model import Params

# the original module IS the definition of the experiment. importing it, rather
# than copying from it, is what makes this runner scientifically equivalent.
import run_experiment_d as original

#: files whose contents determine the NUMBERS. a checkpoint computed under a
#: different version of any of these is not reusable and is rejected on resume.
#: this runner itself is deliberately absent: it schedules work, it does not
#: compute anything, so an orchestration edit must not invalidate finished
#: science. `provenance.py` is absent for the same reason -- it only serialises.
SOURCE_FILES = (
    "scripts/run_experiment_d.py",
    "scripts/proteostasis/__init__.py",
    "scripts/proteostasis/model.py",
    "scripts/proteostasis/simulate.py",
    "scripts/proteostasis/equilibria.py",
    "scripts/proteostasis/sweeps.py",
)

#: seconds between controller polls of its children
POLL_S = 1.0

#: seconds a timed-out background is given to die on SIGTERM before SIGKILL
KILL_GRACE_S = 10.0


# --------------------------------------------------------------------------
# hashing and identity
# --------------------------------------------------------------------------

def sourceHashes(repo: Path = REPO_ROOT) -> dict:
    """sha256 per scientific source file, keyed by repo-relative path."""
    return {rel: hashFile(repo / rel) for rel in SOURCE_FILES}


def combinedSourceHash(repo: Path = REPO_ROOT) -> str:
    return hashObject(sourceHashes(repo))


def sampleMatrix(cfg: dict):
    """regenerate the ORIGINAL LHS matrix.

    identical call sequence to `run_experiment_d.main`, so the returned samples
    are the same objects the original run is evaluating right now.
    """
    base = Params(**cfg.get("base_params", {})).validate()
    ranges = rangesFromConfig(cfg.get("param_ranges"))
    samples = latinHypercube(ranges, cfg["n_backgrounds"], cfg["seed"])
    return base, samples


def runIdentity(cfg: dict, samples: list) -> dict:
    """the four things a checkpoint must agree with to be reusable."""
    return {
        "config_hash": hashObject(cfg),
        "source_hash": combinedSourceHash(),
        "source_files": sourceHashes(),
        "sample_matrix_hash": hashObject(samples),
        "n_backgrounds": int(len(samples)),
    }


# --------------------------------------------------------------------------
# atomic on-disk checkpoints
# --------------------------------------------------------------------------

def _atomicWriteText(path: Path, text: str) -> str:
    """write-to-temp then rename, so a reader never observes a half file."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    with open(tmp, "w") as fh:
        fh.write(text)
        fh.flush()
        os.fsync(fh.fileno())
    os.replace(tmp, path)
    return hashFile(path)


def backgroundDir(outdir: Path, idx: int) -> Path:
    return Path(outdir) / "backgrounds" / f"background_{idx:04d}"


def readCheckpoint(outdir: Path, idx: int, identity: dict, sample: dict) -> dict | None:
    """load a checkpoint only if it is complete, uncorrupted and OURS.

    returns the payload dict, or None with the rejection reason available via
    `checkpointStatus`. every one of the four identity criteria is enforced:
    background index, config hash, source hash, sample hash.
    """
    ok, _reason, payload = checkpointStatus(outdir, idx, identity, sample)
    return payload if ok else None


def checkpointStatus(outdir: Path, idx: int, identity: dict,
                     sample: dict) -> tuple[bool, str, dict | None]:
    d = backgroundDir(outdir, idx)
    done, ck_p = d / "DONE", d / "checkpoint.json"
    rows_p, meta_p = d / "rows.json", d / "meta.json"

    if not done.exists():
        return False, "no completion marker", None
    for p in (ck_p, rows_p, meta_p):
        if not p.exists():
            return False, f"missing {p.name}", None
    try:
        ck = json.loads(ck_p.read_text())
    except Exception as exc:                                      # noqa: BLE001
        return False, f"unreadable checkpoint.json: {type(exc).__name__}", None

    # the marker names the checkpoint it certifies, so a DONE left behind by an
    # earlier run cannot vouch for a checkpoint written later
    if done.read_text().strip() != hashFile(ck_p):
        return False, "completion marker does not match checkpoint.json", None

    if int(ck.get("background_index", -1)) != int(idx):
        return False, "background index mismatch", None
    if ck.get("config_hash") != identity["config_hash"]:
        return False, "config hash mismatch", None
    if ck.get("source_hash") != identity["source_hash"]:
        return False, "source hash mismatch", None
    if ck.get("sample_hash") != hashObject(sample):
        return False, "sample hash mismatch", None

    # payload integrity: the checkpoint records what rows/meta hashed to
    if hashFile(rows_p) != ck.get("rows_sha256"):
        return False, "rows.json hash mismatch", None
    if hashFile(meta_p) != ck.get("meta_sha256"):
        return False, "meta.json hash mismatch", None

    try:
        rows = json.loads(rows_p.read_text())
        meta = json.loads(meta_p.read_text())
    except Exception as exc:                                      # noqa: BLE001
        return False, f"unreadable payload: {type(exc).__name__}", None
    if not isinstance(rows, list) or not isinstance(meta, dict):
        return False, "payload has wrong type", None

    ck = dict(ck, rows=rows, meta=meta)
    return True, "", ck


def writeCheckpoint(outdir: Path, idx: int, identity: dict, sample: dict,
                    rows: list, meta: dict, runtime_s: float) -> Path:
    """persist one background atomically; DONE is written last, and only last."""
    d = backgroundDir(outdir, idx)
    d.mkdir(parents=True, exist_ok=True)

    # canonicalJson round-trips float64 exactly (python repr), unlike the %.12g
    # of writeTable -- the merged table must be bit-comparable to a direct run
    rows_sha = _atomicWriteText(d / "rows.json", canonicalJson(rows) + "\n")
    meta_sha = _atomicWriteText(d / "meta.json", canonicalJson(meta) + "\n")

    ck = {
        "background_index": int(idx),
        "config_hash": identity["config_hash"],
        "source_hash": identity["source_hash"],
        "source_files": identity["source_files"],
        "sample_matrix_hash": identity["sample_matrix_hash"],
        "sample_hash": hashObject(sample),
        "sample": sample,
        "n_rows": int(len(rows)),
        "rows_sha256": rows_sha,
        "meta_sha256": meta_sha,
        "runtime_s": float(runtime_s),
        "python": sys.version.split()[0],
        "numpy": np.__version__,
        "hostname": os.uname().nodename,
    }
    ck_sha = _atomicWriteText(d / "checkpoint.json", json.dumps(ck, indent=2, sort_keys=True) + "\n")
    _atomicWriteText(d / "DONE", ck_sha + "\n")
    return d


def unresolvedPath(outdir: Path, idx: int) -> Path:
    return Path(outdir) / "unresolved" / f"background_{idx:04d}.json"


def writeUnresolved(outdir: Path, idx: int, status: str, detail: dict) -> None:
    """record a background that produced no checkpoint, and why.

    `status` is 'unresolved_timeout' or 'failed_process'. neither is a numerical
    result; both are excluded from interactions.tsv and every summary.
    """
    rec = dict(background_index=int(idx), status=status, **detail)
    _atomicWriteText(unresolvedPath(outdir, idx), json.dumps(rec, indent=2, sort_keys=True) + "\n")


# --------------------------------------------------------------------------
# worker mode -- exactly one background, in this process
# --------------------------------------------------------------------------

def runOneBackground(cfg: dict, outdir: Path, idx: int, force: bool = False) -> int:
    base, samples = sampleMatrix(cfg)
    if not (0 <= idx < len(samples)):
        raise SystemExit(f"background index {idx} outside 0..{len(samples)-1}")
    identity = runIdentity(cfg, samples)
    sample = samples[idx]

    if not force:
        ok, reason, _ = checkpointStatus(outdir, idx, identity, sample)
        if ok:
            print(f"[D/ckpt {idx:04d}] valid checkpoint present, skipping", flush=True)
            return 0
        if reason != "no completion marker":
            print(f"[D/ckpt {idx:04d}] rejecting checkpoint: {reason}", flush=True)

    # this is the ONLY place the science is invoked, and it is the original's
    # own function, given the original's own context dict
    original._CTX.update(cfg=cfg, base=base)
    t0 = time.time()
    rows, meta = original._backgroundTask((idx, sample))
    dt = time.time() - t0

    writeCheckpoint(outdir, idx, identity, sample, rows, meta, dt)
    # a fresh success supersedes any earlier unresolved marker for this index
    unresolvedPath(outdir, idx).unlink(missing_ok=True)
    print(f"[D/ckpt {idx:04d}] usable={meta.get('usable')} rows={len(rows)} "
          f"error={meta.get('error') or '-'} in {dt:.1f}s", flush=True)
    return 0


# --------------------------------------------------------------------------
# merge -- original table and summary semantics over completed checkpoints
# --------------------------------------------------------------------------

def mergeCheckpoints(cfg: dict, outdir: Path, identity: dict, samples: list,
                     extra: dict | None = None) -> dict:
    """assemble interactions.tsv / backgrounds.tsv / summary.json.

    ascending background order reproduces `Pool.map`'s order-preserving output,
    so the merged tables are ordered exactly as the original would have written
    them for the same set of backgrounds.
    """
    outdir = Path(outdir)
    rows, metas = [], []
    completed, rejected, missing = [], {}, []

    for idx in range(len(samples)):
        ok, reason, ck = checkpointStatus(outdir, idx, identity, samples[idx])
        if ok:
            rows += ck["rows"]
            metas.append(ck["meta"])
            completed.append(idx)
        elif reason == "no completion marker":
            missing.append(idx)
        else:
            rejected[idx] = reason
            missing.append(idx)

    unresolved = {}
    udir = outdir / "unresolved"
    if udir.is_dir():
        for p in sorted(udir.glob("background_*.json")):
            try:
                rec = json.loads(p.read_text())
            except Exception:                                     # noqa: BLE001
                continue
            i = int(rec.get("background_index", -1))
            if i in completed:            # superseded by a later success
                continue
            unresolved[i] = rec

    timeouts = sorted(i for i, r in unresolved.items() if r.get("status") == "unresolved_timeout")
    failures = sorted(i for i, r in unresolved.items() if r.get("status") == "failed_process")

    df = pd.DataFrame(rows)
    dfm = pd.DataFrame(metas)
    outputs = {
        "interactions.tsv": writeTable(df, outdir / "interactions.tsv"),
        "backgrounds.tsv": writeTable(dfm, outdir / "backgrounds.tsv"),
    }

    # --- identical to run_experiment_d.main, field for field -----------------
    summary = dict(
        experiment="D",
        n_backgrounds=int(len(dfm)),
        n_usable_backgrounds=int(dfm["usable"].sum()) if "usable" in dfm else 0,
        n_errors=int((dfm["error"] != "").sum()) if "error" in dfm else 0,
        unusable_reasons=(dfm.loc[~dfm["usable"], "reason"].value_counts().to_dict()
                          if "reason" in dfm and (~dfm["usable"]).any() else {}),
        n_cells=int(len(df)),
        burden_threshold=cfg["burden_threshold"],
        by_pair={label: original._pairSummary(df[df["pair"] == label])
                 for label, _, _ in original.PAIRS} if len(df) else {},
        overall=original._pairSummary(df) if len(df) else {},
    )
    # --- recovery bookkeeping, kept OUT of the original fields ---------------
    summary["recovery"] = dict(
        runner="run_experiment_d_checkpointed.py",
        n_backgrounds_requested=int(len(samples)),
        n_backgrounds_completed=len(completed),
        n_unresolved_timeout=len(timeouts),
        unresolved_timeout_backgrounds=timeouts,
        n_failed_process=len(failures),
        failed_process_backgrounds=failures,
        n_missing_checkpoint=len([i for i in missing if i not in unresolved]),
        rejected_checkpoints=rejected,
        config_hash=identity["config_hash"],
        source_hash=identity["source_hash"],
        source_files=identity["source_files"],
        sample_matrix_hash=identity["sample_matrix_hash"],
        timeout_semantics=("a timed-out background is unresolved, not a numerical "
                           "failure: it contributes no rows, is excluded from every "
                           "interaction summary, and is not counted in n_errors"),
        **(extra or {}),
    )

    _atomicWriteText(outdir / "summary.json", canonicalJson(summary) + "\n")
    outputs["summary.json"] = hashFile(outdir / "summary.json")
    writeProvenance(outdir, cfg, outputs,
                    extra=dict(runner="run_experiment_d_checkpointed.py",
                               recovery=summary["recovery"]))
    return summary


# --------------------------------------------------------------------------
# controller -- isolated subprocesses, bounded concurrency, hard wall timeout
# --------------------------------------------------------------------------

def _childEnv(threads: int = 1) -> dict:
    """one BLAS thread per background: 60 backgrounds x N threads oversubscribes."""
    env = dict(os.environ)
    for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[k] = str(threads)
    return env


def _killChildGroup(proc: subprocess.Popen, grace: float = KILL_GRACE_S) -> str:
    """terminate ONLY the recovery child's own process group.

    children are spawned with `start_new_session=True`, so each one leads its own
    session and process group. the identity check below refuses to signal
    anything whose group is not exactly that child -- in particular it can never
    reach the controller's group or any unrelated run on this host.
    """
    try:
        pgid = os.getpgid(proc.pid)
    except ProcessLookupError:
        return "already gone"
    if pgid != proc.pid:
        # not a session leader: refuse to signal a group we do not own
        proc.kill()
        return f"killed pid only (pgid {pgid} != pid {proc.pid})"

    os.killpg(pgid, signal.SIGTERM)
    deadline = time.time() + grace
    while time.time() < deadline:
        if proc.poll() is not None:
            return "SIGTERM"
        time.sleep(0.2)
    try:
        os.killpg(pgid, signal.SIGKILL)
    except ProcessLookupError:
        return "SIGTERM"
    try:
        proc.wait(timeout=grace)
    except subprocess.TimeoutExpired:
        pass
    return "SIGKILL"


def _writeProgress(outdir: Path, state: dict) -> None:
    """flushed on every state change, so an observer always sees current truth."""
    _atomicWriteText(Path(outdir) / "progress.json",
                     json.dumps(state, indent=2, sort_keys=True) + "\n")


def controller(cfg: dict, config_path: str, outdir: Path, concurrency: int,
               timeout_s: float, threads: int, only: list | None = None) -> dict:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(exist_ok=True)

    base, samples = sampleMatrix(cfg)
    identity = runIdentity(cfg, samples)
    n = len(samples)
    # `only` restricts WHICH backgrounds are scheduled; it never changes the LHS
    # matrix, so a selected background is the same draw it would be in a full run
    selected = sorted(set(range(n) if only is None else only))
    for idx in selected:
        if not 0 <= idx < n:
            raise SystemExit(f"background index {idx} outside 0..{n-1}")

    manifest = dict(identity, experiment="D", runner="run_experiment_d_checkpointed.py",
                    config_path=str(Path(config_path).resolve()),
                    outdir=str(outdir.resolve()), concurrency=int(concurrency),
                    background_timeout_s=float(timeout_s), blas_threads=int(threads),
                    selected_backgrounds=selected, started_epoch=time.time(),
                    hostname=os.uname().nodename, python=sys.executable)
    _atomicWriteText(outdir / "run_manifest.json",
                     json.dumps(manifest, indent=2, sort_keys=True) + "\n")

    resumed, pending = [], []
    for idx in selected:
        ok, _r, _c = checkpointStatus(outdir, idx, identity, samples[idx])
        (resumed if ok else pending).append(idx)

    state = dict(n_backgrounds=n, n_selected=len(selected), selected=selected,
                 concurrency=int(concurrency),
                 background_timeout_s=float(timeout_s), resumed=resumed,
                 completed=[], unresolved_timeout=[], failed_process=[],
                 running=[], pending=list(pending), updated_epoch=time.time())
    _writeProgress(outdir, state)
    print(f"[D/ckpt] {len(selected)} of {n} backgrounds selected, "
          f"{len(resumed)} resumed from checkpoint, "
          f"{len(pending)} to run, concurrency {concurrency}, "
          f"timeout {timeout_s:.0f}s/background", flush=True)

    queue, running = deque(pending), {}
    env = _childEnv(threads)
    t0 = time.time()

    def _touch():
        state["running"] = sorted(running)
        state["pending"] = list(queue)
        state["updated_epoch"] = time.time()
        _writeProgress(outdir, state)

    while queue or running:
        while queue and len(running) < concurrency:
            idx = queue.popleft()
            # n_backgrounds is passed explicitly: it sizes the LHS matrix, so a
            # child left to re-read it from the config would draw a DIFFERENT
            # sample whenever the controller was given --n-backgrounds
            cmd = [sys.executable, str(Path(__file__).resolve()),
                   "--config", str(config_path), "--outdir", str(outdir),
                   "--n-backgrounds", str(cfg["n_backgrounds"]),
                   "--background-index", str(idx)]
            fo = open(outdir / "logs" / f"background_{idx:04d}.out", "w")
            fe = open(outdir / "logs" / f"background_{idx:04d}.err", "w")
            proc = subprocess.Popen(cmd, stdout=fo, stderr=fe, env=env,
                                    cwd=str(REPO_ROOT), start_new_session=True)
            running[idx] = dict(proc=proc, started=time.time(), fo=fo, fe=fe)
            _touch()

        time.sleep(POLL_S)

        for idx in sorted(running):
            job = running[idx]
            rc, elapsed = job["proc"].poll(), time.time() - job["started"]

            if rc is None:
                if elapsed <= timeout_s:
                    continue
                how = _killChildGroup(job["proc"])
                writeUnresolved(outdir, idx, "unresolved_timeout",
                                dict(elapsed_s=round(elapsed, 1),
                                     timeout_s=float(timeout_s), killed_by=how,
                                     note="not a numerical failure; excluded from all summaries"))
                state["unresolved_timeout"].append(idx)
                print(f"[D/ckpt {idx:04d}] UNRESOLVED after {elapsed:.0f}s "
                      f"(wall limit {timeout_s:.0f}s, {how})", flush=True)
            else:
                ok, reason, _c = checkpointStatus(outdir, idx, identity, samples[idx])
                if ok:
                    state["completed"].append(idx)
                    print(f"[D/ckpt {idx:04d}] done in {elapsed:.0f}s "
                          f"({len(state['completed'])+len(resumed)}/{len(selected)})", flush=True)
                else:
                    writeUnresolved(outdir, idx, "failed_process",
                                    dict(returncode=int(rc), elapsed_s=round(elapsed, 1),
                                         reason=reason,
                                         stderr_tail=(outdir / "logs" / f"background_{idx:04d}.err")
                                         .read_text()[-2000:]))
                    state["failed_process"].append(idx)
                    print(f"[D/ckpt {idx:04d}] FAILED rc={rc} ({reason})", flush=True)

            job["fo"].close()
            job["fe"].close()
            del running[idx]
            _touch()

    state["running"], state["pending"] = [], []
    state["finished_epoch"] = time.time()
    _writeProgress(outdir, state)

    summary = mergeCheckpoints(cfg, outdir, identity, samples,
                               extra=dict(concurrency=int(concurrency),
                                          background_timeout_s=float(timeout_s),
                                          blas_threads=int(threads),
                                          selected_backgrounds=selected,
                                          wall_s=round(time.time() - t0, 1),
                                          n_resumed_from_checkpoint=len(resumed)))
    return summary


# --------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--config", required=True)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--n-backgrounds", type=int, default=None)
    ap.add_argument("--background-index", type=int, default=None,
                    help="worker mode: compute exactly this background and exit")
    ap.add_argument("--concurrency", type=int, default=None,
                    help="controller mode: simultaneous background subprocesses")
    ap.add_argument("--background-timeout", type=float, default=3600.0,
                    help="controller mode: hard wall limit per background, seconds")
    ap.add_argument("--threads", type=int, default=1,
                    help="BLAS threads per background subprocess")
    ap.add_argument("--only", default=None,
                    help="controller mode: comma-separated background indices to "
                         "schedule (the LHS matrix is unchanged; this selects which "
                         "of its draws are evaluated)")
    ap.add_argument("--merge-only", action="store_true",
                    help="merge existing checkpoints; compute nothing")
    ap.add_argument("--force", action="store_true",
                    help="worker mode: recompute even if a valid checkpoint exists")
    args = ap.parse_args(argv)

    cfg = loadConfig(args.config)
    if args.n_backgrounds:
        cfg["n_backgrounds"] = args.n_backgrounds
    outdir = Path(args.outdir or cfg["outdir"])
    outdir.mkdir(parents=True, exist_ok=True)

    if args.background_index is not None:
        return runOneBackground(cfg, outdir, args.background_index, force=args.force)

    if args.merge_only:
        _base, samples = sampleMatrix(cfg)
        summary = mergeCheckpoints(cfg, outdir, runIdentity(cfg, samples), samples,
                                   extra=dict(merge_only=True))
        print(canonicalJson(summary))
        return 0

    concurrency = max(1, int(args.concurrency or cfg.get("workers") or 1))
    only = ([int(x) for x in args.only.split(",") if x.strip() != ""]
            if args.only else None)
    t0 = time.time()
    summary = controller(cfg, args.config, outdir, concurrency,
                         args.background_timeout, args.threads, only=only)
    print(canonicalJson(summary))
    rec = summary["recovery"]
    print(f"[D/ckpt] done in {time.time()-t0:.1f}s -> {outdir} "
          f"({rec['n_backgrounds_completed']} completed, "
          f"{rec['n_unresolved_timeout']} unresolved_timeout, "
          f"{rec['n_failed_process']} failed)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
