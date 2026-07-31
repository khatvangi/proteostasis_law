# Phase 1 jobs

Single clean submission of the Phase 1 computational falsification run for the
Proteostasis Law. This file records what was actually launched on the host, not
what any result means. **No scientific claim is made here.**

- Run timestamp: `20260731T162946-0500`
- Host: `boron` (Linux 5.14.0-611.16.1.el9_7.x86_64, AlmaLinux, 64 CPUs)
- User: `kiran` (uid 1000)
- Run root: `results/phase1/run_20260731T162946-0500/`
- Host-level ops directory: `results/phase1/run_20260731T162946-0500/ops/`
- Operational manifest: `results/phase1/ops/clean_submission_manifest.json`

## Gate: canonical test evidence (checked before launch)

The launch was gated on the repaired canonical suite. Both conditions had to
hold or the run would not have started.

- `results/phase1/ops/full_tests_after_repair.exitcode` = `0`
- `results/phase1/ops/full_tests_after_repair.log` reports
  `69 passed, 19 warnings in 66.36s (0:01:06)`, from `collected 69 items`
- No `failed` or `error` token anywhere in the log
- Python 3.12.11, pytest 8.4.2

The 19 warnings are the pre-existing multiprocessing `fork()` DeprecationWarning
and are not failures. The suite was repaired earlier by correcting the
saddle-direction and event-resolved slow-mode recovery tests; theory, configs,
and manuscript were not modified.

## Environment

- Python (absolute): `/home/kiran/miniforge3/bin/python` -> 3.12.11
- numpy 2.2.6, scipy 1.14.0
- Thread limits applied to every job and **verified inside each running process**
  via `/proc/<pid>/environ`:
  `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`,
  `NUMEXPR_NUM_THREADS=1`
- Working directory for all units: repository root

## Launch mechanism

`systemd-run --user` transient units. Chosen over `nohup setsid` because the
user manager is running with `Linger=yes`, so units survive session exit, and
because each unit is a supervised cgroup: worker descendants are tracked
automatically and the exit status is recorded rather than inferred from a PID.
`RemainAfterExit=yes` keeps each unit queryable after it finishes.

Scheduler note: Slurm was not used. A previous probe recorded
`Error creating slurm stream socket: Operation not permitted`, and the configs'
`resources.scheduler: "slurm"` fields are resource *requests*, not a routing
instruction. No scientific config was altered to accommodate this.

`--workers` is deliberately **not** passed on the command line. Each script
resolves `args.workers or cfg.get("workers")`, so the unchanged config file is
the single source of truth for parallelism (B=8, C=16, D=8).

## Jobs

| Exp | systemd unit | ExecMainPID | Config SHA-256 | Output |
|-----|--------------|-------------|----------------|--------|
| A | `proteostasis-phase1-20260731T162946-0500-A.service` | 1269390 | `81d4d05b…6fd692` | `results/phase1/run_20260731T162946-0500/A` |
| B | `proteostasis-phase1-20260731T162946-0500-B.service` | 1269401 | `5802e1a3…c096b76` | `results/phase1/run_20260731T162946-0500/B` |
| C | `proteostasis-phase1-20260731T162946-0500-C.service` | 1269445 | `9c561bd2…c4b8fcc` | `results/phase1/run_20260731T162946-0500/C` |
| D | `proteostasis-phase1-20260731T162946-0500-D.service` | 1269455 | `fcc39896…88d7169` | `results/phase1/run_20260731T162946-0500/D` |

Full config hashes (verified unchanged against
`results/phase1/ops/config_hashes.sha256` immediately before launch):

```
81d4d05b53a65355379616b4b64d616ff52f9feb5b4edce1603c7d14456fd692  configs/phase1/experiment_a.json
5802e1a34c613edd3bd825dae55c3cf8f63b9c772c104b13e37534667c096b76  configs/phase1/experiment_b.json
9c561bd2aee2f9de8f734a21bb00aa1f69eca79db8e5db45c08ac40ecf4b8fcc  configs/phase1/experiment_c.json
fcc398966baf8f3c5a739772c044f464396ad3155afa91dc4460292a188d7169  configs/phase1/experiment_d.json
```

### A — baseline sanity and invariant-domain tests

- Config request: 4 CPUs, 8 GB. Single process, one numerical thread.
- Logs: `results/phase1/run_20260731T162946-0500/ops/experiment_A.{stdout,stderr}.log`

```bash
systemd-run --user --unit=proteostasis-phase1-20260731T162946-0500-A \
  --description='Proteostasis Law Phase 1 experiment A' \
  --property=WorkingDirectory=/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory \
  --property=RemainAfterExit=yes \
  --property=StandardOutput=file:/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/run_20260731T162946-0500/ops/experiment_A.stdout.log \
  --property=StandardError=file:/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/run_20260731T162946-0500/ops/experiment_A.stderr.log \
  --setenv=OMP_NUM_THREADS=1 --setenv=OPENBLAS_NUM_THREADS=1 \
  --setenv=MKL_NUM_THREADS=1 --setenv=NUMEXPR_NUM_THREADS=1 \
  /home/kiran/miniforge3/bin/python /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/scripts/run_experiment_a.py \
    --config /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/configs/phase1/experiment_a.json \
    --outdir /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/run_20260731T162946-0500/A
```

### B — stability maps over (j, nu, chi)

- Config request: 8 CPUs, 8 GB. 8 workers, one numerical thread per worker.
- Logs: `results/phase1/run_20260731T162946-0500/ops/experiment_B.{stdout,stderr}.log`

```bash
systemd-run --user --unit=proteostasis-phase1-20260731T162946-0500-B \
  --description='Proteostasis Law Phase 1 experiment B' \
  --property=WorkingDirectory=/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory \
  --property=RemainAfterExit=yes \
  --property=StandardOutput=file:/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/run_20260731T162946-0500/ops/experiment_B.stdout.log \
  --property=StandardError=file:/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/run_20260731T162946-0500/ops/experiment_B.stderr.log \
  --setenv=OMP_NUM_THREADS=1 --setenv=OPENBLAS_NUM_THREADS=1 \
  --setenv=MKL_NUM_THREADS=1 --setenv=NUMEXPR_NUM_THREADS=1 \
  /home/kiran/miniforge3/bin/python /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/scripts/run_experiment_b.py \
    --config /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/configs/phase1/experiment_b.json \
    --outdir /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/run_20260731T162946-0500/B
```

### C — global parameter sweep (5000 LHS draws, 16 parameters)

- Config request: 16 CPUs, 16 GB. 16 workers, one numerical thread per worker.
- Logs: `results/phase1/run_20260731T162946-0500/ops/experiment_C.{stdout,stderr}.log`

```bash
systemd-run --user --unit=proteostasis-phase1-20260731T162946-0500-C \
  --description='Proteostasis Law Phase 1 experiment C' \
  --property=WorkingDirectory=/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory \
  --property=RemainAfterExit=yes \
  --property=StandardOutput=file:/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/run_20260731T162946-0500/ops/experiment_C.stdout.log \
  --property=StandardError=file:/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/run_20260731T162946-0500/ops/experiment_C.stderr.log \
  --setenv=OMP_NUM_THREADS=1 --setenv=OPENBLAS_NUM_THREADS=1 \
  --setenv=MKL_NUM_THREADS=1 --setenv=NUMEXPR_NUM_THREADS=1 \
  /home/kiran/miniforge3/bin/python /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/scripts/run_experiment_c.py \
    --config /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/configs/phase1/experiment_c.json \
    --outdir /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/run_20260731T162946-0500/C
```

### D — perturbation interaction study

- Config request: 8 CPUs, 8 GB. 8 workers, one numerical thread per worker.
- Logs: `results/phase1/run_20260731T162946-0500/ops/experiment_D.{stdout,stderr}.log`

```bash
systemd-run --user --unit=proteostasis-phase1-20260731T162946-0500-D \
  --description='Proteostasis Law Phase 1 experiment D' \
  --property=WorkingDirectory=/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory \
  --property=RemainAfterExit=yes \
  --property=StandardOutput=file:/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/run_20260731T162946-0500/ops/experiment_D.stdout.log \
  --property=StandardError=file:/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/run_20260731T162946-0500/ops/experiment_D.stderr.log \
  --setenv=OMP_NUM_THREADS=1 --setenv=OPENBLAS_NUM_THREADS=1 \
  --setenv=MKL_NUM_THREADS=1 --setenv=NUMEXPR_NUM_THREADS=1 \
  /home/kiran/miniforge3/bin/python /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/scripts/run_experiment_d.py \
    --config /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/configs/phase1/experiment_d.json \
    --outdir /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/run_20260731T162946-0500/D
```

## Verification at +20 s

Log: `results/phase1/run_20260731T162946-0500/ops/verification_20s.log`

| Exp | ActiveState | SubState | Result | stderr bytes |
|-----|-------------|----------|--------|--------------|
| A | active | running | success | 0 |
| B | active | running | success | 0 |
| C | active | running | success | 0 |
| D | active | running | success | 0 |

An argv-exact `/proc` census returned **36** live canonical processes:
4 parents + 8 (B) + 16 (C) + 8 (D) workers. This count matches the configured
worker totals exactly and confirms **exactly one instance of each experiment**.
Each parent's PGID equals its own PID; every worker's PGID equals its parent's
PID, so the four jobs occupy four disjoint process groups. Per-unit
`TasksCurrent` was A=2, B=13, C=21, D=13.

Progress banners at +20 s (stderr empty for all four):

```
[A] 25 parameter sets x 12 states / A1 done / A2 done / A3-A6 done
[B] 325 (nu,chi) cells x 60 influx values, 8 workers
[C] 5000 LHS draws over 16 parameters, 16 workers
[D] 60 backgrounds x 3 pairs x 6x6 cells, 8 workers
```

## Duplicate cleanup record

Record: `results/phase1/ops/duplicate_cleanup_record_20260731T162755-0500.log`
Census: `results/phase1/ops/host_process_census_20260731T162755-0500.log`

Two earlier duplicate launch attempts recorded PIDs 743-746 and 749-752 (plus an
initial 741-744), all pointed at the shared `results/phase1/{A,B,C,D}`
directories.

**No process was terminated, because none of them was alive on this host.**
Every recorded PID was resolved against `/proc` before any action:

| PID | Recorded as | Host state now |
|-----|-------------|----------------|
| 741, 742, 744, 745, 746, 751, 752 | experiments A-D | no `/proc` entry — already gone |
| 743 | C (attempt 0) / A (attempt 1) | **kernel thread** `kworker/26:1H-kblockd`, PPid 2, state `I (idle)` |
| 749 | A (attempt 2) | **kernel thread** `kworker/40:1H-kblockd`, PPid 2, state `I (idle)` |
| 750 | B (attempt 2) | **kernel thread** `kworker/39:1H-kblockd`, PPid 2, state `I (idle)` |

The duplicates ran under a sandboxed session with its own PID namespace, so
those PIDs never referred to host processes. Three have since been reused by
kernel block-layer workers. Killing by recorded PID would have signalled
unrelated kernel objects, so identity was verified via `/proc/<pid>/status`
(`Name`, `PPid`, zero-length `cmdline`) before acting. No signal was sent.

Independent confirmations that nothing survived:

- `pgrep -af 'run_experiment_[a-d]\.py'` — no matches
- argv-exact `/proc` scan (`results/phase1/ops/scan_live_experiments.sh`) — 0
- `systemctl --user list-units` matching `phase1|experiment|proteostasis` — none
- python processes with PPID 1 — only `firewalld` and a `uvicorn` service, both unrelated

Unrelated long-running jobs were identified and deliberately left untouched:
two `run_md_charmm.py` CHARMM MD jobs (PIDs 815702 and 969170, ~18 h and ~8 h
elapsed).

## Quarantine

`results/phase1/{A,B,C,D}` were **moved**, not deleted, to:

```
results/phase1/quarantine_preclean_20260731T162755-0500/{A,B,C,D}
```

All four were empty — the duplicate attempts never flushed any output — so no
scientific evidence existed to lose. Existing ops logs were left in place in
`results/phase1/ops/`. A `QUARANTINE_README.txt` records the move.

## Commit and push status

- `.git` is writable from this host session. The earlier
  `Unable to create '.git/index.lock': Read-only file system` failure was an
  artifact of the sandboxed session, not a host condition.
- Committed: canonical scripts, tests, configs, `JOBS.md`, `STATUS.md`.
- Not committed: `results/**`. The repository's `.gitignore` ignores `results/`,
  which covers the ops logs and `clean_submission_manifest.json`. There is no
  repository convention tracking manifests, so nothing under `results/` was
  force-added.
- **Not pushed.** The repository has no remote configured at all
  (`git remote -v` is empty, and `master` has no upstream), so no push target
  exists and remote safety could not be established.
