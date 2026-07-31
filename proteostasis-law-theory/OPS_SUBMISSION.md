# Phase 1 operational submission

- Submission timestamp: `2026-07-31T16:12:20-05:00`
- Verification timestamp: `2026-07-31T16:12:36-05:00` (15 seconds after launch)
- Host: `boron` (`hostname -f` returned `localhost`)
- Canonical repository HEAD: `73e1c0ab341530f7b2d369c459008674897f6287`
- Repository: `/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory`
- Submission mode: detached local CPU jobs using `nohup setsid`; Slurm controller unavailable
- Python: `/home/kiran/miniforge3/bin/python` (`Python 3.12.11`)
- NumPy: `2.2.6`
- SciPy: `1.14.0`
- Host logical CPUs: `64`

## Duplicate check

Before submission, no `OPS_SUBMISSION.md`, phase-1 operational manifest, matching live `run_experiment_[a-d].py` process, or pre-existing phase-1 job record was found. No existing job was duplicated.

An initial detached launch at approximately `2026-07-31T16:11:55-05:00` was immediately reaped by the execution sandbox: none of the four PIDs survived and all stdout/stderr files remained empty. It produced no experiment output. The jobs below are the single live submissions, relaunched in a persistent supervisory shell so the detached processes remain available in the execution namespace.

## Fast test suite

Exact command:

```bash
/home/kiran/miniforge3/bin/python -m pytest -q
```

Result: exit code `1`; `67 passed, 2 failed, 19 warnings in 69.81s (0:01:09)`.

Failures:

- `tests/test_equilibria.py::TestReachability::testAboveTheUnstableStateTheSystemEscapes`: expected `blowup`, observed `timeout`.
- `tests/test_simulate.py::TestRecoveryTime::testRecoveryTimeGrowsTowardTheFold`: recovery times were `[25.012506253126563, 25.012506253126563, 25.012506253126563]`, not strictly increasing.

No scientific code was changed. Full output is in `results/phase1/ops/fast_tests.log`; the numeric exit status is in `results/phase1/ops/fast_tests.exitcode`.

## Scheduler detection

`sbatch`, `squeue`, and `sinfo` are installed under `/usr/bin`. Read-only `sinfo` and `squeue` probes were each bounded by `timeout 8s`; both repeatedly returned `Error creating slurm stream socket: Operation not permitted` and timed out with exit code `124`. Scheduler submission was therefore unavailable. Full probe output is in `results/phase1/ops/scheduler_probe.log`.

## Jobs

All commands set numerical-library thread counts to one per worker: `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1`.

### Experiment A

- PID: `741`
- Config: `configs/phase1/experiment_a.json`
- Config SHA-256: `81d4d05b53a65355379616b4b64d616ff52f9feb5b4edce1603c7d14456fd692`
- Resources: config requests 4 CPUs and 8 GB; local launcher is a single process with numerical-library threads limited to 1.
- Output directory: `results/phase1/A`
- Verification: alive after 15 seconds, process state `Rsl`, elapsed `00:15`; stderr empty.

```bash
nohup setsid env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 /home/kiran/miniforge3/bin/python /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/scripts/run_experiment_a.py --config /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/configs/phase1/experiment_a.json --outdir /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/A > /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/ops/experiment_A.stdout.log 2> /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/ops/experiment_A.stderr.log < /dev/null &
```

PID file: `results/phase1/ops/experiment_A.pid`.

### Experiment B

- PID: `742`
- Config: `configs/phase1/experiment_b.json`
- Config SHA-256: `5802e1a34c613edd3bd825dae55c3cf8f63b9c772c104b13e37534667c096b76`
- Resources: 8 workers/CPUs and 8 GB requested; launched with `--workers 8` and one numerical-library thread per worker.
- Output directory: `results/phase1/B`
- Verification: alive after 15 seconds, process state `Ssl`, elapsed `00:15`; stderr empty.

```bash
nohup setsid env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 /home/kiran/miniforge3/bin/python /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/scripts/run_experiment_b.py --config /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/configs/phase1/experiment_b.json --outdir /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/B --workers 8 > /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/ops/experiment_B.stdout.log 2> /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/ops/experiment_B.stderr.log < /dev/null &
```

PID file: `results/phase1/ops/experiment_B.pid`.

### Experiment C

- PID: `743`
- Config: `configs/phase1/experiment_c.json`
- Config SHA-256: `9c561bd2aee2f9de8f734a21bb00aa1f69eca79db8e5db45c08ac40ecf4b8fcc`
- Resources: 16 workers/CPUs and 16 GB requested; launched with `--workers 16` and one numerical-library thread per worker.
- Output directory: `results/phase1/C`
- Verification: alive after 15 seconds, process state `Ssl`, elapsed `00:15`; stderr empty.

```bash
nohup setsid env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 /home/kiran/miniforge3/bin/python /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/scripts/run_experiment_c.py --config /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/configs/phase1/experiment_c.json --outdir /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/C --workers 16 > /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/ops/experiment_C.stdout.log 2> /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/ops/experiment_C.stderr.log < /dev/null &
```

PID file: `results/phase1/ops/experiment_C.pid`.

### Experiment D

- PID: `744`
- Config: `configs/phase1/experiment_d.json`
- Config SHA-256: `fcc398966baf8f3c5a739772c044f464396ad3155afa91dc4460292a188d7169`
- Resources: 8 workers/CPUs and 8 GB requested; launched with `--workers 8` and one numerical-library thread per worker.
- Output directory: `results/phase1/D`
- Verification: alive after 15 seconds, process state `Ssl`, elapsed `00:16`; stderr empty.

```bash
nohup setsid env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 /home/kiran/miniforge3/bin/python /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/scripts/run_experiment_d.py --config /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/configs/phase1/experiment_d.json --outdir /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/D --workers 8 > /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/ops/experiment_D.stdout.log 2> /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/ops/experiment_D.stderr.log < /dev/null &
```

PID file: `results/phase1/ops/experiment_D.pid`.

The consolidated liveness evidence is in `results/phase1/ops/verification_15s.log`. Package/host details and config hashes are in `environment.log` and `config_hashes.sha256` in the same directory.

## Git status after operational recording

No commit was created. `git status --short` reported the following. `OPS_SUBMISSION.md` is the operational file created by this submission; the scientific/config/test paths were pre-existing and were not modified. Operational result paths are ignored by the repository:

```text
?? OPS_SUBMISSION.md
?? configs/
?? scripts/proteostasis/
?? scripts/run_experiment_a.py
?? scripts/run_experiment_b.py
?? scripts/run_experiment_c.py
?? scripts/run_experiment_d.py
?? tests/_context.py
?? tests/test_equilibria.py
?? tests/test_model.py
?? tests/test_reproducibility.py
?? tests/test_simulate.py
?? tests/test_units.py
```
