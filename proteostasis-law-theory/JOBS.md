# Phase 1 jobs

## Gate and provenance

- Full-suite command: `/home/kiran/miniforge3/bin/python -m pytest`
- Result: `69 passed, 19 warnings in 60.62s (0:01:00)`, exit code 0.
- Coverage retained: dimensions, positivity/invariant domain, conserved-pool
  and total-resource balances, analytic Jacobian correctness, fold finding,
  reachability, reproducibility, and recovery dynamics.
- Launch method: separate `nohup setsid` processes held by a supervisory shell.
- Scheduler probe: `timeout 3s sinfo` and `timeout 3s squeue`; both exited 124
  after `Error creating slurm stream socket: Operation not permitted`.
- Absolute Python: `/home/kiran/miniforge3/bin/python` (Python 3.12.11).
- Repository HEAD visible at launch:
  `73e1c0ab341530f7b2d369c459008674897f6287` (working tree contains the
  canonical repair because `.git` is read-only and cannot accept a commit).
- Numerical-library settings for every job: `OMP_NUM_THREADS=1`,
  `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `NUMEXPR_NUM_THREADS=1`.
- 15-second verification:
  `results/phase1/ops/verification_15s_20260731T162300-0500.log`.

## A

- PID: 743; verified `Rsl` and alive at 15 seconds.
- Resources: config request 4 CPUs, 8 GB; one top-level process, numerical
  libraries restricted to one thread.
- Config: `configs/phase1/experiment_a.json`.
- SHA-256: `81d4d05b53a65355379616b4b64d616ff52f9feb5b4edce1603c7d14456fd692`.
- Output: `results/phase1/A`.
- Logs: `results/phase1/ops/experiment_A_20260731T162300-0500.{stdout,stderr}.log`.
- PID file: `results/phase1/ops/experiment_A_20260731T162300-0500.pid`.

```bash
nohup setsid env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 /home/kiran/miniforge3/bin/python /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/scripts/run_experiment_a.py --config /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/configs/phase1/experiment_a.json --outdir /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/A
```

## B

- PID: 744; verified `Ssl` and alive at 15 seconds.
- Resources: 8 workers/CPUs, 8 GB; one numerical-library thread per worker.
- Config: `configs/phase1/experiment_b.json`.
- SHA-256: `5802e1a34c613edd3bd825dae55c3cf8f63b9c772c104b13e37534667c096b76`.
- Output: `results/phase1/B`.
- Logs: `results/phase1/ops/experiment_B_20260731T162300-0500.{stdout,stderr}.log`.
- PID file: `results/phase1/ops/experiment_B_20260731T162300-0500.pid`.

```bash
nohup setsid env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 /home/kiran/miniforge3/bin/python /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/scripts/run_experiment_b.py --config /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/configs/phase1/experiment_b.json --outdir /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/B --workers 8
```

## C

- PID: 745; verified `Ssl` and alive at 15 seconds.
- Resources: 16 workers/CPUs, 16 GB; one numerical-library thread per worker.
- Config: `configs/phase1/experiment_c.json`.
- SHA-256: `9c561bd2aee2f9de8f734a21bb00aa1f69eca79db8e5db45c08ac40ecf4b8fcc`.
- Output: `results/phase1/C`.
- Logs: `results/phase1/ops/experiment_C_20260731T162300-0500.{stdout,stderr}.log`.
- PID file: `results/phase1/ops/experiment_C_20260731T162300-0500.pid`.

```bash
nohup setsid env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 /home/kiran/miniforge3/bin/python /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/scripts/run_experiment_c.py --config /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/configs/phase1/experiment_c.json --outdir /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/C --workers 16
```

## D

- PID: 746; verified `Ssl` and alive at 15 seconds.
- Resources: 8 workers/CPUs, 8 GB; one numerical-library thread per worker.
- Config: `configs/phase1/experiment_d.json`.
- SHA-256: `fcc398966baf8f3c5a739772c044f464396ad3155afa91dc4460292a188d7169`.
- Output: `results/phase1/D`.
- Logs: `results/phase1/ops/experiment_D_20260731T162300-0500.{stdout,stderr}.log`.
- PID file: `results/phase1/ops/experiment_D_20260731T162300-0500.pid`.

```bash
nohup setsid env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 /home/kiran/miniforge3/bin/python /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/scripts/run_experiment_d.py --config /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/configs/phase1/experiment_d.json --outdir /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/D --workers 8
```

## Commit status

Attempted command:

```bash
git add configs scripts tests && git commit -m "Repair saddle reachability and recovery timing tests"
```

It failed before staging or committing with
`fatal: Unable to create '.git/index.lock': Read-only file system`.  No push
was attempted.  A writable Git metadata mount is required to produce the
requested commit hash.
