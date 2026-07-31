# Phase 2 jobs — operational record

Operational record only. **No scientific conclusion appears in this file.** The
jobs below were launched and verified alive; nothing here asserts that any
Phase 2 prediction has been tested, confirmed, or falsified.

Read `theory/MATCHED_IMPLEMENTATION_PROTOCOL.md` before interpreting any output
from these jobs — in particular section 6, which states the conditions under
which a percentage may be compared at all.

---

## Phase 2A matched equivalence benchmark — submitted 2026-07-31

Run root (boron): `results/phase2/matched_20260731T175912-0500/`

```
matched_20260731T175912-0500/
  smoke/       rep1, rep2 (workers 1), rep3 (workers 4) -- n = 8 determinism check
  boron/       full matched factorial, n = 2000       <- boron job output
  nitrogen/    reserved for the retrieved free-limit counterpart
  evidence/    test logs, T0 report, smoke and launch verification
  logs/        boron stdout / stderr
```

`results/` is gitignored, so the run root is **not** in the repository. The
submission manifest is `results/phase2/matched_20260731T175912-0500/SUBMISSION_MANIFEST.json`.

### Gate evidence recorded before launch

| gate | outcome |
|---|---|
| Phase 2 targeted suite | `80 passed` (`python -m pytest tests/phase2 -q`), exit 0 |
| full canonical suite | `149 passed, 23 warnings` (`python -m pytest tests -q`), exit 0 |
| T0 epsilon -> 0 identity | **PASSED**, 10/10 checks, exit 0 |
| smoke matched benchmark | **PASSED** — schema, determinism, 0 numerical failures, exact label agreement at epsilon = 1e-6 |

The 23 warnings are the pre-existing Python multiprocessing/fork deprecation
warning, unchanged from Phase 1.

---

### Job 1 — boron, full matched factorial

- **Host:** boron
- **Launcher:** `systemd-run --user` transient unit (survives the shell)
- **Unit:** `phase2-matched-boron-20260731T180846.service`
- **MainPID:** `1446897`
- **Started:** Fri 2026-07-31 18:08:46 CDT
- **Resources:** 16 worker processes, **one numerical-library thread per worker**
  (`OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=MKL_NUM_THREADS=NUMEXPR_NUM_THREADS=VECLIB_MAXIMUM_THREADS=1`,
  verified by reading `/proc/1446897/environ` of the live process)
- **Cells:** `{boron, free}` x 7 epsilon rungs x `{nitrogen, boron}` protocols = **28 cells**
- **Output:** `results/phase2/matched_20260731T175912-0500/boron/`
- **Logs:** `results/phase2/matched_20260731T175912-0500/logs/boron_matched.{stdout,stderr}.log`

```bash
systemd-run --user --unit=phase2-matched-boron-20260731T180846 \
  --description="Phase 2A matched equivalence benchmark (boron, n=2000, seed 20260731)" \
  --property=WorkingDirectory=/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory \
  --property=StandardOutput=append:/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase2/matched_20260731T175912-0500/logs/boron_matched.stdout.log \
  --property=StandardError=append:/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase2/matched_20260731T175912-0500/logs/boron_matched.stderr.log \
  --setenv=OMP_NUM_THREADS=1 --setenv=OPENBLAS_NUM_THREADS=1 \
  --setenv=MKL_NUM_THREADS=1 --setenv=NUMEXPR_NUM_THREADS=1 \
  --setenv=VECLIB_MAXIMUM_THREADS=1 \
  /home/kiran/miniforge3/bin/python \
  /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/scripts/phase2/run_matched_benchmark.py \
    --outdir /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase2/matched_20260731T175912-0500/boron \
    --n-samples 2000 --seed 20260731 --workers 16 \
    --model-forms boron,free --protocols nitrogen,boron \
    --epsilons 1e-6,1e-3,1e-2,1e-1,0.3,1.0,2.0 \
    --label matched_boron_full
```

**Verification at +32 s** (`evidence/boron_verification.txt`):
`ActiveState=active`, `SubState=running`, `Result=success`, parent state `Ssl`,
exactly 16 worker children, thread environment confirmed inside the live
process, stderr **0 bytes**, exactly one benchmark instance on the host
(17 processes = 1 parent + 16 workers; no other `run_matched_benchmark.py`).

---

### Job 2 — nitrogen, free-limit counterpart

Nitrogen runs the **free arm only**. The free-limit model form is the
`epsilon = 0` face of the boron family (protocol document section 2), and the
nitrogen host must never import the boron `proteostasis` package — asserted
structurally by
`tests/phase2/test_benchmark_design.py::TestFreeArmNeedsNoBoronCode`.

- **Host:** nitrogen (10.147.17.116)
- **Scheduler:** Slurm, partition `main`
- **Job ID:** `4`
- **State at +40 s:** `RUNNING` (`JobState=RUNNING`, `Reason=None`)
- **Resources:** `NumCPUs=1`, `NumTasks=1`, `CPUs/Task=1`, `TRES=cpu=1,mem=4G,node=1`
- **Live PID on nitrogen:** `1566011`
- **Workspace (new, unique):**
  `/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory-nitrogen-check/matched_20260731T180944/`
- **WorkDir:** as above
- **Command:** `.../matched_20260731T180944/matched_free_limit.sbatch`
- **StdOut:** `.../matched_20260731T180944/logs/matched_free.4.out`
- **StdErr:** `.../matched_20260731T180944/logs/matched_free.4.err` (0 bytes at verification)
- **Cells:** `{free}` x 7 epsilon rungs x `{nitrogen, boron}` protocols = **14 cells**

`--workers 1` is deliberate. The script's parallelism *is* deterministic — the
smoke check confirmed `payload_sha256` is identical at `workers = 1` and
`workers = 4` — but a single core removes the question entirely and matches the
1 CPU / 4 GB allocation.

Payload executed under Slurm:

```bash
python3 run_bench.py \
  --outdir .../matched_20260731T180944/out \
  --n-samples 2000 --seed 20260731 --workers 1 \
  --model-forms free --protocols nitrogen,boron \
  --epsilons 1e-6,1e-3,1e-2,1e-1,0.3,1.0,2.0 \
  --label matched_free_limit_nitrogen
```

`run_bench.py` is a four-line shim that puts the staged directory on `sys.path`
and delegates to `phase2/run_matched_benchmark.py` via `runpy`. The benchmark
source itself was **not** modified for nitrogen.

#### What was copied, and proof it is the same code

Only the modules the free arm needs were staged. `boron_continuation.py` and the
`proteostasis` package were deliberately **not** copied.

SHA-256, boron side and nitrogen side, byte-identical:

| file | sha256 |
|---|---|
| `__init__.py` | `bf3b4bda585753cbc839e80e43af223bebe19526edc778e2f617e0264d53fcf2` |
| `lhs.py` | `f9dc9f2456118074d41ba6663b76e9ed9bc95dc6512a615922b490734c167e58` |
| `mapping.py` | `a413536da5ae6a599253ff633deed38f95dc9957d66d972110becd087fe65621` |
| `models.py` | `61322d9f25ab79f15a2399a41e64d7837567d2d7a87e3311fbb523e7ae3e3153` |
| `nitrogen_limit.py` | `6c6cdb25b11ea03edb36ef792a8a96df870276c78b9d0866f1393bb56d2074c7` |
| `protocols.py` | `27f637fd4c824b6c8cc7d25efa36b8d2ed9f551dc64ac7ec942043b234e9f448` |
| `run_matched_benchmark.py` | `5631e52db5f27acc26f12abc7416bb33568bebcefeb1f325809897103c1e1092` |

Job scripts: `matched_free_limit.sbatch`
`e4745a6666622f72c90a4da8dabe00a40ac11746e81c2d6cc36c6e944f9de1e2`,
`run_bench.py`
`4382602b05780786630b07bc14cfe3f93169b3bac91ab5be0fe8b2b94e65708f`.

#### Proof both hosts score the same parameter draws

This is the condition without which no cross-host cell may be compared at all.
The n = 2000, seed 20260731 sample matrix was exported from boron and staged
alongside the source. On nitrogen, under **numpy 1.26.4** (boron runs a different
numpy), the matrix was regenerated and checked:

```
regenerated sha256  = a937f4bc68faa9bb404bfa81c190a168f19a866d5402954a691e9b2d650d85fa
matches pinned      = True     (phase2/lhs.py::SAMPLE_HASHES[2000])
bit-identical to staged matrix from boron = True
```

---

## (!) Not the matched result

The pre-existing **50,000-sample independent nitrogen sweep is not the matched
counterpart** and must not be reported as the matched result. It uses a
different sample matrix, so it cannot be differenced against any boron cell —
the two would be scoring different parameter draws. It remains valid only as a
free-limit result in its own right. See protocol document section 6.

## (!) Percentages

**Percentages are not comparable outside the matched benchmark.** The five
conditions are listed in protocol document section 6. The benchmark's own
manifest repeats the constraint:

> no percentage in this output may be interpreted until every matched cell has
> completed.

Both jobs above were **still running** when this record was written. Nothing in
their partial output has been interpreted.

## Concurrent work on this host — do not disturb

- Phase 1 experiment D (`run_experiment_d.py`, PID `1269455` + 8 workers) is
  still running into `results/phase1/run_20260731T162946-0500/D`. Untouched.
- A separate Phase 2 audit wrote into
  `results/phase2/audit_20260731T172711-0500/`. Untouched, and gitignored, so it
  cannot enter a commit. Its `c02` stage had already recorded `done in 1132.4s`
  before this session's work began; no process for the reported parent PID
  `1361326` was present on boron at that time (see `HISTORY.md`).
