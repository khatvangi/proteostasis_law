# Phase 1 experiment D — checkpointed recovery run

Operational record for the checkpointed re-run of experiment D. Written at launch;
the numbers below are launch-time facts, not results. No scientific conclusion is
drawn here — the run was still in progress when this file was committed.

## Why this run exists

`scripts/run_experiment_d.py` evaluates all 60 backgrounds through a single
`Pool.map`. `Pool.map` returns only when every task has returned, and the script
writes nothing before that point. A small number of stiff parameter draws
therefore holds the entire experiment hostage and leaves **zero** partial output.

Observed on the live run at 22:11 on 2026-07-31: 5 h 40 min wall, 9 h 26 min CPU
across 9 processes (one parent, eight workers), and an output directory that was
still completely empty. Measured per-background cost in this recovery: the median
background finishes in ~25 s, and only backgrounds **14, 19, 37 and 44** exceed
45 s. Four pathological draws out of sixty were blocking the other fifty-six.

## The original run is untouched and still live

| | |
|---|---|
| Unit | `proteostasis-phase1-20260731T162946-0500-D.service` |
| MainPID | `1269455` |
| State at recovery launch | `active (running)`, `NRestarts=0` |
| Started | Fri 2026-07-31 16:30:01 CDT |
| Output directory | `results/phase1/run_20260731T162946-0500/D/` — still empty, mtime 16:29 |

No signal was sent to that unit or to any of its PIDs. It was not stopped,
reniced, traced, or attached to. Its output directory was not created, written,
moved, or read from. The recovery writes to a different root entirely.

The recovery's own children are spawned with `start_new_session=True`, so each
leads its own process group, and the timeout path signals a group only after
asserting `pgid == child pid`. It is structurally incapable of reaching the
original run's process group, the controller's own group, or any unrelated job on
this host. The concurrent Phase 2 closure audit under
`results/phase2/closure_20260731T220024-0500` was likewise not touched.

## Commit

| | |
|---|---|
| Commit | `ee64a3ff9b024b54c523294674cbdc401f6c5d9b` (`ee64a3f`) |
| Branch | `master` |
| Message | Add checkpointed, timeout-bounded recovery runner for experiment D |
| Contents | `scripts/run_experiment_d_checkpointed.py`, `tests/test_experiment_d_checkpointed.py` |

No results were committed — `results/` is covered by `.gitignore`. No config,
theory, manuscript or existing script was modified. `configs/phase1/experiment_d.json`
is byte-unchanged and is the same file the original run is using.

## Unit and command

| | |
|---|---|
| Unit | `proteostasis-phase1-D-checkpointed-20260731T223225-0500.service` |
| MainPID (controller) | `1622784` |
| Started | Fri 2026-07-31 22:32:25 CDT |
| Host | boron |
| Output root | `results/phase1/D_checkpointed_20260731T223225-0500/` (gitignored via `results/`) |

```
systemd-run --user \
  --unit=proteostasis-phase1-D-checkpointed-20260731T223225-0500 \
  --description="Proteostasis Law Phase 1 experiment D checkpointed recovery (60 backgrounds, 16 concurrent, 3600 s/background)" \
  -p RemainAfterExit=yes \
  -p WorkingDirectory=/storage/kiran-stuff/proteostasis_law/proteostasis-law-theory \
  -p Environment="OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1" \
  -p StandardOutput=file:.../controller.out \
  -p StandardError=file:.../controller.err \
  /home/kiran/miniforge3/bin/python \
    /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/scripts/run_experiment_d_checkpointed.py \
    --config /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/configs/phase1/experiment_d.json \
    --outdir /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory/results/phase1/D_checkpointed_20260731T223225-0500 \
    --concurrency 16 --background-timeout 3600 --threads 1
```

16 concurrent background subprocesses, one BLAS thread each (60 backgrounds x
multi-threaded BLAS would oversubscribe the host, which also runs the original D
and the Phase 2 matched benchmark).

## Identity hashes

Recorded in `run_manifest.json` and in every checkpoint. A checkpoint is reused on
resume only if its background index, config hash, source hash **and** sample hash
all match, and its recorded payload hashes verify.

| | |
|---|---|
| `config_hash` | `18280a039eea9d1cac54a18f09e865b56da69feadd510b80ef15cb0555ab29a5` |
| `source_hash` | `809b3cd7a312e1bfd1e63f5760eccce021cbcf8010a185548c2619daf3b205a9` |
| `sample_matrix_hash` | `73b6a58d07462256e08c0e77cdd30fe5e59f450792ffac2ca4d96362fd1b7975` |

`source_hash` is over the files that determine the numbers:

```
scripts/run_experiment_d.py          9ff94b8f6e676587e6153087cb247ac3174918d556feaa077ca8c7347f938c67
scripts/proteostasis/__init__.py     2ddbe26173697aec44b7188ddc7588713a17401c67ce635e462596e3fc8b074a
scripts/proteostasis/model.py        bf3b642f147315f51454ce6304d6e6152c51b483a8eaee5d971b2f75c650ff0c
scripts/proteostasis/simulate.py     d4d0cdeed0b5d4adc6589940dbc80978cf338d37d74a16bcbadf2160bd5f4883
scripts/proteostasis/equilibria.py   c00bb65eee6ac126c09a2e4d38354819a5876f447ff96a77c8c731ec11c3ebd3
scripts/proteostasis/sweeps.py       d3cd0649e39905d20a49af2ad1ec0544c8861e3635b861ef9acf1d43d72631d3
```

The recovery runner itself is deliberately **excluded** from `source_hash`: it
schedules work and does not compute anything, so an orchestration edit must not
invalidate finished science. A test enforces the exclusion.

Environment: python 3.12.11, numpy 2.2.6, scipy 1.14.0, pandas 2.3.3.

## What is preserved

The science is not reimplemented. The runner imports `run_experiment_d` and calls
that module's own `_backgroundTask`, `_pairSummary` and `PAIRS`. A test greps the
runner's source and fails if `_score`, `_perturb`, `_backgroundTask` or
`_pairSummary` is ever redefined there.

Preserved verbatim: the equations; the LHS seed (20260731) and the 60-background
sample matrix; the baseline fraction 0.4; the three perturbation pairs and their
axes; burden levels `[1.0, 1.2, 1.4, 1.7, 2.0, 2.5]`; capacity levels
`[1.0, 0.9, 0.8, 0.7, 0.6, 0.5]`; `burden_threshold` 1.0; `burden_censor` 1000.0;
`t_end` 50000.0; `n_out` 200; `rtol` 1e-8; `atol` 1e-11; the additive,
multiplicative and Bliss nulls; synthetic-collapse recording; and the summary
logic including the restriction to genuinely double perturbations.

Changed: only scheduling, checkpointing and time bounds.

Checkpoint rows are stored as JSON, not TSV, because `writeTable`'s `%.12g` would
truncate float64. Formatting happens once, at merge, so the merged
`interactions.tsv` is byte-identical to a direct write of the same rows.

## Timeout semantics

**A timed-out background is `unresolved_timeout`. It is not a numerical failure.**

A background exceeding the 3600 s wall limit is recorded in
`unresolved/background_XXXX.json` with status `unresolved_timeout`. It:

- contributes **no rows** to `interactions.tsv`;
- appears in **no** interaction summary, per-pair or overall;
- is **not** counted in `n_errors`, which counts only genuine numerical failures
  raised inside `_backgroundTask` on backgrounds that ran to completion;
- is **not** counted in `unusable_reasons`, which records model-level exclusions
  such as "no viable state at j_lo";
- is reported separately under `summary.json → recovery.n_unresolved_timeout` and
  `recovery.unresolved_timeout_backgrounds`.

The honest reading of an unresolved background is *"this draw was not evaluated
within the budget"* — a statement about the budget, not about the model. Any
downstream use of this run must report the unresolved count alongside the
interaction statistics and must not describe those backgrounds as failures.

A process that exits non-zero without producing a valid checkpoint is a distinct
category, `failed_process`, also reported separately.

## Output layout

```
results/phase1/D_checkpointed_20260731T223225-0500/
  run_manifest.json                    identity hashes, concurrency, timeout
  progress.json                        rewritten atomically on every state change
  controller.out / controller.err
  logs/background_XXXX.{out,err}
  backgrounds/background_XXXX/
      rows.json       factorial rows, exact float repr
      meta.json       background metadata incl. every p_<param>
      checkpoint.json index, hashes, the sample itself, runtime
      DONE            written last; contains the sha256 of checkpoint.json
  unresolved/background_XXXX.json      unresolved_timeout / failed_process
  interactions.tsv / backgrounds.tsv / summary.json / provenance.json   (at merge)
```

Resume is safe at any point: `--merge-only` re-merges without computing, and
rerunning the controller recomputes only backgrounds without a valid checkpoint.

## Tests

`tests/test_experiment_d_checkpointed.py` — 38 tests, all passing.

- **Sample matrix.** The regenerated LHS matrix equals the original construction
  element for element, under exact float equality, for all 60 backgrounds.
- **Numerical equivalence.** For four measured-fast backgrounds — 48 and 56 (full
  108-row factorial) and 5 and 35 (early-return path) — the checkpointed output
  is compared against a direct `run_experiment_d._backgroundTask` call via
  canonical JSON, which is bitwise equality of every float. A companion test
  asserts the two usable backgrounds really do produce all 108 cells, so the
  equivalence check cannot pass vacuously.
- **Merged output.** The merged `interactions.tsv` is compared **byte for byte**
  against a direct `writeTable` of the same rows, and the merged summary is
  compared against the original `run_experiment_d.main` summary logic recomputed
  independently from the direct results.
- **Resume.** After a second controller pass, every checkpoint file is unchanged
  in both mtime and content — proving the work was skipped, not recomputed.
- **Rejection.** Eight modes are each rejected with a specific reason: tampered
  rows, truncated rows, tampered meta, unparseable checkpoint, config-hash
  mismatch, source-hash mismatch, sample-hash mismatch, background-index
  mismatch, missing completion marker, stale completion marker, missing payload.
  A control test confirms an uncorrupted copy is accepted, and a further test
  confirms a rejected checkpoint is recomputed rather than trusted.
- **Timeout.** A real controller run with a 0.5 s limit exercises the real kill
  path: all backgrounds are classified `unresolved_timeout`, `n_errors` is 0,
  `n_cells` is 0, `overall` and `by_pair` are empty, no `DONE` is left behind,
  and the unresolved markers state that this is not a numerical failure.
- **Selection.** `--only` restricts scheduling without perturbing the sample
  matrix; an out-of-range index is rejected.

Suite results:

```
python -m pytest tests/test_experiment_d_checkpointed.py -q   ->  38 passed
python -m pytest tests -q                                     -> 187 passed
```

187 = 149 pre-existing + 38 new. No pre-existing test was modified and none
regressed.

## Smoke run

Backgrounds 5, 48 and 56, concurrency 3, in a temporary root, before launch.

| Check | Result |
|---|---|
| Fresh run | 3/3 completed in 18.1 s, 0 unresolved, 0 failed |
| Rerun into same root | 3/3 resumed in 0.1 s; all 12 checkpoint files identical in mtime and content |
| Independent fresh run | `interactions.tsv` and `backgrounds.tsv` byte-identical to the first run |
| Per-background payloads | `rows.json` and `meta.json` sha256 equal across the two independent runs for all three backgrounds |
| Summary | scientific fields identical; `recovery` block identical except `wall_s` and `n_resumed_from_checkpoint` |
| stderr | empty everywhere |

`interactions.tsv` sha256 `0fa5c6f5731f51fba09a730751d4c6d1a09eb94af4c51cb7717cbc8cedaac13a`
(216 cells, 2 usable backgrounds), `backgrounds.tsv` sha256
`904afb82a660c1c5f809202be4aa825cc4783963bb2d335a04de7b8eec715563`, both
reproduced identically by the independent run.

## Launch verification (+40 s)

| Check | Result |
|---|---|
| Recovery unit | `active (running)`, `NRestarts=0`, controller uptime 40 s |
| Exactly one recovery unit | 1 unit matching `D-checkpointed` |
| Exactly one process tree | 17 processes — 1 controller + 16 background workers — all 17 inside the unit's cgroup |
| Checkpoints appearing | 19 backgrounds complete with `DONE`; 16 running, 25 pending |
| Unresolved / failed so far | 0 / 0 |
| `progress.json` | present, `concurrency=16`, `background_timeout_s=3600.0`, `n_backgrounds=60` |
| stderr | `controller.err` 0 bytes; 0 non-empty per-background stderr files |
| Original D | still `active (running)`, MainPID `1269455`, `NRestarts=0`, output dir still empty |

Nineteen backgrounds were checkpointed in the first forty seconds. The original
run has produced no output in six hours.

## Reading the result when it finishes

`summary.json` keeps the original field names with their original meanings over
the backgrounds that resolved, and puts all recovery bookkeeping under a separate
`recovery` key. Note that `n_backgrounds` is the number of rows in
`backgrounds.tsv`, i.e. resolved backgrounds — `recovery.n_backgrounds_requested`
is 60. Report `recovery.n_unresolved_timeout` alongside any interaction statistic
taken from this run.
