# History

## 2026-07-31

**Q:** Host-level operational cleanup and a single clean Phase 1 submission on boron:
find and gracefully kill any live `run_experiment_[a-d].py` processes from the two
duplicate launch attempts (recorded PIDs 743-746 and 749-752), quarantine the existing
`results/phase1/{A,B,C,D}`, confirm the repaired test suite passed, launch exactly one
clean run of A-D, verify at 20 s, write a manifest plus JOBS.md/STATUS.md, and commit.
Do not touch scientific code, tests, configs, theory, or manuscript.

**A:** No live experiment processes existed on the host — nothing was killed. The
recorded PIDs came from a sandboxed session with its own PID namespace; on boron,
743/749/750 have since been reused by kernel worker threads (`kworker/*-kblockd`,
PPid 2), so killing by recorded PID would have signalled unrelated kernel objects.
Identity was verified via `/proc/<pid>/status` before any action. The empty A-D
directories were moved (not deleted) to
`results/phase1/quarantine_preclean_20260731T162755-0500/`.

The test gate passed (exit code 0, `69 passed`), so one instance each of A-D was
launched as `systemd-run --user` transient units into
`results/phase1/run_20260731T162946-0500/`, with thread limits verified inside each
running process. At +20 s all four were active with empty stderr and exactly 36 live
processes (4 parents + 32 workers), matching the configs. Experiment A has since
finished successfully (11/11 internal-consistency checks); B, C, D still running.

Committed as `a17dfaf` on `master`. Not pushed — the repository has no remote
configured at all. Results and ops logs are excluded by the `results/` gitignore rule,
so the manifest was written but not committed.

**Note for future sessions:** a second agent session was concurrently active in this
repo during this work. Re-check for duplicate launches before starting any new run.

**Q:** Complete the unfinished operational submission for Phase 2A: run the Phase 2
and full pytest suites and fix only real defects, run a smoke matched benchmark,
write `theory/MATCHED_IMPLEMENTATION_PROTOCOL.md` / `PHASE2_JOBS.md` / STATUS
updates, create a timestamped matched run root, launch the full boron n=2000
matched benchmark as a persistent systemd unit, submit the free-limit counterpart
on nitrogen via Slurm, write a submission manifest, and commit. Do not touch the
c02 audit or Phase 1 D.

**A:** Suites passed with no pre-existing failures: `tests/phase2` 76 passed,
full `tests` 145 passed. T0 was rerun and reproduced its headline exactly
(RHS 4.800103e-06, Jacobian 8.263311e-06 at epsilon 1e-6; slopes 1.0023 / 1.0028),
10/10 checks, exit 0.

The smoke run exposed one genuine defect. Every cell summary recorded only
`tsv_sha256`, a byte hash of a file whose rows each carry a wall-clock `seconds`
column, so that hash can never match between boron's free arm and nitrogen's free
arm even when every computed value is identical — which is exactly the cross-host
comparison the benchmark exists to make. Diagnosed by column-wise diff across
three replicates: `seconds` was the *only* differing column, including between
`workers=1` and `workers=4`. Fix was additive — a `payload_sha256` over the result
columns plus `payload_excludes` — with four new tests, one of which asserts that
`seconds` is the only non-deterministic column so the exclusion set cannot drift
silently. No scientific criterion, T0 threshold, or mapping was touched. Suites
after the fix: 80 and 149 passed.

Smoke then passed on all criteria: schema equal to `ROW_FIELDS`, `payload_sha256`
identical across a repeat run and across worker counts, zero numerical failures,
and 8/8 exact label *and* admissibility agreement between the boron and free arms
at epsilon 1e-6 under both protocols.

`theory/MATCHED_IMPLEMENTATION_PROTOCOL.md` did not exist even though `mapping.py`
and `t0_equivalence.py` cite its sections 3 and 5 by number — those were dangling
references to the derivation licensing the tolerances. It was written with the
term-by-term flux identity verified against `model.py`, the a priori tolerance
derivations, and explicit statements that nitrogen is the epsilon = 0 limit and
that percentages are not comparable outside the matched benchmark.

Boron runs the full 28-cell factorial as unit
`phase2-matched-boron-20260731T180846.service` (MainPID 1446897, 16 workers, one
thread each, stderr empty, exactly one instance). Nitrogen runs the 14-cell
free-limit counterpart as Slurm job 4, RUNNING at 1 CPU / 4 GB. Both hosts were
proven to regenerate a bit-identical n=2000 seed-20260731 sample matrix under
different numpy builds (2.2.6 vs 1.26.4) — without that, no cross-host cell could
be compared. Run root `results/phase2/matched_20260731T175912-0500/`, manifest
`SUBMISSION_MANIFEST.json`. No result was interpreted; both jobs were still
running at commit time.

**Note:** the reported audit parent PID 1361326 with 24 workers had no matching
process on boron; its `c02` stage had already written `done in 1132.4s` before this
session began. Its directory was left untouched and is gitignored. Phase 1 D
(PID 1269455) was still running and was not signalled.

**Q:** Operational continuation of the "already-running" Phase 2B multistability audit:
verify the live c02 process (parent PID 1361326, `--workers 24 --n-controls 200`, started
17:35:30), create a launcher + manifest under the gitignored audit root only, run it as a
persistent `systemd-run --user` unit that waits for c02 to exit, gates on completeness,
then runs c03/c04/c05 and d01. No tracked edits, no commits, no interference with the
matched-benchmark task.

**A:** The stated premise was wrong on two points, both verified before acting. PID 1361326
does not exist on boron (`/proc/1361326` absent), and a full `/proc` sweep matching on the
string `c02_dense_multistability` in the command line found no live process at all — the
same PID-namespace mismatch recorded earlier today. More importantly, **c02 had already
finished**: its log reads `304 draws (68 candidate / 36 zero_stable / 200 control), 24
workers` and `done in 1132.4s`, and 17:35:32 + 1132.4 s = 17:54:24, exactly the c02 output
mtime. So the existing JSON/TSV are not a smoke run — they are the full run's own output
(304 records, 0 errors, 608 points).

That inverted the `smoke_pre_full/` step: moving those files would have destroyed the
completed result and broken c03/c04, which both open `c02_dense_multistability.json` from
the audit root (`c03_basins.py:238`, `c04_taxonomy.py:164`). The continuation therefore
took the documented "unsafe -> leave in place, validate mtime and record count" branch,
with the decision and its reasoning recorded in the manifest. The move logic is still
implemented and still checks live open FDs; it simply did not fire.

The completeness gate was not weakened. All 15 checks pass, with record and error counts
recomputed from the `records` array rather than read out of c02's own summary line.

Launched as `proteostasis-phase2-cont-20260731T180145-0500.service`, main PID 1424672,
running `ops/run_continuation.sh -> ops/continuation.py` under
`results/phase2/audit_20260731T172711-0500/ops/`. At +105 s: active/running, 25 c03
processes (1 parent + 24 workers), 67 (draw, evaluation point) pairs in flight, empty
stderr, exactly one such unit. `git status` byte-identical to session start; nothing
tracked was touched, no commit, no signal sent to any process, Phase 1 D and the
matched-benchmark run untouched.

**Q:** Recover Phase 1 experiment D without disturbing the live run: add a checkpointed,
timeout-bounded runner that is scientifically equivalent to `run_experiment_d.py`, prove
the equivalence with tests, smoke it, commit source/tests/docs only, then launch the full
60-background recovery as a persistent systemd unit with 16 concurrent subprocesses, one
BLAS thread each, and a 3600 s per-background wall limit.

**A:** The diagnosis held up under measurement. `Pool.map` returns only when every task
returns and writes nothing before that, so the live run had spent 5 h 40 min wall and
9 h 26 min CPU with a completely empty output directory. Timing every background
individually showed the median finishes in ~25 s and only four draws — 14, 19, 37, 44 —
exceed 45 s. Four pathological backgrounds out of sixty were blocking the other fifty-six.

The equivalence was made structural rather than asserted: the new runner imports
`run_experiment_d` and calls that module's own `_backgroundTask`, `_pairSummary` and
`PAIRS`, so there is no second copy of the model to drift. A test greps the runner's
source and fails if any of those is ever redefined there. Checkpoint rows are stored as
JSON, not TSV, because `writeTable`'s `%.12g` truncates float64 — formatting happens once
at merge, so the merged `interactions.tsv` is byte-identical to a direct write.

One real defect was found and fixed while writing the tests: the controller did not
forward `--n-backgrounds` to its children, so a reduced-size run would have handed each
child a different LHS matrix. The sample-hash gate would have caught it as a rejected
checkpoint rather than producing wrong numbers, but the flag was silently broken.

38 new tests, 187 total (149 pre-existing, none modified, none regressed). The strongest
are byte-for-byte equality of the merged table against a direct `writeTable`, and
float-for-float equality of the checkpointed path against a direct `_backgroundTask` call
on four measured-fast backgrounds. Timeout is tested through the real kill path, not a
mock. Smoke on backgrounds 5/48/56 reproduced identical artifact hashes across two
independent runs and resumed in 0.1 s versus 18.1 s with every checkpoint file untouched.

Committed as `ee64a3f` (source + tests) and a follow-up docs commit for
`PHASE1_D_CHECKPOINTED.md`. Launched as
`proteostasis-phase1-D-checkpointed-20260731T223225-0500.service`, MainPID 1622784, into
`results/phase1/D_checkpointed_20260731T223225-0500/`. At +40 s: active, exactly one unit,
1 controller + 16 workers all inside the unit's cgroup, 19 backgrounds already
checkpointed, stderr empty everywhere. The original unit
`proteostasis-phase1-20260731T162946-0500-D.service` (MainPID 1269455) was never
signalled and is still running with its output directory still empty.

**Note on timeout semantics:** a background that exceeds the wall limit is recorded as
`unresolved_timeout` — it contributes no rows, enters no interaction summary, and is not
counted in `n_errors`. It is a statement about the budget, not about the model, and must
not be reported as a numerical failure.

**Note on the duplicate-launch check:** the first pre-launch scan matched its own shell,
because the `-c` string of the scanning bash contains the pattern being searched for.
This is the exact pitfall documented in `results/phase1/ops/scan_live_experiments.sh`.
The check was redone as an argv-element test rather than a substring match.

## 2026-08-01

**Q:** Finish only the tracked documentation, tests and commit for the already-computed
experiment D closure. Do not recompute the 10,000-bootstrap analysis unless a validation
file fails, and do not alter raw Phase 1 results. Write `EXPERIMENT_D_FINAL.md` and
`PHASE2_CLOSURE_FINAL.md`, add a compact tracked validator, update STATUS.md, run the
suites, and commit as `Close Phase 2 interaction analysis`.

**A:** The bootstrap was not recomputed — `validation.json` reports 36/36 with
`all_passed: true`, so the existing output stood. Committed on `master` as
`Close Phase 2 interaction analysis`; full suite **249 passed**,
`scripts/phase2/check_d_closure.py` **177/177**.

A previous timed-out pass had already drafted both documents. Rather than trust them,
every descriptive claim was re-derived from the float64 checkpoint payloads, and that
found three real defects.

The largest was in the sensitivity section of *both* documents: "the chaperone-only
additive claim fails to reach a majority of backgrounds under every column" is false.
`conditional` is 0.6087 and the `usable` bounds are [0.5833, 0.6250] — both above one
half. Only the `requested` column fails to reach a majority (≤ 0.5000). The claim
actually fails because the grouped bootstrap CI for the majority fraction, [0.457,
0.739], straddles one half, and the three sensitivity columns are arithmetic
reweightings of the same 28 backgrounds rather than a second test. Conflating the two
would have implied the sensitivity analysis was doing inferential work it cannot do.

Second, `9.5e-10` as the "worst absolute discrepancy" of the `%.12g` TSV round-trip does
not reproduce from any column. The actual figures are 4.9e-12 *relative* (which is what
the validation check asserts) and, because burdens run to the censor at 1000, 4.8e-09
absolute in `excess_multiplicative` and 5.0e-10 in `excess_additive`.

Third, the draft overstated the influx-vs-nascent contrast. `nascent × total_capacity`
*is* supported on all three nulls, so "the synergy requires damage influx" and "the
interaction is driven by damage influx, not folding-machinery occupancy" were too strong.
The contrast that survives is the categorical one — synthetic collapse 13/46 against
43/46 — plus an order-of-magnitude smaller graded excess on the escape-free subset
(0.0028 against 0.0348). Scoped accordingly in the refined law, in §5.5, and in P4.

Everything else in the drafts verified exactly: the 470 lethal-single cells with Bliss
excess never negative, `mult_pred ≥ add_pred` in 100 % of both-singles-damaging cells,
`censored_12` identical to the blowup set, 40 of the 80 anti-additive cells carrying a
censored single, the 31/15 protective split with 25 and 3 additive majorities, and
Spearman −0.60. A numeric diff over the carried-forward §1–§3 confirmed zero drift from
the pending note, so every negative result survived intact; two caveats the draft had
trimmed (the weak markers' CIs touching 0.5, and the c05 design note) were restored.

The new tracked validator `scripts/phase2/check_d_closure.py` exists because the
evidence is gitignored while the prose is not, so nothing would otherwise notice the two
drifting apart. It re-checks the pinned source and artifact hashes, the counts, the nine
majority counts and verdicts, the collapse rates, the arithmetic of all 13 sensitivity
bounds, and that both documents still state those numbers. It carries no copy of the
estimator. Three negative controls were run to prove it is not vacuous, and the
absent-results branch returns `None` and exits 0 with an explicit `SKIP` rather than
reporting a success it did not verify — that branch is itself tested.

One test-gate defect was fixed: `TestDFinalOutputsAgreeWithRecomputation` skipped on
`D_final/d_final_results.json` but its `setUpClass` globs checkpoints out of the Phase 1
run root, so it would have errored rather than skipped had `results/phase1/` been cleaned
independently.

**Note:** no claim ledger exists in this repository (only `theory/SCOPE_AND_NONCLAIMS.md`),
so none was created. `OPS_SUBMISSION.md` was left uncommitted — it is an untracked
artefact of an earlier session and outside this task's scope.
