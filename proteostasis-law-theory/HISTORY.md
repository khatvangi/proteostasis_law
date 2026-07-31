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
