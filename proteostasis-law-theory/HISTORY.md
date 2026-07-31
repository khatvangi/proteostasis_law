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
