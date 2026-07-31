# Status

## Phase 1 canonical tests repaired; experiments running

On 2026-07-31 the two invalid numerical tests were repaired without changing
the theory, experiment configurations, or manuscript.  The full canonical
suite passes: `69 passed, 19 warnings in 66.36s` under Python 3.12.11.  The
warnings are the existing Python multiprocessing/fork deprecation warning.

Phase 1 experiments A-D are running.  One instance of each was launched on
2026-07-31 at 16:29:46 -0500 on `boron` with the unchanged Phase 1 configs,
as `systemd-run --user` transient units, writing to the clean run root
`results/phase1/run_20260731T162946-0500/`.  Two earlier duplicate launch
attempts left no live host processes; their empty output directories were
moved to `results/phase1/quarantine_preclean_20260731T162755-0500/` rather
than deleted.  Exact commands, unit names, PIDs, resources, config hashes,
logs, and output locations are in `JOBS.md`.

**No result has been produced, inspected, or validated.**  The jobs are
running; nothing in this file asserts that any Phase 1 prediction has been
tested, confirmed, or falsified.  Read `JOBS.md` for the operational record
and wait for the runs to complete before drawing any scientific conclusion.

The repaired scripts, tests, and configs are committed on branch `master`.
The repository has no remote configured, so nothing was pushed.

## Phase 0: theory reconstruction completed

The canonical law, Pareto geometry, site-resolved damage influx, conserved-resource dynamics, scope, falsifiable predictions, empirical program, and a substantive first-pass manuscript are now defined consistently. Discredited components are quarantined in the rejection ledger.

## Next scientific tasks, in priority order

1. Prove or numerically delimit existence and local stability regions for the two-state conserved-resource model.
2. Specify operational damage thresholds and viability readouts for an initial organism and condition.
3. Design a perturbation matrix that independently varies burden composition, translation strategy, chaperone allocation, and degradation capacity.
4. Build uncertainty-aware estimators for site-resolved substitution probabilities and damage weights.
5. Test proximity-to-boundary scaling and burden-capacity interactions in preregistered held-out conditions.
6. Compare multi-proxy predictions against scalar mistranslation, growth rate, and expression-only baselines.
7. Curate and insert primary literature citations.
8. Only after validation, consider comparative or evolutionary extensions.
