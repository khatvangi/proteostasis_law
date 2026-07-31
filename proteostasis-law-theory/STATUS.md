# Status

## Phase 1 canonical tests repaired; experiments running

On 2026-07-31 the two invalid numerical tests were repaired without changing
the theory, experiment configurations, or manuscript.  The full canonical
suite passes: `69 passed, 19 warnings in 60.62s` under Python 3.12.11.  The
warnings are the existing Python multiprocessing/fork deprecation warning.

Experiments A-D were launched separately with the unchanged Phase 1 configs
after the passing suite.  Slurm was unavailable (`sinfo` and `squeue` could not
create the controller stream socket), so the fallback is supervised
`nohup setsid`.  Current canonical PIDs are A=743, B=744, C=745, and D=746;
all four were alive after 15 seconds with empty stderr.  Exact commands,
resources, config hashes, logs, and output locations are in `JOBS.md`.

The scientific changes could not be committed in this execution environment:
`.git` is mounted read-only, and `git commit` failed while creating
`.git/index.lock`.  The unchanged repository HEAD is
`73e1c0ab341530f7b2d369c459008674897f6287`; it is not a commit containing
these repairs.  Nothing was pushed.

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
