# Status

## Phase 2A matched equivalence benchmark submitted

On 2026-07-31 the Phase 2A matched benchmark was gated and launched on both
hosts.  The targeted Phase 2 suite passes (`80 passed`, exit 0) and the full
canonical suite passes (`149 passed, 23 warnings`, exit 0) under Python 3.12.11.
The warnings are the pre-existing multiprocessing/fork deprecation warning.

T0, the `epsilon -> 0` right-hand-side and Jacobian identity test, **passed** all
10 checks with exit 0: RHS relative discrepancy `4.800103e-06` and Jacobian
`8.263311e-06` at `epsilon = 1e-6`, log-log slopes `1.0023` and `1.0028`.  The
observed exponent of 1 is what section 3 of
`theory/MATCHED_IMPLEMENTATION_PROTOCOL.md` requires.  T0 passing licenses the
downstream label comparison; it is not itself a scientific result.

One code defect was found and fixed, and it was a defect *for the matched
benchmark specifically*: each cell recorded only `tsv_sha256`, a byte hash of a
file whose every row carries a wall-clock `seconds` column.  That hash can never
match between boron's free arm and nitrogen's free arm even when every computed
value is identical, which is precisely the cross-host comparison the benchmark
exists to make.  A `payload_sha256` over the result columns was added alongside
it, and four tests now pin the property, including one asserting that `seconds`
is the *only* non-deterministic column.  No scientific criterion, threshold, or
mapping was changed; the T0 tolerances and the epsilon ladder are untouched.

A smoke benchmark passed before any full job was launched: correct schema,
`payload_sha256` identical across a repeat run and across `workers = 1` versus
`workers = 4`, zero numerical failures, and exact label and admissibility
agreement between the boron and free arms at `epsilon = 1e-6` under both
protocols.

Two jobs are **running**.  Boron runs the full 28-cell factorial as a
`systemd-run --user` transient unit; nitrogen runs the 14-cell free-limit
counterpart under Slurm at 1 CPU / 4 GB.  Both use n = 2000, seed 20260731, and
one numerical-library thread per worker.  Both hosts were proven to regenerate a
bit-identical sample matrix under different numpy builds.  Exact commands, unit
and job IDs, PIDs, source hashes, and verification evidence are in
`PHASE2_JOBS.md`.

**No Phase 2 result has been produced, inspected, or validated.**  The jobs were
still running when this was written and their partial output has not been
interpreted.

Two constraints govern any later reading of that output, and both are stated in
full in `theory/MATCHED_IMPLEMENTATION_PROTOCOL.md` section 6:

- **nitrogen's model is the `epsilon = 0` limit** of the boron family, not a
  different model.  Setting every substrate-binding affinity to zero turns
  boron's equations into nitrogen's term for term; the sequestration ladder
  makes that limit a single scannable scalar.
- **percentages are not comparable outside the matched benchmark.**  They require
  the same sample matrix, the same root protocol, the same named admissibility
  criterion, a stated epsilon rung, and every matched cell complete.  The older
  50,000-sample independent nitrogen sweep is **not** the matched counterpart and
  must not be reported as the matched result.

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
