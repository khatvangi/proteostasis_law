# Status

## Phase 3: the fold is derived, not sampled

On 2026-08-02 the collapse boundary was characterised analytically. Because `j`
enters `du/dt` and nowhere else, the aggregate nullcline is a fixed curve, and
mass balance plus one determinant-preserving row operation gives the identity
`det J = -(grad R x grad G)` for total removal `R`. A saddle-node is therefore
exactly a constrained critical point of `R` on the nullcline, so

```
j_crit = R(u*, a*)          exact, no continuation sweep
```

**The fold is the constrained maximum of total removal on the aggregate
nullcline.** Statement, proof, verification and limits are in
`theory/FOLD_THEOREM.md`; decisions D007-D009 record what changed.

Verified against the Phase 1 run root: the identity holds to a median relative
error of **1.436e-07** (the central-difference floor), the parallelism residual
correlates with the recorded leading eigenvalue at **+0.9987** — showing the
residual is bracket tolerance and not a failure of the identity — and the 2x2
solver reproduces the continuation-derived folds to **6.652e-07**. `phi` rebuilds
from first principles at all **2884** folds, median absolute error **1.3e-13**.

**What sets phi.** At collapse the machinery runs at `s_ref` **0.175**, `s_u`
**0.155**, `s_a` **0.056** — roughly 6-18 % of V_max. Collapse is *not* capacity
exhaustion; superlinear nucleation overtakes sublinear saturating removal long
before removal saturates. Counterfactually the saturation deficit accounts for
**35.8 %** of the shortfall against **12.6 %** for sequestration. This supplies
the mechanism behind P2 and shows `removalCeiling`, while a correct bound, is
about 13x too loose to predict anything alone.

**One earlier reading withdrawn (D009).** A within-draw comparison over the
allocation `chi` alone gave a between/within variance ratio of 9.6, read as `phi`
being a load-invariant network constant. The properly nested design — kinetic
draws crossed with a `(nu, chi)` grid, 446 folds solved — gives **5.9**, with
per-network spreads up to **13.6x** when both load coordinates sweep, overlapping
the **8.86x** between-network spread. `phi` is network-characteristic but not
load-invariant. Holding one load coordinate gives 1.59-1.80x.

**No empirical claim.** No organism data entered Phase 3 either. The sharpened
prediction — collapse at `s_a` far below saturation, against the
capacity-exhaustion alternative — is an untested empirical hypothesis, and it
requires growth rate to be held fixed because growth rate sets `nu`.

**Growth dilution, and the coupling the theory had not stated.** Cell division
was absent from the model, and for most proteins dilution outpaces proteolysis.
`scripts/phase3/dilution.py` adds it without touching the frozen upstream model.
The theorem survives for any dilution law — `mu -> 0` reduces to the frozen model
exactly (0.0), the identity holds at 1.2e-10 (constant `mu`) and 4.7e-10
(burden-dependent), and for constant `mu` the diluted Jacobian is exactly
`J - mu.I`, so the saddle-node condition says `mu` is an eigenvalue of the
undiluted Jacobian.

The **removal ceiling does not survive**: `R_tot = R + mu.(u+a)` is unbounded in
burden, so A8 is false once cells divide. Under *constant* dilution the fold is
destroyed above a threshold — continuation gives `j_crit` rising 0.1542 → 0.2456
as `mu` goes 0 → 0.08 while `a*` diverges 0.265 → 0.750, and by `mu = 0.10` no
fold exists. Growth feedback restores it at every rate tested.

> **The collapse boundary of a dividing cell exists because burden slows growth.**

An experiment that lets growth rate float is therefore not holding disposal
capacity fixed, because growth rate is part of disposal.

**One gloss withdrawn (D011).** The proven statement is constrained *critical
point*, not *maximum*. `R` rises monotonically along the first-root branch and
the solved fold has `dG/da > 0`, so the critical point sits past the nullcline's
turning point. No number changes — everything solves `det J = 0` directly.

**Division makes the system bistable, and uniqueness fails (D012, D013).**
Enumerating the whole nullcline rather than one branch settles the uniqueness
question with a split verdict. Without dilution the curve closes (152 lower, 152
upper) and carries exactly **one** critical point. With dilution at `mu = 0.04`
there is a **second** genuine saddle-node at `u = 0.1585, a = 1.9835`
(`|G| = 2.8e-17`, `|det J| = 1.6e-17`, eigenvalues `-1.083, 0`) whose critical
influx **0.1585 lies below** the first at 0.1950.

They are the two folds of an S-curve, because dilution makes the high-burden
state bounded rather than divergent. Sweeping influx up from zero burden and back
down: the healthy branch survives to `j = 0.194` and jumps at 0.196 to a bounded
state (`u = 0.079, a = 3.94`); the damaged branch survives back to `j = 0.160`
and recovers at 0.158. The bistable window **[0.160, 0.194]** lies inside the two
computed critical points. **Losing proteostasis under division is a transition to
a persistent aggregate-laden state, not a divergence, and recovery requires
lowering influx below the level that triggered it.** This does not reinstate
`notes/REJECTED_COMPONENTS.md` item 7 — it comes from a model feature neither
that claim nor Phase 2 §2.1 contained — and it was found at one parameter point,
not surveyed.

**A margin that survives division (D014).** `j_crit = C_enz . phi_enz / (1-delta)`
exactly (identity to 1.6e-16), with `phi_enz` the enzymatic capacity in use at
collapse and `delta` the share of disposal done by division. `phi_enz` is nearly
invariant to dilution (0.125–0.134) while `delta` carries the variation, so
division multiplies tolerable influx by `1/(1-delta)` without changing the
enzymatic condition. Across 25 draws, 23 lose their boundary under constant
dilution; thresholds span 3.3 decades and `delta` at the threshold has median
0.35. Conclusions hold under both the hyperbolic and the linear
(proteome-partitioning) growth-burden form, though form-dependent values differ.

A working manuscript for the whole phase 3 result is
`manuscript/COLLAPSE_BOUNDARY.md`.

Reproduce with `python scripts/phase3/fold_theorem.py`,
`python scripts/phase3/dilution.py` and
`python scripts/phase3/boundary_structure.py`; asserted by `tests/phase3/`
(42 checks, of which 35 are model-level and run on a clean checkout).

## Phase 1 experiment D closed; Phase 2 synthesis final

On 2026-08-01 experiment D was closed scientifically and integrated into the
Phase 2 synthesis.  The two documents are `EXPERIMENT_D_FINAL.md` and
`PHASE2_CLOSURE_FINAL.md`; the latter supersedes the working note
`PHASE2_CLOSURE_PENDING_D.md` under the gitignored closure directory and
carries every one of its negative results forward unchanged.

The run analysed is the checkpointed recovery run
`results/phase1/D_checkpointed_20260731T223225-0500/`.  **60 backgrounds
requested, 58 completed, 46 usable, 12 model-unusable (`no viable state at
j_lo`), 2 unresolved timeouts (backgrounds 19 and 37), zero numerical errors,
zero process failures.**  An unresolved background exceeded the 3600 s wall
limit and was not evaluated within the budget.  It is **not a failure**, is not
counted in `n_errors`, contributes no rows, and must never be reported as one.

Validation passes **36/36**: live source, config and Latin-hypercube hashes match
what the run recorded; all 58 checkpoints pass the runner's own identity gate
with zero rejections; re-merging them reproduces `interactions.tsv` and
`backgrounds.tsv` **byte for byte**; all three nulls and all four excess columns
recompute exactly; and `_pairSummary` recomputed from the shipped TSV reproduces
every field of `summary.json`.  `46 usable x 3 pairs x 36 = 4968` cells, of which
`46 x 3 x 25 = 3450` are genuinely double.

**The raw summary's headline fractions were cell-level descriptives and are not
inference.**  Their 3450 cells are 46 parameter draws x 3 pairs x 25 correlated
cells.  Every estimate in the closure is formed with the background as the unit,
with 95 % CIs from 10,000 replicates resampling whole backgrounds (seed
20260801); the only p-value anywhere is an exact binomial sign test over
backgrounds.

Supported at the background level, on all three nulls: `influx x total_capacity`
and `nascent x total_capacity`.  **Not** supported on the primary set under the
additive or Bliss null: `influx x chaperone_only` — and its multiplicative pass
is not corroboration, because the multiplicative null is the weaker test wherever
a single perturbation is protective, which it is in 32.6 % of that pair's cells.
The most robust result is categorical and null-free: **synthetic collapse in
43 of 46 backgrounds** (95 % CI [0.848, 1.000]), at least 71.7 % of all sixty
requested draws under the most conservative denominator.

Two audit findings travel with those numbers.  `multiplicative_median_excess` in
`summary.json` is a **log-scale** median under a name that states no scale — a
labelling defect, not a sign error.  And the Bliss result is **not** internally
inconsistent: worse means *negative* excess on the survival scale, so a negative
median and a majority above 0.5 are the same statement.

One new negative result: **chaperone-only knockdown is not universally a
burden.**  In 15 of 46 usable backgrounds it *lowers* total burden in all 25
double cells, which follows from the model's own form (chaperone binding
sequesters substrate away from the protease).  Rescue capacity is therefore not a
scalar, and the law is restated accordingly.

**No empirical claim was made.**  No organism data entered Phase 1 or Phase 2.
The six falsifiable predictions the computational work now justifies are listed
in `PHASE2_CLOSURE_FINAL.md` §7 and are labelled empirical hypotheses, untested.

Analysis code is tracked (`scripts/phase2/d_final.py`), the run identity is
pinned (`scripts/phase2/D_RUN_HASHES.json`), and both are asserted by
`tests/phase2/test_d_final.py`.  Detailed outputs are gitignored under
`results/phase2/closure_20260731T220024-0500/D_final/`, so
`scripts/phase2/check_d_closure.py` is the tracked bridge across that boundary:
**177 checks** over the pinned hashes, the counts, the nine majority counts and
verdicts, the collapse rates and every sensitivity bound, plus the requirement
that both closure documents still state them.  Without the run root it prints an
explicit `SKIP` and exits 0 rather than passing silently.  Full suite:
**249 passed** under Python 3.12.11.

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

**Superseded by `PHASE2_CLOSURE_FINAL.md` §8.**  The list below is the Phase 0
ordering and is kept as a record of what was planned before Phase 1 and Phase 2
ran; items 1 and 3 are done, and the current open list is §8 of the closure.


1. Prove or numerically delimit existence and local stability regions for the two-state conserved-resource model.
2. Specify operational damage thresholds and viability readouts for an initial organism and condition.
3. Design a perturbation matrix that independently varies burden composition, translation strategy, chaperone allocation, and degradation capacity.
4. Build uncertainty-aware estimators for site-resolved substitution probabilities and damage weights.
5. Test proximity-to-boundary scaling and burden-capacity interactions in preregistered held-out conditions.
6. Compare multi-proxy predictions against scalar mistranslation, growth rate, and expression-only baselines.
7. Curate and insert primary literature citations.
8. Only after validation, consider comparative or evolutionary extensions.
