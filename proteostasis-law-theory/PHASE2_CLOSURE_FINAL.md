# Phase 2 scientific synthesis — final, with experiment D integrated

This supersedes the working note
`results/phase2/closure_20260731T220024-0500/PHASE2_CLOSURE_PENDING_D.md`, which
was written while Phase 1 experiment D was still running and excluded it. Every
result in that note is carried forward unchanged, including all of its negative
results. What is new is §5 (experiment D) and the consequences in §6–§8.

Sources: Phase 1 A/B/C run root `results/phase1/run_20260731T162946-0500/`;
Phase 1 D checkpointed recovery run
`results/phase1/D_checkpointed_20260731T223225-0500/`; Phase 2A matched
benchmark `results/phase2/matched_20260731T175912-0500/` plus its nitrogen
counterpart; Phase 2B audit `results/phase2/audit_20260731T172711-0500/`; and
`results/phase2/closure_20260731T220024-0500/` including its `D_final/`
subdirectory. All result directories are gitignored.

## How claims are labelled

| Label | Meaning |
|---|---|
| **Mathematical** | An identity, invariant, or theorem-level consequence. Verified symbolically or to machine precision; does not depend on how the parameter space was sampled. |
| **Computational** | A property of a finite sample of a specific model under a specific numerical protocol. Depends on the sampling design and the solver; a percentage here is a property of that design, not of biology. |
| **Empirical hypothesis** | A statement about organisms. **Nothing in Phase 1 or Phase 2 has tested one.** No organism data was used anywhere in either phase, experiment D included. |

That last row is still the single most important caveat in this document. Every
number below is a statement about a model, not about a cell.

---

## 1. What the theory has survived

### 1.1 Internal mathematical consistency — Mathematical

Experiment A, 11/11 checks passed, `all_passed: true`:

| check | metric | tolerance |
|---|---|---|
| A1 dimensional consistency (+ negative control rejecting bad units) | 0.0 | 0.0 |
| A2 nondimensionalisation is an exact rescaling | 1.31e-13 | 1e-06 |
| A3 free-pool solution unique (least = greatest fixed point of the monotone map) | 3.13e-13 | 1e-10 |
| A4 pool-closure mass balance | 9.18e-14 | 1e-10 |
| A5 ODE mass balance `d(u+a)/dt = influx − refold − deg_u − deg_a` | 3.75e-16 | 1e-11 |
| A6 analytic Jacobian vs central differences | 8.31e-08 | 1e-05 |
| A7 nonnegative orthant invariant (300 trajectories, 0 negative terminations) | 0.0 | 0.0 |
| A8 no bounded state above the removal ceiling `j > c_tot + (ρ_U+ρ_A)·d_tot` | 0.0 | 0.0 |
| A9 blind and continuation equilibria agree (14 influx values) | 3.93e-14 | 1e-04 |
| A10 deterministic rerun | 0.0 | 0.0 |

A3 and A8 are the load-bearing ones. A3 makes the conserved-pool closure
well-posed (the free-pool map has a *unique* fixed point, so "the" equilibrium
is not an artefact of which root the solver happened to find). A8 is a genuine
theorem-level statement: above the removal ceiling no bounded state exists at
all, which is the hard edge the law's viability condition refers to.

### 1.2 The ε → 0 limit identity — Mathematical

T0, 10/10 checks, exit 0. The free (nitrogen) model form is the ε = 0 face of
the boron family, not a different model:

```
RHS      relative discrepancy at ε = 1e-6 : 4.800103e-06 , log-log slope 1.0023
Jacobian relative discrepancy at ε = 1e-6 : 8.263311e-06 , log-log slope 1.0028
```

Slope 1 is first-order convergence, which is what the term-by-term flux identity
in `theory/MATCHED_IMPLEMENTATION_PROTOCOL.md` §3 requires. This licenses
comparing labels across the two forms; it is not itself a result.

The model-form gap is not merely O(ε), it is quantitatively the
sequestered-substrate fraction the mapping predicts in closed form,
`1 − 1/(1 + ε·σ₀)` with `σ₀ = 1/(1+n_load) + 1`. Row by row on 2000 samples,
median absolute deviation of observed from predicted:

| ε | 1e-6 | 1e-3 | 1e-2 | 0.1 | 0.3 | 1 | 2 |
|---|---|---|---|---|---|---|---|
| predicted median gap | 1.3e-06 | 1.333e-03 | 0.013171 | 0.117832 | 0.286113 | 0.572056 | 0.727878 |
| observed median gap | 1.3e-06 | 1.333e-03 | 0.013154 | 0.116966 | 0.283485 | 0.566242 | 0.721294 |
| median abs deviation | 9.6e-15 | 9.6e-09 | 9.2e-07 | 7.1e-05 | 3.9e-04 | 1.3e-03 | 1.9e-03 |

The deviation grows with ε because the closed form is evaluated at u = a = 0 and
is a low-burden approximation — a property of the evaluation point, not a defect
in the mapping.

### 1.3 Cross-host implementation equivalence — Computational

Two independent implementations, two library stacks (python 3.12.11 / numpy
2.2.6 / scipy 1.14.0 against python 3.10.12 / numpy 1.26.4 / scipy 1.15.2), same
7 source files by SHA-256, same sample matrix by SHA-256:

- **28,000 / 28,000 rows agree on the operational label.** Diagonal confusion
  matrix in every one of 14 cells; zero `numerical_failure` anywhere.
- **112,000 / 112,000 admissibility booleans agree** (four criteria × 28,000).
- **Zero one-sided missing roots** — never once did one host find an equilibrium
  the other did not.
- Equilibria, burdens, eigenvalues and residuals agree to **≤ 2.44e-14
  absolute**, against a nearest decision boundary of **1.04e-04** — a margin of
  10¹⁰–10¹².

The payload hashes differ on 14/14 cells, and this is *not* a discrepancy: six of
the seven sampled-input columns pass through `np.exp(np.log(...))`, which the two
numpy builds round differently in the last bit on ~22 % of rows. `n_load`, the
only linear-transform coordinate, is bit-identical on all 28,000 rows — the
control that localises the divergence to the inverse transform. Full argument in
`MATCHED_CROSSHOST_AUDIT.md`.

### 1.4 Solver protocol is not driving any conclusion — Computational

Within boron's 2×2 factorial, the two root protocols (dense 9×9 log multistart +
Radau confirmation vs four fixed linear guesses + DOP853 return check) find **the
same root to ~1 ULP** (median relative difference ≈ 2.5e-16) at every ε and in
both model forms. Label disagreements: 0–3 rows per 2000 (≤ 0.15 %), all of them
*existence* misses. No outcome falls below 99 % agreement anywhere.

This matters because the multistability results below rest on root-finding. The
factorial shows the answer is not an artefact of the search strategy.

### 1.5 Fold structure and critical slowing — Computational

Experiment B (325 cells, structured grid):

- a fold exists in **325/325 cells** (`frac_cells_with_fold = 1.0`);
- median fold-to-ceiling ratio **0.1378** — collapse happens at ~14 % of the hard
  removal ceiling, i.e. the *dynamical* boundary is well inside the
  *thermodynamic* one;
- the two-state → zero-state pattern holds in **98.77 %** of cells; a third state
  appears in only 0.92 %;
- nascent-load monotonicity: **1.0** (325/325);
- critical slowing down, n = 2272: eigenvalue rate vs distance-to-fold Spearman
  **ρ = +0.893**, recovery time vs distance **ρ = −0.899**, both p = 0.0.

Experiment C (5000 Latin-hypercube draws over 16 parameters, deliberately wide
ranges chosen to find where the predictions *fail*):

| prediction | fraction | n |
|---|---|---|
| C1 fold exists | **57.68 %** | 2884 / 5000 |
| C2 collapse below the removal ceiling | **100 %** | 2884 / 2884 |
| C3 nascent load shrinks the viable window | **97.11 %** | 2760 / 2842 |
| C4 second stable attractor (raw) | 2.36 % | 68 / 2884 |
| C4 second attractor confirmed | 2.15 % | 62 / 2884 |
| C5 critical slowing down | **91.83 %** | 2023 / 2203 |
| C6 interior allocation optimum | **15.94 %** | 433 / 2716 |

C2 at 100 % is the strongest survival: in every one of 2884 draws with a fold,
the fold sits strictly below the removal ceiling. That is the geometric content
of the law — viability is lost at a *dynamical* boundary before the resource
budget is exhausted.

C5's ρ ≈ 0.9 in B and 91.8 % in C is the one prediction with an obvious
experimental analogue (recovery time lengthens near the boundary), and it is the
most promising route to an eventual empirical test.

---

## 2. What failed, or came out much weaker than advertised

Nothing in this section changed when experiment D closed.

### 2.1 Multistability is rare and was over-counted — Computational

This is the clearest negative result of Phase 2, and it went against the theory's
own preliminary count.

Phase 1 C reported 68 draws with a second stable equilibrium and 62 "confirmed".
The Phase 2B dense audit (c02 → c03 basins → c04 taxonomy) re-audited **every**
one of the 68 candidates exhaustively and reclassified them:

| class | candidates (n=68) | controls (n=200) | zero-stable (n=36) |
|---|---|---|---|
| confirmed multistable | **57** | 1 | 0 |
| locally stable, no basin | 6 | 1 | 0 |
| root-finder artifact | **5** | 0 | 0 |
| single attractor | 0 | 198 | 36 |
| threshold only | 0 | 0 | 0 |
| numerical failure | 0 | 0 | 0 |

So **62 → 57**: five of the reported second states were root-finder artifacts and
six were locally stable with no reachable basin — stable to a linearisation but
not reachable from the biologically relevant initial set, which is exactly the
condition the law's viability criterion requires.

Headline rate: **57/5000 = 1.14 % of all draws**, or **57/2884 = 1.98 % of
fold-evaluable draws**, 95 % CI **[1.50 %, 2.55 %]**.

The control arm is what makes this honest in both directions. 200 draws from the
single-attractor pool were audited blind; 1 turned out to be multistable, a
false-negative rate of 1 % (95 % CI [0.12 %, 3.57 %]), implying **3.4 to 99.1
missed multistable draws** in the single-attractor pool. So 57 is a *lower* bound
with a wide upper uncertainty, and the control was a seeded sample rather than a
census — that asymmetry must be carried forward, not dropped.

**Conclusion**: bistability is a real but marginal feature of this model, not a
generic one. Any part of the theory that leaned on hysteresis or on "two
coexisting proteostasis regimes" as a typical outcome is not supported. The law
does not need bistability — the fold in C1/C2 is a saddle-node collapse, and that
is a different and much better-supported claim.

### 2.2 Interior allocation optimum largely fails — Computational

C6 holds in only **15.94 %** (433/2716) of draws, and experiment B's structured
grid found an interior optimum in **0 %** of its 325 cells
(`allocation_interior_optimum_frac = 0.0`). The prediction that there is a
generically interior optimal chaperone/degradation split is **not supported** and
the corollary is **withdrawn**. In most of the sampled box the optimum is at a
boundary of the allocation simplex.

Experiment D adds an independent reason to distrust any interior-optimum
intuition: §5.7 shows that reducing chaperone alone *lowers* total burden in a
third of viable backgrounds, so "more chaperone is better up to a point" is not
even locally true across the sampled box.

### 2.3 The fold itself is not universal — Computational

C1 holds in 57.68 % of draws. Nearly **half** the sampled parameter box has no
fold at all in the scanned influx range — those draws simply have no collapse
boundary to be near. The law's geometry applies where a fold exists; it is not a
property of every parameterisation of the model. The 42.32 % without a fold are
not a failure of the model, but they *are* a limit on how universally the
geometric picture can be stated.

Experiment D independently reproduces a version of this limit from a different
direction: 12 of its 58 completed backgrounds (20.7 %) had **no viable state at
all** at the low end of the influx scan, so no baseline existed to perturb.

### 2.4 What did **not** fail but is often misreported

The two model forms differ by up to **72 % in total-coordinate burden** at ε = 2.
That is not a failure — it is the predicted sequestered fraction (§1.2), and in
*free* coordinates the two forms agree to ~1 ULP at every rung. But it means a
burden percentage is meaningless without stating which coordinate it is in.
Nothing in `adm_*_free` ever disagrees; the only disagreements are in
`adm_*_total`, at ≤ 1 % of rows, and they are **strictly one-way** (the boron
form is the one that becomes inadmissible, 0 rows in the reverse direction across
all 28,000 comparisons).

The one real dynamical difference: the boron form relaxes **more slowly** in
essentially every row (1751–1974 of 1753–1974), monotonically in ε, up to a
**3.6× slowdown at ε = 2**, with no eigenvalue sign change. Sequestration buffers
the rate without moving the free-coordinate fixed point.

---

## 3. The c05 parameter analysis — descriptive, not causal

`c05_parameters.py` asked which of the 16 sampled parameters mark the draws that
the audit confirmed as multistable. Base rate in its merged label: **58/2884 =
2.011 %**.

Rank-AUC effect sizes, `audit_confirmed` label, sorted by distance from 0.5:

| feature | rank AUC | 95 % CI | CI excludes 0.5 | Cohen's d | median (positive) | median (negative) |
|---|---|---|---|---|---|---|
| **log10_kappa_da** | **0.7400** | [0.6812, 0.7959] | yes | **0.849** | 0.4398 | −0.6147 |
| log10_rho_U | 0.3620 | [0.3003, 0.4319] | yes | −0.477 | −0.4169 | −0.1053 |
| log10_kappa_du | 0.3771 | [0.3049, 0.4577] | yes | −0.425 | −1.1088 | −0.5005 |
| log10_rho_A | 0.6220 | [0.5508, 0.6925] | yes | 0.428 | 0.1258 | −0.5623 |
| m | 0.4051 | [0.3402, 0.4715] | yes | −0.328 | 2.1338 | 2.2616 |
| log10_alpha_n | 0.4059 | [0.3439, 0.4708] | yes | −0.311 | −1.2471 | −0.6994 |
| log10_kappa_ref | 0.4170 | [0.3327, 0.5028] | no | −0.325 | −0.8030 | −0.2870 |
| log10_kappa_cu | 0.4285 | [0.3489, 0.5067] | no | −0.246 | −1.0332 | −0.5045 |

Cross-validated permutation importance agrees on the ordering: `log10_kappa_da`
**0.0817 ± 0.0291**, an order of magnitude above the next terms (`log10_alpha_n`
0.0250 ± 0.0140, `log10_kappa_du` 0.0189 ± 0.0140, `log10_kappa_a`
0.0129 ± 0.0074).

**The finding**: confirmed multistability is enriched primarily by the
**aggregate-degradation binding and clearance parameters** — `kappa_da`
(aggregate–protease binding) dominant, with `kappa_du`, `rho_A` and `rho_U` next,
and nucleation-related terms (`alpha_n`, `m`) contributing a weaker signal. High
`kappa_da` and high `rho_A`, low `rho_U`, low `kappa_du`, low `alpha_n` mark the
multistable corner of the box.

**This is descriptive, and the four reasons must travel with it:**

1. **It is a statement about location in the sampled box, not about causation.**
   The Latin hypercube varies all 16 parameters *independently*. That is exactly
   the design that licenses "multistability sits here" and exactly the design
   that forbids "this parameter causes multistability".
   `c05_parameters.py`'s own design note says the same thing.
2. **The predictive signal is weak in absolute terms.** Logistic-L2 ROC AUC
   0.786 ± 0.014, average precision **0.097** against a base rate of **0.020** —
   a ~4.8× lift, but the great majority of high-`kappa_da` draws are *not*
   multistable. `hist_gbdt` is no better (ROC AUC 0.755 ± 0.008, AP 0.089).
3. **The positive class is 58 draws.** Every CI above is built on that, and the
   weaker markers' intervals nearly touch 0.5 (`kappa_ref` and `kappa_cu` do not
   exclude it at all).
4. **It inherits §2.1's asymmetry.** The audit's control arm implies 3.4–99.1
   further multistable draws sitting unlabelled in the single-attractor pool. If
   those are not distributed like the 57, the markers shift.

A causal claim would require re-running the model along `kappa_da` with the other
fifteen parameters held fixed. That was not done and is the obvious follow-up.
The same caveat applies verbatim to experiment D's `kappa_cu` marker (§5.7).

---

## 4. Phase 2A matched benchmark and the ε ladder

Carried forward unchanged from §1.2–§1.4 above. The one operational rule that
governs any later reading: **percentages are not comparable outside the matched
benchmark.** They require the same sample matrix, the same root protocol, the
same named admissibility criterion, a stated ε rung, and every matched cell
complete. The older 50,000-sample independent nitrogen sweep is **not** the
matched counterpart.

---

## 5. Experiment D — the burden × capacity interaction

Full detail in `EXPERIMENT_D_FINAL.md`; only the load-bearing content is repeated
here. Analysis code: tracked `scripts/phase2/d_final.py`, tested by
`tests/phase2/test_d_final.py`; outputs under
`results/phase2/closure_20260731T220024-0500/D_final/`.

### 5.1 What completed, and one word that must not drift

60 backgrounds requested, **58 completed**, of which **46 usable** and **12
model-unusable** (`no viable state at j_lo`, all 12). **2 unresolved timeouts**
(backgrounds 19 and 37). **Zero numerical errors, zero process failures.**

An unresolved background exceeded the 3600 s wall limit and was not evaluated
within the budget. It is **not a failure** — it contributes no rows, enters no
summary, and is not counted in `n_errors`. That is a statement about the budget,
not about the model, and it is not reported otherwise anywhere in this closure.

`46 × 3 pairs × 36 factorial cells = 4968` cells written; removing the 11 cells
per (background, pair) with a factor at 1.0 leaves `46 × 3 × 25 = 3450` genuinely
double cells. Validation: **36/36 checks pass**, including byte-identical
re-merge of the shipped tables from the 58 checkpoints and exact recomputation of
all three nulls.

### 5.2 The unit of inference was wrong in the raw summary — Computational

The run's own headline numbers — 0.7965 worse than additive, 0.8325 worse than
multiplicative, 0.7183 worse than Bliss, 0.2206 synthetic-collapse rate — are
**cell-level descriptives over 3450 cells that are not 3450 independent
observations**. They are 46 parameter draws, each contributing 25 cells per pair
that share a parameter vector, a fold, a baseline steady state and a pair of
single-perturbation scores. Every estimate in this closure is formed with the
background as the unit; intervals are percentile 95 % CIs from 10,000 replicates
resampling whole backgrounds, seed 20260801. The only p-value computed anywhere
is an exact binomial sign test over backgrounds.

No number is pooled across the three pairs: the same 46 backgrounds appear in all
three, so the run's "overall" row averages three correlated designs answering
three different questions.

### 5.3 The Bliss sign question — Mathematical

`bliss_frac_worse_than_null = 0.7183` with `bliss_median_excess = −0.0054` is
**not a defect and not a contradiction**. The three nulls sit on three scales,
and `_pairSummary` carries a per-null direction: worse means excess **> 0** for
additive (burden scale) and multiplicative (log-burden scale), and excess **< 0**
for Bliss (survival scale, where less survival is worse). A negative median and a
majority of negative excesses are the same statement. The apparent paradox comes
only from carrying the additive convention across to Bliss.

One genuine **labelling** defect was found: `multiplicative_median_excess` is the
median of `log_excess_multiplicative`, i.e. a log-scale quantity reported under a
name that states no scale, next to an additive field that is on the burden scale.
Direction is unaffected; comparability is not. Every number in this closure
states its scale.

### 5.4 Supported interaction claims — Computational

A claim enters only when three things hold in the same direction on the primary
set: the background-level median excess points the worse-than-null way, its
grouped bootstrap 95 % CI excludes zero, and the CI for the fraction of the 46
backgrounds with a within-background majority lies strictly above 0.5.

| pair | additive | multiplicative | Bliss |
|---|---|---|---|
| **influx × total capacity** | **supported** — Δ 0.3652 [0.221, 350.6], majority 46/46 [1.00, 1.00] | **supported** — Δ 0.3348 [0.232, 0.546], majority 44/46 [0.891, 1.000] | **supported** — Δ −0.0644 [−0.114, −0.045], majority 41/46 [0.804, 0.978] |
| **influx × chaperone only** | **not supported** — Δ 0.00356 [**−0.00022**, 0.0339], majority 28/46 [**0.457**, 0.739] | supported — Δ 0.0270 [0.0061, 0.0557], majority 35/46 [0.630, 0.870] — **but see below** | **not supported** — Δ −0.00406 [−0.0253, **0.0000**], majority 28/46 [**0.457**, 0.739] |
| **nascent × total capacity** | **supported** — Δ 0.00197 [0.00069, 0.00406], majority 37/46 [0.674, 0.913] | **supported** — Δ 0.00859 [0.00494, 0.01449], majority 38/46 [0.717, 0.935] | **supported** — Δ −0.00109 [−0.00319, −0.00033], majority 37/46 [0.674, 0.913] |

**The chaperone-only multiplicative pass is not corroboration.** With
`x = burden_1/burden_0`, `y = burden_2/burden_0`, the identity
`mult_pred − add_pred = burden_0·(x−1)·(y−1)` means the multiplicative null is
the *stricter* one only when both singles push the same way, and the *weaker* one
where a single is protective. Chaperone-only knockdown is protective in 32.6 % of
that pair's cells (§5.7), so multiplicative passing while additive fails is
exactly what an arithmetic artefact looks like. On the subset where both singles
are damaging, `mult_pred ≥ add_pred` holds in 100 % of cells and **both** nulls
pass on the same 31 of 46 backgrounds — which is the honest form of the claim.

**Additive and multiplicative pass but Bliss does not** for
`influx × chaperone_only`, stated plainly rather than smoothed over. Bliss cannot
resolve it there for two structural reasons: `survival = clip(1 − burden/H, 0, 1)`
is floored at zero, so 12.2 % of that pair's cells have a Bliss prediction of
exactly 0 and cannot register synergy at all; and the survival scale saturates,
so it cannot distinguish "dead" from "much more dead".

**Magnitudes are not claimed for `influx × total_capacity` additive.** The
direction is unanimous across all 46 backgrounds, but the excess is dominated by
escapes (43.9 % of that pair's doubles escaped and were censored at the ceiling),
which is why the CI reaches 350.6. On the escape-free, both-singles-damaging
subset the same quantity is 0.0348 [0.0244, 0.0430] on 43 backgrounds.

### 5.5 Synthetic collapse is the more robust result — Computational

Synthetic collapse — both singles survivable, the combination not — is a
categorical interaction defined without reference to any null. It is therefore
untouched by the scale, tie and censoring problems above: a cell in which both
singles are viable and the double escapes is a collapse under any reading.

| pair | backgrounds with ≥ 1 collapse | 95 % CI | cell rate | 95 % CI |
|---|---|---|---|---|
| **any pair** | **43 / 46 = 0.9348** | [0.848, 1.000] | 0.2206 | [0.190, 0.250] |
| influx × total capacity | 43 / 46 = 0.9348 | [0.848, 1.000] | 0.3965 | [0.348, 0.443] |
| influx × chaperone only | 36 / 46 = 0.7826 | [0.652, 0.891] | 0.2235 | [0.179, 0.268] |
| nascent × total capacity | 13 / 46 = 0.2826 | [0.152, 0.413] | 0.0417 | [0.021, 0.066] |

Only 3 of 46 backgrounds show none. The gradient across pairs is itself the
finding: collapse is common when a genuine damage influx meets a reduced
conserved rescue pool, less common under chaperone-only knockdown, and **rare
when the added load carries no damage influx at all** (nascent). The
*categorical* collapse is driven by damage influx meeting reduced clearance
rather than by occupancy of the folding machinery; the graded interaction is
still supported for the nascent pair on all three nulls (§5.4), just far
smaller.

### 5.6 Sensitivity to the 14 backgrounds with no interaction rows

Nothing is imputed. `conditional` is over the 46 usable backgrounds — the
population the experiment can speak about. `usable bounds` lets the 2 unresolved
backgrounds have been usable and gone either way. `requested bounds` counts all
60 draws, treating the 12 structurally unusable ones as backgrounds in which the
property was not demonstrated.

| quantity | conditional | usable bounds | requested bounds |
|---|---|---|---|
| ≥ 1 synthetic collapse, any pair | 0.9348 | [0.8958, 0.9375] | [0.7167, 0.7500] |
| ≥ 1 collapse, influx × total capacity | 0.9348 | [0.8958, 0.9375] | [0.7167, 0.7500] |
| ≥ 1 collapse, nascent × total capacity | 0.2826 | [0.2708, 0.3125] | [0.2167, 0.2500] |
| majority worse than additive, influx × total capacity | 1.0000 | [0.9583, 1.0000] | [0.7667, 0.8000] |
| majority worse than additive, influx × chaperone only | 0.6087 | [0.5833, 0.6250] | [0.4667, 0.5000] |
| majority worse than additive, nascent × total capacity | 0.8043 | [0.7708, 0.8125] | [0.6167, 0.6500] |

Two asymmetries travel with these. The 2 unresolved backgrounds are **not a
random pair** — they are the two that exceeded a 3600 s limit where the median
background finishes in ~25 s, and slow integration here is associated with
stiffness near a boundary, so if anything they are enriched for extreme
behaviour. And the `requested` column mixes two populations, because the 12
unusable draws have no baseline steady state at all; "no interaction
demonstrated" there is a statement about viability geometry, not about
interaction.

Under every column, synthetic collapse in `influx × total_capacity` reaches at
least **71.7 % of all sixty requested draws**.

No column rescues the chaperone-only additive claim, and the reason must be
stated precisely rather than as "it fails everywhere". Its `conditional` and
`usable` values *are* above one half (0.6087 and [0.5833, 0.6250]); only the
most conservative `requested` column fails to reach a majority at all
(≤ 0.5000). What fails is the grouped bootstrap CI for that majority fraction on
the primary set, [0.457, 0.739], which straddles one half. **These columns are
arithmetic reweightings of the same 28 backgrounds, not a second test.**

### 5.7 A new negative result: chaperone-only knockdown is not universally a burden — Computational

The 46 usable backgrounds split cleanly and without intermediates: in **31**,
chaperone-only knockdown raises burden in all 25 double cells; in **15** it
*lowers* burden in all 25. Of the 31, 25 show an additive majority worse; of the
15, only 3. Spearman correlation between a background's protective fraction and
its additive `frac_worse` is −0.60.

The direction follows from the model's own form. State variables are **total**
pools; free substrate is `uf = u / (1 + cf/κ_cu + df/κ_du)` and degradation acts
on free substrate, `ρ_U · df · uf/(κ_u + uf)`. Chaperone binding sequesters
substrate away from the protease, so where binding is tight, lowering `c_tot`
releases substrate to degradation and lowers total burden.

The descriptive marker agrees: `κ_cu`, the chaperone–soluble-substrate binding
constant, is by far the strongest discriminator (rank AUC 0.187; median 0.162 in
the protective group against 0.950 in the damaging group — low `κ_cu` is *tight*
binding). The next markers are far weaker (`κ_ref` 0.712, `α_n` 0.669). This is
**descriptive, not causal**, for exactly the reason §3 gives; the controlled test
is a one-at-a-time scan along `κ_cu`, and it was not done.

**Consequence for the theory**: "reducing chaperone capacity increases burden"
is *not* a general property of this model, and the law must not be stated as
though rescue capacity were a single scalar whose reduction always hurts.
Chaperone and degradation capacity are not interchangeable, and a proportional
knockdown of both (`rescue_total`) behaves qualitatively differently from a
chaperone-only knockdown of the same nominal size.

---

## 6. The refined canonical law

The Phase 0 statement in `theory/LAW_STATEMENT.md` survives. Phase 1/2 justify
tightening four things and dropping one.

> **Proteostasis Law (refined).** For environment `e`, a translation strategy `u`
> is viable only if the coupled finite-resource proteostasis dynamics
> `dx/dt = F(x; u, e, θ)` possess an **attracting** bounded invariant state that
> is **reached from the biologically relevant initial set** and is contained in
> the admissible region `D = {x : h_j(x) < H_j}`.
>
> Where such a state exists and a fold exists, viability is lost at a
> **saddle-node collapse strictly inside the removal ceiling**, not at the
> ceiling itself, and approach to that collapse is accompanied by **critical
> slowing down**.
>
> Damage influx and **conserved rescue capacity** are **not separable**: raising
> influx while reducing the conserved rescue pool drives the system past that
> collapse **supra-additively**, and does so at combinations where each
> perturbation alone is survivable. The coupling is **much stronger when the added
> load carries a damage influx**: a nascent-load increase with no misfolding
> influx does still interact supra-additively with capacity reduction, but it
> produces synthetic collapse in only 13 of 46 backgrounds against 43 of 46, and
> a median per-background excess an order of magnitude smaller on the escape-free
> subset (0.0028 against 0.0348). And rescue capacity is **not a scalar** —
> chaperone and degradation capacity are not interchangeable.

What changed and why:

1. **"Attracting and reached from the initial set" is load-bearing, not
   decorative.** The c04 audit found 6 of 68 candidates *locally stable with no
   basin*. A linearised-stability criterion would have counted them; the
   reachability criterion correctly does not.
2. **The collapse boundary is dynamical, not thermodynamic.** C2 = 100 %
   (2884/2884), median fold-to-ceiling ratio 0.1378. A8 establishes the ceiling
   as a hard outer bound; the fold is where viability actually ends.
3. **Critical slowing down is promoted toward the central falsifiable
   prediction.** B: ρ = +0.893 / −0.899, n = 2272, p = 0.0. C: 91.83 % of 2203.
   It is the prediction with the most direct experimental readout.
4. **The burden × capacity coupling is added, in the restricted form the
   evidence supports** (new). Supported on all three nulls for
   `influx × total_capacity` and `nascent × total_capacity`; **not** supported on
   the primary set for `influx × chaperone_only` under the additive or Bliss
   null. The categorical statement — synthetic collapse at survivable singles —
   is the robust part, at 43/46 backgrounds.
5. **Rescue capacity is explicitly not a scalar** (new, §5.7). A third of viable
   backgrounds *improve* under chaperone-only knockdown.
6. **Bistability/hysteresis is demoted out of the law.** 1.14 % of draws (§2.1).
7. **Scope is narrowed explicitly**: the geometric statement applies where a fold
   exists, which is 57.68 % of the sampled box (§2.3), and where a viable
   baseline exists at all, which failed in 12 of D's 58 completed draws. Outside
   that region the law reduces to the bare existence-of-attractor condition.
8. **The interior-optimum corollary remains withdrawn** (§2.2: 15.94 % in C, 0 %
   in B, and now §5.7 against it from a second direction).

Claim status of the refined law:

- **Mathematical**: the definition, the necessity of an attracting reachable
  bounded state, the removal-ceiling bound (A8), the ε → 0 identity, and the
  null-ordering identity `mult_pred − add_pred = burden_0·(x−1)·(y−1)`.
- **Computational**: fold-below-ceiling universality, critical slowing, the
  nascent-load direction, the rarity of multistability, the ε-ladder behaviour,
  the burden × capacity supra-additivity and synthetic collapse, and the
  protective chaperone-knockdown regime — all conditional on this model, these
  ranges, this protocol.
- **Empirical hypothesis**: **everything about organisms.** Untested. No organism
  data entered Phase 1 or Phase 2.

---

## 7. The empirical predictions now justified

These are **Empirical hypotheses** and none has been tested. They are stated in
the sharpest falsifiable form the computational work supports, with the model
quantity each corresponds to.

**P1 — Critical slowing near a proteostasis collapse boundary.** Recovery time
from a small proteostatic displacement lengthens as a cell is driven toward
collapse, and does so monotonically in the distance to the boundary. *Model
support*: B ρ = −0.899 (n = 2272), C 91.83 % of 2203. *Falsified by*: recovery
time flat or decreasing with proximity in a system that does subsequently
collapse.

**P2 — Collapse before capacity saturation.** Viability is lost while chaperone
and protease capacity are still measurably unsaturated, not when they are
exhausted. *Model support*: C2 100 % of 2884; median fold at 13.8 % of the
removal ceiling. *Falsified by*: collapse coinciding with saturation of the
rescue machinery.

**P3 — Synthetic lethality of misfolding load with total rescue-capacity
reduction.** There exist doses of a misfolding-inducing insult and of a combined
chaperone + protease reduction that are each individually tolerated and jointly
lethal. *Model support*: synthetic collapse in 43/46 backgrounds; 39.7 % of
`influx × total_capacity` double cells; ≥ 71.7 % of backgrounds even counting
every unusable and unresolved draw against it. *Falsified by*: a broad dose grid
in which no such combination exists.

**P4 — Synthetic lethality is driven by damage influx, not by folding-machinery
occupancy.** Increasing nascent-chain load that carries no misfolding influx
produces a much weaker interaction with capacity reduction than a genuine
misfolding insult at the same nominal load multipliers — weaker, but not absent:
the graded interaction is still supported on all three nulls. *Model support*:
synthetic collapse 4.2 % of cells and 13/46
backgrounds for `nascent × total_capacity`, against 39.7 % and 43/46 for
`influx × total_capacity`; the median background shows none. *Falsified by*:
equal synthetic lethality from a non-misfolding translational load.

**P5 — Chaperone and protease capacity are not interchangeable, and chaperone
knockdown alone can be protective.** In systems where chaperone–substrate binding
is tight, reducing chaperone availability *lowers* aggregate misfolded burden by
releasing substrate to degradation; a proportional reduction of both arms does
not. *Model support*: 15 of 46 backgrounds protective in all 25 cells, marked by
low `κ_cu` (rank AUC 0.187); `influx × chaperone_only` fails the additive and
Bliss tests on the primary set where `influx × total_capacity` passes all three.
*Falsified by*: chaperone knockdown monotonically increasing aggregate burden
across binding regimes. **This is the sharpest new prediction from experiment D,
and it is the one most likely to be wrong** — it rests on a model in which the
burden coordinate counts chaperone-bound substrate.

**P6 — Interaction scale matters and must be prespecified.** Whether a
burden × capacity interaction reads as synergistic depends on the null and on
whether either single perturbation is protective. *Model support*: the identity
in §5.4, and 24.4 % of double cells with exactly one protective single.
*Consequence for design*: an experiment testing P3–P5 must prespecify its null,
report whether each single perturbation raised or lowered the readout, and treat
the categorical synthetic-lethality endpoint as primary.

---

## 8. What is still open

1. **Any empirical contact whatsoever.** Phase 2 closed the computational gates.
   It did not open an empirical one. P1–P6 are the entry points.
2. **Causal tests of the two markers** — a matched one-at-a-time scan along
   `kappa_da` (§3) and along `kappa_cu` (§5.7), not a hypercube.
3. **A census rather than a sample of the single-attractor pool**, to replace the
   3.4–99.1 missed-multistable interval with a number.
4. **Backgrounds 19 and 37**, which exceeded the 3600 s budget. Resolving them
   requires a longer wall limit or a stiffer-tolerant protocol, and would replace
   the sensitivity bounds in §5.6 with point estimates over 48 usable draws.
5. **Steady-state completeness in experiment D.** 250 of 3450 double cells
   (7.2 %) had not reached a steady state at `t_end = 5e4`; their burden is a
   transient. A longer horizon would say whether the interaction estimates move.
6. **A direct test of whether the burden coordinate should count machinery-bound
   substrate.** P5 depends on it, and §1.2 already shows the total- and
   free-coordinate readings diverge by up to 72 %.

## 9. Reproduction

```
# experiment D closure (this document's §5)
cd results/phase2/closure_20260731T220024-0500/D_final
python validate_d_run.py            # 36 checks; exit 1 if any fails
python run_d_final.py               # background-level inference, seed 20260801
python write_validation_report.py

# earlier phase 2 closure work (§1-§3)
cd results/phase2/closure_20260731T220024-0500
python closure_audit.py
python validate_outputs.py

# the tracked assertions
python scripts/phase2/check_d_closure.py   # 177 checks: hashes, counts, headlines
python -m pytest tests -q                  # 249 tests
```

All read-only with respect to every Phase 1 and Phase 2 raw result. `results/` is
gitignored, so none of those directories is in the repository; the statistics,
the run identity and the assertions are tracked in `scripts/phase2/d_final.py`,
`scripts/phase2/D_RUN_HASHES.json`, `scripts/phase2/check_d_closure.py` and
`tests/phase2/test_d_final.py`.

`check_d_closure.py` is what keeps §5 honest across the gitignore boundary: it
holds the counts, the nine majority counts and verdicts, the collapse rates and
every sensitivity bound in this document to the shipped `D_final/` output, and it
fails if either closure document stops stating them. Absent the run root it
prints an explicit `SKIP` and exits 0, so a clean checkout neither fails nor
reports a success it did not verify.
