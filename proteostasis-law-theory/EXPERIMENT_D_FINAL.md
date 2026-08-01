# Experiment D — final result

Phase 1 experiment D asked whether the model actually produces the
burden × capacity interaction the theory advertises (`theory/PREDICTIONS.md` #1),
and against which null. This document closes it.

The run analysed is the checkpointed recovery run
`results/phase1/D_checkpointed_20260731T223225-0500/` (launched at commit
`ee64a3f`, operational record `850726c`, described in
`PHASE1_D_CHECKPOINTED.md`). Its identity is pinned in
`scripts/phase2/D_RUN_HASHES.json` and asserted on every test run.
The detailed outputs are gitignored under
`results/phase2/closure_20260731T220024-0500/D_final/`
(`validate_d_run.py`, `run_d_final.py`, `write_validation_report.py`,
`validation.json`, `d_final_results.json`, five TSVs, `VALIDATION_REPORT.md`).
The statistics themselves live in tracked `scripts/phase2/d_final.py` and are
tested by `tests/phase2/test_d_final.py`.

## Claim labels

| Label | Meaning |
|---|---|
| **Mathematical** | An identity or arithmetic consequence of the definitions. Does not depend on the sample. |
| **Computational** | A property of a finite sample of this model under this numerical protocol. Depends on the sampling design, the integrator and the censoring rule. |
| **Empirical hypothesis** | A statement about organisms. **Nothing here tests one.** No organism data entered experiment D. |

---

## 1. What completed

| quantity | value |
|---|---|
| backgrounds requested | 60 |
| backgrounds completed | 58 |
| usable (viable baseline strictly below threshold) | **46** |
| model-unusable — `no viable state at j_lo`, all 12 | 12 |
| **unresolved timeouts — backgrounds 19 and 37** | **2** |
| numerical errors (`n_errors`) | 0 |
| process failures | 0 |
| cells written | 4968 |
| genuinely double cells summarised | 3450 |

**The two unresolved backgrounds are not failures.** A background that exceeded
the 3600 s wall limit was not evaluated within the budget. It contributes no
rows, enters no summary, and is not counted in `n_errors`. That is a statement
about the budget, not about the model, and it is never reported otherwise here.

The 12 unusable backgrounds failed a different and *structural* test: at the low
end of the influx scan the model admits no viable state at all, so no baseline
steady state exists from which to perturb. There is no interaction to measure in
those draws, not an interaction that came out null.

## 2. The 4968 cells — Mathematical

The design is a full 6 × 6 factorial per (background, pair): six burden
multipliers `(1.0, 1.2, 1.4, 1.7, 2.0, 2.5)` × six capacity multipliers
`(1.0, 0.9, 0.8, 0.7, 0.6, 0.5)`, over three pairs.

```
46 usable backgrounds × 3 pairs × 36 cells = 4968      cells written
46 usable backgrounds × 3 pairs × 25 cells = 3450      genuinely double cells
```

Exactly one level on each axis is the unperturbed 1.0, so the 11 cells per
(background, pair) with `burden_factor == 1` **or** `capacity_factor == 1` are
single perturbations, which every null reproduces by construction. Removing
them leaves 5 × 5 = 25. Verified cell by cell: every one of the 138
(background, pair) groups has exactly 36 cells and exactly 25 doubles.

## 3. Validation — Computational

`validate_d_run.py`: **36 / 36 checks pass.**

The load-bearing ones:

- the live source files, config and Latin-hypercube matrix hash to exactly what
  the run recorded (`source_hash 809b3cd7…`, `config_hash 18280a03…`,
  `sample_matrix_hash 73b6a58d…`);
- all 58 background checkpoints pass the runner's own gate — DONE marker
  matching `checkpoint.json`, payload hashes matching, and index / config /
  source / sample identity all matching; **zero rejected checkpoints**;
- re-merging those 58 checkpoints with the original `writeTable` reproduces the
  shipped `interactions.tsv` and `backgrounds.tsv` **byte for byte**;
- the three null predictions and all four excess columns recompute **exactly**,
  bit for bit, from `burden_0/1/2` and `survival_0/1/2`;
- `_pairSummary` recomputed from the shipped TSV reproduces every field of
  `summary.json`.

One precision fact matters downstream. The excess columns are differences of
nearly equal burdens — catastrophic cancellation — so the `%.12g` TSV cannot
round-trip them. The discrepancy is at the twelfth significant digit by
construction (worst *relative* deviation over the four excess columns 4.9e-12,
which is what the `tsv_reproduces_checkpoint_values_to_12_significant_digits`
check asserts), but because the burdens themselves
run up to the censor at 1000, the worst *absolute* deviation is 4.8e-09 in
`excess_multiplicative` and 5.0e-10 in `excess_additive`. All inference here is
therefore computed from the float64 checkpoint payloads rather than the TSV. This
changes no sign: zero cells flip direction between the TSV and the checkpoints
under any of the three nulls.

## 4. The nulls, re-audited from source

### 4.1 The sign conventions — Mathematical

The three nulls are stated on three different scales, so "worse than null" does
**not** have the same sign for all three. `run_experiment_d._pairSummary`
encodes this explicitly as a per-null direction:

| null | excess column | scale | worse when |
|---|---|---|---|
| additive | `burden_12 − (burden_1 + burden_2 − burden_0)` | burden | excess **> 0** |
| multiplicative | `log burden_12 − log(burden_1·burden_2/burden_0)` | log burden | excess **> 0** |
| bliss | `survival_12 − survival_1·survival_2/survival_0` | survival | excess **< 0** |

**So `bliss_frac_worse_than_null = 0.7183` together with
`bliss_median_excess = −0.0054` is not a contradiction — it is the same
statement twice.** On the survival scale, less survival is worse, so a negative
median excess and a majority of negative excesses agree. The appearance of a
paradox comes only from carrying the additive convention ("worse = positive")
across to Bliss. The code is right; the field names invite the misreading.

Verified structurally: `tests/phase2/test_d_final.py` reads the source of
`_pairSummary` and fails if either side's direction triples ever change.

### 4.2 One real labelling defect — Mathematical

`summary.json`'s `multiplicative_median_excess` (0.02964 overall) is the median
of `log_excess_multiplicative`, i.e. a **log-scale** median. The linear-scale
column `excess_multiplicative` exists and is not summarised anywhere. The name
carries no scale, and the additive field next to it *is* on the burden scale.
This is a labelling defect, not a sign error — the direction is monotone either
way — but a reader comparing 0.0110 (additive, burden units) with 0.0296
(multiplicative, log units) is comparing two different things. Every number in
this document states its scale.

### 4.3 Three limitations that bound what the nulls can show — Computational

**(a) Bliss is floored, and 13.6 % of double cells cannot show synergy at all.**
`survival = clip(1 − burden/H, 0, 1)` has a hard floor at zero. Whenever either
single perturbation is already lethal, the Bliss prediction is exactly 0 and an
equally lethal double gives excess exactly 0. Every one of the 456 exact Bliss
ties (13.2 %) is such a cell, and across the 470 cells with a lethal single the
Bliss excess is **never negative**. Those cells are structurally incapable of
registering Bliss synergy, so `bliss_frac_worse_than_null` is diluted downward.
Restricted to cells where both singles leave nonzero survival, the cell-level
fraction rises from 0.7183 to **0.8315**.

**(b) The multiplicative null is not uniformly the stricter one.** With
`x = burden_1/burden_0` and `y = burden_2/burden_0`,
`mult_pred − add_pred = burden_0·(x−1)·(y−1)`. So multiplicative demands more
than additive exactly when both singles push burden the same way — and demands
*less* when one single is protective, where it is correspondingly easier to
exceed. In **24.4 %** of double cells exactly one single is protective
(32.6 % for `influx × chaperone_only`, 31.3 % for `nascent × total_capacity`).
This is the whole reason the multiplicative fraction (0.8096) exceeds the
additive one (0.6765) for the chaperone-only pair. **A higher `frac_worse` under
the multiplicative null is not automatically stronger evidence**, and for that
pair it is weaker evidence. Restricted to cells where both singles are damaging,
`mult_pred ≥ add_pred` holds in 100 % of cells, as the identity requires.

**(c) Censoring is escape, and it is conservative — except in one direction.**
`censored_12` is exactly the `blowup` set (verified: zero exceptions). A censored
double did not have a large finite burden; it escaped, and 1000 was recorded in
place of an unbounded value. Recording the ceiling therefore *understates*
supra-additivity in the 23.7 % of doubles that escaped. The exception runs the
other way: whenever a single escaped, the double escaped too (zero exceptions),
and then `excess_additive = burden_12 − (1000 + burden_2 − burden_0) ≤ 0` is
recorded as "better than additive" for a purely arithmetic reason. That accounts
for **40 of the 80** anti-additive cells in `influx × total_capacity`. Both
directions are reported; neither is corrected.

**(d) 250 of 3450 doubles (7.2 %) were not at a steady state** when integration
stopped at `t_end = 5e4` (`status_12 == "timeout"`: final relative rate still
≥ 1e-6). Their burden is a transient value, not an equilibrium. Overall
30.99 % of doubles are non-`converged` (819 escaped + 250 still moving); for
`influx × total_capacity` it is 54.1 %.

---

## 5. Inference at the background level — Computational

**The 3450 double cells are not 3450 independent observations.** They are 46
parameter draws, each contributing 25 cells per pair that share one parameter
vector, one fold, one baseline steady state and one pair of single-perturbation
scores. The independent unit is the background. Every estimate below is formed
at that level; the only p-value computed anywhere is an exact binomial sign test
over backgrounds. Intervals are percentile 95 % CIs from 10,000 replicates
resampling **whole backgrounds** with replacement, seed 20260801.

Cell-level fractions are still reported — they are the run's own headline
numbers — but with intervals that resample backgrounds, which is roughly √25
wider than pretending the cells were independent.

### 5.1 Primary analysis set (all 25 double cells per background)

`Δ` is the background-level median of each background's own median excess.
`maj` is the fraction of the 46 backgrounds whose own majority of double cells is
worse than null.

| pair | null | Δ (95 % CI) | maj (95 % CI) | cell frac worse (95 % CI) | sign test |
|---|---|---|---|---|---|
| influx × total capacity | additive | **0.3652** [0.2213, 350.6] | **1.000** [1.000, 1.000] | 0.9304 [0.892, 0.963] | 46/0, p = 2.8e-14 |
| influx × total capacity | multiplicative (log) | **0.3348** [0.2319, 0.5460] | **0.9565** [0.891, 1.000] | 0.9078 [0.863, 0.946] | 44/2, p = 3.1e-11 |
| influx × total capacity | bliss | **−0.0644** [−0.1136, −0.0451] | **0.8913** [0.804, 0.978] | 0.7896 [0.709, 0.862] | 41/0, p = 9.1e-13 |
| influx × chaperone only | additive | 0.00356 [**−0.00022**, 0.0339] | 0.6087 [**0.457**, 0.739] | 0.6765 [0.571, 0.780] | 28/18, p = 0.184 |
| influx × chaperone only | multiplicative (log) | **0.0270** [0.0061, 0.0557] | **0.7609** [0.630, 0.870] | 0.8096 [0.712, 0.897] | 35/11, p = 5.4e-4 |
| influx × chaperone only | bliss | −0.00406 [−0.0253, **0.0000**] | 0.6087 [**0.457**, 0.739] | 0.6304 [0.512, 0.748] | 28/13, p = 0.028 |
| nascent × total capacity | additive | **0.00197** [0.00069, 0.00406] | **0.8043** [0.674, 0.913] | 0.7826 [0.665, 0.887] | 37/9, p = 4.1e-5 |
| nascent × total capacity | multiplicative (log) | **0.00859** [0.00494, 0.01449] | **0.8261** [0.717, 0.935] | 0.7800 [0.671, 0.876] | 38/8, p = 9.2e-6 |
| nascent × total capacity | bliss | **−0.00109** [−0.00319, −0.00033] | **0.8043** [0.674, 0.913] | 0.7348 [0.626, 0.835] | 37/8, p = 1.5e-5 |

Bold marks the entries that decide a verdict; the entries in bold **inside**
brackets are the ones that fail.

No number is pooled across the three pairs. The same 46 backgrounds appear in
all three, so the run's "overall" row averages three correlated designs that
answer three different questions; it is descriptive only.

### 5.2 Distribution across backgrounds, not just its median

The medians above hide real structure, and one case is not noise:

- `influx × total_capacity`, additive: per-background median excess ranges from
  2.2e-4 to 999.6 across the 46 backgrounds; the upper quartile is 998.8,
  because in a quarter of backgrounds the median double cell escaped. This is
  why the CI upper limit is 350.6 — the *direction* is unanimous (46/46) but the
  *magnitude* is meaningless as a single number. On the escape-free subset the
  same quantity is 0.0348 [0.0244, 0.0430].
- `influx × chaperone_only`, additive: the per-background fraction worse is
  sharply **bimodal** — quartiles 0.200 / 0.980 / 1.000, with 23 backgrounds at
  exactly 1.0 and a large cluster at 0.200. That bimodality has a cause, in §7.

### 5.3 Prespecified robustness sets

Two subsets condition only on the **single** perturbation scores, never on the
double outcome, so neither can select for or against an interaction:
`both_singles_damaging` (both singles raise burden) and `bliss_informative`
(both singles leave nonzero survival). A third, `uncensored`, conditions on the
double outcome and is **conservative by construction** — it removes the escapes,
which are the most strongly supra-additive cells. `clean` is the intersection of
the first and third.

| pair | null | set | n backgrounds | Δ (95 % CI) | maj (95 % CI) | passes |
|---|---|---|---|---|---|---|
| influx × total capacity | additive | clean | 43 | 0.0348 [0.0244, 0.0430] | 0.953 [0.884, 1.000] | yes |
| influx × total capacity | multiplicative | clean | 43 | 0.0742 [0.0572, 0.0836] | 0.930 [0.837, 1.000] | yes |
| influx × total capacity | bliss | informative | 45 | −0.1841 [−0.2778, −0.1356] | 0.956 [0.889, 1.000] | yes |
| influx × chaperone only | additive | clean | 31 | 0.0122 [0.0054, 0.0314] | 0.839 [0.710, 0.968] | yes |
| influx × chaperone only | multiplicative | clean | 31 | 0.0353 [0.0089, 0.0467] | 0.774 [0.613, 0.903] | yes |
| influx × chaperone only | bliss | informative | 45 | −0.00772 [−0.0349, −0.00017] | 0.667 [0.533, 0.800] | yes (marginally) |
| nascent × total capacity | additive | clean | 30 | 0.00276 [0.00122, 0.00551] | 0.833 [0.700, 0.967] | yes |
| nascent × total capacity | multiplicative | clean | 30 | 0.00791 [0.00494, 0.01707] | 0.833 [0.700, 0.967] | yes |
| nascent × total capacity | bliss | informative | 46 | −0.00223 [−0.00351, −0.00063] | 0.826 [0.717, 0.935] | yes |

The surviving background counts (43, 31, 30) are reported because these subsets
are **not** random subsamples of the 46.

---

## 6. Synthetic collapse — Computational

Synthetic collapse is a categorical interaction defined without reference to any
null: both single perturbations survivable, the combination not
(`viable_1 ∧ viable_2 ∧ ¬viable_12`). Verified to be exactly that definition,
cell by cell.

| pair | backgrounds with ≥ 1 collapse | 95 % CI | median per-background collapse rate | 95 % CI | cell rate | 95 % CI |
|---|---|---|---|---|---|---|
| **any pair** | **43 / 46 = 0.9348** | [0.848, 1.000] | 0.2267 | [0.207, 0.267] | 0.2206 | [0.190, 0.250] |
| influx × total capacity | 43 / 46 = 0.9348 | [0.848, 1.000] | 0.4000 | [0.360, 0.480] | 0.3965 | [0.348, 0.443] |
| influx × chaperone only | 36 / 46 = 0.7826 | [0.652, 0.891] | 0.2000 | [0.200, 0.241] | 0.2235 | [0.179, 0.268] |
| nascent × total capacity | 13 / 46 = 0.2826 | [0.152, 0.413] | 0.0000 | [0.000, 0.000] | 0.0417 | [0.021, 0.066] |

Only 3 of 46 backgrounds show no collapse anywhere. This is the most robust
result in experiment D: it needs no null, no scale choice, no tie convention,
and it is unaffected by the censoring and Bliss-floor problems in §4.3 — a cell
in which both singles are viable and the double escapes is a collapse under any
reading.

The three pairs are not equally susceptible. Collapse is common when total
rescue capacity is knocked down against a real damage influx, appreciably less
common under chaperone-only knockdown, and **rare when the added burden carries
no damage influx at all** (`nascent`, 4.2 % of cells, and a median background
shows none). That last contrast is informative: it says the collapse is driven
by damage influx meeting reduced clearance, not by occupancy of the folding
machinery per se.

---

## 7. Chaperone-only knockdown is protective in a third of backgrounds — Computational

The bimodality in §5.2 is not noise. The 46 usable backgrounds split cleanly:

- in **31**, chaperone-only knockdown raises burden in **all 25** double cells;
- in **15**, it *lowers* burden in **all 25** double cells.

Of the 31, 25 show an additive-majority-worse; of the 15, only 3.
Spearman correlation between a background's protective fraction and its additive
`frac_worse` is −0.60.

This direction is a consequence of the model's own form, not an anomaly. In
`scripts/proteostasis/model.py` the state variables are **total** pools and free
substrate is `uf = u / (1 + cf/κ_cu + df/κ_du)`, while degradation acts on free
substrate, `ρ_U · df · uf/(κ_u + uf)`. Chaperone binding therefore *sequesters
substrate away from the protease*. Where binding is tight, lowering `c_tot`
releases substrate to degradation and lowers total burden.

The descriptive marker agrees: `κ_cu`, the chaperone–soluble-substrate binding
constant, is by far the strongest discriminator (rank AUC 0.187; median 0.162 in
the protective group against 0.950 in the damaging group — low `κ_cu` is *tight*
binding). The next markers are much weaker (`κ_ref` 0.712, `α_n` 0.669).

**This marker is descriptive, not causal**, for the same reason the `κ_da`
marker in the Phase 2 closure is: the Latin hypercube varies all 16 parameters
independently, which licenses "the protective backgrounds sit here" and forbids
"this parameter causes it". The controlled test is a one-at-a-time scan along
`κ_cu` with the other fifteen parameters held fixed. That was not done.

---

## 8. Sensitivity to the 14 backgrounds with no interaction rows

Nothing is imputed. Three denominators answer three different questions.

`conditional` is over the 46 usable backgrounds — the population the experiment
can speak about. `usable bounds` allows the 2 unresolved backgrounds to have
been usable and to have gone either way. `requested bounds` counts all 60 draws,
treating the 12 structurally unusable ones as backgrounds in which the property
was not demonstrated.

| quantity | k / 46 | conditional | usable bounds | requested bounds |
|---|---|---|---|---|
| ≥ 1 synthetic collapse, any pair | 43 | **0.9348** | [0.8958, 0.9375] | [0.7167, 0.7500] |
| ≥ 1 collapse, influx × total capacity | 43 | 0.9348 | [0.8958, 0.9375] | [0.7167, 0.7500] |
| ≥ 1 collapse, influx × chaperone only | 36 | 0.7826 | [0.7500, 0.7917] | [0.6000, 0.6333] |
| ≥ 1 collapse, nascent × total capacity | 13 | 0.2826 | [0.2708, 0.3125] | [0.2167, 0.2500] |
| majority worse than additive, influx × total cap. | 46 | 1.0000 | [0.9583, 1.0000] | [0.7667, 0.8000] |
| majority worse than multiplicative, influx × total cap. | 44 | 0.9565 | [0.9167, 0.9583] | [0.7333, 0.7667] |
| majority worse than bliss, influx × total cap. | 41 | 0.8913 | [0.8542, 0.8958] | [0.6833, 0.7167] |
| majority worse than additive, influx × chap. only | 28 | 0.6087 | [0.5833, 0.6250] | [0.4667, 0.5000] |
| majority worse than multiplicative, influx × chap. only | 35 | 0.7609 | [0.7292, 0.7708] | [0.5833, 0.6167] |
| majority worse than bliss, influx × chap. only | 28 | 0.6087 | [0.5833, 0.6250] | [0.4667, 0.5000] |
| majority worse than additive, nascent × total cap. | 37 | 0.8043 | [0.7708, 0.8125] | [0.6167, 0.6500] |
| majority worse than multiplicative, nascent × total cap. | 38 | 0.8261 | [0.7917, 0.8333] | [0.6333, 0.6667] |
| majority worse than bliss, nascent × total cap. | 37 | 0.8043 | [0.7708, 0.8125] | [0.6167, 0.6500] |

Two asymmetries must travel with these bounds.

1. **The 2 unresolved backgrounds are not a random pair.** They are the two that
   exceeded a 3600 s wall limit whose median background finishes in ~25 s. Slow
   integration in this model is associated with stiffness near a boundary, so if
   anything they are enriched for the extreme behaviour, not representative of
   it. The bounds are honest ranges, not a probability statement.
2. **The `requested` column mixes two populations.** The 12 unusable draws have
   no baseline steady state, so "no interaction demonstrated" there is a
   statement about the model's viability geometry, not about interaction. That
   column is the most conservative possible reading and is included for
   completeness, not as the estimate.

Under every column, one conclusion is unchanged: synthetic collapse in
`influx × total_capacity` reaches at least **71.7 %** of *all sixty requested
draws* even in the worst case.

No column rescues the chaperone-only additive claim either, but the reason has to
be stated precisely. Its `conditional` and `usable` values *are* above one half
(0.6087 and [0.5833, 0.6250]); only the most conservative `requested` column
fails to reach a majority at all (≤ **0.5000**). What actually fails is the
grouped bootstrap CI for that majority fraction on the primary set,
[0.457, 0.739], which straddles one half. **These three columns are arithmetic
reweightings of the same 28 backgrounds, not a second test**, and none of them
can repair a CI that straddles the decision boundary.

---

## 9. Verdicts

A claim of non-additivity is entered only when three things hold in the **same
direction** on the primary set: the background-level median excess points the
worse-than-null way, its grouped bootstrap 95 % CI excludes zero, and the CI for
the fraction of backgrounds with a within-background majority lies strictly above
0.5. Any one of the three alone is descriptive.

| pair | null | verdict |
|---|---|---|
| influx × total capacity | additive | **supra-additive — supported** |
| influx × total capacity | multiplicative | **supra-multiplicative — supported** |
| influx × total capacity | bliss | **worse than Bliss independence — supported** |
| influx × chaperone only | additive | **not supported** on the primary set; supported only conditional on both singles being damaging (31 of 46 backgrounds) |
| influx × chaperone only | multiplicative | supported on the primary set, **but see the caveat below** |
| influx × chaperone only | bliss | **not supported** on the primary set (CI upper limit is exactly 0.0; majority CI [0.457, 0.739] straddles a half); supported marginally on the informative subset |
| nascent × total capacity | additive | **supra-additive — supported** |
| nascent × total capacity | multiplicative | **supra-multiplicative — supported** |
| nascent × total capacity | bliss | **worse than Bliss independence — supported** |

**The caveat on `influx × chaperone_only` multiplicative.** It passes while
additive fails, in exactly the pair where a single perturbation is protective in
32.6 % of cells and the multiplicative prediction therefore sits *below* the
additive one (§4.3b). That is what an arithmetic artefact looks like, and the
multiplicative pass must not be reported as independent corroboration of a
supra-additive claim for that pair. On the `clean` subset, where the ordering
identity forces `mult_pred ≥ add_pred`, both nulls pass on the same 31
backgrounds — which is the honest form of the claim.

**Where additive and multiplicative pass but Bliss does not:** that is
`influx × chaperone_only` on the primary set, stated as such above and not
smoothed over. Bliss disagrees there for two reasons that are both visible in
§4.3: its floor makes 12.2 % of that pair's cells structurally unable to show
synergy, and the survival scale saturates once burden exceeds the threshold, so
it cannot distinguish "dead" from "much more dead".

**Synthetic collapse is the more robust result**, and it is a *different* kind of
claim — categorical, null-free, and unaffected by every scale and censoring
problem above. 43 of 46 backgrounds show at least one, and the pattern across
pairs (0.93 / 0.78 / 0.28) is itself informative.

---

## 10. What experiment D does not establish

- **Nothing about organisms.** Every statement above is Computational or
  Mathematical. No organism data entered experiment D. Label: **Empirical
  hypothesis — untested.**
- **No causal parameter claim.** §7's `κ_cu` marker is a location in the sampled
  box, not a cause.
- **No claim about the 12 unusable draws.** They have no baseline; the
  experiment is silent about them.
- **No claim about backgrounds 19 and 37.** They were not evaluated within the
  budget.
- **No magnitude claim for `influx × total_capacity` additive excess.** The
  direction is unanimous; the size is dominated by escapes and is reported only
  on the escape-free subset.
- **No claim that "worse than null" means synergy where a single is protective.**
  In 24.4 % of double cells one single perturbation lowers burden, and there the
  phrase does not carry its intended meaning.

## 11. Reproduction

```
cd results/phase2/closure_20260731T220024-0500/D_final
python validate_d_run.py            # 36 checks; exit 1 if any fails
python run_d_final.py               # background-level inference, seed 20260801
python write_validation_report.py   # renders VALIDATION_REPORT.md
```

Read-only with respect to the experiment D run root. The statistics are in
tracked `scripts/phase2/d_final.py`; the run identity is in tracked
`scripts/phase2/D_RUN_HASHES.json`; both are asserted by
`python -m pytest tests/phase2/test_d_final.py`.

Because the outputs above are gitignored, the bridge between them and this
document is itself tracked:

```
python scripts/phase2/check_d_closure.py   # 177 checks; hashes, counts, headlines
```

It re-checks the pinned source and artifact hashes, the counts in §1, the
majority counts and verdicts in §5.1, the collapse rates in §6, the arithmetic of
every sensitivity bound in §8, and that this document and
`PHASE2_CLOSURE_FINAL.md` still state those numbers. On a checkout without the
run root it prints an explicit `SKIP` line and exits 0 rather than passing
silently.
