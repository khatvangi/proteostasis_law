# Gate 4 proposal — first empirical contact for the fold theorem

**Status: proposal only. Nothing here has been executed, and no outcome value in
any repository has been read in writing it.**

This is the theory side's specification of what an empirical test of
`theory/FOLD_THEOREM.md` would have to look like. It is deliberately written as a
*gate proposal* rather than an analysis, for the reason in §1.

---

## 1. Why this is a proposal and not an analysis

The empirical data that could test this sits in a sibling repository,
`../proteostasis-law-empirical-phase3`, which operates a **preregistered outcome
firewall**: predictors are built and frozen before any outcome value is read, and
a machine-readable record (`analysis/OUTCOME_FIREWALL.json`) attests to which
numbers have been looked at.

That firewall is an asset, and an ad-hoc look would destroy it. A theory that
arrived after the data, and then went looking, would produce exactly the kind of
result this whole project exists to avoid. **So the correct move is to specify the
test completely, in advance, and only then execute it.**

Three constraints follow and are not negotiable:

1. No outcome value is read while this document is being written or revised.
2. The hypothesis, the estimator, the kill criteria and the analysis set are
   fixed and hash-frozen before execution.
3. If the theory's prediction and the data disagree, that is recorded as a
   failure of the theory, not repaired by re-specification.

## 2. What the theory predicts

> **Superseded in part by §9-10.** H1 below is retained as the original framing;
> §9 shows it is not sharp enough to test, and §10 replaces it with H2.

From `theory/FOLD_THEOREM.md`, Consequences 2 and 4:

> **H1 (far-from-saturation).** At the point where viability is lost, the
> protein quality-control machinery is operating at a small fraction of its
> maximum rate — of order 0.05–0.4 in the aggregate-clearance arm — rather than
> at saturation.

The alternative it is being tested against is the standard one:

> **H0 (capacity exhaustion).** Viability is lost when quality-control flux
> saturates, i.e. the saturation fraction approaches 1.

These differ in the *direction* of a measurable quantity, not merely in a
parameter value, which is what makes the contrast worth running.

**The observable is a ratio** — flux relative to that arm's own maximum — so the
test does not require absolute copies-per-cell calibration. That is the single
most important design property, because absolute quantification is the step most
likely to fail.

## 3. The design constraint that the dilution work imposes

From Consequence 5 and D010: **growth rate is part of disposal.** In a dividing
cell, dilution is a removal channel that for most proteins outpaces proteolysis.

Therefore an experiment in which growth rate is allowed to float has *not* held
disposal capacity fixed, and its perturbation is confounded with an uncontrolled
change in the dominant removal route. Any Gate 4 design must either

- hold growth rate fixed (chemostat / turbidostat, or matched-rate sampling), or
- measure growth rate per condition and enter it as a disposal term, not as a
  covariate to be regressed away.

This constraint is a consequence of the model, not a methodological preference,
and it is the main reason a naive re-use of existing batch-culture data would not
constitute a test.

## 4. What the sibling repository already has, and what it does not

Recorded from that repository's own documentation, without reading outcomes:

| Available | Bearing on Gate 4 |
|---|---|
| Pranjic et al. 2024 — aggregation and DnaK enrichment under Ile mistranslation (norvaline) | the perturbation axis exists |
| Stikeleather et al. 2026 — measured per-codon substitution rates | the error axis exists |
| Structural and chaperone-engagement predictors, frozen at Gate 2A | per-protein covariates exist |

**What is missing for H1 is the saturation fraction itself.** Nothing in that
data set measures quality-control *flux relative to its maximum*. Aggregation
level and chaperone enrichment are burden proxies, not saturation states.

This is the honest conclusion of the exercise: **the existing data cannot test
the theory's central prediction**, however it is analysed. It can test weaker,
directional statements, and those should not be dressed up as the main test.

## 5. Two tiers, and they must not be conflated

**Tier A — executable now, weak.** Directional predictions that the frozen
Gate 2A predictors can address, e.g. whether site-resolved proteostasis demand
outperforms scalar error burden. This is the question that repository already
preregistered, and it is *not* a test of the fold theorem.

**Tier B — the actual test, requires new measurement.** H1 versus H0 needs
quality-control flux against its own maximum, measured across a mistranslation
dose series at fixed growth rate. Minimum requirements:

1. a dose series spanning survivable to non-survivable misfolding load;
2. growth rate held fixed or measured per dose;
3. a clearance-flux readout with a determinable maximum — a titratable reporter
   substrate whose degradation rate can be driven to saturation in the same
   cells, so the ratio is internal;
4. an independent viability endpoint, prespecified.

Without (3) there is no test, and (3) is the part that does not exist in any
staged data set. **It does, however, exist as a method** — see §8, which tests
that question and resolves it positively, at the cost of one substitution
(§8.1).

## 6. Kill criteria (superseded by §10.3)

> Retained as a record of the original design. The operative criteria are §10.3.

- **K1.** If the saturation fraction at the viability boundary exceeds 0.8, H1 is
  rejected and the far-from-saturation claim is retracted.
- **K2.** If the boundary cannot be located because no dose is non-survivable at
  fixed growth rate, the gate is VOID, not negative.
- **K3.** If growth rate varies by more than a prespecified tolerance across the
  dose series, the gate is VOID under §3 regardless of outcome.

## 7. What this proposal deliberately does not do

- It does not analyse, load, or reference any outcome value.
- It does not propose modifying the sibling repository, whose preregistration is
  frozen and whose gates are its own.
- It does not claim that the fold theorem is supported by anything. **No organism
  data has entered Phase 1, 2 or 3.** The one measured constant now used
  (`scripts/phase3/calibration.py`) is a growth-burden slope from the literature,
  in yeast, and it calibrates a model input rather than testing a model output.

## 8. Feasibility of the missing instrument — TESTED, and it exists

§5 item (3) was the blocker: a clearance-flux readout whose maximum is
determinable in the same cells. A literature test was run via PubMed, and the
instrument exists and is routine.

**Proteolytic queueing.** Jadhav et al. (2025) *ACS Synth Biol* 14:1062-1071,
[doi:10.1021/acssynbio.4c00612](https://doi.org/10.1021/acssynbio.4c00612),
PMID 40106229: the ClpXP-SspB complex is held at low abundance, and the
rate-limiting step "leads to proteolytic queueing, where the proteins form
waiting lines, and their overall degradation rate is slowed." SsrA-tagged
substrates are the standard handle, and synthetic circuits are routinely built on
this saturation. Critically for our purpose, they **overexpressed each component
in turn** and localised the bottleneck to **ClpX** — the ATPase — rather than
ClpP or SspB.

**A saturation-state signature.** Ogle & Mather (2016) *Phys Biol* 13:025002,
[doi:10.1088/1478-3975/13/2/025002](https://doi.org/10.1088/1478-3975/13/2/025002),
PMID 27042892: proteolytic machinery "is limited in capacity and can lead to a
bottleneck ... whereby many proteins compete ('queue') for proteolytic
resources," and the resulting inter-substrate correlations are "strongest **near
the queueing theoretic point of balance**."

That supplies all three requirements:

| requirement | how it is met |
|---|---|
| titratable substrate | inducible SsrA- (or LAA-) tagged fluorescent reporter |
| reachable saturation | proteolytic queueing, routinely engineered |
| maximum determinable **in the same cells** | relieve the queue by ClpX overexpression; and independently, the queueing-correlation peak locates the balance point |

**Gate 4 is therefore executable in principle.** The blocker identified in §5 is
removed.

### 8.1 The substitution this forces, and it must not be glossed

The accessible arm is **ClpXP-mediated degradation of soluble substrate**. In the
model that is the `s_u` term, not `s_a`.

`theory/FOLD_THEOREM.md` Consequence 2 gives medians `s_u = 0.155` and
`s_a = 0.056`. H1 is a claim about the machinery being far from saturation, and
`s_u` tests it — but it is the **less extreme** of the two arms, so this is a
weaker version of the prediction than the headline number suggests. Aggregate
clearance involves disaggregation (ClpB/DnaK) upstream of proteolysis and has no
comparable titratable handle located here.

H1 must therefore be restated for execution:

> **H1'.** At the point where viability is lost, ClpXP-mediated degradation of a
> reporter substrate is operating at a small fraction of its own maximum — of
> order 0.1-0.3 — rather than at saturation.

with K1 rewritten against `s_u` and its threshold reset accordingly. Quoting the
`s_a` figure while measuring `s_u` would be a bait-and-switch.

### 8.2 Two cautions that survive

- **The chaperone arm remains unmeasurable here.** No equivalent titratable,
  saturation-reachable handle for DnaK/GroEL flux was located. The folding side
  of the theory stays untested.
- **The reporter is a proxy for the endogenous burden.** Queueing measured on a
  synthetic substrate reports the state of the protease, which is what H1' is
  about — but it assumes the synthetic substrate and the mistranslation-induced
  burden compete for the same bottleneck. Given that ClpX is the bottleneck and
  handles misfolded protein generally, that is plausible and is itself testable
  by the correlation signature, but it is an assumption and belongs in the
  preregistration as one.

## 9. Setting K1 killed the design — the saturation test is underpowered

Deriving K1's threshold from the model, rather than choosing it, is what exposed
this. `scripts/phase3/gate4_prediction.py` computes the predicted distribution of
`s_u` at the collapse boundary, in the regime the experiment would actually sit
in (dividing cells, calibrated growth law):

| `s_u` at the boundary | p1 | p5 | p25 | **p50** | p75 | p95 | p99 |
|---|---|---|---|---|---|---|---|
| non-dividing, n=2884 | 0.0001 | 0.0010 | 0.0300 | **0.155** | 0.507 | 0.877 | 0.966 |
| **dividing, calibrated, n=269** | 0.0020 | 0.0043 | 0.0310 | **0.119** | 0.316 | 0.835 | **0.897** |

**The prediction is not sharp.** It spans nearly the whole interval. 17.8 % of
dividing draws exceed `s_u = 0.5`, and the p99 is 0.897 — against H0's prediction
of `s_u -> 1`, that leaves a separation of **0.103**.

The "6-18 % of V_max" figure in `theory/FOLD_THEOREM.md` is a **median**, and the
distribution behind it covers almost everything. A single measurement of the
saturation fraction at the boundary therefore cannot discriminate H1' from H0
except in an extreme tail, and the limit is the theory's own parameter
uncertainty — the Michaelis constants that set saturation were never measured —
not the assay.

**This design is rejected as the primary test.** Preregistering it would have
produced an underpowered study with a foreseeable null.

## 10. Restructured Gate 4

### 10.1 Primary outcome — critical slowing, a directional prediction

The robust prediction is not a magnitude but a **sign**. From
`PHASE2_CLOSURE_FINAL.md` §1.5: recovery time lengthens as the boundary is
approached, holding in **91.83 %** of parameter draws in experiment C, with
Spearman ρ = **−0.899** (n = 2272) in experiment B.

Directional predictions survive parameter uncertainty in a way magnitude
predictions do not, and this one is measurable with the same apparatus: perturb
proteostasis with a pulse, follow the reporter's return, and repeat at graded
distances from the boundary. It is also a **within-experiment relative**
comparison, so it cancels most calibration.

The sharpest available form is not "recovery time increases" but a
**parameter-free exponent**. At a saddle-node one eigenvalue passes through zero
and the normal form gives

```
tau  =  1/|lambda|  ~  (j_crit - j)^(-1/2)
```

The 1/2 follows from the bifurcation type, not from any rate constant, which is
exactly why it escapes the parameter uncertainty that killed §9's design.

**Experimentally convenient restatement:** `tau^-2` is then *linear* in dose, and
its x-intercept locates the boundary. That needs no prior knowledge of `j_crit`,
which is not measurable in advance.

> **H2 (primary).** At fixed growth rate, `tau^-2` measured across a
> mistranslation dose series is linear in dose with negative slope, and its
> x-intercept lies within the series or just beyond it.
>
> **H2-null.** `tau` stays bounded as dose increases through a series that does
> subsequently lose viability, or `tau^-2` is not linear in dose.

### 10.2 Secondary outcome — the saturation fraction, reported not tested

`s_u` at the boundary is still worth measuring and reporting, as a descriptive
quantity and as a constraint on the Michaelis constants. It is **not** a
hypothesis test, and no claim of support for the fold theorem may rest on it.

### 10.3 Kill criteria, revised

- **K1 (primary).** If the `tau^-2` versus dose regression does not have a
  negative slope with its CI excluding zero, across a series that does lose
  viability, H2 is rejected and the critical-slowing prediction is retracted.
  The fitted exponent is reported alongside; a value far from 1/2 with a good
  linear fit indicates a different bifurcation, not a failed measurement.
- **K2.** If no dose is non-survivable at fixed growth rate, the gate is VOID.
- **K3.** If growth rate varies beyond the prespecified tolerance across the dose
  series, the gate is VOID under §3 regardless of outcome.
- **K4 (descriptive, non-testing).** `s_u` at the boundary is reported with its
  interval. A value above 0.897 — the model's dividing p99 — is recorded as
  inconsistent with the model, but on its own it neither confirms nor rejects,
  because §9 shows the predicted distribution is too wide to support a test.

### 10.4 Prerequisite — DISCHARGED

Critical slowing had been established only *without* dilution, and dilution
changes both the boundary location and the Jacobian (`J - mu.I` for constant
`mu`, D010), so the eigenvalue governing recovery time is directly affected. That
had to be re-established before freezing. `scripts/phase3/gate4_slowing.py` does
it, fitting `log|lambda|` against `log(j_crit - j)` by continuation outward from
the exactly-known fold state:

| regime | networks | slope median | IQR | r² median |
|---|---|---|---|---|
| no dilution | 4 | 0.5134 | 0.501–0.539 | 0.9997 |
| calibrated `mu0 = 0.05` | 8 | 0.5080 | 0.497–0.513 | 0.9999 |
| calibrated `mu0 = 0.10` | 10 | 0.5054 | 0.499–0.511 | **1.0000** |

Overall slope median **0.5077**, r² median **1.0000**, with **86.4 %** of
networks within 0.05 of the predicted 0.5.

**The exponent is unchanged by dilution and the fits are essentially exact.**
Contrast §9: `s_u` spans nearly [0,1] across the same box, while this spans
0.497–0.513. The prerequisite is discharged and H2 stands as the primary outcome.

Two honest limits. Only 22 of 42 attempted ladders converged, so the sample is
small and selected for convergence. And the 1/2 exponent is **generic to
saddle-nodes**: confirming it supports "the boundary is a saddle-node", which is
what the fold theorem asserts the boundary *is*, but it does not by itself
select the two-state model over any other model with a fold. The alternative it
genuinely discriminates against is a smooth decline with no bifurcation, where
`tau` stays bounded.

## 12. The exponent is a ruler, not a test — and here is the test

§10.4 leaves an honest problem. The 1/2 exponent is robust *because* it is
generic to saddle-nodes, and that is precisely why confirming it selects this
model over no alternative at all. §9 leaves the mirror problem: the saturation
fraction discriminates in principle but cannot be measured against a prediction
spanning [0,1].

A useful test needs both properties. The resolution is to stop treating the
exponent as the hypothesis:

> **The exponent is the instrument.** `tau^-2` regressed on dose gives `j_crit`
> as an x-intercept with a confidence interval, without needing `j_crit` known in
> advance. That is a *ruler for locating the boundary*.
>
> **The test is whether the boundary moves** when a perturbation is applied that
> the competing models treat differently.

### 12.1 The perturbation, and why it discriminates

The model's distinctive structural claim is that rescue capacity is **conserved
and shared**: free chaperone is what remains after ordinary nascent chains,
damaged monomer and aggregate have all taken their share,

```
cf = c_tot / (1 + nu + uf/kappa_cu + af/kappa_ca)
```

so the nascent load `nu` consumes capacity while contributing **no damage influx
of its own**. Hence:

> **H3 (discriminating).** Raising the load of *perfectly-folding* protein lowers
> the tolerable mistranslation dose, even though that protein causes no damage.
>
> - **shared-capacity model (this one):** `j_crit` falls as `nu` rises.
> - **independent-handling model:** no shift.

The perturbation is orthogonal to the outcome — gratuitous expression of a
well-folding protein adds folding load without adding misfolding load — so a
positive result cannot be explained by the perturbation itself causing damage.

### 12.2 Effect size, measured rather than assumed

Phase 1 recorded this as C3 in 97.11 % of draws, but only as yes/no.
`scripts/phase3/gate4_discriminating.py` measures the size, which is what decides
executability. Over a **100x** nascent-load ladder (`nu` 0.05 → 5):

| regime | n | direction correct | monotone | fold drop in `j_crit`, median (IQR) |
|---|---|---|---|---|
| no dilution | 21 | 21/21 | 21/21 | 1.294 (1.055–1.721) |
| calibrated `mu0 = 0.05` | 23 | 23/23 | 23/23 | 1.257 (1.077–1.524) |
| calibrated `mu0 = 0.10` | 24 | 23/24 | 23/24 | 1.169 (1.073–1.458) |

**Direction correct in 67 of 68 (98.5 %)**, median shift **1.22x**, p90 3.51,
and 25 % of networks above 1.5x.

The ladder was chosen by sweeping candidates: 30x gives 1.12x, 100x gives 1.24x
with direction still correct in 100 % of converged networks, and wider ladders
(400x, 5000x) buy more effect at the cost of direction consistency (95 %, 94 %)
as some networks lose the fold entirely.

**So the direction is near-universal and the magnitude is modest.** A ~22 %
shift in boundary location is resolvable by a regression intercept but is not
resolvable by eye, and power depends on the intercept CI. That is a genuine
constraint on the design, not a formality.

### 12.3 The confound this design must defeat

Gratuitous protein expression is the standard way to raise `nu`. It is also the
standard way to *lower growth rate* — that is what Scott et al. (2010) used it
for. Growth rate is part of disposal (D010, §3).

**In batch culture the perturbation moves both the thing being varied and the
thing being measured, and the result is uninterpretable in either direction.**

The design is therefore valid only at **externally fixed growth rate** — a
chemostat or turbidostat, where dilution rate is set by the operator and
gratuitous expression changes `nu` without changing `mu`. This is not a
refinement to be traded away for convenience; a batch version of this experiment
cannot test H3, and would produce a shift under both competing models.

### 12.4 Revised criteria

- **K1 (instrument).** `tau^-2` versus dose must be linear with negative slope,
  CI excluding zero, in every arm. If not, the boundary cannot be located and the
  gate is VOID rather than negative.
- **K5 (primary, discriminating).** The fitted `j_crit` x-intercept must be lower
  in the high-`nu` arm than the low-`nu` arm, with non-overlapping CIs. Failure
  rejects H3 and with it the shared-capacity claim.
- **K6.** If measured growth rate differs between arms beyond the prespecified
  tolerance, the gate is VOID under §12.3 regardless of outcome.
- `s_u` remains reported and non-testing (§10.2).

## 13. Remaining specification before freezing

The instrument is settled (§10.4) and the discriminating test is specified
(§12). What is still unspecified, and must be fixed before any hash-freeze:

1. **Dose series.** Number of levels, spacing, and the criterion for "just beyond
   the boundary". Should be geometric in dose, since the ladder in §10.4 is.
2. **Growth-rate tolerance.** The numeric bound above which K3 voids the gate.
3. **Recovery-time estimator.** Pulse magnitude, the fitting window for the
   return, and the rule for rejecting a trace as non-recovering.
4. **Viability endpoint.** Prespecified and independent of the recovery readout.
5. **Nascent-load calibration.** How the gratuitous-expression level maps onto
   the `nu` ladder, and the evidence that the expressed protein really is
   well-folding rather than adding misfolding load of its own — on which the
   whole orthogonality argument in §12.1 rests.
6. **Growth-rate tolerance between arms** (K6), which is a different and tighter
   quantity than the within-series tolerance in K3.

None of these is a scientific choice the theory constrains, and all four are
places where post-hoc flexibility could manufacture a result — which is why they
are named here rather than left to execution.

The `s_a` arm, and with it the sharpest form of the original prediction, stays
untestable until a disaggregation-flux handle exists.
