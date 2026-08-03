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

## 2. What the theory now predicts, sharply

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

## 6. Kill criteria, to be fixed before execution

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

## 9. Recommended next action

Execution is now blocked on specification, not on method. Before any data is
touched:

1. restate H1 as H1' above and reset K1's threshold against `s_u`;
2. fix the dose series, the growth-rate control (§3), and the estimator;
3. hash-freeze this document, then run.

The `s_a` arm, and with it the sharpest form of the prediction, stays untestable
until a disaggregation-flux handle exists. That limitation should appear in the
manuscript rather than be left implicit.
