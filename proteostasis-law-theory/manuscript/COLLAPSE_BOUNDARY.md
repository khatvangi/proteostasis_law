# Protein quality control fails long before it is saturated, and in dividing cells it fails only because damage slows growth

**Status: working draft. No organism data was used anywhere in this work.**
Every quantitative statement below is either a mathematical identity or a
property of a specific model over a specific parameter range. Nothing here has
been tested against a cell. See *Claim labels* before quoting any number.

---

## Abstract

Cells continuously produce misfolded protein and continuously dispose of it, by
refolding, by proteolysis, and — in dividing cells — by dilution into daughter
cells. Viability requires that disposal keep pace. The intuitive failure mode is
saturation: the cell collapses when its quality-control machinery is working flat
out and still cannot keep up.

We show that in a minimal finite-resource model of this system, collapse occurs
while the machinery is running at roughly **6–18 % of its maximum rate**, far
from saturation, and we derive why. Because damage influx enters only the
soluble-burden equation, and because mass balance makes the sum of the two
governing equations exactly *influx minus total removal*, the Jacobian
determinant is identically the cross product of the gradients of total removal
and of the aggregate nullcline. The saddle-node condition is therefore a
constrained criticality condition, and the critical influx equals total removal
evaluated there. The collapse boundary becomes a two-equation algebraic system
containing no reference to the load, replacing a numerical continuation sweep and
reproducing it to seven significant figures.

Decomposing the resulting margin shows that the shortfall between the collapse
boundary and the naive capacity bound is dominated by enzyme saturation state
(35.8 % of the shortfall) rather than by sequestration of machinery on substrate
(12.6 %). The mechanism is a race: aggregation is superlinear in burden while
enzymatic removal saturates, so aggregation wins well before removal is exhausted.

Extending the model to dividing cells changes the picture qualitatively. Dilution
is an unbounded disposal channel, and the analytic capacity bound does not survive
it. Under dilution at a burden-independent rate, the collapse point migrates to
ever larger aggregate burden and then ceases to exist: a cell that could divide at
full rate regardless of damage would have no discontinuous transition at all.
Restoring the physiological coupling — burden slows growth — restores a boundary
at every rate tested. **In a dividing cell the point of no return exists because
damage slows growth**, not in spite of it.

Division also changes what crossing the boundary means. Because dilution bounds
the high-burden state, the system becomes bistable: the constrained critical point
is unique without division and is not with it, the two points being the folds of
an S-shaped curve. Loss of the healthy branch is a jump to a persistent,
aggregate-laden but bounded state, and recovery requires lowering damage influx
below the level that triggered the switch. Finally, splitting critical influx into
its enzymatic and dilution shares gives an exact dimensionless decomposition that
survives division, in which the enzymatic term is nearly invariant and the
dilution share carries the variation.

The theorem survives the extension unchanged, for any dilution law.

---

## Introduction

### The problem, without jargon

A cell is a factory that never stops building proteins: chains of amino acids
that must fold into precise shapes to do their jobs. The process is error-prone.
Occasionally the wrong amino acid is inserted; more often, a correctly built chain
simply fails to fold. Misfolded protein is sticky. Left alone it clumps into
aggregates that are useless to the cell and that interfere with everything else.

Cells defend themselves with two systems. **Chaperones** bind misfolded chains and
give them another chance to fold correctly. **Proteases** destroy them outright
and recycle the amino acids. Both are made of protein, both are present in finite
amounts, and both are shared among all the cell's clients at once — including the
ordinary, correctly folding new chains that need chaperone help simply because
they are still being made.

The question is how much damage this system can absorb before it fails.

### The intuitive answer

The natural model is a bucket with a hole. Damage pours in at some rate; the
quality-control machinery drains it at some maximum rate. The cell survives while
inflow is below the drain rate and fails when it exceeds it. On this picture the
failure point is a **capacity** limit, and a cell approaching death should show
its chaperones and proteases working at full tilt.

This picture yields a clean bound, and the bound is genuinely correct as far as
it goes: total removal cannot exceed the sum of the machinery's maximum rates, so
an influx above that admits no steady state at all. In the model studied here
that bound is `c_tot + (rho_U + rho_A)·d_tot`, and it is provable.

### Why the intuitive answer is the wrong quantity

Earlier phases of this project found numerically that collapse never happens
anywhere near that bound — the median collapse influx was about 14 % of it — but
located the boundary only by simulation, cell by cell in a large parameter sweep.
That is evidence that a cliff exists in the cases examined. It is not an
explanation of where the cliff is or why.

This paper derives it. The derivation is short, requires no new data, and turns
out to be robust to a substantial extension of the model.

---

## The model

The state is two burden variables: `u`, soluble misfolded monomer, and `a`,
aggregated material, both carried in monomer-equivalent units so that transfer
between them conserves mass. Damage arrives at rate `j`. Misfolded monomer is
refolded by chaperones, degraded by proteases, and lost to aggregation by
nucleation (`alpha_n·u_f^m`, with `m > 1`) and by growth onto existing aggregate.
Aggregate is redissolved by chaperones and cleared by proteases.

Two features matter for what follows.

**Resources are conserved and shared.** Free chaperone is what is left after
ordinary nascent chains, misfolded monomer, and aggregate have taken their share.
The nascent-chain load `nu` contributes no damage at all; it acts purely by
occupying chaperone. Free protease is likewise depleted by bound substrate. These
are solved as simultaneous binding equilibria, never by substituting a total pool
into a free-resource formula.

**Removal saturates; aggregation does not.** Every enzymatic term is
Michaelis-like and approaches a ceiling. Nucleation is superlinear. This asymmetry
is the engine of everything below.

---

## Results

### 1. The collapse boundary is a constrained criticality condition

Write `G(u,a) = da/dt` and let `R(u,a)` be total removal — refolding plus soluble
degradation plus aggregate clearance.

Two structural facts:

**(i)** The influx `j` appears in `du/dt` and in no other equation. The aggregate
nullcline `{G = 0}` is therefore a curve in burden space that does not move when
the load changes.

**(ii)** The internal transfer between `u` and `a` cancels between the two
equations, so `du/dt + da/dt = j − R` exactly.

Applying the determinant-preserving row operation `row1 → row1 + row2` to the
Jacobian gives

```
det J  =  −( R_u G_a − R_a G_u )  =  −( ∇R × ∇G ) .
```

A saddle-node requires `det J = 0`, which is exactly the condition that `∇R` be
parallel to `∇G` — the Lagrange condition for a constrained critical point of
total removal on the aggregate nullcline. With (ii), the critical value is

```
j_crit  =  R(u*, a*) .
```

Neither defining equation contains `j`, so the collapse boundary is a two-equation
algebraic system in the burden coordinates alone.

**Verification.** Against the earlier phase's recorded folds, the identity holds
to a median relative error of 1.4 × 10⁻⁷, which is the finite-difference floor
rather than a property of the identity. The residual in the parallelism condition
is not zero at those states and should not be: they are bracketed approximations
whose leading eigenvalue is about −2 × 10⁻⁴. That residual correlates with the
recorded eigenvalue at **r = +0.9987**, and the one state bracketed to
`eig = −4.2 × 10⁻⁹` has a parallelism residual of 3.8 × 10⁻⁸ — three orders
tighter. Solving the two equations directly reproduces the continuation-derived
boundary to 6.7 × 10⁻⁷.

**What is *not* claimed.** The proven statement is *critical point*, not
*maximum*. Total removal rises monotonically along the branch reached by taking
the first aggregate root, and the solved boundary state lies past that branch's
turning point. Whether the relevant critical point is a maximum over the whole
nullcline is open.

### 2. Collapse occurs far from saturation, and saturation state is why

Writing `R = c_f·s_ref + rho_U·d_f·s_u + rho_A·d_f·s_a` against the capacity bound
`c_tot + rho_U·d_tot + rho_A·d_tot`, the margin between them has exactly two
sources: machinery **sequestered** on substrate (`c_f < c_tot`, `d_f < d_tot`) and
Michaelis factors **below saturation** (`s < 1`). Separating them by counterfactual
across 2884 independently sampled collapse points, and reconstructing the recorded
margin from first principles at every one (median error 1.3 × 10⁻¹³):

| at the collapse point | median |
|---|---|
| margin below the capacity bound | 0.0769 |
| refolding saturation `s_ref` | 0.175 |
| soluble degradation saturation `s_u` | 0.155 |
| aggregate clearance saturation `s_a` | 0.056 |
| share of shortfall from **saturation state** | **35.8 %** |
| share of shortfall from **sequestration** | 12.6 % |

The machinery is running at roughly 6–18 % of maximum when viability is lost.

The mechanism is the asymmetry noted in the model section. Aggregation accelerates
with burden; enzymatic removal decelerates. An accelerating process racing a
decelerating one overtakes it early, and it does so while the decelerating process
still has most of its headroom unused. The capacity bound is correct but roughly
13× too loose to predict anything on its own.

### 3. Growth dilution removes the bound, and the boundary needs feedback to exist

Dividing cells dispose of material by dilution, and for most bacterial proteins
dilution outpaces proteolysis. Adding it — `du/dt − mu·u`, `da/dt − mu·a` —
leaves the binding equilibria untouched.

**The theorem survives, for any dilution law.** Influx still enters one equation
only, and the transfer still cancels, giving `du/dt + da/dt = j − R_tot` with
`R_tot = R + mu(u,a)·(u+a)`, so `det J = −( ∇R_tot × ∇G_dil )` again. Verified at
1.2 × 10⁻¹⁰ for constant dilution and 4.7 × 10⁻¹⁰ for burden-dependent dilution;
the `mu → 0` limit reproduces the undiluted model exactly. For constant `mu` the
diluted Jacobian is exactly `J − mu·I`, so the collapse condition states that the
growth rate is an eigenvalue of the undiluted Jacobian.

**The capacity bound does not survive.** `R_tot` contains `mu·(u+a)`, unbounded in
burden. Continuation in `mu`, each solve seeded from the previous and every
accepted point satisfying `|G| < 10⁻¹⁷` and `|det J| < 10⁻¹⁷`:

| `mu` | 0 | 0.02 | 0.04 | 0.06 | 0.08 | 0.10 |
|---|---|---|---|---|---|---|
| `j_crit` | 0.1542 | 0.1738 | 0.1950 | 0.2186 | 0.2456 | **none** |
| `a*` | 0.265 | 0.318 | 0.389 | 0.496 | **0.750** | — |

The collapse state migrates to ever larger aggregate burden and then escapes. Past
that point the low-burden branch no longer terminates: burden rises smoothly with
influx, and a cell dividing at a rate independent of its damage can dilute its way
out. Real cells slow down when damaged, and adding that coupling restores a
boundary at every rate tested — under either of two functional forms for the
coupling, the hyperbolic one and the linear one implied by proteome partitioning.

> **In a dividing cell the point of no return exists because damage slows growth.**
> Dilution alone is unbounded disposal; only its throttling by the burden it
> disposes of makes viability finite.

### 4. Under division the transition is to a persistent damaged state, with hysteresis

Division changes not only *whether* there is a boundary but *what crossing it
means*. Because dilution removes material without limit, the high-burden state is
bounded rather than divergent — and the system becomes bistable.

At `mu = 0.04` there are two constrained critical points, not one: the expected
boundary at `j = 0.1950`, and a second genuine saddle-node at `j = 0.1585`
(eigenvalues −1.083 and 0). They are the two folds of an S-shaped curve. Sweeping
influx up from zero burden and then back down:

| | |
|---|---|
| healthy branch survives up to | `j = 0.194`, tips at 0.196 |
| state reached after tipping | `u = 0.079, a = 3.94` — **bounded, not divergent** |
| damaged branch survives back down to | `j = 0.160`, recovers at 0.158 |
| bistable window | **[0.160, 0.194]**, lying inside the two computed folds |

So in a dividing cell, losing proteostasis is not a blow-up. It is a switch into a
persistent, aggregate-laden but still viable state — and **getting back out
requires reducing the damage load well below the level that caused the switch**.
That asymmetry is a direct, quantitative prediction, and it is the kind of
irreversibility that a purely capacity-based account has no way to produce.

Two cautions. First, this arises from a model feature — cell division — that
earlier phases of this project did not contain, and it does not reinstate an
earlier bistability claim that was rejected on other grounds for a different
model. Second, it was found at one parameter point and has **not** been surveyed;
how often division makes such a system bistable is unknown.

### 5. A margin that survives division

The margin used above divides by the enzymatic capacity bound, which Result 3
shows is not a bound once cells divide. Splitting critical influx into its two
disposal routes repairs this exactly:

```
j_crit = C_enz · φ_enz / (1 − δ)
```

where `φ_enz` is the fraction of enzymatic capacity in use at collapse and `δ` is
the share of disposal done by division. Both are dimensionless fractions; the
identity is algebra and holds to 1.6 × 10⁻¹⁶.

The informative result is that **`φ_enz` barely moves** — 0.125 to 0.134, a ±4 %
band, as dilution rises from zero to the threshold — while `δ` grows from 0 to
0.39 and carries all the variation. Division does not change the enzymatic
condition for collapse; it multiplies the tolerable influx by `1/(1 − δ)`. The
escape in Result 3 is `δ` approaching 1.

Across 25 parameter draws, 23 lose their boundary under constant dilution, at
thresholds spanning 3.3 decades, with `δ` at the threshold having median 0.35 —
the boundary is typically lost once division does about a third of the disposal.

---

## Discussion

### What would distinguish this from the intuitive picture experimentally

The two accounts make opposite predictions about a measurable quantity:

- **Capacity exhaustion:** viability is lost when quality-control flux saturates.
- **This work:** viability is lost while that flux is a small fraction of maximum.

Usefully, the discriminating observable is a **ratio** — how close the machinery
is to its own maximum rate — rather than an absolute count of molecules per cell.
Ratios cancel most of the calibration burden that makes absolute proteome
quantification hard.

### A design constraint that follows from the dilution result

If growth rate is part of disposal, then an experiment in which growth rate is
allowed to change has not held disposal capacity fixed. Perturbation designs that
compare damaged and undamaged cells without controlling growth are measuring a
composite of the intended perturbation and an uncontrolled change in the dominant
disposal channel. This constraint is a consequence of the model, not a
methodological preference.

### Relation to the surrounding framework

This paper concerns the stability filter only. It says nothing about *why*
translation operates where it does, nothing about the evolutionary origin of any
feature of the genetic code, and nothing about a trade-off surface — the six
objectives named in the surrounding framework are not computed anywhere in this
work. The feasibility constraint derived here is the object that such a surface
would be intersected with, but that intersection has not been constructed.

---

## Limitations

1. **No organism data.** None, at any stage. Every number is a statement about
   equations.
2. **One model.** Two burden states, no gene regulation, no compartmentation, no
   oxidative damage or trafficking.
3. **No growth-burden relation has been measured.** Two functional forms are
   used and the qualitative conclusions agree under both, but neither is fitted
   to data, and form-dependent quantities differ between them.
4. **Parameter ranges were chosen adversarially wide**, to find where predictions
   fail rather than to represent any organism. Since saturation state dominates
   the margin and the Michaelis constants set saturation, the spread across
   parameter draws is partly a property of the sampling design.
5. **Bistability was found at one parameter point, not surveyed.** How often
   division makes such a system bistable is unknown, and it is the obvious next
   question.
6. **Uniqueness is settled only in part.** It holds without division and fails
   with it. Which critical point is operative is established at one parameter
   point by the sweep, not in general; whether the relevant one is a maximum
   remains open.

---

## Claim labels

Following the project convention:

| Label | Content |
|---|---|
| **Mathematical** | the determinant identity; the constrained-criticality characterisation; the capacity bound; the `mu → 0` reduction; the `J − mu·I` form |
| **Computational** | every percentage, margin and threshold — properties of this model over these ranges under this numerical protocol |
| **Empirical hypothesis** | everything about organisms. **Untested.** |

---

## Methods and reproduction

The model, the derivation and every number are in `theory/FOLD_THEOREM.md`.
Analysis code is tracked:

```
python scripts/phase3/fold_theorem.py        # theorem, margin, nested design
python scripts/phase3/dilution.py            # growth dilution, boundary survival
python scripts/phase3/boundary_structure.py  # uniqueness, bistability, thresholds
python -m pytest tests/phase3 -q             # 42 checks
```

Raw result directories are gitignored; the modules print an explicit skip and the
artefact-dependent tests skip when they are absent, while the model-level tests —
including all dilution tests — run on a clean checkout and are the ones that pin
the mathematics. Decisions D007–D014 record what changed and what was withdrawn.
