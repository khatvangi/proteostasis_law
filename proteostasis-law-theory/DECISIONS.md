# Decision Log

## D001 — Law is a feasibility restriction

The law is defined as membership in a feasible subset of a translation Pareto surface, determined by bounded stable proteostasis dynamics below damage thresholds.

## D002 — Damage influx remains site resolved

The canonical input retains protein abundance/translation flux, codon and condition, substituted residue identity, site criticality, and a damage function. A scalar mistranslation rate may appear only as a labeled coarse approximation.

## D003 — Minimal dynamics have two burden states

The baseline system contains misfolded monomer `U` and aggregate `A`, refolding, degradation, nucleation/growth, disaggregation, and aggregate clearance.

## D004 — Rescue pools are conserved

Free chaperone and degradation machinery are obtained from explicit total-pool balances. Chaperone occupancy includes ordinary nascent-chain folding as well as damaged substrate. Total substrate is never inserted into a free-resource binding formula as though it were free substrate.

## D005 — Claim classes are labeled

The project distinguishes definitions and theorem-level statements, conditional model consequences, falsifiable empirical hypotheses, and speculative evolutionary implications.

## D006 — Historical provenance is isolated

The manuscript contains no draft history. Recovered and rejected prior components are documented only under `notes/`.

## D007 — The fold is characterised analytically, not sampled

Because `j` enters `du/dt` only, the aggregate nullcline `{da/dt = 0}` is a
fixed curve, and mass balance gives `du/dt + da/dt = j - R` for total removal
`R`. A determinant-preserving row operation then yields the identity
`det J = -(grad R x grad G)`. A saddle-node therefore occurs exactly at a
constrained critical point of `R` on the nullcline, so `j_crit = R(u*,a*)`
exactly. This is theorem-level and replaces the Phase 1/2 practice of locating
the fold only by continuation. It also reduces fold location to a 2x2 root solve
with no sweep in `j`. Stated and verified in `theory/FOLD_THEOREM.md`.

## D008 — phi is reported with its two deficits separated

`phi = j_crit / removalCeiling` may not be reported as a bare number. It must be
accompanied by the counterfactual split into the availability deficit (rescue
machinery sequestered, `cf < c_tot`, `df < d_tot`) and the saturation deficit
(Michaelis factors below their asymptote). Across all 2884 Phase 1 folds the
saturation deficit dominates, 35.8 % against 12.6 % of the shortfall, and the
machinery at collapse runs at `s_ref` 0.175, `s_u` 0.155, `s_a` 0.056. The
mechanism claim is that superlinear nucleation overtakes sublinear saturating
removal well before removal saturates — not that capacity is exhausted.

The magnitude of `phi` and the existence of the turning point are separate
questions and must not be conflated: the counterfactuals address the magnitude
only, while the turning point requires the sequestration/aggregation runaway.

## D009 — The load-invariance reading of phi is withdrawn

An earlier decomposition of the Phase 1 hypercube gave a between/within variance
ratio of 9.6 and was read as `phi` being a load-invariant property of a network.
That comparison varied only the allocation `chi` within a draw. The properly
nested design (`nestedInvariance`: kinetic draws crossed with a `(nu, chi)` load
grid) gives 5.9, with per-network spreads up to 13.6x when both load coordinates
sweep together, overlapping the 8.86x between-network spread.

`phi` is therefore network-characteristic but **not** load-invariant, and the
"material constant" analogy is not to be reintroduced. Holding one load
coordinate fixed does give ~1.6-1.8x, which is what makes the empirical contrast
in `theory/FOLD_THEOREM.md` testable — but only under fixed growth rate, since
growth rate sets `nu`.

## D010 — Growth dilution is part of disposal, and the ceiling does not survive it

Cell division dilutes every species, and for most proteins dilution outpaces
proteolysis. `scripts/phase3/dilution.py` adds `-mu(u,a).u` and `-mu(u,a).a`
without modifying the frozen upstream model, since binding is unaffected.

D007's theorem survives for ANY dilution law: `j` still enters `du/dt` only and
the internal transfer still cancels, so `du/dt + da/dt = j - R_tot` with
`R_tot = R + mu.(u+a)`, and `det J = -(grad R_tot x grad G_dil)` again. For
constant `mu` the diluted Jacobian is exactly `J - mu.I`, so the saddle-node
condition states that `mu` is an eigenvalue of the undiluted Jacobian.

The removal ceiling does NOT survive. `R_tot` is unbounded in burden, so A8 is
false once cells divide, and under *constant* dilution the fold is destroyed
above a threshold (at base parameters, present at `mu = 0.08` with `a*` diverging
to 0.750, absent by `mu = 0.10`). A fold reappears at every rate tested once
growth slows with burden.

**The collapse boundary of a dividing cell exists because burden slows growth.**
Consequence for experimental design: a perturbation that lets growth rate float
is not holding disposal capacity fixed, because growth rate is part of disposal.

## D011 — The "constrained maximum" gloss is withdrawn

D007 proves that the saddle-node is a constrained CRITICAL POINT of total removal
on the aggregate nullcline. An earlier informal gloss called it the constrained
maximum. That is withdrawn: `R` rises monotonically along the branch reached by
taking the first root in `a`, and the solved fold state has `dG/da > 0`, so the
critical point lies past the nullcline's turning point rather than on that
branch. `lowerNullclineA` is a bracketing heuristic only and must not be
described as identifying the fold's branch. Whether the critical point is a
maximum over the whole curve is open. No number changes, because every result
solves `det J = 0` directly.

## D012 — Uniqueness holds without division and fails with it

Enumerating the whole aggregate nullcline, rather than the branch a bracketing
heuristic happens to reach, settles the question D007 left open. Without dilution
the curve closes (152 lower, 152 upper points) and carries exactly one sign
change of `det J`: the constrained critical point is unique.

With dilution it is not. At `mu = 0.04` a second, distinct critical point exists
at `u = 0.1585, a = 1.9835` with `|G| = 2.8e-17`, `|det J| = 1.6e-17` and
eigenvalues `(-1.083, 0)`, and its critical influx (0.1585) is BELOW the first
(0.1950). Any claim of a unique collapse boundary must therefore be qualified by
whether division is included.

## D013 — Under division, loss of viability is a transition, not a divergence

**(!) CAVEAT AT THE CLAIM, not in a footnote.** Every number in this entry was
computed under CONSTANT dilution (`k_mu = inf`), a regime this project has since
rejected on measured grounds. Constant dilution means growth rate cannot respond
to burden, so the regime predicts EXACTLY ZERO growth-rate loss at any aggregate
load. That is not a small idealisation: it contradicts the only dosage-resolved
measurement the project holds (D015 — Geiler-Samerotte 2011, 3.2% growth-rate
reduction at under 0.1% of protein misfolded) and it contradicts the observation
the first post-diction tried to explain (D026 — Lindner 2008, over 30%
reproductive loss). **The regime that produces the bistability is the same regime
that gets the measured quantity wrong.** Under the physiological growth laws the
result does not survive as stated: linear arrest gives no bounded high-burden
state at all, and hyperbolic feedback is monostable in four of six settings
tested (D026). D015 already withdrew the claim that bistability follows from cell
division per se and made it form-dependent; this caveat records the sharper form
— the numbers below are quoted from a regime that is unphysical in a specific,
measured way, and they may not be quoted without it.

D010 said constant dilution "destroys the fold" and an earlier phrasing said
there is then "no collapse boundary at all". Both overstate it.

With dilution, `G -> -infinity` at large aggregate whenever `alpha_g.u_f < mu`,
so the high-burden state is BOUNDED. The two critical points of D012 are the two
folds of an S-shaped curve. Sweeping influx up from zero burden and back down at
`mu = 0.04`: the low branch survives to `j = 0.194` and jumps at 0.196 to a
bounded state (`u = 0.079, a = 3.94`); the high branch survives back down to
`j = 0.160` and drops at 0.158. The bistable window [0.160, 0.194] lies inside
the two computed critical points 0.158496 and 0.195047.

So losing the low-burden branch is a transition to a persistent, bounded,
aggregate-laden state, and recovery requires lowering influx BELOW the level that
triggered the transition. Above the dilution threshold what disappears is the
discontinuity, not viability.

This does not reinstate `notes/REJECTED_COMPONENTS.md` item 7, which rejects
bistability attributed to the old one-variable ODE, nor does it overturn Phase 2
§2.1, which demoted bistability in the NON-dividing model. Bistability here comes
from a model feature neither analysis contained. It was found at the base
parameters, has NOT been surveyed across the box, and those base parameters carry
the constant-dilution defect stated at the head of this entry.

## D014 — The margin is reported as (phi_enz, delta) once division is included

`phi = j_crit / removalCeiling` divides by a quantity that stops bounding
`j_crit` under division. The replacement is exact algebra:

    j_crit = C_enz . phi_enz / (1 - delta)

with `phi_enz` the enzymatic capacity in use at collapse and `delta` the share of
disposal done by division; both dimensionless and in [0,1). The identity holds to
1.6e-16. `phi_enz` is nearly invariant to dilution (0.125-0.134 across
`mu = 0` to 0.08) while `delta` carries the variation, so division multiplies the
tolerable influx by `1/(1 - delta)` without changing the enzymatic condition.
Every `phi` reported before this decision is `phi_enz` at `delta = 0`.

Across 25 draws with a boundary at `mu = 0`, 23 lose it under constant dilution;
the threshold spans 3.3 decades and `delta` at the threshold has median 0.35.

## D015 — The growth-burden relation is calibrated to a measurement, and it retracts D013's generality

`dilution.py` and `boundary_structure.py` used guessed functional forms. A real
dosage-resolved measurement was located and is now the anchor:

  Geiler-Samerotte KA et al. (2011) PNAS 108(2):680-685,
  doi:10.1073/pnas.1017570108, PMID 21187411, retrieved via PubMed:
  "a 3.2% growth rate reduction when misfolded YFP represents less than 0.1% of
  total cellular protein."

That gives a slope of 32 per unit misfolded proteome fraction — a LOWER bound,
since the abundance figure is an upper bound — and linear arrest at a misfolded
fraction of 0.03125, an UPPER bound. It is YEAST; no equivalent dosage-resolved
bacterial measurement was found.

Two conversions are needed to use it and NEITHER was obtained: the chaperone plus
protease share of the proteome, and the refolding turnover that sets model time.
They are therefore SWEPT, not assumed, and no result is quoted at a single
calibrated point.

Three outcomes:

1. **The prior guesses were in the right regime.** Calibrated `k_mu` runs 0.31 to
   6.25 over proteome shares of 0.10 down to 0.005; the guessed 0.5 and 2.0 both
   lie inside it.
2. **A collapse boundary survives in 30 of 30 calibrated cells.** The
   "constant dilution destroys the boundary" pathology of D010 does not occur at
   any calibrated setting — it is an artifact of omitting the measured coupling.
   `phi_enz` stays in 0.072-0.147 throughout.
3. **D013's bistability does NOT survive, and is form-dependent.** Under the
   linear (measured-shape) law the down-sweep fails to settle at every influx
   tested — a runaway, not a second attractor — while the hyperbolic form settles
   on both branches with a window of (0.17, 0.19). The reason is mechanical:
   linear arrest sets `mu = 0` exactly beyond `k_mu`, switching dilution off and
   removing the bound on the high-burden state; the hyperbolic form only
   approaches zero, so dilution keeps bounding it.

**The measurement cannot adjudicate**, because it constrains only the slope at
below 0.1% misfolded — roughly three decades below the arrest burden where the
two forms diverge. Whether proteostasis collapse in a dividing cell is a
reversible switch or an irreversible runaway therefore turns on an unmeasured
property: whether growth arrest under misfolding burden is complete or merely
asymptotic. D013 must be quoted with its growth law attached, and the claim that
bistability is a consequence of cell division per se is withdrawn.

This also identifies the highest-value measurement the theory currently asks for.

## D016 — Sweeps must verify that reported states are equilibria

`hysteresisSweep` originally returned whatever the integrator reached at
`t_end`. Under the linear law that produced apparent high-burden "branches" at
`a` of order 10^4 that were still growing between sweep steps, which would have
been reported as bistability. It now checks each endpoint against the vector
field and exposes `settled_up`, `settled_down` and `all_settled`, and only influx
values where BOTH branches settled can enter a bistable window. A large state is
not an attractor.


## D017 — The Pareto surface is computed, and the optimum binds the constraint

`theory/PARETO_GEOMETRY.md` asserted a trade-off surface cut by the proteostasis
condition; no script ever computed one. `scripts/phase3/pareto.py` supplies the
minimal version: strategy `(alpha, R)` for accuracy and quality-control
investment, both paid out of one proteome, cut by the derived constraint
`j(alpha,R) < j_crit(R)`. Cost forms are stipulated, not fitted.

The throughput optimum lies exactly on the feasibility boundary,
`j/j_crit = 1.000000`, which follows from throughput being strictly decreasing in
both coordinates. A grid gave 0.8975 and that was discretisation error, not an
interior optimum — grid optima must not be reported as interior without tracing
the boundary.

Along the non-dominated front `j/j_crit` runs 0.227 to 0.965, so the framework's
"strategies sit near the boundary" implication holds only at the
throughput-maximising end. A deterministic maximiser has zero margin, so any
observed margin requires a mechanism outside this layer.

## D018 — Empirical contact is specified as a gate, never as a look

`empirical/GATE4_PROPOSAL.md` specifies what a test of the fold theorem requires,
without reading any outcome value. The sibling empirical repository operates a
preregistered outcome firewall; an ad-hoc analysis would destroy it, and a theory
that arrived after the data and then went looking is exactly what this project
exists to avoid.

The proposal records one conclusion that matters more than the design: **the
staged data cannot test the central prediction.** H1 concerns quality-control flux
relative to its own maximum, and nothing in those data sets measures a saturation
state — aggregation level and chaperone enrichment are burden proxies. Tier A
(directional predictor comparisons) is executable and is NOT a test of the
theorem; Tier B requires a clearance-flux readout with an internally determinable
maximum, which does not currently exist in any staged data.


## D019 — The Gate 4 instrument exists, and it forces an arm substitution

D018 recorded that the staged data cannot test H1 because nothing measures a
saturation state, and identified the blocker: a clearance-flux readout whose
maximum is determinable in the same cells. That question was tested against the
literature and resolves POSITIVELY.

Proteolytic queueing at ClpXP is a routine, engineered phenomenon. Jadhav et al.
(2025) ACS Synth Biol 14:1062-1071, doi:10.1021/acssynbio.4c00612, overexpressed
each component of the ClpXP-SspB complex in turn and localised the queueing
bottleneck to ClpX rather than ClpP or SspB. Ogle & Mather (2016) Phys Biol
13:025002, doi:10.1088/1478-3975/13/2/025002, show that inter-substrate
correlations peak near the queueing point of balance, which is an internal
saturation-state signature. Together these give a titratable substrate, reachable
saturation, and a maximum determinable in the same cells.

**The substitution this forces must be stated, not glossed.** The accessible arm
is ClpXP degradation of soluble substrate, which is `s_u` (median 0.155), not the
aggregate-clearance arm `s_a` (median 0.056). H1 is therefore restated as H1'
against `s_u`, with K1's threshold reset. Quoting the `s_a` figure while measuring
`s_u` would be a bait-and-switch, and `s_u` is the LESS extreme arm, so the
executable test is weaker than the headline number suggests.

**The chaperone arm stays untestable.** No comparable titratable,
saturation-reachable handle for DnaK/GroEL flux was located, so the folding side
of the theory remains without an instrument. That limitation belongs in the
manuscript.


## D020 — The saturation-fraction test is underpowered by the theory's own spread

Deriving Gate 4's K1 threshold from the model, instead of choosing it, killed the
design. `scripts/phase3/gate4_prediction.py` computes the predicted `s_u` at the
collapse boundary in the regime an experiment sits in (dividing cells, calibrated
growth law): median 0.119 but p95 0.835 and **p99 0.897**, with 17.8 % of draws
above 0.5. Against H0's `s_u -> 1` that leaves a separation of 0.103.

The "6-18 % of V_max" figure in `theory/FOLD_THEOREM.md` is a MEDIAN, and the
distribution behind it covers nearly the whole interval. A single measurement of
the saturation fraction cannot discriminate H1' from H0, and the limit is the
theory's own parameter uncertainty -- the Michaelis constants were never measured
-- not the assay. `s_u` is demoted to a reported descriptive quantity; no claim
of support may rest on it.

## D021 — Gate 4's primary outcome is the critical-slowing exponent

The robust prediction is a parameter-free exponent, not a magnitude. At a
saddle-node one eigenvalue passes through zero, so
`tau = 1/|lambda| ~ (j_crit - j)^(-1/2)`, and the 1/2 follows from the
bifurcation type rather than any rate constant.

`scripts/phase3/gate4_slowing.py` discharges the §10.4 prerequisite by
continuation outward from the exactly-known fold state: slope median **0.5077**,
r2 median **1.0000**, and 86.4 % of networks within 0.05 of 0.5. It is unchanged
by dilution (0.5134 undiluted, 0.5080 and 0.5054 at calibrated `mu0` of 0.05 and
0.10). Compare D020: `s_u` spans nearly [0,1] over the same box while this spans
0.497-0.513.

Restated for execution as `tau^-2` linear in dose, whose x-intercept locates the
boundary -- this needs no advance knowledge of `j_crit`, which is not measurable
beforehand.

Two limits travel with it. Only 22 of 42 attempted ladders converged, so the
sample is small and selected for convergence. And the 1/2 exponent is GENERIC to
saddle-nodes: confirming it supports "the boundary is a saddle-node" but does not
select the two-state model over any other model with a fold. What it genuinely
discriminates against is a smooth decline with no bifurcation.


## D022 — The exponent is the instrument; the discriminating test is the nascent-load shift

D021 left an honest gap. The 1/2 exponent is robust BECAUSE it is generic to
saddle-nodes, which is exactly why confirming it selects this model over no
alternative. D020 left the mirror gap: the saturation fraction discriminates but
cannot be measured against a prediction spanning [0,1].

The resolution is to stop treating the exponent as the hypothesis. `tau^-2`
regressed on dose yields `j_crit` as an x-intercept with a CI, needing no advance
knowledge of `j_crit` -- that is a RULER. The test is whether the boundary MOVES
under a perturbation the competing models treat differently.

**H3.** Raising the load of perfectly-folding protein lowers the tolerable
mistranslation dose, though that protein causes no damage. This follows from
rescue capacity being conserved and shared: `nu` enters only the denominator of
the free-chaperone balance, consuming capacity while contributing no influx. An
independent-handling model predicts no shift.

`scripts/phase3/gate4_discriminating.py` measures the effect size, which Phase 1
recorded only as yes/no (C3, 97.11 %). Over a 100x nascent-load ladder: direction
correct in **67 of 68 (98.5 %)**, monotone in 98.5 %, median shift **1.22x**,
p90 3.51, and 25 % of networks above 1.5x. The ladder was chosen by sweep -- 30x
gives 1.12x, 100x gives 1.24x with direction correct in 100 % of converged
networks, and 400x/5000x buy effect at the cost of direction consistency
(95 %, 94 %).

Direction is near-universal; magnitude is modest. A ~22 % shift is resolvable by
a regression intercept and not by eye, so power depends on the intercept CI.

## D023 — H3 is invalid in batch culture, and that is not a tradeable detail

Gratuitous protein expression is the standard way to raise `nu`. It is also the
standard way to lower growth rate -- Scott et al. (2010) used it precisely for
that. Growth rate is part of disposal (D010).

So in batch culture the perturbation moves both the variable and the readout, and
the result is uninterpretable in either direction: a shift would appear under
BOTH the shared-capacity and the independent-handling model. H3 is testable only
at externally fixed growth rate, in a chemostat or turbidostat where dilution
rate is set by the operator.

K6 voids the gate if measured growth rate differs between arms beyond tolerance.
This is recorded as a hard design requirement rather than a preference, and is
pinned by test so it cannot be edited away.


## D024 — The fold theorem generalises to n states

The objection that the theorem is "exact about a toy" is answered. With state
`x = (u, a, c, ...)`, influx entering only `du/dt`, and mass balance still giving
`du/dt + da/dt = j - R(x)`, the same row operation yields

    det J = -det[ grad R ; grad G ; grad C ]

so a saddle-node is a constrained critical point of R on the intersection of the
non-influx nullclines. Verified on a three-state system with sigma-32-style
chaperone control: relative error 0.000e+00 unregulated, 2.5e-11 and 2.9e-11
regulated, and `sigma0 -> 0` reproduces the frozen two-state model exactly.

The theorem is therefore a structural property of the model class, not of the
two-state reduction, and extending the model does not require re-deriving the
boundary.

## D025 — Regulation does not rescue the predictions; the theory is structural, not predictive

The hypothesis was that the quantitative predictions are weak because the model
lacks regulation, and that a controlled cell sits where its controller puts it
rather than where its kinetics do, collapsing the spread that made `s_u`
untestable (D020).

**Refuted.** The p5-p95 width of `s_u` at collapse goes from 0.8904 unregulated
to 0.9677 regulated. It widens.

One tentative observation survives and cuts against the paper's headline: the
regulated median `s_u` is 0.323 against 0.169, so control pushes the collapse
point CLOSER to saturation, partially toward the capacity-exhaustion picture the
theory argues against. Only 14 of 30 regulated networks converged against 24
unregulated; not to be quoted without a larger sample.

Two attempts to sharpen the quantitative predictions have now failed --
calibration (D015) and regulation -- while the structural core has survived every
extension, including dilution and an added controlled state (D024). The honest
position is that this is a STRUCTURAL theory: it says exactly where the boundary
is given the parameters, and that this holds for any model in the class. It does
not predict a number without measured parameters, and claims must be pitched
accordingly.


## D026 — The aging/rejuvenation post-diction FAILS, and points at spatial sequestration

The first post-diction attempted was Lindner et al. (2008) PNAS 105:3076-3081,
doi:10.1073/pnas.0708931105, PMID 18287048: E. coli under NON-STRESSED growth
accumulate aggregates in old-pole cells, losing ">30% of reproductive ability",
while new-pole progeny "exhibit rejuvenation".

**The logical point stands and is worth keeping.** Rejuvenation is only a
coherent category in a BISTABLE system; in a monostable one a daughter with less
aggregate simply relaxes back to the single attractor, so there is nothing to be
rejuvenated into. Two attractors and a separatrix are what make inheritance of a
low-burden state possible at all.

**The model does not supply that bistability.**

- CONSTANT dilution is bistable, with a 12.6-32.3 fold aggregate ratio and a
  35-82% shed fraction needed to escape. But `k_mu = inf` means growth cannot
  respond to burden, so it predicts ZERO reproductive loss -- contradicting the
  >30% that is the actual measurement.
- HYPERBOLIC FEEDBACK, the physiologically appropriate law, is MONOSTABLE in four
  of six settings tested. Where bistability appears the predicted reproductive
  loss is 48-95% (median 51%), more severe than reported.
- LINEAR arrest gives no bounded high state, so no bistability and no
  rejuvenation.

Reporting this as a success would have required quoting the constant-dilution
regime, which is precisely the regime that gets the measured quantity wrong. The
tests pin the negative so it cannot drift into a positive.

**What the failure points at.** The observation is not merely that a high-burden
state exists; it is that aggregates are SPATIALLY SEQUESTERED at a pole and
segregated asymmetrically. This model is well-mixed. Sequestration into an
inclusion body removes aggregate from the reactive pool -- changing the kinetics,
not the bookkeeping -- and is a candidate mechanism for a stable high-burden state
that does not require the growth law to do the work.

Spatial sequestration is therefore the next mechanism to add, identified by a
failed post-diction rather than a confirmed one. It is also the first candidate
in this project that arrived from an observation rather than from inspecting the
model.

## D027 — Antecedent check A1: the machinery may damage itself; the identity survives, the shortcut does not

The derivation rests on `j` entering `du/dt` and nowhere else. In a cell it does
not: chaperones and proteases are themselves translated at the same per-codon
error rate, so raising `j` degrades the machinery that clears the damage. If that
coupling broke `det J = -(grad R x grad G)`, the theorem would hold for a model
class that EXCLUDES cells. `scripts/phase3/capacity_self_damage.py` tests it with
`C_enz(load) = C_0/(1 + eps.load)` applied to both pools, `eps = 1/j_ref` swept
over four decades and NOT tuned to any biological value. `eps = 0` reproduces the
frozen model exactly (asserted at 0.0).

Two modes, and the second is the one that can actually break anything.
INFLUX (`load = j`) is the stated worry, but `j` is a parameter, so it changes the
equations without changing their dependence on the state. BURDEN
(`load = u + a`) puts capacity inside `grad R` and `grad G` themselves.

**1. The identity survives, in both modes, at machine precision.** Gradient-
normalised residual: floor 2.2e-14 at `eps = 0`; worst median 6.4e-14 (influx)
and 4.6e-13 (burden) at `eps = 100`, where capacity is down to 16.7% and 1.8% of
nominal. Log-log slopes in `eps` are 0.12 and 0.31 on quantities that never leave
the 1e-14 to 5e-13 band. Outcome 1 of the three: slope ~0 at the differencing
floor.

**There is no corrected algebraic form, and that is the finding.** The row
operation needs only (a) `j` additive in `du/dt` and absent from `da/dt`, and
(b) `du/dt + da/dt = j - R` exact. How the PARAMETERS depend on `j`, or on the
state, is irrelevant to either. The gradients are taken at fixed `j`, so
`det J = -(grad R x grad G)` is unchanged.

**2. The numerical check is weaker than it looks, and the reason is worth
recording.** `du/dt = j - R - G` holds pointwise, and the central-difference
operator is linear, so it reproduces the row operation exactly at ANY step size:
the identity carries no truncation term. Measured slope in `h` is -0.97
(burden, `eps = 100`) and -0.94 (`eps = 10`) with NO V-shaped minimum over four
decades -- pure roundoff, falling as `h` grows, exactly as predicted. Refining the
step from 1e-6 to 1e-2 drops the residual 5.6e-9 to 4.6e-13. So the check
certifies the IMPLEMENTATION preserves mass balance; the analytic argument
carries the theorem. Reporting the ladder at the repo's habitual `h = 1e-6` would
have shown a spurious rise and invited a false alarm.

**3. What the coupling DOES destroy is the shortcut.** The second half of fact
(i) was that `{G = 0}` is a fixed curve. Under self-damage it moves with the
load, so `j_crit = R(u*,a*)` stops being an evaluation and becomes a
self-consistency condition. Fold-finding grows from two equations in `(u,a)` to
three in `(u,a,j)`: `G = 0`, `det J = 0`, `R = j`. A real loss -- of the
algorithm, not the theorem.

**4. Direction: self-damage lowers the boundary substantially, and does NOT
steepen the approach.** Median `j_crit` ratio against the frozen fold falls
0.999 / 0.990 / 0.925 / 0.640 / 0.322 across the influx ladder and
0.995 / 0.949 / 0.740 / 0.324 / 0.131 across the burden ladder -- a 3x to 7.6x
reduction in tolerable influx. The critical-slowing exponent does not move:
paired over networks with both values, median damaged -0.4763 vs frozen -0.4813,
Wilcoxon p = 0.312, n = 19. **No new prediction.** Collapse under self-damage is
still a generic saddle-node; it just happens sooner.

**5. One exact new result, and it is loose.** In influx mode every removal flux
carries `1/(1 + eps.j)`, so the frozen model's necessary condition `j <= C_0`
becomes `eps.j^2 + j - C_0 <= 0`, i.e.

    j <= ( sqrt(1 + 4.eps.C_0) - 1 ) / (2.eps)   ->   sqrt(C_0/eps) for large eps

A LINEAR capacity ceiling becomes a SQUARE-ROOT one: doubling the machinery buys
only `sqrt(2)` in tolerable error rate once the machinery is itself error-prone.
Exact, and never violated -- but never binding either, with `j_crit/j_max` median
0.039-0.186 and max 0.623. Recorded as a necessary condition, not as the
boundary.

**6. Reported as a limitation: fold recovery is incomplete at large `eps`.**
Continuation with intermediate rungs recovers folds that a direct solve loses
(burden `eps = 100`: 0 of 7 direct, 2 of 7 continued), but recovery counts are
NON-MONOTONE in `eps` -- influx 7/5/6/6/4, burden 6/7/7/4/2. A genuine loss of
the boundary would be monotone, so this is continuation failure and is NOT
reported as folds disappearing. Where a fold does solve it is a real saddle-node:
`sin(grad R, grad G)` at the solved state is below 2.0e-9 throughout.

**Consequence for the manuscript.** The antecedent "influx enters one equation"
must be stated as what it actually requires -- that total influx be
state-independent and that mass balance count all outflow -- rather than as
independence of influx and capacity, which the theorem does not need. Whether
capacity is a function of the error rate is immaterial to the identity, and
material to where the boundary sits.
