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

**(!) Empirical motivation withdrawn (D031).** The bistability of D012/D013 had
one pull towards biology: rejuvenation appeared to require two attractors. It does
not. Lindner's old-pole cell inherits a physical inclusion body at every division,
a continuously renewed perturbation rather than an attractor, and a MONOSTABLE
model reproduces the measured lineage difference. D012 and D013 are therefore
model properties under a growth law the project has disqualified, and in the
manuscript they appear as such. They are not a biological claim about how
proteostasis is lost.


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
the first post-diction tried to explain (D026 — Lindner 2008, a 1.2-1.8%
aggregate-attributable growth deficit; see D028 for why this is not the "over
30%" an earlier reading of that abstract reported). **The regime that produces the bistability is the same regime
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
accumulate aggregates in old-pole cells, while new-pole progeny "exhibit
rejuvenation".

**(!) The number quoted in this entry was misread and is corrected in D028.** The
abstract's ">30% of the loss of reproductive ability" is a SHARE of the aging
effect, not the effect; the full text gives the effect as 3.95 +/- 0.5% and the
aggregate share as ~30-40%, so the measured quantity is **1.2-1.8%**, not 30%.
This entry's verdict is unaffected and in fact strengthened: the losses it called
"more severe than reported" are too severe by 27x to 92x, not by a hair.

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

## D028 — PREREGISTRATION: spatial sequestration, written before the run

This entry is committed BEFORE `scripts/phase3/sequestration.py` is run, per
`notes/POSTDICTION_PROTOCOL.md` rule 2. The reason is specific rather than
ceremonial: the first Lindner pass produced exactly the tidy quantitative match
it was hoping for, under a regime later rejected. That is a
researcher-degrees-of-freedom hazard, and the fix is procedural.

### The observation

Lindner AB, Madden R, Demarez A, Stewart EJ, Taddei F (2008) PNAS
105(8):3076-3081, doi:10.1073/pnas.0708931105, PMID 18287048, retrieved via
PubMed. Verbatim from the abstract:

  "This accretion is associated with >30% of the loss of reproductive ability
  (aging) in these cells relative to the new-pole progeny, devoid of parental
  inclusion bodies, that exhibit rejuvenation."

Condition: E. coli, NON-STRESSED growth, time-lapse lineage tracking with
IbpA-tagged inclusion bodies.

### The corresponding model quantity

`1 - mu_high / mu_low`, the growth rate of the aggregate-laden attractor against
the growth rate of the low-aggregate attractor. This is a ratio BETWEEN TWO
LINEAGES, matching the paper's "relative to the new-pole progeny". D026 already
computed this quantity correctly and it is carried over unchanged.

Bistability is a PRECONDITION, not the prediction. Rejuvenation is only a
coherent category if two attractors and a separatrix exist; in a monostable
system a daughter with less aggregate relaxes straight back.

### The band, both edges — DATA-DERIVED (amended before the run)

**The original band [0.30, 0.60] was wrong by a factor of about thirty, and the
error was a misreading of the abstract.** ">30% **of** the loss of reproductive
ability" is a SHARE of the aging effect, not the aging effect. The full text
(PMC2268587, via PubMed) settles it:

  "[Delta(GR_old - GR_new)]mean/GR_mean = -3.95 +/- 0.5%"

  Table 1, Population 1 (332 cell pairs), units 1e-2 min^-1:
      GR_old = 3.54 +/- 0.02,  GR_new = 3.69 +/- 0.02   -> deficit 4.07%

  "the fraction of the growth rate decrease (aging) associated with the presence
   of the aggregate, Agg/(Agg + Pole) is ~30-40%"

So the total old-vs-new-pole growth deficit is 3.95 +/- 0.5%, and the
AGGREGATE-ATTRIBUTABLE part of it is 30-40% of that. The model has no pole-age
mechanism other than aggregate, so `1 - mu_high/mu_low` maps to the
aggregate-attributable share, and the band is their product:

    MATCH  if  0.0104 <= (1 - mu_high/mu_low) <= 0.0178

    lower = 0.30 x (0.0395 - 0.0050) = 0.01035
    upper = 0.40 x (0.0395 + 0.0050) = 0.01780

Both edges are now data-derived, which is what rule 2 asks for. The paper's
denominator is `GR_mean` rather than `GR_new`; converting costs at most 2%
relative (GR_mean/GR_new = 3.61/3.69 = 0.978), far inside the +/-0.5 percentage
point standard error, so it is absorbed rather than propagated.

**Disclosure, because it matters here.** This amendment was made AFTER a first
sequestration scan had been run and inspected. It is recorded rather than hidden,
and the direction is checkable: under the old band the in-band cells carried
losses 0.482 and 0.508, and under the new band **every one of them falls out**.
The re-anchoring is strictly MORE stringent and cannot have been chosen to
manufacture a pass. No cell moved into the band.

**Consequence for D026.** Its verdict that a predicted 48-95% loss was "more
severe than reported" is now correct by a data-derived margin — too severe by a
factor of 30 to 60. `notes/POSTDICTION_PROTOCOL.md` §5 previously withdrew that
verdict on the strength of the misreading; that withdrawal is itself withdrawn.

### The required regime

The claim counts ONLY under a physiological growth law: calibrated hyperbolic
(D015) or linear arrest. **A match that appears only under constant dilution
(`k_mu = inf`) is a FAILURE and will be reported as one**, per rule 4, because
that regime predicts zero reproductive loss by construction and so contradicts
the very number being explained.

### The model form, fixed in advance

Three states `(u, a_r, a_s)`. Aggregate splits into REACTIVE `a_r` and
SEQUESTERED `a_s`:

    da_s/dt = k_seq . a_r^q  -  k_rel . a_s        ( - mu.a_s under dilution )

with `q = 1` tested first and `q > 1` second. `a_s` does NOT appear in the
chaperone or protease resource denominators and does NOT nucleate. That is the
entire mechanism: sequestration changes KINETICS, not bookkeeping. `k_seq = 0`
must reduce to the two-state model exactly.

D024 (the n-state generalisation) applies directly, so the boundary condition
needs no re-derivation — but it must be VERIFIED that it applies, by checking
`det J = -det[grad R; grad G; grad C]` on the extended system, BEFORE any number
from this model is interpreted.

### What would falsify

- **Monostable under the physiological law** -> sequestration does not supply
  the precondition; the mechanism is wrong or insufficient. Third failed
  post-diction.
- **Bistable under the physiological law but the loss falls OUTSIDE
  [0.30, 0.60]** -> a SECOND FAILED POST-DICTION, stated here in advance as a
  failure rather than as a partial success. It would name the next missing
  mechanism, and the direction of the miss identifies which: a loss below 0.30
  says the model under-couples burden to growth; above 0.60 says the sequestered
  state is too costly and something is protecting the cell that the model lacks.
- **Bistable and in band ONLY under constant dilution** -> failure, per rule 4.

### What is NOT claimed either way

Nothing here bears on the fold theorem, which is D007/D024 and holds independent
of how many aggregate states there are. A failure of this post-diction is a
failure of a MECHANISM, not of the identity.

## D029 — The sequestration post-diction FAILS. The model overpredicts the cost of an aggregate-laden cell by 27x to 92x

Run of `scripts/phase3/sequestration.py` against the preregistered criterion in
D028. Reported failure-first, per `notes/POSTDICTION_PROTOCOL.md`.

### Gate first: D024 holds on the extended system

Checked before any other number was looked at, as D028 required. Over 144
combinations of sequestration setting, growth law, state and growth-cost
convention: `det J = -det[grad R; grad G; grad C]` with median relative error
1.5e-12 and maximum 4.7e-11, none above 1e-6. `k_seq = 0` reproduces the
two-state field to 0.00e+00 in every component, diluted and undiluted. The
extension is legitimate and nothing downstream is void.

### The verdict: FAIL, on both criteria, in both regimes

|                          | `a_s` costs growth | `a_s` free |
|---|---|---|
| qualified cells          | 384 | 384 |
| settled                  | 128 | 249 |
| bistable                 | 4   | 15  |
| bistable with `k_seq>0`  | 1   | 12  |
| bistable in the control  | 3   | 3   |
| **in band**              | **0** | **0** |

Every bistable cell's predicted reproductive loss is 0.482, 0.508, 0.954 or
1.000, against a measured band of [0.0104, 0.0178]. The miss is **27x to 92x**,
in the direction of excess severity, and it is uniform — there is no marginal
case.

### Two things the run established that are not the verdict

**1. Sequestration does supply the precondition, in one regime.** With `a_s`
exempt from the growth cost, bistable cells with the mechanism ON outnumber the
control 12 to 3. So the mechanism is not inert: draining aggregate into an inert
compartment can create a second attractor under linear arrest. It was the right
kind of idea. The resulting high state is simply always catastrophic — every one
of those 12 cells sits at loss 1.000, complete arrest.

**2. Under the ORIGINAL band, the literal criterion would have PASSED, carried
entirely by the control arm.** Two cells scored in [0.30, 0.60] and both had
`k_seq = 0` — the two-state model, mechanism switched off, reproducing D026's
numbers exactly. D028 fixed the band, the growth law and the falsifier but never
required the mechanism under test to be ON, so a control cell could have carried
a "pass". That omission is now protocol rule 6, and `verdict()` reports
`mechanism_passes` alongside D028's literal `passes` rather than quietly
replacing it. The band correction made the point moot; the defect was real
regardless.

### What the failure names, quantitatively

D028 stated in advance that a miss above the band means "the sequestered state is
too costly and something is protecting the cell that the model lacks". That is the
direction, and the size can be computed. Inverting the hyperbolic law for the
burden that WOULD produce a measured-size loss:

| `k_mu` | `B_high` actual | `B_high` needed | ratio |
|---|---|---|---|
| 1.0 | 35.55 | 0.140-0.149 | **239x-254x** |
| 2.0 | 3.596 | 0.399-0.417 | 8.6x-9.0x |
| 2.0 | 4.146 | 0.531-0.550 | 7.5x-7.8x |

**The model's bistable high state is 7.5x to 254x more aggregate-laden than the
cell Lindner actually measured.**

**(!) These severity figures UNDERSTATE the miss, and D031 established why.** For
the twelve mechanism-on cells the reported loss of `1.000` is not a measurement
but the LINEAR-ARREST LAW CLAMPING: `mu = mu0.max(0, 1 - (u+a)/k_mu)` returns
exactly zero, and those states sit 3.4x to 43x past arrest. The true severity
there is unbounded. Quote "27x to 92x" only as a lower bound on the miss. A real old-pole E. coli carries a visible
inclusion body at a cost of 1.2-1.8% of growth rate. Every high state this model
can produce is a cell that has essentially stopped dividing.

That is a sharper and more useful statement than "the post-diction failed": the
coupling from aggregate burden to growth rate, or the burden at which the model
places its second attractor, is wrong by one to two orders of magnitude.

### The next mechanism, and its provenance

The observation points at it, not the model. Lindner's cell holds a SMALL,
STABLE, spatially localised deposit and keeps dividing at 96-99% of normal. This
model's high attractor exists only where aggregate has run away far enough for
saturating removal plus dilution to catch it, which is intrinsically a
large-burden state. Sequestration as specified in D028 does not fix that, because
`k_seq . a_r^q` is unbounded: it moves aggregate to a different compartment
without limiting how much there is.

The candidate is therefore a **size-limited deposit** — saturating sequestration
into a finite number of foci, `k_seq . a_r/(K_seq + a_r)`, rather than an
unbounded sink. A deposit that cannot exceed a fixed size gives a high state
whose burden is set by the deposit capacity rather than by where the removal
curve bends, which is the only way this model can place a second attractor at 1%
growth cost instead of 100%.

Stated as the next candidate, not as a result. It has not been run.

### Count

Three collisions with published observation so far, three informative failures:
regulation (D025), sequestration-as-reservoir (D026 named it, D029 tested it),
and now the magnitude of the high state. None has confirmed the theory. The fold
theorem (D007, D024) is untouched by any of them — a mechanism failing is not the
identity failing.

## D030 — PREREGISTRATION: the framing test

Written before `scripts/phase3/asymmetric_division.py` was run.

**A monostable two-state diluted model under the calibrated hyperbolic law
(D015), with aggregate partitioned asymmetrically at division in ratio
`(f, 1-f)` with `f > 0.5` and no second attractor anywhere, scores
`Delta(GR_old - GR_new)/GR_mean` inside the FIXED band [0.0104, 0.0178] of D028,
with the mechanism required ON — `f = 0.5` control cells cannot carry the pass
(protocol rule 6).**

The band is imported from `postdiction_aging.aggregateAttributableLoss()` and is
not re-derived. Both outcomes are results and are reported as such:

- **IN BAND** -> the bistability requirement was a category error. D026's
  surviving claim, that rejuvenation is only coherent in a bistable system, is
  WITHDRAWN: the old-pole cell is not sitting in a second basin, it inherits a
  physical inclusion body at every division, which is a continuously renewed
  perturbation and not an attractor. The three prior failures were then scoring a
  quantity the observation never required.
- **NOT IN BAND** -> no route tried produces a mildly-burdened stable state, and
  that is a structural limitation of the model class, stated as a finding rather
  than as a fourth failed post-diction.

## D031 — The framing test PASSES. Bistability was never required, and D026's surviving claim is withdrawn

Run of `scripts/phase3/asymmetric_division.py` against D030, which was committed
before the run.

### Result: IN BAND

| | |
|---|---|
| cells | 728 |
| settled | 709 |
| **in band, mechanism ON (`f > 0.5`)** | **43** |
| in band, control (`f = 0.5`) | **0** |
| `f` values scoring in band | 0.60 to 0.99 |

The control is exact rather than approximate: across all 66 `f = 0.5` cells the
aging effect is **0.0** with standard deviation **0.0**. Symmetric partitioning
in half the volume leaves concentration unchanged, which is precisely what the
`-mu.x` dilution term already encodes, so the control is the plain diluted model
and returns identically zero. Rule 6 is satisfied by construction and not by a
margin.

The effect is **monotone in `f`** in every setting, rising smoothly from zero. It
is a continuous one-parameter family crossing a narrow band, not a sliver: the
band spans a factor of 1.7, so a monotone curve crosses it in one or two rungs of
an eleven-point ladder, which is geometry rather than fragility. At `mu0 = 0.1`,
**8 of 18** `(p_qc, j)` settings contain an in-band `f`; in the rest the whole
ladder either under- or over-shoots.

### What this means

**The bistability requirement was a category error, and D026's surviving claim is
WITHDRAWN.** That claim was: "rejuvenation is only a coherent category in a
bistable system, since in a monostable one a daughter with less aggregate relaxes
straight back to the single attractor." It assumed the old-pole cell is SITTING in
a second basin. It is not. It inherits a physical inclusion body at every
division. That is a **continuously renewed perturbation, not an attractor**, and a
monostable system under a renewed perturbation has a stationary offset with no
separatrix anywhere.

A MONOSTABLE two-state diluted model, with no sequestration, no second attractor
and no additional state variable, reproduces the measured lineage difference. The
only addition is that division partitions aggregate asymmetrically.

**Three prior failures were scoring a quantity the observation never required.**
D026, and D029's two regimes, all searched for a second attractor whose burden
would have to be 7.5x to 254x higher than the measurement permits. The
measurement never called for one. The old-pole cell is mildly burdened because it
keeps receiving a little more than its share, not because it has fallen into a
basin.

This does not rescue the model class on the point D029 established. It relocates
it: what the model cannot do is place a STABLE ATTRACTOR at a 1% growth cost. It
was never asked to.

### The audit flag: settled

Twelve bistable sequestration cells produced four distinct loss values. The draws
are **independent** — 12 distinct `(mu0, k_mu, k_seq, q, j)` tuples for 12 cells,
nothing coarse-grained in the sampling. The identical `1.000` is the LINEAR-ARREST
LAW SATURATING: `mu = mu0.max(0, 1 - (u+a)/k_mu)` returns exactly `0.0` whenever
burden exceeds `k_mu`, and every one of those high states sits **3.4x to 43x past
arrest** (burden 1.71 to 21.51 against `k_mu = 0.5`). So `loss = 1 - 0/mu_low` is
identically 1.000 in floating point.

The consequence corrects D029 in the conservative direction: for those twelve
cells `1.000` is a **clamp, not a measurement**, and the true severity is
unbounded. D029's "27x to 92x too severe" therefore UNDERSTATES the miss for the
mechanism-on cells rather than overstating it. The four distinct values decompose
cleanly as three hyperbolic control cells plus one saturated value shared by all
twelve.

## D032 — The framing test made parameter-free. It is well-posed and currently unmeasurable

D031 is an accommodation, not a post-diction, and the entry stands corrected. A
monotone one-parameter family rising from exactly zero crosses any band below its
maximum, so "43 cells in band, `f` from 0.60 to 0.99" shows the curve is tall
enough, not that the model predicts 1.0 to 1.8 percent. Worse, there were FOUR
free knobs, not the two an initial reading suggested: `f`, `mu0`, `p_qc` and
`j/j_crit`. Four parameters, one target number, zero residual degrees of freedom.

### Two of the four cancel, analytically

At small burden the hyperbolic law gives `aging = (B_old - B_new)/k_mu`, and since
`k_mu = ARREST_PROTEOME_FRACTION/p_qc` while a model burden of 1 equals a proteome
fraction of `p_qc`,

    (B_old - B_new)/k_mu  ==  32 . (B_old - B_new) as a proteome fraction

**identically** — verified in code, `np.allclose` exact. `p_qc` cancels and `mu0`
never appears. Both enter only through what the stationary aggregate load IS, so
taking that load from measurement removes them. `j/j_crit` goes the same way, for
the same reason. The prediction collapses onto ONE measurable number.

### The fourth is pinned by data, and not where a guess would have put it

**(!) SUPERSEDED BY D033.** The argument below sets `f = 1` by treating the
visible focus as the whole of `a`. It is not: `a` is TOTAL aggregate and the focus
is one object inside it. The correct statement is `f_eff = 0.5 + 0.5.beta` for
focus share `beta`, and the requirement below is the `beta = 1` corner of a
`beta`-indexed interval. The prediction is NOT parameter-free.

Lindner's full text, via PubMed: "52.3% have no inclusion body, 46.5% of the cells
contain only one inclusion body, and only 1.2% carry two inclusion bodies
immediately after cell division", and the new-pole progeny are "devoid of parental
inclusion bodies". The inclusion body is a SINGLE INDIVISIBLE OBJECT. It goes
entirely to one daughter. So `f = 1`, not 0.6 and not 0.95, and the continuous-`f`
model is a mean-field stand-in for an all-or-nothing partition. That is a
modelling mismatch and is stated rather than smoothed over.

### What the model then requires

At `f = 1` the daughters start at `2.a_end` and `0`, so `aging = 64 . a_end` as a
proteome fraction in the linear regime, times a damping factor of **0.4386**
measured by the machinery (the daughter's load relaxes during its own generation,
so the time-averaged difference is smaller than the initial one). Inverting the
fixed band:

    old-pole aggregate = 3.69e-04 to 6.34e-04 of the proteome
                       = 0.0369% to 0.0634%

**This is the whole prediction, with no free parameters.** It is falsifiable by a
single measurement.

### And that measurement does not exist at the required precision

Nearest available, via PubMed: Tomoyasu T, Mogk A, Langen H, Goloubinoff P,
Bukau B (2001) Mol Microbiol 40(2):397-413,
doi:10.1046/j.1365-2958.2001.02383.x, PMID 11309122 — "In DeltarpoH mutants ...
5-10% and 20-30% of total protein aggregated at 30 degrees C and 42 degrees C
respectively", while "In rpoH+ cells, DnaK depletion did not lead to protein
aggregation at 30 degrees C".

The wild-type value is UNDETECTED, which is a bound and not a measurement. The
model's requirement sits **79x to 271x below** the chaperone-crippled `rpoH`-null
figure and below the detection threshold of the standard aggregate-isolation
assay. Consistent, therefore, and uninformative.

### Verdict

**The framing test's structural conclusion stands and does not depend on any of
this.** Bistability is not required to explain Lindner: the control at `f = 0.5`
is exactly zero for an algebraic reason, and a monostable model produces a
stationary lineage offset. That is a negative claim and it needs no fit.

**The quantitative match does not stand as a post-diction.** It is an
accommodation, now converted into a well-posed prediction that current data cannot
score. The paper reports it that way: a falsifiable requirement of
**0.037% to 0.063% of the proteome in the old-pole inclusion body**, together with
the fact that the wild-type aggregate fraction has been reported only as
undetectable. Measuring it to a precision of about 0.01% decides the question.

That is the correct place for it — a stated experimental target, not a claimed
success.

## D033 — beta: the focus is not the whole of `a`, so the requirement is an interval, not a number

D031 and D032 set `f = 1` on the grounds that the inclusion body is indivisible.
The premise is right and the inference was not. The model's `a` is TOTAL
aggregate; the IbpA-marked focus Lindner tracks is one object inside that pool,
and the diffuse and small oligomeric species outside it partition roughly evenly.
**D031 assumed the focus was all of `a`, and never said so.**

### The correction is a reparameterisation, not a new mechanism

Let `beta` = (aggregate mass in the visible focus)/(total aggregate). A fraction
`beta` goes entirely to the old-pole daughter; `1 - beta` splits evenly. Each
daughter has half the volume, so the concentration multipliers are

    old : 2.[beta.a + 0.5(1-beta).a] = (1 + beta).a
    new : 2.[         0.5(1-beta).a] = (1 - beta).a

which is exactly the scalar rule with `2.f_eff = 1 + beta`, i.e.

    f_eff = 0.5 + 0.5.beta

**No new dynamics, and that is the point.** `f` was never a free knob to be
pinned; it is DETERMINED by `beta`, a physical quantity nobody measured. Both
reductions are exact and asserted at 0.0: `beta = 1` gives `f_eff = 1` and
reproduces D031; `beta = 0` gives `f_eff = 0.5`, the control, at which the aging
effect is identically zero.

### Damping recomputed per beta

D032 carried a single damping of 0.4386 measured under `f = 1`. Under
two-compartment partitioning the new-pole daughter is no longer aggregate-free, so
its own relaxation differs. Recomputed: **0.346 to 0.355** across
`beta = 0.145` to `1.0` — the beta dependence is weak, under 3%, but the value is
25% below 0.4386 because that figure came from a wider `(mu0, j)` grid. The
requirement inherits that as a systematic and is quoted with it.

### The beta-indexed requirement

`B_old - B_new = (1+beta).a - (1-beta).a = 2.beta.a`, so the requirement is
proportional to `1/beta` and rises without bound as the focus share falls. LOWER
`beta` needs MORE aggregate for the same lineage difference.

| beta | f_eff | damping | required old-pole aggregate (% of proteome) | ratio to rpoH-null 5-10% |
|---|---|---|---|---|
| 1.00 | 1.000 | 0.346 | 0.0467 - 0.0803 | 62x - 214x |
| 0.75 | 0.875 | 0.355 | 0.0607 - 0.1044 | 48x - 165x |
| 0.50 | 0.750 | 0.355 | 0.0911 - 0.1567 | 32x - 110x |
| 0.25 | 0.625 | 0.355 | 0.1824 - 0.3137 | 16x - 55x |
| 0.145 | 0.573 | 0.355 | 0.3146 - 0.5411 | **9.2x - 31.8x** |

### What the literature bounds, which is very little

Searched via PubMed and checked against full text, not abstracts. **No source
bounds `beta` for the unstressed condition Lindner studied.**

- **Winkler J et al. (2010) EMBO J 29(5):910-923,
  doi:10.1038/emboj.2009.412, PMID 20094032** is the right paper and does NOT
  report the share. It reports the two ingredients separately, under HEAT STRESS:
  "1.5-3% of total cytosolic E. coli proteins, corresponding to approximately
  17500-33000 molecules, aggregate on heat stress"; "individual aggregates are
  composed of approximately 2400-16500 protein molecules"; "62% of all
  heat-treated cells (n=200) exhibited two fluorescent foci at both poles".
  Combining these gives `beta` in **0.145 to 1.0** IF two foci are assumed.
  **That combination is our arithmetic, not their measurement, and it is the wrong
  condition** — heat stress, not the unstressed growth Lindner used.

  **(!) And the assumption does not transfer.** The 62% two-foci figure is from
  heat-treated cells. In the unstressed cells Lindner measured, 46.5% carry ONE
  focus and only 1.2% carry two, so the one-focus arithmetic applies and gives
  `beta` in **0.073 to 0.94** — a floor half as large and nearly the whole unit
  interval. The 0.145 row of the table above is therefore an ILLUSTRATIVE
  pessimistic case, not a bound. **No defensible lower bound on `beta` exists**,
  and the requirement is quoted as a family spanning roughly 0.047% to 0.5% of the
  proteome over the plausible range, with the direction stated rather than an
  endpoint claimed.
- **Lindner et al. 2008** quantifies focus fluorescence and focus counts (46.5% of
  cells carry one focus) but not the focus share of total aggregate.
- **Coquel A-S et al. (2013) PLoS Comput Biol 9(4):e1003038,
  doi:10.1371/journal.pcbi.1003038, PMID 23633942** measures aggregate diffusion
  constants and sizes and the nucleoid-crowding mechanism. It does not report a
  mass share, and is therefore NOT cited for one.

### Verdict

**The prediction is not parameter-free, and D032's paragraph overstated it.** It
is a `beta`-indexed interval, and `beta` has NO defensible lower bound. The
requirement spans roughly

    0.047%  (beta = 1, the D031 assumption)
    ~0.5%   (beta ~ 0.15, an illustrative pessimistic case, NOT a bound)

and rises without limit as `beta` falls further. The direction is the useful
part: **every departure from `beta = 1` moves the requirement toward
detectability**, narrowing the gap to the only measured aggregate load from
62x-214x down to roughly 10x. So the correction improves the prediction's
experimental reach rather than damaging it. No endpoint is claimed.

**The measurement that closes this** is the focus share of total aggregate in
unstressed E. coli — quantitative fluorescence of the focus against total
aggregated protein from a sedimentation assay in the same cells. With `beta`
pinned, the requirement becomes a single interval and the prediction becomes a
test.

The structural conclusion of D031 is untouched: it depends on the control being
exactly zero, which holds at `beta = 0` for the same algebraic reason.

## D034 — Figure 1, and a structural finding it forced: this model class has no infinitesimal-homeostasis point

Panel (b) of Figure 1 was specified to mark a VERTICAL tangent (the saddle-node)
and a HORIZONTAL tangent (infinitesimal homeostasis). Differentiating the
equilibrium condition gives

    da*/dj = G_u / det J

so the vertical tangent is `det J = 0`, which is the theorem, and the horizontal
tangent is `G_u = 0`. **The horizontal tangent does not exist**, and the reason is
structural rather than a property of the base parameters.

`G = v_nuc + v_grow - v_dis - v_degA`. Nucleation and aggregate growth rise with
`u` directly. Raising `u` also sequesters chaperone and protease, which LOWERS
disaggregation and aggregate clearance — and both enter `G` with a minus sign, so
those contributions are positive too. All four terms push the same way, so
`G_u > 0` and `da*/dj` never vanishes.

Measured: `G_u > 0` at all 305 points of the traced nullcline at the base point
(minimum 2.08e-03), and across 40 draws from the Phase 1 box **0 draws** contain
any point with `G_u < 0` (smallest value seen anywhere 6.67e-12). This is a
statement about the model class, not about one parameter point, and §3.1 should
say so rather than positioning the result as "homeostasis sits elsewhere on the
branch".

**(!) AMENDED: panel (b) needs no schematic.** The reasoning above stopped one
derivative short. Cramer on the OTHER column gives `du*/dj = -G_a/det J`, so the
SOLUBLE coordinate has a horizontal tangent where `G_a = 0`. That point exists, is
unique, and is generic: `det J = R_a.G_u = 2.027e-03` there, not zero. It is the
nullcline's own turning point, and it is the same locus at which root-finding in
`a` at fixed `u` loses the curve -- the numerical difficulty and the horizontal
tangent are one thing. Panel (b) therefore plots `u*(j)` and `a*(j)` together,
with the turn at `j_turn = 0.15409` and the fold at `j_crit = 0.15424`, all
computed, and the inset is a x182 ZOOM rather than a drawing.

**It is not a second prediction, and is not offered as one.** `j_turn/j_crit` is
0.9990 at the base point. Across 30 draws the locator placed the turn on the
stable branch in only 6 cases, with median ratio 0.9963, min 0.3513 and max 1.4032
(above `j_crit`, hence off the accessible branch). It is branch geometry, not
something an experiment could separate from the boundary.

### Two things the figure work corrected

**Root-finding loses the curve at its turn.** Tracing `{G = 0}` by finding roots
in `a` at each fixed `u` returns two disconnected pieces with the fold sitting in
the GAP between them, because the two roots merge at the turning point. The figure
traces the nullcline as a contour instead, which follows it through the turn. The
first rendering of Figure 1 showed the fold floating off its own nullcline.

**`rel_err` must not appear in a caption.** `determinantIdentity`'s `rel_err`
divides by `max(|det J|, |cross|)`, and BOTH vanish at a saddle-node, so at the
plotted fold it is 0/0 and returns exactly 1.0 no matter how well the identity
holds. The caption quotes `sin(grad R, grad G) = 3.5e-10` instead, which is
scale-stable there. A test asserts the caption does not contain the phrase
"relative error".

### What is deliberately not figured

Recorded so the omissions are choices rather than oversights:

- **The §8.2 empirical failures.** The 7.5x-254x miss is a table. Four failed
  post-dictions rendered as a figure would give them visual weight the argument
  does not want.
- **The square-root capacity ceiling (D027 §5).** Exact, but never binding in the
  tested range, and least numerically checkable in the regime where it would bind.
  A plot would overstate it.


## D035 — Section 5, defect one: the verification statistic degraded toward the thing it verified

`determinantIdentity`'s `rel_err` divides by `max(|det J|, |cross|)`, and BOTH
vanish at a saddle-node. Over the whole load grid of 325 folds:

| metric | median | p90 | max | corr with \|eig\| |
|---|---|---|---|---|
| max-normalised (what §5 reported) | 2.00e-07 | 7.81e-07 | **1.55e-02** | **-0.262** |
| gradient-normalised (D027) | **2.34e-10** | 6.16e-10 | **1.29e-09** | +0.060 |

The negative correlation is the defect stated as a measurement: the error GROWS as
the bracket tightens on the true fold. Splitting at the median eigenvalue, the
tighter half has median 2.79e-07 against 1.38e-07 for the looser half, twice as
bad when closer to the object being verified. The `1.55e-02` tail is what a
referee would find in a paper whose selling point is exactness.

**Correction: change the statistic.** §5 now reports the gradient-normalised
residual and Fig. S1 carries the contrast rather than asserting it.

Found by the caption audit for Figure 1, where the metric returns exactly 1.0 at
an exact fold, which is the one place a caption would quote it.

## D036 — Section 5, defect two: three headline numbers came from a 6% subsample, and all three were optimistic

Separate from D035 and with a different correction. `verifyAgainstRun` takes
`n_identity = 20` states from the load grid's 325, and §5 quoted those values
without saying they were a subsample. Recomputed over ALL 325:

| quantity | 20-state subsample | all 325 | direction |
|---|---|---|---|
| identity residual, median | 1.436e-07 | 2.00e-07 | worse |
| corr(log sin angle, log \|eig\|) | +0.9987 | **+0.9960** | worse |
| \|G\| at fold states, maximum | 8.2e-10 | **1.63e-09** | **2x worse** |

**Every one moved in the flattering direction.** That is what a subsample does
when it is not declared, and it is a separate failure from the normalisation
choice: this one is not fixed by changing the metric, only by not subsampling.

**Correction: use the whole population, and name it.** §5 previously described its
box as "2884 fold states" while the identity rows were computed on experiment B,
a different population entirely. The two are now distinguished at every row:

- **load grid** — 325 folds, nascent occupancy against rescue allocation at fixed
  kinetics. The identity, parallelism, `|G|` and solver rows.
- **kinetic box** — 5000 Latin-hypercube draws, 2884 admitting a fold. The `phi`
  rebuild and the saturation-fraction distributions.

A referee comparing "2884 folds" in the header to a table computed over 325 would
have assumed the worse thing.

### Forensics: nothing selected those 20, and only one of the three was biased

"Three independent quantities, all flattering" looks like a selection process. It
is not. The draw is `b.sample(n=20, random_state=1)` -- uniform, seeded, with no
filter on bracket quality or convergence. Over 4000 redraws of 20 from the 325:

| statistic | P(subsample flatters) | reading |
|---|---|---|
| median residual | 0.501 | coin flip |
| correlation | 0.438 | near chance |
| **maximum** | **0.939** | **deterministic** |

`1 - 20/325 = 0.938`. A maximum over a subset can never exceed the maximum over
the set, so the `|G|` bound was guaranteed to understate; the other two are
ordinary sampling noise, and both landing flattering has probability about 0.22.

So the bias IS bounded and IS explained: one structural underestimate plus two
coin flips, not an unknown selection. The general rule this yields is sharper
than "do not subsample": **never subsample an extremum.** Medians and
correlations on a subsample are noisy but unbiased; maxima and minima are biased
by construction.

**Neither defect was reachable by the test suite.** A test asserting that
1.436e-07 is reproduced passes whether or not the metric degrades and whether or
not the sample is complete. See `notes/VERIFICATION_RULE.md`.

### Every other site carrying these values is corrected or marked

Grepping the three numbers across tracked text found them in `STATUS.md`,
`theory/FOLD_THEOREM.md`, `manuscript/COLLAPSE_BOUNDARY.md` and in this log.
`STATUS.md` is corrected in place. `theory/FOLD_THEOREM.md` carries a superseded
banner above the affected table rather than a rewrite, since the numbers there
were correct for what was computed at the time. `COLLAPSE_BOUNDARY.md` already
carries a whole-document superseded banner. Entries in this log are left as the
historical record they are. A corrected value in one document and the original in
another is a lineage split inside one repository, which is the failure this
project has already paid for once.


## D037 — The extremum rule is an audit criterion, and applying it found a fifth affected number

D036's rule was recorded and then under-applied. "Never subsample an extremum" is
not a note about one table row: the paper reports maxima and worst-cases in at
least eight places, each of which needs its population named and its completeness
checked.

Auditing them found that `direct solver against continuation sweep, maximum
relative error` was ALSO a 20-state draw (`random_state=3`). Full load grid:
**7.56e-07**, against the 6.652e-07 previously reported. Same mechanism, same
direction, and it survived four earlier corrections because the rule had been
written down but not swept.

### Every §5 extremum, with its population and a p99

| quantity | population | median | p99 | max |
|---|---|---|---|---|
| identity residual | load grid, 325, complete | 2.34e-10 | 9.67e-10 | 1.29e-09 |
| `\|G\|` at fold states | load grid, 325, complete | 4.95e-14 | 9.40e-10 | 1.63e-09 |
| solver vs continuation | load grid, 325, complete | 3.03e-07 | 7.20e-07 | 7.56e-07 |
| `phi` rebuild | kinetic box, 2884, complete | 1.3e-13 | 2.98e-09 | 7.25e-09 |

### Why the p99 and not only the max

**A maximum grows with the size of the population it is drawn from.** The `phi`
maximum of 7.25e-09 over 2884 draws and the identity maximum of 1.29e-09 over 325
are not comparable as stated, and neither is comparable to whatever a reader gets
on a rerun of a different size. A maximum is therefore a weak verification
statistic even when computed correctly. A p99 is stable under resampling and
carries the same bounding claim, so every extremum doing bounding work is now
reported beside one.

### theory/FOLD_THEOREM.md is LIVE, and is corrected in place

The README describes `theory/` as where the theorem is stated, so that document is
current rather than historical. A banner above a live table is exactly how a future
session reads the table and not the banner, which is the mechanism behind the P3
caption divergence. The table is therefore corrected in place, with a note
recording what the earlier values were and why they moved.

### Pin the property, never the token

Two tests asserting `1.436e-07` was ABSENT from the manuscript failed once §5
began stating its own corrections explicitly. A string-absence test penalises the
correction discipline the project runs on. They now assert the value never appears
as a CURRENT claim -- absent from the results table, present only inside the
superseding sentence. Recorded in `notes/VERIFICATION_RULE.md`.


## D038 — The extremum sweep is exhaustive, not partial

D037 swept §5. Five further extrema sit outside it and are now checked against the
same criterion, so the sweep can be described as exhaustive rather than partial.

| quantity | population | complete? | action |
|---|---|---|---|
| n-state identity max 4.7e-11 | 144-point grid | yes, full enumeration | population named |
| D027 worst medians 6.4e-14 / 4.6e-13 | 20 networks from the kinetic box | subsample | named; a MEDIAN is unbiased under subsampling |
| `j_crit/j_max` max 0.623 | 8 drawn networks | subsample | relabelled a **lower bound** |
| `sin(grad R, grad G)` "below 2.0e-9 throughout" | folds solved in that sweep | yes | "throughout" replaced by an explicit completeness statement |
| per-network spread 13.6x | 10 draws | subsample | relabelled largest-observed |

Two were genuine repeats of D036: a maximum over 8 draws and one over 10, both
reported as maxima. Neither changes a conclusion, and each now says so where it
appears -- the square-root ceiling is not binding at any value below 1, and a
larger sample can only widen an overlap. The medians at `eps = 100` needed only
their population named, since subsampling makes a median noisy but not biased.

"below X throughout" is a maximum wearing a completeness claim, and is the form
most likely to survive an audit unnoticed. It is now phrased as what it is.

### Section 5's apparatus stops here

§5 now carries a four-row table with median, p99 and max per row, prose on why
maxima are incomparable across populations, and one paragraph recording that
corrections were made. Every element earns its place: the p99 is a methodological
improvement, the populations prevent a specific misreading, and the corrections
paragraph buys credibility that §8 will spend. **No further defensive apparatus is
to be added to §5 unless a new defect forces it**, and the corrections paragraph
does not grow with each find. A sixth affected number is corrected in place and
covered by the existing paragraph. Past a point, self-correction reads as anxiety
rather than rigour, and §5 is close to that point.


## D039 — The saturation screen is not applied, because the data does not support one

Figure 2 was specified to screen draws collapsing at `s_a` near 0.003, per a
standing limitation. Building it showed the screen is not defensible.

A screen requires the low-`s_a` draws to form a distinct cluster. On a log axis the
distribution runs smoothly from 1e-5 to 1 with no gap, and the median of the
survivors slides continuously with the floor:

| floor | 1e-4 | 5e-4 | 1e-3 | 2e-3 | 3e-3 | 5e-3 | 1e-2 | 2e-2 |
|---|---|---|---|---|---|---|---|---|
| median `s_a` | 0.090 | 0.124 | 0.162 | 0.194 | 0.227 | 0.267 | 0.309 | 0.355 |

Any floor is a free parameter moving a load-bearing number by a factor of four.
The first version of the figure used 5e-3 and reported a median `s_a` of 0.267
against the 0.056 in §6 -- a 5x divergence between figure and text, created by the
screen alone.

**No screen, and no exclusion at all.** Dropping even the 47 numerically-zero
draws moves the medians to 0.1858 / 0.1602 / 0.0603 against 0.175 / 0.155 / 0.056
in the text. Introducing a second population to fix a cosmetic issue is how the
325-against-2884 confusion arose. The complete population reproduces §6 exactly,
and the figure carries an inset showing the sensitivity so the decision is visible
rather than asserted.

### A sixth affected number, corrected in place

Per D038, corrected without extending §5's list. `s_u` p5-p95 width 0.890 was
quoted in §9 as a general property; it is the REGULATION experiment's value over
its own 30 networks. The kinetic box's 2884 give **0.876**. Both sites now name
their population. Nothing material changes.
