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

| metric | median | p90 | max | corr with \|eig\|, log-log | raw Pearson |
|---|---|---|---|---|---|
| max-normalised (what §5 reported) | 2.00e-07 | 7.81e-07 | **1.54e-02** | **-0.835** | -0.221 |
| gradient-normalised (D027) | **2.34e-10** | 6.16e-10 | **1.29e-09** | +0.041 | +0.073 |

**Four entries corrected by D041.** This table first read `1.55e-02`, `-0.262` and
`+0.060`, and gave the median split below as 2.79e-07 against 1.38e-07. Those four
were computed in-session and never by committed code; none reproduces. The values
above come from `scripts/figures/fig_identity.py:captionNumbers`, over the same complete
325. Every correction moves in the same direction, which is that the defect is
LARGER than this entry claimed when it diagnosed it.

The negative correlation is the defect stated as a measurement: the error GROWS as
the bracket tightens on the true fold. Splitting at the median eigenvalue, the
tighter half has median 3.13e-07 against 1.46e-07 for the looser half, 2.15x worse
when closer to the object being verified. The `1.54e-02` tail is what a referee
would find in a paper whose selling point is exactness.

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


## D040 — A limitation is an observation; the screen was a pre-authorised free parameter

D039 removed Figure 2's screen. The screen was not a figure-authoring mistake. It
came from a standing limitation in `theory/FOLD_THEOREM.md`, written months
earlier, which read "these should be screened rather than averaged into a median."
The defect was latent in the repository, recorded as a limitation, and would have
propagated into the paper the first time anyone acted on it — which is what
happened.

The observation and the instruction are separable, and only the first survives.
Those draws are marginal; the exclusion does not follow. Both claims are now
stated separately at the site, with the tested reason the screen fails.

### Sweep for other limitations phrased as directions

`theory/`, `notes/` and `STATUS.md`, all limitation, caveat and open-question
text. Classification: observation (safe), restriction (safe — adds a requirement
rather than removing data), direction (a free parameter with prior
authorisation).

| site | form | verdict |
|---|---|---|
| `FOLD_THEOREM` — "these should be screened rather than averaged" | direction | **rewritten** |
| `FOLD_THEOREM` — "should not be quoted as a result without a larger sample" | restriction | kept |
| `FOLD_THEOREM` — "should not be quoted without stating which law produced them" | restriction | kept |
| `FOLD_THEOREM` — "every `phi` ... should be read as `phi_enz` at `delta = 0`" | direction in grammar | kept — an algebraic identity, changes no number |
| `FOLD_THEOREM` — "claims should be pitched accordingly" | direction about rhetoric | kept — touches no data |
| `PREDICTIONS.md` — six "should" statements | hypotheses about the world | kept |
| `DYNAMICAL_SYSTEM` / `MATCHED_IMPLEMENTATION_PROTOCOL` "must" statements | restrictions | kept |

One in the dangerous form, and it is the one already found. Recorded in
`notes/VERIFICATION_RULE.md` as an audit form so the next such bullet is caught
when written rather than when acted on.

### D038's relabelling had reached the manuscript only

Found by this sweep, and it is the lineage-split failure `VERIFICATION_RULE.md`
warns about, occurring inside the correction that introduced the warning. D038
relabelled `j_crit/j_max` max 0.623 as a lower bound and the 13.6x per-network
spread as largest-observed. Both were corrected in `manuscript/bmb_v4.md`. Neither
was corrected in `STATUS.md` or `theory/FOLD_THEOREM.md`, which still read "up to
13.6x", "as much as 13.6x" and "max 0.623" — three sites carrying the superseded
form against a corrected manuscript. Now corrected at all sites.

The 13.6x position is stated rather than left to be reconstructed under review:
**the design was not sized to bound the per-network spread, and we cannot bound
it.** That is a weak claim and it is the correct one. The overlap with the 8.86x
between-network spread survives regardless, because more sampling only widens an
overlap.

### The sixth number was a different failure from the other five

Recorded because the checklist did not distinguish them. Five were incomplete
populations. The sixth — `s_u` p5–p95 0.890 attributed to the kinetic box when it
belongs to the regulation experiment's 30 networks — is a complete population
under the wrong sentence. "Is the population complete?" cannot catch it, because
it is. `VERIFICATION_RULE.md` now asks both questions separately.


## D041 — Section 5's normalisation contrast had no generator, and does not reproduce

Found while building Figure S1, by the ordinary discipline of reconciling a
figure's numbers against the text before writing its caption. §5 and D035 report
that the max-normalised residual correlates with the recorded eigenvalue at
**-0.262** and the gradient-normalised one at **+0.060**. Neither reproduces.

Checked exhaustively on the complete 325-state load grid: two residual forms x
{raw, log10} x four proximity variables x {Pearson, Spearman, Kendall} — 48
combinations — and separately on the 20-state draw that caused D036, recomputed
from the run root through `fold_theorem.determinantIdentity` itself. The published
pair appears in none of them.

| statistic | published | log-log | raw Pearson |
|---|---|---|---|
| max-normalised vs \|eig\| | -0.262 | **-0.835** | -0.221 |
| gradient-normalised vs \|eig\| | +0.060 | **+0.041** | +0.073 |

The medians and p90s in D035's table DO reproduce exactly, which is what made this
hard to see: a wrong population would have moved every entry, and only some moved.
The correlations, the max (1.55e-02 against a true 1.5448e-02) and the median split
(2.79e-07 / 1.38e-07 against 3.13e-07 / 1.46e-07) came from a session computation
that was never committed. This is precisely the failure `CLAUDE.md` names — a
number whose lineage runs to a summary rather than to code — occurring inside the
entry that diagnosed the metric defect.

**The claim survives and strengthens.** The point of the number is that the
max-normalised statistic degrades toward the object it verifies. At -0.835 it
degrades far more decisively than at -0.262. Every correction here moves against
the earlier reporting, so nothing was flattered in the direction of the argument.

**The correction is structural, not a retyped number.** `scripts/figures/fig_identity.py`
now owns every quantity §5's normalisation paragraph and Fig. S1's caption use, in
one function, `captionNumbers`, which prints them all; a test asserts the
manuscript against that function rather than against a stored value. The
convention is stated in the paper: log-log throughout, matching the +0.9960
statistic, with the raw Pearson given beside it so a reader who correlates the
columns directly meets their own number instead of a discrepancy.

### A seventh affected number, and a form the extremum sweeps could not see

§5 read "the single state bracketed to eigenvalue -4.2e-9 has sin(angle) =
3.8e-8". That state is the **minimum of the same 20-state draw**. Over the full
325 the tightest bracket is |eig| = 8.10e-10 with sin = 7.75e-09, and three states
are tighter than 4.2e-9 — so "the single state" was also false as a description.

D037 and D038 swept every maximum, minimum, worst-case and "throughout" and both
missed this, because it is an extremum wearing a **definite description**. "The
single state where...", "the one draw that...", "the tightest/closest/last..." are
extrema in prose form and carry the same `1 - k/N` bias. The audit criterion is
extended to cover them in `notes/VERIFICATION_RULE.md`.

The direction is the usual one: a subsample understates tightness, so the true
tightest bracket is 5x tighter and its parallelism error 5x smaller than reported.
The corrected number is better for the claim than the published one.

### Two smaller repairs at the same site

- §5 said "Three of these values were previously reported from a 20-state random
  subsample" and then listed **four**. D037 added the solver maximum to the list
  without updating the count. Now five, listing five.
- Fig. S1 supplies a fact the text did not have: the log-log slope of sin(angle)
  against |eig| is **1.00**. The parallelism error is *proportional* to bracket
  looseness, not merely increasing in it, which is a sharper statement of "bracket
  tolerance" than the correlation alone. Added to §5.

### The generator that caused D036 was still armed

`verifyAgainstRun(n_identity: int = 20)` still defaulted to a 20-of-325 draw. §5's
numbers had been re-sourced elsewhere, so nothing was wrong in the paper — but the
function whose printout a future session would read still silently subsampled, and
still led with the retracted max-normalised residual. Correcting the numbers while
leaving the mechanism armed is fixing the symptom.

Default is now `None`, meaning the complete population; an integer still
subsamples, deliberately, and the printout says `(SUBSAMPLE)` when it does. The
gradient-normalised residual is the reported statistic and the max-normalised one
is printed beneath it labelled RETRACTED, so the retraction stays visible rather
than being quietly deleted. No test depended on the old default.

One run now reproduces §5 end to end: median 2.340e-10, max 1.292e-09,
corr +0.9960, solver max 7.564e-07, phi over 2884. Previously no single script
printed those five.


## D042 — The figures are wired into the prose, and two divergences surfaced doing it

Five main figures and one supplementary now exist. Attaching them to the text is
where the remaining defects were, because a figure and a paragraph reporting the
same quantity are two lineages for one number.

### Numbering follows first mention, so the new figure is Figure 2

The dilution decomposition supports Corollary 3 in §4.2, which precedes §6 and §7.
By order of first mention it is therefore **Figure 2**, and the saturation, front
and beta figures shift to 3, 4 and 5. It was specified as "Figure 5 (optional
sixth)" on build order, and build order is not mention order — the same confusion
that had the beta figure built before the two that precede it.

**Filenames are now semantic, not positional.** `fig_dilution.py`, not `fig5.py`;
`saturation.tsv`, not `fig2_saturation.tsv`. Each script declares its number in one
`FIGURE` constant and writes `figures/figN.*` from it. A positional scheme would
have to be renamed again if this figure drops on page count, which is a live
possibility, and the data files carried stale numbers into the sha256 manifest. The
manifest hashes are unchanged by the rename, which confirms only names moved.

### The prose reads the figures rather than pointing at them

Added where the text was silent: §3.1 now states Fig. 1b's two ordered
singularities with values (`j_turn = 0.154090`, `j_crit = 0.154239`,
`det J = R_a·G_u = 2.02749e-03` at the turn) and, in one sentence, that the
root-finder's failure and the horizontal tangent are the same locus rather than two
problems. §6 states the no-screen decision with the sensitivity ladder in the prose
where a reader meets the medians, not only in the caption.

**Duplication resolved per case, recorded here.** §6's table keeps the medians and
the Fig. 3 caption no longer restates them; the caption keeps the p5-p95 widths,
which the table lacks, and §9 now points at the panel instead of restating 0.876.
The regulation widths stay in prose because they belong to an experiment no figure
carries — the distinction the sixth affected number turned on. §5's table and the
S1 caption both carry +0.9960 because the panel renders it; a test asserts both
against `captionNumbers` rather than against each other.

### Divergence one: the beta figure used a single damping, the table used four

`fig_beta.DAMPING = 0.35` against a table built from the per-beta measurement. Every
row of the figure sat ~1.5% off §8.4, with no cause but convenience. The damping is
weakly beta-dependent (0.3462 at beta = 1, 0.3544-0.3552 below) and is now stored at
full precision in `asymmetric_division.DAMPING_BY_BETA` with `dampingAtBeta`
interpolating between measured points. Storing it at 3 decimals was tried first and
lost the agreement with D033 in the fourth digit, which is how the precision
requirement was discovered rather than assumed.

The stored table is a cache — `dampingForBeta` takes about a minute per beta — and
a cached number with no live check is the D041 defect. `tests/phase3/
test_asymmetric_division.py` recomputes it at both ends and asserts the cache.

### Divergence two: a claim quoted at the marked beta, not over the plotted range

Writing the Fig. 5 callout, the draft said the requirement never comes within a
factor of fifteen of the only measured aggregate load. That is true at the three
marked beta and **false at the left edge of the figure's own x-axis**: the range
plotted runs to beta = 0.05, where the gap closes to **3.19x**. Corrected to
"between 3x and 214x below", both computed by the figure script and asserted.

A smaller one at the same site: the closest marked ratio is 15.94x, which the table
correctly rounds to 16x, but the prose had promoted the rounded figure to a claim
of "16x or more". It now says "at least fifteenfold", and a test bounds the value on
both sides so neither the rounding nor the claim can drift.

### What the wiring test pins

`tests/figures/test_figure_wiring.py`: numbering follows first mention; every main
figure is cited in prose BEFORE its embed; embed paths match their labels; each
caption number agrees with the text number for the same quantity and both come from
the generator; and every figure is 84 or 174 mm wide and under 234 mm tall.

### A second lineage split, same shape as D038's

Found by sweeping for stale values after the rename. D037 corrected the solver
maximum from the 20-state draw's **6.652e-07** to **7.564e-07** over all 325, and
D041 corrected the tightest bracket. Both reached the manuscript. Neither reached
`STATUS.md` or `theory/FOLD_THEOREM.md`, which still carried 6.652e-07 in two
places — and `FOLD_THEOREM.md`'s table sits directly under a banner reading "the
values below are full-population", so the document asserted its own correctness
while one row contradicted it. The `phi_enz` band was stale in the same two files
(0.125-0.134 against a computed 0.1245-0.1343).

This is the second occurrence of the identical failure in three decisions. The
manuscript is corrected first because it is what gets read; the supporting
documents are corrected only if someone remembers. **The rule to carry forward: a
correction is not applied until `grep` for the old value across the whole repository
returns nothing outside a supersession note.** D038 said this and it was still not
done, which is why it is now written as a step rather than an observation.


## D043 — The submission PDF is generated from the markdown, not hand-ported

`manuscript/bmb_v4.{tex,pdf}` and `manuscript/bmb_v4_supplementary.{tex,pdf}` are
built by `scripts/manuscript/to_latex.py`. A hand-made `.tex` would be a second
lineage for every number in the paper, which is the defect D041 found in the
opposite direction; the markdown stays the single source and the `.tex` files
carry a do-not-edit banner.

**The supplementary is a separate document.** Journals take it as its own file, and
a supplementary figure inside the main PDF is numbered as part of the paper.
`\thefigure` is redefined to `S\arabic{figure}` there, so the two sequences cannot
collide.

### Every span is math unless declared code

276 backtick spans; exactly TWO are literal code. A regex deciding which is which
would silently mangle the paper, so the rule is inverted: math by default, and a
span that LOOKS code-shaped but is not in `CODE_SPANS` raises. A new code span
fails the build rather than rendering as mathematics.

### Four silent corruptions, all found by rendering rather than by exit codes

1. **`j_crit` in math is `j_c` followed by `rit`.** Multi-character subscripts must
   be braced. Unbraced, it typesets as a different symbol and no error is raised.
   This appears dozens of times in the paper.
2. **An unescaped `%` in a caption centred the second half of the document.**
   Figure 5's caption is the only one quoting percentages. `%` opens a LaTeX
   comment, which swallowed the closing brace of `\caption{...}`; pandoc could not
   find `\end{figure}`, escaped the entire raw block into `\textbackslash
   begin\{figure\}`, and the orphaned `\centering` centred every paragraph from
   section 10 to the end. The build reported "figures 6" throughout. The converter
   now asserts that no raw block was escaped and that the `\includegraphics` count
   matches the figures converted.
3. **`width=\linewidth` threw away the entire width discipline.** The figures are
   built at exactly 84 or 174 mm because that is what the journal specifies and
   every font size in them was chosen for that width. Scaling an 84 mm panel to the
   160 mm text block enlarged its 7 pt labels by 1.9x. Now `max width`, which
   shrinks only what would overflow.
4. **Pandoc will not parse `$...$` as math if a space precedes the closing
   delimiter** — it escapes both dollars instead. The Greek substitution left a
   trailing space, so `μ` rendered as a literal `$\mu$`.

Two smaller ones: the markdown numbers its own headings and LaTeX numbers them
again, giving "0.1 1. Introduction", so the manual number is stripped **after
asserting that LaTeX's counter will produce the same number** — the prose
cross-references those numbers, and a silent shift would break all of them. And
the byline printed twice, because `\maketitle` and the markdown front matter both
supply it.

### What this cost, and the pattern

Every one of these produced a valid PDF. Not one raised an error. They were found
by rendering pages to PNG and reading them, which is the same method that found
the fold sitting off its own nullcline (D034), the `0/0` caption metric (D035) and
Figure 2's screen (D039). **Four of the last seven real defects in this project
were found by looking at output, and none by a passing test.** The tests added
here pin what rendering revealed; they could not have found it first.

`tests/manuscript/test_to_latex.py`, 15 checks.


## D044 — The class swap happened before the audit, and it changed the answer

D043's visual read validated the output of `article`. Three of its four
corruptions were class-interacting, so that read certified an artifact nobody will
submit. `sn-jnl.cls` (Springer Nature, December 2024, LPPL) is vendored at
`manuscript/springer/` and the build now targets it.

**The reference-style option is a trap and the template arms it.** Springer ships
`sn-article.tex` with `sn-mathphys-num` UNCOMMENTED and `sn-mathphys-ay` commented
out one line below. BMB is author-year, and this manuscript's reference list and
every in-text citation are author-year. Taking the template's default would have
silently renumbered the list and left the citations pointing at nothing. The build
declares `sn-mathphys-ay` and a test asserts `-num` never appears.

**The text block is 130.7 mm, not 160.** Measured, not assumed:
`textwidth = 372 pt`. So the 84 mm figures sit at natural size, and **the 174 mm
figures do not fit** and are reduced to 75%. That is inherent to a one-column
submission format and cannot be designed away, but at 75% the smallest text in
those two figures fell to 3.0-3.2 pt.

The fix was not to rescale. Those figures carried 4.0-4.6 pt annotations, which is
**already below Springer's minimum at their own final width** — the class swap
exposed a defect that predated it. The font floor in the two wide figures is now
6 pt, and Figure 3 was re-laid out around the larger text: taller panel, the three
median/width labels aligned at a fixed x clear of the `s = 1` rule, and the
sensitivity inset moved into the empty band below the last strip.

### The equations were verbatim, not mathematics

The markdown writes displays as fenced code blocks. Pandoc renders those as
`verbatim`, so **every displayed equation in the paper was typeset in typewriter
face with ASCII operators**, and the widest ran past the right margin. Ten blocks
now have declared layouts in `DISPLAY_MATH`, keyed by their exact source so that
editing an equation invalidates its key and stops the build. Layout is a judgement;
the coupling is deliberate.

### A caption bug that survived a full visual read

`(u*, a*)` inside a caption became `\emph{, a}`: the markdown-emphasis rule added
in D043 matched from the first star to the second and deleted both starred
equilibria. It was present in the `article` build and **the 21-page read did not
catch it** — which is the honest limit of reading as a method. Backtick spans are
now protected before emphasis conversion.

### Counts before pages

`figures 6` appeared in the log on every run of the broken build while one figure
was mangled. The count was wrong, persistently, and nothing compared it against
anything. **A count in a log that nobody asserts against is a comment, not a
check.** The build now asserts sections, main figures, supplementary figures,
tables, display equations, code spans and stripped sections against `EXPECTED`,
that `figure`/`longtable`/`equation*`/`align*` are balanced, that no `\caption{`
is unclosed, and that page counts fall in a stated range. Most of what reading
caught this round is now caught in milliseconds.

### `### Reference verification` no longer reaches the reader

It records which references were checked against PubMed and what was corrected.
That is repository provenance, not part of the paper. Stripped from both PDFs by
`INTERNAL_SECTIONS`, kept in the markdown.

### Two things decided now rather than under deadline

**The generated-`.tex` discipline has an expiry.** It holds through review, where
the markdown is edited and the `.tex` regenerated. **It breaks at proofs**, where
Springer returns a typeset file and expects corrections against that. The handover
is declared here: at the moment proofs arrive, `manuscript/bmb_v4.tex` stops being
a build artefact and becomes the artifact; `to_latex.py` is retired for that
manuscript; and `bmb_v4.md` is banner-marked superseded rather than left as a
stale second lineage. The rule that a number must be recomputable from a script
still holds — it just points at the proof file from then on.

**Figure 5 was included by default, not by decision.** It was specified as
droppable on page count. The main text is 23 pages, which is comfortable for BMB,
so it stays — but this records that no one chose to keep it. If anything must come
out later, it is the one already marked droppable, and nothing else should be cut
before it.


## D045 — Block 1: the converse needs hypotheses, and three of my own numbers were artefacts first

External review, major revision, core sound. Block 1 is the part that gates the
title, so it ran before any rewriting.

### 1a. The converse, with its hypotheses, verified

Theorem 1 proves saddle-node => grad R parallel to grad G. It does NOT prove the
converse: `det J = 0` is equally consistent with a transcritical or pitchfork
bifurcation, a cusp, a double-zero eigenvalue, or a curve of equilibria. The word
"exactly" in part (3) claimed an iff that was never established.

`scripts/phase3/genericity.py` verifies the standard conditions at solved fold
states: (G1) `|grad G| != 0`, (G2) `tr J != 0`, (G3) `d2R/ds2 != 0` along the
constraint, (G4) parameter transversality `w . dF/dj != 0`.

| population | evaluated | violations |
|---|---|---|
| load grid | 325 of 325 | **0** — minima 0.106 / 0.303 / 0.0929 / 0.341 |
| kinetic box | 2767 of 2884 | **2** — `draw115`, `draw1700` |

(G4) is the one that could have failed on its merits: `dF/dj = (1, 0)` exactly, so
transversality reduces to `w_1 != 0` and nothing in the model forces it. It is
bounded away from zero by 0.34 at worst on the load grid.

The two failures fail (G1), (G3) and (G4) TOGETHER, at `a = 1.2e-11` and
`1.3e-15`. The aggregate pool has collapsed to numerical zero, `grad G -> 0`, and
every condition defined on the constraint degenerates with it. That is the edge of
the model's domain, not a counterexample, and is reported as such. 117 of 2884
networks have no solvable fold state and carry no result.

**`d2R/ds2 < 0` at all 325 load-grid states.** The fold is a constrained LOCAL
maximum of removal. Remark 3 withdrew the *global* maximum claim over the closed
curve and that withdrawal stands, but the local statement is now measured, and
Remark 3 should say which of the two it concedes.

### 1b. Hopf exclusion, load grid complete

Tracing `tr J` along the stable low-burden branch from small `j` to `j_crit`,
contour-traced so it survives the nullcline's turn:

- `max(tr J)` median **-0.338**, largest **-0.243**
- networks reaching `tr J >= 0` anywhere: **0 of 325**
- singular points per network: median 1, max 1, `>1` in **0 of 325**

So on this population the fold is the first loss of stability, and Corollary 1's
"no continuation sweep" needs no branch-identification caveat. The load grid
varies load at fixed kinetics and cannot speak for the network ensemble; the
kinetic box decides the title and is still running.

### Three of my own numbers were method artefacts before they were results

**The 105 Bogdanov-Takens states did not exist.** Testing at RECORDED kinetic-box
states gave 105 of 2884 with `tr J ~ 0`, which reads as 3.6% of networks at a
double-zero eigenvalue. Those states have `det J` up to 3.6 — they are not folds.
`|det J|` at recorded kinetic-box states has median 2.8e-05 but **p99 2.4 and max
52**, against a load-grid maximum of 8.9e-05. At solved states (G2) has zero
violations. Reporting the 105 would have manufactured a serious finding out of
loose bracketing.

**`foldSolve` cost 25% of the ensemble, non-randomly.** It seeds from a coarse
scan over `lowerNullclineA` — a first-root heuristic whose own docstring says it
does not identify the branch the fold lives on — and failed for 718 of 2884.
Seeding Newton from the recorded state, which experiment C's continuation already
placed on the right branch, recovers 2767 and runs ten times faster. "No
violations among the survivors" over a 75% subset selected by solver convergence
would have been a selected result on precisely the geometry under test.

**Half the section 6 shift was population, not correction.** Rebuilding at solved
states moves `s_ref` 0.175->0.180, `s_u` 0.155->0.159, `s_a` 0.056->0.050, `phi`
0.0769->0.0825. Paired on the same 2761 networks the per-network median relative
difference is 2.5e-04, 2.6e-04, 1.4e-03 and **4.1e-07**. The correction is carried
by a tail of 120-200 networks; the rest is the 123 networks that do not re-solve
leaving the population. Section 6's mechanism is untouched — but the numbers are
quoted to three decimals and `phi` 0.0769 -> 0.0825 turns "about thirteenfold too
loose" into twelvefold. Only **1828 of 2884** recorded states were within
`|det J| <= 1e-4` to begin with.

**What connects all three:** each was caught by reporting a count that could have
been quietly omitted — states not evaluable, solver failures, population size.
This is the `figures 6` pattern from D044 in a different domain, and it argues the
count assertions added to the LaTeX build should extend to ensemble runs:
population in, population evaluated, population excluded, asserted rather than
printed.

### D034's failure recurred, in a third guise

Computing `d2R/ds2` by fixing `u` and solving `G = 0` in `a` failed one-sided at
exactly 12 of 325 — the states where the nullcline runs near-vertical. Correcting
along `grad G` by arclength evaluates all 325. Root-finding in `a` at fixed `u`
has now cost this project three separate defects.

### 64 cores

The branch trace is one independent field evaluation per network and was serial:
146.6 s against 11.1 s on 64 networks, **13.3x**. Parallel and serial outputs
compare EQUAL, checked rather than assumed. Behind a `workers` argument with the
serial path preserved so the equivalence check stays runnable.

## D046 — the Hopf is real, the title changes, and folds have two orientations

### The verification, and why it needed one

`branchToFold` reported `tr J >= 0` on the stable branch for 101 of 2765
kinetic-box networks. Three of the four Block 1 numbers before it were method
artefacts, and this one had two candidate mechanisms with identical signatures:
the `u + a` ordering heuristic walking the wrong arc, and contour vertices that
do not lie on `{G = 0}`.

Both were removed rather than argued about. Ordering replaced by the maximal
contiguous `det J > 0` run adjacent to the solved fold — connectivity decides,
nothing guesses. Every vertex reprojected onto `{G = 0}` along `grad G` before
its Jacobian is taken, so each point is an exact equilibrium (`|G| <= 2.4e-14`;
with `j := R`, `du/dt = -G = 0` too). **The count went UP, 101 -> 108.** Four sit
on branches whose `j`-maximum is not the fold and are excluded as
multiplicity-ambiguous, leaving **104 of 2766**.

### The control is what makes the dynamics check mean anything

| | escaped | grew >10x | median max\|delta\|/delta_0 |
|---|---|---|---|
| 104 crossers | 102 | 104 | 9851 |
| 205 stable controls | **0** | **0** | **1.00** |

No distributional overlap. The integration also recovers the eigenvalue from
nothing the contour touched: growth exponent within 5% of `max Re lambda` in 90
of 104, period within 5% of `pi/omega` in 47 of 49.

**A test without a control measures nothing.** The first version ran only on the
crossers, and would have been unfalsifiable had the integrator drifted.

### Four defects in the verifier, all found by running the whole population

- fixed horizon: `lambda_max` spans 5.8e-05 to 3.6, so `t = 50` is 170 e-folds
  for one network and 0.005 for another
- growth read at the ENDPOINT: a slowly growing spiral can finish nearer the
  equilibrium than it started; three were scored "did not grow" on that alone
- whole-window exponent fit: wrong for REAL pairs (median 0.11) because the
  subdominant mode contaminates the early window. Two estimators now — peaks for
  complex pairs (identical phase, unbiased, 2.9e-04), late half for real (1.3e-02)
- **the predicted period was off by exactly 2x.** `|delta|^2` oscillates at
  `2 omega`, so peaks come every `pi/omega`. The measurement disagreed twice and
  was right both times; the prediction was checked third.

Also: the horizon must be CAPPED. `25/lambda` reaches 4.3e+05, and a horizon that
long is spent past the escape, where `solveFreePools` can run to its 2e5-iteration
cap on every RHS call. Stopping at escape is both the cheap fix and the correct
one — leaving the linear neighbourhood answers the question.

### The title

Pre-committed rule fires: **"An Exact Fold Condition for Mass-Balanced Models of
Protein Quality Control"**, and the paper says fold, not collapse threshold.
Multiplicity alone would have forced it.

### The load grid's cleanliness proves less than it looks

0 of 325, and that is a control, not an absence of Hopf in general. The load grid
FIXES kinetics and sweeps load; the crossers are distinguished by KINETIC
parameters (`kappa_a` 11.7x lower, `rho_A` 7.5x higher, `alpha_g/alpha_n` 3.7 ->
33.6). So the load grid is one point in the space where the distinction lives.
The defensible pair, and the form the manuscript uses: an oscillatory instability
precedes the fold in 3.76% of a deliberately wide kinetic ensemble, concentrated
in an identifiable region; and the base kinetic set lies outside that region,
verified over its complete load grid. The 3.76% is a property of the sampling box
and is NOT the headline — the parameter-corner characterisation is.

### Remark 3 becomes a classification, not an exception list

`j_crit = R(u*, a*)`, so `sign(d2R/ds2)` is `sign(d2j/ds2)`. At a constrained
MAXIMUM the equilibrium pair is destroyed as `j` rises: collapse. At a MINIMUM it
is BORN: the lower turning point of a hysteresis loop, not a threshold.

`d2R/ds2 < 0` at 325 of 325 load-grid folds and 2739 of 2765 kinetic-box folds.
At the other 26 the solver returned a birth-oriented point, and **all 26 of those
networks carry a collapse-oriented candidate too, at strictly higher influx in
26 of 26.** They are hysteresis loops whose lower turning point was reached first.
7 of a 153-network control sample carry the same pair. Orientation computed twice
by independent routes: agree 54/54 and 162/167.

Consequence: solving `det J = 0` returns candidates without distinguishing them,
so orientation must be COMPUTED, never assumed from the fact that a solver
converged.

### Three degeneracies, all now instantiated or excluded

`det H = 0` homeostasis (structurally excluded, Section 3.1), `det J = 0` with
`tr J != 0` the fold (computed), `tr J = 0` with `det J > 0` the Hopf (computed).
Section 3.1 says so in one sentence.

## D047 — v5 adopted as the submission source; four defect classes, one structural

### The beta margin: a structural fault, not a lapse

Section 8.3's table stopped at `beta = 0.25` and showed 16x as the closest
approach to the measured `rpoH` load. Figure 5 plotted to `beta = 0.05`, where
the requirement comes within **3.19x**. Prose written from the table therefore
overstated the separation FIVEFOLD, in the paper's only falsifiable prediction.

The same error was corrected in v4 earlier and came back in a clean rewrite. That
is the signature of a structural fault: **two ranges for one quantity will always
let prose land in the wrong place.** The fix is not the sentence. `fig_beta.py`
now owns `BETA_LO/BETA_HI` and `tableRows()`; the manuscript table is generated
from it, the figure sweeps the same endpoints, and the caption echoes the table's
own endpoints rather than restating a span. The table gains its `beta = 0.05` row
(0.913-1.570% of proteome, 3.2x-11x).

### Figure 3 and Section 6 now share one population

`saturation.tsv` is rebuilt at RE-SOLVED fold states, so the figure's medians
(0.180 / 0.159 / 0.049) are the table's. A figure and a table reporting one
quantity over two populations is the P3 divergence; a caption noting the
discrepancy is a note, not a fix. Population 2884 -> 2767, the 117 that do not
re-solve stated. The D039 floor ladder moves 0.0899->0.0763 at 1e-4 and
0.3551->0.3978 at 2e-2, a factor of five rather than four. **The no-screen
decision is unchanged** -- it rests on the absence of a gap and on the median
sliding continuously, neither of which depends on the specific floor values.

### Typography is not structure, and only rendering finds it

Four defects reached the PDF with every count assertion green, a zero-overfull
LaTeX log, and balanced-environment checks passing:

- a `longtable` sharing a page with a top float displaced the page FOLIO into the
  body text -- page 6 printed its "6" over the "u" of "input-output" and carried
  no number at its foot. Short tables are now `table` floats.
- markdown's `\|` table escape reached LaTeX as a DOUBLE bar, so |tr J|,
  |d2R/ds2| and |w_1| -- all scalars -- were set as norms.
- `tr J`, and later `Re(lambda)`, set as products of italic letters. LaTeX has no
  `\tr`, so these belong in `_OTHER_OPS` (`\operatorname{}`), not `_BUILTIN_OPS`.
- `+-i sqrt(det J)` reached the PDF as `\sqrt\{\det J\}` -- an EMPTY radical
  followed by a literal set -- because the sqrt conversion ran before the
  brace-escaping step and its own braces were escaped as set delimiters.

A programmatic scan of all 273 inline math spans for these four shapes now
returns zero, but the scan was written FROM the rendered read, not before it.

### Tests: deleted, not rewritten

16 assertions pinned v4 prose in sections v5 restructured. They were deleted
rather than rewritten, because rewriting means pinning text that has not settled.
What survives is the class that is about correctness rather than wording:
populations named and sized, no maximum without a p99, no figure script reading
the gitignored run root, caption numbers recomputed from their generator.

**Sections are indexed by HEADING TEXT, never by number.** Three renumberings
have now broken assertions; a section's words survive a renumber and its number
does not.

### Token-pinning recurred three times in one session

Once in tests written immediately after the lesson was recorded (`### 3.2`
indexing), once banning a string the manuscript uses correctly, and once on the
word `phase1RunDir` appearing in a COMMENT EXPLAINING ITS OWN REMOVAL. The last
is D042 exactly. That test now parses the AST and asks whether the module CALLS
the function.

### Figure 4 could not have been rebuilt from a clean checkout

`fig_hopf.py` read the gitignored run root directly. The reduction moved into
`build_figure_data.py`, the one script permitted to; `fig4.pdf` hashes IDENTICAL
afterwards, which is how we know the path changed and the figure did not.

### Reference verification (recorded here; the appendix stays removed)

- **Jin J, Rempala GA (2026)** *J Math Biol* 92(3):35,
  doi:10.1007/s00285-026-02352-y. v5 cited "Jin J, Yu PY" -- **the second author
  does not exist on that paper**. The arXiv listing (2407.11248) gives Rempala
  too, so this was not an arXiv-versus-published discrepancy; it looks like a
  conflation with Polly Y. Yu, who works in the same area. Corrected, volume and
  article number added.
- **Xu L, Qu Z (2012)** *PLoS ONE* 7(4):e34616, doi:10.1371/journal.pone.0034616.
  v5 had "Xu Y". Claim verified verbatim against the abstract.
- **Goldbeter A (2013)** *FEBS Lett* 587(17):2778-84,
  doi:10.1016/j.febslet.2013.07.031. Citation correct.
- **Krishna S, Jensen MH, Sneppen K (2006)** *PNAS* 103(29):10840-5,
  doi:10.1073/pnas.0604085103. Citation correct; the abstract supports the claim
  directly ("depending mainly on the saturation of the active degradation rate
  of IkappaB").
- **Cotton MW, Goriely A, Klenerman D, Meisl G (2025)** arXiv:2511.18893.
  Authors and title correct.

**Two claims remain unverified and must not be certified from an abstract:**
Cotton's "eight aggregation and removal mechanisms" and the two named sufficient
conditions (the abstract says only "a wide range"); and the Goldbeter claim that
saturable degradation lowers the cooperativity required for oscillation and
destabilises the steady state as the Michaelis constant falls. Both need the full
text. See `abstracts-are-not-sources-for-numbers`.

## D048 — provenance, and two claims closed

### Cotton et al. 2025: both claims VERIFIED from the full text

Table I lists eight models -- Constant Removal, MM Removal, Mass Saturating,
Combination Removal, Overlapping Aggregates, Conserved Protein, No
Self-Replication, Fragmenting. Section IV C states the two sufficient conditions
explicitly: a self-replicating aggregation mechanism with a monomer
concentration-dependent rate, and a limited-capacity aggregate removal mechanism.
Sections 1 and 3.3 characterise both correctly.

Also confirmed from the full text: their monomer concentration enters primary
nucleation, secondary nucleation and elongation in BOTH moment equations, so H1
fails for their system. Section 3.3's statement stands.

### Goldbeter 2013: claim NARROWED, not verified

The cooperativity clause came through a secondary summary rather than the
article. It is not load-bearing for Section 7 and is removed. What remains is
what the NF-kappaB and ubiquitination sources support independently: saturable
degradation drives a stable steady state unstable as the Michaelis constant of
the degrading step falls.

### OPS_SUBMISSION.md documents a run that produced NOTHING

It records a submission at 16:12:20 writing to `results/phase1/{A,B,C,D}`. Those
four directories were EMPTY -- the launch was reaped by the sandbox -- and were
quarantined at 16:28:56. The run that produced every number in the paper is
`results/phase1/run_20260731T162946-0500`, launched 17 minutes later.

Committing it as "the provenance" would have attributed the results to a run that
generated nothing. It is committed as history; `data/PROVENANCE.md` carries the
identity of the run actually used, extracted from the per-experiment
`provenance.json` records rather than from the operational narrative.

**The run is not single-commit.** A ran at `73e1c0ab`, B and C at `a17dfafd`, D
at `850726c7`, all with a dirty working tree. `OPS_SUBMISSION.md` asserts one
canonical HEAD; that is true only of A. B and C sharing a commit is the case that
matters for comparing the load grid against the kinetic box.

### Read success from the operation, never from a later step

`git push … | tail -3 && echo "pushed"` printed `pushed` while the push failed:
the branch was `master`, not `main`, and the `&&` chained off `tail`. A false
success was reported to the user and acted on. Same shape as a metric that
returns 1.0 exactly where it is quoted -- structurally unable to report failure
at the point failure matters. Written into notes/VERIFICATION_RULE.md.

---

## D049 — The model statement was not the model: §2.1 printed total concentrations

Phase A of the merged review order. The finding blocks Sections 2, 3.3, 6 and 7,
and it is the most damaging kind because nothing in the paper could have caught
it: every number in the manuscript is recomputed from `model.py`, and `model.py`
was right. Only the *printed* model was wrong.

### What was wrong

§2.1 printed `v_ref = k_ref C_f u/(K_ref + u)` — a Michaelis factor in **total**
substrate against free enzyme. `model.py` solves it in **free** substrate
throughout, per decision D004, which forbids exactly the substitution §2.1 made.
The manuscript also never said that `u` and `a` are total pools, free plus
machinery-bound, so a reader had no way to detect the mismatch. At the base
parameter set the two forms of `v_ref` differ by more than 20% at the fold state.

§2.1 now prints the four conservation laws in full and every rate law in free
concentrations. `tests/theory/test_manuscript_model.py` implements §2.1 again
from the printed equations, by naive fixed-point iteration rather than the
safeguarded Newton in `model.py`, and asserts the two right-hand sides agree over
a 676-point grid and over 40 random parameter sets. A third assertion checks that
the *old* printed form does **not** reproduce `rhs`, so the guard cannot pass for
the wrong reason.

### The second defect: the catalytic layer is stipulated, not derived

For `C + U <-> CU -> C + native` the flux is `k_cat[CU]`, which the closure
already makes `k_cat C_f U_f/K_CU`. The coded form saturates a *second* time, in
the same substrate, against a *second* constant. That asserts a mechanism.

`theory/RATE_LAWS.md` writes the reaction network under which the coded form is
exact — nonproductive sequestration counted in all four conservation laws, then
an independent productive cycle on the unheld pool — and then measures its two
load-bearing assumptions at the 2767 solved fold states of the kinetic box:

| assumption | measure | result |
|---|---|---|
| cycles sharing a pool have independent sites | `s_ref + s_dis` | median 0.455, **exceeds 1 in 334 of 2767** (max 1.97) |
| productive complexes hold negligible substrate | `(C_f s_ref + D_f s_u)/u` | median **0.204**, exceeds 0.10 in 1863 of 2767, max 22.7 |

The first could be true — GroEL is a 14-mer — but it binds in a tenth of the box.
The second is not an approximation but a contradiction of a conservation law: the
classical `E_tot << S_tot + K_M` condition cannot hold in a model whose subject is
machinery being scarce relative to burden.

**Verdict: option (b). The rate laws are phenomenological closures.** Theorem 1
and its corollaries are untouched — they need only H1–H3, and mass balance is
exact by construction for either closure. The quantitative layer is not, which is
what D051 measures.

## D050 — Lemma 0: H3 was an assumption about implicitly defined functions

`R` and `G` are not formulae in `(u, a)`. They are defined through a nonlinear
algebraic system the manuscript never stated as such, and nothing established
that the system has a solution, that it is unique, or that it is differentiable.
H3 asserted the conclusion.

All three are provable, and the proof is elementary — `theory/LEMMA0_BINDING.md`:

* **Existence.** The closure map is monotone increasing and carries
  `[0, C_tot] × [0, D_tot]` into itself; Knaster–Tarski.
* **Uniqueness.** It is strictly sublinear, because scaling the argument by
  `λ < 1` shrinks each denominator's occupancy terms but not its additive
  constant. A monotone strictly sublinear map has at most one positive fixed
  point.
* **Regularity.** All four partials of the free substrates in the free pools are
  negative, so the binding Jacobian is a Z-matrix; contracting it against
  `(C_f, D_f)` and using that each free substrate is homogeneous of degree `−1`
  in `(1, C_f, D_f)` gives row values `C_f(1+ν+…)` and `D_f(1+…)`, both positive.
  A Z-matrix with a positive vector of positive image is a nonsingular M-matrix.

The last step yields a sharper statement than well-posedness:

> **det of the binding Jacobian ≥ 1 + ν, everywhere.**

Verified over 24,300 states spanning four decades of every dissociation constant:
smallest observed margin `1.1×10⁻⁶`, largest off-diagonal `−8.0×10⁻⁹`, smallest
row-contraction margin `3.2×10⁻⁹` — all three attained in the vanishing-burden
limit where the bounds are equalities. The proof does not rest on the check.

Two consequences beyond H3. The runtime singularity guard in `model.py:jacobian`
**cannot fire**; it is kept as a defensive assertion and annotated as such.
And inverse-positivity gives `∂u_f/∂u > 0` and `∂a_f/∂u ≥ 0` — raising soluble
burden titrates machinery away and *raises* free aggregate. That second
inequality is what kills §3.3's sign argument (task B6).

## D051 — The twelvefold ceiling factor is closure-dependent. It is 2.8× under the alternative.

Task A3, required because D049 landed as option (b). Same kinetic box, same
draws, same seeding from the same recorded states; the two pipelines differ in
exactly one field of `Params`. The alternative closure is
`v_ref = C_f U_f/K_CU = [CU]` — catalytic flux proportional to bound complex,
which removes both stipulations at once, since there is no second saturation and
cycles sharing a pool compete automatically. The removal ceiling
`c_tot + (ρ_U + ρ_A)d_tot` is unchanged, because every complex concentration is
bounded by its own pool total, so `φ` means the same thing in both.

| quantity, median at the fold | michaelis (published) | bound complex |
|---|---|---|
| `φ = j_crit/ceiling` | 0.0825 | **0.353** |
| `1/φ`, the ceiling overestimate | **12.1×** | **2.8×** |
| refolding utilisation `v/V_max` | 0.0363 | 0.256 |
| soluble degradation utilisation | 0.0584 | 0.330 |
| aggregate clearance utilisation | 0.0190 | 0.0658 |
| networks yielding a solvable fold state | 2767 of 2884 | 2027 of 2884 |

Paired over the 1991 networks solved under both, the median ratio is **3.81×**
and 1516 of 1991 lie outside a factor of two. The work order's threshold was a
factor of two. It is exceeded by the median alone.

**The twelvefold claim is a property of the closure, not of the model class**, and
must say so where it is made — abstract, §6 headline, and Discussion.

Note what is *not* closure-dependent: the identity, the converse, the `n`-state
form, the dilution corollaries, and Corollary 4. Those follow from H1–H3.

### A third finding, forced by the same computation

"The fold occurs while the machinery runs at roughly 5–18% of `V_max`" is not
what `s_ref = 0.180`, `s_u = 0.159`, `s_a = 0.049` measure. Those are Michaelis
saturation factors. Utilisation of nominal capacity is `(c_f/c_tot)·s`, which is
**0.036 / 0.058 / 0.019** — about a fifth of the quoted figures. `φ = 0.0825` is
the honest utilisation measure and `1/0.0825 = 12.1` is where "twelvefold" comes
from, so the abstract currently reports two different quantities as one, and its
own arithmetic does not close: `1/0.18` to `1/0.05` is 5.6× to 20×, not 12×.
This is task E2(a), and it is confirmed rather than merely alleged.

## D052 — Theorem 1 was three claims with different hypotheses wearing one label

Task B1. Two defects, one of them a contradiction the paper carried in plain
sight.

**H1 contradicted §4.3.** H1 said `j` "does not appear in `da/dt`". §4.3 lets
`C_tot` depend on `j`, which puts `j` in `da/dt`, yet §3.2 said transversality
reduces to `w₁ ≠ 0` *because* `∂F/∂j = (1,0)` exactly. Both could not hold. The
identity survives because the gradients are taken with respect to the state at
fixed `j`, not because H1 holds. H1 is now split at the point of use: **(H1a)**
total influx is state-independent, which is what puts mass balance in the form H2
needs; **(H1b)** `j` is absent from `da/dt`, so the nullcline does not move.
Only (H1a) enters the identity.

**"`det J = 0` holds exactly at constrained critical points" is false where
`∇G = 0`.** The determinant vanishes automatically there, `{G = 0}` need not be a
manifold, and the Lagrange reading is vacuous — as the paper's own two
exceptions demonstrate. Regularity is a hypothesis of the forward direction, not
only of the converse. The claim that it is needed only for the converse is gone.

The replacement is four statements, collectively Theorem 1:

| | hypotheses | content |
|---|---|---|
| 1.1 | (H1a), H2, H3 | `det J = −(∇R × ∇G)` |
| 1.2 | + `∇G ≠ 0` | `r'(s) = det J/‖∇G‖` **pointwise** |
| 1.3 | + (H1b) | candidates solve `{G = 0, det J = 0}` in the state alone |
| 1.4 | `∂G/∂j ≠ 0` | three equations in `(u, a, j)` |

**Theorem 1.2 is the strengthening and it costs three lines.** `det J` *is* the
derivative of total removal along the aggregate nullcline, up to the positive
factor `‖∇G‖` — at every point, not only where it vanishes. The paper's slogan
was a gloss on a vanishing condition; it is an exact identity. Two things follow
immediately: the equilibrium branch is the level-set inverse of one scalar
function of one variable, and `det J = 0` is `r'(s) = 0`.

## D053 — Theorem 2 was asserted. It is now proved, and (G3) turns out to be Sotomayor's condition.

Task B2. Numerical verification of (G1)–(G4) is evidence that the conditions
hold. It is not a proof that they imply the normal form, and the paper had no
proof.

The proof is half a page and needs no citation, because of Step 4. With
`w = (1, 1+λ)` the left null vector and `v` the right one,

```
w·D²F(v,v) = −D²R(v,v) + λD²G(v,v)
r''(0)     =  D²R(v,v) − λD²G(v,v)     [from ∇G·γ'' = −D²G(v,v)]
⟹ r''(0) = −w·D²F(v,v)
```

so **(G3) is not merely equivalent to the classical nondegeneracy condition of
the saddle-node but equal to it up to sign**. Step 3 is what makes the comparison
legitimate: the second row of `Jv = 0` reads `∇G·v = 0`, so the right null vector
is tangent to the nullcline and the centre manifold and the arclength reduction
are the same object.

The orientation classification of §3.4 becomes a corollary of Step 1 rather than
a separate observation, and Table 1 stops carrying the weight of the converse.

## D054 — (G4) follows from (G1). The genericity list is three conditions, not four.

Task B3. §3.2 said "(G4) is not automatic … no structural feature of the model
forces it." That was wrong, and the paper's own numbers were the evidence.

**Proof.** `∂F/∂j = e₁`, so (G4) is `w₁ ≠ 0`. With `L = I + e₁e₂ᵀ`, `LJ` has first
row `−∇R` and the rest the non-influx constraint gradients, and `L` fixes the
first coordinate of a left-multiplying row vector. Put `z = wL⁻¹`; then
`z(LJ) = 0`, and `z₁ = 0` would make the remaining components annihilate
full-rank gradients. So `z₁ ≠ 0` and `w₁ = z₁ ≠ 0`. **Every state dimension.**

The margin has a closed form: `|w₁|/‖w‖ = (1 + (1+λ)²)^(−1/2)`, small exactly
when `|λ| → ∞`, i.e. when `∇G → 0`. **The co-occurrence of the (G1) and (G4)
near-failures is a theorem, not an observation about two networks.**

Verified by two independent routes — SVD of the analytic Jacobian against
central-difference gradients — at all 325 load-grid and 2765 of 2767 kinetic-box
states, median relative difference 8.5×10⁻¹¹, maximum 1.8×10⁻⁵. The two
disagreements are exactly Table 1's two exceptions, at `‖∇G‖` of 1.04×10⁻⁸ and
3.99×10⁻⁸ against a smallest value of 1.41×10⁻⁶ elsewhere: there the right-hand
side is differencing noise. The closed form fails only where it has nothing to
evaluate, which is confirmation rather than qualification.

### B7, from the same computation: Table 1 reported `|tr J|` and hid the classification

At a two-dimensional fold the eigenvalues are `0` and `tr J`. A negative trace is
a stable node colliding with a saddle — the fold terminates a stable branch. A
positive trace is an unstable node colliding with a saddle, and **the fold is not
a stability boundary at all**. Table 1 now reports the signed range.

Load grid: all 325 negative, −0.365 to −0.303. Kinetic box: −626 to **+3.41**,
with **61 of 2767 positive**.

All 61 are among §7's networks that lose stability before the fold — as they must
be, since a positive trace at the fold puts the branch maximum above zero, so
that count is a consistency check between two independent computations and not a
finding. The finding is the converse, and the work order's expectation does not
survive it: it supposed the 104 crossers "have `tr J > 0` at their folds by
construction". **47 of the 108 crossers have `tr J < 0` at the fold** — the branch
destabilises and recovers before reaching it. An oscillatory excursion below
`j_crit` does not imply the fold itself is unstable.

## D055 — Corollary 2's converse is false above two dimensions, and §4.1 carried a 0 = 0 check

Task B4. "The converse requires, in addition to (G2)–(G4), full rank" uses
`tr J ≠ 0` for two jobs that only coincide in the plane. Writing
`p(λ) = λⁿ + … + c₁λ + c₀`, a zero eigenvalue is `c₀ = 0` and is algebraically
simple exactly when `c₁ ≠ 0`, where `c₁` is up to sign the sum of the
`(n−1)×(n−1)` principal minors; partial hyperbolicity is separate. For `n = 2`,
`c₁ = −tr J` and the remaining eigenvalue is `tr J`, which is why the error was
invisible. As stated, the converse was false for the two three-state systems §4.1
verifies.

The corrected conditions hold at every fold state of both: `|c₁| ≥ 9.61×10⁻³`
(regulated), `≥ 0.132` (strong controller), `≥ 2.13×10⁻²` (sequestered, 46 fold
states from 63 seeded attempts), with the rest of the spectrum off the imaginary
axis by at least 1.62×10⁻².

**A second finding, nearly reported as a result.** The `σ₀ = 0` and `k_seq = 0`
settings make a row of the Jacobian identically zero, so `det J = 0` everywhere
and a fold solve converges to an arbitrary point of the nullcline. Its 24
"fold states" are not fold states. Excluding them was not fastidiousness — the
first run included them and produced a `|c₁|` margin of 1.5×10⁻⁴ that measured
nothing.

The same structural zero makes §4.1's unregulated identity check vacuous:
`∇C = 0`, so both sides are exactly zero. **"median relative error 0.000×10⁰"
was `0 = 0`** — the defect the work order flags for Table 2 row 1, present in a
second place. §4.1 now says what that limit does and does not test.

## D056 — Transversality under load-dependent capacity, in closed form

Task B5. The handoff said "use the full `F_j`". It evaluates:
`F_j = (1 − R_j − G_j, G_j)` and `w = (1, 1+λ)` give

```
w·F_j = 1 − R_j + λG_j,    so transversality is    R_j − λG_j ≠ 1.
```

That is the exact statement that D054's automatic transversality is lost
precisely when capacity feeds back on the influx.

Measured at solved self-damaged fold states, 99 converging from 200 attempts. In
**burden mode** the capacity factor contains no `j`, so `R_j = G_j = 0` exactly
and `w·F_j = 1` at all 46 states — (G4) returning for the reason it held in §3.2.
In **influx mode** the margin *rises* with the coupling, 1.000 at `ε = 0.01` to
1.38 at `ε = 100` over 53 states. Minimum over the whole population: 1.0. The
feedback moves transversality away from failure in this model; that it could move
the other way is what the condition is for.

## D057 — C0a: 104 and 108 were never in conflict. Neither is superseded.

The work order's hypothesis was that 104 came from `genericity.py:branchToFold` and
108 from the re-solved fold states. It is wrong. Both come from the same file,
`data/computed/hopf_refined_kinetic_box.tsv`, produced by one computation. They
are different quantities and neither was labelled:

| number | what it is |
|---|---|
| 2767 | fold states solved over the kinetic box |
| 2766 | of those whose low-burden branch traces — **1 trace failure** |
| 108 | of those reaching `tr J = 0` before the fold |
| 104 | of *those*, whose branch's influx maximum is the fold |

The 108 → 104 step is a filter on meaning, not on quality. Where the branch's
`j`-maximum is not the fold, a second singular point terminates the run the walk
returned, so "the crossing precedes the fold" concerns a different branch's
endpoint. Those 4 were excluded from the integration and from §7's headline.
`104/2766 = 3.76%`, `108/2766 = 3.90%`. Both are now labelled wherever they
appear, in §7, Fig. 4's caption and Remark 3.

## D058 — C0b: the crossers are two populations, and one of them has a window

Splitting by the sign of `tr J` at the fold (task B7):

- **61 terminal** — one crossing, never recovers; the fold is reached by an
  already-unstable branch.
- **47 window** — `tr J` returns below zero before the fold, so the branch carries
  an **instability window strictly interior to it**.

The parity check earned its place. A terminal branch must cross an odd number of
times and a window branch an even number ≥ 2. All 61 terminal networks cross
exactly once. Only **38 of the 47** window networks cross exactly twice; the other
nine are branches the trace does not resolve, and they are reported rather than
assigned: 2 are the multiplicity-ambiguous cases, 2 leave the window through the
high-burden side, 5 stop short of low burden. Their single traced crossing is a
*downward* one, which is the same fact seen from the other side and is why all 7
of the truncated ones show negative transversality at their only crossing.

Window geometry over the 38 resolved: band width 0.0179 to 0.658 of `j_crit`,
median **0.199**, opening at median 0.64 and closing at 0.95.

**The escape statistic cannot answer the question it was asked.** The work order's
decision rule was that if the window group escapes at the terminal group's rate,
both crossings are subcritical. It does escape at the same rate — 42 of 43 against
60 of 61 — and the inference does not follow. `d_escape` is 10⁻² of the state
scale against a 10⁻⁶ perturbation, so "escaped" means a ten-thousandfold
amplification, which growth onto a stable limit cycle produces exactly as
reliably as divergence. The test measures departure from a small neighbourhood.
Criticality needed `l₁`, and `l₁` says the opposite (D059).

## D059 — C2 and C3: 120 of 145 crossings are supercritical. "Sustained" was right.

`hopf_check.py` said a `tr J = 0` crossing with `det J > 0` is "a Hopf point by
definition, not an inference". That is false — it is necessary, not sufficient —
and it is corrected in the docstring in the same commit.

Computed per crossing, not per network, because a window network has two and they
have opposite transversality sign by construction.

**C2.** All 146 crossings transversal, `|d(tr J)/dj| ≥ 0.293`. Sign positive at
all 61 terminal crossings and at the 40 opening crossings; negative at all 38
closing ones. The 7 remaining negatives are the truncated branches of D058.

**C3.** `l₁` by the planar formula (Guckenheimer and Holmes 3.4.11), evaluable at
145 of 146 with the sign stable across a factor of six in step size at every one.
The differentiator was validated first against the two textbook normal forms,
recovering `a = −1` and `a = +1` exactly; a convention error here reverses the
paper's conclusion, so the self-test is asserted before any model result is
computed.

| | supercritical | subcritical |
|---|---|---|
| terminal, single crossing | 39 | 21 |
| window, opening crossing | 46 | 1 |
| window, closing crossing | 35 | 3 |
| **all** | **120** | **25** |

So a **stable limit cycle occupies the window**: born as the influx enters it,
destroyed as it leaves. Burden oscillates over a band of influx and then stops
oscillating as the influx rises further, before the fold. That licenses
"sustained" for the window group — the word the earlier review ordered deleted,
deleted on the strength of an escape statistic that could not see the difference.
For the terminal group about a third of crossings are subcritical and lose
stability hard, with no local attractor past the crossing, so the word must not
be used there.

## D060 — C4: three degeneracies, three different networks, no common centre

| degeneracy | condition | closest approach | network |
|---|---|---|---|
| Bogdanov–Takens | `tr J = 0` at a fold | `\|tr J\| = 8.38×10⁻⁴` | draw861 |
| Hopf-pair annihilation | `max(tr J) → 0` on a window branch | `4.56×10⁻⁴` | draw2863 |
| Bautin | `l₁ = 0` at a crossing | `5.65×10⁻⁵` | draw820 |

They are reached at three different networks, so they are three boundaries and
not one organising centre; the manuscript says so rather than presenting a single
attractor of the geometry.

One clarification falls out. `det J = 0` holds at *every* fold, so the distance to
a Bogdanov–Takens point just *is* the (G2) margin — Table 1 row 2 has a geometric
meaning it was not given, and it is exactly where Theorem 2's converse fails,
since (G2) is the condition that breaks there. Closing crossings sit an order of
magnitude nearer the Bautin boundary than opening ones, median `|l₁|` 0.051
against 0.858.

## D061 — Goldbeter 2013: not verifiable here; claim narrowed to the abstract's own sentence

The article has no PubMed Central record, so full text is unavailable through the
tooling in this environment and the claim cannot be verified as instructed. What
the record does support, from the abstract of Goldbeter A (2013) *FEBS Lett*
587(17):2778–84, doi:10.1016/j.febslet.2013.07.031, retrieved via PubMed: the
oscillations in question arise "as a result of feedback regulation" — allosteric
product activation — and "the Michaelian nature of the product decay term markedly
influences the period, amplitude and waveform of the oscillations".

That is a modulation claim about an already-oscillating system, not the
destabilisation claim §7 attributed to it. The destabilisation claim now rests on
Krishna et al. 2006 and Xu and Qu 2012, which support it independently, and
Goldbeter is cited for what its own abstract states. This is a citation of a
qualitative statement made unambiguously in an abstract, not of a number.

## D062 — B8: three theorem statements were formally wrong, and one of them was mine

An external review of the Phase A–C state found six defects. Three were errors in
statements I had written in this pass, which is the useful kind of finding: the
mathematics was right and the statement of it was not.

### The centre-manifold identification was false (B8.1)

§3.2 said "It is the nullcline itself … the centre manifold and the arclength
reduction are the same object." `{G = 0}` is **not invariant**. Along the flow

```
dG/dt = ∇G·(du/dt, da/dt) = G_u·(j − R),
```

which vanishes only at the equilibrium, so trajectories cross the curve rather
than run along it. Confirmed numerically before removal: `dG/dt` is `10⁻⁴` to
`10⁻⁶` at nullcline points either side of the base fold and `10⁻¹⁸` at the fold
itself. A centre manifold is invariant; this is not one.

What survives is the tangency at `x*`, which is all the proof ever needed. The
proof is restated in four steps and cites Sotomayor for the normal form while
keeping `r''(0) = −w·D²F(v,v)` as the reason no citation is needed for the
**nondegeneracy condition**. The correct global statement — that the curve of
equilibria lies in the extended centre manifold of the `dj/dt = 0` suspension —
is given in its place.

### (H1b) was used unnamed in three statements (B8.2)

Theorem 1.1 correctly needs (H1a) and H2 only. But Theorem 1.2's corollary that
the branch is a level-set inverse needs a *fixed* `r`, Theorem 2's proof needs
`F_j = e₁`, and the transversality corollary opens with it. All three now name
it; the level-set statement moved to Theorem 1.3, where (H1b) is a hypothesis.

### Corollary 2 did not cover its own second verification (B8.3)

It assumed `du/dt + da/dt = j − R`, while the sequestered three-state system it
verifies satisfies `du/dt + da_r/dt + da_s/dt = j − R`. Generalised to a constant
**mass vector** `m` with `m₁ = 1` and `mᵀF = j − R`; the row operation
`row₁ → row₁ + Σ_{k≥2} m_k row_k` preserves the determinant *because* `m₁ = 1`,
and yields `−∇R`. The regulated system is `m = (1,1,0)` and the sequestered
`m = (1,1,1)`, so both verifications are now instances of the corollary. Full
rank moved from a remark after the equivalence to a hypothesis of it — where the
non-influx gradients are dependent the determinant vanishes for every `∇R` and
the Lagrange reading is empty, exactly as `∇G = 0` empties it in the plane.
Partial hyperbolicity added to the converse summary.

### B8.5: a coefficient was being reported as a margin

`|w·F_j|` under `w₁ = 1` proves nonvanishing but a left null vector rescales
freely, so its size means nothing. The scale-invariant cosine
`|w·F_j|/(‖w‖‖F_j‖)` is a different number: minimum **0.0070** against a
coefficient that never drops below 1.0. Burden mode's exact 1.0 at all 46 states
is the §3.2 corollary reappearing — with no `j` in the capacity factor,
`R_j = G_j = 0` — and is now flagged as such rather than presented as data.

### B8.6: a test pinned the claim B3 disproved

`testTransversalityIsTheConditionWithNoStructuralGuarantee` asserted that the
manuscript flags (G4) as having no structural guarantee, and passed for as long
as the manuscript was wrong. Replaced by three tests pinning the result against
its generator. The audit then found **two more** prose-pinning tests that failed
on *correct* v5 text — both rewritten to pin the property rather than the
sentence. The rule earns its place again: a test that asserts the manuscript
contains a sentence will pass while the sentence is wrong.

### Two items of the review's own were stale, and one was wrong

`.zenodo.json` was already at 1.2.0 with a consistent description as of
`aa97b09`; no change was made. And the handoff's parenthetical that the ORCID
"was set as plain text" was false — it was absent from the PDF entirely. It now
renders, on the **fifth** placement: `\equalcont` takes arbitrary text, needs no
vendored graphic, and prints a dagger footnote under the corresponding-author
line with no label of its own. Verified by rendering page 1 and reading it. The
four earlier failures are kept in the source comment because each is a distinct
failure mode, and one of them — `\footnotetext` — succeeded at the build while
putting nothing on the page.

## D063 — C5d: the small transversality cosine is (G1) again, not a new caution

The smallest scale-invariant transversality margins are not scattered. With
`w = (1, 1+λ)` the cosine falls as `|λ|` rises, and `|λ|` rises as `‖∇G‖` falls
at fixed `‖∇R‖`. Over the 99 self-damaged fold states, log cosine correlates with
log `‖∇G‖` at **+0.90** and with log `|λ|` at **−0.93**, and the single smallest
cosine is the single largest `|λ|` — the same state, `‖∇G‖ = 2.7×10⁻³`.

So Table 1's two exceptions, the two failures of the closed form in D054, and
this margin are **one geometric fact reported three times**: the constraint
losing regularity. §4.3 says it once and drops the reassurance that the margin
is never small, which was true only of the coefficient.
