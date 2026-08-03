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
parameters and has NOT been surveyed across the box.

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
