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
