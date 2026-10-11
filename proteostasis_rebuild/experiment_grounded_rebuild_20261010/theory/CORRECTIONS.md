# Corrections, 2026-10-10

An independent physics audit (`../independent_review/AUDIT_A6_A7.md`, `results_a6_a7.json`,
18/18 tests, re-run here and still passing) showed that the overload section (A6) and the
maintenance section (A7) of this theory were wrong. A6 was not a small slip: its central claim
was false. This file says what was wrong, what replaces it, and how it was checked. Every
coefficient is hypothetical and dimensionless. Nothing here is an E. coli estimate, and none of it
is biological validation.

## The errors

**E1. Chaperone conservation could go negative (A1.1).** The old law,
M^tot = M^free + ΣC + βΣA with constant β, made each aggregated chain hold β units however much
aggregate there was. Once βΣA > M^tot − ΣC, free chaperone is negative, which is physically
impossible.
*Fix:* sequestration is capped at what the pool can supply. B = min(βΣA, M^tot − ΣC), so
M^free = max(M^tot − ΣC − βΣA, 0). β is now a binding capacity, not a fixed occupancy.
*Evidence:* in the full one-client model, the capped form keeps every pool ≥ 0. The uncapped form,
in the same instance, drives bound chaperone to −0.119 and fails both positivity checks. I ran this
as a mutation of `checks.py`; it is not a standing check.

**E2. The A6 "runaway" and "two steady states or none" were artefacts of E1.** The old A6 called
the right side of dA/dt convex on [0, a/b). The pole at a/b, which drove the runaway, sits
a0/(k_c β) beyond A_max = M^tot/β, in the region of negative free chaperone. With the pool capped,
the dynamics above A_max are linear: dA/dt = k_a L/(k_a + μ) − (μ + δ)A. Whenever μ + δ > 0 that
gives a finite, stable high state. There is always at least one steady state. "None" never happens,
and nothing diverges.

**E3. The fold needed a condition A6 never stated.** The saddle-node lies in the physical region
iff k_c M^tot > k_a + μ. The old check `bistable_below_none_above` used the lumped numbers a = 3,
b = 1, k_a = 1, μ + δ = 1 without saying how a splits into k_c M^tot and k_a + μ.
- In the split where a fold exists, the "unstable threshold" it counted at 0.9 L_crit sits at
  A = 1.974 > A_max = 1.75. It is unphysical, and the true count there is 1 state, not 2.
- In another split with the same four numbers, there is no fold at any load.

So that check passed without testing A6's claim. **It has been removed.**

**E4. The maintenance condition (A7.1) was incomplete.** When the required native level depends
on growth rate, the quota boundary moves. The condition on that face must include
−(dN^req/dμ)(dμ/dt). The old form omitted it. In the counterexample, the old condition holds with
zero slack and the state still leaves the viable set at once. The old form was also written "for
every x ∈ ∂K". It applies only where a quota is active, and on its own it is necessary, not
sufficient.

**E5. Native protein was never degraded (A1.3).** There was no N → D flux. Every quota therefore
missed a k_dN·N^req term. The slow-growth limit was wrong as well: it is (k_u + k_dN)N^req(0), not
k_u N^req. An aggregate-degradation flux k_dA·A → D was also added, so that δ in A6 has a
meaning in the full model.

**E6. "Graceful failure" was asserted, not derived.** A2 said synthesis scales with μ. A6 said
performance falls smoothly "until the quota cannot be met even as μ → 0". With strictly
growth-proportional synthesis, growth rates near zero are never viable. The viable interval of
growth rates can also vanish abruptly at a positive μ* = √((k_u + k_dN)n0/n1). Smooth decline is
now stated as conditional: it needs basal synthesis y s⁰ ≥ (k_u + k_dN)n0 and a stable growth
closure, neither of which is established.

**E7. Approximations in the reduction were not declared.** A6 dropped the client-bound term ΣC
from the pool and treated disaggregation as removal. Both are now stated explicitly, and the fold
condition in E3 is labelled as belonging to the reduced model only. Neither effect has been
tested, by this correction or by the audit.

## The corrected result (EQUATIONS A6, Theorem A6)

In the reduced model, dA/dt = k_a L/(k_c max(M^tot − βA, 0) + k_a + μ) − (μ + δ)A, with β > 0.
- If μ + δ > 0: there are 1 to 3 steady states at every load, and all trajectories stay bounded.
  Above L_sat = (μ+δ)M^tot(k_a+μ)/(β k_a) there is a stable high state
  A_hi = k_a L/((k_a+μ)(μ+δ)).
- If k_c M^tot ≤ k_a + μ: there is one globally stable state at every load. No fold, no
  hysteresis.
- If k_c M^tot > k_a + μ: there are three states (stable low, unstable middle, stable high) exactly
  for L_sat < L < L_crit = (μ+δ)(k_c M^tot + k_a + μ)²/(4 k_c β k_a). Outside that window there is
  one state. The window's relative width is ((k_c M^tot − k_a − μ)/(k_c M^tot + k_a + μ))².
- If μ + δ = 0: aggregate grows without bound. This is the only unbounded case.

Also new: **Proposition A1**, exact forward invariance of the nonnegative state set with the
capped pool. **(A3.1)** is now P ≥ (μ + k_u + k_dN)N^req(μ). **(A7.1)** is now
P − (μ + k_u + k_dN)N^req − (dN^req/dμ)(dμ/dt) ≥ 0 on active quota faces.

## Test verdict

`python checks.py` gives **22/22 pass**. The old check `bistable_below_none_above` is removed.
These checks were added:

| check | what it shows |
|---|---|
| a6_symbolic_fold_and_cap | sympy: A_c, L_crit, L_sat, A_max − A_c, loop width, pole position, A_hi |
| a6_capped_root_counts_three_state_window | 1 / 3 / 1 states at 0.9 / 0.99 / 1.1 L_crit; stable/unstable/stable; high root > A_max; brute-force sign changes agree |
| a6_no_fold_without_capture_condition | same lumped numbers, condition reversed: 1 state at all 300 loads |
| a6_high_state_finite_and_attracting | from A = 0 at 1.1 L_crit, A → 1.98 = A_hi; free chaperone never < 0 |
| a6_hysteresis_loop_between_Lsat_and_Lcrit | up-jump 2.263 (L_crit 2.25), down-jump 2.192 (L_sat 2.1875), step 0.014; upper branch finite |
| a6_unbounded_only_without_removal | μ = δ = 0: A(200) = 898 (513 × A_max) |
| a6_random_topology_local_2000 | 0 violations, including "3 states ⇔ condition and window" |
| a6_matches_independent_audit | numbers agree with the auditor to 10⁻¹²; the auditor's 20,000-draw sweep has 0 violations (1,796 three-state cases). It is cross-read, not re-run |
| a1_conservation_and_positivity | full A1 model with cap active: mass-balance drift < 2·10⁻¹³; no pool < 0; C ≤ M^tot |
| a1_positivity_from_boundary | boundary start: minimum −2·10⁻¹³ (round-off; the check fails below −10⁻¹⁰) |
| a3_native_quota_includes_degradation | P = (μ + k_u + k_dN)N at steady state; the old quota is short by the degradation flux |
| a7_face_condition_symbolic | sympy derivation of the moving-face condition |
| a7_moving_quota_counterexample_and_fix | old condition: zero slack, state leaves V (gap −1.88); corrected condition met on the whole face: stays in V |
| a7_feasible_growth_interval | y_crit and μ* symbolic; interval vanishes at μ* > 0; matches the auditor at k_dN = 0; k_dN raises μ* |

The eight checks of sections 1–4 (HtpG sign rule, GroEL sign, saturation, feasible-not-sufficient,
deletion, Little's law) are unchanged and still pass.

## What did not change, and what is still owed

- §§1–4 of the theory and A2, A4, A5 and A8 stand. Only wording was changed, for consistency with
  the cap and the basal synthesis term.
- **Hysteresis prediction (P3): survives as a conditional prediction, in narrower form.** The loop
  is between two finite states and needs one more unmeasured condition. The window may be narrow.
  No experiment has been run.
- Not derived: the fold condition with the client-bound term kept, with disaggregated chains
  recycled to U, with several pools, with σ32 regulation, or with a growth closure μ(A). There is
  no calibrated E. coli parameter set. The earlier draft's ×25 headroom stays withdrawn.
