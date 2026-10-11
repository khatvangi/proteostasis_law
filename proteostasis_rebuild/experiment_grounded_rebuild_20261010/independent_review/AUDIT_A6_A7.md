# Independent audit of EQUATIONS.md, sections A6 and A7

Audited file: `theory/EQUATIONS.md` (not edited). Scripts: `audit_a6_a7.py` (writes
`results_a6_a7.json`), tests: `python -m unittest test_audit_a6_a7 -v` (18/18 pass, ~30 s).
Every coefficient below is hypothetical and dimensionless. Nothing here is an E. coli estimate,
and none of it is biological validation. These results are about the mathematics only.

## Bottom line

1. **The A6 "runaway" is an artefact of letting the free chaperone pool go negative.** In the
   physical model (free chaperone ≥ 0), any positive dilution or aggregate removal (μ + δ > 0)
   gives a **finite, stable high-aggregate state** for every load. There is always at least one
   steady state. "Two steady states or none" is false.
2. **The saddle-node survives, but only under a condition A6 never states:** k_c·M_tot > k_a + μ.
   Without it there is no fold at all, only one steady state at every load. The A6 bistability
   claim therefore depends on one more unmeasured inequality besides β > 0.
3. What the fold really produces is **hysteresis between a low and a titrated high state**, inside
   a load window [L_sat, L_crit]. It is not divergence.
4. **A7.1 is not the correct invariance condition** when the required native level depends on growth
   rate. It misses the term (dN^req/dμ)·dμ/dt. I built a case where A7.1 holds exactly and the
   state still leaves the viable set at once.
5. **No native-protein degradation exists in A1.3.** Every quota (A3.1, A7.1, and the μ → 0 limit)
   leaves out a k_dN·N^req term.

## A6.1 Exact physical domain of the written saddle-node

Notation: a0 = k_a + μ, a = k_c M_tot + a0, b = k_c β, m = μ + δ.

- Free chaperone M_free = M_tot − βA is nonnegative only for **A ≤ A_max = M_tot/β**.
- A6 says the right side is convex "on [0, a/b)". But a/b − A_max = a0/(k_c β) > 0. So the
  whole stretch (A_max, a/b), including the pole that drives the "runaway", has negative free
  chaperone. [sympy, `pole_ok`]
- The saddle-node A_c = a/(2b), L_crit = m a²/(4 b k_a) is reproduced exactly. But

      A_max − A_c = (k_c M_tot − k_a − μ) / (2 k_c β)            [sympy, `A_gap_ok`]

  so **the fold lies in the physical domain iff k_c M_tot ≥ k_a + μ**. In words: chaperone
  capture at full pool has to beat the spontaneous exits (aggregation plus dilution) combined.
- `checks.py` tests a = 3, b = 1, k_a = 1, μ+δ = 1 and never says how a splits into k_c M_tot
  and k_a + μ. That split decides the answer:
  - P_FOLD (k_c M_tot = 1.75, k_a+μ = 1.25): A_max = 1.75. The "unstable threshold" that
    checks.py counts at 0.9 L_crit sits at A = 1.974 > A_max, so it is unphysical.
  - P_NOFOLD (k_c M_tot = 1.10, k_a+μ = 1.90): has the same four numbers as checks.py and gives
    no fold at any load (300 loads tested).

  So `bistable_below_none_above` passing does not test what A6 claims.

## A6.2 Saturated pool M_free = max(M_tot − βA, 0)

For A > A_max the clearance term vanishes: U = L/a0 and dA/dt = k_a L/a0 − mA. That is linear
with slope −m.

**Proof of the topology.** Steady states are solutions of H(A) = L, where
H(A) = (m/k_a)·A·(k_c·max(M_tot − βA, 0) + a0) is the load needed to hold aggregate A. H is
continuous, H(0) = 0 and H → ∞, so **at least one steady state exists for every L whenever m > 0**.
- If k_c M_tot ≤ a0, H increases everywhere. There is a unique, stable steady state at every load,
  and no fold.
- If k_c M_tot > a0, H rises to L_crit at A_c, falls to L_sat = H(A_max) = m M_tot a0/(β k_a), then
  rises linearly. For L in (L_sat, L_crit) there are three states: low stable, middle unstable, and
  high stable at A_max or above. Otherwise there is exactly one.
- L_crit − L_sat = m(k_c M_tot − a0)²/(4 k_c β k_a) ≥ 0. [sympy, `L_gap_ok`]

**The high state** is A_hi = k_a L / ((k_a+μ)(μ+δ)). That is exactly the aggregate level of a cell
with no free chaperone. It is finite and linear in L.

Numerical confirmation (results_a6_a7.json):

| check | result |
|---|---|
| P_FOLD steady states at 0.9 / 0.99 / 1.1 L_crit | 1 / 3 / 1 (checks.py asserts 2 / – / 0) |
| closed-form count vs brute-force sign changes on a 200k grid | agree in all 10 cases |
| ODE from A=0 at 1.1 L_crit, saturated | settles at A = 1.98 = A_hi |
| same ODE, written formula | hits the pole at t ≈ 8, M_free = −1.25 |
| quasi-static sweep up / down | jumps up at L ≈ L_crit (2.263 vs 2.25), down at L ≈ L_sat (2.192 vs 2.1875), within grid step |
| 20,000 random positive parameter sets, L from 0.01 to 10 × L_crit | 0 violations of: ≥1 root, ≤3 roots, #stable = #unstable + 1, unique high state above L_crit, no fold when k_c M_tot ≤ k_a+μ (1,796 bistable cases) |

**When does a real runaway remain?** Only when μ + δ = 0 on the titrated branch. Then dA/dt =
k_a L/a0 > 0 forever (tested: A reaches > 500·A_max). I also made disaggregation depend on free
chaperone, δ = δ_deg + δ_dis·M_free/M_tot. The high state stays finite, k_a L/((k_a+μ)(μ+δ_deg)),
as long as μ + δ_deg > 0 (tested for δ_deg = 0 and 0.05). A runaway would also need a separately
specified coupling that drives μ → 0 as A grows, with no aggregate degradation. A6 mentions that
growth-dilution feedback but does not model it, and this audit does not either.

**What in the old A6 text must change:**
- "two of them … or none": false. There is always ≥ 1 steady state.
- "Above L_crit the aggregate grows until chaperone is fully titrated": true up to titration. It
  then stops at A_hi. It does not diverge.
- "a large enough perturbation … produces the same runaway": no. It produces a jump to the
  finite high state (history dependence survives, as hysteresis).
- The bistability claim must add k_c M_tot > k_a + μ as a stated, unmeasured condition.

**Upstream inconsistencies the reduction inherits (not fixed here):**
- (A1.1) with constant β is itself impossible once βΣA > M_tot. Sequestration has to saturate in
  the conservation law, not only in A6.
- A6 writes M_free = M_tot − βA and drops the client-bound term ΣC that (A1.1) includes. I did not
  test whether keeping ΣC changes the fold condition.

## A7 Maintenance

**Moving boundary.** With φ = N_i − N_i^req(μ(x)), Nagumo's condition on the face φ = 0 is

    P_i − (μ + k_u,i + k_dN,i) N_i^req − (dN_i^req/dμ)·(dμ/dt) ≥ 0       [sympy, `phidot_face`]

A7.1 keeps only the first two terms. That is enough only if N^req does not depend on μ, or if μ is
stationary on the boundary. Counterexample (`movingQuotaCounterexample`): N^req = 1 + 2μ, and μ
relaxes from 0.5 to 1.0. A7.1 holds with zero slack at t = 0, yet dφ/dt = −1.0 and φ(0.1) = −0.10,
so the state leaves V immediately. Also, A7.1 is written "for every x ∈ ∂K". It should apply only at
points of K where N_i = N_i^req. At other boundary points of K the condition is K's own tangent
cone, and A7 leaves that unspecified. So A7.1 is a necessary condition on one face of K, not a
sufficient one for invariance.

**No native degradation.** (A1.3) has no N → D flux. Regulated proteolysis of folded proteins is
therefore absent. Adding k_dN changes the steady quota to P_i ≥ (μ + k_u + k_dN)N^req. It also
changes the μ → 0 limit in A3 to (k_u + k_dN)N^req(0), not k_u N^req, and the A3 limit should use
N^req(0) in any case.

**Viability boundary when N^req depends on μ.** This is an illustration under two assumptions, not
results: synthesis proportional to growth (s = s0 μ, as A2 states), productive fraction y, and
N^req = n0 + n1 μ. The quota then reads n1μ² + (n0 + n1k_u − y s0)μ + k_u n0 ≤ 0.
- With strict s ∝ μ and k_u n0 > 0, **μ near 0 is never viable** (lower root > 0). A6's "graceful
  failure … until the quota cannot be met even as μ → 0" assumes μ → 0 is the last viable point,
  but under A2's own scaling it is not. A basal synthesis term restores low-μ viability only if
  y·s_basal ≥ k_u n0.
- As y falls, the feasible interval [μ_lo, μ_hi] shrinks and **vanishes at a positive growth rate**
  μ* = √(k_u n0/n1), at y_crit = (√n0 + √(n1k_u))²/s0 [sympy]. In the example, the feasible
  interval is [0.055, 0.916] at 1.5 y_crit, [0.180, 0.278] at 1.01 y_crit, and empty at 0.99 y_crit.
  So even with β = 0, feasibility can disappear abruptly at a nonzero μ. This is a statement about
  the feasible set only. Whether the realized dynamics jumps depends on the μ(N) closure, which
  EQUATIONS.md does not specify.

## Limits of this audit
Only the reduced one-variable A6 model and toy A7 models were tested. Regulation (σ32 induction),
ΣC in the pool, seeded aggregation and an explicit μ(A) were not included. None of β, k_c, k_dN,
n1 or the fold condition has a measured value here.
