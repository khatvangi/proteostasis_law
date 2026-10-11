# Theorems and proofs

Statements, hypotheses and proofs for the spine in [MASTER_ARGUMENT.md](MASTER_ARGUMENT.md).
Notation and each equation's physical meaning are in [FORMAL_MODEL.md](FORMAL_MODEL.md).
Every claim marked **[test: name]** is checked by `test_laws.py`, and every quoted number comes
from `results.json`, which `python laws.py` writes. All toy parameters are dimensionless
placeholders. None is an E. coli value, and none was fitted to anything.

Reproduce:

    cd experiment_grounded_rebuild_20261010/laws
    python laws.py                        # writes results.json
    python -m unittest test_laws -v       # 32 tests, about 90 s

---

## T0. Conservation and positive invariance (reused, not new)

The state is (U_i, N_i, A_i, C_{m,i}, D_i) for each client i and machine m. The machine pool
is capped (`../theory/EQUATIONS.md` A1.1a):

    M_m^free = max(M_m^tot − Σ_i C_{m,i} − β_m Σ_i A_i, 0).

Client mass obeys

    d/dt (U_i + N_i + A_i + Σ_m C_{m,i} + D_i) = s_i − μ (U_i + N_i + A_i + Σ_m C_{m,i}).      (T0.1)

The nonnegative set with Σ_i C_{m,i} ≤ M_m^tot is forward invariant (Proposition A1 there,
proved by Nagumo's theorem). Everything below assumes these two facts.
**[test: every ODE in this folder carries a mass-balance counter; residuals ≤ 1.3·10⁻¹² in the race
model and ≤ 10⁻⁸ in the two-client model, and no pool goes below −10⁻¹⁰]**

---

## T1. Native replacement quota (necessary)

Write λ_i ≡ μ + k_u,i + k_dN,i. The native balance is

    dN_i/dt = P_i − λ_i N_i,    P_i = σ_i s_i + k_f,i U_i + Σ_m κ^N_{m,i} C_{m,i}.            (T1.1)

**T1(a), steady state.** A steady state with N_i ≥ N_i^req(μ) exists only if
P_i ≥ λ_i N_i^req(μ).

**T1(b), moving face.** Take a closed set K ⊆ V = {N_i ≥ N_i^req(μ)}. If K is forward invariant
and x ∈ ∂K with N_i = N_i^req(μ), then

    P_i − λ_i N_i^req(μ) − (dN_i^req/dμ) · μ̇ ≥ 0.                                           (T1.2)

*Proof.* (a) Set (T1.1) to zero and use N_i ≥ N_i^req. (b) φ = N_i − N_i^req(μ(t)) is C¹ along
solutions, and φ = 0 on the face. If φ̇ < 0 there, the solution leaves V, and with it K. □

Without the μ̇ term the condition is wrong whenever N^req depends on μ.
**[test: test_moving_quota_face_term_is_needed]** Take N^req = 1 + 2μ, with μ relaxing from 0.5 to 1
at rate 1 and k_u + k_dN = 0.1. The static condition holds with zero slack at P = 1.2, yet the gap
falls to −1.91. The face condition requires P = 3.3, and the gap then stays ≥ 0. This agrees with
`../theory/EQUATIONS.md` A7.

---

## T2. Occupancy and the productive-flux ceiling (necessary)

**T2(a), Little's law.** At steady state C_{m,i} = J_{m,i} τ_{m,i}, with
τ = 1/(κ^N + κ^U + κ^A + κ^D + μ). Hence Σ_i J_{m,i} τ_{m,i} ≤ M_m^tot − B_m (EQUATIONS A4).

**T2(b), ceiling at every instant.** Write κ^N_{m,i} = p^N_{m,i}/τ_{m,i}, where p^N is the
probability that a bound chain is released native. Then, along every solution,

    σ_i s_i  ≤  P_i(t)  ≤  σ_i s_i + k_f,i U_i(t) + Σ_m (p^N_{m,i}/τ_{m,i}) · (M_m^tot(t) − B_m(t) − Σ_{j≠i} C_{m,j}(t)).   (T2.1)

**T2(c), inventory bound.** Write Z_i = U_i + N_i + A_i + Σ_m C_{m,i}. If s_i ≤ s̄_i and μ ≥ μ_lo > 0,
then U_i(t) ≤ Z_i(t) ≤ Z_i(0) e^{−μ_lo t} + (s̄_i/μ_lo)(1 − e^{−μ_lo t}).

*Proof.* (b) follows from (T1.1), C ≥ 0 and (A1.1a). (c) By (T0.1), dZ/dt ≤ s̄ − μ_lo Z, then
comparison. □

**Corollaries.**
- **Co-translational floor.** If σ_i s_i ≥ λ_i N_i^req, client i's quota cannot fail whatever the
  machines do. Machine-limited failure therefore needs σ_i s_i < λ_i N_i^req.
  **[test: test_cotranslational_floor_excludes_failure]**
- **Obligate clients.** If k_f,i ≈ 0 under the condition, as for GroEL-obligate clients (Kerner
  2005; Fujiwara 2010, where only about 60% of the Class III set is truly obligate), the ceiling
  contains only measurable or boundable quantities: σ s, p^N/τ and M^tot(t).
- The bound concerns *productive* flux. Total synthesis s_i is not a bound on P_i once the
  inventory term k_f U is included, because a backlog of U built earlier can fold later.

---

## T3. Time-to-depletion window (bounds from bounds)

**Hypotheses on [0, T].** P(t) ∈ [P_lo, P_hi], λ(t) ∈ [λ_lo, λ_hi] with λ_lo > 0, and
N^req(t) ∈ [R_lo, R_hi]. N(0) = N_0 > R_hi. The failure time is t_f = inf{t : N(t) < N^req(t)}.

**T3(a), envelope.**

    N_0 e^{−λ_hi t} + (P_lo/λ_hi)(1 − e^{−λ_hi t})  ≤  N(t)  ≤  N_0 e^{−λ_lo t} + (P_hi/λ_lo)(1 − e^{−λ_lo t}).   (T3.1)

**T3(b), guaranteed failure.** If P_hi < λ_lo R_lo, then

    t_f ≤ t_up = λ_lo⁻¹ ln[(N_0 − P_hi/λ_lo) / (R_lo − P_hi/λ_lo)],                          (T3.2)

provided the hypotheses hold on [0, t_up].

**T3(c), guaranteed survival.** If P_lo ≥ λ_hi R_hi, then t_f = ∞ on [0, T]. Otherwise

    t_f ≥ t_lo = λ_hi⁻¹ ln[(N_0 − P_lo/λ_hi) / (R_hi − P_lo/λ_hi)].                          (T3.3)

**T3(d), transient overload.** Suppose an overload lasts T_ov, and the bounds describe the
overload. If T_ov < t_lo, the quota holds throughout it. If T_ov > t_up, the quota is violated
during it. Between the two, the bounds do not decide.

*Proof.* N ≥ 0 by T0, so −λN ≤ −λ_lo N, and dN/dt ≤ P_hi − λ_lo N. Let y solve y' = P_hi − λ_lo y
with y(0) = N_0. Then (N − y)' ≤ −λ_lo (N − y), and Grönwall gives N ≤ y. The lower bound is
symmetric, using λ ≤ λ_hi and P ≥ P_lo ≥ 0. When P_hi/λ_lo < R_lo < N_0, y decreases
monotonically and crosses R_lo once, at t_up. For t > t_up, N ≤ y < R_lo ≤ N^req(t). Likewise,
for t < t_lo, N ≥ y_lo > R_hi ≥ N^req(t). □
**[test: test_symbolic_solution_and_crossing (sympy); test_failure_time_inside_window_random (150
draws with P, λ and N^req oscillating inside their bounds, every simulated t_f ∈ [t_lo, t_up]);
test_upper_bound_is_attained_by_extreme_inputs (t_up is exact when P ≡ P_hi and λ ≡ λ_lo)]**

**What has to be measured.** Each bound needs only quantities that can be measured or bounded:

| quantity | what fixes it | status in E. coli |
|---|---|---|
| N_0 | native level before the overload (abundance × specific activity) | measurable per protein |
| N^req | lowest native level that supports the growth rate (promoter or CRISPRi titration) | measurable; few done |
| λ_lo | μ + k_dN (k_u ≥ 0 is dropped), from growth plus pulse-chase of the native protein | measurable |
| λ_hi | also needs an upper bound on k_u under the condition | rarely measured |
| P_hi | T2.1, needs σ s, p^N/τ and M^tot(t) | τ, p^N not measured in vivo |
| P_lo | σ s (co-translational floor); P_lo = 0 is always valid | σ unmeasured |

With only P_lo = 0 the survival bound reduces to t_lo ≥ λ_hi⁻¹ ln(N_0/R_hi), the pure-decay floor.
The guaranteed-failure bound has no machine-free form, because it needs P_hi.

---

## T4. Route-qualified linear feasibility, with a Farkas certificate

**Objects.** A route r belongs to one client i(r). It is a way a non-native chain, having left
U_i, can go on. Route data:

- a_{m,r} ≥ 0 [time]: expected time machine m is held per chain sent down r. It is
  (expected visits) × τ, and a serial route has a_{m,r} > 0 for several m at once.
- ρ_r [–]: probability that r ends native.
- b_r [–]: probability that r returns the chain to U. ρ_r + b_r ≤ 1.
- an optional ceiling x_r ≤ c_r [conc/time]. For a departure out of U with hazard k_r(x),
  c_r = k̄_r Ū, where k̄_r ≥ k_r and Ū bounds U.

Route flux x_r ≥ 0 [conc/time]. Pools M_m [conc]. Entry e_i = (1−σ_i)s_i + k_u,i N_i^req
[conc/time]. Deficit d_i = λ_i N_i^req − σ_i s_i [conc/time]. The system F is

    (cap)    Σ_r a_{m,r} x_r ≤ M_m^tot − B_m        for every machine m      [conc]
    (sup)    Σ_{r∈i} (1 − b_r) x_r ≤ e_i              for every client i       [conc/time]
    (dem)    Σ_{r∈i} ρ_r x_r ≥ d_i                    for every essential i    [conc/time]
    (ceil)   x_r ≤ c_r                                                         [conc/time]

**T4(a), necessity.** Assume the kinetic model (A1) has a steady state with N_i ≥ N_i^req for every
i ∈ E, and that some machine-free loss route exists (dilution of U at rate μ > 0 always is one).
Take as the route library the one-step routes "leave U_i by option o", evaluated at that steady
state. Then F is feasible.

*Proof.* At the steady state every hazard is constant. A chain leaving U_i by binding machine m
stays bound for an exponential time with mean τ_{m,i}, then is released native (p^N), back to U
(p^U), or lost. So a_{m,r} = τ, ρ = p^N, b = p^U. Spontaneous folding has a = 0, ρ = 1. Loss
options have ρ = 0. Put x_o = hazard_o · U_i^*.
(cap): Σ_o a x = Σ_i τ J = Σ_i C ≤ M − B, by Little and (A1.1a).
(sup): departures from U equal arrivals into U, so Σ(1−b)x = (1−σ)s + k_u N_i (the
disaggregation return is part of the aggregation route's b).
(dem): the native balance gives Σ ρx = λN − σs ≥ d_i.
(ceil): x_o = k_o(x^*) U^* ≤ k̄_o Ū.
This holds at N = N^*. To reduce to N^req, lower the route flux by a total (1−b)-weighted amount
k_u(N^* − N^req). The native loss is at most that amount, since ρ ≤ 1 − b, and it is at most
λ(N^* − N^req), since k_u ≤ λ. Capacity use only falls. □
**[test: test_realized_steady_state_is_lp_feasible, 200 random finite-pool instances]**

For serial **hand-off** routes (TF or DnaK → GroEL) the residence time in the first machine
depends on the free GroEL pool. The necessary LP then uses the pool-independent *lower* bound
a ≥ 1/(κ_tot + k_h M_G^tot + μ). Shrinking coefficients is a relaxation, so it stays necessary.

**T4(b), Farkas alternative.** Write F as Gx ≤ h, x ≥ 0. Exactly one of the following holds:
(i) F is feasible; (ii) there are multipliers y ≥ 0 with Gᵀy ≥ 0 and hᵀy < 0. In physical form,
(ii) asks for

- z_i ≥ 0 per essential client [dimensionless: value per native chain],
- u_i ≥ 0 per client [value per entering chain],
- w_m ≥ 0 per machine [value per (machine · unit time), units 1/time],
- c'_r ≥ 0 per ceiling,

such that

    for every route r of client i:  Σ_m w_m a_{m,r} + u_i (1 − b_r) + c'_r ≥ z_i ρ_r          (T4.1)
    Σ_i z_i d_i  >  Σ_m w_m (M_m − B_m) + Σ_i u_i e_i + Σ_r c'_r c_r.                         (T4.2)

Every term of (T4.2) has units of conc/time. Every term of (T4.1) is dimensionless, because
w_m [1/time] × a_{m,r} [time]. (T4.1) says no route turns resources into native chains more cheaply
than z values them. (T4.2) says the valued demand exceeds the valued resources. Together they
prove infeasibility. *Proof:* the Farkas lemma for {Gx ≤ h, x ≥ 0}. □
**[test: test_farkas_alternative_random. 300 random serial and parallel instances; exactly one
alternative holds each time; feasible points satisfy every row, and certificates satisfy (T4.1–2)
independently of the solver. Presolve is disabled because of the earlier HiGHS false-infeasible
episode. test_time_unit_invariance: rescaling time leaves feasibility and margin unchanged]**

**T4(c), why a "Hall" test is wrong here.** A Hall-type test compares the demand of any client set
S with the capacity of the machines S touches, Σ_{m∈N(S)} M_m/a_m. That treats the machines in
N(S) as substitutes. A serial route needs all of them at once.
- *Serial counterexample* **[test: test_serial_route_breaks_hall_union]**: one client, one
  route A→B, a_A = a_B = 1, M_A = M_B = 1, d = 1.5. The union test passes (1.5 ≤ 2). F is infeasible,
  and an exact rational certificate is z = 1, w_A = 1, w_B = 0, u = 0. (T4.1) gives 1 − 1 = 0 ≥ 0, and
  (T4.2) gives 3/2 − 1 = 1/2 > 0. The solver's certificate puts the weight on B, which is equally
  valid.
- *Parallel counterexample to the opposite rule* **[test: test_parallel_route_breaks_bottleneck_min]**:
  "each client is limited by its slowest machine" rejects two parallel routes (A or B, d = 1.5),
  which F accepts (x_A = 1, x_B = 0.5).
- *Where a Hall/Gale condition is exact* **[test: test_gale_cut_condition_exact_for_single_machine_routes,
  200 instances]**: every route uses exactly one machine, ρ = 1, b = 0, and the service time
  depends only on the machine. This is a bipartite supply–demand flow, and Gale's theorem applies.
  Serial routes, yields below 1, returns, and client-specific τ all break it.

**T4(d), global pooling is a strict relaxation.** Summing the (cap) rows gives a valid inequality,
so qualified feasibility implies pooled feasibility. The converse is false. In the instance, client
g can use only machine G and client k only K, with M_G = 0.5, a_G = 2, d_g = 1. The pooled system
is feasible; the qualified one is infeasible. **[test: test_pooling_is_a_strict_relaxation, plus
200 random checks of the implication]**

**T4(e), max-load margin and its sensitivity.** Θ* = max{θ : F with d replaced by θd is
feasible}, and F is feasible iff Θ* ≥ 1. Θ* is nondecreasing in every M_m.
**[test: test_margin_monotone_in_every_pool]** At a nondegenerate optimum, LP sensitivity gives
∂Θ*/∂M_m = w_m* and ∂Θ*/∂a_{m,r} = −w_m* x_r*, where w_m* is the capacity dual. A machine with
slack has w_m* = 0. **[test: test_margin_sensitivity_binding_vs_slack. Finite difference −0.75000
against −w_G x_r = −0.75 on the binding machine, and exactly 0 on the slack one]**

**T4(f), what F does not do.** F lets the split x be chosen. The cell does not choose it: rate
competition sets it (T5). F is therefore only necessary. Without ceilings it says nothing about
non-obligate clients: if a route with a = 0 and ρ = 1 (spontaneous folding) is qualified, every
demand up to the entry flux is "feasible". Capacity statements bite only through the ceilings
k̄ Ū, or for clients with no machine-free productive route.

---

## T5. Separation 1: capacity is monotone, realized yield is not

**Setting.** One client, entry flux e into U, first-order loss rate k_l from U (unrecovered
aggregation, protease, dilution), spontaneous folding k_f, and machines m with binding hazard
θ_m = k_on,m M_m^free and one-step outcomes p^N_m, p^U_m. The per-entry native yield is

    Y = (k_f + Σ_m θ_m p^N_m) / (k_f + k_l + Σ_m θ_m (1 − p^U_m)).                            (T5.1)

Define the **committed success** of route m, ρ̂_m = p^N_m/(1 − p^U_m). It is the probability that
a chain bound to m ends native rather than lost, once returns are discounted.

**Theorem S1.**
(a) ∂Y/∂θ_m = [(1 − p^U_m)/Q] · (ρ̂_m − Y), with Q the denominator of (T5.1). Raising route m's
    hazard raises realized yield iff its committed success exceeds the current yield.
(b) With one machine of finite pool M^tot and no sequestration, dM^free/dM^tot > 0 at the unique
    steady state. Hence sign(dY/dM^tot) = sign(ρ̂ − Y).
(c) Θ*(M) is nondecreasing in every M_m (T4e).
(d) So whenever ρ̂_m < Y, adding capacity of machine m strictly lowers realized yield while the
    capacity margin cannot fall.
(e) Quota failure is possible at arbitrarily low occupancy. For every ε > 0 there is an instance
    with every machine below occupancy ε and Y < ε.

*Proof.* (a) Differentiate (T5.1): ∂Y/∂θ_m = [p^N_m Q − (1 − p^U_m) P]/Q² with P the numerator,
which is the stated form. (b) The steady state solves e = (k_f + k_l)U + (1 − p^U) C/τ with
C = M^tot − M^free and C = k_on M^free U τ. Suppose M^tot rises and M^free does not. Then C rises,
so U must fall to keep e fixed. But M^free = M^tot/(1 + k_on τ U) then rises, a contradiction.
Since θ = k_on M^free and Y depends on M^tot only through θ, (a) gives the sign. (c) The feasible
set of F only grows with M. (e) Take k_f = 1, k_l = 0, a productive machine with ρ̂ = 1, and a
disposal machine with ρ̂ = ε³ and capture rate ε⁻² per unit pool, both pools ε⁻¹ (laws.py
`findLowOccupancyFailure`). □
**[test: test_sign_law_symbolic (sympy); test_sign_law_transfers_to_finite_pool_single_machine
(500 draws, every evaluable sign agrees); test_s1_capacity_up_realized_down;
test_quota_failure_at_vanishing_occupancy]**

**Instance** (`results.json` → S1_comparison). G is productive (p^N = 0.9) and H is
disposal-dominated (ρ̂_H = 0.1). As M_H goes 0 → 0.5 → 2 → 5, realized Y goes
0.898 → 0.692 → 0.392 → 0.234, while Θ* stays at 1.632. With ε = 0.01: Y = 1.0·10⁻⁴ at maximum
occupancy 1.0·10⁻⁶.

**Fails with two pools.** The transfer in (b) does **not** extend to two or more finite pools. In
400 random two-machine draws, sign(dY/dM_1^tot) differs from sign(ρ̂_1 − Y) in **21**
(`results.json` → S1_two_machine_sign_exceptions). The reason is that adding pool 1 moves U,
which moves the *other* machine's free pool. Rule (a), at fixed free hazards, still holds exactly.
Rule (b) is a single-pool statement. With shared pools the realized sign has to be computed from
the full steady state, as in T9(c).

**What it means experimentally.** A growth or yield gain on deleting a route (ΔhtpG in ΔdnaKJ at
37/39 °C; ΔhslV in ΔdnaKJ at 39 °C; Fauvet 2021) cannot be explained by capacity, because removing
a route never raises Θ*. Within this model it says that in that background the route's committed
success is below the background's realized yield. Deleting HslUV raises the HtpG route's ρ̂, since
disposal becomes return (p^U rises, p^N/(1 − p^U) rises), so it acts in the same direction as
deleting HtpG. An optimizing model (FoldME-type) cannot produce such a gain without an added cost
term.

---

## T6. Separation 2: the induction race

**Hypotheses.**
(H1) The client's productive flux satisfies P(t) ≤ p_0 + κ C_m(t). Here p_0 is the machine-free
     part (σ s + k_f Ū, which is σ s for an obligate client), κ = p^N_m/τ_m, and C_m(t) is the
     client bound to m.
(H2) C_m ≤ M_m(t) (T0).
(H3) Machine synthesis is bounded but arbitrarily regulated: dM_m/dt = b_0 + b_1 u(t) − γ M_m, with
     u(t) ∈ [0, 1] measurable and otherwise arbitrary. This includes every σ32 law, every gain and
     every delay. γ is dilution plus machine turnover.
(H4) After the shock at t = 0, λ is constant (or bounded below by λ) and N^req is constant.

**Theorem S2.** Write M_c = (b_0 + b_1)/γ and M_0 = M_m(0). For λ ≠ γ:

(a) M_m(t) ≤ M̄(t) = M_c + (M_0 − M_c) e^{−γt}.

(b) N(t) ≤ N̄(t) = N_0 e^{−λt} + [(p_0 + κM_c)/λ](1 − e^{−λt}) + κ(M_0 − M_c)(e^{−γt} − e^{−λt})/(λ − γ).   (T6.1)

(c) If min_{t≥0} N̄(t) < N^req, the quota is violated at some t > 0 under **every** regulation law
    u(·). This holds even when the static, full-induction condition p_0 + κM_c ≥ λN^req holds.

(d) Hence neither steady-state capacity feasibility (T1, T4) nor a locally stable equilibrium inside
    V is sufficient for maintenance through a transient. Sufficiency needs a forward-invariant set
    that contains the post-shock state (T7).

*Proof.* (a) Comparison with u ≡ 1. (b) dN/dt ≤ p_0 + κ M̄(t) − λN by (H1–H4) and (a). Integrate
the linear comparison equation. sympy confirms that (T6.1) solves it with N̄(0) = N_0. (c) N ≤ N̄
(comparison), and N̄ < N^req at some t. (d) The instance below has both properties and still
violates. □
**[test: test_bound_symbolic; test_machine_ceiling_any_regulation (50 random piecewise-constant
controls, worst excess of M over M̄ is 1.8·10⁻¹⁵); test_trajectory_below_bound_and_conserved;
test_static_feasible_yet_guaranteed_transient_violation; test_post_shock_equilibrium_locally_stable;
test_violation_survives_maximal_instant_induction]**

**Instance** (`results.json` → S2_race). This is a mechanistic closed loop: an obligate client
(k_f = 0, σ = 0.2), one productive machine with σ32-type synthesis repressed by free machine,
and a heat-like step in native unfolding k_u from 0.02 to 1.0, so λ = 1.085, against γ = μ = 0.035.
N_0 = 28.57 and N^req = 12. The static full-induction steady state of the bound is
(p_0 + κM_c)/λ = 36.4 > 12. The closed loop settles at N^* = 18.99 > 12, a locally stable
equilibrium (largest real eigenvalue −0.065). The bound reaches 6.60 at t = 3.0, so violation is
guaranteed by t = 1.04. The simulated N first falls below 12 at t = 1.03, with minimum 6.54. With
forced maximal induction from t = 0 (u ≡ 1, zero delay), N still falls to 6.60. Mass-balance error
is 8·10⁻¹³.

**Physical reading.** What decides the race is how fast regulated synthesis can raise the
machine pool (b_0 + b_1 − γM_0 early on, and 1/γ for the approach to M_c), compared with how fast
native protein is lost (λ). Feedback gain and delay cannot beat (T6.1). They can only lose
further ground to it. For a heat step, λ jumps through k_u while γ stays near μ. So the bound
predicts transient loss of obligate-client function after a shift, even when the induced steady
state is fully adequate. The size of that loss is set by measurable rates, and the bound is
falsifiable (PREDICTIONS P2).

---

## T7. Maintenance: necessary versus sufficient

Definitions as in `../theory/EQUATIONS.md` A7. V = {N_i ≥ N_i^req(μ), i ∈ E}. W is the
disturbance class. "Maintainable from x_0" means there is a closed K ⊆ V containing x_0 that is
forward invariant for every admissible w(·). With a regulator u (σ32), "controlled-maintainable"
means some admissible law makes K invariant: x_0 lies in the robust viability kernel of V.

**Necessary conditions.** Each must hold. None is sufficient.
- N1: static quota P ≥ λN^req at any equilibrium in V (T1a).
- N2: occupancy ΣJτ ≤ M − B (T2a).
- N3: route feasibility F ≠ ∅, equivalently Θ* ≥ 1, with no Farkas certificate (T4).
- N4: the face condition (T1.2) on every active quota face of K.
- N5: for every admissible disturbance, the race bound N̄_w(t) ≥ N^req for all t (T6). Otherwise
  even the best regulation fails.

**Sufficient conditions, and what each needs.**
- S-floor: P_lo ≥ λ_hi R_hi with P_lo = σs, or more generally T3(c). This needs no knowledge of the
  machines, and it is rarely met by a chaperone-dependent client.
- S-overload: T_ov < t_lo (T3d). This needs a lower bound on productive flux over the overload.
  Any machine contribution to P_lo needs a lower bound on bound client C(t), and so a lower bound on
  U(t) and M^free(t). That is a statement about an **invariant set** of the realized kinetics, not
  about capacity.
- S-Nagumo: a candidate K with the tangent condition on every face (quota faces with (T1.2) for the
  worst w, plus faces bounding U, C, M from below). This is the only general sufficient route, and
  it **requires** an explicit invariant set.
- S-basin: a locally stable equilibrium in int V whose basin contains every post-disturbance state.
  This **requires** an estimate of the basin, such as a Lyapunov sublevel set inside V. T6 shows the
  equilibrium alone is not enough: the S2 instance's equilibrium is stable and inside V, yet the
  trajectory from the pre-shock state leaves V.

So no condition on capacity (pools, occupancy, LP margin) is sufficient. Every sufficient
condition is a statement about the realized dynamics: an invariant set, a basin, or a trajectory
bound.

---

## T8. σ32: one input, several machines, and what free-machine sensing cannot see

**T8(a), the capacity set is up-closed.** If F is feasible at pools M, it is feasible at every M'
≥ M. Suppose a single regulator drives every machine monotonically, M_m(u) nondecreasing in u ∈ [0, 1].
Then some u gives F ≠ ∅ iff u = 1 does. **[test: test_capacity_set_up_closed]** So "one input cannot
set two capacities" (feedback/MINIMAL_DYNAMICS §6) is **not** a capacity restriction: for
feasibility, more of everything is never worse. It becomes one only through realized kinetics, as
in (b), or through costs (T10).

**T8(b), realized yield can be non-monotone in regulator output.** The σ32 regulon raises proteases
(ClpXP, Lon, HslUV, FtsH) as well as chaperones (Zhang 2014; Fauvet 2021 for ΔdnaKJ). By S1, adding
a disposal-dominated machine can lower Y. In the instance (`results.json` → L6_coregulated),
M_G = 0.1 + 2u and M_H = 2u: Y(0) = 0.231, Y(0.31) = 0.316, Y(1) = 0.255.
**[test: test_realized_yield_has_interior_optimum]**

**T8(c), free-machine sensing is blind to routing failure.** A regulator repressed by free machine
(σ32 by free DnaK/J and GroEL/S; Tomoyasu 1998, Guisbert 2004) responds to occupancy. By S1(e),
quota failure can coexist with free fractions above 1 − ε. Such a failure produces near-basal
induction. This is immediate from S1(e) and the sensing assumption, which is a property of the
model loop, not a measured gain.

**Reused, not re-derived here** (feedback_review/REVISIONS_REQUIRED.md): the linearised delayed
loop is stable at every delay iff G ≤ 1, and G > 1 gives a finite τ*. "C ≤ C_max" holds only for
C(0) ≤ C_max at constant μ. The collapsed growth-mode state depends on assuming chaperone synthesis
∝ μ.

---

## T9. Codon-conditioned perturbations

Synonymous codon choice z_i enters through σ_i (the co-translational native fraction) and through
the birth distribution π⁰ over non-native states (free U, TF-held, DnaK-held), for example by
elongation pausing. Misincorporation frequency ε creates variant species with their own route data.
These are different parameters, and proteomic error tables constrain only ε.

**T9(a), single client, fixed pools.** Y_tot = σ + (1 − σ) Σ_k π⁰_k Y_k, where Y_k is the realized
yield starting from state k. Then ∂Y_tot/∂σ = 1 − Σ π⁰_k Y_k ≥ 0. Moving birth mass from k to k'
changes Y_tot by (1 − σ)(Y_k' − Y_k)Δπ, which has the sign of Y_k' − Y_k. **[test:
test_cotranslational_symbolic]** "Synonymous change X lowers chaperone demand" is therefore a
statement about σ or π⁰, each of which must be measured. It is not a consequence of lower error.

**T9(b), demand effects matter only on binding machines.** By T4(e), a change Δa_{m,r} moves the
margin by −w_m* x_r* Δa to first order. That is zero for any machine with slack. A codon effect on
growth should therefore be **conditional on the capacity regime** of the machine it loads.

**T9(c), cross-client sign under shared pools.** Raising σ of client 1 can lower client 2's yield.
With σ_1 going 0 → 0.3, Y_1 rises 0.164 → 0.405 and Y_2 falls 0.647 → 0.626. The mechanism is that
the free disposal pool H rises 0.924 → 1.120, and ρ̂_{H,2} = 0.05 < Y_2, so S1 applies to client 2.
The fixed point agrees with direct ODE integration to 10⁻⁶. **[test:
test_two_client_cross_effect_and_ode_agreement]** A codon-optimized gene can therefore lower
another protein's yield. This is a two-client consequence of S1, not a general rule.

**T9(d), mistranslation load bound.** For a protein carrying n_c codons of class c with per-codon
error ε_c, the variant fraction is 1 − Π(1 − ε_c)^{n_c} ≤ Σ_c n_c ε_c. The extra machine-time load
is ≤ s Σ_c n_c ε_c · a_var. To move Θ* it must be comparable to the binding machine's slack.
**[test: test_mistranslation_first_order_bound]** The only synonymous per-codon contrast that
survived the calibrated test is serine AGC: pooled detected rate 3.1·10⁻⁴ against
1.8–2.6·10⁻⁵ for TCT/TCC, 96% Ser→Asn (`../../paper_reassessment_20261010/codon_bridge/CODON_RESULT.md`).
That rate is detection-conditional, and there is no folding readout. The bound says an AGC-derived
variant load matters only if the Asn variant's machine time per chain exceeds the native's by a
factor of order slack/(s · n_AGC · 3.1·10⁻⁴). No measurement of that factor exists.

---

## T10. Energy: an optional constraint, not an axiom

An energy row Σ_r E_r x_r ≤ E_max can be added. E_r includes ATP per machine cycle and, for routes
that lose chains, the resynthesis cost (1 − ρ_r)·c_syn. Fauvet 2021 cite an estimate that
resynthesis costs 100–1000× more ATP than repair; that number is not used here. With the row
absent, the only objective is the margin Θ*. With it present, (Θ*, E_max) has a front.
**[test: test_front_exists_only_with_energy_row]** In the instance, Θ* = 0.30, 0.60, 0.87, 1.14,
1.41, 1.74, 1.74 at E_max = 0.5, 1, 2, 3, 4, 6, 10, against 1.74 without the row. Multiobjective
language is legitimate only once such a row, or another objective, is stated explicitly.
