# Equations and proofs (appendix to SCIENTIFIC_THEORY.md)

Every sign, threshold and counterexample marked **[checked]** is recomputed by `checks.py`
(22/22 pass, output in `checks.json`). Parameter values in `checks.py` are illustrative and
dimensionless. None is an E. coli estimate.

Sections A1, A3, A6 and A7 were corrected on 2026-10-10 after an independent audit
(`../independent_review/AUDIT_A6_A7.md`). What was wrong, and why, is in
[CORRECTIONS.md](CORRECTIONS.md).

---

## A1. State variables, conservation and positivity

For each client species i (a protein, or a class of proteins treated as one):

| symbol | meaning |
|---|---|
| U_i | soluble non-native chains: released nascent chains not yet folded, misfolded, unfolded natives |
| N_i | native, functional chains |
| A_i | chains in aggregates |
| C_{m,i} | chains bound to machine m as clients (one variable per machine, per client) |
| D_i | cumulative chains removed by degradation (a counter, not a pool) |

For each machine m (TF, DnaK/DnaJ/GrpE, GroEL/GroES, HtpG, ClpB, Lon, ClpXP, HslUV, …):

    M_m^tot = M_m^free + Σ_i C_{m,i} + B_m                                         (A1.1)

B_m is the number of machine units held on aggregates. Sequestration is capped by the pool:

    B_m = min( β_m Σ_i A_i ,  M_m^tot − Σ_i C_{m,i} ),
    so  M_m^free = max( M_m^tot − Σ_i C_{m,i} − β_m Σ_i A_i , 0 ).                   (A1.1a)

β_m is a binding *capacity*: the number of units one aggregated chain can hold when chaperone is
in excess. It is not the number it holds. Once the aggregate sites exceed what is left after
client binding, B_m equals that remainder and further aggregate holds nothing more. β_m = 0 means
no sequestration. (A1.1a) is the fast, tight-binding limit, with client-bound units not displaced
by aggregate. Both are assumptions. A slower mass-action version, with B_m as a state and its own
on/off rates, would also keep 0 ≤ B_m ≤ M_m^tot − Σ_i C_{m,i}, but its overload structure has not
been derived here.

The earlier form of (A1.1), M^tot = M^free + ΣC + βΣA with constant β, is impossible once
βΣA > M^tot − ΣC: it forces M^free < 0. That error produced the false "runaway" in the earlier A6
(CORRECTIONS.md, E1).

Cochaperones obey their own version of (A1.1). An effective rate that needs a cochaperone (GroES
for GroEL, DnaJ/GrpE for DnaK) is a function of the cochaperone's free pool too.

Process fluxes (per volume per time, all ≥ 0):

| process | from → to | flux |
|---|---|---|
| birth (synthesis) | – → chains | s_i(t) = elongation-limited synthesis rate |
| co-translational native production | birth → N | σ_i s_i |
| post-translational entry | birth → U | (1−σ_i) s_i |
| spontaneous folding | U → N | k_f,i U_i |
| misfolding / unfolding of native | N → U | k_u,i N_i |
| **degradation of native protein** | N → D | k_dN,i N_i (regulated proteolysis of folded chains) |
| machine binding | U → C_m | r^on_{m,i} M_m^free U_i |
| productive release | C_m → N | κ^N_{m,i} C_{m,i} |
| unproductive release | C_m → U | κ^U_{m,i} C_{m,i} |
| release into aggregate | C_m → A | κ^A_{m,i} C_{m,i} |
| hand-off to degradation | C_m → D | κ^D_{m,i} C_{m,i} |
| aggregation | U → A | k_a,i U_i (or a nucleated form; see A6) |
| disaggregation | A → U | v^dis_i(A, M_ClpB^free, M_DnaK^free), with v^dis_i = 0 when A_i = 0 |
| **degradation of aggregate** | A → D | k_dA,i A_i |
| protease degradation of U | U → D | Σ_p r^on_{p,i} M_p^free U_i (proteases are machines p with κ^D = 1) |
| dilution by growth | every pool | μ × pool |

Balance equations:

    dU_i/dt = (1−σ_i)s_i + k_u,i N_i + Σ_m κ^U_{m,i} C_{m,i} + v^dis_i
              − [k_f,i + k_a,i + Σ_m r^on_{m,i} M_m^free + μ] U_i                     (A1.2)
    dN_i/dt = σ_i s_i + k_f,i U_i + Σ_m κ^N_{m,i} C_{m,i} − (k_u,i + k_dN,i + μ) N_i   (A1.3)
    dA_i/dt = k_a,i U_i + Σ_m κ^A_{m,i} C_{m,i} − v^dis_i − (k_dA,i + μ) A_i           (A1.4)
    dC_{m,i}/dt = r^on_{m,i} M_m^free U_i − (κ^N + κ^U + κ^A + κ^D + μ)_{m,i} C_{m,i}    (A1.5)
    dD_i/dt = Σ_m κ^D_{m,i} C_{m,i} + k_dN,i N_i + k_dA,i A_i                          (A1.6)

Summing (A1.2–A1.6) gives the conservation law

    d/dt (U_i + N_i + A_i + Σ_m C_{m,i} + D_i) = s_i − μ (U_i + N_i + A_i + Σ_m C_{m,i}).   (A1.7)

Every synthesized chain is either diluted in a pool or counted in D_i, whichever pool it was
degraded from. The bound units B_m are machine, not client, so they do not enter (A1.7). When
aggregate is removed, (A1.1a) releases the units it held at once.
**[checked: a1_conservation_and_positivity. One client, one chaperone, first-order protease,
native and aggregate degradation, cap active (M^tot − C − βA = −0.11 at the end). Drift of
(U+N+A+C+D) − ∫(s − μ·pools) dt stays below 2·10⁻¹³.]**

**Proposition A1 (exact positive invariance).** Assume s_i ≥ 0, all rate constants ≥ 0, every
v^dis_i ≥ 0 and locally Lipschitz with v^dis_i = 0 at A_i = 0, and r^on, κ locally Lipschitz in the
state. Then the closed set

    Ω = { all U, N, A, C, D ≥ 0  and  Σ_i C_{m,i} ≤ M_m^tot for every m }

is forward invariant. Along every solution M_m^free ≥ 0 and 0 ≤ B_m ≤ M_m^tot − Σ_i C_{m,i}.

*Proof.* The right side is locally Lipschitz, since max(·, 0) is Lipschitz, so solutions are unique.
Ω is a polyhedron. Its tangent cone at a point is the set of directions that keep every active
inequality satisfied, so it suffices that on each face the field does not point outward. At
U_i = 0, every negative term of (A1.2) is proportional to U_i and vanishes, and the rest are ≥ 0.
The same holds for N_i in (A1.3), for A_i in (A1.4) (v^dis_i vanishes), and for C_{m,i} in (A1.5)
(binding is ≥ 0). D_i never decreases. On the face Σ_i C_{m,i} = M_m^tot, (A1.1a) gives
M_m^free = 0, so binding stops and d(Σ_i C_{m,i})/dt = −Σ_i(κ + μ)C_{m,i} ≤ 0. Nagumo's theorem
then gives forward invariance. The bounds on M^free and B follow from (A1.1a). □
**[checked: a1_conservation_and_positivity, a1_positivity_from_boundary. Starting from all pools at
zero except a large aggregate, with no synthesis, the minimum over pools and times is −2·10⁻¹³,
which is integrator round-off. The tolerance is −10⁻¹⁰.]**

The cap is what makes this true. Put the uncapped M^free = M^tot − C − βA into the same instance
and binding turns negative. C then falls to −0.12 and both checks fail (CORRECTIONS.md, E1).

## A2. State-dependent effective rates

Every rate is a function, not a constant:

    r_{k,r} = r_{k,r}(T, c, z_k, C_free)

- T is temperature. Folding, unfolding, aggregation and binding rates all have different
  temperature dependences. No monotone "temperature = load multiplier" is assumed; the
  Kanemori plating table already rules that out (EXPERIMENT_CONSTRAINTS row K99).
- c is the cochaperone and nucleotide state (GroES, DnaJ/GrpE, ATP/ADP).
- z_k is client-intrinsic. It includes the sequence and **the elongation-time profile set by
  synonymous codon choice**.
- C_free is the vector of free machine pools from (A1.1, A1.1a).

The synonymous-codon variable z_k enters only through σ_k and the partition among
k_f, k_a and r^on_{m,k} at the moment of release. That is a timing effect. It is a different
parameter from the substitution frequency ε_k that mass spectrometry measures. ε_k creates a
separate variant species k′ with its own rates. The Landerer and Stikeleather tables constrain ε
(standing frequency, detection-conditional). They contain no information about σ or the partition.

Synthesis rate s_k enters (A1.2) linearly as load. Write it as s_k = s_k⁰ + s_k¹ μ, a basal part
plus a growth-proportional part, so faster growth raises load and dilution together. Strict
proportionality (s⁰ = 0) has a consequence that is easy to miss: at slow growth the native quota
cannot be met (A3.2, A7).

## A3. The productive essential quota (necessary condition 1)

Let E be the set of essential species. Each needs a standing native level N_i ≥ N_i^req(μ) to
support growth rate μ. At steady state, with μ fixed, (A1.3) gives

    P_i ≡ σ_i s_i + k_f,i U_i + Σ_m κ^N_{m,i} C_{m,i} = (μ + k_u,i + k_dN,i) N_i.

So a steady state with N_i ≥ N_i^req requires the **native quota**

    P_i ≥ (μ + k_u,i + k_dN,i) N_i^req(μ)      for all i ∈ E.                         (A3.1)

The quota has three parts: dilution μN^req, repair of spontaneous unfolding k_u N^req, and
replacement of regulated native proteolysis k_dN N^req. The earlier form left out the third.
**[checked: a3_native_quota_includes_degradation. At the steady state of the A1 instance,
P = 1.5680 = (μ + k_u + k_dN)N. The old (μ + k_u)N = 1.2763 is short by exactly the
native-degradation flux.]** (A3.1) is a steady-state condition. When μ moves, the relevant
condition is the moving-boundary condition (A7.1).

Consequences:
- **The quota does not vanish at slow growth.** As μ → 0 it tends to (k_u,i + k_dN,i) N_i^req(0).
  Note the argument N^req(0), not N^req(μ). This matches the Santra 2017 model statement that very
  slow growth still needs chaperones. Here it follows from (A1.3) without fitting anything.
- **Reaching slow growth is not automatic.** Write synthesis as s = s⁰ + s¹μ (basal plus
  growth-proportional, A2), with productive fraction y (P = y s), N^req = n0 + n1 μ and
  k = k_u + k_dN. The quota becomes

      y (s⁰ + s¹ μ) ≥ (μ + k)(n0 + n1 μ).                                             (A3.2)

  If synthesis is strictly proportional to growth (s⁰ = 0), the left side → 0 as μ → 0 while the
  right side → k n0 > 0. So growth rates near zero are **never** feasible. Low-μ feasibility needs
  y s⁰ ≥ k n0. With s⁰ = 0 the feasible growth rates form an interval [μ_lo, μ_hi]. The interval
  shrinks as y falls and vanishes at a **positive** growth rate

      μ* = √(k n0 / n1),   at   y_crit = (√n0 + √(n1 k))² / s¹.

  **[checked: a7_feasible_growth_interval, symbolic and numeric. Example: s¹ = 10, n0 = 1, n1 = 2,
  k_u = 0.1, k_dN = 0. The interval is [0.055, 0.916] at 1.5 y_crit, [0.180, 0.278] at 1.01 y_crit,
  and empty at 0.99 y_crit. These are the auditor's numbers. Adding k_dN = 0.05 moves μ* from 0.224
  to 0.274. With k = 0.15 and s⁰ at 0.9 / 1.1 × k n0 / y, the low end is 0.0066 / 0.]** This describes
  the feasible set only. Whether the realized dynamics jump depends on the closure μ(x), which is
  not specified here.
- The quota bounds productive flux, not total synthesis. Synthesis minus productive flux must
  leave through aggregation, degradation of non-native, native or aggregated chains, or dilution
  (A1.7).

## A4. Machine occupancy (necessary condition 2): Little's law

At steady state, (A1.5) gives C_{m,i} = J_{m,i} τ_{m,i}, where J_{m,i} = r^on M^free U_i is the binding
flux and τ_{m,i} = 1/(κ^N+κ^U+κ^A+κ^D+μ)_{m,i} is the mean residence (service) time. Substituting
into (A1.1), with B_m the capped aggregate-held units of (A1.1a):

    Σ_i J_{m,i} τ_{m,i} = M_m^tot − M_m^free − B_m ≤ M_m^tot.                     (A4.1)

This is an identity of the steady state, not an assumption **[checked: littles_law_occupancy,
ΣJτ = 0.65249218 against 0.65249218 bound]**. Repeated binding cycles count as separate visits
(J counts visits, not clients). What is and is not measured:

- M_m^tot is measured (Schmidt 2016, BW25113, glucose, exponential).
- J_{m,i} and τ_{m,i} are **not measured in vivo for any E. coli machine–client pair**. Kerner's
  ">75% of GroEL capacity" is a reported occupancy share, an estimate of the left side of (A4.1)
  divided by M^tot for one client class. It is not J and not τ.

So (A4.1) is a constraint with a known right side and an unknown left side. It cannot be turned
into a number for the margin without in-vivo turnover data.

## A5. Feasible region versus realized routing. Necessary ≠ sufficient.

**Capacity-feasible region.** For fixed condition θ = (T, μ, s, z), with τ evaluated at θ:

    F(θ) = { J ≥ 0 : flux balance (A1.7) at steady state,
                      quota (A3.1) for i ∈ E,
                      occupancy (A4.1) for every m,
                      each J_{m,i} only on qualified routes (κ^N_{m,i} > 0 for productive use) }.

**Theorem 1 (necessity).** If (A1.2–A1.6) has a steady state with N_i ≥ N_i^req for all i ∈ E, then
F(θ) ≠ ∅.
*Proof.* The steady-state fluxes satisfy A3 and A4 by construction. □

**Realized routing.** The kinetic system does not choose J. From state U_i the fraction routed to
exit r is

    π_{i,r} = rate_r / Σ_{r′} rate_{r′},  with rates k_f, k_a, r^on_{m,i} M_m^free, protease rates, μ.

The realized point J^real(θ) is fixed by rate competition and free pools. It is not an optimum
over F.

**Theorem 2 (non-sufficiency).** F(θ) ≠ ∅ does not imply J^real(θ) satisfies (A3.1).
*Counterexample* **[checked: feasible_not_sufficient]**. One client class, s = 1, quota 0.8.
Capacity is ample (GroEL τ = 2, pool 5; HtpG τ = 1, pool 5). The best routing over F sends
everything productive and reaches native flux 1.0. With realized binding rates GroEL 1.0,
HtpG 6.0, spontaneous 0.6, loss 0.2, and HtpG disposal probability q = 0.9, realized native
flux = 0.28. Neither machine is more than 16% occupied. The quota fails with no capacity
shortage at all. □

**Corollary (a deletion can only help kinetically)** **[checked:
deletion_can_help_only_kinetically]**. Removing a route can only shrink F, so the optimum over F
cannot rise (1.0 → 1.0). In the same instance, removing HtpG raises realized native flux from 0.28
to 0.89. Any model in which the cell attains the best point of F therefore cannot reproduce an
observed growth improvement on deleting a chaperone (Fauvet 2021, ΔhtpG in ΔdnaKJ), unless it
adds an explicit cost term to the objective. A kinetic-routing model reproduces it with no extra
term.

A third gap between feasibility and maintenance comes from the three-state window of A6, when that
window exists (only if β > 0 and, in the reduced model, k_c M^tot > k_a + μ). There F can be
non-empty and a good steady state can exist while the cell sits at the finite high-aggregate state.

## A6. Overload in a reduced one-variable model: smooth response, kink, or a bounded hysteresis loop

**The reduction, and what it leaves out.** Take one aggregation-prone load L (the rate at which
non-native chains are produced) and one clearance chaperone with pool M^tot, capture rate k_c per
free unit, and aggregate capacity β. Hold U in quasi-steady state:

    U = L / (k_c M^free + k_a + μ),     dA/dt = k_a U − (μ + δ) A.

Three approximations are made, and none is tested here.
1. **The bound-client term is omitted.** The pool is M^free = M^tot − B, not M^tot − ΣC − B as in
   (A1.1). Keeping ΣC lowers the free pool at a given A and changes the fold condition below by an
   amount that has not been derived. The independent audit did not test this either.
2. **δ is net aggregate removal.** Captured chains are taken to leave the load, either folded or
   degraded. Disaggregated chains are taken not to return to U. So δ is aggregate degradation k_dA
   plus whatever disaggregation is not followed by re-aggregation. In the full A1 model,
   disaggregated chains re-enter U and add to the load.
3. There is one pool, one client class, no regulation, and no seeded aggregation.

With the capped pool (A1.1a), M^free(A) = max(M^tot − βA, 0) and

    dA/dt = k_a L / (k_c max(M^tot − βA, 0) + k_a + μ) − (μ + δ) A.                   (A6.1)

Notation: a0 = k_a + μ, a = k_c M^tot + a0, b = k_c β, m = μ + δ. A_max = M^tot/β is the
aggregate level at which the whole pool is held.

(A6.1) is piecewise. On [0, A_max] it is k_aL/(a − bA) − mA. Above A_max the chaperone term is
gone and it is k_aL/a0 − mA, which is linear. The earlier A6 used the first formula for every A.
That formula has a pole at a/b = A_max + a0/(k_c β). The pole lies where free chaperone would be
negative, and the earlier "runaway" was the approach to that pole (CORRECTIONS.md, E2).

**No sequestration (β = 0).** There is one steady state, A* = k_a L/(a m). It is linear in L and
globally stable for m > 0 **[checked]**.

**Single bottleneck (β = 0, saturable client binding).** Take one pool and one client class, with
client residence time τ. Then M^free = M^tot/(1 + k_c τ U), and the binding flux
k_c M^tot U/(1 + k_c τ U) increases with U. The steady state is still unique. Throughput caps at
M^tot/τ, and the excess spills to aggregation and degradation, which gives a kink in the
yield-versus-load curve. This is an argument, not a check in `checks.py`. With several pools
sharing clients, uniqueness is not proven. Which machine saturates first depends on τ, which is
unmeasured, so abundance alone cannot name the bottleneck.

**Theorem A6 (capped sequestration, β > 0).** Let H(A) = (m/k_a)·A·(k_c max(M^tot − βA, 0) + a0)
be the load needed to hold aggregate A. Steady states are the solutions of H(A) = L.
Define

    L_crit = m a² / (4 b k_a)  at  A_c = a/(2b),      L_sat = H(A_max) = m M^tot a0 / (β k_a).

Assume first that m = μ + δ > 0.

(a) **Existence and boundedness.** For every L ≥ 0 there is at least one steady state and at most
    three. Every trajectory with A(0) ≥ 0 stays bounded and nonnegative.
(b) **High state.** If L > L_sat, there is a steady state on the titrated piece,

        A_hi = k_a L / ((k_a + μ)(μ + δ))  >  A_max.                                   (A6.2)

    It is finite, linear in L, and stable with eigenvalue −m. It is the aggregate level of a cell
    whose clearance chaperone is entirely held.
(c) **No fold.** If k_c M^tot ≤ k_a + μ, H is nondecreasing. There is exactly one steady state at
    every load, it is globally stable on A ≥ 0, and there is no hysteresis.
(d) **Fold and hysteresis.** If k_c M^tot > k_a + μ, then A_c < A_max and L_sat ≤ L_crit. For
    L_sat < L < L_crit there are three steady states: a stable low state (A < A_c), an unstable
    middle state (A_c < A < A_max) and the stable high state A_hi. Below L_sat only the low state
    exists, and above L_crit only A_hi. The window is

        L_crit − L_sat = m (k_c M^tot − k_a − μ)² / (4 k_c β k_a),
        (L_crit − L_sat)/L_crit = ((k_c M^tot − k_a − μ)/(k_c M^tot + k_a + μ))².        (A6.3)

(e) **No removal.** If m = 0 (no growth dilution and no net aggregate removal) and L > 0, there is
    no steady state. Then dA/dt ≥ k_a L/a > 0 at every A, and A grows without bound.

*Proof.* Write the right side of (A6.1) as f(A) = [k_a/(k_c M^free(A) + a0)]·(L − H(A)). The
bracket is positive, so f has the sign of L − H. H is continuous, H(0) = 0, and for m > 0
H → ∞, so (a) follows. Boundedness holds because f < 0 once H(A) > L. On [0, A_max],
H = (m/k_a)(aA − bA²), which peaks at A_c. Above A_max, H = (m a0/k_a)A is increasing.
A_max − A_c = (k_c M^tot − a0)/(2k_cβ), so the peak lies inside the physical piece iff
k_c M^tot > a0. That gives (c) and (d). In one dimension, a crossing where H increases is stable
and one where H decreases is unstable, so the states alternate stable, unstable, stable. (b) is
the root of the linear piece, and (e) holds because the bracket is at least k_a/a. □

**[checked: a6_symbolic_fold_and_cap. Sympy confirms A_c, L_crit, L_sat, A_max − A_c, the window
(A6.3), the pole position and A_hi.
a6_capped_root_counts_three_state_window. Set P_FOLD has k_c M^tot = 1.75 > k_a + μ = 1.25. It gives
1, 3 and 1 states at 0.9, 0.99 and 1.1 L_crit (L_sat = 2.1875, L_crit = 2.25). At 0.99 L_crit the
states are 1.35 (stable), 1.65 (unstable) and 1.782 (stable, > A_max = 1.75). A brute-force count of
sign changes agrees.
a6_no_fold_without_capture_condition. Set P_NOFOLD has the same a, b, k_a and m but
k_c M^tot = 1.10 < 1.90. It has one state at all 300 loads tested.
a6_high_state_finite_and_attracting. From A = 0 at 1.1 L_crit, A → 1.98 = A_hi, and free
chaperone is never below 0.
a6_hysteresis_loop_between_Lsat_and_Lcrit. In a quasi-static sweep, the jump up comes at 2.263
(L_crit 2.25) and the jump down at 2.192 (L_sat 2.1875), with grid step 0.014.
a6_unbounded_only_without_removal: with μ = δ = 0, A(200) = 898, which is 513 times A_max.
a6_random_topology_local_2000: 0 violations, including "three states ⇔ k_c M^tot > k_a + μ and
L_sat < L < L_crit".
a6_matches_independent_audit. The same numbers come out to 10⁻¹², and the auditor's own
20,000-draw sweep records 0 violations in 1,796 three-state cases.]**

**What the theorem says in plain terms.**
- **At positive dilution or removal there is no runaway.** Above L_crit the aggregate jumps to a
  finite high state proportional to load. Inside the window, a large enough perturbation sends the
  cell to that same finite state, and it stays there until the load drops below L_sat. History
  dependence survives, but as hysteresis between two finite states, not as divergence.
- **Bistability needs two conditions, not one.** They are β > 0 and k_c M^tot > k_a + μ: at full
  pool, chaperone capture of a non-native chain must outpace spontaneous aggregation plus dilution.
  Neither is measured. The second is a property of this **reduced** model only. With ΣC kept, with
  disaggregated chains returning to U, or with several pools, the fold condition is different and
  has not been derived.
- **The loop may be narrow.** By (A6.3) its relative width is ((x − 1)/(x + 1))² with
  x = k_c M^tot/(k_a + μ). At x = 1.4 it is 2.8 % of L_crit. A coarse load ladder can step straight
  over it.
- **When does unbounded accumulation actually occur?** Only when m = 0, that is, with no growth
  dilution and no net aggregate removal. Making disaggregation depend on free chaperone,
  δ = δ_deg + δ_dis M^free/M^tot, does not change this. The high state is then
  k_a L/((k_a+μ)(μ+δ_deg)), finite whenever μ + δ_deg > 0. The auditor tested δ_deg = 0 and 0.05;
  this was not re-run here. A true runaway would also need a separately specified coupling that
  drives μ → 0 as A grows, with δ_deg = 0. No such coupling is modeled.
- **Other feedbacks.** Seeded aggregation, with k_a(A) increasing, and growth-dilution feedback
  μ(A) may produce folds of their own. Their structure with a capped pool has not been derived, and
  each would need its own measured coupling.

None of this defines a death threshold. "Failure" means only that no state meets (A3.1). A_hi is
finite. Whether it violates the quota depends on N^req and on how aggregate affects μ and P, and
neither is specified. L_crit and L_sat depend on β and k_c, which have not been measured in vivo, so
no collapse threshold is stated, universal or otherwise.

## A7. Maintenance, stated formally

Setting. The state is x ∈ Ω (Proposition A1). Growth rate is either a state function μ = g(x),
locally Lipschitz and C¹ where differentiated, or an input. If it is an input, μ is appended to the
state with its own dynamics dμ/dt = h(μ, w), so the system stays autonomous. Disturbances are
w(t) ∈ W, with W compact and measurable in t. They cover bounded perturbations of s, T, k_a, k_u,
k_dN and machine levels. The field f(x, w) is locally Lipschitz in x, uniformly in w, and
continuous in w. N_i^req(μ) is C¹.

Viability set:

    V = { x ∈ Ω : φ_i(x) ≡ N_i − N_i^req(μ(x)) ≥ 0  for i ∈ E }.

A condition with disturbance class W is **maintainable from x₀** if there is a closed set K ⊆ V
with x₀ ∈ K that is forward invariant for every admissible w(·). Under the assumptions above, K is
forward invariant for every admissible disturbance iff f(x, w) lies in the tangent cone T_K(x) for
every x ∈ ∂K and every w ∈ W. This is Nagumo's condition, applied to every w. Necessity follows by
taking w constant. Where K is cut out by finitely many C¹ inequalities with linearly independent
active gradients, T_K(x) is the set of directions v with ∇ψ_j·v ≥ 0 for every active j.

**The quota face.** At every point of ∂K where a quota is active (N_i = N_i^req(μ)):

    dφ_i/dt = P_i(x, w) − (μ + k_u,i + k_dN,i) N_i^req(μ) − (dN_i^req/dμ)·(dμ/dt) ≥ 0,      (A7.1)

with dμ/dt = ∇g(x)·f(x, w) or h(μ, w). The earlier A7.1 left out k_dN and the last term. It is
right only if N^req does not depend on μ, or if μ is stationary on the face (CORRECTIONS.md, E4).
**[checked: a7_face_condition_symbolic.
a7_moving_quota_counterexample_and_fix. Take N^req = 1 + 2μ with μ relaxing from 0.5 to 1. The old
condition holds with zero slack at t = 0, but the corrected left side is −1.0, and the state leaves
V (minimum gap −1.88). Raise P to 3.3, the value that satisfies (A7.1) at every face point
μ ∈ [0.5, 1]. The gap then stays ≥ 0.]**

(A7.1) applies only at boundary points of K where a quota is active. Ω's own faces are handled by
Proposition A1. Any further faces chosen for K, such as bounds on A or on μ, carry their own
tangency conditions, and A7 does not specify those in general. So (A7.1) is **necessary** for
invariance of any K whose boundary meets a quota face. Sufficiency needs the tangent-cone condition
on all of ∂K. Corners, where several faces are active, need all of them at once.

Weaker local version: a steady state x* ∈ int V that is asymptotically stable, with the perturbed
state inside its basin. Inside the A6 window both the low and the high state are asymptotically
stable. Whether the high state lies in V is not determined by A6.

What this does and does not say:
- A disturbance is maintainable only if (A7.1) can hold somewhere. If a disturbance raises k_u,
  k_dN or k_a so far that no state meets it, no invariant subset of V exists and maintenance fails.
  The theory does not assert that every disturbance is maintainable.
- **Smooth decline to slow growth is not guaranteed.** Under (A3.2), if synthesis is strictly
  proportional to growth, growth rates near zero lie outside V for every y. As y falls, the
  feasible growth interval closes at a positive rate μ*. So even with β = 0, feasibility can
  vanish abruptly at nonzero growth. A smooth decline toward μ → 0 needs basal synthesis with
  y s⁰ ≥ (k_u + k_dN) n0 and a closure μ(x) whose slow-growth state is stable. Neither is
  established here.
- F(θ) ≠ ∅ is necessary for any equilibrium in V (Theorem 1). It is not sufficient (Theorem 2,
  and the A6 window).
- Regulation (σ32-driven chaperone induction, titrated by free DnaK) makes M^tot a state variable
  with negative feedback. It is not included in A6 or here. It must be included before any
  quantitative claim.

## A8. The solved competing-fates example (GroEL productive, HtpG disposal, spontaneous folding)

From U: spontaneous folding k_f; GroEL binding k_G G with outcomes N (p_G), back to U (r_G),
to A (a_G = 1−p_G−r_G); HtpG binding k_H H with outcomes N (1−q), degradation (q); aggregation
k_a; disaggregation k_dis; dilution μ. Machines are unsaturated. Y = μN/s is the native fraction of
what is synthesized. Exact solution in `checks.py`.

**Result 1, the HtpG sign rule (a_G = 0)** **[checked symbolically, ratio positive over 20,000 draws]**:

    sign(∂Y/∂H) = sign[ (1−q)(k_a^eff + μ) − q (k_f + p_G k_G G) ],   k_a^eff = k_a μ/(k_dis+μ).

In words, more HtpG raises native yield only when its disposal odds q/(1−q) are smaller than the
ratio of the client's loss rate (unrecovered aggregation plus dilution in the non-native state) to
its productive alternatives (spontaneous folding plus productive GroEL). Both signs occur:
∂Y/∂H = −0.046 with strong alternatives, +0.103 when aggregation dominates **[checked]**.
Efficient disaggregation (large k_dis) lowers k_a^eff and so moves HtpG toward harm.

**Result 2** **[checked]**: ∂A/∂H < 0 always. HtpG always lowers the aggregate burden. Whether
it helps growth therefore depends on whether growth is limited by native function (Y) or by
aggregate toxicity (A). That is a separate, measurable assumption, and the theory does not choose.

**Result 3** **[checked]**: with a_G = 0, ∂Y/∂G ∝ k_G p_G (k_dis+μ)(k_a μ + k_dis μ + μ² + q k_H H(k_dis+μ)) > 0,
so GroEL is monotone beneficial. Counterexample: with an unproductive exit to aggregate (a_G = 0.6,
p_G = 0.1) and a fast spontaneous folder (k_f = 5), ∂Y/∂G = −0.058. So "chaperones are beneficial"
is a parameter-dependent statement even for GroEL. A further cost not modeled here is the
proteome share spent on making GroEL (the FoldME reallocation argument).

**Result 4** **[checked: saturation_preserves_sign]**: with a finite HtpG pool and a single client,
H^free increases with H^tot. So the sign of ∂Y/∂H^tot equals the sign in Result 1. In the harm
regime, yield falls monotonically over a 400-fold range of H^tot. With several clients sharing a
pool, this sign transfer is not proven.
