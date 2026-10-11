# Equations and proofs (appendix to SCIENTIFIC_THEORY.md)

Every sign, threshold and counterexample marked **[checked]** is recomputed by `checks.py`
(10/10 pass, output in `checks.json`). Parameter values in `checks.py` are illustrative and
dimensionless. None is an E. coli estimate.

---

## A1. State variables and conservation

For each client species i (a protein, or a class of proteins treated as one):

| symbol | meaning |
|---|---|
| U_i | soluble non-native chains: released nascent chains not yet folded, misfolded, unfolded natives |
| N_i | native, functional chains |
| A_i | chains in aggregates |
| C_{m,i} | chains bound to machine m (one variable per machine, per client) |
| D_i | cumulative chains removed by degradation (a counter, not a pool) |

For each machine m (TF, DnaK/DnaJ/GrpE, GroEL/GroES, HtpG, ClpB, Lon, ClpXP, HslUV, …):

    M_m^tot = M_m^free + Σ_i C_{m,i} + β_m Σ_i A_i                       (A1.1)

β_m is the number of machine units held per aggregated chain. β_m = 0 means no sequestration.
It is the one parameter that switches on the positive feedback in A6. Cochaperones obey their own
version of (A1.1). An effective rate that needs a cochaperone (GroES for GroEL, DnaJ/GrpE for
DnaK) is a function of the cochaperone's free pool too.

Process fluxes (per volume per time, all ≥ 0):

| process | from → to | flux |
|---|---|---|
| birth (synthesis) | – → chains | s_i(t) = elongation-limited synthesis rate |
| co-translational native production | birth → N | σ_i s_i |
| post-translational entry | birth → U | (1−σ_i) s_i |
| spontaneous folding | U → N | k_f,i U_i |
| misfolding / unfolding of native | N → U | k_u,i N_i |
| machine binding | U → C_m | r^on_{m,i} M_m^free U_i |
| productive release | C_m → N | κ^N_{m,i} C_{m,i} |
| unproductive release | C_m → U | κ^U_{m,i} C_{m,i} |
| release into aggregate | C_m → A | κ^A_{m,i} C_{m,i} |
| hand-off to degradation | C_m → D | κ^D_{m,i} C_{m,i} |
| aggregation | U → A | k_a,i U_i (or a nucleated form; see A6) |
| disaggregation | A → U | v^dis_i(A, M_ClpB^free, M_DnaK^free) |
| protease degradation of U | U → D | Σ_p r^on_{p,i} M_p^free U_i (proteases are machines p with κ^D = 1) |
| dilution by growth | every pool | μ × pool |

Balance equations:

    dU_i/dt = (1−σ_i)s_i + k_u,i N_i + Σ_m κ^U_{m,i} C_{m,i} + v^dis_i
              − [k_f,i + k_a,i + Σ_m r^on_{m,i} M_m^free + μ] U_i                (A1.2)
    dN_i/dt = σ_i s_i + k_f,i U_i + Σ_m κ^N_{m,i} C_{m,i} − (k_u,i + μ) N_i      (A1.3)
    dA_i/dt = k_a,i U_i + Σ_m κ^A_{m,i} C_{m,i} − v^dis_i − μ A_i                (A1.4)
    dC_{m,i}/dt = r^on_{m,i} M_m^free U_i − (κ^N + κ^U + κ^A + κ^D + μ)_{m,i} C_{m,i}   (A1.5)
    dD_i/dt = Σ_m κ^D_{m,i} C_{m,i}                                             (A1.6)

Summing (A1.2–A1.6) gives the conservation law

    d/dt (U_i + N_i + A_i + Σ_m C_{m,i} + D_i) = s_i − μ (U_i + N_i + A_i + Σ_m C_{m,i}).   (A1.7)

Every synthesized chain ends up either diluted in a pool or counted in D_i. No flux is
unaccounted for.

## A2. State-dependent effective rates

Every rate is a function, not a constant:

    r_{k,r} = r_{k,r}(T, c, z_k, C_free)

- T is temperature. Folding, unfolding, aggregation and binding rates all have different
  temperature dependences. No monotone "temperature = load multiplier" is assumed; the
  Kanemori plating table already rules that out (EXPERIMENT_CONSTRAINTS row K99).
- c is the cochaperone and nucleotide state (GroES, DnaJ/GrpE, ATP/ADP).
- z_k is client-intrinsic. It includes the sequence and **the elongation-time profile set by
  synonymous codon choice**.
- C_free is the vector of free machine pools from (A1.1).

The synonymous-codon variable z_k enters only through σ_k and the partition among
k_f, k_a and r^on_{m,k} at the moment of release. That is a timing effect. It is a different
parameter from the substitution frequency ε_k that mass spectrometry measures. ε_k creates a
separate variant species k′ with its own rates. The Landerer and Stikeleather tables constrain ε
(standing frequency, detection-conditional). They contain no information about σ or the partition.

Synthesis rate s_k enters (A1.2) linearly as load. In exponential growth s_k scales with μ, so
faster growth raises load and dilution together.

## A3. The productive essential quota (necessary condition 1)

Let E be the set of essential species. Each needs a standing native level N_i ≥ N_i^req(μ) to
support growth rate μ. At steady state, (A1.3) gives

    P_i ≡ σ_i s_i + k_f,i U_i + Σ_m κ^N_{m,i} C_{m,i} = (μ + k_u,i) N_i.

So a steady state with N_i ≥ N_i^req requires

    P_i ≥ (μ + k_u,i) N_i^req      for all i ∈ E.                                (A3.1)

Two consequences follow directly.
- As μ → 0 the quota does not vanish. It tends to k_u,i N_i^req, the maintenance flux needed to
  repair spontaneous unfolding. This matches the Santra 2017 model statement that very slow growth
  still needs chaperones. Here it follows from (A1.3) without fitting anything.
- The quota bounds productive flux. It does not bound total synthesis. Synthesis minus productive
  flux must leave the cell through aggregation, degradation or dilution of non-native pools (A1.7).

## A4. Machine occupancy (necessary condition 2): Little's law

At steady state, (A1.5) gives C_{m,i} = J_{m,i} τ_{m,i}, where J_{m,i} = r^on M^free U_i is the binding
flux and τ_{m,i} = 1/(κ^N+κ^U+κ^A+κ^D+μ)_{m,i} is the mean residence (service) time. Substituting
into (A1.1):

    Σ_i J_{m,i} τ_{m,i} = M_m^tot − M_m^free − β_m Σ A_i ≤ M_m^tot.             (A4.1)

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

A third gap between feasibility and maintenance comes from bistability (A6): F can be non-empty
and a good steady state can exist while the cell sits in the basin of a bad one.

## A6. Overload: graceful failure, single bottleneck, or runaway

Take one aggregation-prone load L (synthesis of non-native chains). DnaK-like clearance has
free pool M^free = M^tot − βA. Holding chaperone binding in quasi-steady state,
U = L/(k_c M^free + k_a + μ), and

    dA/dt = k_a L / (a − bA) − (μ + δ) A,   a = k_c M^tot + k_a + μ,  b = k_c β.     (A6.1)

δ lumps disaggregation and aggregate degradation.

**No feedback (β = 0).** There is a unique steady state, A* = k_a L / (a(μ+δ)), linear in L
**[checked]**. Damage rises smoothly with load and there is no threshold. If growth is coupled to
native function (μ = μ(N)), falling μ lowers synthesis (s ∝ μ) and dilution together. This gives a
stable slower-growth state. Call it **graceful failure**: performance falls continuously until a
quota in (A3.1) cannot be met even as μ → 0.

**Single bottleneck (β = 0, saturable machine).** With finite pools, the machine with the smallest
M^tot/(Σ J τ) saturates first (M^free → 0). Its throughput caps at M^tot/τ, and excess client
flux spills to the competing exits (aggregation, degradation). The yield-versus-load curve then has
a kink. It is still continuous and still has a unique steady state. Which machine binds is set by
τ, which is unmeasured. Abundance alone cannot name the bottleneck.

**Runaway (β > 0, specified positive feedback).** The right side of (A6.1) is convex in A on
[0, a/b). Steady states are the intersections of a convex increasing curve with the line
(μ+δ)A. There are two of them (lower stable, upper unstable threshold) or none. A saddle-node
occurs at

    A_c = a / (2b),   L_crit = (μ+δ) a² / (4 b k_a)                                 (A6.2)

**[checked: sequestration_saddle_node; 2 steady states at 0.9 L_crit, 0 at 1.1 L_crit]**.
Above L_crit the aggregate grows until chaperone is fully titrated. Below it, a large enough
perturbation crosses the unstable branch and produces the same runaway. That is a history
dependence a no-feedback model cannot produce.

Two other specified feedbacks give the same structure:
- seeded aggregation, k_a = k_a(A) increasing;
- growth-dilution feedback, where aggregates lower μ and lower μ removes less aggregate.
Each needs its own measured coupling.

None of this defines a death threshold. "Failure" means only that no steady state or invariant set
meets (A3.1). Growth arrest, filamentation or death need further state equations.

L_crit depends on β, which has not been measured in vivo. A universal collapse threshold
cannot be stated.

## A7. Maintenance, stated formally

Viability set V = { x ≥ 0 : N_i ≥ N_i^req(μ(x)) for i ∈ E }.

A condition θ, with disturbance class W (bounded perturbations of s, T, k_a, machine levels),
is **maintainable from x₀** if there is a closed set K ⊆ V with x₀ ∈ K that is forward-invariant
for every admissible disturbance trajectory. By Nagumo's theorem (for locally Lipschitz
dynamics), a closed K is forward-invariant iff on its boundary the vector field lies in the
tangent cone. For the quota face N_i = N_i^req, that means

    P_i(x, w) ≥ (μ + k_u,i) N_i^req    for every boundary state x ∈ ∂K and every w ∈ W.   (A7.1)

Weaker, local version: a steady state x* ∈ int V that is asymptotically stable, with the
perturbed state inside its basin.

What this does and does not say:
- A disturbance is maintainable only if (A7.1) can hold somewhere. If a disturbance raises k_u or
  k_a so far that no state meets (A7.1), no invariant subset of V exists and maintenance fails. The
  theory does not assert that every disturbance is maintainable.
- F(θ) ≠ ∅ is necessary for any equilibrium in V (Theorem 1). It is not sufficient (Theorem 2, A6).
- Regulation (σ32-driven chaperone induction, titrated by free DnaK) makes M^tot a state variable
  with negative feedback. It enlarges K and is not included in A6. Including it is required before
  any quantitative claim.

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
