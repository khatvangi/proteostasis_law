# Formal model: one card per equation

Each card gives the equation, its physical meaning, its assumptions, the quantities a
measurement would need, its limiting cases, what would falsify it, and how it relates to the
three published models:

- **FoldEco**: Powers, Powers & Gierasch 2012, Cell Rep 1:265, doi:10.1016/j.celrep.2012.02.011
- **Santra**: Santra, Farrell & Dill 2017, PNAS 114:E2654, doi:10.1073/pnas.1620646114
- **FoldME**: Chen et al. 2017, PNAS 114:11548, doi:10.1073/pnas.1705524114

What those models contain is taken from `../theory/PRIOR_ART.md`, which read the main texts; their
SI parts were not read. Proofs are in [THEOREMS.md](THEOREMS.md).

Symbols: U, N, A, C_{m,i}, D are non-native, native, aggregated, machine-bound and degraded chains
of client i. M_m^tot and M_m^free are machine pools. s_i is synthesis and σ_i the co-translational
native fraction. μ is growth rate, k_f folding, k_u unfolding of native, k_dN and k_dA degradation of
native and aggregate. κ^X_{m,i} are release rates from machine m to X ∈ {N, U, A, D}, and
τ_{m,i} = 1/(Σκ + μ) is residence time. p^X = κ^X τ. λ_i = μ + k_u,i + k_dN,i.

---

### L0. Conservation of chains and machines

    d/dt(U+N+A+ΣC+D)_i = s_i − μ(U+N+A+ΣC)_i ;   M_m^free = max(M_m^tot − Σ_i C_{m,i} − β_m Σ_i A_i, 0)

- **Meaning.** Every chain made is somewhere or has been degraded. Machine units are free, holding
  a client, or held on aggregate.
- **Assumptions.** Well-mixed cell. Sequestration is fast and capped (β is a capacity).
  Co-chaperones obey their own version.
- **Testable.** Nothing directly; it is bookkeeping. Pulse-chase mass closure
  (total = soluble + insoluble + degraded) tests the accounting.
- **Limits.** β = 0: no sequestration. μ = 0 and no degradation: chains accumulate without bound.
- **Falsifier.** A pulse-chase in which labelled chains are not accounted for by the compartments,
  which would mean a missing state such as membrane insertion or secretion.
- **Prior art.** The same accounting underlies FoldEco and Santra. FoldEco has no dilution.
  FoldME has no aggregate or degraded state. The cap on M^free is a correction this project made
  (`../theory/CORRECTIONS.md` E1).

### L1. Native replacement quota (T1)

    dN_i/dt = P_i − λ_i N_i ;   steady: P_i ≥ λ_i N_i^req(μ) ;   face: P_i − λ_i N_i^req − (dN_i^req/dμ) μ̇ ≥ 0

- **Meaning.** Native protein is lost by dilution, unfolding and regulated proteolysis. Productive
  flux must replace it, and must also chase a moving target when growth changes.
- **Assumptions.** N^req(μ) exists, is C¹, and is client-specific. That is a hypothesis about how
  growth depends on the protein, not a law.
- **Testable.** N^req(μ) by promoter or CRISPRi titration at several growth rates. k_dN by
  pulse-chase of the folded protein. μ.
- **Limits.** μ → 0: the quota tends to (k_u + k_dN)N^req(0) > 0. Strictly growth-proportional
  synthesis then cannot meet it (`../theory/EQUATIONS.md` A3.2). If N^req is constant, the face term
  vanishes.
- **Falsifier.** A protein held below its measured N^req at steady state while growth continues
  at that μ, which would mean N^req is not a requirement. Or growth that persists with P measured
  below λN^req, which would mean a missing productive route.
- **Prior art.** Santra reaches a non-vanishing slow-growth chaperone need by simulation; here it
  follows from the balance. Neither FoldEco nor Santra has native degradation in the quota or the
  moving face term. FoldME imposes demand through growth optimization.

### L2. Machine occupancy and productive-flux ceiling (T2)

    Σ_i J_{m,i} τ_{m,i} ≤ M_m^tot − B_m ;   σ_i s_i ≤ P_i ≤ σ_i s_i + k_f,i U_i + Σ_m (p^N_{m,i}/τ_{m,i}) M_m^avail,i

- **Meaning.** A machine unit holding a client cannot hold another. Productive output through a
  machine is at most (success probability / residence time) × units available.
- **Assumptions.** One client per unit (GroEL ring, DnaK monomer: a stoichiometry convention).
  Visits are counted separately.
- **Testable.** M^tot (Schmidt 2016, measured). τ and p^N in vivo for a machine–client pair:
  **not measured** for any E. coli pair (`EXPERIMENT_CONSTRAINTS` NONE-1). σ s, by ribosome
  profiling for s; σ not measured.
- **Limits.** τ → 0: the machine never limits. σ s ≥ λN^req: the client never needs a machine.
- **Falsifier.** A measured occupancy share above 1, which would mean a wrong stoichiometry
  convention, or productive flux above the ceiling, which would mean an unmodelled route.
- **Prior art.** Santra reports GroEL near 100% occupancy at 40-min doubling, a model output.
  FoldME's coupling constraints (folding flux ≤ k_cat × chaperone) are the same kind of
  inequality, used as optimization constraints. Here it is a necessary condition with its
  measurable parts named.

### L3. Time-to-depletion window (T3)

    N_0 e^{−λ_hi t} + (P_lo/λ_hi)(1−e^{−λ_hi t}) ≤ N(t) ≤ N_0 e^{−λ_lo t} + (P_hi/λ_lo)(1−e^{−λ_lo t})
    t_up = λ_lo⁻¹ ln[(N_0 − P_hi/λ_lo)/(N^req − P_hi/λ_lo)]   (if P_hi < λ_lo N^req)
    t_lo = λ_hi⁻¹ ln[(N_0 − P_lo/λ_hi)/(N^req − P_lo/λ_hi)]   (if P_lo < λ_hi N^req)

- **Meaning.** Bounds on productive flux and loss rate give a time window in which a
  depletion-driven quota failure must occur, or cannot yet have occurred.
- **Assumptions.** Bounds hold over the window. N^req is bounded within [R_lo, R_hi].
- **Testable.** N_0, N^req, μ, k_dN (measurable). P_hi needs T2 (τ, p^N unmeasured). Without
  P_lo, the floor t_lo ≥ λ_hi⁻¹ ln(N_0/N^req) still holds.
- **Limits.** P_hi → 0: pure decay, t_up = λ_lo⁻¹ ln(N_0/N^req). P_hi → λ_lo N^req: t_up → ∞.
- **Falsifier.** A native-activity time course that crosses N^req outside [t_lo, t_up] with the
  bounds independently measured. That would mean an unmodelled productive or loss route. A
  conditional GroE-depletion strain (Fujiwara 2010 used one) gives exactly this experiment for
  obligate clients.
- **Prior art.** None of the three states time-to-failure bounds. FoldEco and Santra simulate
  trajectories at given parameters. FoldME is steady state.

### L4. Route-qualified feasibility with a Farkas certificate (T4)

    Σ_r a_{m,r} x_r ≤ M_m − B_m ;  Σ_{r∈i}(1−b_r)x_r ≤ e_i ;  Σ_{r∈i} ρ_r x_r ≥ d_i ;  x_r ≤ c_r ;  x ≥ 0
    infeasible ⇔ ∃ (z, u, w, c') ≥ 0:  Σ_m w_m a_{m,r} + u_i(1−b_r) + c'_r ≥ z_i ρ_r  ∀r,
                                         Σ z_i d_i > Σ w_m(M_m−B_m) + Σ u_i e_i + Σ c'_r c_r

- **Meaning.** Is there any way, routing chains through qualified routes with their own
  machine-time costs and yields, to cover every essential deficit? If not, there is a pricing of
  machine-seconds that proves it.
- **Assumptions.** Linear route data a, ρ, b at the condition. A machine-free loss route exists.
  For hand-off routes, lower bounds on a. The split x is free, which is what makes the condition
  only necessary.
- **Units.** a [time], x, e, d, c [conc/time], M [conc], w [1/time], z and u dimensionless.
- **Testable.** Route qualification (which machines can serve which client, for example the
  Kerner obligate class or the Vorderwülbecke direct GroEL route). τ and visit counts (unmeasured
  in vivo). M (measured).
- **Limits.** Single-machine routes with ρ = 1 and machine-specific τ: reduces to Gale's
  supply–demand cut condition (tested). With serial routes no cut condition is exact. Global pooling
  is a strict relaxation.
- **Falsifier.** A strain that maintains growth while a Farkas certificate built from measured
  M, τ and route qualification proves infeasibility. That would mean an unrecognized route or
  wrong τ.
- **Prior art.** FoldME's chaperone-coupling constraints are the closest relative (LP/ME family).
  What differs: route qualification with serial multi-machine coefficients, yields and returns;
  the certificate and its units; and use as a necessary condition rather than as an optimization
  the cell is assumed to solve. FoldEco and Santra do not formulate feasibility. The earlier
  project route model used Hall-style reasoning; T4(c) shows why that fails for serial routes.

### L5. Realized competing-fates kinetics (T5)

    dU/dt = (1−σ)s + k_u N + Σ_m κ^U C_m + v_dis − (k_f + k_a + Σ_m k_on,m M_m^free + k_p + μ) U
    dC_m/dt = k_on,m M_m^free U − (κ^N + κ^U + κ^A + κ^D + μ)_m C_m
    per-entry yield at fixed free pools: Y = (k_f + Σ θ_m p^N_m)/(k_f + k_l + Σ θ_m(1−p^U_m)),  θ_m = k_on,m M_m^free

- **Meaning.** The cell does not choose routes. At the branch point, hazards compete, and free
  pools set the binding hazards.
- **Assumptions.** Mass-action hazards. Each bound state's outcome probabilities are fixed by
  machine–client constants. Linear aggregation in the closed form (`../kinetics/` treats second
  order).
- **Testable.** Committed success ρ̂_m = p^N/(1 − p^U) per machine–client pair: refolding yield
  per binding event. In vitro values exist for model substrates; in vivo values do not.
- **Limits.** θ → 0: Y → k_f/(k_f + k_l). θ_m → ∞: Y → ρ̂_m.
- **Falsifier.** A route whose overexpression raises yield in a background where its measured
  ρ̂ is below the measured Y (contradicts S1a at fixed other pools).
- **Prior art.** Kinetic competition is FoldEco's and Santra's core, and is **not new**. New here:
  the exact sign rule with the committed-success threshold, and its formal separation from
  capacity (S1). Neither model includes HtpG or HslUV.

### L6. σ32 free-machine feedback (T8; dynamics in ../feedback)

    dM_m/dt = b_0,m + b_1,m h(F(t−τ)) − γ M_m ,  F = free chaperone signal,  h decreasing, 0 ≤ h ≤ 1

- **Meaning.** Free DnaK/J and GroEL/S repress σ32 (activity and FtsH-mediated stability). Load
  titrates them and raises regulon output after a delay.
- **Assumptions.** One lumped regulator. Synthesis bounded by b_0 + b_1. Monotone regulon.
- **Testable.** Regulon induction rate and ceiling (pulse labelling after a shift). Loop gain G
  and delay τ: not measured (`../feedback/EXPERIMENTAL_CONTROL.md`).
- **Limits.** b_1 → 0: no feedback. τ → 0: stable at any gain. G ≤ 1 ⇔ stable at every delay
  (feedback_review).
- **Falsifier.** Chaperone accumulation faster than the bound in T6(a) allows, given measured
  synthesis. Or strong induction while free chaperone is high (T8c).
- **Prior art.** The loop itself is established (Straus 1987, Tomoyasu 1998, Guisbert 2004,
  El-Samad 2005), and none of FoldEco, Santra or FoldME models it. New here: the capacity set is
  up-closed, so co-regulation is not a feasibility constraint (T8a), but realized yield can be
  non-monotone in regulon output (T8b), and free-machine sensing is blind to routing failure (T8c).

### L7. Controlled robust invariant viable set (T6, T7)

    K ⊆ V closed, f(x,u,w) ∈ T_K(x) ∀x ∈ ∂K, ∀w ∈ W (some admissible u) ;
    race bound: N(t) ≤ N̄(t) = N_0 e^{−λt} + (p_0+κM_c)/λ (1−e^{−λt}) + κ(M_0−M_c)(e^{−γt}−e^{−λt})/(λ−γ)

- **Meaning.** Maintenance is a property of trajectories: the cell must stay in V through
  disturbances, under the regulation it actually has. The race bound caps native recovery after a
  shock by the fastest machine induction any regulator could produce.
- **Assumptions.** H1–H4 of T6. Bounded disturbance class.
- **Testable.** λ after the shift, M_0, b_0 + b_1, γ, p^N/τ, σ s, N_0, N^req. Then
  min N̄ < N^req predicts a transient loss of function.
- **Limits.** λ ≪ γ: induction keeps up and the dip vanishes. M_0 ≥ M_c: no induction needed.
- **Falsifier.** A measured native time course that rises above N̄(t) computed from measured
  inputs.
- **Prior art.** Santra and FoldEco simulate shifts at fitted parameters. None states
  viability-kernel conditions or a regulator-independent transient bound.

### L8. Codon-conditioned perturbations (T9)

    Y_tot = σ + (1−σ) Σ_k π⁰_k Y_k ;  ∂Θ*/∂a_{m,r} = −w_m* x_r* ;  variant fraction ≤ Σ_c n_c ε_c

- **Meaning.** Synonymous codons act on demand and routing: co-translational folding σ, the
  birth distribution π⁰ over held states, and route coefficients. Misincorporation adds variant
  species. Effects on growth are conditional on whether the loaded machine binds.
- **Assumptions.** Codon effects enter through σ, π⁰ and a, not through new routes. ε is per-codon
  and independent across codons.
- **Testable.** σ and π⁰ per variant (co-translational folding assays, ribosome-bound chaperone
  profiling). ε (Landerer; detection-conditional). Machine capacity regime, by titration.
- **Limits.** Slack machine: no first-order effect. σ → 1: no machine demand.
- **Falsifier.** A synonymous variant whose growth effect does not depend on the titrated level
  of the machine it loads (contradicts T9b), or a codon-driven yield change whose sign disagrees
  with the measured change in σ or π⁰ (contradicts T9a).
- **Prior art.** FoldEco mimics translation rate by about a 2× change in k_f. None of the three
  treats synonymous codons or misincorporation.

### L9. Optional energy row (T10)

    Σ_r E_r x_r ≤ E_max ,  E_r = ATP per route + (1−ρ_r) c_syn

- **Meaning.** Routes cost ATP, and lossy routes cost resynthesis.
- **Assumptions.** Linear costs. A budget exists. This is an extra constraint, not an axiom.
- **Testable.** ATP per chaperone cycle (in vitro). Resynthesis cost per chain (bookkeeping).
  The budget: not defined by any measurement.
- **Limits.** E_max large: the row is inactive and the margin is the only objective.
- **Falsifier.** None until a budget is measured; it is optional.
- **Prior art.** FoldME's proteome reallocation is a cost of the same kind. Its biomass objective
  makes the trade-off central. Here it is optional and explicit.
