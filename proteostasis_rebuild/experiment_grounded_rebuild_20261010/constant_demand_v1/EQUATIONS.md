# Equations, units and proofs (constant demand, steady state)

Numbering matches `model.py`. Nothing here varies in time except N(t) relaxing to steady state;
demand, environment and parameters are constant.

## Units

| symbol | meaning | units |
|---|---|---|
| D | required function, fixed by the environment | fu/min (fu = function unit, e.g. substrate molecules converted) |
| N | **functional** copies of the client protein (active, folded, cofactor-loaded) | molecules/cell |
| a | function delivered per functional molecule under the fixed condition | fu/(molecule·min) |
| Q = aN | available functional capacity | fu/min |
| P | supply of *new functional* molecules | molecules/min |
| λ | loss of function per molecule: dilution + degradation + inactivation | 1/min |
| s | chains synthesized | chains/min |
| J_r | chains entering processing route r | chains/min |
| y_r | fraction of chains on route r that end functional | dimensionless |
| c_{m,r} | machine m occupied per chain on route r (dwell time × units held) | machine·min/chain |
| M_m | machine m complexes present | complexes/cell |
| λ_m | loss rate of machine m | 1/min |
| ℓ | synthesis cost of one chain or one machine complex | aa |
| R | synthesis budget shared by clients and machines | aa/min |
| ε | molecules inactivated per fu delivered (mechanism-based loss) | molecule/fu |

## Scalar model

    Q = a N                                                         (E1)
    dN/dt = P − λ N                                                 (E2)
    N* = P/λ,   N(t) = N* + (N0 − N*) e^{−λt}                        (E3)
    Q* = a P / λ                                                    (E4)
    N_req = D/a,   P_req = λ D / a                                  (E5)

**Theorem (E6).** Assume (i) P is a constant chosen anywhere in [0, Pmax], independent of N;
(ii) λ > 0 and a > 0 are constants; (iii) nothing about the route that produces P depends on N or
on Q. Then a steady state with Q* ≥ D exists **if and only if**

    D ≤ a Pmax / λ.                                                 (E6)

*Proof.* By (E3) the steady state is unique, N* = P/λ, and globally attracting because λ > 0.
Q* = aP/λ is increasing in P. Necessity: if D > aPmax/λ then every admissible P gives
Q* ≤ aPmax/λ < D. Sufficiency: P = P_req = λD/a lies in [0, Pmax] and gives Q* = D. ∎

Checked against a brute-force grid over [0, Pmax] in 500 random draws (`test_condition_is_necessary_and_sufficient`).

**Installed, sustainable, potential (E7).**

    Q_inst(t) = a N(t)      Q_sus = a P_real / λ      Q_pot = a Pmax / λ,     Q_sus ≤ Q_pot

**Time until installed capacity stops covering demand (E8).** If N0 > N_req > N* = P_real/λ:

    t_short = ln[(N0 − N*)/(N_req − N*)] / λ

**Service-dependent loss (E9).** Only if each fu delivered destroys ε molecules does demand cost
protein directly. With Q ≥ D, loss = λN + εD, so

    P_req = D (λ/a + ε),    condition  D ≤ Pmax / (λ/a + ε).

ε = 0 for ordinary catalysis: the enzyme is not used up by its reaction.

**Route feedback breaks sufficiency (E10).** If P = P(N) ≤ Pmax, then N* = P(N*)/λ ≤ Pmax/λ, so
(E6) stays **necessary**. It is no longer sufficient. Example P(N) = Pmax − βN gives
N* = Pmax/(λ + β), and any D in (aPmax/(λ+β), aPmax/λ] passes (E6) and is not met.

**Machine self-maintenance (E11).** A machine pool is held, not stored:

    dM_m/dt = P_m − λ_m M_m   ⇒   holding M_m costs ℓ_m λ_m M_m aa/min, indefinitely.

## Where Pmax comes from: finite, machine-specific processing

    P = Σ_r y_r J_r,     Σ_r J_r ≤ s                                (E12)
    Σ_r c_{m,r} J_r ≤ M_m        for every machine m                 (E13)

(E13) is Little's law: mean machine units occupied = arrival rate × mean time held. A serial route
(K then G) has c_{K,r} > 0 and c_{G,r} > 0: one chain consumes both. Machines are not
interchangeable: each m has its own row.

**Pmax** = max (E12) subject to (E13), J ≥ 0, with routing treated as freely allocable. This is
an upper bound, not a prediction: nothing in the cell solves this LP.
**P_real** = Σ y_r J_r at the J_r actually produced by kinetic competition. Because the realized
J is one feasible point, P_real ≤ Pmax (100 random draws, `test_realized_never_exceeds_potential`).

Kinetic routing used in the example: a new chain commits to the first of four competing events,
hazards k_f (fold), k_a (aggregate), b_K K_free, b_G G_free, with free pools solving (E13) as an
equality for the occupied part, M_free = M_tot − Σ_r c_{m,r} J_r.

## Many clients, many functions

    A N ≥ D,   N_i = P_i/λ_i  (P = Λ N, Λ = diag λ_i),   P_i = Σ_r y_{ir} J_{ir}     (E14)

A_{ji} [fu_j/(molecule_i·min)] is the function j delivered per functional molecule of protein i.
Write the decision vector x = (J_{ir} ≥ 0, M_m ≥ 0). Then F = A Λ^{-1} Y and the problem is

    find x ≥ 0 with  F x ≥ D  and  G x ≤ h                          (E15)

where the rows of G x ≤ h are (E13) for each machine, per-client synthesis caps if given, and

    dM_m/dt = P_m − λ_m M_m = 0                                     (E16)
    Σ_{i,r} ℓ_i J_{ir} + Σ_m ℓ_m λ_m M_m ≤ R                        (E17)

(E17) charges failed chains too, and charges machines for their own turnover.

**Farkas alternative (E18).** Exactly one of the following holds:

1. some x ≥ 0 satisfies F x ≥ D and G x ≤ h;
2. some prices u ≥ 0 (per unit of each resource row) and w ≥ 0 (per unit of each demanded
   function) satisfy  Gᵀu ≥ Fᵀw  and  hᵀu < Dᵀw.

*Proof.* (E15) is {x ≥ 0 : [G; −F] x ≤ [h; −D]}. Farkas' lemma (inequality form) says this is
empty iff some y = (u, w) ≥ 0 has [G; −F]ᵀ y ≥ 0 and [h; −D]ᵀ y < 0, which is statement 2. ∎

Reading: every route and every machine is "worth" (via w) no more than the resources it uses
(via u), yet the demand is worth more than all resources together. Units: w_j in 1/(fu_j/min),
u in 1/(resource unit); both sides dimensionless; only the ratios of prices matter.
A certificate is checked by hand (`checkCertificate`), never trusted from the solver.

**Scalar reduction.** One client, one function: (E15) is feasible iff D ≤ (a/λ) · Pmax, with
Pmax the LP value of (E12–E13). That is (E6) again (`test_lp_reduces_to_scalar_condition`).

**Why not Hall's condition.** Hall-type subset counting ("for every set S of clients, the
demand of S fits in the machines S touches") is exact when each unit of flux uses exactly one
unit of exactly one machine. Serial routes break that. Counterexample (`test_hall_counting_...`):
K = G = 1 complex, dwell 1 min; client 1 must go K→G, client 2 needs G; required fluxes 1 and
0.5. Every subset passes the count ({1}: 1 ≤ 2, {2}: 0.5 ≤ 1, {1,2}: 1.5 ≤ 2), but G must carry
1.5 > 1. Infeasible, with a verified certificate. We therefore claim nothing Hall-like.

## Assumed vs derived

| ASSUMED | DERIVED from the assumptions |
|---|---|
| D constant and externally fixed | N_req = D/a, P_req = λD/a (E5) |
| Q = aN, a constant (no substrate or crowding change with N) | unique, attracting steady state (E3) |
| first-order loss λN with λ > 0 | D ≤ aPmax/λ is necessary and sufficient under (i)–(iii) (E6) |
| P constant, controllable on [0, Pmax], independent of N | with feedback P(N): necessary only (E10) |
| no service-dependent inactivation (ε = 0) unless stated | turnover, not service, sets maintenance cost (E5, E9) |
| route yields y and machine costs c constant | occupancy bound (E13) is exact for averages (Little's law) |
| machines fold without assistance, lose at λ_m | machine maintenance consumes budget, can decide feasibility (E17) |
| linear resource use | Farkas certificate is exact for (E15) (E18) |
| — | Hall counting is not sufficient with serial routes |
| — | P_real ≤ Pmax; equal machine counts do not fix P_real |

Not assumed: that the cell optimizes, that chaperones help, or that any machine is saturated.
