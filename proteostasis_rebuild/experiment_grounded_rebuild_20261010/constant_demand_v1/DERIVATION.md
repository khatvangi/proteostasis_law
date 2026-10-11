# 1. Holding a fixed function against turnover

*Draft first section for a theoretical proteostasis paper. Equation numbers refer to
EQUATIONS.md. Every number below is HYPOTHETICAL unless stated otherwise and comes from
`python model.py` (saved in `example_output.txt`).*

## 1.1 Three levels under a constant demand

We start from the simplest case: the environment is constant and requires a fixed amount of one
function per unit time, D (fu/min, e.g. substrate molecules converted per minute per cell).
Three quantities have to be kept apart.

1. **Required function, D.** Set by the environment, not by the cell.
2. **Available functional capacity, Q = aN.** N counts *functional* copies of the protein that
   performs the function, and a is what one functional copy delivers per minute under this
   condition. Q has the same units as D.
3. **Maintenance capacity.** Functional copies are lost at rate λN (dilution by growth,
   degradation, inactivation), so they must be replaced: dN/dt = P − λN, with P the supply of
   new *functional* molecules per minute.

## 1.2 The scalar condition

At steady state N* = P/λ and Q* = aP/λ. The demand is met when Q* ≥ D, i.e. when the supply is at
least P_req = λD/a. If P can be set anywhere in [0, Pmax], independently of N, and λ > 0, then a
steady state meeting the demand exists if and only if

    D ≤ a Pmax / λ.                                                 (E6)

The proof is two lines (EQUATIONS.md): Q* is increasing in P, and P_req is admissible exactly
when (E6) holds. The theorem says an admissible supply *exists*. It does not say the cell
chooses it.

The condition holds only under its three assumptions. If the supply depends on N, for example
because the product competes for its own folding machinery, (E6) stays necessary but is no longer
sufficient. With P(N) = Pmax − βN the steady state is aPmax/(λ + β), so demands just below
aPmax/λ pass (E6) and still go unmet (E10).

## 1.3 What costs protein: turnover, not service

An enzyme is not consumed by its reaction. Delivering D costs the cell nothing per reaction; the
cost is replacing copies that are lost, λN_req = λD/a molecules per minute. Doubling λ doubles the
cost. If λ were zero the installed copies would serve forever without new synthesis. Demand enters
the cost directly only when the reaction destroys the catalyst (ε molecules per fu), in which case
P_req = D(λ/a + ε) (E9).

Two distinctions follow.

- **Stationary is not static.** A constant N* is a flux of λN* molecules per minute passing
  through the pool. The same holds for every chaperone and protease pool: holding M_m costs
  λ_m M_m complexes per minute, indefinitely (E11).
- **Installed is not sustainable.** Installed capacity Q_inst(t) = aN(t) is what is present now.
  Sustainable capacity Q_sus = aP_real/λ is what the supply can hold. A cell carrying a surplus
  from an earlier condition can meet the demand for a while and then fail. In the hypothetical
  example below, N0 = 5000 copies against N_req = 3000 and a sustainable N* = 2184 crosses below
  N_req after 61.9 min (E8). Measuring abundance once cannot tell these apart.

## 1.4 Where Pmax comes from

The supply of functional protein passes through finite, machine-specific processing. A new chain
takes one of several routes r (spontaneous folding, DnaK-like holding, GroEL-like folding, DnaK
then GroEL in series). Route r turns a fraction y_r of its chains into functional protein and
occupies machine m for c_{m,r} machine·min per chain. Machines are not interchangeable, so each has
its own occupancy limit, Σ_r c_{m,r} J_r ≤ M_m (E13; Little's law). A serial route draws on both
machines it passes through. The machines themselves must be resynthesized at λ_m M_m, which uses
the same synthesis budget as the clients (E17).

Pmax is the largest productive flux these constraints allow *if routing could be allocated
freely*, a linear program (E12–E13). It is a ceiling, not a prediction. What the cell actually
gets, P_real, is set by kinetic competition for each new chain, and P_real ≤ Pmax always.

## 1.5 Same chaperone counts, different outcome (hypothetical)

One client, synthesized at 100 chains/min, with D = 1.8×10⁶ fu/min, a = 600 fu/(molecule·min) and
λ = 0.02/min. Then N_req = 3000 copies and P_req = 60 molecules/min. The cell has 200 K and 40 G
complexes. Left alone, a chain folds with yield 0.3. Via G the yield is 0.9, and via K (then G
for chains K cannot fold) it is 0.95.

| | Pmax (LP) | P_real | K occupied | G occupied | meets D? |
|---|---|---|---|---|---|
| good routing (b_K = 0.05, b_G = 0.02 per molecule·min) | 95.0 | 89.0 | 2.1% | 24.1% | yes |
| poor routing (b_K = 0.001, b_G = 0.002) | 95.0 | 43.7 | 0.4% | 6.8% | **no** |

The two cells have the same machine counts and the same Pmax: the potential capacity
aPmax/λ = 2.85×10⁶ fu/min exceeds D in both. In the poor case, 78.5% of chains are never captured
by either machine and fold spontaneously at low yield. The shortfall happens with both machines
more than 93% free, so it is not a saturation effect. Doubling both pools raises P_real only to
52.6, still short. A count of chaperone molecules, or a capacity bound built from counts, cannot
predict whether this demand is met.

## 1.6 Many clients and functions

With several proteins and functions, the requirement is A N ≥ D, where A_{ji} is the function j
delivered per functional copy of protein i. At steady state N_i = P_i/λ_i, and each P_i is a sum of
productive route fluxes that share the machines and the synthesis budget. Feasibility is a linear
system (E15). By Farkas' lemma it is infeasible exactly when there exist nonnegative prices on
resources and on demands under which every route and machine is worth no more than what it
consumes, while the total demand is worth more than all the resources (E18). Such a certificate
can be checked by hand. For one client and one function the system reduces to (E6).

We do not use a Hall-type subset count ("each group of clients fits in the machines it uses").
That test is exact only when each unit of flux occupies exactly one machine. A serial route
violates this, and a two-client example passes every subset count but is infeasible (EQUATIONS.md).

In a hypothetical two-client example (E1 needs G; E2 goes K→G; synthesis budget R = 60,000 aa/min;
machines lost at 0.02/min), the demand D0 = (1.8×10⁶, 0.9×10⁶) fu/min is feasible and 1.5·D0 is not
(verified certificate). At R = 48,000 aa/min, D0 is infeasible when machines must be maintained
and feasible if they were free. Machine upkeep alone decides the verdict.

## 1.7 Experimental mapping

A chemostat is the experimental form of this section: the dilution rate fixes λ's growth part and
the environment, so the demand is constant. Take an enzyme E that the medium makes essential, put
E under a titratable promoter, and choose E to be chaperone-dependent. (Candidate clients should
come from a retrieved GroEL-client study; none is cited here.)

| quantity | measure | gives |
|---|---|---|
| D | steady flux through E's reaction per cell (feed − residual substrate × dilution, cell count) | fu/min |
| N (functional) | E activity in lysate ÷ specific activity of pure native E; separately total E copies by targeted MS | functional copies, functional fraction |
| a | per-copy activity at in-vivo substrate level (permeabilized cells), ÷ functional copies | fu/(molecule·min) |
| λ | dilution rate + decay of E *activity* after a translation block (not abundance decay) | 1/min |
| P_real | λ · N_func at steady state; cross-check with pulse-labelled synthesis × functional fraction | molecules/min |
| M_m | DnaK, GroEL, protease complexes by quantitative proteomics | complexes/cell |

**Prediction.** As the promoter is titrated down, the culture keeps up at the set dilution rate
until a·N_func* falls to D, and washes out or shifts state below it. N_req = D/a is computed before
the titration. Washout at a clearly different N_func refutes Q = aN with constant a for this enzyme.
**Routing test.** Two strains with equal measured machine counts (e.g. a client variant with
altered chaperone recognition) should differ in functional fraction and in where they fail.

**What cannot be measured yet.** Dwell times c_{m,r} in vivo are unknown for any E. coli
machine-client pair, so the occupancy side of (E13), and therefore Pmax, cannot be computed from
data. The experiment measures P_real and Q_sus. Pmax remains a model ceiling.

**What fails without functional-activity measurement.** Proteomics counts all copies, including
misfolded, apo and aggregated ones, so a·N_total overstates Q. Abundance decay misses inactivation,
which understates λ. Synthesis rates count chains, not functional molecules, so yields are
invisible. Without activity, a cell with poor routing looks the same as one with good routing
(same counts, same synthesis), and the condition cannot be tested.

## 1.8 Scope

Assumed: constant D; Q linear in N; first-order loss; constant route yields and machine costs;
machines that fold unaided; linear resource use. Derived: the exact scalar condition and its
failure under feedback; turnover as the sole maintenance cost for ordinary catalysis; P_real ≤
Pmax with no fixed relation to machine counts; the Farkas feasibility test and the failure of
Hall counting. Not assumed: that the cell optimizes, that chaperones help, that any machine
saturates. Not treated here: demand that changes in time, growth-rate feedback on λ, and
aggregate sequestration of machines. Those come after the steady state.
