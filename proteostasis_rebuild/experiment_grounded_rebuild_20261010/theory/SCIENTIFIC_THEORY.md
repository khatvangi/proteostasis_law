# Maintaining a functional proteome with finite, machine-specific capacity and condition-dependent routing

A minimal theory for E. coli. Argument first; the algebra is in [EQUATIONS.md](EQUATIONS.md),
the sources in [PRIOR_ART.md](PRIOR_ART.md), the data constraints in
`EXPERIMENT_CONSTRAINTS.tsv`, what was dropped in [CLAIM_DECISIONS.md](CLAIM_DECISIONS.md), and
what would kill it in [FALSIFICATION_REGISTER.md](FALSIFICATION_REGISTER.md). Every signed or
numeric statement labelled **[checked]** is recomputed by `checks.py` (22/22 pass). Sections
2, 5 and 6 were corrected after an independent audit; see [CORRECTIONS.md](CORRECTIONS.md). No parameter
here is an E. coli estimate.

## 1. The claim

A cell keeps its essential proteins folded by sending newly made chains down competing kinetic
routes — folding alone, folding on a machine, aggregating, being degraded, being diluted away.
Each machine is a finite pool, and each chain that occupies one holds it for some time. Two
conditions must therefore hold in any steady state:

- **the productive quota** — for every essential protein, the flux arriving in the native state
  must cover what dilution and spontaneous unfolding remove;
- **machine occupancy** — for every machine, the sum over clients of (arrival flux × holding time)
  cannot exceed the pool.

Both are consequences of bookkeeping, not hypotheses. The scientific content of the theory is
what they do *not* determine, and that is the part this project is about:

**Neither condition says the cell will meet the quota, because the cell does not choose its
routing.** Routing is set by rate competition at the branch point. The set of flux patterns the
capacity permits — call it the feasible region — can be wide open while the pattern actually
realized misses the quota entirely. Capacity and routing are two different failure modes, and
distinguishing them is the whole point.

Three consequences follow, and they are the testable part.

1. **Increasing a chaperone can be harmful.** A machine whose bound clients mostly exit to
   degradation competes with productive routes for the same chains. Whether adding more of it
   raises or lowers native yield is fixed by a sign rule with measurable ingredients (§4).
2. **Overload need not look like collapse.** Without a specified positive feedback, the aggregate
   response to overload is smooth, or kinks where one machine saturates. A history-dependent switch
   needs a named feedback, such as aggregates holding chaperone. In the model analyzed it also needs
   chaperone capture that outpaces spontaneous aggregation. Even then the switch is between two
   *finite* states: a hysteresis loop, not a runaway. This is a prediction about *which* feedback
   and *which* condition, not a claim that cells have collapse thresholds.
3. **Maintenance is conditional.** There is no universal law here. For a given condition and a
   given class of disturbances there either is or is not a viable set the cell can stay inside.
   Some disturbances are not maintainable at all, and the theory says so rather than assuming
   otherwise.

What the theory does **not** assume: that the cell solves an optimization problem; that chaperones
are beneficial; that temperature is a load multiplier; that deleting capacity can only hurt; that
there is a death threshold.

## 2. Bookkeeping, and the two necessary conditions

For each protein species, chains exist as soluble non-native (U), native (N), aggregated (A),
machine-bound (C), or removed (D). Chains are born by synthesis and enter N co-translationally or
by folding. They leave N by unfolding **or by regulated proteolysis of the folded protein**. They
enter A by aggregation or by unproductive release from a machine, and leave A by disaggregation or
degradation. They enter D by proteolysis from any of U, C, N or A, and every pool is diluted by
growth. The balance equations (EQUATIONS A1) conserve chains exactly. Whatever is synthesized is
either diluted in some pool or degraded from one **[checked: drift below 2·10⁻¹³]**.

Machines obey their own conservation law, and here the earlier version was physically wrong. It
had each aggregated chain hold a fixed β units of chaperone however much aggregate there was. Past
the point where aggregate sites outnumber the pool, that forces the free chaperone negative. In the
corrected law, β is a *capacity*. Aggregates hold β units per chain only while chaperone is
available. Once the pool left after client binding is used up, they hold that remainder and
nothing more, and the free pool is exactly zero, never negative. With that cap, the model has a
property that can be proved rather than hoped for. Starting from nonnegative amounts, every pool
stays nonnegative for all time, and bound chaperone never exceeds the pool (EQUATIONS,
Proposition A1) **[checked; with the uncapped form, the same instance drives bound chaperone to
−0.12]**.

Setting dN/dt = 0 gives the **native quota**. Native flux must cover dilution, repair of
spontaneous unfolding and replacement of degraded native protein: (μ + k_u + k_dN)N^req. Three
things follow. First, the quota does not vanish as growth slows. It tends to
(k_u + k_dN)N^req(0), the same conclusion Santra 2017 reached by simulation, here from the balance
equation alone. Second, reaching slow growth is not automatic. If synthesis is strictly
proportional to growth rate, then near zero growth the productive flux goes to zero while the
quota does not, so very slow growth is never viable. Smooth survival down to μ → 0 needs a basal
synthesis term large enough to pay the maintenance quota. Without one, the range of viable growth
rates can close abruptly at a *positive* growth rate as productive yield falls **[checked]**.
Third, the quota bounds *productive* flux only. The rest of synthesis must leave through
aggregation, degradation or dilution.

Setting dC/dt = 0 gives **occupancy**: the bound fraction of a machine equals Σ(flux × service
time). This is Little's law, and it is exact **[checked: 0.65249218 against 0.65249218]**. It is
worth being blunt about what that buys. The right side, the pool, is measured (Schmidt 2016,
BW25113, glucose, exponential). The left side is not: **no in-vivo arrival flux or service time
exists for any E. coli machine–client pair.** Kerner's ">75% of GroEL capacity" is an occupancy
share for one client class, which estimates the left side divided by the right, not either
factor. So occupancy is a constraint with a known bound and an unknown load. It cannot be turned
into a margin, and this theory does not report one. The previous draft's ×25 margin is withdrawn
for a different reason as well (CLAIM_DECISIONS).

Rates are functions, not constants: r = r(T, cochaperone and nucleotide state, client identity,
free machine pools). Free pools are themselves state variables through the conservation of each
machine. Temperature enters every rate separately (folding, unfolding, aggregation, binding), so
no monotone temperature-to-load map is assumed. Kanemori's protease plating data, non-monotone in
temperature, already rules that map out.

Two different codon variables appear, and conflating them was an error in the earlier work. Codon
choice sets *when* a chain is released and how fast it elongates. That changes the co-translational
partition and the competition at the branch point. Mass spectrometry measures something else: the
standing frequency of amino-acid substitutions. A substitution makes a different protein with its
own rates; a synonymous change does not. The proteomic error tables constrain only the second
variable. Nothing in them bears on folding, and the theory does not let them.

## 3. Necessary, not sufficient — and why that is the useful part

Feasibility is necessary: any steady state that keeps the essential proteins folded yields fluxes
obeying the quota and occupancy, so the feasible region is non-empty (EQUATIONS, Theorem 1).

It is not sufficient. A worked instance **[checked: feasible_not_sufficient]**: one client class,
quota 0.8 of synthesis, ample pools. The best routing the capacity permits delivers native flux
1.0. The routing the rates actually produce delivers 0.28, with neither machine above 16%
occupancy. The quota fails with no capacity shortage anywhere. A capacity audit would have
declared the cell healthy.

This gap has a sharp experimental consequence **[checked:
deletion_can_help_only_kinetically]**. Removing a route can only shrink the feasible region, so
the best attainable point cannot improve: 1.0 before, 1.0 after. In the same instance, removing the
route raises the *realized* native flux from 0.28 to 0.89. So a model in which the cell attains the
best point of its feasible region — the FoldME family — cannot reproduce a growth improvement on
deleting a chaperone unless a cost term is added, whereas a kinetic-routing model reproduces it
with nothing added. Fauvet's ΔhtpG suppression of ΔdnaKJ is exactly that observation. The caveat is
real: FoldME does carry a proteome-cost term, so this discriminates the plain optimizing model and
needs care against the costed one.

There is a third way maintenance can fail with a non-empty feasible region. Where the overload
model has a three-state window, a good steady state can exist while the cell sits at the finite
high-aggregate state (§5).

## 4. The solved example: GroEL productive, HtpG disposal, spontaneous folding

Take one client with three ways out of the non-native state — fold alone, bind GroEL (which mostly
returns it native), bind HtpG (which mostly hands it to a protease) — plus aggregation,
disaggregation and dilution. The steady state is solvable in closed form (EQUATIONS A8).

**The sign rule.** More HtpG raises native yield only if

    (1−q)(k_a^eff + μ)  >  q (k_f + p_G k_G G),

where q is the probability that an HtpG-bound chain is destroyed and k_a^eff is the aggregation
rate discounted by the fraction that disaggregation recovers **[checked symbolically; ratio to the
exact derivative positive over 20,000 random parameter draws]**. Read plainly: a disposal route
helps when the chain's alternative is being lost anyway, and hurts when the chain had good
productive alternatives. Both signs occur in the same model: ∂Y/∂H = −0.046 when folding and GroEL
are strong, +0.103 when aggregation dominates **[checked]**.

Three things sharpen this.

- **More HtpG always lowers the aggregate burden** **[checked: ∂A/∂H < 0 over 20,000 draws]**. So
  the same intervention can reduce aggregation and reduce function at once. Which one growth
  follows is a separate, measurable question, and the theory deliberately leaves it open.
- **Better disaggregation makes disposal worse.** Efficient ClpB lowers k_a^eff, so losing ClpB
  moves HtpG toward being helpful. That is a signed, testable interaction.
- **Saturation does not change the sign.** With a finite HtpG pool, free HtpG rises with total
  HtpG, so the sign of the effect of total HtpG matches the free-HtpG rule; yield falls
  monotonically over a 400-fold range in the harm regime **[checked]**. With several clients
  sharing one pool this transfer is not proven, and is listed as an open point.

**GroEL is not automatically safe either.** With all GroEL cycles productive, ∂Y/∂G > 0
analytically. But give GroEL an unproductive exit to aggregate and a client that folds fast alone,
and ∂Y/∂G = −0.058 **[checked]**. So "chaperones are beneficial" is a statement about parameters,
not a property of chaperones — and that is before the proteome cost of making them, which FoldME
treats and this model does not.

## 5. What overload does

Three regimes are possible. Which one occurs depends on a feedback that has to be named and on
one further condition. All of this is worked out in a **reduced one-variable model**: one load,
one clearance chaperone, chaperone held by clients left out, net aggregate removal lumped into δ
(EQUATIONS A6). Its conditions belong to that model, not to the full network.

**Smooth response.** With no feedback, the aggregate steady state is unique, stable and linear in
load **[checked]**. Damage rises smoothly with load. This does *not* by itself mean that growth
declines gracefully. Whether coupling growth to native function gives a stable slower-growing
state depends on how growth rate depends on the proteome (the closure μ(x)), which is not
specified. If synthesis is strictly proportional to growth rate, the viable range of growth rates
can close abruptly at a nonzero rate (§2). So "graceful failure" is a possibility conditional on
basal synthesis and a stable closure, not a theorem.

**Single bottleneck.** With one finite pool and one client class, the steady state stays unique.
The machine's throughput caps at pool/service-time, and the excess spills to aggregation and
degradation, so the yield-versus-load curve kinks. With several pools sharing clients, uniqueness
has not been proven. Which machine binds first is set by service time, which is unmeasured. So
abundance alone cannot name the bottleneck, and claims that it can are not supported. That
includes the earlier network work, where ClpB appeared binding largely because its pool is small.

**A bounded hysteresis loop, only with sequestration, and only if capture is fast enough.** Let
aggregates hold chaperone, up to β units per aggregated chain, with total sequestration capped by
the pool. Clearance then weakens as aggregate grows until the pool is fully held. After that it
cannot weaken further, and aggregate is removed only by dilution and net removal at rate μ + δ.
The exact result (EQUATIONS, Theorem A6):

- If μ + δ > 0 there is **always** at least one steady state, and every trajectory stays bounded.
  A **finite, stable high-aggregate state** A_hi = k_a L/((k_a + μ)(μ + δ)) exists at every load
  above L_sat. It is the aggregate level of a cell whose clearance chaperone is entirely held.
- **Coexistence needs a second condition.** It requires k_c M^tot > k_a + μ: at full pool,
  chaperone capture of a non-native chain must outpace spontaneous aggregation plus dilution.
  When it holds, the load window L_sat < L < L_crit has three states: low (stable), middle
  (unstable) and high (stable). Below the window only the low state exists, and above it only the
  high one. When it fails, there is one state at every load and no fold at all
  **[checked: 1 / 3 / 1 states at 0.9 / 0.99 / 1.1 L_crit in one parameter set; 1 state at all 300
  loads in a set with the same four lumped numbers but the condition reversed]**.
- Ramp the load up and the aggregate jumps to the high state at L_crit. Ramp it down and the
  aggregate stays high until L_sat. That loop is the signature **[checked: jumps at 2.263 and 2.192
  against 2.25 and 2.1875]**. Its relative width is ((x − 1)/(x + 1))², with x = k_c M^tot/(k_a + μ).
  The loop can be narrow: 2.8 % of L_crit at x = 1.4.
- **Unbounded accumulation happens only with no removal at all** (μ + δ = 0). A runaway would
  need a further, separately specified coupling that drives growth to zero as aggregate builds,
  with no aggregate degradation. That coupling is not modeled.

The earlier version of this section said there were "two states or none" and that the cell "runs
away" above L_crit. Both claims are withdrawn. They came from letting free chaperone go negative
(CORRECTIONS.md).

Two other feedbacks might give their own folds: seeded aggregation, and aggregates slowing growth
so that dilution removes less aggregate. Neither has been analyzed with a capped pool, and each
needs its own measurement.

No death cutoff is imposed anywhere. "Failure" means only that no state satisfies the quota. The
high state is finite, and whether it breaks the quota depends on things not specified here.
Whether failure means arrest, filamentation or death needs more equations and more data. Neither
β nor k_c has been measured in vivo, so **no collapse threshold is claimed, and certainly not a
universal one or several.**

## 6. When is maintenance possible?

Define the viable set as the states in which every essential protein is at or above its required
native level. A condition and a class of disturbances is *maintainable* if some closed subset of
the viable set stays forward-invariant under every admissible disturbance. For a smooth enough
model this holds exactly when, at every boundary point of that subset and for every disturbance,
the dynamics do not point outward (Nagumo's condition; EQUATIONS A7). The weaker, local version is
an asymptotically stable steady state inside the viable set, with the disturbed state in its
basin.

On a face where a quota is just met, "not pointing outward" means more than productive flux
covering the quota. The required level itself moves when growth rate moves. So the condition is

    productive flux  ≥  (μ + k_u + k_dN) N^req(μ)  +  (dN^req/dμ) · (dμ/dt).

The earlier statement left out native degradation and the last term. It fails in exactly the
case that matters, a shift in growth rate. The worked counterexample meets the old condition with
zero slack and still leaves the viable set at once **[checked]**. With the derivative term
included, a productive flux large enough to meet the condition along the whole face keeps the
state inside **[checked]**. This condition is necessary on quota faces. It is sufficient only
together with the tangency conditions on every other face of the chosen set.

This is deliberately conditional. If a disturbance raises unfolding, native degradation or
aggregation enough that no state satisfies the quota, no invariant subset exists and maintenance
fails. The theory asserts no universal robustness. In particular it does **not** assert that
decline is always smooth. Under growth-proportional synthesis, viability can be lost at a positive
growth rate (§2). Inside the A6 window both the low and the high state are locally stable, and
which one the cell occupies depends on its history.

Vorderwülbecke and Genevaux make the conditional nature concrete. Δtig ΔdnaK is maintainable at
23 °C but not at 30 °C, and maintainable again at 30 °C with GroEL/ES overproduced: the same
genotype, three conditions, two answers. Still missing, and owed before any quantitative claim:
σ32-driven chaperone induction makes the pools state variables with negative feedback, and that
enlarges the viable set.

## 7. Honest position relative to prior art

The finite-capacity chaperone network is **not new**. FoldEco already models synthesis, folding,
misfolding, nucleated aggregation, TF, DnaK/DnaJ/GrpE, GroEL/ES, ClpB and Lon by mass action with
finite pools, and already has degradation competing with folding. Santra 2017 adds growth coupling,
client sorting by misfolded-state properties, and occupancy near saturation at fast growth. FoldME
adds genome-scale temperature-dependent folding and the ribosome-versus-chaperone proteome
trade-off. Anyone presenting a capacity network as original here would be wrong, and the earlier
drafts in this project came close to doing so.

What is left that looks genuinely testable: the explicit separation of a necessary occupancy bound
from realized routing and the signed experimental consequence of that gap (§3); an unfavourable
disposal route with an analytic sign rule and the ClpB and protease interactions it predicts (§4);
and a hysteresis loop tied to one named feedback coefficient and one stated capture condition,
with a finite upper branch, rather than an asserted collapse (§5). Of
these, the first is likely a presentation advance rather than new physics, and all three still owe
a targeted citation search before any novelty is claimed in a manuscript.

## 8. Predictions

Three, pre-registered in [FALSIFICATION_REGISTER.md](FALSIFICATION_REGISTER.md) with their
discriminating data and their minimum readouts. In brief: the HtpG disposal sign should flip with
genetic background in a direction the rule specifies; genetic-interaction signs should follow
serial versus parallel route structure; and a hysteresis loop under a graded load up-and-down
sweep, with a finite high-aggregate branch, should appear only if chaperone sequestration by
aggregates is real *and* capture outpaces spontaneous aggregation plus dilution.

**No manuscript is to be drafted from this.** The theory is a hypothesis under test; the register
lists what must survive first.
