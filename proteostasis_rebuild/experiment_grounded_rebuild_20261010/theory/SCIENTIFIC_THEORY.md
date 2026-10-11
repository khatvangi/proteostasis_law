# Maintaining a functional proteome with finite, machine-specific capacity and condition-dependent routing

A minimal theory for E. coli. Argument first; the algebra is in [EQUATIONS.md](EQUATIONS.md),
the sources in [PRIOR_ART.md](PRIOR_ART.md), the data constraints in
`EXPERIMENT_CONSTRAINTS.tsv`, what was dropped in [CLAIM_DECISIONS.md](CLAIM_DECISIONS.md), and
what would kill it in [FALSIFICATION_REGISTER.md](FALSIFICATION_REGISTER.md). Every signed or
numeric statement labelled **[checked]** is recomputed by `checks.py` (10/10 pass). No parameter
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
2. **Overload need not look like collapse.** Without a specified positive feedback, overload
   produces either a smooth decline or a kink where one machine saturates. A sharp, history-dependent
   collapse requires a named feedback — aggregates holding chaperone, or seeded aggregation, or
   growth-dilution coupling. This is a prediction about *which* feedback, not a claim that cells
   have collapse thresholds.
3. **Maintenance is conditional.** There is no universal law here. For a given condition and a
   given class of disturbances there either is or is not a viable set the cell can stay inside.
   Some disturbances are not maintainable at all, and the theory says so rather than assuming
   otherwise.

What the theory does **not** assume: that the cell solves an optimization problem; that chaperones
are beneficial; that temperature is a load multiplier; that deleting capacity can only hurt; that
there is a death threshold.

## 2. Bookkeeping, and the two necessary conditions

For each protein species, chains exist as soluble non-native (U), native (N), aggregated (A),
machine-bound (C), or removed (D). Chains are born by synthesis, enter N co-translationally or via
folding, leave N by unfolding, enter A by aggregation or by unproductive release from a machine,
leave A by disaggregation, enter D by proteolysis, and every pool is diluted by growth. The
balance equations (EQUATIONS A1) conserve chains exactly: what is synthesized is either diluted in
some pool or degraded.

Setting dN/dt = 0 gives the **productive essential quota**: native flux must equal (μ + k_u) N.
Two things fall out. The quota does not vanish as growth slows — maintenance against spontaneous
unfolding remains — which is the same conclusion Santra 2017 reached from simulation, here from
the balance equation alone. And the quota bounds *productive* flux only; the rest of synthesis must
leave through aggregation, degradation or dilution.

Setting dC/dt = 0 gives **occupancy**: the bound fraction of a machine equals Σ(flux × service
time). This is Little's law, and it is exact **[checked: 0.65249218 against 0.65249218]**. It is
worth being blunt about what that buys. The right side — the pool — is measured (Schmidt 2016,
BW25113, glucose, exponential). The left side is not: **no in-vivo arrival flux or service time
exists for any E. coli machine–client pair.** Kerner's ">75% of GroEL capacity" is an occupancy
share for one client class, which is an estimate of the left side divided by the right, not of
either factor. So occupancy is a constraint with a known bound and an unknown load. It cannot be
turned into a margin, and this theory does not report one. The previous draft's ×25 margin is
withdrawn for a different reason as well (CLAIM_DECISIONS).

Rates are functions, not constants: r = r(T, cochaperone and nucleotide state, client identity,
free machine pools). Free pools are themselves state variables through the conservation of each
machine. Temperature enters every rate separately — folding, unfolding, aggregation, binding — so
no monotone temperature-to-load map is assumed. Kanemori's protease plating data, non-monotone in
temperature, already rules that map out.

Two different codon variables appear, and conflating them was an error in the earlier work. Codon
choice sets *when* a chain is released and how fast it elongates, which changes the co-translational
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

There is a third way maintenance can fail with a non-empty feasible region: a good steady state can
exist while the cell sits in another basin (§5).

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

Three regimes, and which one occurs depends on a feedback that must be named.

**Graceful failure.** With no feedback, the aggregate steady state is unique and linear in load
**[checked]**. Damage rises smoothly; coupling growth to native function gives a stable
slower-growing state. Performance degrades continuously until some quota cannot be met even at
zero growth.

**Single bottleneck.** With finite pools, the machine with the smallest pool-to-load ratio
saturates first, its throughput caps at pool/service-time, and the excess spills to aggregation and
degradation. The yield-versus-load curve kinks but stays continuous and single-valued. Note which
machine binds is set by service time, which is unmeasured — so abundance alone cannot name the
bottleneck, and claims that it can (including in the earlier network work, where ClpB appeared
binding largely because its pool is small) are not supported.

**Runaway, only with specified positive feedback.** Let each aggregated chain hold β machine units.
Then aggregate clearance weakens as aggregate grows, and the steady-state equation becomes convex:
either two states (a low stable one and an unstable threshold) or none. The saddle-node is at
L_crit = (μ+δ)a²/(4 b k_a) with b = k_c β **[checked: two states at 0.9 L_crit, none at 1.1]**.
Above L_crit the cell runs away; below it, a large enough perturbation crosses the threshold and
runs away too. That history dependence is the signature, and it is absent when β = 0.

Two other feedbacks give the same structure and need their own measurements: seeded aggregation,
and aggregates slowing growth so that dilution removes less aggregate.

No death cutoff is imposed anywhere. "Failure" means only that no state satisfies the quota.
Whether that means arrest, filamentation or death needs more equations and more data. And since β
has never been measured in vivo, **no collapse threshold is claimed, and certainly not a universal
one or several.**

## 6. When is maintenance possible?

Define the viable set as the states where every essential protein is above its required native
level. A condition and a class of disturbances is *maintainable* if some closed subset of the
viable set is forward-invariant under every admissible disturbance — formally, if on the boundary
of that subset the productive flux still covers the quota (Nagumo's condition; EQUATIONS A7). The
weaker local version is an asymptotically stable steady state inside the viable set, with the
disturbed state in its basin.

This is deliberately a conditional statement. If a disturbance raises unfolding or aggregation
enough that no state satisfies the quota, no invariant subset exists and maintenance fails. The
theory asserts no universal robustness. Vorderwülbecke and Genevaux make this concrete: Δtig ΔdnaK
is maintainable at 23 °C, not at 30 °C, and maintainable again at 30 °C with GroEL/ES
overproduced. The same genotype, three conditions, two answers. Also missing from §5 and owed
before any quantitative claim: σ32-driven chaperone induction makes the pools state variables with
negative feedback, which enlarges the viable set.

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
and bifurcation structure tied to one named feedback coefficient rather than asserted (§5). Of
these, the first is likely a presentation advance rather than new physics, and all three still owe
a targeted citation search before any novelty is claimed in a manuscript.

## 8. Predictions

Three, pre-registered in [FALSIFICATION_REGISTER.md](FALSIFICATION_REGISTER.md) with their
discriminating data and their minimum readouts. In brief: the HtpG disposal sign should flip with
genetic background in a direction the rule specifies; genetic-interaction signs should follow
serial versus parallel route structure; and hysteresis under a graded load up-and-down sweep should
appear only if chaperone sequestration by aggregates is real.

**No manuscript is to be drafted from this.** The theory is a hypothesis under test; the register
lists what must survive first.
