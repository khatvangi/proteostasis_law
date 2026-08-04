# Status

## Phase 3: the fold is derived, not sampled

On 2026-08-02 the collapse boundary was characterised analytically. Because `j`
enters `du/dt` and nowhere else, the aggregate nullcline is a fixed curve, and
mass balance plus one determinant-preserving row operation gives the identity
`det J = -(grad R x grad G)` for total removal `R`. A saddle-node is therefore
exactly a constrained critical point of `R` on the nullcline, so

```
j_crit = R(u*, a*)          exact, no continuation sweep
```

**The collapse boundary is where total removal stops responding to burden along
the aggregate nullcline.** (The proven statement is constrained *critical point*;
an earlier "constrained maximum" gloss is withdrawn — see D011 below.) Statement,
proof, verification and limits are in `theory/FOLD_THEOREM.md`; decisions
D007-D016 record what changed.

Verified against the Phase 1 run root: the identity holds to a median relative
error of **1.436e-07** (the central-difference floor), the parallelism residual
correlates with the recorded leading eigenvalue at **+0.9987** — showing the
residual is bracket tolerance and not a failure of the identity — and the 2x2
solver reproduces the continuation-derived folds to **6.652e-07**. `phi` rebuilds
from first principles at all **2884** folds, median absolute error **1.3e-13**.

**What sets phi.** At collapse the machinery runs at `s_ref` **0.175**, `s_u`
**0.155**, `s_a` **0.056** — roughly 6-18 % of V_max. Collapse is *not* capacity
exhaustion; superlinear nucleation overtakes sublinear saturating removal long
before removal saturates. Counterfactually the saturation deficit accounts for
**35.8 %** of the shortfall against **12.6 %** for sequestration. This supplies
the mechanism behind P2 and shows `removalCeiling`, while a correct bound, is
about 13x too loose to predict anything alone.

**One earlier reading withdrawn (D009).** A within-draw comparison over the
allocation `chi` alone gave a between/within variance ratio of 9.6, read as `phi`
being a load-invariant network constant. The properly nested design — kinetic
draws crossed with a `(nu, chi)` grid, 446 folds solved — gives **5.9**, with
per-network spreads up to **13.6x** when both load coordinates sweep, overlapping
the **8.86x** between-network spread. `phi` is network-characteristic but not
load-invariant. Holding one load coordinate gives 1.59-1.80x.

**No empirical claim.** No organism data entered Phase 3 either. The sharpened
prediction — collapse at `s_a` far below saturation, against the
capacity-exhaustion alternative — is an untested empirical hypothesis, and it
requires growth rate to be held fixed because growth rate sets `nu`.

**Growth dilution, and the coupling the theory had not stated.** Cell division
was absent from the model, and for most proteins dilution outpaces proteolysis.
`scripts/phase3/dilution.py` adds it without touching the frozen upstream model.
The theorem survives for any dilution law — `mu -> 0` reduces to the frozen model
exactly (0.0), the identity holds at 1.2e-10 (constant `mu`) and 4.7e-10
(burden-dependent), and for constant `mu` the diluted Jacobian is exactly
`J - mu.I`, so the saddle-node condition says `mu` is an eigenvalue of the
undiluted Jacobian.

The **removal ceiling does not survive**: `R_tot = R + mu.(u+a)` is unbounded in
burden, so A8 is false once cells divide. Under *constant* dilution the fold is
destroyed above a threshold — continuation gives `j_crit` rising 0.1542 → 0.2456
as `mu` goes 0 → 0.08 while `a*` diverges 0.265 → 0.750, and by `mu = 0.10` no
fold exists. Growth feedback restores it at every rate tested.

> **The collapse boundary of a dividing cell exists because burden slows growth.**

An experiment that lets growth rate float is therefore not holding disposal
capacity fixed, because growth rate is part of disposal.

**One gloss withdrawn (D011).** The proven statement is constrained *critical
point*, not *maximum*. `R` rises monotonically along the first-root branch and
the solved fold has `dG/da > 0`, so the critical point sits past the nullcline's
turning point. No number changes — everything solves `det J = 0` directly.

**Division makes the system bistable, and uniqueness fails (D012, D013).**
Enumerating the whole nullcline rather than one branch settles the uniqueness
question with a split verdict. Without dilution the curve closes (152 lower, 152
upper) and carries exactly **one** critical point. With dilution at `mu = 0.04`
there is a **second** genuine saddle-node at `u = 0.1585, a = 1.9835`
(`|G| = 2.8e-17`, `|det J| = 1.6e-17`, eigenvalues `-1.083, 0`) whose critical
influx **0.1585 lies below** the first at 0.1950.

They are the two folds of an S-curve, because dilution makes the high-burden
state bounded rather than divergent. Sweeping influx up from zero burden and back
down: the healthy branch survives to `j = 0.194` and jumps at 0.196 to a bounded
state (`u = 0.079, a = 3.94`); the damaged branch survives back to `j = 0.160`
and recovers at 0.158. The bistable window **[0.160, 0.194]** lies inside the two
computed critical points. **Losing proteostasis under division is a transition to
a persistent aggregate-laden state, not a divergence, and recovery requires
lowering influx below the level that triggered it.** This does not reinstate
`notes/REJECTED_COMPONENTS.md` item 7 — it comes from a model feature neither
that claim nor Phase 2 §2.1 contained — and it was found at one parameter point,
not surveyed.

**That parameter point is unphysical, in a specific and measured way.** It uses
constant dilution (`k_mu = inf`), under which growth rate cannot respond to
burden, so the same regime predicts **exactly zero** growth-rate loss at any
aggregate load. That contradicts the one dosage-resolved measurement the project
holds (D015: 3.2 % loss at <0.1 % misfolded) and the observation the first
post-diction tried to explain (D026/D028: a 1.2–1.8 % aggregate-attributable
growth deficit — not the ">30 %" an earlier reading of that abstract took it to
be). The regime that
produces the bistability is the regime that gets the measured quantity wrong.
Under the physiological laws it does not survive as stated — linear arrest gives
no bounded high-burden state, hyperbolic feedback is monostable in four of six
settings. Pinned by test so it cannot be edited away.

**A margin that survives division (D014).** `j_crit = C_enz . phi_enz / (1-delta)`
exactly (identity to 1.6e-16), with `phi_enz` the enzymatic capacity in use at
collapse and `delta` the share of disposal done by division. `phi_enz` is nearly
invariant to dilution (0.125–0.134) while `delta` carries the variation, so
division multiplies tolerable influx by `1/(1-delta)` without changing the
enzymatic condition. Across 25 draws, 23 lose their boundary under constant
dilution; thresholds span 3.3 decades and `delta` at the threshold has median
0.35. The `(phi_enz, delta)` decomposition holds under both the hyperbolic and the
linear growth-burden form; the bistability above does NOT — see D015 immediately
below, which is a retraction, not a confirmation.

**Calibrated against a measured growth-burden relation (D015).** The dilution
laws above were guesses. A real dosage-resolved measurement is now the anchor:
Geiler-Samerotte et al. 2011 PNAS, doi:10.1073/pnas.1017570108 (PMID 21187411),
via PubMed — 3.2 % growth-rate reduction at <0.1 % of total cellular protein
misfolded, in **yeast**. Slope 32 per unit proteome fraction (a lower bound),
linear arrest at 0.03125 (an upper bound). The two conversions needed to use it —
the quality-control proteome share and the refolding turnover — were **not**
obtained and are swept, not assumed.

Three results. (1) The prior guesses were in the right regime: calibrated `k_mu`
runs 0.31–6.25 and the guessed 0.5 and 2.0 both lie inside. (2) **A boundary
survives in 30 of 30 calibrated cells**, `phi_enz` confined to 0.072–0.147 — so
the "constant dilution destroys the boundary" behaviour is an artifact of
omitting the measured coupling, not a prediction. (3) **The bistability above
does NOT survive and is form-dependent**: under the linear (measured-shape) law
the down-branch fails to settle at 0/8 influx values — a runaway at `a ~ 10^4`,
not an attractor — while the hyperbolic form settles on both branches with window
(0.17, 0.19). Linear arrest sets `mu = 0` exactly beyond `k_mu`, switching
dilution off; the hyperbolic form only approaches zero.

**The measurement cannot adjudicate**, constraining only the slope three decades
below the arrest burden. So whether collapse in a dividing cell is a reversible
switch or an irreversible runaway turns on an unmeasured property — is growth
arrest under misfolding burden complete or asymptotic? That is now the
highest-value measurement the theory asks for, and the bistability claim is
quoted with its growth law attached rather than as a consequence of division.

**The Pareto layer now exists (D017).** `theory/PARETO_GEOMETRY.md` asserted a
trade-off surface cut by the proteostasis condition, and no script had ever
computed one. `scripts/phase3/pareto.py` supplies the minimal version: strategy
`(alpha, R)` for accuracy and quality-control investment, both paid out of one
proteome, cut by the derived constraint. **The throughput optimum sits exactly on
the feasibility boundary** (`j/j_crit = 1.000000`) — a grid gave 0.8975 and that
was discretisation. But along the non-dominated front `j/j_crit` runs
**0.227–0.965**, so "strategies sit near the boundary" holds only at the
throughput-maximising end. A deterministic maximiser has zero margin, so any
observed margin needs a mechanism outside this layer — worth noting against the
superseded envelope-paper's order-of-magnitude margin claim.

**Empirical contact is specified, not attempted (D018).**
`empirical/GATE4_PROPOSAL.md` states what a test of the theorem requires, with no
outcome value read. Its load-bearing conclusion: **the staged data cannot test the
central prediction**, because H1 concerns quality-control flux relative to its own
maximum and nothing staged measures a saturation state. Directional predictor
comparisons remain executable but are not a test of the theorem.

**Arrest-shape literature search returned no usable measurement.** The question
D015 identified as highest-value — is growth arrest under misfolding burden
complete or asymptotic — was searched for via PubMed. Samhita et al. 2025
(doi:10.1093/molbev/msaf312) establish that 10–50x wild-type mistranslation
remains viable in E. coli but measure costs only in the weak regime; Melnikov
et al. 2020 (doi:10.1073/pnas.2003132117) report norvaline intolerance
qualitatively. No dose-resolved growth-versus-quantified-burden curve into the
arrest regime was found, and two relevant full texts were not retrievable through
that route, so a targeted PDF read could still settle it.

**The Gate 4 instrument was tested for and exists (D019).** The blocker D018
identified — a clearance-flux readout with an internally determinable maximum —
resolves positively. Proteolytic queueing at ClpXP is routine and engineered:
Jadhav et al. 2025 (doi:10.1021/acssynbio.4c00612) localise the queueing
bottleneck to **ClpX** by overexpressing each component in turn, and Ogle & Mather
2016 (doi:10.1088/1478-3975/13/2/025002) show inter-substrate correlations peak
at the queueing balance point — an internal saturation-state signature. Titratable
substrate, reachable saturation, and a maximum determinable in the same cells.

**But it forces an arm substitution.** The accessible arm is ClpXP degradation of
soluble substrate, `s_u` (median 0.155), not aggregate clearance `s_a` (0.056).
H1 is restated as H1' against `s_u` with K1 reset. `s_u` is the **less extreme**
arm, so the executable test is weaker than the headline figure. The chaperone arm
has no comparable handle and stays untestable.

**Preregistering Gate 4 killed its own design, then produced a better one
(D020, D021).** Deriving K1's threshold from the model rather than choosing it
showed the saturation test is underpowered: predicted `s_u` at the boundary in
dividing cells has median 0.119 but **p99 0.897**, only 0.103 from H0's
`s_u -> 1`, with 17.8 % of draws above 0.5. The "6-18 % of V_max" headline is a
**median** over a distribution covering nearly the whole interval, and the limit
is the theory's own parameter uncertainty, not the assay. `s_u` is demoted to a
reported descriptive quantity.

The replacement primary outcome is a **parameter-free exponent**:
`tau ~ (j_crit - j)^(-1/2)`, since one eigenvalue passes through zero at a
saddle-node. Fitted by continuation outward from the exactly-known fold state,
the slope is **0.5077** median with **r2 = 1.0000**, unchanged by dilution
(0.5134 undiluted; 0.5080 and 0.5054 at calibrated `mu0` 0.05 and 0.10), and
86.4 % of networks fall within 0.05 of 0.5. Against `s_u`'s [0,1] spread over the
same box, this spans 0.497-0.513. Restated for execution as **`tau^-2` linear in
dose**, whose x-intercept locates the boundary without needing `j_crit` in
advance.

Limits carried: 22 of 42 ladders converged, so the sample is small and selected;
and the 1/2 exponent is generic to saddle-nodes, so it supports "the boundary is
a saddle-node" rather than selecting this model. What remains unspecified before
any hash-freeze is listed in `empirical/GATE4_PROPOSAL.md` §11.

**The exponent is a ruler, not a test — and the test now exists (D022, D023).**
D021's exponent is robust *because* it is generic to saddle-nodes, so confirming
it selects this model over no alternative; D020's saturation fraction
discriminates but cannot be measured. The resolution: `tau^-2` regressed on dose
locates `j_crit` as an x-intercept with a CI — an instrument — and the test is
whether that boundary **moves**.

**H3:** raising the load of *perfectly-folding* protein lowers the tolerable
mistranslation dose, though it causes no damage. This follows from rescue
capacity being shared — `nu` enters only the denominator of the free-chaperone
balance. An independent-handling model predicts no shift. Over a 100x
nascent-load ladder: direction correct in **67 of 68 (98.5 %)**, monotone in
98.5 %, median shift **1.22x**, p90 3.51. Ladder chosen by sweep; wider ladders
buy effect at the cost of direction consistency.

**And it is invalid in batch culture.** Gratuitous expression raises `nu` and
lowers `mu`, and `mu` is disposal — so a batch experiment would show a shift
under *both* competing models. H3 requires externally fixed growth rate
(chemostat/turbidostat). K6 voids the gate on between-arm growth differences.
Pinned by test so it cannot be edited away.

**The theorem generalises to n states (D024).** The "exact about a toy"
objection is answered: with state `(u,a,c,...)`, influx in one equation and mass
balance intact, `det J = -det[grad R; grad G; grad C]`, so a saddle-node is a
constrained critical point of R on the intersection of the non-influx nullclines.
Verified on a three-state system under sigma-32-style chaperone control —
relative error **0.000e+00** unregulated, ~2.5e-11 regulated, and `sigma0 -> 0`
reproduces the frozen model exactly. Extending the model no longer requires
re-deriving the boundary.

**Regulation does NOT rescue the predictions (D025).** The hypothesis was that a
controlled cell sits where its controller puts it, collapsing the spread that
made `s_u` untestable. **Refuted** — p5-p95 width goes 0.8904 -> **0.9677**, it
widens. One tentative observation survives and cuts against the headline: the
regulated median `s_u` is 0.323 vs 0.169, so control pushes collapse CLOSER to
saturation, partway toward the capacity-exhaustion picture. 14 of 30 converged;
not quotable yet.

**Standing of the theory.** Two attempts to sharpen the quantitative predictions
have failed (calibration D015, regulation D025) while the structural core has
survived every extension (dilution, an added controlled state). This is a
**structural theory, not a predictive one**: it says exactly where the boundary is
given the parameters and that this holds for the whole model class; it does not
predict a number without measured parameters. The defensible claim is that the
naive capacity bound is wrong by ~an order of magnitude and the boundary is a
computable constrained critical point — not that collapse occurs at any
particular fraction of V_max.

**First post-diction attempted, and it FAILED (D026).** Lindner et al. 2008
(doi:10.1073/pnas.0708931105) report that unstressed E. coli accumulate
aggregates in old-pole cells, losing >30% of reproductive ability, while new-pole
progeny rejuvenate. The logical point stands — **rejuvenation is only coherent in
a bistable system**, since in a monostable one a daughter with less aggregate
relaxes straight back. But the model does not supply that bistability where it
matters: constant dilution is bistable yet predicts **zero** reproductive loss by
construction; hyperbolic feedback is monostable in four of six settings and where
bistable predicts 48-95% loss, more severe than observed; linear arrest gives no
bounded high state at all.

Reporting a success would have meant quoting the constant-dilution regime — the
one that gets the measured quantity wrong. The negative is pinned by test.

**The failure names the next mechanism: spatial sequestration.** The observation
is about aggregates being localised at a pole and segregated asymmetrically; this
model is well-mixed. Sequestration removes aggregate from the reactive pool,
changing kinetics rather than bookkeeping, and could produce a stable high-burden
state without making the growth law do the work. It is the first candidate
mechanism in this project to arrive from an observation rather than from
inspecting the model.

**(!) THE FRAMING TEST PASSES: bistability was never required (D031).** Before
adding a fourth mechanism, the requirement itself was tested. Lindner's old-pole
cell is not sitting in a second basin. It inherits a physical inclusion body at
every division, which is a **continuously renewed perturbation, not an
attractor**, and a monostable system under a renewed perturbation has a
stationary offset with no separatrix anywhere.

A MONOSTABLE two-state diluted model under the calibrated hyperbolic law, with no
sequestration, no second attractor and no extra state variable, reproduces the
measured lineage difference. The only addition is asymmetric partitioning of
aggregate at division. Across 728 cells: **43 in band with the mechanism on, 0 in
the control**, `f` from 0.60 to 0.99. The control is exact rather than close — all
66 `f = 0.5` cells give an aging effect of **0.0 with standard deviation 0.0**,
because symmetric partitioning in half the volume leaves concentration unchanged,
which is what the dilution term already encodes. The effect is monotone in `f` in
every setting.

**(!) The quantitative match is an ACCOMMODATION, not a post-diction (D032).**
A monotone family rising from exactly zero crosses any band below its maximum, so
"43 cells in band" shows the curve is tall enough, not that the model predicts
1.0-1.8%. There were FOUR free knobs, not two: `f`, `mu0`, `p_qc`, `j/j_crit`.

Removing them: `p_qc` and `mu0` **cancel analytically** — `(B_old - B_new)/k_mu`
is identically `32 x (B_old - B_new)` as a proteome fraction, verified exact in
code — and `f` is pinned by the full text at **1, not 0.6-0.99**, because the
inclusion body is a single indivisible object ("46.5% of the cells contain only
one", new-pole progeny "devoid of parental inclusion bodies"). The continuous-`f`
model is a mean-field stand-in for an all-or-nothing partition, and that is a
stated mismatch.

What remains is one falsifiable number: the model requires the old-pole aggregate
to be **0.037% to 0.063% of the proteome**. The nearest measurement (Tomoyasu et
al. 2001 Mol Microbiol, doi:10.1046/j.1365-2958.2001.02383.x, via PubMed) reports
5-10% in `rpoH`-null cells at 30 C and wild-type aggregation as **undetected** —
a bound, not a value, sitting 79x to 271x above the requirement. So the test is
**well-posed and currently unmeasurable**, and is reported as an experimental
target rather than a claimed success. Measuring the wild-type aggregate fraction
to about 0.01% decides it.

**The structural conclusion does not depend on any of that**, because it rests on
a control that is exactly zero for an algebraic reason rather than on a fit.

**D026's surviving claim is WITHDRAWN.** "Rejuvenation is only a coherent
category in a bistable system" assumed the old-pole cell occupies a second basin.
It does not. The three prior failures were scoring a quantity the observation
never required: they searched for an attractor whose burden would have to be 7.5x
to 254x above what the measurement permits, and the measurement never called for
one.

This does not rescue the model class on D029's point; it relocates it. What the
model cannot do is place a stable ATTRACTOR at a 1% growth cost. It was never
asked to.

**Audit flag settled.** Twelve bistable sequestration cells gave four distinct
loss values. The draws are independent (12 distinct parameter tuples for 12
cells). The identical `1.000` is the linear-arrest law saturating —
`mu = mu0.max(0, 1 - (u+a)/k_mu)` returns exactly zero and every one of those
states sits 3.4x to 43x past arrest. So `1.000` is a **clamp, not a measurement**,
and D029's "27x to 92x too severe" understates the miss rather than overstating
it.

**Running count of collisions with existing observation: three.** Regulation
(D025), sequestration-as-reservoir (D026 named it, D029 tested it), and the
magnitude of the aggregate-laden state (D029). None confirmed the theory; each
named something the model lacks. That is the correct count and it is not to be
inflated by upgrading a qualitative match to a confirmation.

**(!) The Lindner number was misread, and the correction is a factor of thirty
(D028).** ">30% **of** the loss of reproductive ability" is a SHARE of the aging
effect, not the effect. The full text (PMC2268587, via PubMed) gives
`[Delta(GR_old - GR_new)]mean/GR_mean = -3.95 +/- 0.5%` and the
aggregate-attributable share as `Agg/(Agg + Pole) ~ 30-40%`, so the measured
quantity is **1.2-1.8%**. The misreading survived two full analysis cycles and
was written into a protocol as a worked example before the full text stopped it.
An abstract is not a source for a number.

**The sequestration post-diction FAILED, and it is the sharpest failure yet
(D029).** D024 was verified on the extended three-state system first (median
1.5e-12, max 4.7e-11; `k_seq = 0` reduces to the two-state field at 0.00e+00), so
nothing downstream is void. Across 384 qualified cells in each of two regimes:
**0 in band**. Every bistable cell predicts a reproductive loss of 0.482, 0.508,
0.954 or 1.000 against a measured [0.0104, 0.0178] — **27x to 92x too severe**,
uniformly, with no marginal case.

Sequestration is not inert: with sequestered aggregate exempt from the growth
cost, bistable cells with the mechanism ON outnumber the control 12 to 3. It was
the right kind of idea. The high state it creates is simply always complete
arrest.

**Inverting the growth law gives the size of the miss.** The burden that would
produce a measured-size loss is 7.5x to 254x below where the model puts its
second attractor. A real old-pole E. coli carries a visible inclusion body and
still divides at 96-99% of normal; every high state this model can reach is a
cell that has essentially stopped. The next candidate — from the observation, not
from the equations — is a **size-limited deposit**, saturating sequestration into
a finite number of foci rather than an unbounded sink, so that the high state's
burden is set by deposit capacity rather than by where the removal curve bends.
Not yet run.

**One process defect, recorded rather than patched.** Under the original band the
literal criterion would have PASSED, carried entirely by `k_seq = 0` control
cells. D028 fixed the band, the growth law and the falsifier but never required
the mechanism under test to be switched on. That is now protocol rule 6, and
`verdict()` reports `mechanism_passes` beside the literal `passes` instead of
replacing it.

**Antecedent check A1 — the machinery damages itself, and the theorem survives
(D027).** The derivation assumes influx and clearance capacity are independent.
In a cell they are not: chaperones and proteases are themselves translated at the
per-codon error rate that produces the damage. Making capacity error-dependent as
`C_enz = C_0/(1 + eps.load)` over four decades of `eps` — in both a parametric
mode (`load = j`) and a state-dependent mode (`load = u + a`, the one that puts
capacity inside both gradients) — leaves `det J = -(grad R x grad G)` at machine
precision: floor 2.2e-14, worst median 4.6e-13 where capacity is down to **1.8 %**
of nominal.

**There is no corrected form, and that is the result.** The row operation needs
only that `j` be additive in `du/dt` and absent from `da/dt`, and that
`du/dt + da/dt = j - R` be exact. How the parameters depend on `j`, or on the
state, is irrelevant to both. The antecedent must therefore be stated as what it
requires — state-independent total influx, and mass balance counting all outflow
— not as independence of influx and capacity, which the theorem never needed.

The check is also weaker than it looks, and that is recorded rather than glossed:
since `du/dt = j - R - G` holds pointwise and differencing is linear, the identity
carries **no truncation term**, so the residual can only be roundoff. Measured
slope in the step size is **-0.97**, with no V-shaped minimum over four decades —
that prediction confirmed. Run at the repo's habitual `h = 1e-6` the ladder shows
a spurious rise; the analytic argument, not the ladder, carries the theorem.

**What the coupling does destroy is the shortcut.** `{G = 0}` stops being a fixed
curve, so `j_crit = R(u*,a*)` becomes a self-consistency condition and fold-finding
grows from two equations in `(u,a)` to three in `(u,a,j)`.

**Direction, reported independently: the boundary moves a long way, the exponent
does not move at all.** Median `j_crit` falls to **0.32x** (influx) and **0.13x**
(burden) of the frozen value at the top of the ladder. The critical-slowing
exponent is unchanged — paired median -0.4763 damaged against -0.4813 frozen,
Wilcoxon p = 0.312, n = 19. So self-damage makes collapse happen **sooner, not
steeper**: it is still a generic saddle-node. **No new prediction.**

One exact new necessary condition falls out: in influx mode `j <= C_0` becomes
`j <= (sqrt(1 + 4.eps.C_0) - 1)/(2.eps)`, tending to `sqrt(C_0/eps)`. A linear
capacity ceiling becomes a **square-root** one — doubling the machinery buys only
`sqrt(2)` in tolerable error rate once the machinery is itself error-prone. It is
never violated and never binding (`j_crit/j_max` median 0.039–0.186, max 0.623),
so it is recorded as a bound, not as the boundary. Fold recovery at large `eps` is
incomplete and **non-monotone** in `eps`, which identifies it as continuation
failure; folds are not reported as disappearing.

A working manuscript for the whole phase 3 result is
`manuscript/COLLAPSE_BOUNDARY.md`.

Reproduce with `python scripts/phase3/fold_theorem.py`,
`python scripts/phase3/dilution.py` and
`python scripts/phase3/boundary_structure.py` and
`python scripts/phase3/calibration.py` and `python scripts/phase3/pareto.py`;
asserted by `tests/phase3/` (63 checks, of which 56 are model-level and run on a
clean checkout).

## Phase 1 experiment D closed; Phase 2 synthesis final

On 2026-08-01 experiment D was closed scientifically and integrated into the
Phase 2 synthesis.  The two documents are `EXPERIMENT_D_FINAL.md` and
`PHASE2_CLOSURE_FINAL.md`; the latter supersedes the working note
`PHASE2_CLOSURE_PENDING_D.md` under the gitignored closure directory and
carries every one of its negative results forward unchanged.

The run analysed is the checkpointed recovery run
`results/phase1/D_checkpointed_20260731T223225-0500/`.  **60 backgrounds
requested, 58 completed, 46 usable, 12 model-unusable (`no viable state at
j_lo`), 2 unresolved timeouts (backgrounds 19 and 37), zero numerical errors,
zero process failures.**  An unresolved background exceeded the 3600 s wall
limit and was not evaluated within the budget.  It is **not a failure**, is not
counted in `n_errors`, contributes no rows, and must never be reported as one.

Validation passes **36/36**: live source, config and Latin-hypercube hashes match
what the run recorded; all 58 checkpoints pass the runner's own identity gate
with zero rejections; re-merging them reproduces `interactions.tsv` and
`backgrounds.tsv` **byte for byte**; all three nulls and all four excess columns
recompute exactly; and `_pairSummary` recomputed from the shipped TSV reproduces
every field of `summary.json`.  `46 usable x 3 pairs x 36 = 4968` cells, of which
`46 x 3 x 25 = 3450` are genuinely double.

**The raw summary's headline fractions were cell-level descriptives and are not
inference.**  Their 3450 cells are 46 parameter draws x 3 pairs x 25 correlated
cells.  Every estimate in the closure is formed with the background as the unit,
with 95 % CIs from 10,000 replicates resampling whole backgrounds (seed
20260801); the only p-value anywhere is an exact binomial sign test over
backgrounds.

Supported at the background level, on all three nulls: `influx x total_capacity`
and `nascent x total_capacity`.  **Not** supported on the primary set under the
additive or Bliss null: `influx x chaperone_only` — and its multiplicative pass
is not corroboration, because the multiplicative null is the weaker test wherever
a single perturbation is protective, which it is in 32.6 % of that pair's cells.
The most robust result is categorical and null-free: **synthetic collapse in
43 of 46 backgrounds** (95 % CI [0.848, 1.000]), at least 71.7 % of all sixty
requested draws under the most conservative denominator.

Two audit findings travel with those numbers.  `multiplicative_median_excess` in
`summary.json` is a **log-scale** median under a name that states no scale — a
labelling defect, not a sign error.  And the Bliss result is **not** internally
inconsistent: worse means *negative* excess on the survival scale, so a negative
median and a majority above 0.5 are the same statement.

One new negative result: **chaperone-only knockdown is not universally a
burden.**  In 15 of 46 usable backgrounds it *lowers* total burden in all 25
double cells, which follows from the model's own form (chaperone binding
sequesters substrate away from the protease).  Rescue capacity is therefore not a
scalar, and the law is restated accordingly.

**No empirical claim was made.**  No organism data entered Phase 1 or Phase 2.
The six falsifiable predictions the computational work now justifies are listed
in `PHASE2_CLOSURE_FINAL.md` §7 and are labelled empirical hypotheses, untested.

Analysis code is tracked (`scripts/phase2/d_final.py`), the run identity is
pinned (`scripts/phase2/D_RUN_HASHES.json`), and both are asserted by
`tests/phase2/test_d_final.py`.  Detailed outputs are gitignored under
`results/phase2/closure_20260731T220024-0500/D_final/`, so
`scripts/phase2/check_d_closure.py` is the tracked bridge across that boundary:
**177 checks** over the pinned hashes, the counts, the nine majority counts and
verdicts, the collapse rates and every sensitivity bound, plus the requirement
that both closure documents still state them.  Without the run root it prints an
explicit `SKIP` and exits 0 rather than passing silently.  Full suite:
**249 passed** under Python 3.12.11.

## Phase 2A matched equivalence benchmark submitted

On 2026-07-31 the Phase 2A matched benchmark was gated and launched on both
hosts.  The targeted Phase 2 suite passes (`80 passed`, exit 0) and the full
canonical suite passes (`149 passed, 23 warnings`, exit 0) under Python 3.12.11.
The warnings are the pre-existing multiprocessing/fork deprecation warning.

T0, the `epsilon -> 0` right-hand-side and Jacobian identity test, **passed** all
10 checks with exit 0: RHS relative discrepancy `4.800103e-06` and Jacobian
`8.263311e-06` at `epsilon = 1e-6`, log-log slopes `1.0023` and `1.0028`.  The
observed exponent of 1 is what section 3 of
`theory/MATCHED_IMPLEMENTATION_PROTOCOL.md` requires.  T0 passing licenses the
downstream label comparison; it is not itself a scientific result.

One code defect was found and fixed, and it was a defect *for the matched
benchmark specifically*: each cell recorded only `tsv_sha256`, a byte hash of a
file whose every row carries a wall-clock `seconds` column.  That hash can never
match between boron's free arm and nitrogen's free arm even when every computed
value is identical, which is precisely the cross-host comparison the benchmark
exists to make.  A `payload_sha256` over the result columns was added alongside
it, and four tests now pin the property, including one asserting that `seconds`
is the *only* non-deterministic column.  No scientific criterion, threshold, or
mapping was changed; the T0 tolerances and the epsilon ladder are untouched.

A smoke benchmark passed before any full job was launched: correct schema,
`payload_sha256` identical across a repeat run and across `workers = 1` versus
`workers = 4`, zero numerical failures, and exact label and admissibility
agreement between the boron and free arms at `epsilon = 1e-6` under both
protocols.

Two jobs are **running**.  Boron runs the full 28-cell factorial as a
`systemd-run --user` transient unit; nitrogen runs the 14-cell free-limit
counterpart under Slurm at 1 CPU / 4 GB.  Both use n = 2000, seed 20260731, and
one numerical-library thread per worker.  Both hosts were proven to regenerate a
bit-identical sample matrix under different numpy builds.  Exact commands, unit
and job IDs, PIDs, source hashes, and verification evidence are in
`PHASE2_JOBS.md`.

**No Phase 2 result has been produced, inspected, or validated.**  The jobs were
still running when this was written and their partial output has not been
interpreted.

Two constraints govern any later reading of that output, and both are stated in
full in `theory/MATCHED_IMPLEMENTATION_PROTOCOL.md` section 6:

- **nitrogen's model is the `epsilon = 0` limit** of the boron family, not a
  different model.  Setting every substrate-binding affinity to zero turns
  boron's equations into nitrogen's term for term; the sequestration ladder
  makes that limit a single scannable scalar.
- **percentages are not comparable outside the matched benchmark.**  They require
  the same sample matrix, the same root protocol, the same named admissibility
  criterion, a stated epsilon rung, and every matched cell complete.  The older
  50,000-sample independent nitrogen sweep is **not** the matched counterpart and
  must not be reported as the matched result.

## Phase 1 canonical tests repaired; experiments running

On 2026-07-31 the two invalid numerical tests were repaired without changing
the theory, experiment configurations, or manuscript.  The full canonical
suite passes: `69 passed, 19 warnings in 66.36s` under Python 3.12.11.  The
warnings are the existing Python multiprocessing/fork deprecation warning.

Phase 1 experiments A-D are running.  One instance of each was launched on
2026-07-31 at 16:29:46 -0500 on `boron` with the unchanged Phase 1 configs,
as `systemd-run --user` transient units, writing to the clean run root
`results/phase1/run_20260731T162946-0500/`.  Two earlier duplicate launch
attempts left no live host processes; their empty output directories were
moved to `results/phase1/quarantine_preclean_20260731T162755-0500/` rather
than deleted.  Exact commands, unit names, PIDs, resources, config hashes,
logs, and output locations are in `JOBS.md`.

**No result has been produced, inspected, or validated.**  The jobs are
running; nothing in this file asserts that any Phase 1 prediction has been
tested, confirmed, or falsified.  Read `JOBS.md` for the operational record
and wait for the runs to complete before drawing any scientific conclusion.

The repaired scripts, tests, and configs are committed on branch `master`.
The repository has no remote configured, so nothing was pushed.

## Phase 0: theory reconstruction completed

The canonical law, Pareto geometry, site-resolved damage influx, conserved-resource dynamics, scope, falsifiable predictions, empirical program, and a substantive first-pass manuscript are now defined consistently. Discredited components are quarantined in the rejection ledger.

## Next scientific tasks, in priority order

**Superseded by `PHASE2_CLOSURE_FINAL.md` §8.**  The list below is the Phase 0
ordering and is kept as a record of what was planned before Phase 1 and Phase 2
ran; items 1 and 3 are done, and the current open list is §8 of the closure.


1. Prove or numerically delimit existence and local stability regions for the two-state conserved-resource model.
2. Specify operational damage thresholds and viability readouts for an initial organism and condition.
3. Design a perturbation matrix that independently varies burden composition, translation strategy, chaperone allocation, and degradation capacity.
4. Build uncertainty-aware estimators for site-resolved substitution probabilities and damage weights.
5. Test proximity-to-boundary scaling and burden-capacity interactions in preregistered held-out conditions.
6. Compare multi-proxy predictions against scalar mistranslation, growth rate, and expression-only baselines.
7. Curate and insert primary literature citations.
8. Only after validation, consider comparative or evolutionary extensions.
