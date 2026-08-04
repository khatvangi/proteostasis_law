# An Exact Collapse Threshold for Conserved-Resource Models of Protein Quality Control

**Kiran Boggavarapu**

Department of Chemistry and Physics, McNeese State University, Lake Charles, LA 70609, USA
kiran@mcneese.edu

**MSC** 92C40 (primary); 37G10, 34C23, 92C42 (secondary)

**Keywords** proteostasis; saddle-node bifurcation; conserved resource; collapse threshold; growth dilution; clearance network

---

## Abstract

A cell clearing misfolded protein with a finite pool of chaperones and proteases tolerates damage influx only up to a threshold, beyond which burden runs away. Where that threshold sits is normally found by numerical continuation, one parameter set at a time. We show it need not be. For any clearance network in which the damage influx is state-independent and mass balance accounts for all outflow, the Jacobian determinant factors as `det J = −(∇R × ∇G)`, where `R` is total removal flux and `G` the aggregate nullcline field. A saddle-node therefore occurs exactly at a constrained critical point of removal on the aggregate nullcline, and the critical influx is `j_crit = R(u*, a*)` with no sweep required. The identity generalises to arbitrary state dimension, survives growth dilution, and survives the case in which the clearance machinery is itself damaged by the influx it clears. Three consequences follow. Collapse occurs while the machinery runs at roughly 6–18% of maximum velocity, so the naive capacity ceiling overestimates the tolerable influx by about thirteenfold. Under cell division the ceiling fails outright, and the boundary of a dividing cell exists only because burden slows growth; the critical influx then decomposes exactly as `j_crit = C_enz·φ_enz/(1−δ)`, with the enzymatic factor nearly invariant to division rate. When the machinery damages itself, a linear capacity ceiling becomes a square-root one. Against four published observations the underlying two-state model fails in a single direction: it cannot place a stable attractor at low burden. One of those failures resolves not by adding mechanism but by removing a requirement the data never imposed, which yields a falsifiable prediction: the aggregate load of an aging *Escherichia coli* old-pole cell must lie between 0.047% and roughly 0.5% of total protein, depending on how much of that aggregate sits in the visible polar focus, and in every case below the current detection bound.

---

## 1. Introduction

Protein quality control is a resource-limited clearance problem. Chaperones refold damaged monomers, proteases degrade them, and both draw on pools that are finite and shared. Damage arrives continuously from translation error, thermal stress, and oxidative modification. When influx exceeds what the machinery can dispose of, aggregate accumulates and the cell loses viability. That much is not controversial and is the shared premise of essentially every quantitative model of proteostasis.

What is not settled is where the threshold sits, and models built on this premise answer it numerically. One fixes a parameter set, integrates, continues the stable branch in the influx, and records where it terminates. The answer holds for those parameters. Changing the network requires the whole computation again, and the result carries no statement about why the boundary is where it is.

Two intuitions fill that gap and both are wrong. The first is that collapse occurs when clearance saturates, so the boundary should sit near the maximum removal flux. The second is that adding cell division only helps, since dilution is an additional disposal channel. We show below that collapse occurs far below saturation, and that division removes the bound that made the capacity argument work at all.

The result reported here is that the threshold is algebraic. Two structural facts about clearance models — the influx enters one state equation, and mass balance is exact — force the Jacobian determinant to factor into the cross product of two gradients. A saddle-node is then a constrained critical point of total removal, in the Lagrange sense, on the nullcline of the non-influx equation. The critical influx is the value of removal at that point.

The hypothesis this requires is weaker than it first appears, and identifying its sharp form is part of the contribution. The theorem does not need influx and clearance capacity to be independent. It needs only that total influx be state-independent and that mass balance count all outflow. Capacity may depend on the error rate, on the burden, or on both, and the identity is unaffected — which matters because in a cell the chaperones and proteases are themselves translated by the error-prone apparatus whose product they clear.

Section 2 fixes the model class and states the hypotheses. Section 3 states and proves the theorem. Section 4 gives four corollaries: the *n*-state generalisation, the decomposition under growth dilution, the behaviour under self-damaging capacity, and the square-root capacity ceiling that follows from it. Section 5 reports numerical verification. Sections 6 and 7 use the theorem to characterise the boundary and to cut a strategy space. Section 8 reports what happened when the model was put against published measurements, which is mostly failure, and Section 9 states what the theory predicts and what would falsify it.

## 2. The model class and its hypotheses

### 2.1 State equations

Let `u` denote damaged or misfolded soluble monomer and `a` aggregate burden, in compatible concentration units. The baseline field is

```
du/dt = j − v_ref(u, C_f) − v_degU(u, D_f) − n(u) − g(u,a) + v_dis(a, C_f)
da/dt =                                       n(u) + g(u,a) − v_dis(a, C_f) − v_degA(a, D_f)
```

with rate laws

```
v_ref  = k_ref C_f u/(K_ref + u)      refolding
v_degU = k_U D_f u/(K_U + u)          soluble degradation
n      = k_n u^m,  m > 1              nucleation
g      = k_g u a                      aggregate growth
v_dis  = k_dis C_f a/(K_dis + a)      disaggregation
v_degA = k_A D_f a/(K_A + a)          aggregate clearance
```

The nucleation term is superlinear and the clearance terms saturate, and this asymmetry does all the work in Section 6. Both forms are standard. Autocatalytic aggregate formation with an effective order greater than one is the shared structure of nucleation-dependent polymerisation models (Ferrone 1999; Knowles et al. 2009; Cohen et al. 2011), and secondary nucleation, in which existing aggregate surface catalyses the formation of new nuclei, supplies the mechanism (Cohen et al. 2013; Meisl et al. 2016). Michaelis-Menten clearance follows from enzyme-limited refolding and proteolysis. We take `m > 1` as a modelling assumption rather than a fitted exponent; the theorem does not depend on its value.

Chaperone and protease pools are conserved and shared. Under rapid-equilibrium binding with free substrate concentrations,

```
C_f = C_tot /(1 + N_f/K_N + u_f/K_CU + a_f/K_CA + O_f/K_CO)
```

and similarly for `D_f`. Nascent chains `N` consume chaperone capacity without contributing damage influx, which is the mechanism by which ordinary translational load can shift the boundary.

Write the total removal flux and the aggregate field as

```
R(u,a) = v_ref + v_degU + v_degA        G(u,a) = da/dt
```

`R` is what leaves the damaged compartment; `G` is the non-influx equation. Both are used throughout with these meanings.

### 2.2 Hypotheses

**H1 (state-independent influx).** The damage influx `j` is a parameter. It appears additively in `du/dt` and does not appear in `da/dt`.

**H2 (exact mass balance).** The internal transfer between `u` and `a` cancels between the two equations, so `du/dt + da/dt = j − R` holds identically.

**H3 (smoothness).** `R` and `G` are continuously differentiable on the region of interest.

**H4 (deterministic, well-mixed, finite-dimensional).** The system is an ODE in concentrations.

*Remark 1.* H1 constrains the influx, not the capacity. Nothing forbids `C_tot` or `D_tot` from depending on `j` or on the state. The gradients in Section 3 are taken at fixed `j`, so parameter dependence on the influx cannot enter them, and the row operation needs only that `j` be additive in one equation and absent from the other. This distinction is not cosmetic: chaperones and proteases are translated by the same error-prone machinery that generates the damage, so a cell violates the stronger reading of H1 and satisfies the weaker one. Corollary 3 makes this quantitative.

*Remark 2.* H4 is false in cells in two known respects. Aggregates in *E. coli* are spatially segregated into polar inclusion bodies rather than well mixed, and molecule numbers per damaged species are small. Both push in the same direction: noise-induced escape precedes the deterministic fold, so a measured boundary should sit below `j_crit` and be smeared across a population rather than sharp. We return to this in Section 9.

## 3. The theorem

**Theorem 1.** *Under H1–H3, for the system of Section 2.1:*

1. *equilibria are exactly the points of the curve `{G = 0}` at which `R = j`, and that curve does not depend on `j`;*
2. `det J = −(∇R × ∇G)` *identically;*
3. *a saddle-node bifurcation therefore occurs exactly at a constrained critical point of `R` restricted to `{G = 0}`, and the critical influx is*

```
j_crit = R(u*, a*).
```

*Proof.* Two structural facts do the work.

By H1, `j` appears in `du/dt` and nowhere in `da/dt`. The aggregate nullcline `{G = 0}` is therefore a fixed curve in burden space: changing the load slides the state along it without reshaping it. An equilibrium requires `G = 0` and `du/dt = 0`, and by H2 the second condition is `R = j` on that curve. This gives (1).

By H2, `du/dt + da/dt = j − R` exactly. Apply the determinant-preserving row operation `row₁ → row₁ + row₂` to the Jacobian:

```
det J = det [ ∂(du/dt)/∂u   ∂(du/dt)/∂a ]  =  det [ −R_u  −R_a ]
            [ ∂(da/dt)/∂u   ∂(da/dt)/∂a ]        [  G_u   G_a ]

      = −(R_u G_a − R_a G_u) = −(∇R × ∇G).
```

That is (2). A saddle-node requires `det J = 0`, hence `∇R ∥ ∇G`, which is the Lagrange condition for a critical point of `R` subject to `G = 0`. Combining with (1), the critical value is `j_crit = R(u*, a*)`. ∎

Informally: **the collapse boundary is where total removal stops responding to burden along the aggregate nullcline.**

*Remark 3.* The proven statement is *critical point*, not *maximum*. An earlier version of this result asserted a constrained maximum; that gloss is withdrawn. Scanning `R` along the branch reached by taking the first root in `a` shows it rising monotonically to the end of that branch in every case tested, and the solved fold state has `dG/da` between 9.1×10⁻³ and 3.9×10⁻², so the critical point lies past the nullcline's turning point rather than on the first-root branch. Whether it is a maximum over the whole closed curve is not established. No quantity in this paper depends on the distinction, because all of them solve `det J = 0` directly.

**Corollary 1 (fold location is a root solve).** Neither `G = 0` nor `det J = 0` contains `j`. Locating the fold is therefore a two-dimensional root solve rather than a continuation sweep in the influx, which reduces the cost of surveying a parameter box by the length of the sweep.

At the base parameter point the Lagrange condition is visible directly rather than inferred: the removal contours meet the aggregate nullcline tangentially at the solved fold, so `∇R` and `∇G` align there to within a sine of `3.5×10⁻¹⁰`, and the critical influx is simply the removal flux evaluated at that state (Fig. 1a).

![Figure 1](../figures/fig1.pdf)

**Fig. 1** The collapse threshold as a constrained critical point. **(a)** Phase
plane at the base parameter point. The aggregate nullcline `{G = 0}` (dark) is a
fixed curve; grey contours are total removal `R`. At the solved fold
`(u*, a*) = (0.4166, 0.2650)` the gradients `∇R` and `∇G` point the same way,
which is the Lagrange condition of Theorem 1; they are drawn from offset origins
because, being parallel, a common origin hides one behind the other. The sine of
the angle between them is `3.5×10⁻¹⁰`. The critical influx is the removal flux
evaluated there, `j_crit = R(u*, a*) = 0.1542`. **(b)** Both equilibrium
coordinates against influx, traced along the same curve. Cramer's rule on the two
columns gives `du*/dj = −G_a/det J` and `da*/dj = G_u/det J`, so the branch
carries two distinct, ordered singularities. At `j_turn = 0.15409` the soluble
pool has a HORIZONTAL tangent, where `G_a = 0`; this is a regular point of the
branch, since `det J = R_a·G_u = 2.027×10⁻³` there rather than zero, and it lies on
the stable branch at `j_turn/j_crit = 0.9990`. At `j_crit` both coordinates have a
VERTICAL tangent, which is `det J = 0`. The inset is a ×182 zoom on the same
computed data, not a schematic. The aggregate coordinate has no horizontal
tangent anywhere: `G_u > 0` at all 305 traced points (minimum `2.08×10⁻³`), and in 40
draws from the parameter box no point with `G_u < 0` was found, for the reason
given in Section 3.1.

### 3.1 Relation to infinitesimal homeostasis

Theorem 1 sits close to an existing body of work, and the relation needs stating precisely because the hypotheses nearly coincide.

Wang et al. (2021) study input-output networks in which a parameter enters exactly one node's equation, which is our H1, and ask when the input-output function `x_o(I)` has vanishing derivative. That question has a long history in the adaptation literature, where the stronger condition of perfect adaptation requires the derivative to vanish identically over an interval, and where the network topologies capable of achieving it have been characterised numerically and then structurally (Ma et al. 2009; Ferrell 2016; Araujo and Liotta 2018; Khammash 2021). Cramer's rule applied to the implicitly differentiated equilibrium condition gives `x'_o = ±f_{ι,I} det(H)/det(J)`, where the homeostasis matrix `H` is obtained from the Jacobian by deleting its first row and last column. Infinitesimal homeostasis is therefore `det(H) = 0`, and the programme built on that identity factors the determinant combinatorially, associating each irreducible factor with a subnetwork motif of structural or appendage type (Golubitsky and Stewart 2017; Reed et al. 2017; Golubitsky and Wang 2020; Huang and Golubitsky 2022; Madeira and Antoneli 2024; Antoneli et al. 2025; Lin et al. 2026).

The two conditions are complementary degeneracies of the same expression: theirs is the vanishing of the numerator, ours of the denominator. For the two-state system the account can be made completely explicit, and it separates into three cases, because Cramer's rule applied to the two columns gives both

`du*/dj = −G_a / det J`  and  `da*/dj = G_u / det J`.

**Output the aggregate.** Taking `ι = u` and `o = a` is literally the Wang et al. setup: two nodes, one arrow, no regulatory node, so `H` is 1×1 and `det H = G_u`. In their classification this is the degree-one structural type, Haldane homeostasis, which occurs only when the single coupling along the input-output path vanishes. It does not vanish here, and the reason is a sign structure rather than a parameter value. Nucleation and aggregate growth both rise with `u`. Raising `u` also titrates chaperone and protease away from disaggregation and aggregate clearance, and those two fluxes enter `G` with a minus sign, so their contributions to `G_u` are positive as well. All four terms have the same sign. Numerically, `G_u > 0` at all 305 traced points of the nullcline at the base parameter point, minimum 2.08×10⁻³, and across 40 draws from the parameter box no point with `G_u < 0` was found, the smallest value encountered anywhere being 6.7×10⁻¹². **The only homeostasis type available to this network is Haldane, and the sign structure of `G` excludes it.**

**Output the soluble monomer.** Taking `o = u` places the output at the input node, which their framework excludes by requiring the two to be distinct, so what follows is an analogue rather than an instance and we do not call it infinitesimal homeostasis. The derivative condition is nonetheless the natural one, and `du*/dj` vanishes where `G_a = 0`. That point exists, is unique on the branch, and is generic rather than degenerate: at it `det J = R_a·G_u`, which is `2.027×10⁻³` at the base parameter point rather than zero. It is the turning point of the nullcline itself, where `{G = 0}` runs vertical in the burden plane. Physically the soluble pool stops responding to influx there because every further increment is routed into aggregate.

The branch therefore carries **two distinct singularities in a fixed order**, and Fig. 1b traces both on the same computed curve: the horizontal tangent of `u*` at `j_turn = 0.154090`, then the vertical tangent of both coordinates at `j_crit = 0.154239`. The first is a regular point of the branch — `det J = R_a·G_u = 2.02749×10⁻³` there, not zero — so it is generic rather than a degeneracy that a perturbation would remove. The separation is `j_turn/j_crit = 0.9990` at the base parameters. We do not offer it as a testable prediction: across 30 draws the locator placed it on the stable branch in only six cases, with median ratio 0.9963, so it is a feature of the branch geometry rather than something an experiment could separate from the boundary itself.

One numerical consequence is worth stating because it is structural rather than incidental. Tracing `{G = 0}` by solving for `a` at each fixed `u` fails exactly at this locus: the two roots in `a` merge where the nullcline runs vertical, so the traced curve returns as two disconnected pieces with the fold sitting in the gap between them. The root-finder's failure and the horizontal tangent of `u*(j)` are the same point, not two problems, and contour tracing follows the curve through it.

**Neither numerator.** `det J = 0` is the fold, excluded from their setting by the hyperbolicity assumption their derivation requires. Geometrically, infinitesimal homeostasis is a horizontal tangent of the equilibrium branch and a saddle-node is a vertical tangent of both coordinates at once; neither condition implies the other, and on this branch the two occur at different, ordered values of the influx.

What makes our determinant factor is H2, which has no counterpart in the input-output setting. Mass balance supplies the row operation that replaces the influx row by `−∇R`, so `det J` factors into gradients of identified physical fluxes rather than into combinatorial blocks. That is also why `j_crit = R(u*, a*)` is an evaluation rather than a classification: the theorem locates a threshold, where the homeostasis programme enumerates the mechanisms by which a different degeneracy arises.

Whether the combinatorial machinery transfers is open and worth having. A Frobenius-König factorisation (Brualdi and Ryser 1991) of `det J = −det[∇R; ∇G; ∇C]` would classify collapse mechanisms by network topology in the way structural and appendage blocks classify homeostasis mechanisms. Section 6 decomposes the shortfall into saturation and sequestration contributions, which is a two-term physically motivated version of the same question.

The complementary question, when regulation fails rather than when it holds, was raised early in that literature and has been developed less. Nijhout et al. (2014) ask what happens when a homeostatic mechanism is pushed past its range, and Reed et al. (2018) treat homeostasis at equilibria that are not stable. Both predate the determinant machinery and neither is a singularity-theoretic account of the failure. Theorem 1 supplies one for the case in which the failure is a fold.

## 4. Extensions

### 4.1 Arbitrary state dimension

The standing objection to Theorem 1 is that its exactness might be a property of a two-state reduction rather than of clearance networks. It is not.

**Corollary 2.** *Let the state be `x = (u, a, c, …)` with `du/dt` the only equation containing `j`, and let mass balance give `du/dt + da/dt = j − R(x)`. Writing the remaining equations as `G` (aggregate) and `C` (controller and any further states),*

```
det J = −det [ ∇R ; ∇G ; ∇C ]
```

*which vanishes exactly when `∇R` lies in the span of the remaining gradients — the Lagrange condition for a constrained critical point of `R` on the intersection of the non-influx nullclines.*

The proof is the same row operation applied to the first row of an *n* × *n* Jacobian. The practical consequence is that a model in this class can be extended without re-deriving its boundary condition.

We verified Corollary 2 on a three-state system in which the chaperone pool is dynamical under σ32-style control, with synthesis rising as free chaperone falls, which is the actual titration mechanism rather than a phenomenological feedback. Median relative error of the identity is 0.000×10⁰ unregulated and 2.5×10⁻¹¹ to 2.9×10⁻¹¹ regulated, with the `σ₀ → 0` limit reproducing the two-state field exactly. A separate three-state extension splitting aggregate into reactive and sequestered compartments (Section 8.3) gives median 1.5×10⁻¹² and maximum 4.7×10⁻¹¹ over the complete 144-point grid of sequestration settings, growth laws, states and growth-cost conventions.

### 4.2 Growth dilution, and why the capacity bound fails

The model as stated has no cell division, and for most bacterial proteins dilution outpaces proteolysis. Adding it leaves the theorem intact and destroys the naive bound.

Dilution does not affect binding, so the diluted field is `du/dt − μu`, `da/dt − μa`. H1 and H2 both survive, with total removal now

```
R_tot = R + μ(u,a)·(u + a)
```

so `det J = −(∇R_tot × ∇G_dil)` for any dilution law. Verified at median relative error 1.2×10⁻¹⁰ for constant `μ` and 4.7×10⁻¹⁰ for burden-dependent `μ`, with the `μ → 0` reduction exact at 0.0. For constant `μ` the diluted Jacobian is `J − μI`, so the saddle-node condition says exactly that `μ` is an eigenvalue of the undiluted Jacobian.

The removal ceiling does not survive. `R_tot` contains `μ(u+a)`, which is unbounded in burden, so the analytic bound `j_crit ≤ C_enz` is false once cells divide. The consequence is not cosmetic. Under *constant* dilution the fold is destroyed above a threshold: continuing in `μ` at the base parameters gives `j_crit` rising 0.1542 → 0.2456 as `μ` goes 0 → 0.08 while the fold state `a*` diverges 0.265 → 0.750, and by `μ = 0.10` no fold exists. Past that point the low-burden branch no longer terminates and burden rises smoothly with influx. A cell dividing at a burden-independent rate dilutes its way out.

Growth feedback restores the boundary. With `μ(u,a) = μ₀/(1 + (u+a)/k_μ)` a fold exists at every rate tested up to `μ₀ = 0.3`.

> **The collapse boundary of a dividing cell exists because burden slows growth.** Dilution alone is an unbounded disposal channel; only its throttling by the burden it disposes of makes viability finite.

The coupling itself is not in dispute. Expressing protein a cell does not need reduces its growth rate, and the reduction scales with expression level (Dekel and Alon 2005; Shachrai et al. 2010). Coarse proteome partitioning makes the mechanism explicit: raising any sector's share comes at the expense of the ribosomal sector, and the resulting growth law predicts the burden of heterologous expression quantitatively (Scott et al. 2010; Scott et al. 2014; Scott and Hwa 2023). What the theory adds is that this coupling is load-bearing for the existence of the boundary rather than a correction to its position.

This changes what a perturbation experiment must control. Growth rate is part of disposal, so an experiment that lets growth rate float is not holding disposal capacity fixed.

**Corollary 3 (exact decomposition under division).** *Splitting the critical influx into its enzymatic and dilution parts,*

```
j_crit = C_enz · φ_enz /(1 − δ),     φ_enz = R_enz(u*,a*)/C_enz,     δ = R_dil(u*,a*)/j_crit
```

*with both factors dimensionless and in [0,1). The identity is algebraic and holds to 1.6×10⁻¹⁶.*

The two factors behave completely differently under division, and that asymmetry is the content of the corollary. At the base parameters `φ_enz` stays within 0.1245–0.1343 as `μ` runs 0 → 0.08 — a total width of 7.5% of its own mean — while `δ` runs 0 → 0.39 over the same sweep and carries essentially all the variation (Fig. 2). Division therefore does not change the enzymatic condition for collapse; it multiplies the tolerable influx by `1/(1 − δ)`. The escape above is `δ → 1`.

![Figure 2](../figures/fig2.pdf)

**Fig. 2** The enzymatic condition for collapse is insensitive to division; the
dilution share carries the variation. Both factors of Corollary 3 are plotted on
one linear axis against the dilution rate, as a 33-point continuation at the base
parameters under constant dilution — a complete enumeration of the sweep, not a
sample. The shaded band spans the full range of `φ_enz`. Twin axes are avoided
deliberately: rescaling the flat curve to fill the panel would make its variation
look like structure, which would invert the reading. The flatness is not a
flatness of the underlying state — `j_crit` and the fold state `a*` both move
substantially over this same sweep, by the amounts given above.

Across 25 parameter draws having a boundary at `μ = 0`, 23 lose it under constant dilution. The threshold spans 3.3 decades (p10/p50/p90 = 0.0033/0.086/0.328), and `δ` at the threshold has median 0.35: the boundary is typically lost once division does roughly a third of the disposal work.

Calibrating the growth-burden coupling against a measured dosage-response — a 3.2% growth rate reduction when misfolded YFP represents less than 0.1% of total cellular protein (Geiler-Samerotte et al. 2011) — gives slope 32 per unit misfolded proteome fraction and linear arrest at fraction 0.03125, both bounds rather than point values. Under this calibration a boundary survives in 30 of 30 cells with `φ_enz` confined to 0.072–0.147. The destruction of the boundary under constant dilution is therefore an artifact of omitting the measured coupling, not a prediction.

### 4.3 Self-damaging machinery

In a cell, chaperones and proteases are themselves translated by the error-prone apparatus. Raising `j` degrades the machinery that clears the damage, so capacity is a decreasing function of the load. This is the case in which the strong reading of H1 fails.

We modelled it as `C_enz(load) = C_0/(1 + ε·load)` applied to both pools, with `ε` swept over four decades and not tuned to any biological value, in two modes: `load = j` (influx mode) and `load = u + a` (burden mode). The second is the one that can break anything, since it places capacity inside `∇R` and `∇G` themselves.

**The identity survives in both modes at machine precision.** The gradient-normalised residual has floor 2.2×10⁻¹⁴ at `ε = 0`, and worst median 6.4×10⁻¹⁴ (influx) and 4.6×10⁻¹³ (burden) at `ε = 100`, where capacity has fallen to 16.7% and 1.8% of nominal. These are medians over 20 networks drawn from the kinetic box; a median is unbiased under subsampling, unlike a maximum.

There is no corrected algebraic form, and that absence is the finding. The row operation requires only that `j` be additive in `du/dt` and absent from `da/dt`, and that `du/dt + da/dt = j − R` be exact. How the parameters depend on `j` or on the state is irrelevant to either, because the gradients are taken at fixed `j`. Remark 1 is the general statement; this is its verification.

What the coupling destroys is Corollary 1. Under self-damage `{G = 0}` moves with the load, so `j_crit = R(u*,a*)` stops being an evaluation and becomes a self-consistency condition, and fold-finding grows from two equations in `(u,a)` to three in `(u,a,j)`. The theorem survives; the algorithm does not.

Self-damage lowers the boundary substantially without steepening the approach. Median `j_crit` relative to the frozen fold falls 0.999 / 0.990 / 0.925 / 0.640 / 0.322 across the influx ladder and 0.995 / 0.949 / 0.740 / 0.324 / 0.131 across the burden ladder, a threefold to 7.6-fold reduction in tolerable influx. The critical-slowing exponent does not move: paired over networks with both values, median −0.4763 damaged against −0.4813 frozen, Wilcoxon *p* = 0.312, *n* = 19. Collapse under self-damage is still a generic saddle-node. It happens sooner.

**Corollary 4 (square-root capacity ceiling).** *In influx mode every removal flux carries a factor `1/(1 + εj)`, so the necessary condition `j ≤ C_0` of the frozen model becomes `εj² + j − C_0 ≤ 0`, that is*

```
j ≤ ( √(1 + 4εC_0) − 1 )/(2ε)   →   √(C_0/ε)   for large ε.
```

A linear capacity ceiling becomes a square-root one: doubling the machinery buys only √2 in tolerable error rate once the machinery is itself error-prone. The bound is exact and never violated. It is also never binding in the range tested, with `j_crit/j_max` of median 0.039–0.186 and a largest observed value of 0.623 over 8 drawn networks. That figure is a lower bound on the population maximum rather than the maximum itself, since an extremum over a subsample understates by construction; the conclusion does not depend on it, because the bound is not binding at any value below 1. It is recorded as a necessary condition and not as the boundary.

Two limitations attach. Fold recovery is incomplete at large `ε`: continuation with intermediate rungs recovers folds that a direct solve loses, but recovery counts are non-monotone in `ε` (influx 7/5/6/6/4, burden 6/7/7/4/2). A genuine loss of the boundary would be monotone, so this is continuation failure rather than folds disappearing, and where a fold does solve it is a real saddle-node: `sin(∇R, ∇G)` is below 2.0×10⁻⁹ at every one of the folds solved in this sweep, which is a complete enumeration of that sweep rather than a sample from it. The consequence is that Corollary 4 is least numerically checkable in exactly the regime where it would bind.

## 5. Numerical verification

Every quantity below is recomputed by the accompanying code and asserted by its test suite. Two distinct populations are involved and each row names its own, because they differ in what they vary and are therefore not interchangeable. The **load grid** holds kinetics fixed and sweeps the two load coordinates, nascent-chain occupancy against rescue allocation, giving 325 fold states; it answers how the boundary moves with load in one network. The **kinetic box** instead varies the kinetic parameters, 5000 Latin-hypercube draws of which 2884 admit a fold, and answers how the boundary varies across networks. Both are populations of folds, which is why a reader meeting 325 in one row and 2884 in another needs to be told which question each answers. Every entry uses the whole of its population; none is a subsample.

Five of these values were previously reported from a 20-state random subsample of the load grid, under a normalisation discussed below, and all five moved when the full population was used: the median residual from 1.436×10⁻⁷ to 2.00×10⁻⁷ under the old normalisation, the correlation from +0.9987 to +0.9960, the `|G|` maximum from 8.2×10⁻¹⁰ to 1.63×10⁻⁹, the solver maximum from 6.652×10⁻⁷ to 7.56×10⁻⁷, and the tightest bracket in the population from `|λ| = 4.2×10⁻⁹` to `8.10×10⁻¹⁰`. **None of these changes weakens any claim made here.** The identity still sits at the differencing floor, the correlation is still decisive, and `|G|` at 1.63×10⁻⁹ remains four orders below anything that would affect a conclusion. We report the corrections explicitly rather than silently, because a section whose subject is exactness should survive arithmetic done by a reader.

| quantity | population | median | p99 | max |
|---|---|---|---|---|
| `det J` against `−(∇R × ∇G)`, residual | load grid, 325 | 2.34×10⁻¹⁰ | 9.67×10⁻¹⁰ | 1.29×10⁻⁹ |
| \|G\| at recorded fold states | load grid, 325 | 4.95×10⁻¹⁴ | 9.40×10⁻¹⁰ | 1.63×10⁻⁹ |
| direct solver against continuation sweep, relative | load grid, 325 | 3.03×10⁻⁷ | 7.20×10⁻⁷ | 7.56×10⁻⁷ |
| `φ` rebuilt from first principles | kinetic box, 2884 | 1.3×10⁻¹³ | 2.98×10⁻⁹ | 7.25×10⁻⁹ |
| correlation of log sin(angle) with log \|eigenvalue\| | load grid, 325 | +0.9960 | — | — |
| *n*-state identity, regulated three-state | base point | 2.5×10⁻¹¹ | — | — |
| dilution identity, burden-dependent `μ` | base point | 4.7×10⁻¹⁰ | — | — |

Each maximum is reported beside a p99 deliberately. A maximum grows with the size of the population it is drawn from, so the `φ` maximum of 7.25×10⁻⁹ over 2884 draws and the identity maximum of 1.29×10⁻⁹ over 325 are not comparable as stated, and neither is comparable to whatever a reader obtains on a rerun of a different size. The p99 is stable under resampling and carries the same bounding claim. Where a number below is doing bounding work, the p99 is the one to read.

The identity residual sits at the differencing floor. The parallelism residual is not zero at the recorded states and should not be: those states are bracketed approximations whose leading eigenvalue is about −2×10⁻⁴ rather than 0. The +0.9960 correlation between parallelism error and recorded eigenvalue shows the residual is bracket tolerance rather than failure of the identity, and the log–log slope is 1.00, so sin(angle) is proportional to bracket looseness rather than merely increasing with it. The tightest-bracketed state in the population has `|λ| = 8.10×10⁻¹⁰` and `sin(angle) = 7.75×10⁻⁹`.

The residual is normalised by `|∇R||∇G|` rather than by `max(|det J|, |cross|)`. The distinction matters, and not only for presentation: both `det J` and the cross product vanish at a saddle-node, so the second normalisation is `0/0` at an exact fold and returns 1 regardless of how well the identity holds. Measured on the same 325 states, its residual correlates with the recorded eigenvalue at −0.835 in log–log, meaning it degrades as the bracket tightens on the true fold, and its maximum reaches 1.54×10⁻². The gradient normalisation correlates at +0.041 and its maximum is 1.29×10⁻⁹. Stated without a correlation, splitting the population at the median eigenvalue: the max-normalised residual has median 3.13×10⁻⁷ on the tighter-bracketed half against 1.46×10⁻⁷ on the looser, 2.15 times worse where it matters most. Correlations are reported in log–log throughout, matching the +0.9960 statistic above, and the raw Pearson values are −0.221 and +0.073, same signs and same conclusion. An earlier draft reported this pair as −0.262 and +0.060; those values were computed outside the repository and do not reproduce under any definition, and they understated the effect they were cited for. A verification statistic that worsens as it approaches the object being verified is the wrong statistic, and Fig. S1 in the supplementary material reports the contrast.

One caution on the self-damage check of Section 4.3. Because `du/dt = j − R − G` holds pointwise and the central-difference operator is linear, the check reproduces the row operation exactly at any step size and carries no truncation term. Measured slopes in the step size `h` are −0.97 and −0.94 with no V-shaped minimum over four decades, which is pure roundoff falling as `h` grows. The check therefore certifies that the implementation preserves mass balance; the analytic argument carries the theorem. Reporting the ladder at the customary `h = 10⁻⁶` would have shown a spurious rise and invited a false alarm.

## 6. Where the boundary sits

The theorem locates the fold. It does not by itself say how far below the capacity ceiling that is. Writing `R = c_f s_ref + ρ_U d_f s_u + ρ_A d_f s_a` against `ceiling = c_tot + ρ_U d_tot + ρ_A d_tot`, and defining `φ = j_crit/ceiling`, the saturation fractions at the fold answer it.

| at the collapse point | median |
|---|---|
| `φ` observed | 0.0769 |
| refolding saturation `s_ref` | 0.175 |
| soluble degradation saturation `s_u` | 0.155 |
| aggregate clearance saturation `s_a` | 0.056 |
| shortfall recovered by removing saturation | 35.8% |
| shortfall recovered by removing sequestration | 12.6% |

**Collapse occurs while the clearance machinery runs at roughly 6–18% of `V_max`.** It is not capacity exhaustion. Superlinear nucleation overtakes sublinear saturating removal long before removal maxes out, which makes the removal ceiling a correct bound that is about thirteenfold too loose to predict anything on its own.

Two questions must not be merged here. Saturation dominates the *magnitude* of `φ`; the *existence* of a turning point requires the aggregation runaway that drives the free pools down at high burden. The counterfactuals above answer the first only.

The medians above are computed over the whole kinetic box and no draw is excluded from them. This is worth stating because a subset invites exclusion: some draws collapse at `s_a` near 10⁻³, where aggregation is fast enough that the low-burden branch barely exists, and 47 sit below 10⁻⁹. Screening them would be defensible only if they formed a distinct group, and they do not — the distribution runs smoothly across five decades with no gap, and the median of the survivors slides continuously with whatever floor is imposed, from 0.090 at a floor of 10⁻⁴ to 0.355 at 2×10⁻² (Fig. 3, inset). Any floor is therefore a free parameter that moves a load-bearing number by a factor of four. A suspicion that some draws are marginal does not by itself license a threshold, and the complete 2884 reproduce the medians in the table exactly.

![Figure 3](../figures/fig3.pdf)

**Fig. 3** Saturation of the clearance machinery at the collapse boundary, over
all 2884 folds of the kinetic box, with no screen and no exclusion. Shaded strips
are the full distributions, bars the p5–p95 span, open circles the medians
tabulated above, against the dashed line at `s = 1` that capacity exhaustion would
predict. The spread is drawn as prominently as the centre because both halves of
the claim matter: collapse occurs far below saturation, and the p5–p95 widths of
0.881, 0.876 and 0.863 are wide enough that a single measurement discriminates
weakly. The inset is the sensitivity that rules out a screen, plotting the median
`s_a` of the survivors against the floor imposed on them.

The dispersion of these quantities is large enough to matter, and two attempts to reduce it failed. A properly nested design crossing 10 kinetic draws with a 7×7 load grid gives a between-to-within variance ratio for `φ` of 5.9, with per-network spread reaching 13.6× when both load coordinates sweep, overlapping the 8.86× between-network spread. Both figures are largest-observed values over the 10 draws of that design and so understate their populations; the overlap they establish is unaffected, since a larger sample can only widen it. `φ` is therefore network-characteristic but not load-invariant, and an earlier reading of it as a load-invariant material constant is withdrawn. Adding σ32-style regulation, on the hypothesis that a controlled cell sits where its controller puts it rather than where its raw kinetics do, *widens* the p5–p95 width of `s_u` from 0.890 to 0.968 rather than narrowing it, both measured on that experiment's own 30 networks rather than on the kinetic box. The regulated median `s_u` of 0.323 against 0.169 unregulated points toward the capacity-exhaustion picture rather than away from it, though only 14 of 30 regulated networks converged against 24 unregulated, so it should not be quoted as a result without a larger sample.

The defensible claim from this section is that the naive capacity bound is wrong by roughly an order of magnitude and that the boundary is a computable constrained critical point. Not that collapse occurs at any particular fraction of `V_max`.

## 7. The boundary as a constraint on strategy

The theorem supplies a constraint with no free parameters of its own, `j(α,R) < j_crit(R)`, and that makes a trade-off surface computable. Take two strategy coordinates paid out of the same proteome: accuracy investment `α`, which lowers the error rate and hence the influx, and quality-control investment `R`, which raises `j_crit`. Both reduce the ribosomal share and therefore throughput. Cost forms are stipulated, not fitted.

The throughput optimum sits exactly on the feasibility boundary, at `j/j_crit = 1.000000`. This is not numerical coincidence: throughput is strictly decreasing in both coordinates, so any interior point improves by lowering `α` until the constraint binds. A 26×22 grid placed the optimum at `j/j_crit = 0.8975`; tracing the boundary exactly gives 1.0 and 0.89% higher throughput, so the grid figure was discretisation rather than an interior optimum.

The front is not uniformly at the boundary. Along the non-dominated set in (throughput, accuracy), `j/j_crit` runs 0.227 to 0.965, rising monotonically with throughput, so proximity to the boundary is not a property of being on the front but of *where* on it a strategy sits: accuracy-favouring strategies lie far inside the envelope and only the throughput-favouring end presses against it (Fig. 4).

This gives a precise and partly deflationary form to the expectation that evolution concentrates strategies near feasibility boundaries: it holds for throughput maximisers and fails for the rest of the front. A deterministic throughput maximiser has zero safety margin, so any observed margin must come from outside this model — noise, environmental fluctuation, or robustness against parameter uncertainty. The layer cannot adjudicate which, and does not claim to.

![Figure 4](../figures/fig4.pdf)

**Fig. 4** The strategy front, and where along it the feasibility constraint
actually binds. Grey points are all 469 feasible strategies on a 26×22 grid of
(`α`, `R`); coloured points are the 13 non-dominated ones, shaded by `j/j_crit`.
The connecting line is grey rather than colour-mapped because nothing was solved
between the 13 points and colour along a segment would assert values that were
never computed. The star is the throughput optimum, obtained by tracing the
feasibility boundary rather than from the grid, at `j/j_crit = 1.000000`; the best
grid strategy reaches only 0.8975, and that gap is discretisation. The inset gives
the deflationary half the same visual weight as the optimum: `j/j_crit` along the
front spans 0.227 to 0.965, so accuracy-favouring strategies sit far inside the
envelope and only the throughput end presses on it.

## 8. What the model class cannot do

Sections 3 to 5 are mathematics. This section reports what happened when the model was put against published measurements, and the answer is that it failed three times in one direction before a fourth test showed that one of the failures was mis-specified. We report the sequence because the diagnosis is the useful part.

### 8.1 A measurement, misread

The first target was the aging asymmetry in *E. coli*: unstressed cells accumulate aggregates at the old pole and lose reproductive ability relative to new-pole progeny (Lindner et al. 2008). The abstract reports ">30% of the loss of reproductive ability," and the first two analysis cycles took 30% as the measured growth-rate deficit.

That is a misreading, and by a factor of thirty. The full text gives both parts: the lineage difference is `[Δ(GR_old − GR_new)]_mean/GR_mean = −3.95 ± 0.5%`, and the fraction of that deficit attributable to the aggregate rather than to the pole itself is 30–40%. The aggregate-attributable loss is therefore 1.04% to 1.78%, not 30%. All subsequent scoring uses the derived band [0.0104, 0.0178], computed from the two reported quantities rather than quoted from the abstract.

An abstract is not a source for a number.

### 8.2 Three failures in one direction

Three model variants were then scored against that band, each asking whether the model can produce a bistable system whose high-burden attractor carries a 1% growth cost.

Under constant dilution the model is bistable, with a 12.6- to 32.3-fold aggregate ratio between branches. That regime predicts a growth-rate loss of exactly zero, because `μ` cannot respond to burden. The regime producing the bistability is the regime that gets the measured quantity wrong.

Under the physiologically appropriate hyperbolic law the system is monostable in four of six settings tested. Under linear arrest there is no bounded high-burden state at all: `μ` reaches exactly zero beyond `k_μ`, switching dilution off and removing the bound, and the apparent high branch sits at `a ~ 10⁴` and is still growing between sweep steps.

Adding a sequestered compartment — aggregate drained into an inert pool that neither nucleates nor occupies chaperone or protease, which is the natural reading of an inclusion body — is not inert in its effect. It produces 12 bistable cells against the control's 3. But 0 of 384 qualified cells score in band, in either growth-cost regime. Every bistable cell predicts a loss of 0.482, 0.508, 0.954 or 1.000 against a measured 0.0104–0.0178. Inverting the growth law for the burden that would produce a measured-size loss, the model's high state is 7.5 to 254 times more aggregate-laden than the cell that was measured. The reported severity is a lower bound: the 1.000 values are the linear-arrest law clamping at states 3.4 to 43 times past arrest, so the true miss is larger.

A real old-pole *E. coli* carries a visible inclusion body and still divides at 96–99% of normal. Every high state this model can reach is a cell that has essentially stopped.

### 8.3 The requirement was never imposed

The three variants share an assumption that was never checked: that rejuvenation requires bistability. The reasoning was that a daughter inheriting less aggregate has nothing to be rejuvenated *into* unless a second attractor exists.

That is wrong, and the reason is physical. The old-pole cell does not occupy a second basin. It inherits a physical inclusion body at every division, which is a continuously renewed perturbation rather than an attractor. A monostable system under a renewed perturbation has a stationary offset and needs no separatrix anywhere.

Testing it required no new state variable and no new mechanism: the two-state diluted model under the calibrated hyperbolic law, with asymmetric partitioning of aggregate at division in ratio `(f, 1−f)`. Of 728 cells, 709 settled, and 43 scored in band with the mechanism on. The control is exact rather than approximate: all 66 cells at `f = 0.5` give an aging effect of 0.0 with standard deviation 0.0. That falls out of the modelling — symmetric partitioning into half the volume leaves concentration unchanged, which is what the `−μx` term already encodes — so a nonzero control would have signalled an implementation error rather than a small effect.

Bistability is not required to explain the observation. The claim that rejuvenation is coherent only in a bistable system is withdrawn, and the three failures of Section 8.2 were scoring a quantity the observation never called for.

This does not rescue the model class. It relocates the limitation. What the model cannot do is place a stable attractor at a 1% growth cost, and nothing in Lindner et al. asks it to.

### 8.4 The residual test, and why it is currently unmeasurable

The partitioning model as first run had four free parameters against one target number, and a monotone one-parameter family crosses any band below its maximum. Two of them cancel analytically: with `k_μ = 0.03125/p_qc` and a model burden of 1 corresponding to a proteome fraction of `p_qc`, the scored quantity is `32 × (B_old − B_new)` as a proteome fraction exactly, so the quality-control proteome share and `μ₀` enter only through the stationary load.

The partition rule is not a free coefficient, and it is not `f = 1` either. Lindner et al. report that 52.3% of cells have no inclusion body, 46.5% have one, and 1.2% have two, and that new-pole progeny are devoid of parental inclusion bodies. The focus is a single indivisible object and goes entirely to one daughter. But the model's `a` is total aggregate, of which the visible focus is only a share. Writing `β` for the fraction of aggregate held in the focus, the old-pole daughter receives `(1+β)a` and the new-pole daughter `(1−β)a` after division into half the volume, which is exactly the scalar rule at `f_eff = (1+β)/2`. An earlier version of this analysis set `β = 1` without stating it and reported the result as parameter-free. That claim is withdrawn.

At each `β` the aging effect is `64 × a_end` as a proteome fraction times a damping factor measuring the daughter's relaxation during its own generation. That factor is itself `β`-dependent, though weakly — it is 0.355 over most of the range and falls to 0.346 as `β` approaches 1 — and is measured at each `β` rather than carried over, because under two-compartment partitioning the new-pole daughter is no longer aggregate-free and so relaxes differently. Inverting the fixed band gives a requirement that scales as `1/β`:

| `β` | `f_eff` | damping | required old-pole aggregate (% of proteome) | ratio to the Δ*rpoH* load |
|---|---|---|---|---|
| 1.00 | 1.000 | 0.3462 | 0.0467 – 0.0803 | 62× – 214× |
| 0.75 | 0.875 | 0.3552 | 0.0607 – 0.1044 | 48× – 165× |
| 0.50 | 0.750 | 0.3549 | 0.0911 – 0.1567 | 32× – 110× |
| 0.25 | 0.625 | 0.3546 | 0.1824 – 0.3137 | 16× – 55× |

The nearest available measurement does not reach any row of it. Tomoyasu et al. (2001) report 5–10% of total protein aggregated in Δ*rpoH* mutants at 30 °C and no detectable aggregation in *rpoH*⁺ cells at the same temperature, so wild type is a bound rather than a value.

No source we could find bounds `β` under unstressed growth. Winkler et al. (2010) report the ingredients separately and under heat stress, which is not the condition at issue, and combining them requires assumptions about focus number that conflict with the distribution Lindner et al. measured in unstressed cells. We therefore quote the requirement as a `β`-indexed family rather than a number, spanning 0.047% to roughly 0.5% of the proteome over the plausible range, and note the direction: the less of the aggregate that sits in the focus, the more aggregate the account requires, and the closer the prediction moves to the only load anyone has measured. Plotting the family against that measurement shows both how far the two are from meeting and how fast the gap closes (Fig. 5). Over the plotted range the requirement stays between 3× and 214× below the sole reported aggregate load, reaching that closest approach only at `β = 0.05`; at the three focal shares the measurement makes plausible it is at least fifteenfold below. No value of `β` in the plausible range puts the account in contact with an existing number.

![Figure 5](../figures/fig5.pdf)

**Fig. 5** What the account requires of an old-pole cell. The band is the
aggregate load implied by the measured lineage difference, as a fraction of total
protein, against the share `β` of that aggregate held in the visible inclusion
body. It scales as `1/β`, a straight line on these axes, because only the focal
share segregates asymmetrically. Vertical bars mark `β = 1, 0.5, 0.25`, giving
0.05–0.08 %, 0.09–0.16 % and 0.18–0.31 %. The shaded region is the only aggregate
load measured in *E. coli*, 5–10 % of total protein in Δ*rpoH* cells at 30 °C;
wild-type aggregation at the same temperature is reported as undetected, which is
a bound and not a value, so the requirement is consistent with everything measured
and constrained by none of it. No lower limit on `β` is drawn, because none is
defensible: combining Winkler et al.'s heat-stress counts to obtain one requires
assuming two foci per cell, and Lindner et al. report one focus in 46.5 % of
unstressed cells and two in 1.2 %.

One measurement closes it. Quantitative fluorescence of the polar focus against total aggregated protein by sedimentation, in the same unstressed cells, fixes `β` and collapses the table to a single interval.

## 9. Predictions and falsification

Four statements are worth separating by what would refute them.

**Mathematical, and refutable only by error.** The identity `det J = −(∇R × ∇G)`, its *n*-state form, the `μ → 0` and `σ₀ → 0` reductions, the `J − μI` form under constant dilution, the `(φ_enz, δ)` decomposition, and the square-root ceiling of Corollary 4. These follow from H1–H3 and no observation bears on them.

**Structural and untested.** Collapse occurs at a saturation fraction far below 1, against the capacity-exhaustion alternative in which collapse occurs as `s → 1`. The observable is a ratio, so testing it does not require absolute copies-per-cell calibration. Two cautions apply, both from our own results: the predicted saturation fractions are dispersed widely enough across the kinetic box that a single measurement discriminates weakly, as the spans in Fig. 3 show, and regulation widens that dispersion rather than narrowing it.

**The sharper version of the same test.** Since growth rate sets the load coordinate `ν`, and `ν` drift alone contributes ~1.6× against the ~1.8× contrast expected from shifting the chaperone/protease split, an experiment that does not hold growth rate fixed cannot discriminate. Chemostat or turbidostat, not batch culture. Under fixed growth rate, raising the load of *perfectly folding* protein should lower the tolerable damage influx, because nascent chains consume chaperone capacity without contributing influx; a model in which the two loads are handled independently predicts no shift. Direction was correct in 67 of 68 settings over a 100-fold ladder, with median shift 1.22×.

**A ruler, not a test.** Approaching the fold, the relaxation time scales as `τ ~ (j_crit − j)^(−1/2)` (fitted slope 0.5077, *r*² = 1.0000). This is generic to saddle-nodes, so confirming it selects this model over nothing in particular. Its use is as an instrument for locating the boundary in an experiment that then asks whether the boundary *moves*.

**The measurement the theory most wants.** Whether growth arrest under misfolding burden is complete or merely asymptotic. Under the hyperbolic form, dilution keeps bounding the high-burden state and collapse is a reversible switch; under linear arrest, `μ` reaches exactly zero and the high state is a runaway. The available dosage-response constrains the slope only below 0.1% misfolded protein, roughly three decades below the arrest burden where the two forms diverge, so it cannot adjudicate. That single functional form decides whether losing proteostasis in a dividing cell is recoverable.

## 10. Discussion

The threshold of a finite-resource clearance network is not something one has to find by sweeping. Given the removal flux and the aggregate nullcline, it is where their gradients become parallel, and the critical influx is the removal flux evaluated there. That statement requires no fitted parameter, holds for any number of states, survives growth dilution, and survives capacity that degrades with the load it carries.

What it does not do is predict a number without measured parameters. Two attempts to make the quantitative predictions sharp have failed — calibration against a measured growth-burden relation, then the addition of regulation — while the structural core survived every extension attempted. This is a structural theory rather than a predictive one, and claims should be pitched accordingly.

The negative results are the more transferable half. Three mechanisms were added to satisfy a structural requirement that nobody had checked, and the check ran in an afternoon against three cycles of mechanism-building. The requirement dissolved. Before adding a mechanism to reproduce an observation, it is worth asking whether the observation imposes the structure the mechanism is meant to supply.

The nearest prior model is Santra, Dill and de Graff (2019), which drives proteostasis collapse by chaperone titration under cumulative damage and locates a tipping point beyond which replenishment no longer keeps pace with loss. The qualitative claim is the same one made here, and the difference is worth stating precisely rather than as a claim of novelty. That model reaches its tipping point by simulation and fits lifespan data, which this paper does not attempt; against that, the threshold here is a condition on two gradients, so it is located without sweeping, holds at any state dimension, and survives dilution and self-damaging capacity without re-derivation.

The substantive disagreement is about mechanism. Titration accounts place collapse at the point where chaperone capacity is exhausted. Section 6 finds the opposite: at the boundary the machinery runs at 6 to 18 percent of its maximum rate, the saturation deficit accounts for 35.8 percent of the shortfall against 12.6 percent for sequestration, and the capacity bound is about thirteen times too loose to locate anything on its own. Collapse in this model class occurs because superlinear nucleation overtakes sublinear saturating removal well before removal saturates. If that is right, interventions that add capacity buy less than a titration picture predicts, and the two accounts are separable by measuring how far from saturation the clearance machinery runs when a cell loses proteostasis.

Two limitations bound everything above. Deterministic well-mixed dynamics is false in the system we scored against — aggregates are spatially segregated and molecule numbers are small — and both effects place a measured boundary below the deterministic fold and smear it across a population. And the parameter box was chosen adversarially wide in the saturation constants, which set `φ`, so the 8.86× between-network spread is partly a property of the sampling rather than of biology.

The prediction that follows from Section 8.4 is the one a reader can act on. If the aggregate burden of an aging *E. coli* old-pole cell falls outside the interval that Section 8.4 requires at the measured focal share, the account given here of the aging asymmetry is wrong. Two measurements decide it, and neither has been attempted as far as we can determine: the wild-type aggregate fraction under unstressed growth, and the share of that aggregate held in the polar focus.

---

## Supplementary figure

![Figure S1](../figures/figS1.pdf)

**Fig. S1** The parallelism residual is bracket tolerance, not failure of the
identity. Over all 325 folds of the load grid, `sin θ` between `∇R` and `∇G`
tracks the leading eigenvalue recorded at the same state across five decades, with
a log–log correlation of +0.9960 and a slope of 1.00 — the residual is
*proportional* to how loosely the state is bracketed, which is what bracket
tolerance means and what a failing identity would not produce. The circled point
is the tightest bracket in the population, `|λ| = 8.10×10⁻¹⁰` with
`sin θ = 7.75×10⁻⁹`. Normalisation, reported here rather than given a panel of its
own: on these same 325 states the max-normalised residual correlates with `|λ|` at
−0.835 and reaches 1.54×10⁻², while the gradient-normalised residual used
throughout correlates at +0.041 and reaches 1.29×10⁻⁹. The first degrades as the
bracket tightens on the object it is meant to verify, because both `det J` and the
cross product vanish at a saddle-node and their ratio is 0/0 there; the second does
not. Every number in this caption is recomputed by
`scripts/figures/fig_identity.py:captionNumbers` and asserted by the test suite.

---

## Data and code availability

All analysis code, the parameter configurations, and the test suite that asserts every numerical quantity reported here are at https://github.com/khatvangi/proteostasis-law-theory under the MIT licence. Archived version and DOI: https://doi.org/10.5281/zenodo.21794565 (concept DOI, resolving to the latest version). Every figure is regenerated by a script in `scripts/figures/`, each of which reads a reduced array carrying its own provenance — population, count, and whether it is a subsample — so the figures rebuild on a clean checkout without the run root. Detailed run outputs are regenerable from the tracked scripts; without the run root, artefact-dependent scripts print an explicit SKIP and exit 0 rather than passing silently.

## References

Antoneli F, Golubitsky M, Jin J, Stewart I (2025) Homeostasis in input-output networks: structure, classification and applications. *Math Biosci*. arXiv:2405.03861

Araujo RP, Liotta LA (2018) The topological requirements for robust perfect adaptation in networks of any size. *Nat Commun* 9:1757.

Brualdi RA, Ryser HJ (1991) *Combinatorial matrix theory*. Cambridge University Press, Cambridge.

Dear AJ, Meisl G, Michaels TCT, Zimmermann MR, Linse S, Knowles TPJ (2020) The catalytic nature of protein aggregation. *J Chem Phys* 152(4):045101. doi:10.1063/1.5133635

Cohen SIA, Vendruscolo M, Welland ME, Dobson CM, Terentjev EM, Knowles TPJ (2011) Nucleated polymerization with secondary pathways. I. Time evolution of the principal moments. *J Chem Phys* 135:065105.

Cohen SIA, Linse S, Luheshi LM, Hellstrand E, White DA, Rajah L, Otzen DE, Vendruscolo M, Dobson CM, Knowles TPJ (2013) Proliferation of amyloid-β42 aggregates occurs through a secondary nucleation mechanism. *Proc Natl Acad Sci USA* 110:9758–9763.

Dekel E, Alon U (2005) Optimality and evolutionary tuning of the expression level of a protein. *Nature* 436:588–592.

Ferrell JE (2016) Perfect and near-perfect adaptation in cell signaling. *Cell Syst* 2:62–67.

Ferrone F (1999) Analysis of protein aggregation kinetics. *Methods Enzymol* 309:256–274.

Frere S, Slutsky I (2018) Alzheimer's disease: from firing instability to homeostasis network collapse. *Neuron* 97(1):32–58. doi:10.1016/j.neuron.2017.11.028

Geiler-Samerotte KA, Dion MF, Budnik BA, Wang SM, Hartl DL, Drummond DA (2011) Misfolded proteins impose a dosage-dependent fitness cost and trigger a cytosolic unfolded protein response in yeast. *Proc Natl Acad Sci USA* 108(2):680–685. doi:10.1073/pnas.1017570108

Golubitsky M, Stewart I (2006) Nonlinear dynamics of networks: the groupoid formalism. *Bull Amer Math Soc* 43(3):305–364.

Golubitsky M, Stewart I (2017) Homeostasis, singularities and networks. *J Math Biol* 74:387–407. doi:10.1007/s00285-016-1024-2

Golubitsky M, Wang Y (2020) Infinitesimal homeostasis in three-node input-output networks. *J Math Biol* 80:1163–1185. doi:10.1007/s00285-019-01457-x

Huang Z, Golubitsky M (2022) Classification of infinitesimal homeostasis in four-node input-output networks. *J Math Biol* 84:25. doi:10.1007/s00285-022-01727-1

Hipp MS, Kasturi P, Hartl FU (2019) The proteostasis network and its decline in ageing. *Nat Rev Mol Cell Biol* 20(7):421–435. doi:10.1038/s41580-019-0101-y

Khammash M (2021) Perfect adaptation in biology. *Cell Syst* 12:509–521. doi:10.1016/j.cels.2021.05.020

Knowles TPJ, Waudby CA, Devlin GL, Cohen SIA, Aguzzi A, Vendruscolo M, Terentjev EM, Welland ME, Dobson CM (2009) An analytical solution to the kinetics of breakable filament assembly. *Science* 326:1533–1537.

Lin X, Antoneli F, Wang Y (2026) Automated classification of homeostasis structure in input-output networks. *Bull Math Biol*, article 113. arXiv:2603.08882

Lindner AB, Madden R, Demarez A, Stewart EJ, Taddei F (2008) Asymmetric segregation of protein aggregates is associated with cellular aging and rejuvenation. *Proc Natl Acad Sci USA* 105(8):3076–3081. doi:10.1073/pnas.0708931105

Ma W, Trusina A, El-Samad H, Lim WA, Tang C (2009) Defining network topologies that can achieve biochemical adaptation. *Cell* 138:760–773. doi:10.1016/j.cell.2009.06.013

Madeira JLO, Antoneli F (2024) Homeostasis in networks with multiple inputs. *J Math Biol* 89:17. doi:10.1007/s00285-024-02117-5

Meisl G, Kirkegaard JB, Arosio P, Michaels TCT, Vendruscolo M, Dobson CM, Linse S, Knowles TPJ (2016) Molecular mechanisms of protein aggregation from global fitting of kinetic models. *Nat Protoc* 11(2):252–272.

Nijhout HF, Best J, Reed MC (2014) Escape from homeostasis. *Math Biosci* 257:104–110. doi:10.1016/j.mbs.2014.08.015

Reed M, Best J, Golubitsky M, Stewart I, Nijhout HF (2017) Analysis of homeostatic mechanisms in biochemical networks. *Bull Math Biol* 79:2534–2557. doi:10.1007/s11538-017-0340-z

Reed MC, Duncan W, Nijhout HF, Best J, Golubitsky M (2018) Homeostasis despite instability. *Math Biosci* 300:130–137. doi:10.1016/j.mbs.2018.03.025

Santra M, Dill KA, de Graff AMR (2019) Proteostasis collapse is a driver of cell aging and death. *Proc Natl Acad Sci USA* 116(44):22173–22178. doi:10.1073/pnas.1906592116

Schmidt A, Kochanowski K, Vedelaar S, Ahrné E, Volkmer B, Callipo L, Knoops K, Bauer M, Aebersold R, Heinemann M (2016) The quantitative and condition-dependent *Escherichia coli* proteome. *Nat Biotechnol* 34(1):104–110. doi:10.1038/nbt.3418

Scott M, Gunderson CW, Mateescu EM, Zhang Z, Hwa T (2010) Interdependence of cell growth and gene expression: origins and consequences. *Science* 330:1099–1102.

Scott M, Klumpp S, Mateescu EM, Hwa T (2014) Emergence of robust growth laws from optimal regulation of ribosome synthesis. *Mol Syst Biol* 10:747.

Scott M, Hwa T (2023) Shaping bacterial gene expression by physiological and proteome allocation constraints. *Nat Rev Microbiol* 21(5):327–342. doi:10.1038/s41579-022-00818-6

Shachrai I, Zaslaver A, Alon U, Dekel E (2010) Cost of unneeded proteins in *E. coli* is reduced after several generations in exponential growth. *Mol Cell* 38:758–767.

Stewart EJ, Madden R, Paul G, Taddei F (2005) Aging and death in an organism that reproduces by morphologically symmetric division. *PLoS Biol* 3(2):e45. doi:10.1371/journal.pbio.0030045

Tomoyasu T, Mogk A, Langen H, Goloubinoff P, Bukau B (2001) Genetic dissection of the roles of chaperones and proteases in protein folding and degradation in the *Escherichia coli* cytosol. *Mol Microbiol* 40(2):397–413. doi:10.1046/j.1365-2958.2001.02383.x

Wang Y, Huang Z, Antoneli F, Golubitsky M (2021) The structure of infinitesimal homeostasis in input-output networks. *J Math Biol* 82. doi:10.1007/s00285-021-01614-1

Winkler J, Seybert A, König L, Pruggnaller S, Haselmann U, Sourjik V, Weiss M, Frangakis AS, Mogk A, Bukau B (2010) Quantitative and spatio-temporal features of protein aggregation in *Escherichia coli* and consequences on protein quality control and cellular ageing. *EMBO J* 29(5):910–923. doi:10.1038/emboj.2009.412

---

### Reference verification

All thirty-five references above were checked against PubMed records, not against
secondary reference lists. Four corrections resulted.

- *The catalytic nature of protein aggregation* is **Dear AJ** et al., not
  Michaels et al.; Michaels is the third author. Article number 045101, confirmed.
- Frere & Slutsky (2018) runs **97(1):32–58**; the end page was previously unknown.
- Stewart et al. (2005) is *PLoS Biol* **3(2):e45**, disambiguated from another
  2005 *PLoS Biol* paper with a Stewart among its authors.
- Scott & Hwa carries **doi:10.1038/s41579-022-00818-6** and appeared online in
  November 2022; it is cited by its 2023 issue.

Coquel et al. (2013) was verified but is **not cited**: it reports aggregate
diffusion constants and the nucleoid-crowding mechanism, and does not report the
focal mass share that §8.4 would need it for.

Two gaps identified earlier are now closed. Hipp, Kasturi & Hartl (2019) supports
the finite-shared-resource premise of §2.1, and Schmidt et al. (2016) supplies the
quantitative *E. coli* proteome from which the chaperone and protease share used in
§4.2 can be obtained rather than swept. Santra, Dill & de Graff (2019) is added as
the closest prior model of proteostasis collapse: it drives collapse by chaperone
titration under cumulative damage and reaches a tipping point, which is the same
qualitative claim this paper makes exactly rather than numerically.
