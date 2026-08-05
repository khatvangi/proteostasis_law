# An Exact Fold Condition for Mass-Balanced Models of Protein Quality Control

**Kiran Boggavarapu**

Department of Chemistry and Physics, McNeese State University, Lake Charles, LA 70609, USA
kiran@mcneese.edu
ORCID 0000-0003-0751-6459

**MSC** 92C40 (primary); 37G10, 34C23, 92C42 (secondary)

**Keywords** proteostasis; saddle-node bifurcation; conserved resource; clearance network; growth dilution; Hopf bifurcation

---

## Abstract

A cell clearing misfolded protein with finite pools of chaperones and proteases tolerates damage influx only up to a threshold, beyond which burden runs away. Where that threshold sits is normally found by numerical continuation, one parameter set at a time. We show that it is algebraic. For any clearance network in which the damage influx is state-independent and mass balance accounts for all outflow, the Jacobian determinant factors as `det J = −(∇R × ∇G)`, where `R` is total removal flux and `G` the aggregate nullcline field. A fold therefore occurs exactly at a constrained critical point of removal on the aggregate nullcline, with critical influx `j_crit = R(u*, a*)`; the converse holds under four genericity conditions which we state and verify. The identity generalises to arbitrary state dimension, survives growth dilution, and survives clearance machinery that is itself damaged by the influx it clears. Three consequences follow. The fold occurs while the machinery runs at roughly 5–18% of maximum velocity, so the capacity ceiling overestimates the tolerable influx about twelvefold. Under cell division the ceiling ceases to bound, and the critical influx decomposes exactly as `j_crit = C_enz·φ_enz/(1−δ)` with the enzymatic factor nearly invariant to division rate. When the machinery damages itself, a linear capacity ceiling becomes a square-root one. The fold is not always the first loss of stability: in a wide kinetic ensemble, 3.8% of networks lose stability through a Hopf bifurcation first, in a region characterised by sharply saturating aggregate clearance against growth-dominated aggregation. Applied to the aggregation asymmetry of dividing *Escherichia coli*, the framework yields a falsifiable requirement on the old-pole aggregate load as a function of the share held in the visible polar focus.

---

## 1. Introduction

Protein quality control is a resource-limited clearance problem. Chaperones refold damaged monomers, proteases degrade them, and both draw on pools that are finite and shared across substrates (Hipp et al. 2019). Damage arrives continuously from translation error, thermal stress and oxidative modification. When influx exceeds what the machinery can dispose of, aggregate accumulates and the cell loses viability.

Models built on this premise locate the threshold numerically. One fixes a parameter set, integrates, continues the stable branch in the influx and records where it terminates. Cotton et al. (2025) illustrate both the power and the cost of this approach: across eight variants of aggregate formation and removal kinetics they demonstrate that the same bifurcation structure is preserved, with a stable low-aggregate state and an unstable upper branch merging at a critical monomer concentration. The universality is established by exhibiting it case by case. Nothing in that treatment says why it should hold, or what would have to change for it to fail.

We show that for one structural class the threshold is algebraic rather than computational. Two facts about clearance models — the damage influx enters one state equation, and mass balance is exact — force the Jacobian determinant to factor into the cross product of two gradients. A fold is then a constrained critical point of total removal, in the Lagrange sense, on the nullcline of the non-influx equation, and the critical influx is the value of removal at that point. Adding states, adding regulation, adding cell division, or making the clearance capacity itself depend on the load leaves the identity intact.

The hypothesis required is weaker than it first appears. The theorem does not need influx and clearance capacity to be independent. It needs only that total influx be state-independent and that mass balance count all outflow. Capacity may depend on the error rate, on the burden, or on both. This matters because in a cell the chaperones and proteases are themselves translated by the error-prone apparatus whose product they clear.

Two results bound what the fold condition delivers. Solving the fold equations returns the candidate set exactly and without a sweep, but in about 9% of a wide kinetic ensemble that set has more than one element, and identifying which candidate terminates the accessible branch requires information the equations do not carry. And the fold is not always the first loss of stability: in the same ensemble a Hopf bifurcation precedes it in a characterisable region of parameter space.

Section 2 fixes the model class and its hypotheses. Section 3 states and proves the fold condition, gives the genericity conditions for the converse, places the result relative to the homeostasis literature, and treats orientation and multiplicity. Section 4 gives three extensions. Section 5 reports numerical verification. Section 6 characterises where the fold sits relative to the capacity ceiling. Section 7 treats instabilities that precede the fold. Section 8 applies the framework to aggregation asymmetry in dividing *E. coli*. Section 9 states what the theory predicts and what would refute it.

## 2. The model class

### 2.1 State equations

Let `u` denote damaged or misfolded soluble monomer and `a` aggregate burden. Both are **total** pools — free plus machinery-bound — carried in monomer-equivalent units, so that transfer between them conserves mass one for one. The baseline field is

```
du/dt = j − v_ref − v_degU − n − g + v_dis
da/dt =              n + g − v_dis − v_degA
```

Chaperone and protease pools are conserved and shared, and every rate law is written in free concentrations. With 1:1 rapid-equilibrium complexes `CU`, `CA`, `CN`, `DU` and `DA`, the four conservation laws close simultaneously on the free concentrations:

```
u_f = u /(1 + C_f/K_CU + D_f/K_DU)     C_f = C_tot /(1 + ν + u_f/K_CU + a_f/K_CA)
a_f = a /(1 + C_f/K_CA + D_f/K_DA)     D_f = D_tot /(1 + u_f/K_DU + a_f/K_DA)
```

Given rapid equilibrium these four equations are identities rather than approximations, and they are solved jointly at every state; total substrate is never substituted into a formula expecting a free concentration. Here `ν = N_f/K_N` is the occupancy contributed by nascent chains, which consume capacity without contributing damage influx — the mechanism by which ordinary translational load can move the fold.

`R` and `G` are therefore defined implicitly rather than by formulae in `(u, a)`, and their regularity has to be established rather than assumed.

**Lemma 0 (well-posedness of the closure).** *For every `(u, a)` in the nonnegative quadrant the four equations above have exactly one solution, and it is strictly positive. The binding Jacobian in `(C_f, D_f)` is a nonsingular M-matrix with `det ≥ 1 + ν`, so the free concentrations are real-analytic in `(u, a)` on the open quadrant. Consequently `∂u_f/∂u > 0` and `∂a_f/∂u ≥ 0`.*

*Proof.* The map sending `(C_f, D_f)` to the right-hand sides of the two pool equations carries the box `[0, C_tot] × [0, D_tot]` into itself and is monotone increasing, since raising either free pool lowers both free substrates and so lowers both denominators; Knaster–Tarski gives a least and a greatest fixed point. It is also strictly sublinear — scaling the argument by `λ < 1` shrinks the occupancy terms of each denominator by `λ` but leaves the additive constant, so the quotient rises by more than `λ` — and a monotone strictly sublinear map has at most one positive fixed point, by the standard argument on the largest `t` with `x ≥ t y`. For the Jacobian, all four partials of the free substrates in the free pools are negative, so the off-diagonal entries are nonpositive. Contracting against the positive vector `(C_f, D_f)` and using that each free substrate is homogeneous of degree `−1` in `(1, C_f, D_f)` gives row values `C_f(1 + ν + s_u^{-1}u_f/K_CU + s_a^{-1}a_f/K_CA)` and `D_f(1 + s_u^{-1}u_f/K_DU + s_a^{-1}a_f/K_DA)`, both positive, where `s_u` and `s_a` are the two substrate denominators. A matrix with nonpositive off-diagonal entries admitting a positive vector of positive image is a nonsingular M-matrix, whence `det ≥ 1 + ν` and the inverse is nonnegative. The implicit function theorem then applies at every state, and inverse-positivity applied to the derivative of the closure in `u` gives the two monotonicity statements. ∎

The rate laws are

```
v_ref  = k_ref C_f u_f/(K_ref + u_f)      refolding
v_degU = k_U D_f u_f/(K_U + u_f)          soluble degradation
n      = k_n u_f^m,  m > 1                primary nucleation
g      = k_g u_f a_f                      elongation
v_dis  = k_dis C_f a_f/(K_dis + a_f)      disaggregation
v_degA = k_A D_f a_f/(K_A + a_f)          aggregate clearance
```

The two aggregate-forming terms are distinct processes. `n = k_n u_f^m` is primary nucleation with an effective reaction order greater than one, in the sense in which fitted orders are understood as coarse-grained descriptions of multi-step processes (Ferrone 1999; Knowles et al. 2009; Meisl et al. 2017). `g = k_g u_f a_f` is elongation: growth of existing aggregate at a monomer-dependent rate. The model carries no secondary nucleation term, which would be autocatalytic in aggregate surface (Cohen et al. 2013); adding one changes the coefficients below but not the structure, since it enters `G` and not `R`.

The dissociation constants `K_CU` and `K_CA` describe competitive occupancy of the pool; the Michaelis constants `K_ref`, `K_U`, `K_dis` and `K_A` describe a rate-limiting catalytic step acting on the pool not held in any complex. Treating them as successive steps is a stipulation and not a derivation. Written against the bound complex instead, each flux would carry no second saturation, so the coded form asserts nonproductive sequestration followed by independent productive capture, and gives cycles sharing one pool independent Michaelis factors rather than one competitive denominator. Both stipulations bind. Over the kinetic ensemble of Section 5 the two catalytic occupancies drawing on the chaperone pool sum above unity at 334 of 2767 fold states, and the substrate that productive complexes would hold reaches a fifth of the total soluble pool at the median. Section 4.4 reports what moves under the alternative closure. Nothing in Theorem 1 depends on the choice: the hypotheses below hold for either.

Write the total removal flux and the aggregate field as

```
R(u,a) = v_ref + v_degU + v_degA        G(u,a) = da/dt
```

`R` is what leaves the damaged compartment; `G` is the non-influx equation. Note that `n` and `g` appear in `R` nowhere: they transfer material between `u` and `a` and cancel from `du/dt + da/dt`. Aggregation does not oppose removal directly. It reshapes the nullcline, redistributes burden between the two pools, and suppresses effective removal through saturation and through sequestration of the shared machinery. Section 6 separates those two contributions.

### 2.2 Hypotheses

**H1 (state-independent influx).** The damage influx `j` is a parameter appearing additively in `du/dt`. Two consequences are used separately and are worth separating. **(H1a)** Total influx does not depend on the state, which is what puts mass balance in the form H2 requires. **(H1b)** `j` does not appear in `da/dt`, so the curve `{G = 0}` does not move with the load. Only (H1a) enters the identity. Section 4.3 relaxes (H1b) by letting clearance capacity depend on the influx, and Theorem 1.4 covers that case.

**H2 (exact mass balance).** Internal transfer between `u` and `a` cancels between the two equations, so `du/dt + da/dt = j − R` holds identically.

**H3 (smoothness).** `R` and `G` are continuously differentiable on the region of interest. For the model class of Section 2.1 this is a consequence rather than an assumption: Lemma 0 makes the free concentrations real-analytic in the state, and every rate law is rational in them apart from `k_n u_f^m`, which is `C¹` at the origin for `m > 1`. H3 is retained as a hypothesis because Theorem 1 is stated for the wider class in which the fields are given directly.

**H4 (deterministic, well-mixed, finite-dimensional).** The system is an ODE in concentrations.

*Remark 1.* H1 constrains the influx, not the capacity. Nothing forbids `C_tot` or `D_tot` from depending on `j` or on the state. The gradients in Section 3 are taken at fixed `j`, so parameter dependence on the influx cannot enter them; the row operation needs (H1a) and H2 alone, and neither mentions `da/dt`. What capacity feedback costs is (H1b), and with it the fixed nullcline that makes the candidate solve cheap. Corollary 4 treats the consequences of capacity that degrades with load.

*Remark 2.* H4 is false in the system of Section 8 in two respects. Aggregates in *E. coli* are spatially segregated into polar inclusion bodies rather than well mixed, and molecule numbers per damaged species are small. Both effects act in the same direction: noise-induced escape precedes the deterministic fold, so a measured threshold should sit below `j_crit` and be smeared across a population rather than sharp.

## 3. The fold condition

### 3.1 The identity

Four statements are separated because they need different hypotheses; they are referred to collectively as Theorem 1.

**Theorem 1.1 (mass-balance identity).** *Assume (H1a), H2 and H3. Then at fixed `j`*

```
det J = −(∇R × ∇G)
```

*identically, the gradients being taken with respect to the state.*

*Proof.* By H2, `du/dt + da/dt = j − R` exactly. The determinant-preserving row operation `row₁ → row₁ + row₂` replaces the first row by the gradient of `j − R`, and `j` is a parameter, so

```
det J = det [ ∂(du/dt)/∂u   ∂(du/dt)/∂a ]  =  det [ −R_u  −R_a ]
            [ ∂(da/dt)/∂u   ∂(da/dt)/∂a ]        [  G_u   G_a ]

      = −(R_u G_a − R_a G_u) = −(∇R × ∇G). ∎
```

(H1b) is not used, and the identity is indifferent to whether `j` enters `G` or the capacity depends on the load, because the gradients are taken at fixed influx.

**Theorem 1.2 (removal along the nullcline).** *Assume in addition that `∇G ≠ 0`. Parametrise `{G = 0}` by arclength `s` with unit tangent `γ' = (−G_a, G_u)/‖∇G‖`, and put `r(s) = R(γ(s))`. Then*

```
r'(s) = ∇R·γ' = det J /‖∇G‖
```

*at every point of the nullcline, not only where it vanishes.*

`det J` **is** the derivative of total removal along the aggregate nullcline, up to the positive factor `‖∇G‖`. The informal reading — the fold is where total removal stops responding to burden along the nullcline — is therefore an exact pointwise identity rather than a gloss on a vanishing condition. Two things follow at once. The equilibria on the nullcline are the solutions of `r(s) = j`, so the equilibrium branch is the level-set inverse of one scalar function of one variable; and `det J = 0` is `r'(s) = 0`, a constrained critical point of `R` on `{G = 0}`, with critical value `j_crit = R(u*, a*)`.

Regularity of the constraint is a hypothesis of the Lagrange reading and not only of the converse. Where `∇G = 0` the determinant vanishes automatically, `{G = 0}` need not be a manifold, and there is no tangent to differentiate along; the two networks discussed in Section 3.2 are that case, not counterexamples to anything.

**Theorem 1.3 (fixed-nullcline shortcut).** *Assume (H1b), so `∂G/∂j = 0`. Then `{G = 0}` is independent of the influx — changing the load slides the state along it — equilibria are exactly its points with `R = j`, and the fold candidates solve `{G = 0, det J = 0}` in the state alone.*

This is the whole practical content: neither equation contains `j`, so locating candidates is a two-dimensional root solve rather than a continuation sweep.

**Theorem 1.4 (load-dependent capacity).** *When `∂G/∂j ≠ 0` the nullcline moves with the influx and `j_crit = R(u*,a*)` becomes a self-consistency condition; candidates solve `{G = 0, det J = 0, R = j}` in `(u, a, j)`. Theorem 1.1 is unaffected.*

At the base parameter set the alignment is visible rather than inferred: the removal contours meet the aggregate nullcline tangentially at the solved fold, to a sine of 3.54×10⁻¹⁰ (Fig. 1a). The same computed branch separates two singularities that are easily conflated. The soluble coordinate has a horizontal tangent at `j_turn = 0.154090`, where `G_a = 0` and the nullcline runs vertical in the burden plane; that point is regular, since `det J = R_a·G_u = 2.027×10⁻³` there. The fold is at `j_crit = 0.154239`, where both coordinates turn together (Fig. 1b). The two lie in a fixed order and differ by one part in a thousand. This is also why tracing `{G = 0}` by solving for `a` at fixed `u` fails: the two roots in `a` merge exactly at `j_turn`, returning the curve in two disconnected pieces with the fold in the gap between them. The root-finder's failure and the horizontal tangent are one locus, not two problems.

![Figure 1](../figures/fig1.pdf)

**Fig. 1** The fold as a constrained critical point. **(a)** Phase plane at the base parameter set. The aggregate nullcline `{G = 0}` (dark) is fixed under change of influx; grey contours are total removal `R`. At the solved fold `(u*, a*) = (0.4166, 0.2650)` the gradients `∇R` and `∇G` are parallel, which is the Lagrange condition of Theorem 1.2; they are drawn from offset origins because a common origin would hide one behind the other. The sine of the angle between them is 3.54×10⁻¹⁰, and `j_crit = R(u*, a*) = 0.1542`. **(b)** Both equilibrium coordinates against influx along the same curve. Cramer's rule gives `du*/dj = −G_a/det J` and `da*/dj = G_u/det J`, so the branch carries two distinct, ordered singularities. At `j_turn = 0.154090` the soluble coordinate has a horizontal tangent where `G_a = 0`; this is a regular point, `det J = R_a·G_u = 2.027×10⁻³` there rather than zero, and it lies on the stable branch at `j_turn/j_crit = 0.9990`. At `j_crit = 0.154239` both coordinates have a vertical tangent, which is `det J = 0`. The inset is a ×182 zoom on the same computed data, not a schematic. The aggregate coordinate has no horizontal tangent anywhere: `G_u > 0` at all 305 traced points, minimum 2.077×10⁻³.

### 3.2 Genericity, and the converse

Theorem 1.2 gives one direction: a fold satisfies gradient alignment. The converse requires that the degeneracy be a fold rather than a transcritical or pitchfork bifurcation, a cusp, a double-zero eigenvalue, or a curve of equilibria. Four conditions suffice:

- **(G1)** regularity of the constraint, `∇G ≠ 0`;
- **(G2)** a simple zero eigenvalue, `tr J ≠ 0`;
- **(G3)** nondegeneracy, `d²R/ds² ≠ 0` along the nullcline at the critical point.

Parameter transversality, `w·∂F/∂j ≠ 0` for the left null vector `w`, is the fourth condition of the classical statement. It is not independent here, and the corollary below shows it follows from (G1).

**Theorem 2 (converse).** *Let `F` be `C²`. Under (G1)–(G3), a constrained critical point of `R` on `{G = 0}` is a generic saddle-node in `j`, with `j_crit = R(u*, a*)`.*

*Proof.* By (G1) the nullcline is a `C²` curve near the point, so Theorem 1.2 applies: equilibria on it are the solutions of `r(s) = j`, and `det J = 0` is `r'(0) = 0`. With (G3),

```
r(s) = j_crit + ½ r''(0) s² + O(s³),
```

a nondegenerate extremum, so for `j` on one side of `j_crit` there are exactly two equilibria near the point and on the other side none: the pair collides and annihilates. By (G2) the zero eigenvalue is algebraically simple and the remaining eigenvalue `tr J` lies off the imaginary axis, so the centre manifold is one-dimensional. It is the nullcline itself — the second row of `Jv = 0` reads `∇G·v = 0`, so the right null vector is tangent to `{G = 0}` and parallel to `γ'(0)`; the centre manifold and the arclength reduction are the same object. Finally the two nondegeneracy conditions coincide. With `w = (1, 1+λ)` the left null vector and `∇R = λ∇G`,

```
w·D²F(v,v) = −D²R(v,v) + λ D²G(v,v),
```

while differentiating `G(γ(s)) = 0` twice gives `∇G·γ'' = −D²G(v,v)`, so that `r''(0) = D²R(v,v) − λ D²G(v,v)`. Hence

```
r''(0) = −w·D²F(v,v),
```

and (G3) is not merely equivalent to the classical nondegeneracy condition of the saddle-node but equal to it up to sign. The reduced equation on the centre manifold is therefore the saddle-node normal form with quadratic coefficient of sign `−r''(0)`. ∎

Two things come free. The sign of `r''(0)` is the orientation of Section 3.4, so that classification is a corollary rather than a separate observation. And nothing in the proof appeals to numerical verification; Table 1 reports whether the hypotheses hold in the ensembles studied, which is a question about the model rather than about the theorem.

*Remark 3 (what a fold terminates).* Step two of the proof also classifies the collision. Since the eigenvalues at a fold are `0` and `tr J`, a negative trace means a stable node and a saddle merge, so the fold terminates a stable branch and is the collapse threshold. A positive trace means an unstable node and a saddle merge, and the low-burden state has already lost stability somewhere below: the fold is a turning point of the branch but not a boundary of viability. On the load grid every fold has `tr J < 0`. Over the kinetic box 61 of 2767 have `tr J > 0`, and all 61 are among the networks of Section 7 whose branch loses stability before the fold — as they must be, which makes the count a consistency check between two independent computations rather than a finding. The finding is the converse: of the 108 networks that lose stability first, 47 recover it before reaching the fold, so an oscillatory excursion below `j_crit` does not imply that the fold itself is unstable.

**Corollary (transversality is automatic).** *Under (G1), `w·∂F/∂j ≠ 0`.*

*Proof.* Since `∂F/∂j = (1, 0)` exactly, the condition is `w₁ ≠ 0`. Let `L = I + e₁e₂^T`, so that `LJ` is the matrix whose first row is `−∇R` and whose remaining rows are the non-influx constraint gradients, and note that `L` fixes the first coordinate of any row vector multiplying it on the left. Put `z = wL⁻¹`, so `z(LJ) = 0`. If `z₁ = 0` the remaining components annihilate the non-influx constraint gradients, which have full rank by (G1), forcing `z = 0` and hence `w = 0`. So `z₁ ≠ 0`, and `w₁ = z₁ ≠ 0`. ∎

The argument uses only `∂F/∂j = e₁` and full rank, so it holds in every state dimension. It also gives the margin in closed form: normalising `w = (1, 1+λ)` with `∇R = λ∇G`,

```
|w₁|/‖w‖ = (1 + (1+λ)²)^(−1/2),
```

which is small exactly when `|λ|` is large — that is, when `∇G → 0` at fixed `∇R`. The near-failures of the constraint and of transversality are therefore one event and not two, which is what Table 1's two exceptions are.

Computing the two sides of the closed form by independent routes — the left from a singular value decomposition of the analytic Jacobian, the right from central-difference gradients — they agree at all 325 load-grid states and at 2765 of the 2767 kinetic-box states, at median relative difference 8.5×10⁻¹¹ and maximum 1.8×10⁻⁵. The two that disagree are precisely the two exceptions of Table 1, where `‖∇G‖` is 1.04×10⁻⁸ and 3.99×10⁻⁸ against a smallest value of 1.41×10⁻⁶ everywhere else: there the gradient ratio is differencing noise and the closed form has nothing to evaluate.

**Table 1** The genericity conditions of Theorem 2, evaluated at every solved fold state in each ensemble. The first three rows are the hypotheses; each entry is the minimum of that condition's margin over the whole population, so a value bounded away from zero means the condition holds throughout. The fourth row is a consequence of the first rather than an independent condition, and takes the closed form given above wherever the constraint is regular; its two exceptions are the same two networks as (G1)'s, which is forced. The last row reports a sign rather than a margin, because `\|tr J\|` conceals the classification: the eigenvalues at a fold are `0` and `tr J`, so a negative trace is a stable node colliding with a saddle and the fold terminates a stable branch, while a positive trace is an unstable node colliding with a saddle and the fold is not a stability boundary at all.

| condition | load grid (325) | kinetic box (2767) |
|---|---|---|
| (G1) `\|∇G\|` | min 0.106 | min 1.04×10⁻⁸, 2 exceptions |
| (G2) `\|tr J\|` | min 0.303 | min 8.38×10⁻⁴ |
| (G3) `\|d²R/ds²\|` | min 0.0929 | min 9.21×10⁻⁶ |
| consequence of (G1), `\|w₁\|/\|w\|` | min 0.341 | min 3.34×10⁻⁹, 2 exceptions |
| sign of `tr J` | all negative, −0.365 to −0.303 | −626 to +3.41; positive at 61 |

The two exceptions are the same two networks, and they fail (G1) and (G3) for one reason, taking the transversality margin down with them: both sit at aggregate burden of order 10⁻¹¹ and 10⁻¹⁵, where the aggregate compartment is empty, `∇G → 0`, the constraint ceases to be a manifold, and every condition defined on it degenerates. At both, `tr J = −∇R` exactly, the signature of a state where `G` has gone flat. These are the boundary of the model's domain of definition rather than counterexamples to Theorem 2. A further 117 of 2884 sampled networks carry no solvable state satisfying `{G = 0, det J = 0}` and are excluded from all conditions above.

Both computations use arclength along the nullcline rather than the `u`-coordinate. The nullcline turns vertical where `G_a = 0`, which is a coordinate singularity of the `u`-parametrisation, and every quantity evaluated on the curve is affected there.

### 3.3 Relation to infinitesimal homeostasis

The hypotheses of Theorem 1 nearly coincide with those of an existing programme, and the relation is worth stating precisely.

Wang et al. (2021) study input-output networks in which a parameter enters exactly one node's equation — the condition H1 restates — and ask when the input-output function `x_o(I)` has vanishing derivative. That question has a long history in the adaptation literature, where the stronger requirement of perfect adaptation demands the derivative vanish identically over an interval and the capable network topologies have been characterised numerically and then structurally (Ma et al. 2009; Ferrell 2016; Araujo and Liotta 2018; Khammash 2021). Cramer's rule applied to the implicitly differentiated equilibrium condition gives `x'_o = ±f_{ι,I} det(H)/det(J)`, where the homeostasis matrix `H` is the Jacobian with its first row and last column deleted. Infinitesimal homeostasis is `det(H) = 0`, and the programme built on that identity applies the Frobenius-König theorem to factor the determinant into irreducible blocks, each associated with a subnetwork motif of structural or appendage type (Golubitsky and Stewart 2017; Reed et al. 2017; Golubitsky and Wang 2020; Huang and Golubitsky 2022; Madeira and Antoneli 2024; Antoneli et al. 2025; Lin et al. 2026). The framework has since been extended to mass-action systems obeying conservation laws, using reduced Jacobian and reduced homeostasis matrices that incorporate the conserved quantities (Jin and Rempala 2026).

The two conditions are complementary degeneracies of the same expression. Theirs is the vanishing of the numerator; ours of the denominator. Their derivation presupposes `det J ≠ 0`, since the equilibrium is assumed hyperbolic, so the case Theorem 1 characterises is excluded from their setting by hypothesis. Geometrically, infinitesimal homeostasis is a horizontal tangent of the equilibrium branch in the input-output plane, where the output stops responding to the input; a fold is a vertical tangent, where the branch turns and the implicit function theorem fails. Neither condition implies the other.

Three degeneracies of the same equilibrium branch therefore arise: `det H = 0` is homeostasis, `det J = 0` with `tr J ≠ 0` is the fold, and `tr J = 0` with `det J > 0` is a Hopf bifurcation. This paper exhibits computed instances of the second and third (Sections 5 and 7) and a structural exclusion of the first. Taking the aggregate as output, the network is two nodes with a single arrow, `H` is 1×1, `det H = G_u`, and the only homeostasis type available is the degree-one Haldane form. It does not occur: `G_u > 0` throughout, because primary nucleation and elongation both rise with `u`, and raising `u` also titrates chaperone and protease away from disaggregation and aggregate clearance, both of which enter `G` with a negative sign. All four contributions carry the same sign. Across 40 kinetic draws no traced point has `G_u < 0`, the smallest value anywhere being 6.67×10⁻¹².

What makes our determinant factor differently is H2, which has no counterpart in the input-output setting. Mass balance supplies the row operation that replaces the influx row by `−∇R`, so `det J` factors into gradients of identified physical fluxes rather than into combinatorial blocks. Whether the combinatorial machinery transfers to `det J = −det[∇R; ∇G; ∇C]` is open; a Frobenius-König factorisation (Brualdi and Ryser 1991) of that determinant would classify fold mechanisms by network topology in the way structural and appendage blocks classify homeostasis mechanisms, and the decomposition of Section 6 is a two-term, physically motivated version of the same question.

The complementary question — when regulation fails rather than when it holds — was raised early in that literature and developed less. Nijhout et al. (2014) ask what happens when a homeostatic mechanism is pushed past its range, and Reed et al. (2018) treat homeostasis at equilibria that are not stable. Neither is a singularity-theoretic account of the failure. Theorem 1 supplies one for the case in which the failure is a fold.

The clearance-limited aggregation literature approaches the same object from the modelling side. Thompson et al. (2021) analyse the stability of a minimal aggregate-removal model; Brennan et al. (2023) treat clearance in a spatial setting; Cotton et al. (2025) generalise to a phase-plane framework across eight aggregation and removal mechanisms, extracting two sufficient conditions for bistability and seeding — self-replicating aggregation with a monomer-dependent rate, and limited-capacity removal. Their input parameter is the monomer concentration, which enters every term of both moment equations, so H1 fails for their system and Theorem 1 does not apply to it as stated. The model classes are adjacent rather than nested. What the present result adds within its own class is a derivation of why the bifurcation structure is preserved under extension, rather than a demonstration that it is.

### 3.4 Orientation and multiplicity

Two properties of the candidate set follow from the theorem and bound what it delivers.

**Orientation.** The sign of `d²R/ds²` distinguishes two kinds of fold. Where it is negative, `R` has a constrained local maximum, `R = j` has solutions below `R(x*)` and none above, and a pair of equilibria is destroyed as the influx rises: a collapse-oriented fold. Where it is positive, the pair is created as the influx rises: a birth-oriented fold, which is the lower turning point of a hysteresis loop rather than a collapse point. On the load grid `d²R/ds² < 0` at all 325 folds. Over the kinetic box, 26 of 2765 are birth-oriented, and all 26 carry a collapse-oriented candidate at strictly higher influx; 7 of 153 ordinary folds carry the same pair. Hysteresis is therefore a feature of the model class rather than of the 26. The operative consequence is that solving `det J = 0` returns candidates without distinguishing them, so orientation must be computed and not inferred from a solver having converged.

**Multiplicity.** Theorem 1.3 returns candidates, not a candidate.

**Corollary 1.** *Solving `{G = 0, det J = 0}` returns the candidate set exactly and without a sweep. Identifying which candidate terminates the accessible low-burden branch requires branch information in 9.1% of the kinetic ensemble (252 of 2765 networks, up to five candidates).*

## 4. Extensions

### 4.1 Arbitrary state dimension

**Corollary 2.** *Let the state be `x = (u, a, c, …)` with `du/dt` the only equation containing `j`, and mass balance giving `du/dt + da/dt = j − R(x)`. Writing the remaining equations as `G` and `C`,*

```
det J = −det [ ∇R ; ∇G ; ∇C ]
```

*which vanishes exactly when `∇R` lies in the span of the remaining gradients — the Lagrange condition for a constrained critical point of `R` on the intersection of the non-influx nullclines.*

The proof is the same row operation on the first row of an `n × n` Jacobian. The practical consequence is that a model in this class can be extended without re-deriving its fold condition.

The converse needs restating above two dimensions, because (G2) does two jobs at once in the plane and neither survives. Write the characteristic polynomial as `p(λ) = λⁿ + … + c₁λ + c₀`, so that `c₀ = ±det J` and `c₁` is, up to sign, the sum of the `(n−1)×(n−1)` principal minors. A zero eigenvalue is `c₀ = 0`; it is algebraically simple exactly when `c₁ ≠ 0`; and partial hyperbolicity is the separate requirement that no other eigenvalue lie on the imaginary axis. For `n = 2` both reduce to `tr J ≠ 0`, since `c₁ = −tr J` and the remaining eigenvalue is `tr J` — which is why the distinction is invisible in the plane and why (G2) may be stated as the trace condition there. In general the converse requires (G1) with the non-influx constraint gradients of full rank, `c₁ ≠ 0`, and (G3).

We verify Corollary 2 on two three-state systems. In the first the chaperone pool is dynamical under σ32-style control, with synthesis rising as free chaperone falls; the identity holds at median relative error 2.5×10⁻¹¹ to 2.9×10⁻¹¹, and the `σ₀ → 0` limit reproduces the two-state field exactly. That limit is not itself a test of the identity: with synthesis off, `dc/dt` vanishes identically, `∇C = 0`, and both sides of the identity are exactly zero, so the reduction is checked against the two-state field rather than against the determinant. In the second system aggregate is split into reactive and sequestered compartments: median 1.5×10⁻¹², maximum 4.7×10⁻¹¹.

The corrected conditions hold at every fold state of both. Over the regulated system `|c₁| ≥ 9.61×10⁻³` at the standard controller and `≥ 0.132` at the stronger one; over the sequestered system `|c₁| ≥ 2.13×10⁻²` across 46 fold states solved from 63 seeded attempts. In all three the remaining spectrum stays off the imaginary axis, the closest approach being `|Re λ| = 1.62×10⁻²`. The `k_seq = 0` and `σ₀ = 0` settings are excluded from these margins rather than reported as passing cases: both make a row of the Jacobian identically zero, so `det J = 0` everywhere and the states a fold solve returns for them are not folds.

### 4.2 Growth dilution

Adding cell division leaves the theorem intact and removes the bound that made the capacity argument work.

Dilution does not affect binding, so the diluted field is `du/dt − μu`, `da/dt − μa`. H1 and H2 both survive, with total removal

```
R_tot = R + μ(u,a)·(u + a)
```

so `det J = −(∇R_tot × ∇G_dil)` for any dilution law: verified at median relative error 1.2×10⁻¹⁰ for constant `μ` and 4.7×10⁻¹⁰ for burden-dependent `μ`, with the `μ → 0` reduction exact. For constant `μ` the diluted Jacobian is `J − μI`, so the fold condition says exactly that `μ` is an eigenvalue of the undiluted Jacobian.

The removal ceiling does not survive. `R_tot` contains `μ(u+a)`, unbounded in burden, so the bound `j_crit ≤ C_enz` fails once cells divide. At sufficiently large constant dilution the fold is destroyed outright: continuing in `μ` at the base parameters gives `j_crit` rising 0.1542 → 0.2456 as `μ` goes 0 → 0.08 while the fold state `a*` diverges 0.265 → 0.750, and by `μ = 0.10` no fold exists. Past that point the low-burden branch no longer terminates and burden rises smoothly with influx. Across 25 parameter draws having a fold at `μ = 0`, 23 lose it under constant dilution, at thresholds spanning 3.3 decades (p10/p50/p90 = 0.0033/0.086/0.328).

Burden-dependent growth restores it. With `μ(u,a) = μ₀/(1 + (u+a)/k_μ)` a fold exists at every rate tested up to `μ₀ = 0.3`. The physiological content is that burden-dependent growth reduction preserves the fold in a regime where constant dilution destroys it, not that division alone is incompatible with a fold.

The coupling itself is well established. Expressing protein a cell does not need reduces its growth rate, and the reduction scales with expression level (Dekel and Alon 2005; Shachrai et al. 2010). Coarse proteome partitioning makes the mechanism explicit: raising any sector's share comes at the expense of the ribosomal sector, and the resulting growth law predicts the burden of heterologous expression quantitatively (Scott et al. 2010, 2014; Scott and Hwa 2023). Calibrating against a measured dosage-response — a 3.2% growth rate reduction when misfolded protein represents less than 0.1% of total cellular protein, measured in yeast (Geiler-Samerotte et al. 2011) and used here as a scaling assumption rather than a bacterial calibration — gives slope 32 per unit misfolded proteome fraction and linear arrest at fraction 0.03125, both bounds rather than point values. Under this calibration a fold survives in 30 of 30 networks with `φ_enz` confined to 0.072–0.147.

Four properties must be distinguished and are not interchangeable: existence of a local fold; boundedness of total burden; genuine runaway; and crossing a viability threshold. The results above concern the first two.

**Corollary 3 (decomposition under division).** *Splitting the critical influx into enzymatic and dilution parts,*

```
j_crit = C_enz · φ_enz /(1 − δ),     φ_enz = R_enz(u*,a*)/C_enz,     δ = R_dil(u*,a*)/j_crit
```

*with both factors dimensionless and in [0,1). The identity is algebraic and holds to 1.6×10⁻¹⁶.*

At the base parameters `φ_enz` stays within 0.125–0.134 as `μ` runs 0 → 0.08, a ±4% band, while `δ` runs 0 → 0.39 and carries essentially all the variation. Division does not change the enzymatic condition for the fold; it multiplies the tolerable influx by `1/(1 − δ)`. Loss of the fold is `δ → 1`, and across the draws above the median `δ` at the threshold is 0.35: the fold is typically lost once division does about a third of the disposal work. Since growth rate is part of disposal, an experiment that lets growth rate float is not holding disposal capacity fixed. The two factors are drawn on one axis in Fig. 2 because the contrast is the content: `φ_enz` is flat to within 7.5% of its own mean across the sweep while `δ` traverses most of its range.

![Figure 2](../figures/fig2.pdf)

**Fig. 2** Under division the enzymatic condition does not move; dilution does. A 33-point continuation at the base parameter set under constant dilution, plotted on one shared linear axis because the claim is a contrast between the two curves. Over `μ` from 0 to 0.080, `φ_enz` stays within 0.1245–0.1343, a full width of 7.5% of its own mean, while `δ` runs 0 to 0.3915 and carries essentially all the variation; `j_crit` rises 0.1542 to 0.2456 and `a*` from 0.265 to 0.750. The band on `φ_enz` is a complete enumeration over the sweep, not an extremum over a subsample. The Corollary 3 identity holds to p99 1.74×10⁻¹⁶ and maximum 1.77×10⁻¹⁶ across the sweep.

### 4.3 Self-damaging machinery

In a cell, chaperones and proteases are translated by the error-prone apparatus, so raising `j` degrades the machinery that clears the damage.

We model this as `C_enz(load) = C_0/(1 + ε·load)` applied to both pools, with `ε` swept over four decades in two modes: `load = j` and `load = u + a`. The second places capacity inside `∇R` and `∇G` themselves.

The identity holds in both modes at machine precision: gradient-normalised residual with floor 2.2×10⁻¹⁴ at `ε = 0`, and worst median 6.4×10⁻¹⁴ (influx mode) and 4.6×10⁻¹³ (burden mode) at `ε = 100`, where capacity has fallen to 16.7% and 1.8% of nominal. No correction to the algebraic form is required, for the reason given in Remark 1: the row operation needs only that `j` be additive in `du/dt` and absent from `da/dt`, and that `du/dt + da/dt = j − R` be exact.

What the coupling removes is the fixed nullcline. Under self-damage `{G = 0}` moves with the load: (H1b) fails and Theorem 1.4 replaces Theorem 1.3, so the candidate solve grows from two equations in `(u,a)` to three in `(u,a,j)`. Corollary 1's saving applies to the frozen-capacity case.

It also removes the automatic transversality of Section 3.2, and puts a checkable scalar in its place. With `F_j = (1 − R_j − G_j, G_j)` and `w = (1, 1+λ)`,

```
w·F_j = 1 − R_j + λ G_j,
```

so transversality is `R_j − λG_j ≠ 1`. This is the exact statement that (G4) ceases to be free precisely when capacity feeds back on the influx. Measured at solved self-damaged fold states — 99 converging from 200 attempts across both ladders — the margin is never small. In burden mode the capacity factor contains no `j`, so `R_j` and `G_j` vanish identically and `w·F_j = 1` exactly at all 46 states, which is (G4) reappearing for the same reason it held in Section 3.2. In influx mode the margin rises with the coupling rather than falling, from 1.000 at `ε = 0.01` to 1.38 at `ε = 100` over 53 states. The feedback moves transversality away from failure in this model; that it could move the other way is what the condition is for.

Self-damage lowers the fold substantially without steepening the approach. Median `j_crit` relative to the frozen fold falls 0.999 / 0.990 / 0.925 / 0.640 / 0.322 across the influx ladder and 0.995 / 0.949 / 0.740 / 0.324 / 0.131 across the burden ladder. The critical-slowing exponent shows no detected shift: paired over the 19 networks carrying both values, median −0.4763 damaged against −0.4813 frozen, with a paired-difference distribution centred near zero (Wilcoxon `p` = 0.312).

**Corollary 4 (square-root capacity ceiling).** *In influx mode every removal flux carries a factor `1/(1 + εj)`, so the necessary condition `j ≤ C_0` of the frozen model becomes `εj² + j − C_0 ≤ 0`, that is*

```
j ≤ ( √(1 + 4εC_0) − 1 )/(2ε)   →   √(C_0/ε)   for large ε.
```

A linear capacity ceiling becomes a square-root one: doubling the machinery buys only √2 in tolerable error rate once the machinery is itself error-prone. The bound is exact and is not binding in the range examined, with `j_crit/j_max` of median 0.039–0.186 and largest observed 0.623 over 8 draws, so it is a necessary condition and not the fold. Continuation with intermediate rungs recovers folds that a direct solve loses at large `ε`, with recovery counts that do not decrease monotonically in `ε` (influx 7/5/6/6/4, burden 6/7/7/4/2); a genuine loss of the fold would be monotone, so this is continuation failure rather than folds disappearing. Where a fold does solve it satisfies `sin(∇R, ∇G) < 2.0×10⁻⁹` throughout that sweep. The consequence is that Corollary 4 is least checkable numerically in the regime where it would bind.

### 4.4 Robustness to the catalytic closure

Section 2.1 stipulates the catalytic step rather than deriving it. Nothing in Section 3 or in Sections 4.1 to 4.3 depends on that choice — the identity, the converse, the `n`-state form and the dilution corollaries follow from H1–H3, which hold for any closure preserving mass balance. The quantitative results of Sections 6 and 7 do depend on it, and this section reports by how much.

The alternative is the mechanistically standard one: catalytic flux proportional to the bound complex, so that `v_ref = k_ref C_f u_f/K_CU` and likewise for the other three removal terms. It removes both stipulations at once. There is no second saturation, and cycles drawing on one pool now compete for it automatically, because the complexes they form appear in the same conservation law. The removal ceiling is unchanged at `C_tot + (ρ_U + ρ_A)D_tot`, since every complex concentration is bounded by its own pool total, so `φ` means the same thing under both and the comparison is well posed. Everything else is held fixed: the same draws, the same seeding from the same recorded states, one field of the parameter set changed.

Two things move, and by similar factors. Median `φ` at the fold rises from 0.0825 to 0.353, so the ceiling overestimates the tolerable influx twelvefold under the published rate laws and 2.8-fold under the alternative; paired over the 1991 networks yielding a fold state under both, the median ratio is 3.81 and 1516 of 1991 lie outside a factor of two. And the incidence of stability loss before the fold falls from 108 of 2766 traced networks to 18 of 2027. Fewer networks admit a solvable fold state at all under the alternative, 2027 against 2767 of the same 2884 draws, and that difference is reported rather than absorbed.

Neither movement is mysterious. The utilisation of nominal capacity at the fold is what `φ` aggregates, and removing a saturating factor from every removal flux raises it. The oscillatory region of Section 7 is characterised by `κ_a`, the Michaelis constant of aggregate clearance, which the alternative closure does not contain; what survives there is the weaker statement that the destabilising mechanism is saturation of clearance, whether that saturation enters through a Michaelis step or through competition for a finite pool.

The consequence for what this paper claims is stated where the claims are made. The ceiling overestimate is a property of the closure as much as of the model class, and the abstract, Section 6 and Section 10 say so. The identity and its corollaries are not.

## 5. Numerical verification

Two ensembles support the results. The **load grid** holds kinetics fixed at the base parameter set and sweeps nascent occupancy and rescue allocation, giving 325 folds; it describes how the fold moves within one network. The **kinetic box** is a Latin-hypercube ensemble of 5000 draws from a stipulated parameter box, of which 2884 admit a fold and 2767 yield a solvable fold state; it describes how the fold varies across networks. Every quantity below is recomputed from tracked generators and asserted by the accompanying test suite.

**Table 2** Numerical verification of the identity and the solver. Each row names its own ensemble because the two are not interchangeable, and each maximum is reported beside a p99 because a maximum grows with the size of the population it is drawn from while a p99 does not.

| quantity | ensemble | median | p99 | max |
|---|---|---|---|---|
| identity residual, gradient-normalised | load grid, 325 | 2.34×10⁻¹⁰ | 9.67×10⁻¹⁰ | 1.29×10⁻⁹ |
| `\|G\|` at fold states | load grid, 325 | 4.95×10⁻¹⁴ | 9.40×10⁻¹⁰ | 1.63×10⁻⁹ |
| direct solver against continuation | load grid, 325 | 3.03×10⁻⁷ | 7.20×10⁻⁷ | 7.56×10⁻⁷ |
| `φ` rebuilt from first principles | kinetic box, 2884 | 1.3×10⁻¹³ | 2.98×10⁻⁹ | 7.25×10⁻⁹ |

Maxima over ensembles of different size are not directly comparable, and neither is comparable to a maximum computed at a third size; the p99 column is given for that reason, and at `n = 325` it rests on a small number of upper-tail observations.

The identity residual sits at the central-difference floor. The parallelism residual is not zero at recorded states and should not be: those are bracketed approximations whose leading eigenvalue is of order 10⁻⁴ rather than 0. The correlation between log parallelism error and log leading eigenvalue is +0.9960 over the load grid, showing the residual to be bracket tolerance rather than failure of the identity, and the tightest-bracketed state in the complete population has `|λ|` = 8.10×10⁻¹⁰ with `sin(angle)` = 7.75×10⁻⁹ (Fig. S1).

One caution attaches to the self-damage check of Section 4.3. Because `du/dt = j − R − G` holds pointwise and the central-difference operator is linear, that check reproduces the row operation exactly at any step size and carries no truncation term; measured slopes in the step size are −0.97 and −0.94 with no minimum over four decades, which is roundoff alone. It certifies that the implementation preserves mass balance. The analytic argument carries the result.

## 6. Where the fold sits

Writing `R = c_f s_ref + ρ_U d_f s_u + ρ_A d_f s_a` against `ceiling = c_tot + ρ_U d_tot + ρ_A d_tot`, and `φ = j_crit/ceiling`, the saturation fractions at the fold answer how far below the ceiling it lies. The decomposition is evaluated at re-solved fold states satisfying `{G = 0, det J = 0}` rather than where continuation stopped; of the 2884 draws admitting a fold, 2767 yield such a state and the remaining 117 are excluded and carry no entry. Evaluated instead at the recorded continuation states, `φ` reads 0.0806 and the ceiling factor thirteenfold rather than twelvefold.

**Table 3** Where the fold sits relative to the capacity ceiling, at re-solved fold states over the 2767 networks that admit one. The saturation fractions are Michaelis factors, so a value of 1 would mean the machinery is running at maximum velocity.

| at the fold, kinetic box, 2767 solved states | median |
|---|---|
| `φ` | 0.0825 |
| refolding saturation `s_ref` | 0.180 |
| soluble degradation saturation `s_u` | 0.159 |
| aggregate clearance saturation `s_a` | 0.049 |
| shortfall recovered by removing saturation | 35.9% |
| shortfall recovered by removing sequestration | 12.6% |

**The fold occurs while the clearance machinery runs at roughly 5–18% of `V_max`**, so the capacity ceiling is a correct bound that overestimates the tolerable influx about twelvefold. That clearance saturation eventually fails to balance aggregation is the qualitative content of the limited-capacity condition identified by Cotton et al. (2025); the decomposition above quantifies how far short of saturation the failure occurs, and splits it.

Two questions must not be merged. Saturation dominates the *magnitude* of `φ`; the *existence* of a turning point requires the aggregation runaway that drives the free pools down at high burden. The counterfactuals answer the first only.

Dispersion is large. A nested design crossing 10 kinetic draws with a 7×7 load grid gives a between-to-within variance ratio for `φ` of 5.9, with per-network spread as large as 13.6× observed when both load coordinates sweep, against an 8.86× between-network spread — so `φ` is network-characteristic but not load-invariant. The saturation box was chosen deliberately wide, so part of the between-network spread is a property of the sampling. Adding σ32-style regulation widens the p5–p95 width of `s_u` from 0.890 to 0.968 over the 30 networks of that experiment rather than narrowing it, with regulated median `s_u` 0.323 against 0.169 unregulated; only 14 of 30 regulated networks converged against 24 unregulated, so that comparison is indicative. Over the kinetic box the p5–p95 width of `s_u` is 0.867.

Fig. 3 shows why the dispersion, rather than the median, is what a reader should take from this section: the p5–p95 spans cover most of the unit interval, so a single measured saturation fraction discriminates weakly. Its inset gives the reason no screening floor is applied to the low-`s_a` tail — the distribution runs smoothly across five decades with no gap, and the median slides continuously from 0.076 to 0.398 with whatever floor is imposed, so any floor is a free parameter moving a load-bearing number fivefold.

The claim this section supports is that the capacity bound overestimates by roughly an order of magnitude and that the fold is a computable constrained critical point — not that the fold occurs at any particular fraction of `V_max`.

![Figure 3](../figures/fig3.pdf)

**Fig. 3** Saturation of the clearance machinery at the fold, over the 2767 kinetic-box networks that yield a solvable fold state, with none excluded. The distributions are evaluated at the same re-solved states as the table in Section 6, so figure and table describe one population; the p5–p95 widths are 0.877 for `s_ref`, 0.867 for `s_u` and 0.892 for `s_a`. The point the figure carries is the width, not the centre: the spans are wide enough that a single measurement discriminates weakly against the capacity-exhaustion alternative. **Inset:** median `s_a` against an imposed lower screening floor, from 0.076 at a floor of 10⁻⁴ to 0.398 at 2×10⁻². The distribution runs smoothly across five decades with no gap, so any floor is a free parameter that moves a load-bearing median by a factor of five. No floor is applied anywhere in this paper.

## 7. Instabilities preceding the fold

The fold is not always the first loss of stability of the low-burden branch. Tracking `tr J` along that branch from small `j` up to `j_crit`, the load grid shows no crossing in any of its 325 networks, with `max(tr J)` of median −0.338 and largest −0.243. Over the kinetic box, 104 of 2766 networks — 3.8% — reach `tr J = 0` before the fold, with `det J > 0` at the crossing in all of them (minimum 1.5×10⁻⁶), so the eigenvalues there are `±i√det J`. The first crossing occurs at a median of 0.83 of the way to `j_crit`, with deciles 0.40 / 0.61 / 0.83 / 0.95 / 0.995.

A crossing network and a non-crossing one are traced together in Fig. 4a: the first rises through zero at about 0.82 of the way to `j_crit` and stays positive to the fold, while the second approaches the fold from below without ever reaching it.

Integration confirms the classification independently of the branch computation. Perturbing each of the 104 crossing networks at its own `tr J` maximum, 104 grow by more than tenfold and 102 escape, with median `max|δ|/δ₀` of 9851; 205 non-crossing networks perturbed at their own `tr J` maximum give 0 escapes and median 1.00, with no overlap between the distributions (Fig. 4b). All 309 states are exact equilibria. The integration also recovers the eigenvalue: growth exponent within 5% of `max Re(λ)` in 90 of 104, and oscillation period within 5% of `π/ω` in 47 of 49.

The crossings occupy a coherent region of parameter space (Fig. 4c). Relative to non-crossing networks, aggregate clearance is faster (`ρ_A` 7.5× higher) and saturates more sharply (`κ_a` 11.7× lower), and aggregate formation is growth-dominated rather than nucleation-dominated (`α_g/α_n` shifting from 3.7 to 33.6, with `α_g` alone only 1.67× higher). That configuration is a recognised destabilising one. Saturable rather than linear degradation of a reaction product drives a stable steady state unstable as the Michaelis constant of the degrading step falls (Goldbeter 2013); the same mechanism has been identified as the determinant of oscillation in NF-κB signalling, where strong saturation of IκB degradation destabilises the fixed point and weakening it restores stability (Krishna et al. 2006), and Michaelis-Menten degradation kinetics promotes oscillation across positive-feedback, negative-feedback and combined architectures (Xu and Qu 2012). Here the autocatalytic partner is elongation.

This is a distinct mechanism from the oscillation question in the regulated heat-shock response, which turns on transcriptional feedback and delay (El-Samad et al. 2005; Zheng et al. 2016; Sivéry et al. 2016). The networks above carry no regulation.

The incidence figure is a property of the stipulated parameter box and is not offered as a prediction. The load grid's cleanliness establishes that the base kinetic parameter set lies outside the oscillatory region, not that the region is small: the load grid holds kinetics fixed, and the crossing networks are distinguished by kinetic parameters.

![Figure 4](../figures/fig4.pdf)

**Fig. 4** The fold is not always the first loss of stability. **(a)** `tr J` along the low-burden branch against `j/j_crit`, on a symmetric-log axis, for one crossing network and one non-crossing network. Both exemplars are chosen by rule rather than by inspection: the crossing case is the network nearest the median crossing position among those whose branch is resolved below 0.3 `j_crit`, and the non-crossing case is the network nearest the population median of `max(tr J)`. **(b)** Amplification of a perturbation applied at each network's own `tr J` maximum, integrating the full nonlinear system. The 104 crossing networks all grow by more than tenfold and 102 leave the linear neighbourhood, at a median amplification of 9851; the 205 non-crossing controls give 0 escapes and a median of 1.00. The distributions do not overlap. This panel shares no failure mode with the branch computation in (a). **(c)** Where the crossing networks sit in the sampled `(κ_a, ρ_A)` plane. Relative to the rest they have `ρ_A` 7.5× higher and `κ_a` 11.7× lower: clearance that is fast and sharply saturating.

## 8. Aggregation asymmetry in dividing *E. coli*

Unstressed *E. coli* accumulate aggregates at the old pole and lose reproductive capacity relative to new-pole progeny (Stewart et al. 2005; Lindner et al. 2008). Two features of the model class bear on this.

### 8.1 What the model class cannot do

No variant of the model examined places a stable attractor at low burden. Under constant dilution the system is bistable, with a 12.6- to 32.3-fold aggregate ratio between branches, but that regime predicts a growth-rate loss of exactly zero because `μ` cannot respond to burden. Under the calibrated hyperbolic law the system is monostable in four of six settings tested. Under linear arrest there is no bounded high-burden state: `μ` reaches exactly zero beyond `k_μ`, switching dilution off and removing the bound. Adding a sequestered compartment — aggregate drained into an inert pool that neither nucleates nor occupies chaperone or protease — produces 12 bistable cells against a control's 3, but none of 384 qualified cells predicts a loss within the measured range, every bistable cell predicting 0.482, 0.508, 0.954 or 1.000. Inverting the growth law, the high state is at least 7.5 to 254 times more aggregate-laden than the measured cell, and the upper figure is a lower bound because the 1.000 values are the linear-arrest law clamping at states 3.4 to 43 times past arrest. A real old-pole cell carries a visible inclusion body and divides at 96–99% of normal.

### 8.2 The observation does not require bistability

The old-pole cell does not occupy a second basin. It inherits a physical inclusion body at each division, which is a continuously renewed perturbation rather than an attractor, and a monostable system under a renewed perturbation has a stationary offset with no separatrix.

The two-state diluted model under the calibrated hyperbolic law, with asymmetric partitioning of aggregate at division and no additional state variable, reproduces the measured lineage difference. Of 728 parameter settings, 709 settle, and 43 score within the measured range with the partitioning active. The control is exact: all 66 settings at symmetric partitioning give a lineage difference of 0.0 with standard deviation 0.0, which follows algebraically, since symmetric partitioning into half the volume leaves concentration unchanged and is what the `−μx` term already encodes.

### 8.3 A requirement on the aggregate load

The scored quantity is Lindner et al.'s `Δ(GR_old − GR_new)/GR_mean`, whose aggregate-attributable part is 1.04–1.78%, obtained from the reported lineage difference of −3.95 ± 0.5% and the reported aggregate-attributable share of 30–40%.

Partitioning is set by one physical quantity. The visible focus is indivisible and passes entirely to one daughter, but the model's `a` is total aggregate, of which the focus is a share `β`. The old-pole daughter then receives `(1+β)a` and the new-pole daughter `(1−β)a` after division into half the volume, which is the scalar rule at `f_eff = (1+β)/2`. With `k_μ = 0.03125/p_qc`, the scored quantity is `32 × (B_old − B_new)` as a proteome fraction exactly, so the quality-control proteome share and the growth rate enter only through the stationary load. Inverting gives a requirement scaling as `1/β`:

**Table 4** The old-pole aggregate load required by the account of Section 8.3, as a function of the share `β` of aggregate held in the visible polar focus, over the same range Fig. 5 plots. The final column compares the requirement with the only aggregate load that has been measured, 5–10% of total protein in Δ*rpoH* at 30 °C; the columns pair the extremes, so they give the widest and narrowest separations consistent with both intervals.

| `β` | `f_eff` | required old-pole aggregate (% of proteome) | ratio to the Δ*rpoH* load |
|---|---|---|---|
| `1.00` | 1.000 | 0.047 – 0.080 | 62× – 214× |
| `0.75` | 0.875 | 0.061 – 0.104 | 48× – 165× |
| `0.50` | 0.750 | 0.091 – 0.157 | 32× – 110× |
| `0.25` | 0.625 | 0.182 – 0.314 | 16× – 55× |
| `0.05` | 0.525 | 0.913 – 1.570 | 3.2× – 11× |

The damping factor for the daughter's relaxation during its own generation is 0.346–0.355, with weak `β` dependence.

The table spans the same `β` range the figure plots. At `β = 1.00`, where the whole aggregate is in the focus, the requirement sits 62× to 214× above the Δ*rpoH* load; at `β = 0.05` it is 3.2× to 11× above it. Reporting only the rows down to `β = 0.25` would put the nearest approach at 16× and overstate the separation fivefold.

No source we could find bounds `β` under unstressed growth. Winkler et al. (2010) report aggregate mass and focus size separately, and under heat stress, which is neither the condition at issue nor a focus-share measurement. The requirement is therefore a `β`-indexed family (Fig. 5), and the direction is the informative part: the less of the aggregate held in the focus, the more aggregate the account requires, and the closer the requirement moves to the only load anyone has measured. Tomoyasu et al. (2001) report 5–10% of total protein aggregated in Δ*rpoH* mutants at 30 °C and no detectable aggregation in *rpoH*⁺ cells at the same temperature; wild type is a bound rather than a value, and the reported mutant figures are not a limit of detection for the assay.

![Figure 5](../figures/fig5.pdf)

**Fig. 5** The old-pole aggregate load required by the account of Section 8.3, as a function of the focus share `β`, over the same range as the table there. The requirement scales as `1/β`: at `β = 1.00` it is 0.0467–0.0803% of the proteome and at `β = 0.05` it is 0.9126–1.5696%. The interval at each `β` comes from the aggregate-attributable part of the measured lineage difference, 1.04–1.78%. The shaded band is the only aggregate load anyone has measured, 5–10% of total protein in Δ*rpoH* at 30 °C; no measurement bounds `β`, so no vertical line is drawn.

Two measurements close the question: the wild-type aggregate fraction under unstressed growth, and the share of that aggregate held in the polar focus, by quantitative fluorescence of the focus against total aggregated protein by sedimentation in the same cells.

## 9. Predictions

**Mathematical.** The identity `det J = −(∇R × ∇G)`, its `n`-state form, the converse under (G1)–(G4), the `μ → 0` and `σ₀ → 0` reductions, the `J − μI` form under constant dilution, the `(φ_enz, δ)` decomposition, and Corollary 4 follow from H1–H3 and no observation bears on them.

**The fold sits far below saturation.** The alternative is capacity exhaustion, in which the fold occurs as `s → 1`. The observable is a ratio, so testing it does not require absolute copies-per-cell calibration. Two cautions: the dispersion in `s_u` is wide (p5–p95 width 0.867 over the kinetic box), so a single measurement discriminates weakly, and adding regulation widens rather than narrows it.

**The fold moves under nascent load.** Since nascent chains consume chaperone capacity without contributing influx, raising the load of perfectly folding protein should lower the tolerable damage influx; a model handling the two loads independently predicts no shift. Direction was correct in 67 of 68 settings over a 100-fold ladder, with median shift 1.22×. Growth rate sets the load coordinate and contributes about 1.6× on its own against the ~1.8× contrast expected, so the experiment requires externally fixed growth rate: chemostat or turbidostat, not batch culture.

**Oscillation before collapse.** In networks with fast, sharply saturating aggregate clearance against growth-dominated aggregation, damped oscillations should appear in aggregate burden as the damage influx rises, becoming sustained before any runaway. The predicted period is `π/ω` with `ω = √det J` at the crossing, and the observable signature is that the amplitude grows while the mean burden remains low.

**A ruler, not a test.** Approaching the fold, relaxation time scales as `τ ~ (j_crit − j)^(−1/2)` (fitted slope 0.5077, `r²` = 1.0000). This is generic to folds, so confirming it selects this model over nothing in particular. Its use is as an instrument for locating the fold in an experiment that then asks whether the fold moves.

**The measurement the theory most wants.** Whether growth arrest under misfolding burden is complete or asymptotic. Under the hyperbolic form, dilution keeps bounding the high-burden state and loss of proteostasis is a reversible switch; under linear arrest, `μ` reaches exactly zero and the high state is a runaway. The available dosage-response constrains the slope only below 0.1% misfolded protein, about 1.5 decades below the arrest burden where the two forms diverge, so it cannot adjudicate. That single functional form decides whether losing proteostasis in a dividing cell is recoverable.

## 10. Discussion

The fold of a mass-balanced clearance network is not something one has to find by sweeping. Given the removal flux and the aggregate nullcline, it is where their gradients become parallel, and the critical influx is the removal flux evaluated there. The statement requires no fitted parameter, holds for any number of states, survives growth dilution, and survives capacity that degrades with the load it carries. Two things bound what it delivers: the candidate set can have more than one element, and the fold is not always the first loss of stability.

What it does not do is predict a number without measured parameters. The dispersion in `φ` across the kinetic ensemble is large enough, and wide enough after calibration and after adding regulation, that the quantitative predictions above are ratio comparisons rather than point predictions. This is a structural result rather than a predictive one.

Santra et al. (2019) reach a tipping point in a proteostasis model by chaperone titration under cumulative damage, which is the same qualitative claim about how collapse arises. The present contribution is the exact flux-geometric condition for where it occurs, and the classification of what else can happen first.

The Hopf result identifies a regime worth looking for. Oscillatory dynamics of aggregate burden in a clearance network requires no regulatory circuit — only fast, sharply saturating clearance against autocatalytic growth, which is the classical saturable-degradation destabilisation acting on an aggregation system. Whether real quality-control networks occupy that region is open, and Section 9 states the signature.

The prediction a reader can act on is Section 8.3. If the aggregate burden of an aging *E. coli* old-pole cell falls outside the interval required at the measured focal share, the account of the aging asymmetry given here is wrong. Neither of the two measurements needed has been attempted as far as we can determine.

---

## Data and code availability

All analysis code, parameter configurations, and the test suite asserting every numerical quantity reported here are at https://github.com/khatvangi/proteostasis-law-theory under the MIT licence, archived at https://doi.org/10.5281/zenodo.21794565 (concept DOI, resolving to the latest version).

## Supplementary material

![Figure S1](../figures/figS1.pdf)

**Fig. S1** Parallelism residual against leading eigenvalue, over all 325 load-grid folds with no subsample. The correlation is +0.9960 in log-log with slope 1.000, so `sin(angle)` is proportional to bracket looseness rather than merely increasing with it; the tightest-bracketed state has `|λ|` = 8.10×10⁻¹⁰ and `sin(angle)` = 7.75×10⁻⁹. The figure also carries the normalisation contrast: normalising by `max(|det J|, |cross|)` is `0/0` at an exact fold, and that residual degrades as the bracket tightens (median 3.133×10⁻⁷ on the tighter half against 1.455×10⁻⁷ on the looser, 2.15× worse where it matters most, log-log correlation −0.835, maximum 1.54×10⁻²), whereas the gradient normalisation used throughout correlates at +0.041 with a maximum of 1.29×10⁻⁹.

![Figure S2](../figures/figS2.pdf)

**Fig. S2** Strategy space cut by the fold condition: 469 feasible strategies, 13 non-dominated, with `j/j_crit` running 0.2271 to 0.9652 along the front (median 0.5010). The exact throughput optimum sits at `j/j_crit` = 1.000000 exactly, on the feasibility boundary; the best grid strategy reaches 0.897488, the gap being discretisation of the grid rather than an interior optimum. Supplementary Section S1 develops a strategy-space consequence of the fold condition: cutting a two-coordinate trade-off surface by the constraint `j(α,R) < j_crit(R)` places the throughput optimum exactly on the feasibility boundary, while along the non-dominated front `j/j_crit` runs 0.227 to 0.965.

---

## References

Antoneli F, Golubitsky M, Jin J, Stewart I (2025) Homeostasis in input-output networks: structure, classification and applications. *Math Biosci*. arXiv:2405.03861

Araujo RP, Liotta LA (2018) The topological requirements for robust perfect adaptation in networks of any size. *Nat Commun* 9:1757.

Brennan GS, Thompson TB, Oliveri H, Rognes ME, Goriely A (2023) The role of clearance in neurodegenerative diseases. *SIAM J Appl Math* S172–S198. doi:10.1137/22M1487801

Brualdi RA, Ryser HJ (1991) *Combinatorial matrix theory*. Cambridge University Press, Cambridge.

Cohen SIA, Linse S, Luheshi LM, Hellstrand E, White DA, Rajah L, Otzen DE, Vendruscolo M, Dobson CM, Knowles TPJ (2013) Proliferation of amyloid-β42 aggregates occurs through a secondary nucleation mechanism. *Proc Natl Acad Sci USA* 110:9758–9763.

Cotton MW, Goriely A, Klenerman D, Meisl G (2025) A universal phase-plane model for in vivo protein aggregation. arXiv:2511.18893

Dekel E, Alon U (2005) Optimality and evolutionary tuning of the expression level of a protein. *Nature* 436:588–592.

El-Samad H, Kurata H, Doyle JC, Gross CA, Khammash M (2005) Surviving heat shock: control strategies for robustness and performance. *Proc Natl Acad Sci USA* 102:2736–2741.

Ferrell JE (2016) Perfect and near-perfect adaptation in cell signaling. *Cell Syst* 2:62–67.

Ferrone F (1999) Analysis of protein aggregation kinetics. *Methods Enzymol* 309:256–274.

Geiler-Samerotte KA, Dion MF, Budnik BA, Wang SM, Hartl DL, Drummond DA (2011) Misfolded proteins impose a dosage-dependent fitness cost and trigger a cytosolic unfolded protein response in yeast. *Proc Natl Acad Sci USA* 108:680–685.

Goldbeter A (2013) Oscillatory enzyme reactions and Michaelis-Menten kinetics. *FEBS Lett* 587:2778–2784. doi:10.1016/j.febslet.2013.07.031

Golubitsky M, Stewart I (2017) Homeostasis, singularities and networks. *J Math Biol* 74:387–407.

Golubitsky M, Wang Y (2020) Infinitesimal homeostasis in three-node input-output networks. *J Math Biol* 80:1163–1185.

Hipp MS, Kasturi P, Hartl FU (2019) The proteostasis network and its decline in ageing. *Nat Rev Mol Cell Biol* 20:421–435.

Huang Z, Golubitsky M (2022) Classification of infinitesimal homeostasis in four-node input-output networks. *J Math Biol* 84:25.

Jin J, Rempala GA (2026) Infinitesimal homeostasis in mass-action systems. *J Math Biol* 92(3):35. doi:10.1007/s00285-026-02352-y

Khammash M (2021) Perfect adaptation in biology. *Cell Syst* 12:509–521.

Knowles TPJ, Waudby CA, Devlin GL, Cohen SIA, Aguzzi A, Vendruscolo M, Terentjev EM, Welland ME, Dobson CM (2009) An analytical solution to the kinetics of breakable filament assembly. *Science* 326:1533–1537.

Krishna S, Jensen MH, Sneppen K (2006) Minimal model of spiky oscillations in NF-κB signaling. *Proc Natl Acad Sci USA* 103:10840–10845. doi:10.1073/pnas.0604085103

Lin X, Antoneli F, Wang Y (2026) Automated classification of homeostasis structure in input-output networks. *Bull Math Biol*, article 113.

Lindner AB, Madden R, Demarez A, Stewart EJ, Taddei F (2008) Asymmetric segregation of protein aggregates is associated with cellular aging and rejuvenation. *Proc Natl Acad Sci USA* 105:3076–3081.

Ma W, Trusina A, El-Samad H, Lim WA, Tang C (2009) Defining network topologies that can achieve biochemical adaptation. *Cell* 138:760–773.

Madeira JLO, Antoneli F (2024) Homeostasis in networks with multiple inputs. *J Math Biol* 89:17.

Meisl G, Michaels TCT, Linse S, Knowles TPJ (2017) Scaling behaviour and rate-determining steps in filamentous self-assembly. *Chem Sci* 8:7087–7097.

Nijhout HF, Best J, Reed MC (2014) Escape from homeostasis. *Math Biosci* 257:104–110.

Reed M, Best J, Golubitsky M, Stewart I, Nijhout HF (2017) Analysis of homeostatic mechanisms in biochemical networks. *Bull Math Biol* 79:2534–2557.

Reed MC, Duncan W, Nijhout HF, Best J, Golubitsky M (2018) Homeostasis despite instability. *Math Biosci* 300:130–137.

Santra M, Dill KA, de Graff AMR (2019) Proteostasis collapse is a driver of cell aging and death. *Proc Natl Acad Sci USA* 116:22173–22178.

Scott M, Gunderson CW, Mateescu EM, Zhang Z, Hwa T (2010) Interdependence of cell growth and gene expression: origins and consequences. *Science* 330:1099–1102.

Scott M, Klumpp S, Mateescu EM, Hwa T (2014) Emergence of robust growth laws from optimal regulation of ribosome synthesis. *Mol Syst Biol* 10:747.

Scott M, Hwa T (2023) Shaping bacterial gene expression by physiological and proteome allocation constraints. *Nat Rev Microbiol* 21:327–342.

Shachrai I, Zaslaver A, Alon U, Dekel E (2010) Cost of unneeded proteins in *E. coli* is reduced after several generations in exponential growth. *Mol Cell* 38:758–767.

Sivéry A, Courtade E, Thommen Q (2016) A minimal titration model of the mammalian dynamical heat shock response. *Phys Biol* 13:066008.

Stewart EJ, Madden R, Paul G, Taddei F (2005) Aging and death in an organism that reproduces by morphologically symmetric division. *PLoS Biol* 3:e45.

Thompson TB, Meisl G, Knowles TPJ, Goriely A (2021) The role of clearance mechanisms in the kinetics of pathological protein aggregation involved in neurodegenerative diseases. *J Chem Phys* 154:125101.

Tomoyasu T, Mogk A, Langen H, Goloubinoff P, Bukau B (2001) Genetic dissection of the roles of chaperones and proteases in protein folding and degradation in the *Escherichia coli* cytosol. *Mol Microbiol* 40:397–413.

Wang Y, Huang Z, Antoneli F, Golubitsky M (2021) The structure of infinitesimal homeostasis in input-output networks. *J Math Biol* 82:62.

Winkler J, Seybert A, König L, Pruggnaller S, Haselmann U, Sourjik V, Weiss M, Frangakis AS, Mogk A, Bukau B (2010) Quantitative and spatio-temporal features of protein aggregation in *Escherichia coli* and consequences on protein quality control and cellular ageing. *EMBO J* 29:910–923.

Xu L, Qu Z (2012) Roles of protein ubiquitination and degradation kinetics in biological oscillations. *PLoS ONE* 7:e34616. doi:10.1371/journal.pone.0034616

Zheng X, Krakowiak J, Patel N, Beyzavi A, Ezike J, Khalil AS, Pincus D (2016) Dynamic control of Hsf1 during heat shock by a chaperone switch and phosphorylation. *eLife* 5:e18638.
