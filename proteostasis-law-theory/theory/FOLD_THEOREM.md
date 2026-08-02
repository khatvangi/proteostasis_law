# The Fold Theorem

Phase 3. This file states and proves where the collapse boundary is, which
Phase 1 and Phase 2 located only numerically. Every number below is recomputed
by `scripts/phase3/fold_theorem.py` from the Phase 1 run root and asserted by
`tests/phase3/test_fold_theorem.py`.

## Statement

Let `G(u,a) = da/dt` and let `R(u,a)` be the total removal flux

```
R = v_ref + v_degU + v_degA .
```

**Theorem.** For the conserved-resource model of `theory/DYNAMICAL_SYSTEM.md`:

1. equilibria are exactly the points of the curve `{G = 0}` at which `R = j`,
   and that curve does not depend on `j`;
2. `det J = -( grad R x grad G )` identically;
3. therefore a saddle-node occurs exactly at a **constrained critical point of
   `R` restricted to `{G = 0}`**, and

```
        j_crit = R(u*, a*)                      (exact)
        phi    = j_crit / removalCeiling
```

Informally: **the collapse boundary is where total removal stops responding to
burden along the aggregate nullcline.**

**Precision note, added after the dilution work below.** The proven statement is
*critical point*, not *maximum*; an earlier version of this line said
"constrained maximum" and that gloss is **withdrawn**. Scanning `R` along the
branch reached by taking the first root in `a` shows it rising monotonically to
the end of that branch in every case tested, and the solved fold state has
`dG/da > 0` (9.1e-03 to 3.9e-02, i.e. not at the nullcline's turning point). The
critical point therefore lies *past the turn*, not on the first-root branch.
`lowerNullclineA` in the phase 3 module is a bracketing heuristic only and does
not identify which branch the fold lives on. That the critical point is a
maximum over the whole curve is **not established**; the theorem and every
number below are unaffected, since all of them use `det J = 0` directly.

## Proof

Two structural facts do all the work.

*The influx enters one equation only.* `j` appears in `du/dt` and nowhere in
`da/dt`. The aggregate nullcline is therefore a fixed curve in burden space;
changing the load slides the state along it but does not reshape it. This gives
(1) once combined with the next step.

*Mass balance.* The internal `u <-> a` transfer (nucleation + growth −
disaggregation) cancels between the two equations, so

```
du/dt + da/dt = j - R
```

exactly. This is the identity already tested as A5 (residual 3.75e-16).

Now apply the determinant-preserving row operation `row1 -> row1 + row2` to the
Jacobian of the equilibrium system:

```
det J = det [ d(du/dt)/du   d(du/dt)/da ]  =  det [ -R_u  -R_a ]
            [ d(da/dt)/du   d(da/dt)/da ]         [  G_u   G_a ]

      = -( R_u G_a - R_a G_u ) = -( grad R x grad G ) .
```

That is (2). A saddle-node requires `det J = 0`, i.e. `grad R` parallel to
`grad G`, which is precisely the Lagrange condition for a critical point of `R`
subject to `G = 0`. With (1), the critical value is `j_crit = R(u*,a*)`. ∎

## Verification

Against the Phase 1 run root, at the recorded fold states:

| quantity | value |
|---|---|
| `det J` vs `-(grad R x grad G)`, median relative error | **1.436e-07** |
| correlation of `log sin(angle)` with `log \|eigenvalue\|` | **+0.9987** |
| `\|G\|` at recorded fold states, max | 8.2e-10 |
| solver `{G=0, det J=0}` vs the continuation sweep, max relative error | **6.652e-07** |
| `phi` rebuilt from first principles, 2884 folds, median / max error | **1.3e-13 / 7.3e-09** |

The identity residual is at the central-difference floor. The parallelism
residual is *not* zero at the recorded states, and should not be: those states
are bracketed approximations whose leading eigenvalue is about −2e-4, not 0. The
+0.9987 correlation between the parallelism error and the recorded eigenvalue
shows the residual is bracket tolerance rather than a failure of the identity —
the one state bracketed to `eig = -4.2e-9` has `sin(angle) = 3.8e-8`.

## Consequence 1 — the fold is a 2x2 root solve

Neither `G = 0` nor `det J = 0` contains `j`, so fold location needs no
continuation sweep. `foldSolve` in the phase 3 module solves the pair directly
and reproduces the Phase 1 sweep to 6.7e-07. This is what makes the nested
design below affordable.

## Consequence 2 — what actually sets phi

Writing `R = cf.s_ref + rho_U.df.s_u + rho_A.df.s_a` against
`ceiling = c_tot + rho_U.d_tot + rho_A.d_tot`, `phi` falls below 1 for exactly
two reasons, separated by counterfactual over all 2884 Phase 1 folds:

| at the collapse point | median |
|---|---|
| `phi` observed | **0.0769** |
| refolding saturation `s_ref` | **0.175** |
| soluble degradation saturation `s_u` | **0.155** |
| aggregate clearance saturation `s_a` | **0.056** |
| share of the shortfall recovered by removing **saturation** | **35.8 %** |
| share of the shortfall recovered by removing **sequestration** | **12.6 %** |

**Collapse happens while the clearance machinery runs at roughly 6–18 % of
V_max.** It is not capacity exhaustion. Superlinear nucleation (`alpha_n.u_f^m`,
`m > 1`) overtakes sublinear saturating removal long before removal maxes out.
This is the mechanism behind prediction P2 in `PHASE2_CLOSURE_FINAL.md` §7, and
it makes `removalCeiling` a correct bound that is about 13x too loose to predict
anything on its own.

Note the two questions are distinct and must not be merged: **saturation
dominates the magnitude of `phi`**, whereas the existence of a turning point
requires the sequestration/aggregation runaway that drives the free pools down
at high burden. The counterfactuals above answer the first question only.

## Consequence 3 — phi is network-characteristic, but NOT load-invariant

The Phase 1 hypercube varied loads and kinetics together and so cannot separate
them. `nestedInvariance` crosses them explicitly — 10 kinetic draws x a 7x7
`(nu, chi)` load grid, 446 folds solved:

| varying | spread in `phi` |
|---|---|
| `chi` with `nu` held fixed | **1.80x** |
| `nu` with `chi` held fixed | **1.59x** |
| both jointly, worst network | up to 13.6x |
| the kinetic parameters `theta` | **8.86x** |

Variance ratio between/within: **5.9**.

So `phi` is more a property of the network than of its load, but the earlier
reading of it as a load-invariant material constant is **withdrawn**: holding one
load coordinate gives ~1.6–1.8x, but sweeping both compounds to as much as 13.6x
in some networks, which overlaps the between-network spread. Load-sensitivity is
itself network-dependent — two of ten draws show essentially none (1.01x, 1.03x
across `nu`), two show 3.7–4.0x.

## Consequence 4 — the empirical prediction this sharpens

`PHASE2_CLOSURE_FINAL.md` §7 P3 asserts that jointly lethal dose pairs exist,
which almost any nonlinear system satisfies. The theorem replaces it with a
contrast against the standard alternative:

> **Capacity-exhaustion model:** collapse when clearance flux saturates, `s -> 1`.
> **This law:** collapse at `s_a` of order 0.05–0.4, far below saturation, at a
> network-characteristic fraction `phi` of total capacity.

The observable is a **saturation fraction**, i.e. a ratio, so testing it does not
require absolute copies-per-cell calibration — which would otherwise be the
hardest quantity to obtain. Secondary contrast: at **fixed growth rate**,
shifting the chaperone/protease split should move the boundary ~1.8x, whereas
changing the network should move it ~8.9x.

**Design constraint, and it is load-bearing:** growth rate sets `nu`, and `nu`
drift contributes its own ~1.6x which compounds with `chi`. An experiment that
does not hold growth rate fixed degrades the contrast toward the between-network
spread and cannot discriminate.

Claim label: **Empirical hypothesis, untested.** No organism data was used
anywhere in Phase 1, 2 or 3.

## Consequence 5 — growth dilution, and why the boundary needs feedback

The model has no cell division, and in growing cells dilution outpaces
proteolysis for most proteins. `scripts/phase3/dilution.py` adds it without
touching the frozen upstream model: binding is unaffected by dilution, so the
diluted field is `du/dt - mu.u`, `da/dt - mu.a`.

**The theorem survives, for any dilution law.** `j` still enters `du/dt` only, and
the internal transfer still cancels, giving `du/dt + da/dt = j - R_tot` with

```
R_tot = R + mu(u,a).(u + a)          total removal, now including dilution
```

so `det J = -( grad R_tot x grad G_dil )` again. Verified at median relative
error **1.2e-10** (constant `mu`) and **4.7e-10** (burden-dependent `mu`), with
the analytic diluted Jacobian matching central differences to **6.7e-11** and the
`mu -> 0` reduction to the frozen model exact at **0.0**.

For constant `mu` the diluted Jacobian is `J - mu.I`, so the saddle-node
condition says exactly that **`mu` is an eigenvalue of the undiluted Jacobian**.

**The removal ceiling does not survive.** `R_tot` contains `mu.(u+a)`, which is
unbounded in burden, so the analytic bound A8 is false once cells divide. The
consequence is not cosmetic — with *constant* dilution the fold is destroyed
outright above a threshold. Continuation in `mu` at the base parameters, each
solve seeded from the previous and every accepted point satisfying
`|G| < 1e-17` and `|det J| < 1e-17`:

| `mu` | 0 | 0.02 | 0.04 | 0.06 | 0.08 | 0.10 |
|---|---|---|---|---|---|---|
| `j_crit` | 0.1542 | 0.1738 | 0.1950 | 0.2186 | 0.2456 | **none** |
| `a*` | 0.265 | 0.318 | 0.389 | 0.496 | **0.750** | — |

`a*` diverges as the threshold is approached: the fold point runs off to large
aggregate burden and escapes. Past it the low-burden branch no longer terminates,
so there is no discontinuous transition — burden rises smoothly with influx
instead. A cell dividing at a burden-independent rate can dilute its way out.

**Correction to an earlier phrasing.** This was first written as "there is no
collapse boundary at all". That overstates it in two ways, both established in
Consequence 7: below the threshold, losing the low-burden branch is a *jump to a
bounded high-aggregate state*, not divergence; and above it, what disappears is
the discontinuity, not viability.

**Growth feedback restores it.** With `mu(u,a) = mu0/(1 + (u+a)/k_mu)` and
`k_mu = 0.5`, a fold exists at every rate tested up to `mu0 = 0.3`, with `j_crit`
rising from 0.1542 to 0.2689.

> **The collapse boundary of a dividing cell exists because burden slows growth.**
> Dilution alone is an unbounded disposal channel; only its throttling by the
> burden it disposes of makes viability finite.

This is the first place the theory requires a coupling it had not previously
stated, and it changes what must be measured: a perturbation experiment that lets
growth rate float is not holding the disposal capacity fixed, because growth rate
*is* part of disposal.

Claim labels: the identity, the `mu -> 0` reduction and the `J - mu.I` form are
**Mathematical**; every threshold and `j_crit` above is **Computational**, at the
base parameters only — the `mu` threshold was not mapped across the parameter box.

## Consequence 6 — a margin that survives division

The margin in Consequence 2 divides by the enzymatic capacity bound, which
Consequence 5 shows is not an upper bound once cells divide. Splitting the
critical influx into its enzymatic and dilution parts repairs this exactly:

```
j_crit = C_enz . phi_enz / (1 - delta)

    phi_enz = R_enzymatic(u*,a*) / C_enz     enzymatic capacity in use at collapse
    delta   = R_dilution(u*,a*) / j_crit     share of disposal done by division
```

Both are dimensionless and lie in [0,1). The identity is algebra and holds to
**0 – 1.6e-16**. At the base parameters:

| `mu` | 0 | 0.02 | 0.04 | 0.06 | 0.08 |
|---|---|---|---|---|---|
| `j_crit` | 0.1542 | 0.1738 | 0.1950 | 0.2186 | 0.2456 |
| `phi_enz` | 0.1285 | 0.1321 | 0.1341 | 0.1335 | 0.1245 |
| `delta` | 0 | 0.0876 | 0.1751 | 0.2671 | 0.3915 |

**`phi_enz` is nearly invariant to dilution** (0.125–0.134, a ±4 % band) while
`delta` carries essentially all the variation. So division does not change the
enzymatic condition for collapse; it multiplies the tolerable influx by
`1/(1 - delta)`. The escape in Consequence 5 is `delta` approaching 1, which
makes it quantitative rather than merely observed.

Across 25 parameter draws that have a boundary at `mu = 0`, **23 lose it** under
constant dilution. The threshold spans **3.3 decades** (p10/p50/p90 = 0.0033 /
0.086 / 0.328), and `delta` at the threshold has median **0.35** (p10 0.195,
p90 0.847) — the boundary is typically lost once division is doing roughly a
third of the disposal work.

## Consequence 7 — division makes the system bistable, and uniqueness fails

Enumerating the whole nullcline rather than one branch settles the uniqueness
question left open above, with a split verdict.

**Without dilution the critical point is unique**: the curve closes (152 lower
and 152 upper points) and carries exactly one sign change of `det J`.

**With dilution it is not.** At `mu = 0.04` a second, distinct constrained
critical point exists at `u = 0.1585, a = 1.9835`, satisfying `|G| = 2.8e-17` and
`|det J| = 1.6e-17` with eigenvalues `(-1.083, 0)` — a genuine saddle-node — and
its critical influx is **0.1585, below** the 0.1950 of the first.

The two are the folds of an S-shaped curve, and the reason is structural: with
dilution, `G -> -infinity` at large aggregate whenever `alpha_g.u_f < mu`, so the
high-burden state is *bounded* rather than divergent. Sweeping influx up from
zero burden and back down:

| | |
|---|---|
| up-sweep: low branch survives to | `j = 0.194`, jumps at 0.196 |
| attained state after the jump | `u = 0.079, a = 3.94` — **bounded** |
| down-sweep: high branch survives to | `j = 0.160`, drops back at 0.158 |
| bistable window | **[0.160, 0.194]** |
| computed critical points | 0.158496 and 0.195047 — the window lies inside |

So under division, losing viability is **not** divergence. It is a transition to a
persistent, bounded, aggregate-laden state, and recovery requires lowering the
damage influx **below** the level that triggered the transition — hysteresis.

This does not reinstate anything from `notes/REJECTED_COMPONENTS.md`. Item 7 there
rejects bistability attributed to the old one-variable ODE, and Phase 2 §2.1
demoted it in the *non-dividing* model on the evidence that it appeared in 1.14 %
of draws. Bistability here arises from a model feature those analyses did not
contain, and it appears at the base parameters rather than in a rare corner.
Whether it is generic across the box has **not** been established.

## Limits

- **One model.** Two states, no regulation, no compartmentation. Growth dilution
  was the most serious of these and is addressed in Consequences 5–7; the `mu`
  threshold is now mapped across 25 draws rather than one point.
- **No growth-burden relation has been measured.** Two functional *forms* are
  used — hyperbolic (`dilution.Growth`) and the linear shape implied by proteome
  partitioning (`boundary_structure.LinearGrowth`) — and the qualitative
  conclusions agree under both. Neither is fitted to data, and `k_mu` remains
  free. Quantities that depend on the form, notably `j_crit` at a given `mu0`,
  differ between them (0.269 vs 0.327 at `mu0 = 0.3`) and should not be quoted
  without stating which law produced them.
- **The old margin is a zero-dilution quantity**, since its denominator stops
  bounding `j_crit` once cells divide. Consequence 6 replaces it with
  `(phi_enz, delta)`, which is exact and dimensionless under division; every
  `phi` in Consequences 2 and 3 should be read as `phi_enz` at `delta = 0`.
- **Bistability was found at one parameter point, not surveyed.** Consequence 7
  does not establish how often division makes the system bistable across the
  box, and that is the obvious next question.
- **Uniqueness holds without division and fails with it** (Consequence 7).
  Ordering is settled at the base parameters by the hysteresis sweep — the upper
  fold bounds the healthy branch, the lower fold bounds recovery — but not in
  general, and not across the parameter box. Whether the relevant critical point
  is a maximum over the whole curve remains open (precision note above).
- **The kappa ranges were chosen adversarially wide.** Saturation dominates
  `phi`, and the kappas set saturation, so the 8.86x between-network spread is
  partly a property of the sampling box rather than of biology.
- **Some draws collapse at `s_a` near 0.003**, where aggregation is so fast the
  low-burden branch barely exists. These should be screened rather than averaged
  into a median.

## Reproduction

```
python scripts/phase3/fold_theorem.py      # the theorem, the margin, the nested design
python scripts/phase3/dilution.py          # growth dilution and boundary survival
python scripts/phase3/boundary_structure.py  # uniqueness, bistability, thresholds
python -m pytest tests/phase3 -q           # 42 checks asserting all three
```

`results/` is gitignored; without the Phase 1 run root `fold_theorem.py` prints
an explicit `SKIP` and exits 0, and its artefact tests skip. The 22 model-level
tests — including every dilution test, which needs no artefacts — run on a clean
checkout and are the ones that pin the mathematics.
