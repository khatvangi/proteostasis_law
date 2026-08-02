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

Informally: **the fold is the constrained maximum of total removal on the
aggregate nullcline.**

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

## Limits

- **One model.** Two states, no regulation, no compartmentation, and — most
  seriously — **no growth dilution**. In dividing cells, dilution outpaces
  proteolysis for most proteins, so `d_tot` may be the wrong object for
  "disposal". This would move `phi` and is not a caveat that can be waved away.
- **The maximum is verified, not proved unique.** If `R` has several constrained
  critical points on the branch, "the" fold is branch-dependent. The theorem
  identifies critical points; ordering them is open.
- **The kappa ranges were chosen adversarially wide.** Saturation dominates
  `phi`, and the kappas set saturation, so the 8.86x between-network spread is
  partly a property of the sampling box rather than of biology.
- **Some draws collapse at `s_a` near 0.003**, where aggregation is so fast the
  low-burden branch barely exists. These should be screened rather than averaged
  into a median.

## Reproduction

```
python scripts/phase3/fold_theorem.py          # prints every number above
python -m pytest tests/phase3 -q               # asserts them
```

`results/` is gitignored; without the Phase 1 run root the module prints an
explicit `SKIP` and exits 0, and the artefact tests skip. The model-level tests
run on a clean checkout and are the ones that pin the theorem itself.
