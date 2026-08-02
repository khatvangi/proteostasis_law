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

`a*` diverges as the threshold is approached: the fold point runs off to infinite
aggregate burden and escapes. Past it there is no collapse boundary at all,
because a cell dividing at a burden-independent rate can always dilute its way
out.

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

## Limits

- **One model.** Two states, no regulation, no compartmentation. Growth dilution
  was the most serious of these and is now addressed in Consequence 5 — the
  theorem survives it, the removal ceiling does not — but the dilution law there
  is a one-parameter guess, not a measured growth-burden relation, and the `mu`
  threshold was mapped only at the base parameters, not across the box.
- **`phi` itself is now ill-defined under dilution.** Its denominator was the
  enzymatic removal ceiling, which Consequence 5 shows is no longer an upper
  bound on `j_crit`. Every `phi` in Consequences 2 and 3 is therefore a
  zero-dilution quantity. A dividing-cell analogue needs a new denominator, and
  none is proposed here.
- **The critical point is verified, not proved unique.** If `R` has several
  constrained critical points on `{G = 0}`, "the" fold is branch-dependent. The
  theorem identifies critical points; ordering them, and establishing whether the
  relevant one is a maximum over the whole curve, is open — see the precision
  note above.
- **The kappa ranges were chosen adversarially wide.** Saturation dominates
  `phi`, and the kappas set saturation, so the 8.86x between-network spread is
  partly a property of the sampling box rather than of biology.
- **Some draws collapse at `s_a` near 0.003**, where aggregation is so fast the
  low-burden branch barely exists. These should be screened rather than averaged
  into a median.

## Reproduction

```
python scripts/phase3/fold_theorem.py    # the theorem, phi, the nested design
python scripts/phase3/dilution.py        # growth dilution and fold survival
python -m pytest tests/phase3 -q         # 29 checks asserting both
```

`results/` is gitignored; without the Phase 1 run root `fold_theorem.py` prints
an explicit `SKIP` and exits 0, and its artefact tests skip. The 22 model-level
tests — including every dilution test, which needs no artefacts — run on a clean
checkout and are the ones that pin the mathematics.
