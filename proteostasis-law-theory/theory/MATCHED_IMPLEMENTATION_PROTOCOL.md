# Matched implementation protocol — Phase 2A

This document is the specification the Phase 2A code is written against. It is
cited by section number from the code:

- `scripts/phase2/mapping.py` cites **section 3** for the term-by-term identity;
- `scripts/phase2/t0_equivalence.py` cites **section 3** for the O(epsilon)
  claim and **section 5** for its tolerances.

Nothing here is a scientific conclusion. This document defines what is being
compared and under what conditions a comparison is legitimate.

---

## 1. The problem this protocol exists to solve

Two independent implementations of the two-state proteostasis model exist:

| | state variables | resource closure |
|---|---|---|
| **boron** — `scripts/proteostasis/model.py` | TOTAL burdens `(u, a)`, free plus machinery-bound | implicit rapid-equilibrium solve |
| **nitrogen** — `.../proteostasis-law-theory-nitrogen-check/src/proteostasis_model.py` | FREE burdens `(u, a)` | closed form |

They disagree on reported percentages. Before that disagreement can be
attributed to anything scientific, it must be established that the two
implementations are being asked the same question. Two separate things can
produce a difference and they are confounded in any naive comparison:

1. **model form** — whether substrate sequestration is represented at all;
2. **root protocol** — the solver, its initial guesses, its residual tolerance,
   and its attractor confirmation.

The matched benchmark is a 2x2 factorial that separates them.

---

## 2. Nitrogen is the epsilon = 0 limit

**Nitrogen's model is not a different model. It is the epsilon = 0 face of the
boron family.**

Boron carries total burdens and divides out two substrate sequestration factors:

```
uf = u / su,   su = 1 + cf/kappa_cu + df/kappa_du
af = a / sa,   sa = 1 + cf/kappa_ca + df/kappa_da
```

Nitrogen asserts `su = sa = 1` identically and keeps only the resource-side
occupancy. Setting every substrate-binding affinity to zero — that is, sending
every `kappa` to infinity — turns boron's equations into nitrogen's, term for
term. This is an algebraic identity, not an approximation or a limit that has to
be argued.

The sequestration ladder makes that limit a scannable single scalar. Every
binding constant is tied to one parameter `epsilon`:

```
kappa_cu = kappa_ca = c_tot / epsilon        c_u = c_a = epsilon / c_tot
kappa_du = kappa_da = d_tot / epsilon        d_u = d_a = epsilon / d_tot
```

so that

```
epsilon = max(c_tot/kappa_cu, d_tot/kappa_du, c_tot/kappa_ca, d_tot/kappa_da)
```

holds exactly, and **both sides use the same resource-side occupancy at every
rung**. Under the ladder `kappa_ca = kappa_cu` and `kappa_da = kappa_du`, so the
two sequestration factors collapse to one scalar:

```
su = sa = 1 + epsilon * sigma,      sigma = cf/c_tot + df/d_tot  in (0, 2]
```

The ladder used by the benchmark is

```
EPSILON_LADDER = (1e-6, 1e-3, 1e-2, 1e-1, 0.3, 1.0, 2.0)
```

`epsilon = 1e-6` is the identity anchor. `epsilon = 0` itself is **refused** as a
boron `Params` — `bindingConstants` raises, because boron's `Params.validate`
requires strictly positive finite kappas. The exact zero is reached by using the
free-limit model form (`phase2/nitrogen_limit.py`), not by passing zero into
boron. That is a deliberate code-level distinction: the anchor is approached in
the boron arm and attained in the free arm.

Boron's own Phase 1 baseline (`kappa_cu = 0.3`, `c_tot = 0.6`) sits at
`epsilon = 2.0` on the chaperone axis and `epsilon = 0.8` on the degradation
axis. The ladder therefore **brackets** the Phase 1 baseline rather than
reproducing it.

### Gauge

Boron's Phase 1 nondimensionalization fixes the concentration scale as
`phi = c_tot + d_tot`, so `c_tot + d_tot == 1` there. The matched benchmark
instead adopts the **nitrogen gauge**: `mu = 1` and `d_tot = D_TOT_GAUGE = 1`.
This is a rescaling of the concentration unit, not a change of dynamics, but it
means **the boron parameter vectors produced here are not in Phase 1's
phi-normalized gauge and must not be compared to Phase 1 parameter tables
coordinate by coordinate.**

---

## 3. The mapping, term by term

`mapping.nitrogenToBoron` in the `mu = 1` gauge. Each line below is an algebraic
identity between the two flux expressions, obtained by writing boron's flux with
`uf = u`, `af = a` (the epsilon = 0 face) and matching coefficients against
nitrogen's.

Boron fluxes (`scripts/proteostasis/model.py::fluxes`) versus nitrogen fluxes
(`scripts/phase2/nitrogen_limit.py::fluxes`):

| flux | boron | nitrogen | forced identity |
|---|---|---|---|
| influx | `j` | `j_u + j_fold + j_other` | `j = q.j / mu` |
| refold | `cf * uf/(kappa_ref + uf)` | `ref_capacity * cf_n * u/(ref_k + u)` | `c_tot = ref_capacity * c_tot_n / mu`; `kappa_ref = ref_k` |
| soluble degradation | `rho_U * df * uf/(kappa_u + uf)` | `deg_u_capacity * df_n * u/(deg_u_k + u)` | `rho_U = deg_u_capacity * d_tot_n / (mu * d_tot)`; `kappa_u = deg_u_k` |
| nucleation | `alpha_n * uf**m` | `nucleation * u**m` | `alpha_n = nucleation / mu`; `m = m` |
| growth | `alpha_g * uf * af` | `growth * u * a` | `alpha_g = growth / mu` |
| disaggregation | `alpha_d * cf * af/(kappa_dis + af)` | `disaggregation * cf_n * a/(disaggregation_k + a)` | `alpha_d = disaggregation * c_tot_n / (mu * c_tot)`; `kappa_dis = disaggregation_k` |
| aggregate clearance | `rho_A * df * af/(kappa_a + af)` | `deg_a_capacity * df_n * a/(deg_a_k + a)` | `rho_A = deg_a_capacity * d_tot_n / (mu * d_tot)`; `kappa_a = deg_a_k` |

Resource denominators:

| | boron | nitrogen | forced identity |
|---|---|---|---|
| chaperone | `cf = c_tot/(1 + nu + uf/kappa_cu + af/kappa_ca)` | `cf_n = c_tot_n/(1 + n_load + o_load + c_u*u + c_a*a)` | `nu = n_load + o_load`; `c_u = 1/kappa_cu`; `c_a = 1/kappa_ca` |
| degradation | `df = d_tot/(1 + uf/kappa_du + af/kappa_da)` | `df_n = d_tot_n/(1 + d_u*u + d_a*a)` | `d_u = 1/kappa_du`; `d_a = 1/kappa_da` |

`bindingConstants` produces `c_u = epsilon/c_tot` and `kappa_cu = c_tot/epsilon`,
which are exact reciprocals by construction. This is asserted rather than
assumed by `tests/phase2/test_mapping.py::testBindingConstantsAreExactReciprocals`,
and the round trip `boron -> nitrogen -> boron` is asserted to be the identity
by `testBoronToNitrogenToBoronIsIdentity`.

`nu` contributes **no damage influx**. It acts purely by consuming chaperone
capacity. `tests/phase2/test_nitrogen_transcription.py::testNascentOccupancyIsNotDamageInflux`
pins this, because reading `n_load` as an influx term would be a silent
falsification of the whole comparison.

### Why the residual discrepancy is exactly first order in epsilon

At finite epsilon the identity above is an approximation, and its error order is
a falsifiable prediction. From `uf = u/(1 + epsilon*sigma) = u*(1 - epsilon*sigma)
+ O(epsilon**2)`:

```
F_free - F_boron = epsilon * sigma * [ u dF/duf + a dF/daf ] + O(epsilon**2)
```

with a coefficient that vanishes nowhere generic. Therefore:

- an observed exponent of **1** confirms the mapping;
- an observed exponent of **2** would mean the leading term cancels — the mapping
  has aligned the wrong pair of models;
- an observed exponent of **0** would mean some term is mismatched independently
  of epsilon, i.e. a straight transcription error.

Only exponent 1 is consistent with this section.

---

## 4. The 2x2 factorial

```
        model form         x        root protocol
        ----------                  -------------
  boron (finite sequestration)      P_BORON     dense 9x9 log multistart,
  free  (epsilon = 0 face)                      hybr in log coordinates,
                                                RELATIVE residual tolerance,
                                                +/-10% kick, 5% relative return,
                                                Radau, t_end = 5e4

                                    P_NITROGEN  four fixed linear guesses,
                                                ABSOLUTE residual tol 1e-7,
                                                +/-0.1% kick, ABSOLUTE 1e-3
                                                return, DOP853, t_end = 150
```

The four cells named in the audit:

| cell | model form | epsilon | protocol | what it isolates |
|---|---|---|---|---|
| 1 | boron | 1e-6 | P_NITROGEN | 1 vs 2 = code equivalence |
| 2 | free | 1e-6 | P_NITROGEN | |
| 3 | boron | 1e-6 | P_BORON | 1 vs 3 = solver effect alone |
| 4 | boron | 2.0 | P_BORON | 3 -> 4 = model form alone |

Every cell is scored on the **same** deterministic 7-D Latin hypercube
(seed `20260731`, n = 2000) with the remaining sixteen coordinates pinned, so a
cell-to-cell difference can only come from the two factors. The sample matrix
SHA-256 is pinned in `phase2/lhs.py::SAMPLE_HASHES` and was verified
bit-identical on boron (scipy 1.14.0 / numpy 2.2.6) and nitrogen (scipy 1.15.2 /
numpy 1.26.4). If the two hosts ever disagreed on the sample matrix, every
cross-host cell comparison would be comparing different parameter draws.

Both protocols are written against a `models.ModelAdapter` and touch the model
only through that interface, so the *same protocol code* runs against both model
forms. If the two arms went through different protocol code the factorial would
no longer separate the factors.

### Admissibility is reported, not folded into the label

Labels are `stable_attractor`, `no_bounded_attractor_operational`, and
`numerical_failure`. They deliberately do **not** encode a threshold split,
because there are four defensible admissibility criteria —
(burden `u + a < H`) versus (componentwise `u < H` and `a < H`), each in TOTAL
and in FREE coordinates — and mixing them silently is exactly the confound this
benchmark exists to remove. All four are computed and written to every row.
`H = 1.0` is shared: boron's `burden_threshold` and nitrogen's
`threshold_u`/`threshold_a` are all 1.0.

`no_bounded_attractor_operational` is an **operational** label. It means no
attractor was found within this finite numerical protocol. **It is not a
mathematical nonexistence result.**

### The continuation, and why it does not break cell 1

Nitrogen's root protocol runs `hybr` in linear coordinates, whose line search
freely probes states with `a < 0`. Nitrogen's model is rational and evaluates
there; boron's `solveFreePools` raises, because the conservation solve is stated
for nonnegative total burdens. Left alone this fails on ~100% of samples and
deletes cell 1.

Three options were considered. Clamping every evaluation at `max(x, 0)` was
**measured and rejected**: applied symmetrically to both arms it changes 247–265
of 400 free-limit results at every epsilon tested, so it is a large uncontrolled
protocol change, not a no-op. The adopted fix is analytic continuation of the
boron field past the boundary (`phase2/boron_continuation.py`), licensed by the
fact that the guard is a domain assertion rather than a singularity.

The discipline that keeps this honest: `solveFreePoolsExtended` **delegates** to
the shipped `model.solveFreePools` whenever `u >= 0` and `a >= 0`. Every initial
guess is positive and every accepted root is nonnegative, so the root, its
residual, its Jacobian and its eigenvalues all come from the shipped boron code
path; the continuation is reached only by intermediate line-search iterates.
Every row records `continuation_evals`, so **a cell can never claim to be a pure
boron code path when it was not.**

---

## 5. T0 and its tolerances

T0 (`scripts/phase2/t0_equivalence.py`) is the mandatory first test of Phase 2A
and the cheapest possible falsification of section 3: no sweep, no root finder,
no integrator. **If T0 fails, every downstream label comparison is void and the
mapping must be re-derived before a single percentage is compared.**

It asserts four things on the 12-state grid taken verbatim from
`configs/phase1/experiment_a.json` (T0 anchors on the shipped config rather than
a table typed into the script, so it cannot silently test a stale parameter set):

1. at `epsilon = 1e-6` the two vector fields agree to relative discrepancy
   below `RHS_ANCHOR_TOL`;
2. the same for the two analytic Jacobians, below `JAC_ANCHOR_TOL`;
3. each model's analytic Jacobian matches its own central-difference Jacobian,
   so a passing comparison cannot be two identical mistakes;
4. the discrepancy scales as O(epsilon) across {1e-6, 1e-3, 1e-2, 1e-1}.

### Tolerances derived a priori

These are derived from the structure, not tuned to observed numbers.

| constant | value | derivation |
|---|---|---|
| `RHS_ANCHOR_TOL` | `1e-5` | Section 3 predicts relative discrepancy `~ epsilon * sigma` with `sigma in (0, 2]`, so at `epsilon = 1e-6` the worst case is `~2e-6`. `1e-5` is the task requirement and sits a factor ~5 above that bound — loose enough not to fail on the predicted first-order term, tight enough that any zeroth-order mismatch fails outright. |
| `JAC_ANCHOR_TOL` | `1e-5` | Same argument. Differentiating the O(epsilon) expansion preserves the order. |
| `SELF_JACOBIAN_TOL` | `1e-6` | Central differencing with `h_rel = 1e-6` has truncation error O(h**2) ~ 1e-12 relative and roundoff O(eps_machine/h) ~ 1e-10, both matrix-scaled. `1e-6` leaves several orders of headroom while still failing on any sign or factor error. |
| `SLOPE_TARGET` | `1.0` | Section 3. |
| `SLOPE_TOL_FULL` | `0.05` | Least squares over all four epsilons. The largest rung (`1e-1`) carries a visible O(epsilon**2) contamination, which biases the four-point fit slightly above 1; 0.05 accommodates that curvature while excluding slope 0 or 2 by a wide margin. |
| `SLOPE_TOL_ASYMPTOTIC` | `0.01` | Two smallest epsilons only, where the O(epsilon**2) term is negligible; the asymptotic slope must therefore be much closer to 1 than the full fit. |
| `R2_MIN` | `0.999` | A power law in log-log is a straight line. Anything materially below this is not a clean power law and the exponent would not be interpretable. |

### Reduction rules

- **RHS**: componentwise `|delta| / max(|reference|, floor)` reduced by max, with
  `floor = j` (the influx scale). Componentwise rather than norm-wise because a
  norm can hide one badly mismatched term behind a large sibling; the floor stops
  a near-zero component from manufacturing a divergent ratio where the field
  genuinely vanishes.
- **Cross-model Jacobian**: componentwise, with the floor set to `1e-6` of the
  largest entry of the reference matrix, so an absolute-zero entry cannot
  manufacture an infinite ratio.
- **Self Jacobian** (analytic vs central difference of the *same* model):
  matrix-scaled, `max|delta| / max|reference|`. Componentwise is not merely
  inconvenient here, it is undefined: at state `(0, 0)` the analytic entry
  `d(da/dt)/du` is exactly zero (`alpha_n * u**m` has zero slope at the origin
  for `m > 1`) while the one-sided boundary stencil returns `alpha_n * h`. No
  relative statement about an exact zero exists. The cost is that an error below
  `tolerance x largest entry` in a small entry is undetected — far below any
  error that could change an eigenvalue sign.

### T0 result on boron, 2026-07-31

Run independently and reproduced during this submission
(`results/phase2/matched_20260731T175912-0500/evidence/t0_equivalence.txt`),
exit code 0, all 10 checks PASS:

| epsilon | `rhs_rel_max` | `jac_rel_max` |
|---|---|---|
| 1e-06 | 4.800103e-06 | 8.263311e-06 |
| 1e-03 | 4.801702e-03 | 8.268664e-03 |
| 1e-02 | 4.816307e-02 | 8.314829e-02 |
| 1e-01 | 4.975588e-01 | 8.623627e-01 |

RHS slope (log-log, 4 pt) = 1.002332, R^2 = 0.99999323, asymptotic = 1.000048.
Jacobian slope = 1.002841, R^2 = 0.99999132, asymptotic = 1.000094.

The exponent is 1, as section 3 requires. **T0 passing licenses the label
comparison; it is not itself a scientific result.**

---

## 6. (!) What may and may not be compared

**Percentages are not comparable outside the matched benchmark.**

Any percentage — fraction viable, fraction multistable, fraction admissible —
produced by a boron run and a nitrogen run is comparable **only** when all of the
following hold:

1. both sides scored the **same sample matrix**: same seed, same n, same pinned
   sixteen coordinates, and a matching `sample_matrix_sha256` in both manifests;
2. both sides used the **same root protocol** (`P_BORON` or `P_NITROGEN`), not
   one each;
3. both sides applied the **same admissibility criterion**, named explicitly out
   of the four reported — burden vs componentwise, total vs free coordinates;
4. the epsilon rung is stated, and boron-arm-vs-free-arm comparisons are made at
   the **same** rung;
5. **every matched cell has completed.** A percentage from a partial run is not
   a percentage.

In particular:

- The pre-existing **50,000-sample independent nitrogen sweep is NOT the matched
  counterpart** and must not be used as the matched result. It uses a different
  sample matrix, so it cannot be differenced against any boron cell — the two
  would be scoring different parameter draws. It remains valid only as a
  free-limit result in its own right.
- Phase 1 percentages are in the phi-normalized gauge and on different parameter
  ranges. They are not comparable to Phase 2A cells.
- A difference between a boron cell and a free cell at `epsilon = 2.0` is a
  **model-form** statement. A difference between the same model form under two
  protocols is a **solver** statement. Reporting either as the other is the
  specific error this benchmark was built to prevent.

The benchmark produces labelled rows. It produces no scientific conclusion.

---

## 7. Determinism and the two hashes

Every cell summary in `manifest.json` carries two hashes, and they mean
different things:

- **`tsv_sha256`** — a byte record of the file as written. It is deliberately
  **not** reproducible: every row carries a per-sample wall-clock `seconds`.
- **`payload_sha256`** — SHA-256 over the result columns only, in `ROW_FIELDS`
  order, excluding the columns named in `payload_excludes` (currently exactly
  `seconds`). **This is the hash that may be compared across reruns, across
  worker counts, and across hosts.**

The scientific columns are bit-identical across repeat runs and across
`workers = 1` versus `workers = 4`;
`tests/phase2/test_benchmark_design.py::TestCellOutputIsReproducible` asserts
that `seconds` is the *only* non-deterministic column, so that if a future change
makes another column wall-clock dependent, the excluded set must be updated
deliberately rather than silently.

---

## 8. Stop-work conditions

Work stops and no percentage is reported if any of these occur:

- T0 fails any of its 10 checks — the mapping is wrong, and every downstream
  label is void;
- the sample-matrix SHA-256 differs between hosts — the two are scoring
  different draws;
- a cell reports `numerical_failure` rows traceable to a coding defect rather
  than to a genuine numerical limit of the protocol;
- `rows_using_continuation` is non-zero in a cell claimed to be a pure boron
  code path under `P_BORON` (it is expected to be non-zero only for the boron
  arm under `P_NITROGEN`);
- any matched cell is incomplete.
