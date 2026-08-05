# A3: what moves when the catalytic closure changes

Deliverable of task A3, required because task A1 landed as option (b) —
`theory/RATE_LAWS.md` establishes that the catalytic layer of the model is a
stipulated closure and measures where its two assumptions bind.

Generator: `scripts/analysis/closure_a3.py` (`--hopf` for the second stage).
Reduction: `data/computed/closure_a3.json`.

## The comparison

| | michaelis (published) | complex (this test) |
|---|---|---|
| refolding | `C_f u_f/(K_ref + u_f)` | `C_f u_f/K_CU = [CU]` |
| soluble degradation | `rho_U D_f u_f/(K_U + u_f)` | `rho_U D_f u_f/K_DU` |
| disaggregation | `alpha_d C_f a_f/(K_dis + a_f)` | `alpha_d C_f a_f/K_CA` |
| aggregate clearance | `rho_A D_f a_f/(K_A + a_f)` | `rho_A D_f a_f/K_DA` |

Held fixed: the same 2884 draws admitting a fold, the same seeding of each fold
solve from the same recorded continuation state, the same solver. The two
pipelines differ in exactly one field of `Params`.

Why the comparison is well posed: the removal ceiling
`c_tot + (rho_U + rho_A) d_tot` bounds removal under BOTH closures, because every
complex concentration is bounded by its own pool total —
`C_f u_f/K_CU <= C_tot` follows from the chaperone conservation law. So
`phi = j_crit/ceiling` measures the same thing in both.

## Results

**A note on which quantity is compared.** `s_ref`, `s_u`, `s_a` are Michaelis
factors and have no counterpart under the complex closure. What both closures
possess is the UTILISATION of nominal capacity, `v/V_max`, which is what "running
at x% of `V_max`" should mean and is what `phi` aggregates. Under the michaelis
closure the two differ by the free-pool fraction: `util_ref = (c_f/c_tot) s_ref`,
which is about a fifth of `s_ref` at the median. Both are given.

| at the fold, kinetic box | michaelis | complex |
|---|---|---|
| networks yielding a solvable fold state | 2767 of 2884 | 2027 of 2884 |
| median `phi` | 0.0825 | 0.3527 |
| ceiling overestimate `1/phi` | **12.1x** | **2.8x** |
| refolding utilisation `v/V_max` | 0.0363 | 0.2562 |
| soluble degradation utilisation | 0.0584 | 0.3298 |
| aggregate clearance utilisation | 0.0190 | 0.0658 |
| free chaperone fraction `c_f/c_tot` | 0.297 | 0.209 |
| Michaelis factor `s_ref` | 0.1803 | — |
| Michaelis factor `s_u` | 0.1590 | — |
| Michaelis factor `s_a` | 0.0495 | — |

Paired over the 1991 networks solved under both closures:

| | value |
|---|---|
| median `phi`, michaelis | 0.0776 |
| median `phi`, complex | 0.3597 |
| median ratio complex/michaelis | **3.81** |
| p5 – p95 of the ratio | 0.220 – 62.3 |
| outside a factor of two | 1516 of 1991 |

### Crossing incidence (stage 2)

Branch re-traced under each closure at that closure's own fold states.

| | michaelis | complex |
|---|---|---|
| fold states | 2767 | 2027 |
| traced | 2766 | 2027 |
| trace failures | 1 | 0 |
| lose stability before the fold | 108 | 18 |
| of those, on branches whose j-maximum is the fold | 104 | 17 |
| incidence | **3.90%** | **0.89%** |
| median `max(tr J)` on the branch | −0.317 | −0.705 |
| median first crossing, as a fraction of `j_crit` | 0.831 | 0.931 |

## Verdict

The work order's threshold was a factor of two in median `phi`. It is exceeded by
the median alone: 0.0825 against 0.3527, a factor of 4.3, and 3.81 paired.

**The twelvefold ceiling factor is closure-dependent and must say so at the
claim** — abstract, Section 6 headline, Section 10. Section 4.4 of the manuscript
carries the comparison.

The Hopf incidence moves in the same direction and by a similar factor. That is
not a coincidence and not a puzzle: the oscillatory region of Section 7 is
characterised by `kappa_a`, the Michaelis constant of aggregate clearance, and
the alternative closure has no such parameter. What survives the closure change
is the qualitative mechanism — saturation of clearance against autocatalytic
growth — rather than the incidence rate, which was already stated to be a
property of the stipulated box rather than a prediction.

**What does NOT move:** Theorem 1.1 through 1.4, Theorem 2 and its corollary,
Corollary 2's `n`-state identity, Corollaries 3 and 4, and the dilution results.
All follow from H1–H3, which hold for either closure, since mass balance is exact
by construction in both — the transfer terms `n`, `g` and `v_dis` appear in
`du/dt` and `da/dt` with opposite signs and cancel identically regardless of how
the removal fluxes are written.
