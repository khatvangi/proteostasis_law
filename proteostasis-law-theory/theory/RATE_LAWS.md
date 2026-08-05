# Rate laws: what is derived and what is stipulated

Deliverable of task A1. This file states, term by term, which parts of the model
in `scripts/proteostasis/model.py` follow from an explicit reaction network and
which are phenomenological closures. The distinction is not cosmetic: Theorem 1
needs only H1–H3, which every version satisfies, but the quantitative
conclusions of Sections 6 and 7 are computed from these particular rate laws and
inherit whatever they assume.

**Summary.** The binding layer is exact — an identity, given rapid equilibrium
and 1:1 complexes. The catalytic layer is a stipulated closure. The two
approximations it makes are measured below and neither is small in the ensemble
the paper reports. Task A3 therefore applies: the quantitative results are
re-run under a mechanistically standard alternative, and any that move are
reported as closure-dependent at the point of claim.

---

## Part 1. The binding layer is exact

### The network

Two conserved machinery pools and three substrate species:

| symbol | species |
|---|---|
| `C` | chaperone pool, total `C_tot` |
| `D` | degradation machinery pool, total `D_tot` |
| `U` | soluble damaged monomer, total `u` |
| `A` | aggregate burden in monomer-equivalent units, total `a` |
| `N` | nascent chain, free concentration `N_f` |

Five 1:1 complexes, each in rapid equilibrium with its own dissociation
constant:

```
C + U <-> CU     K_CU        C + A <-> CA     K_CA        C + N <-> CN     K_N
D + U <-> DU     K_DU        D + A <-> DA     K_DA
```

Rapid equilibrium means `[CU] = C_f U_f / K_CU` and likewise for the rest, where
a subscript `f` denotes a free concentration.

### The closure

The four conservation laws are then identities, not approximations:

```
u     = U_f + [CU] + [DU]                 = U_f (1 + C_f/K_CU + D_f/K_DU)
a     = A_f + [CA] + [DA]                 = A_f (1 + C_f/K_CA + D_f/K_DA)
C_tot = C_f + [CU] + [CA] + [CN]          = C_f (1 + nu + U_f/K_CU + A_f/K_CA)
D_tot = D_f + [DU] + [DA]                 = D_f (1 + U_f/K_DU + A_f/K_DA)
```

with `nu = N_f/K_N` the nascent-chain occupancy, a parameter rather than a state
because nascent chains are supplied by translation and consume chaperone without
contributing damage influx.

These four equations are exactly the closure solved in
`model.py:_bindingResidual` and `model.py:solveFreePools`. Nothing is
linearised, truncated or fitted. Existence, uniqueness and smoothness of the
solution are proved in `theory/LEMMA0_BINDING.md` (task A2), which also removes
smoothness of `R` and `G` from the list of hypotheses and derives it.

Two consequences worth naming, because the manuscript uses both:

* `u` and `a` are **total** pools — free plus machinery-bound. Every rate law
  below is written in free concentrations, and total substrate is never
  substituted into a formula that expects a free one (decision D004).
* The pools are shared. Raising `u` lowers both `C_f` and `D_f`, which raises
  `A_f` at fixed `a`. Section 3.3's sign argument turns on this and is treated
  in task B6.

---

## Part 2. The catalytic layer is stipulated

### What the code computes

```
v_ref  = C_f * U_f/(K_ref  + U_f)        refolding to native
v_degU = rho_U * D_f * U_f/(K_U   + U_f) soluble degradation
v_dis  = alpha_d * C_f * A_f/(K_dis + A_f) disaggregation
v_degA = rho_A * D_f * A_f/(K_A   + A_f) aggregate clearance
```

Each flux is a product of an available-enzyme factor and a Michaelis factor in
that enzyme's free substrate. The manuscript through v5 printed these in **total**
concentrations, which is a different model; §2.1 is rewritten to print the above.

### Why it is not the standard scheme

For `C + U <-> CU -> C + native` the catalytic flux is `k_cat [CU]`, which by the
rapid-equilibrium closure is `k_cat C_f U_f / K_CU` — proportional to free
enzyme and to free substrate, with no saturation beyond that already carried by
the closure. The coded form saturates a **second** time, in the same substrate,
against a **second** constant. `K_CU` and `K_ref` are therefore being asserted to
describe successive steps: nonproductive binding, then a rate-limiting catalytic
capture. That is a claim about mechanism.

The scheme under which the coded form would be exact is:

1. binding to the pool (constants `K_CU`, `K_CA`, `K_N`, `K_DU`, `K_DA`) is
   **nonproductive** sequestration, counted in all four conservation laws;
2. the enzyme not held in those complexes, `C_f` or `D_f`, runs a separate
   Michaelis–Menten cycle on free substrate with constant `K_ref` (or `K_U`,
   `K_dis`, `K_A`), the productive complex being a sub-population of `C_f`;
3. the cycles sharing one pool — refolding and disaggregation on `C`, soluble
   and aggregate degradation on `D` — occupy **independent** sites, so their
   Michaelis factors do not compete for one denominator;
4. the substrate held in the productive complexes is negligible against the
   total substrate pool, so it need not appear in `u` and `a`.

Assumptions 3 and 4 are the content. Both are checkable, and both are measured
below at the states the paper reports.

### The measurement

`scripts/analysis/closure_audit.py`, evaluated at the 2767 re-solved fold states
of the kinetic box (of 2884 draws admitting a fold; 117 yield no solvable state
and carry no entry).

| quantity | median | p99 | max | count violating |
|---|---|---|---|---|
| chaperone catalytic occupancy `s_ref + s_dis` | 0.455 | 1.736 | 1.970 | 334 of 2767 exceed 1 |
| protease catalytic occupancy `s_u + s_a` | 0.487 | 1.781 | 1.962 | 337 of 2767 exceed 1 |
| neglected soluble substrate, `(C_f s_ref + D_f s_u)/u` | 0.204 | 8.70 | 22.7 | 1863 of 2767 exceed 0.10 |
| neglected aggregate substrate, `(C_f s_dis + D_f s_a)/a` | 0.235 | 12.8 | 38.2 | 1831 of 2767 exceed 0.10 |

Read the two halves separately.

**Assumption 3 fails in about 12% of the ensemble.** A single machine processing
two substrates shares one denominator, so its two catalytic occupancies sum to at
most 1. The coded form permits up to 2, and in 334 networks (12.1%) it exceeds 1:
more of the free pool is committed to catalysis than exists. Independent sites
are a real mechanism — GroEL is a 14-mer, ClpB a hexamer, ClpP a 14-mer — so this
is a stipulation that could be true, not an inconsistency. But it is a
stipulation, and it binds in a tenth of the box.

**Assumption 4 fails outright.** At the median state the productive complexes
would hold a fifth as much monomer as the entire soluble pool, and in 1863
networks (67%) more than a tenth. At the upper tail they would hold twenty times
the pool, which is not a small correction to a conservation law but a
contradiction of one. The classical Michaelis–Menten condition
`E_tot << S_tot + K_M` is not satisfied here, and no rescaling makes it so,
because the pools are of the same order by construction: the model is *about*
machinery being scarce relative to burden.

### The verdict

**Option (b) of task A1.** The rate laws are phenomenological closures. They are
a coherent and conventional way to write "flux is proportional to available
enzyme and saturates in substrate", and every structural hypothesis of Theorem 1
holds for them exactly — `j` is additive in one equation and absent from the
other, and mass balance is exact by construction, since the transfer terms `n`,
`g` and `v_dis` appear in `du/dt` and `da/dt` with opposite signs and cancel
identically. Theorem 1 and its corollaries are therefore untouched.

What is not untouched is the quantitative layer: `phi`, the saturation
fractions, the twelvefold ceiling factor, and the Hopf incidence. Those are
computed *from* these rate laws. Task A3 re-runs them under the alternative in
which catalytic flux is proportional to the bound complex, and the comparison
is recorded in Part 3 below.

---

## Part 3. Robustness to the closure (task A3)

See `theory/RATE_LAWS_A3.md`.
