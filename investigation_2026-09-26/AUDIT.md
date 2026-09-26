# Independent scalar and finite-pool validation — 2026-09-26

Scope is proteostasis P1 only. This directory is a new, non-destructive audit;
no pre-existing project file, raw data, manuscript, or git metadata was changed.
The audit reads the upstream source at `../proteostasis-P1/two_pool_ode.py`, the
repair plan at `envelope-paper/PLAN_P1_REPAIR.md`, and the active manuscript
notes. It does not treat old generated outputs as validated evidence.

## Verdict

The scalar redraw has been corrected and independently rechecked. The
stationary-point cubic is `-0.30x^3+0.40x^2+1.70x+5=0`; the prior `+9.60x`
coefficient is rejected. The equilibrium polynomial is
`-chi*x^3+(1-chi)*x^2+(1+rho-lambda)*x-lambda=0`; for `lambda=2` the prior
`rho+1` coefficient is rejected. All plotted positive roots below satisfy
`abs(g(root)-lambda)<1e-9` when evaluated in the original rational function.

| Item | Observation / derivation | Status | Consequence |
|---|---|---|---|
| Two-pool state definition | Source lines 5–6 define `P` and `A` as proteome fractions, while lines 9 and 16 use `M=P Prot_tot` and concentration-based rates. | observed | The state variables are dimensionless fractions, but their physical pools and conversion rates need explicit bookkeeping. |
| `C_free` closure | Source lines 127–129 use `C_tot/(1+M/Kd)` with `M` in µM. This is a ligand-excess/infinite-ligand-style approximation, not finite-pool mass balance. | observed | It violates the stated finite chaperone pool when bound chaperone is non-negligible. The repair plan independently flags this at lines 206–213. |
| Correct finite-pool closure | With `C_T=C_f+C_b`, `M_T=M_f+C_b`, `K_d=C_f M_f/C_b`, exact `C_b` solves `C_b²-(C_T+M_T+K_d)C_b+C_TM_T=0`. | derived | At baseline `C_T=50`, `K_d=1` µM, `M_T=50` µM: `C_f=6.58872`, `C_b=43.41128` µM, versus approximate `C_free=0.98039` µM. |
| Folding rate impact | Source lines 131–133 define `v_fold=k_obs_max C_free/(C_free+K_d)`. | observed + derived | With the same `k_obs_max`, exact/approximate folding rate is 1.75382× at `M_T=50` µM. This is a model correction, not a viability result. |
| Aggregation/drain units | Source lines 146–150 give `drain=k_agg P² Prot_tot`, with units fraction/s if `k_agg` is M⁻¹s⁻¹ and `Prot_tot` is M. | derived | This is dimensionally consistent for a fraction state, but the assumed equality `k_nuc=k_agg` and the mapping from concentration-dependent aggregation to a fractional drain are imposed, not measured here. |
| Inflow amplification `Phi` | Source lines 140–141 define `Phi=1+v_agg/(v_fold+k_deg)` and line 14 sets `J_in=J_bare Phi`. | observed + provenance audit | This is a closure/ansatz in the implementation. The source does not derive it from nascent-chain occupancy, measured feedback, or a mass-action network. |
| Aggregation feedback sign | The scalar manuscript equation (MANUSCRIPT.md lines 117–133; PHYSICS_FRAMEWORK.md lines 50–66) adds `+chi x²` to `dx/dtau`. | observed + derived | This is positive feedback in the reduced scalar model. It is not the same mathematical operation as an aggregation drain that removes material from `P` into `A`. |
| Quasi-steady aggregation pool | Source lines 22–25 eliminate `A` algebraically; lines 169–185 use the resulting `k_clear A` loss. | observed | This is a quasi-steady reduction, not a demonstrated timescale separation. `A_max=0.25` is an operational imposed cap (lines 90–112), not an experimentally established universal viability threshold. |
| Dimensional mapping to codon frequency | Source lines 261–263 convert `J_bare` to `f_codon` using `T_gen`, `N_prot`, `(1-S_avg)`, and `p_baseline`. | observed + imposed | The conversion is algebraically transparent, but its parameter meanings and phase matching are inputs. The repair plan lines 70–90 specifically warn that stationary-phase error rates can mismatch exponential-growth capacity and dilution. |
| Empirical provenance | The upstream script loads an absolute-path proteome table at lines 69–83; Monte Carlo ranges are imposed at lines 392–412. The repair plan lines 13–42 says staged result JSON files were outputs, not generators. | observed | Do not call the old threshold ensemble empirical validation. Its outputs are model-generated conditional on imposed ranges and closures. |
| Functional provenance | The model has no explicit nascent-chain state, ribosome occupancy, substrate classes, or competing chaperone clients; `C_tot` is a single effective pool. | observed | State-dependent regulation and nascent-chain competition are absent, so functional claims about cellular chaperone allocation are not identified by this ODE. |

## Scalar rederivation

The reduced model is

`F(x) = lambda - x - rho*x/(1+x) + chi*x*x = lambda - g(x)`

with

`g(x) = x + rho*x/(1+x) - chi*x*x`.

For `rho=4`, `chi=0.15`:

`g'(x) = 1 + 4/(1+x)^2 - 0.30 x`.

For `x>-1`, stationary points solve the cubic obtained after multiplying by
`(1+x)^2`:

`-0.30 x^3 + 0.40 x^2 + 1.70 x + 5 = 0`.

There is exactly one positive root in this case:

`x_m = 3.890758215117349`, `g(x_m) = 4.80218919587308`.

Thus `lambda_fold = 4.80218919587308` is the maximum admissible positive
horizontal level for the low-load equilibrium in the nonnegative domain. The
positive zero of `g` is `x=9.2645937939` (the other roots are `x=0` and
`x=-3.5979271272`). For `lambda=2`, the positive equilibria are
`x=0.5808674541` and `x=7.9669676044`; direct substitution gives residuals
below `1e-12`. The third algebraic real root is `x=-2.8811683919`, outside
the nonnegative physical domain, and also has residual below `1e-12`. Because
`F'=-g'`, the low root is
locally stable (`g'>0`, `F'<0`) and the high root is unstable (`g'<0`,
`F'>0`). For `lambda>g_max`, there is no positive equilibrium below the
overload crossing; the scalar flow is positive there and the reduced model
runs toward its unbounded/invalid region. This is the plotted overload
behavior, not evidence of a biologically realized state.

A fold requires a simultaneous solution of `F=0` and `F'=0`, equivalently a
turning point of `g` at a level `lambda=g(x)>0`. A sufficient local shape
condition is that `g'` becomes zero and changes sign on the physically allowed
domain. It is not guaranteed by the words “saturating removal plus positive
feedback” alone: parameter values, domain, and signs matter. In particular,
for `chi<=0` the large-x behavior does not produce the same positive-feedback
turning structure, and for `chi>0` a turning point can still lie outside the
allowed domain or at nonpositive `g`.

The `+chi x²` term is a positive-feedback ansatz in the scalar balance. By
contrast, the two-pool source has a drain from `P` to `A` (lines 16, 19–20),
which is a transfer/removal term from the monomeric pool, followed by aggregate
clearance. Calling both “aggregation” would conflate distinct mechanisms.

## Finite-pool comparison

Baseline values are taken from source lines 97–109: `C_T=50 µM`, `K_d=1 µM`,
`Prot_tot=300 µM`; the sweep uses `M_T=0…300 µM`, a feasible nonnegative total
misfolded-ligand range up to the source's total proteome concentration. The
approximation is `C_T/(1+M_T/K_d)`. The exact calculation allows the bound
complex to consume both pools. The rate comparison uses the source's functional
form but reports relative rate with `k_obs_max=1`; choosing another common
`k_obs_max` scales both curves equally.

Selected computed values:

| `M_T` (µM) | approximate `C_f` | exact `C_f` | exact `C_b` | exact/approx. folding rate |
|---:|---:|---:|---:|---:|
| 0 | 50.00000 | 50.00000 | 0.00000 | 1.00000 |
| 10 | 4.54545 | 40.24247 | 9.75753 | 1.19042 |
| 25 | 1.92308 | 25.92839 | 24.07161 | 1.46355 |
| 50 | 0.98039 | 6.58872 | 43.41128 | 1.75382 |
| 100 | 0.49505 | 0.96224 | 49.03776 | 1.48094 |
| 300 | 0.16611 | 0.19905 | 49.80095 | 1.16534 |

The difference peaks around stoichiometric competition, not at the largest
ligand load. The exact/approximate folding-rate ratios at `M_T=0, 10, 50,
300 µM` are `1.00000, 1.19042, 1.75382, 1.16534`. These are numerical consequences of alternate binding closures;
they do not by themselves imply improved or worsened viability.

## Observed / derived / speculative boundary

- **Observed in files:** equation forms, parameter defaults/ranges, absolute
  data path, quasi-steady elimination, operational `A_max`, and manuscript
  language claiming a generic fold.
- **Derived here:** units of the stated rate expressions, scalar stationary
  point and threshold, stability classification, exact finite-pool root, and
  numerical curve differences.
- **Speculative or unestablished:** universal viability threshold, biological
  identification of `Phi`, equivalence of `k_nuc` and `k_agg`, fitted values of
  `chi`/`rho`, general fold behavior in cells, and any experimental support for
  model-generated thresholds.

## Open experimental discriminators

1. Measure free and bound chaperone simultaneously across a calibrated
   misfolded-client titration; test the finite-pool mass balance against the
   approximate closure.
2. Measure `P`, `A`, and fluxes into/out of `A` over time to test whether the
   assumed quasi-steady `A` reduction is valid on the growth/damage timescale.
3. Perturb aggregation nucleation independently of misfolded monomer abundance
   and test whether the inferred monomer-to-aggregate flux is proportional to
   `P²`.
4. Perturb chaperone abundance while measuring folding, degradation, and
   nascent-chain occupancy; distinguish increased free-pool rescue from
   competition by ordinary nascent chains.
5. Measure whether the proposed amplification of damage inflow tracks
   `1+v_agg/(v_fold+k_deg)` or instead follows an independently regulated heat-
   shock / protease response.
6. Repeat burden sweeps in matched exponential and stationary conditions with
   phase-matched dilution, chaperone pool, and error-frequency measurements.
7. Test for hysteresis or a saddle-node signature by controlled up/down burden
   sweeps; a steep response without hysteresis would not establish the scalar
   fold mechanism.
8. Resolve whether any codon-linked phenotype enters through error, folding,
   nascent-chain handling, or regulation using matched synonymous constructs and
   direct proteome/QC readouts.

## Reproducibility

Run `python redraw.py` from this directory. It writes the four requested plot
files and prints the computed scalar and binding values. Run `pytest -q` to
execute the minimal tests. Full equations, assumptions, and units are in
`METHODS.md`.
