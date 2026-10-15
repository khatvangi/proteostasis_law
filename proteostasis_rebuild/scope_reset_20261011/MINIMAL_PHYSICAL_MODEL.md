> **ARCHIVED 2026-10-14 — CLOSED, DO NOT REVIVE. Read /storage/kiran-stuff/proteostasis_law/000_CLOSED_DO_NOT_REVIVE.md first. Numbers in this file include withdrawn claims.**

# Minimal physical accounting

Date: 2026-10-11. **Derived conditional bookkeeping**, not fitted biology. No universal threshold and no measured in-vivo headroom.

## 1. Condition first; subtract the error-free baseline

Fix measured synthesis `s_i` (chains·volume⁻¹·time⁻¹) for each protein, strain, medium, temperature, phase and measured dilution `μ`. Compare error induction with a matched control at those rates. If induction also changes synthesis or damages mature proteins, those changes are separate contributions, not translation-error load.

Let `E=e` denote a specified upstream translation-error pattern, with `e=0` error-free; `π_i(e)` is its probability per newly made chain. Let `Y` be a defined protein outcome within a fixed observation interval, such as entry into non-native handling or degradation. Define `p_iY(e)=P(Y | do(E=e), i, fixed condition)`: the counterfactual outcome if the error pattern were imposed, not an observational association among surviving proteins.

By total probability, under a common specified processing environment:

```
J_Y = Σ_i s_i Σ_e π_i(e) p_iY(e)
J_Y^0 = Σ_i s_i p_iY(0)
ΔJ_error,Y = J_Y − J_Y^0
          = Σ_i s_i Σ_(e≠0) π_i(e) [p_iY(e) − p_iY(0)].
```

For binary error/no-error this is `Σ_i s_i ε_i[P(Y|do(error),i)−P(Y|do(no error),i)]`, with the error mixture specified. Induced-minus-control replaces `π_i(e)` by `π_i^ind(e)−π_i^ctrl(e)` **only if** the conditional outcome probabilities are unchanged; otherwise compare the full two counterfactual expressions. Error-free chains already require folding, chaperone service and degradation. Neither `p_iY(0)=0` nor `p_iY(error)=1` is licensed; the difference need not be positive. Standing substitution frequencies do not directly identify birth probabilities because processing selects survivors.

Separate four levels: **error event → protein outcome → machine-bound occupancy → functional output**. Outcome flux is not workload in machine units. For machine `j`, define `w_ij(e)` as expected cumulative occupied site-time per synthesized chain, counting repeated visits. At stationarity,

```
B_j = Σ_i s_i Σ_e π_i(e) w_ij(e)
ΔB_error,j = Σ_i s_i Σ_(e≠0) π_i(e)[w_ij(e)−w_ij(0)].
```

`B_j` has occupied-site concentration units. Keep a vector over machines; do not add error flux, folding flux, aggregate concentration and occupancy. Background mature-protein maintenance must also occupy the measured pools. These counterfactual differences are **estimands, unmeasured here**; changed competition and feedback alter `w`, so a fixed-environment subtraction is not a universal additive decomposition.

## 2. Exact finite-pool dynamic realization

For clients `i`, use concentrations of free non-native chains `U_i`, aggregated chain equivalents `A_i`, machine-bound chains `C_i`, and native chains `N_i`. One machine pool has total service sites `M`, free sites `F=M−Σ_i C_i`. One bound chain occupies one site; oligomer-to-site conversion must be specified experimentally. `q_i s_i` enters `U_i`, `(1−q_i)s_i` enters `N_i`; error-free `q_i` need not vanish. Index error products separately if their kinetics differ.

Define capture `b_i=k_on,i F U_i`; aggregation transfer `a_i=k_a,i U_i²` in **chain-equivalent flux** units (not aggregate-particle number); release `r_i C_i`; productive exit `f_i C_i`; degradation `d_X,i X_i`; spontaneous folding/unfolding `k_f,i U_i`, `k_u,i N_i`. All constants are nonnegative; `k_on` and `k_a` have inverse-concentration/time units, other rates inverse time.

```
U_i' = q_i s_i + k_u,i N_i + r_i C_i
       − b_i − k_f,i U_i − d_U,i U_i − a_i − μU_i
C_i' = b_i − (r_i+f_i+d_C,i+μ)C_i
A_i' = a_i − (d_A,i+μ)A_i
N_i' = (1−q_i)s_i + k_f,i U_i + f_i C_i
       − (k_u,i+d_N,i+μ)N_i
F'   = α − μF − Σ_i b_i + Σ_i(r_i+f_i+d_C,i)C_i
M'   = α − μM.
```

Thus `M=F+ΣC_i` exactly, with nonnegative states preserved. Fixed `M` requires replacement `α=μM`. This closure assumes no machine decay beyond dilution and releases machines when clients are degraded. Let `P_i=U_i+C_i+A_i+N_i`, `D_i'=Σ_X d_X,i X_i`, and `L_i'=μP_i`; then `(P_i+D_i+L_i)'=s_i`. Aggregation transfers mass from U to A; it creates none. `D,L` are cumulative bookkeeping sinks, not extra concentrations available for processing. No all-client recovery or ClpB route is assumed. Degradation rates are conditional effective rates, not a claim of unlimited protease capacity; quantitative protease competition needs its own measured finite pool.

If the observation window includes transcriptional regulation, minimally replace fixed `α` by `α=k_t m`, `m'=b(S)−γ_m m`, where `S` is measured active σ32, `b` transcription, and `γ_m` transcript loss including dilution. Newly synthesized machine enters F. This is an **unfitted regulatory closure**, not a model predicting σ32. Apply it only to regulated machinery; TF is outside this regulon. Measure assembly delays rather than assuming instantaneous functional induction.

## 3. Conditional steady-state inequalities and limits

At a nonnegative steady state with fixed measured conditions and finite residence times:

```
C_i = b_i τ_i,  τ_i = 1/(r_i+f_i+d_C,i+μ)
Σ_i b_i τ_i ≤ M                         [per machine]
J_native,i = (μ+k_u,i+d_N,i)N_i
J_native,i ≥ (μ+k_u,i+d_N,i)N_req,i      [only if N_i≥N_req,i is required].
```

`N_req,i` must come from an independent functional assay, not a fitted death cutoff. Native conformation need not imply activity for a mistranslated product; measure activity separately. These identities/bounds do **not** establish attainable routing, stability, adequate function, growth, survival, a fold or spare capacity. A transient need not satisfy stationary flux balance. In this minimal closure `P_i'=s_i−μP_i−Σd_X X_i`; constant positive μ bounds client mass, so aggregation alone does not imply runaway. At μ=0, boundedness requires separately established removal.

Source anchors: [finite-pool audit lines 20–30](../../investigation_2026-09-26/AUDIT.md#L20); [theory lines 54–106](../experiment_grounded_rebuild_20261010/theory/SCIENTIFIC_THEORY.md#L54); [state-resolved kinetics lines 30–48](../experiment_grounded_rebuild_20261010/kinetics/RESULTS.md#L30); [σ32 evidence lines 11–38](../experiment_grounded_rebuild_20261010/feedback/EXPERIMENTAL_CONTROL.md#L11). Equations specify a minimal admissible accounting choice; their parameter values and causal workload differences remain conjectural. They do not resolve the Fauvet held-out sign mismatch.
