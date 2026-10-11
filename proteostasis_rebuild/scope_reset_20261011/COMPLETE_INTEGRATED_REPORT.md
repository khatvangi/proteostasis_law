# Complete integrated research report: translation-associated proteostasis burden

11 October 2026. **Scientific reconstruction finished; new biological validation unavailable. Manuscript verdict: theoretical reassessment and experimental constraints / not submission-ready.** All changes are confined to this directory. Earlier branches retain their historical contents, including claims now withdrawn.

## Original question and answer

The canonical *A finite proteostasis envelope for translation*, by Hannah E. Rebbeck, Janak Paudyal, and Boggavarapu Kiran, asks:

> “First, is the envelope real and finite, and where does extant translation sit inside it? Second, does that position have consequences an experiment can reach?”

This is the manuscript's question (canonical line 46), not a recovered quotation of the user's earlier conversation. Its frozen text is [ORIGINAL_FROZEN_MANUSCRIPT.md](ORIGINAL_FROZEN_MANUSCRIPT.md). The revised manuscript preserves the canonical authors and affiliation.

**Answer supported now:** molecular processing pools are finite, and translation-associated clients can engage experimentally alterable processing and feedback. **Answer not established:** where a specified cell sits relative to a functional service limit, how much workload is caused by translation errors rather than ordinary folding, or whether a specified attainable perturbation crosses such a limit. Neither a numerical physiological headroom nor a universal collapse threshold survives. The empirical primary test is blocked by absent matched measurements, not by unfinished algebra.

The appropriate scope is translation-associated processing at measured synthesis and growth in a stated condition. It is not genetic-code origin, a theory of life, entire-cell demand, or a new network model. Synonymous codons remain a secondary intervention whose effects on timing and substitutions require separate measurements.

## Reconstructed physical logic and minimal quantities

A newly synthesized chain may contain an error, but an error is not automatically a misfolded chain. An error-free chain also requires folding and sometimes disposal. A non-native chain may bind a machine repeatedly, be degraded, aggregate, or reach a native conformation. Binding is not productive rescue, and native conformation is not necessarily activity. Stress regulation can change machinery abundance, but induction requires time and does not automatically supply useful throughput.

| Quantity | Meaning | Units |
|---|---|---|
| s_i | protein-specific synthesis | chains / volume / time |
| π_i(e), p_iY(e) | error-pattern probability; outcome probability under the specified intervention | dimensionless |
| J_Y | rate of a defined client outcome | chains / volume / time |
| w_ij(e) | cumulative occupied site-time per chain, including repeat visits | sites × time / chain |
| B_j, M_j, F_j | occupied, total, free machine sites | sites / volume |
| U_i, C_i, A_i, N_i | non-native, bound, aggregated, native chains | chain equivalents / volume |
| μ and first-order rates | dilution and kinetic losses/exits | inverse time |
| k_on, k_a | capture and quadratic chain-transfer aggregation constants | inverse concentration / time, with explicit chain-to-site stoichiometry |

For fixed synthesis and a common processing environment, subtract the error-free baseline:

`ΔJ_error,Y = Σ_i s_i Σ_(e≠0) π_i(e)[p_iY(e)−p_iY(0)]`.

For a stationary machine workload, make the same subtraction using w:

`ΔB_error,j = Σ_i s_i Σ_(e≠0) π_i(e)[w_ij(e)−w_ij(0)]`.

These are unmeasured causal estimands. They can have either sign. An induced-minus-control contrast must compare complete expressions if regulation, competition, synthesis or outcomes change. Standing substitution frequencies among survivors do not directly give error probabilities at birth. Do not add an error flux, folding flux, aggregate concentration and machine occupancy into `B_error+B_fold+B_agg+B_qc`: the units differ and the same chain may be counted repeatedly.

## Derived conditions and assumptions

The manuscript gives a deliberately minimal one-machine realization. Nonnegative capture is `k_on F U`; aggregation `k_a U²` transfers equal chain mass out of U and into A; release, folding, degradation and dilution have explicit exits. New machinery enters F. With one site per bound chain, `M=F+ΣC`, and `M'=α−μM`. A fixed pool requires replacement `α=μM`; finite pools cannot be maintained indefinitely without synthesis at positive dilution.

Summing client states gives `P'=s−μP−Σd_X X`. Including cumulative degradation and dilution counters gives `(P+D+L)'=s`. The boundary vector field points inward, so initially nonnegative states remain nonnegative. For bounded synthesis and constant μ>0, `P'≤s−μP` bounds total client mass. No universal runaway follows from aggregation. At μ=0 an independent removal argument is needed.

At a steady state with finite residence, `C_i=b_i τ_i`, `τ_i=1/(r_i+f_i+d_C,i+μ)` and `Σb_iτ_i≤M`. Arrival counts visits, not distinct proteins. This is a necessary occupancy condition. Productive formation obeys `J_N,i=(μ+k_u,i+d_N,i)N_i`; an independently justified required N gives a corresponding flux inequality. Neither inequality establishes attainable routing, adequate activity, stability, survival or spare capacity. Background clients use the same measured pools.

A response closure `m'=b(S)−γ_m m`, `α=k_t m` describes how measured active σ32 could drive synthesis. It does not predict S. A constant b, blocked translation or unsuccessful assembly prevents automatic useful induction; a disposal route can increase turnover while decreasing native yield. Effective degradation constants do not establish unlimited protease throughput. All illustrative code witnesses are model consequences, not fitted biology.

Corrected earlier models can have two stable finite states separated by an unstable state under specific sequestration assumptions. That conditional result is acknowledged, parked and not extended. The prior constant-demand branch supplies a useful distinction between feasible and realized output; its whole-cell demand scope and numerical examples are excluded. Its warning against a generic Hall subset test with serial routes is retained.

## Experimental evidence and status

| Evidence | What survives | What it cannot establish |
|---|---|---|
| GroEL/ES rescue of the TF/DnaK-deficient background, V04 | Conditional bypass/cooperation; hard universal serial routing fails | Native-state routing at all temperatures; independent validation of a model built around this observation |
| Calloni DnaK/client experiments, C2012 | Redistribution and condition-dependent aggregation patterns | Client flux, residence, productive obligation, or aggregate mass from identity counts |
| σ32 perturbations, T98/G04/L13/Z14 | Regulated machinery and usable buffering for selected recombinant clients | Native gain, assembly delay, endogenous ceiling or quantitative error burden |
| Schmidt abundance measurements, S2016 | Condition-specific amounts, under stated stoichiometry conventions | Processing capacity without in-vivo residence/turnover |
| Reingewertz mDHFR/GroEL, R24 | Specific synonymous variants alter chaperone dependence | Measured timing mediator, universal codon weights or an error-to-folding link |
| Corrected Landerer analysis, LA24 | Narrow pipeline-dependent AGC contrast | General amplitudes, all synonym contrasts, birth rates, folding or workload |
| Fauvet, H2021 | Recorded 30 °C insolubility point estimates and separate higher-temperature growth observations | Statistically established aggregation-sign refutation or 30 °C growth compensation |

The current kinetic model **fails the Fauvet point-estimate sign**: ΔdnaKJ → ΔdnaKJ ΔhtpG changes insoluble protein 20.5% → 18.9%, while the archived model increases aggregate in 200/200 single-client and 199/200 two-client samples. These are simulation counts, not biological replicates. Insolubility includes membrane proteins. No p-value for the total contrast has been established. This is neither successful validation nor a statistically established refutation.

A provenance correction matters: the later novelty archive contains Fauvet Tables S2–S4 and a supplementary PDF. The old statement that no supplement was retrieved is stale. Those tables report per-protein means and SDs across five biological replicates; they do not resolve replicate-level total-fraction uncertainty or cross-protein covariance. The precise blocker is missing total-contrast uncertainty and matched 30 °C synthesis/growth, not absence of every supplement. See [EVIDENCE_RECONCILIATION.md](EVIDENCE_RECONCILIATION.md).

## Kept, narrowed and rejected

**Kept:** conservative finite-pool accounting; ordinary folding baseline; separate error, outcome, occupancy and function levels; machine-specific necessary conditions; condition-dependent cooperation; recorded negative model confrontation.

**Narrowed:** codon evidence to specific client and pipeline observations; stress feedback to demonstrated qualitative mechanisms and engineered interventions; routing results to model consequences under their assumptions.

**Rejected:** ×25 and all replacement physiological margins, generic no-bounded-state collapse, mismatched-unit additive burden, ribosomal synonymous shielding, all nonsynonymous substitutions equaling misfolding, metal-site enrichment, cross-species operational geometry based on corrupt inputs, hard serial routes, all-client ClpB recovery, portable AGC amplitude, starved in-vivo kinetic calibration, and claims that simulations validate biology. The reassessment's ×63–113 arithmetic is not a substitute measured margin.

## Prior art and novelty boundary

FoldEco already supplies finite chaperone kinetics and competing folding/degradation/aggregation; Santra adds growth-dependent folding demand; FoldME constrains folding within a growth model. FoldEcoSlim explicitly includes growth dilution, background-proteome competition and σ32-linked induction. The physiological feedback mechanism also predates this work. The references and access levels are retained in the manuscript and evidence matrix.

No priority claim is made for a capacity network, occupancy law, kinetic routing, feedback, growth coupling, or a combination of these. The surviving work is a transparent reassessment and an experimentally operational definition of the missing increment. Its independent publication novelty is not established merely by correct mathematics.

## Two competing predictions and falsifiers

Use one matched, calibrated error perturbation, orthogonal fidelity rescue, measured synthesis and error spectra, and a prespecified client/machine. Measure occupied site-time and activity early and after induction. Cross error induction with independently controlled functional machinery expression. Prespecify a meaningful workload increment δB and uncertainty from the assay before data collection; none is currently supplied by the records.

**A — material additional processing:** the verified error increment causes ΔB above δB and a machinery-dependent activity effect. Fidelity rescue reverses the error-linked component. A sufficiently precise upper confidence bound below δB falsifies this material-load prediction in the tested condition. Confounding by drug injury or synthesis makes the causal test undecidable.

**B — no material incremental demand on the tested machine:** ΔB remains below δB despite verified extra errors, although reporter or activity changes may occur through other routes. A reproducible lower bound above δB, reversible by fidelity rescue with synthesis and injury controlled, falsifies B. A confidence interval crossing δB decides neither. B does not imply errors are harmless.

Report the factorial interaction `I=y11−y10−y01+y00` on a preselected scale, with uncertainty. Nonzero interaction alone establishes neither saturation nor a threshold; synthetic lethality is not presumed. These are competing interpretations of one bounded experiment, not a new future-work programme.

## Exact blockers, publishing decision and decisive next step

1. No matched causal measurement of per-synthesis error spectrum, ordinary-client baseline and client-specific occupied site-time. This blocks the primary translation-error workload result.
2. No matched in-vivo arrivals/residence and functional machinery calibration. WO-07's matched rows chiefly supply abundance/growth; neither exponential nor stationary bundles support physiological headroom. Stationary synthesis cannot be inferred from `μP` under negative net growth.
3. Fauvet total-fraction replicate uncertainty/covariance and 30 °C growth/synthesis remain unresolved. Per-protein SDs and growth at other temperatures do not close this gap.
4. Original Calloni/Kerner client tables and the Babu interaction table are not authenticated in the available evidence ledger. The historical Babu retrieval failed; missing pairs in a significant-hits-only table would not count as neutral even if retrieved. These are separate retrospective blockers, not prerequisites for inventing new theory.
5. V04 duplicate-publication provenance, raw-spectrum localization control for the AGC contrast, and novelty adequate for an original research submission remain unresolved.

**NO-GO for submitting the original quantitative envelope paper or presenting this revision as validated original empirical research.** The revision is a complete theoretical reassessment and experimental-constraints manuscript. A future submission as a critical reassessment would require editorial fit and resolution of its source/novelty issues, not relabeling conditional examples as biology.

**Decisive next step:** obtain the single matched error-induction/fidelity-rescue × controlled-machinery dataset described above. Before that, this project cannot answer where the cell sits or whether its error-attributable load is materially capacity-limited. Existing-data clarification of Fauvet may reject a kinetic closure but cannot substitute for that measurement. No further simulations are required to finish this reconstruction.

## Completion boundary

The delivered mathematical run passed **20/20 checks**, with zero failures, errors or skips; **15/15 package checks** also passed. Full output is recorded in [CHECK_RESULTS.txt](CHECK_RESULTS.txt) and [DELIVERY_VERIFICATION.txt](DELIVERY_VERIFICATION.txt); the manuscript, evidence matrix and adversarial audit form one complete package. Reconstruction, claim disposition and specification of the decisive test are finished. New wet-lab validation was not performed and remains unavailable. No modifications to earlier branches or the external idea ledger were made; the user required writes only here. An append-ready ledger record is preserved locally.

**STOP / FINISHED SCIENTIFIC SCOPE: bounded reconstruction complete; empirical validation outstanding, not silently claimed complete.**
