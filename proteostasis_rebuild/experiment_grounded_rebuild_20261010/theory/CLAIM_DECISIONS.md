# Claim decisions

What was kept, narrowed, dropped and refused in building this theory, with the reason and the
source. "Earlier work" means the canonical `../envelope-paper/manuscript/MANUSCRIPT.md` and the
`proteostasis_rebuild` work orders and reassessment, all of which are left unmodified.

## Dropped

| claim in earlier work | decision | reason |
|---|---|---|
| usage-weighted per-codon error rate 6.3e-4 as *the* E. coli rate | **dropped** | reproduces the published table only by discarding 3,015 covered zero cells; it is a mean conditional on detection. Zero-inclusive 2.5e-4, pooled 1.4e-4 (`paper_reassessment_20261010/codon_bridge`) |
| ×25 margin, "one order of magnitude inside the envelope" | **dropped** | computed at the rejected evaluation point; and independently, the occupancy condition's left side (flux × service time) has never been measured in vivo, so no margin can be quoted at all |
| "two orders of magnitude" at the bottom of the error window | **dropped** | same reason, plus the window's lower end was already the point the earlier draft rejected |
| all synonymous-codon contrasts except serine AGC | **dropped** | the all-data per-dataset sign test has a 100% family-wise false-hit rate under equal rates at observed depths (mean 6.5 false hits per simulation); proline CCA–CCG has no support |
| a same-synonym substitution error changes the folding outcome of the native protein | **dropped** | no evidence anywhere checked. A substitution makes a different protein; a synonymous change does not change the protein made |
| metal-binding-site codon deployment; cross-species conservation of operational geometry | **stay dropped** | rejected on evidence in earlier rounds; not reintroduced |
| hard serial requirement TF→GroEL or DnaK→GroEL for obligate clients | **dropped as a hard constraint** | contradicted by GroEL/ES overproduction rescue of Δtig ΔdnaK at 30 °C (Vorderwülbecke 2004). Replaced by a conditional, finite direct route |
| all-client radial demand normalization (every aggregate must be productively recovered) | **dropped** | set the demand scale by the smallest pool (ClpB), putting the folding machines in large excess and leaving the central supra-additivity prediction untested. Aggregation and degradation are legitimate outcomes, not failures to be quota'd away |
| ClpB as the binding bottleneck | **dropped** | that was abundance × an assumed aggregate share. Which machine binds depends on service time, which is unmeasured |
| "98.7% single cut" and the multi-stage saturation ordering | **dropped** | artefacts of hand-set route probabilities and a preference ordering (`network_failure_test/REPORT.md`) |
| the cell optimizes / allocates efficiently | **refused** | no evidence licenses it here, and §3 shows an optimizing model makes the wrong call on chaperone-deletion suppression. The word "optimizes" does not appear as a claim |
| universal collapse thresholds, or several of them | **refused** | the threshold exists only with a sequestration coefficient β that has never been measured in vivo |
| a death cutoff | **refused** | failure is defined as no state meeting the quota; what happens then needs more equations |
| temperature as a monotone load multiplier | **refused** | Kanemori's plating data is non-monotone in temperature |
| calibrated E. coli parameters | **refused** | abundances are exponential-phase, kinetics are in vitro, error rates are standing frequencies. No phase-matched set exists |

## Narrowed and kept

| claim | kept as |
|---|---|
| detected mistranslation is organized at the amino-acid level | one depth-calibrated result: serine AGC decoded 12–17× less accurately than TCC/TCT, 96% via Ser→Asn. One amino acid, one growth phase, one pipeline, and a near-isobaric artefact is not excluded |
| codon choice matters for chaperone demand | kept as a *timing* variable (co-translational partition), supported only by Reingewertz 2024 in a heterologous client, where the translation-rate link is computational, not measured. Explicitly separate from substitution frequency |
| proteostasis is a finite-capacity problem | kept, and explicitly credited to FoldEco / Santra / FoldME. Not claimed as original |
| client classes with qualified routes | kept, with route qualification treated as condition-dependent per Vorderwülbecke/Genevaux, and per-protein rather than per-class per Fujiwara |
| stationary-phase error rate 1.8e-3 | kept, reported separately; never pooled with exponential |

## Deliberate non-assumptions

The theory does not assume that every disturbance is maintainable, that chaperones are beneficial,
that deleting capacity can only hurt, that growth follows native function rather than aggregate
toxicity (§4 leaves this open on purpose), or that single-client sign results transfer to shared
saturated pools (listed as open).

## Preserved

Nothing was deleted or rewritten. `../envelope-paper/manuscript/MANUSCRIPT.md`, the work-order
directories, `network_failure_test/` and `paper_reassessment_20261010/` are untouched; this
directory is additive.
