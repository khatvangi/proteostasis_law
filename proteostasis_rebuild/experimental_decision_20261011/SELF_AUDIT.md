# Completion audit

**no wet lab, confirmation not ready.** All templates describe planned collection. No observations, assay qualification, fitted models, simulations, power result or author contact were produced. Structural PASS is not a G0–G6 experimental pass. Every freeze value and signature remains PENDING_LOCK.

The public-data verdict is unchanged: no confirmatory reanalysis from the Pranjic processed tables, no independent holdout, and prior Gate 3 remains VOID. Historical early no-contact text is superseded. Missing concurrent WT factorial cell, separate WT batch, unauthenticated pairing, absent Nva DnaK arm and missing PD mapping cannot be repaired by relabeling. Stage1 already iterated outcomes and wrote previews.

## Physical correction and residual confounds

Total activity and pulse-labeled abundance are separately measurable. Pulse-cohort functional yield is UNIDENTIFIABLE without an independent pulse-cohort activity assay or a validated conditional synthesis-stop preserving physiology. Unlabeled post-pulse synthesis and surviving old reporter contribute activity; error variants may differ in specific activity. No pulse-cohort active-count interval or yield division from total activity is valid here. G6 first requires independent identifiability and specific-activity validation. No synthesis-stop intervention was added. Y, D, Theta, the occupied-site upper-bound statement and sample-size equation are preserved.

- Error source: allele/exposure interactions can affect metabolism beyond mistranslation. Repair and NONE do not remove all interactions. Nva/Val are separate interventions; isomers require same-peptide standards. Off-target spectrum coverage remains incomplete.
- Synthesis: precursor uptake/clearance, recycling, response factors, modifications, DNA variants, localization and detection limits affect nascent estimates. Abundance is not synthesis; surviving pulse cohorts can be selected by degradation.
- DnaK: IP measures association, possibly indirect. Ex vivo exchange, lysis, fixation, antibody specificity, recovery and saturation need validation. Association is not occupied sites, microscopic dwell, throughput or native yield. Rebinding is unresolved; analytical normalization must retain biological DnaK variation.
- ClpB: production cost may change growth, translation, stress and errors; EMPTY removes inducer effects but not production cost. ATP-dependent extract activity is not in vivo capacity. Native ClpB and DnaK partners remain feedback responsive. Disaggregation can transiently increase DnaK association.
- Reporter: the client may perturb physiology or lack measurable DnaK association. Sequence/copy number, transcription, specific activity, matrix interference and substrate excess need validation. One luciferase reporter cannot establish a proteome-wide law or universal sign.
- Losses: aggregate stock is not formation flux; contamination and incomplete recovery affect mass balance. Disappearance can reflect dilution, degradation, recycling, death or extraction. Release and aggregate loss do not establish folding. Residual cohort at 40 min is censored, not extrapolated to zero.
- Feedback: RpoH/RpoS timing is descriptive, not causal mediation. Isotope exchange, medium, aeration, induction and repeated withdrawal can themselves alter physiology.
- Design: three pilot lineage/day blocks confound lineage and day and estimate variance imprecisely. Confirmation must balance lineages and respect clustered terminal sisters and repeated aliquots. Technical splits are not n. Batch confounding, unblinding and outcome-dependent exclusions threaten interpretation.
- Inference: shared calibration uncertainty does not vanish with n. Equivalence requires simultaneous intervals within independently justified margins. Pilot Theta cannot choose delta, dose or the powered effect. Missing/censored times are never zero-filled. Large cancelling control effects do not prove specificity; secondary outcomes cannot rescue primary failure.
- Sources: archived excerpts and metadata are not full supplementary SOP access. S9 and S11 distinguish local citations from full text. Transport failure does not prove global data absence. Only local provenance was checked in this completion.

## Stop and reporting rules

GATE_SPEC.tsv reproduces the corrected protocol rules verbatim. G0 requires fresh confirmation cultures independent of tuning. G1 failure stops a causal-error primary claim. G2 failure restricts interpretation to total perturbation effects. G3 failure stops the primary mechanistic observable; no aggregate substitute or reporter switch. Unresolved freeze, inadequate precision or unaffordable n fails G4 and stops confirmation. The confirmation placeholder is not an allocation or observed replicate.

After earlier gates pass, G5 uses the wider biological versus calibration-inclusive interval: upper 95% CI(Theta)<−delta supports the bounded sign; lower CI>−delta rejects reduction of at least delta; lower CI>0 also falsifies the sign. Equal boundaries or overlap remain inconclusive. G6 permits productive relief only with identifiable increased cohort function and bounded competing losses; otherwise association relief only, no new-native yield claim. Nonsignificant degradation change does not exclude degradation.

No interim significance stopping, endpoint changes, dropping detected substitutions or selective replacement. Repeat failed assays from retained material while blind; irrecoverable cells/blocks remain missing unless the frozen complete-block replacement rule applies. Report every attempted block and conservative censoring bounds. Missing times that prevent a robust G5 decision imply inconclusive.

## Template contract and preservation

The sample sheet has 672 data rows plus a header: four blocks × 24 cultures × seven times. P01–P03 are planned lineage slots, not observed lineage identities. C_BLOCK_REPEAT must be expanded to locked n with fresh unique IDs and balanced verified lineages. All unknowns are PENDING; NA means baseline chase time only. Zeros denote planned NONE added exposure, LOW added IPTG and chase zero, never observations. Proposed 30 C and pulse 20–25 min are design settings, not measurements.

Seven namespaced child metadata groups, five peptide-field groups, activity fields and mass-balance fields are included on each parent row. Expand the peptide slots to long-form records as needed; they are not a one-peptide limit. Preserve parent/child links. There are 4,704 planned child IDs. Separate batch-linked blank, pooled-reference and ATP-omission QC records and a separate PILOT_CAL screening manifest are required before collection; they are not additional factorial rows.

The stdlib checker performs structure/provenance checks only. It checks all 502 PREEXISTING_HASHES.tsv records and protected local documents. Only protocol, analysis plan, README and the append-only local LEDGER_PENDING.md are intentionally updated. All writes stay in this directory. The central ledger is read-only under this task; the correction is staged locally, not claimed as centrally logged.
