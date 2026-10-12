# Read-only status recheck — 2026-10-12

## Verdict and structural check

**49/49 structural PASS reproduced; no wet lab, confirmation not ready.**
Executed `python3 CHECK_PROTOCOL.py > STATUS_RECHECK_20261012.txt 2>&1`; exit status 0.
The new log is byte-for-byte identical to `CHECK_PROTOCOL_RESULTS.txt`: `TOTAL PASS=49 FAIL=0`.
The canonical results and all earlier files were left unchanged. No repairs, experimental analyses, simulations, or new experimental work were performed.
The checker verifies structure/provenance, including 502 preexisting hashes; it does not establish causal identification, assay validity, power, or feasibility.
All G0–G6 statuses remain `PENDING_LOCK`, and freeze values/evidence/signatures remain unlocked (`GATE_SPEC.tsv`; `FREEZE_RECORD_TEMPLATE.tsv`).

## Five high-priority design weaknesses

1. **A coherent factorial estimand exists, but an error-specific causal estimand is not identified.**
   Sources: `ONE_EXPERIMENT_PROTOCOL.md` §§1–2, §5 (D, Theta, G1–G2); `ANALYSIS_PLAN.md`, “Estimation, uncertainty and independent limits,” lines 7, 15–17; `SELF_AUDIT.md`, “Physical correction and residual confounds,” line 11.
   Theta consistently defines a 16-term interaction in normalized reporter association. Isogenic repair and common media remove some differences, but cannot exclude allele-dependent metabolic responses to NVA or their interaction with ClpB induction. G1 establishes an error increment; G2 bounds selected measured changes, not every non-error pathway. The audit already acknowledges this limitation, while “error-attributable” and “causal-error primary” overstate it.
   **Necessary correction:** interpret Theta as the ED/REP × NVA/NONE × induction × cassette intervention contrast. Treat attribution specifically to mistranslation as an additional unverified exclusion assumption, even if G1/G2 pass. This stops the stronger causal-error claim, not the defined intervention contrast; no extra arms are proposed.

2. **ClpB is mechanistically relevant, but neither a selective DnaK-capacity intervention nor a guaranteed negative-sign intervention.**
   Sources: protocol §4 (“Functional ClpB check”), §5 (sign hypothesis, G2/G6); audit lines 13–14 and 26.
   ClpB and DnaK cooperate in disaggregation ([primary study](https://www.nature.com/articles/nsmb.1425)); that supports relevance, not equivalence between extra ClpB and isolated DnaK processing capacity. Extract reactivation measures the coupled system, not direct in vivo ClpB throughput. More delivery of clients to DnaK can increase association, whereas faster release can decrease it; the finite-window integral has no forced sign. Feedback and ClpB production cost remain additional pathways.
   **Necessary correction:** retain the explicitly conditional sign test and describe it as extra-ClpB effects on DnaK-associated reporter. Neither sign tests a universal proteostasis law or identifies throughput; a positive sign falsifies only the selected bounded sign hypothesis. The protocol already states much of this correctly.

3. **The proposed five-minute pulse and rapid chase have not established a distinct cohort in prototrophic MG1655.**
   Sources: protocol §§2–3, §5 (Y, G1/G3); analysis lines 11 and 15; audit line 12.
   Endogenous light amino-acid production, uptake, precursor persistence and recycling can change which molecules are labeled and when. Equal extracellular composition does not prove equal intracellular enrichment or cessation of labeling. Normalizing to labeled input at chase zero cannot repair treatment-dependent cohort composition or continued label entry.
   Prototrophic E. coli labeling is possible ([native SILAC study](https://pmc.ncbi.nlm.nih.gov/articles/PMC3708181/)); that does not validate this strain/medium's five-minute Lys/Arg pulse or two-minute first chase interval.
   **Necessary correction:** enforce the existing G1/G3 stop until condition-relevant enrichment, clearance, temporal resolution and nonperturbing exchange are demonstrated. Do not interpret an unqualified heavy signal as the specified pulse cohort or silently change strain/background.

4. **Luciferase is an unqualified primary client here, and its activity cannot validate pulse-cohort productive relief.**
   Sources: protocol §1 (sole reporter), §4 (“Native function”), §5 (G3/G6); analysis line 17; audit lines 9 and 15.
   Its use in disaggregation assays does not establish quantifiable, preserved DnaK association under this non-heat-shock error perturbation. A successful extract reactivation control therefore cannot qualify the primary in vivo association measurement. Total activity also combines old, pulse-labeled and post-pulse reporter, with potentially different specific activities.
   **Necessary correction:** enforce G3 qualification of this reporter before confirmation; failure stops the primary test without reporter substitution. Keep G6 productive relief/new-native yield unavailable unless its independent cohort-function validation succeeds. G6 failure alone does not invalidate an otherwise qualified association contrast.

5. **The 24-arm design has no demonstrated collection capacity or adequate sample size.**
   Sources: protocol §§1, 3; analysis lines 13, 21, 27–31; audit lines 18, 24 and 32–34.
   Three full pilot blocks require 72 cultures, 504 parent specimens and 3,528 planned assay-child records, before dose screening and QC. Each confirmation block adds 24 cultures, 168 specimens and 1,176 child records; these are workload counts, not independent replicates or necessarily instrument injections. Rapid capture, seven-way aliquoting and repeated withdrawal need feasible staffing, volume and throughput. Three lineage/day-confounded pilot blocks give weak variance information; the 16-term contrast and simultaneous G2 equivalence requirements may demand unaffordable n.
   **Necessary correction:** enforce G4 only after the existing sampling/volume limits, operational capacity, variance upper bound, equivalence precision and affordable block limit are documented. The formula and confirmation placeholder are not a power result. Stop if infeasible; no expansion of the design is recommended.

## Measured-data inspection and scope

Recursive inspection, including hidden files, found no new pilot/sample raw files, populated observation sheet, assay qualification results, or signed freeze evidence in this directory.
`SAMPLE_SHEET_TEMPLATE.tsv` has 672 planned rows, with all 231 checker-designated unknown fields per row still `PENDING`; labels `PILOT` and `C_BLOCK_REPEAT` are plans, not collected samples.
Actual historical measurement values do exist in the four `proteinGroups_*_PREVIEW.tsv` files. `TABLE_INVENTORY.tsv` and `SOURCE_DATA_MANIFEST.tsv` identify them as three-row previews of archived Pranjic processed tables, not new experimental data; `retrieval/` contains retrieval records/errors.
The prior Pranjic Gate 3 VOID verdict and absence of an independent confirmatory dataset therefore remain unchanged (`SELF_AUDIT.md`, opening status).
Only this report and the requested stdout/stderr log were created. All other directories, including the central idea ledger, remained read-only; no central append is claimed.

## Summary

Structural integrity is confirmed, not experimental readiness. The intervention contrast is coherent; error-specific causation and productive relief remain unestablished. Enforce the existing assay/feasibility stops and narrow the causal wording before any confirmation claim. End of check.
