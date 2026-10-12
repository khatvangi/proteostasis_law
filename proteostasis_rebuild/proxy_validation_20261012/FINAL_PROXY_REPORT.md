# Final exploratory proxy report

The specified proxy does **not** show a positive median Nva × ClpB-deletion interaction. Across 1,489 matched complete-positive proteins, the median is **−0.4145 log2 units**, quartiles **−1.2125 and +0.2294**, and **32.50%** are positive. This is a descriptive relative-redistribution result, not evidence that deletion lowers total aggregation. The Nva-versus-Val comparison provides no clear positive Nva-specific location shift under this estimand.

## Status and scope

PROTOCOL.md was written and hashed before this run inspected raw outcome headers or values. The folder already existed with a literature dossier; this run adds an executable archival analysis and preserves those files, with the prior README saved as README_PREEXISTING.md. Existing knowledge and historical analyses make this exploratory, not genuine preregistration confirmation. PROTOCOL_FREEZE.json and OUTCOME_READ_RECEIPT.json document sequence, not third-party attestation. The historical GATE3_STATUS.md remains VOID and unmodified. All writes were confined to this directory; no papers were rewritten and no subagents were used.

## Estimand

Rows flagged reverse, contaminant or site-only were excluded, as were empty/multiple accessions, groups with Number of proteins other than one, and duplicate accessions. Exact accessions were joined across fractions. Primary inclusion required positive finite iBAQ in all 32 measurements: four labels × four genotype/treatment arms × two fractions. Zeros were not imputed.

For the fixed matched set, each sample was divided by the summed iBAQ of that same set in its own fraction. Shares were averaged separately within condition for each fraction; Y = log2(mean insoluble share / mean proteome share). The interaction is [Y(ED Delta clpB,Nva) − Y(ED Delta clpB,blank)] − [Y(ED,Nva) − Y(ED,blank)]. There is no A-D pairing assumption across fractions. Positive values represent a larger Nva-related *relative* insoluble enrichment in the deletion background.

This normalization discards fraction-wide scale: multiplying every aggregate intensity in any sample by any positive factor leaves its shares unchanged. Thus an increase of every protein's aggregate mass could be invisible. Conversely, a few abundant proteins changing can move other proteins' shares. iBAQ shares are not calibrated aggregate mass fractions, and pellet contents can include constitutively insoluble material.

## Primary and predefined specificity results

| analysis | n | q25 | median | q75 | positive_fraction | minimum | maximum |
| --- | --- | --- | --- | --- | --- | --- | --- |
| primary_Nva | 1489 | -1.2125 | -0.4145 | 0.2293 | 0.3251 | -14.4044 | 3.7775 |
| secondary_Val | 1484 | -1.4484 | -0.4188 | 0.5275 | 0.4003 | -4.8273 | 6.0871 |
| common_Nva | 1437 | -1.2253 | -0.4389 | 0.2145 | 0.3194 | -5.5807 | 3.7773 |
| common_Val | 1437 | -1.4487 | -0.4308 | 0.5180 | 0.3987 | -4.8277 | 6.0864 |
| common_Nva_minus_Val | 1437 | -0.4976 | -0.0252 | 0.3517 | 0.4857 | -3.7764 | 3.9679 |

Val uses its own complete-positive set. The common set requires positivity in all 48 measurements, allowing a like-for-like Nva-minus-Val comparison. Its median difference is −0.0252, with 48.57% positive; absence of a large descriptive shift is not an equivalence test. The full-primary minimum (−14.4044) disappears in the more restrictive common set, illustrating sensitivity of individual extremes to coverage; no trimming or outlier exclusion was performed.

## Coverage and QC

| universe | n | eligible_shared | fraction_shared | aggregate_complete | proteome_complete | required_measurements_per_protein |
| --- | --- | --- | --- | --- | --- | --- |
| primary | 1489 | 2221.0000 | 0.6704 | 1787.0000 | 1932.0000 | 32.0000 |
| secondary | 1484 | 2221.0000 | 0.6682 | 1745.0000 | 1986.0000 | 32.0000 |
| common | 1437 | 2221.0000 | 0.6470 | 1721.0000 | 1913.0000 | 48.0000 |

| fraction | stage | n |
| --- | --- | --- |
| aggregates | raw_rows | 2920 |
| aggregates | excluded_reverse_overlapping | 120 |
| aggregates | excluded_contaminant_overlapping | 47 |
| aggregates | excluded_site_only_overlapping | 141 |
| aggregates | excluded_ambiguous_overlapping | 49 |
| aggregates | excluded_duplicate_accession_overlapping | 0 |
| aggregates | excluded_union | 229 |
| aggregates | eligible_unique | 2691 |
| proteomes | raw_rows | 2553 |
| proteomes | excluded_reverse_overlapping | 98 |
| proteomes | excluded_contaminant_overlapping | 35 |
| proteomes | excluded_site_only_overlapping | 127 |
| proteomes | excluded_ambiguous_overlapping | 28 |
| proteomes | excluded_duplicate_accession_overlapping | 0 |
| proteomes | excluded_union | 192 |
| proteomes | eligible_unique | 2361 |

Exclusion categories overlap; excluded_union is the non-overlapping total removed. Eligible sets contain 2,691 aggregate and 2,361 proteome accessions, with 2,221 shared (470 aggregate-only and 140 proteome-only). Requiring complete positivity removes 732 shared accessions from the primary universe. Primary coverage is 55.33% of eligible aggregate accessions and 63.07% of eligible proteome accessions. This selects reliably detected proteins and cannot represent proteins lost or gained below detection.

Retained positive iBAQ signal, relative to all eligible positive measurements in that fraction/sample:

| Universe | Fraction | Minimum | Maximum |
| --- | --- | --- | --- |
| primary | aggregates | 89.83% | 93.75% |
| primary | proteomes | 91.50% | 93.05% |
| secondary | aggregates | 89.92% | 93.66% |
| secondary | proteomes | 91.12% | 92.52% |
| common | aggregates | 89.05% | 92.95% |
| common | proteomes | 90.45% | 92.07% |

Across all 24 ED columns per fraction, no eligible value was missing or nonfinite. Finite nonpositive values occurred in 215–433 eligible aggregate rows and 63–245 eligible proteome rows per sample. SAMPLE_QC.tsv gives each column; DENOMINATOR_QC.tsv gives the actual denominator and retained signal per sample. All normalized sums were one within 1e-12. Both actual archive SHA256 values match historical metadata. GATES.tsv records passing computational gates and the interpretation limitations separately. No post-outcome coverage cutoff was introduced.

## Label-split negative controls

For each treatment, the control is [Y(deletion,AB) − Y(deletion,CD)] − [Y(ED,AB) − Y(ED,CD)]. Both fraction averages are formed separately. These are diagnostics of label structure, not randomized null samples or biological replicate tests.

| analysis | n | q25 | median | q75 | positive_fraction |
| --- | --- | --- | --- | --- | --- |
| primary_AB | 1489 | -1.3444 | -0.4659 | 0.2685 | 0.3499 |
| primary_CD | 1489 | -1.2406 | -0.4425 | 0.2265 | 0.3277 |
| primary_AB_minus_CD | 1489 | -0.6035 | 0.0086 | 0.5620 | 0.5057 |
| negative_control_blank | 1489 | -0.3595 | 0.1172 | 0.5957 | 0.5588 |
| negative_control_Nva | 1489 | -0.2547 | 0.1008 | 0.4631 | 0.5897 |
| negative_control_Val | 1484 | -0.5253 | -0.0257 | 0.3852 | 0.4852 |

The blank and Nva split controls are not exactly centered at zero (medians +0.1172 and +0.1008); their broad distributions caution against treating protein signs as independent evidence. AB-only and CD-only primary medians remain negative. The difference of these primary split estimates has median +0.0086, but this does not certify absence of batch effects. No arbitrary pass/fail tolerance was imposed on the controls.

## Leave-one-label-out robustness

The primary protein universe stays fixed. Eight checks omit a label in only one fraction; four additional checks omit it in both fractions as a sensitivity operation, without asserting pairing.

| analysis | n | median | q25 | q75 | positive_fraction | sign_agreement_with_primary |
| --- | --- | --- | --- | --- | --- | --- |
| LOO_aggregates_A | 1489 | -0.3166 | -1.1443 | 0.3206 | 0.3761 | 0.9046 |
| LOO_aggregates_B | 1489 | -0.4737 | -1.2188 | 0.1767 | 0.3143 | 0.9355 |
| LOO_aggregates_C | 1489 | -0.4626 | -1.2061 | 0.1933 | 0.3049 | 0.9369 |
| LOO_aggregates_D | 1489 | -0.5029 | -1.3991 | 0.2078 | 0.3264 | 0.9382 |
| LOO_proteomes_A | 1489 | -0.4233 | -1.2573 | 0.2173 | 0.3170 | 0.9610 |
| LOO_proteomes_B | 1489 | -0.4817 | -1.2768 | 0.1672 | 0.3076 | 0.9543 |
| LOO_proteomes_C | 1489 | -0.4901 | -1.2668 | 0.1585 | 0.3049 | 0.9570 |
| LOO_proteomes_D | 1489 | -0.3213 | -1.0893 | 0.3348 | 0.3714 | 0.9268 |
| LOO_both_A | 1489 | -0.3265 | -1.1774 | 0.3207 | 0.3606 | 0.9120 |
| LOO_both_B | 1489 | -0.5253 | -1.2804 | 0.1001 | 0.2968 | 0.9167 |
| LOO_both_C | 1489 | -0.5109 | -1.2318 | 0.1431 | 0.2895 | 0.9161 |
| LOO_both_D | 1489 | -0.3996 | -1.2750 | 0.3499 | 0.3593 | 0.8986 |

Every leave-one-label-out median remains negative (range −0.5253 to −0.3166). Per-protein sign agreement is 89.86–96.10%; some individual signs are sensitive even though the median direction persists. These ranges are sensitivity summaries, not confidence intervals. Proteins are correlated response units and are never used as biological n for p-values.

## Independent evidence and prior art

[Ruan 2008](https://doi.org/10.1073/pnas.0809179105) supplies independent evidence of tolerance of misacylation-driven mistranslation and a heat-shock response. It supports physiological buffering, not a quantitative reserve or universal tolerable error rate. Its reporter and selection settings differ from this archive.

[Pranjic 2024](https://doi.org/10.1016/j.ijbiomac.2024.130068) reports proteostasis responses and ClpB-associated counteraction of aggregation in editing-defective IleRS cells challenged with Val or Nva. The archives here come from that same study, so publication claims and this reanalysis are not independent confirmations. A negative median of normalized redistributions does not contradict a bulk aggregation result measured on another scale. Full methods and culture/block linkage were not newly authenticated; header identities justify a bounded descriptive analysis only.

[Martens & Hilser 2025](https://doi.org/10.1101/2025.06.25.661590) is treated as the specified bioRxiv preprint. Synonymous luciferase constructs change protein production while preserving amino acid sequence, and the reported chaperone response is saturable. This supports a production-load mechanism; it does not manipulate or measure mistranslation dose, and it cannot calibrate mistranslation headroom.

[McDonald 2025](https://doi.org/10.1093/nar/gkaf428) distinguishes yeast mistranslation routes: R-to-S affects Hsp70 association, whereas P-to-S is associated with impaired handling/degradation of nascent proteins. Different substitutions and network branches cannot be collapsed into one interchangeable error burden or inferred from a bacterial ClpB contrast.

[Fauvet 2021](https://doi.org/10.3389/fmolb.2021.653073) reports 20.5% insoluble protein in delta dnaKJ and 18.9% in delta dnaKJ delta htpG at 30 C, a −1.6 percentage-point change. The existing local kinetic model instead predicts an increase after HtpG deletion in 200/200 one-client and 199/200 two-client configurations (see ../experiment_grounded_rebuild_20261010/kinetics/RESULTS.md). The model point-estimate direction is wrong for this comparison. The published aggregate totals alone do not establish uncertainty or significance of that difference; neither this caveat nor the current Pranjic proxy rescues the directional mismatch. No model was refit here.

[FoldEco (2012)](https://doi.org/10.1016/j.celrep.2012.02.011) already models protein synthesis, folding, chaperone systems, degradation and aggregation in E. coli. Finite chaperone networks and kinetic competition are prior art. This analysis claims no new discovery of those mechanisms. PROXY_LITERATURE.tsv records primary-source access scope; abstract-only access is not presented as a full methods verification.

## What remains supported and what remains unidentified

**Positive conditional mathematics:** finite resource accounting remains valid under explicit assumptions. If a machine pool has C available sites and each processed client requires mean occupied site-time tau, its long-run processing rate cannot exceed C/tau under that fixed service model. Increasing demand at fixed available sites and service requirement reduces the arithmetic reserve. This is a conditional constraint, not a measured biological law; client routing, sequestration, synthesis, degradation and adaptive responses can change those assumptions. This run does not re-prove all prior nonlinear model claims.

**Biological proxy evidence:** independent studies support proteostasis responses and context-dependent buffering. This reanalysis finds a negative median relative redistribution interaction, robust in median direction to the specified label omissions, with no positive Nva-specific location shift over Val. It does not provide the desired positive proxy confirmation and does not refute all buffering biology.

**Quantitative headroom remains unidentifiable:** these fractions provide neither actual mistranslation burden nor occupied chaperone site-time, available machine capacity, absolute aggregate mass or a calibrated mapping among them. Many absolute loads and capacities can yield the same compositions. There is no numerical headroom estimate, universal threshold or preregistration confirmation to report.

## Reproducibility and stopping

PROXY_ANALYSIS.py was executed on both actual archives. Eight tests passed, including independent scalar/csv reconstruction of every primary Y and interaction, exclusion logic, contrast sign, independent fraction averaging, sample-scale invariance, zero rejection, summary arithmetic and sensitivity/control identities. See PROXY_TEST_RESULTS.txt. SOURCE_HASHES.tsv pins inputs; SHA256SUMS pins delivered files. README.md gives commands. Pending ledger additions are in PROXY_LEDGER_PENDING.md because all paths outside this directory were explicitly read-only. The bounded analysis is complete; no further experiment or paper rewriting was undertaken.
