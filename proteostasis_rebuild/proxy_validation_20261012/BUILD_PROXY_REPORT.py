#!/usr/bin/env python3
"""Render the bounded reanalysis report from delivered result tables."""
from pathlib import Path
import csv, hashlib
P=Path(__file__).resolve().parent
def rows(name):
    with (P/name).open() as h:return list(csv.DictReader(h,delimiter='\t'))
def table(records, columns):
    text='| '+' | '.join(columns)+' |\n| '+' | '.join(['---']*len(columns))+' |\n'
    for r in records:
        vals=[]
        for c in columns:
            v=r[c]
            if c not in ['analysis','n','universe','fraction','stage']:
                try:v=f'{float(v):.4f}'
                except ValueError:pass
            vals.append(v)
        text+='| '+' | '.join(vals)+' |\n'
    return text
s=rows('DISTRIBUTION_SUMMARY.tsv');cov=rows('COVERAGE.tsv');qc=rows('FILTER_QC.tsv');den=rows('DENOMINATOR_QC.tsv')
lit=[
('Ruan 2008','10.1073/pnas.0809179105','https://scholarship.miami.edu/esploro/outputs/journalArticle/Quality-control-despite-mistranslation-caused-by/991031599940802976','Author institution abstract; PMC full text browser challenge','Independent mistranslation tolerance and heat-shock response context; no headroom calibration.'),
('Pranjic 2024','10.1016/j.ijbiomac.2024.130068','https://www.croris.hr/crosbi/publikacija/prilog-casopis/821559','Author institution abstract and original local archived proteinGroups bytes','Same study as the reanalysis; not independent replication; editing-defective IleRS and ClpB context.'),
('Martens & Hilser 2025 preprint','10.1101/2025.06.25.661590','https://pmc.ncbi.nlm.nih.gov/articles/PMC12262517/','Indexed primary article abstract/results/methods; direct open challenged','Synonymous reporter variants alter production with unchanged amino acid sequence; saturable response, not mistranslation dose.'),
('McDonald 2025','10.1093/nar/gkaf428','https://academic.oup.com/nar/article/53/9/gkaf428/8133635','Primary publisher article introduction/results/discussion','Yeast P-to-S and R-to-S challenge distinct branches; not a scalar interchangeable error load.'),
('Fauvet 2021','10.3389/fmolb.2021.653073','https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2021.653073/full','Primary full text, results around Figure 4','20.5% versus 18.9% insoluble in delta dnaKJ versus delta dnaKJ delta htpG at 30 C; opposite local model point-estimate direction.'),
('FoldEco 2012','10.1016/j.celrep.2012.02.011','https://www.sciencedirect.com/science/article/pii/S221112471200068X','Primary publisher summary','Prior art for synthesis/folding/chaperones/degradation/aggregation network modelling.')]
with (P/'PROXY_LITERATURE.tsv').open('w') as h:
    w=csv.writer(h,delimiter='\t');w.writerow(['study','doi','source_url','access_scope','interpretation']);w.writerows(lit)
report='''# Final exploratory proxy report

The specified proxy does **not** show a positive median Nva × ClpB-deletion interaction. Across 1,489 matched complete-positive proteins, the median is **−0.4145 log2 units**, quartiles **−1.2125 and +0.2294**, and **32.50%** are positive. This is a descriptive relative-redistribution result, not evidence that deletion lowers total aggregation. The Nva-versus-Val comparison provides no clear positive Nva-specific location shift under this estimand.

## Status and scope

PROTOCOL.md was written and hashed before this run inspected raw outcome headers or values. The folder already existed with a literature dossier; this run adds an executable archival analysis and preserves those files, with the prior README saved as README_PREEXISTING.md. Existing knowledge and historical analyses make this exploratory, not genuine preregistration confirmation. PROTOCOL_FREEZE.json and OUTCOME_READ_RECEIPT.json document sequence, not third-party attestation. The historical GATE3_STATUS.md remains VOID and unmodified. All writes were confined to this directory; no papers were rewritten and no subagents were used.

## Estimand

Rows flagged reverse, contaminant or site-only were excluded, as were empty/multiple accessions, groups with Number of proteins other than one, and duplicate accessions. Exact accessions were joined across fractions. Primary inclusion required positive finite iBAQ in all 32 measurements: four labels × four genotype/treatment arms × two fractions. Zeros were not imputed.

For the fixed matched set, each sample was divided by the summed iBAQ of that same set in its own fraction. Shares were averaged separately within condition for each fraction; Y = log2(mean insoluble share / mean proteome share). The interaction is [Y(ED Delta clpB,Nva) − Y(ED Delta clpB,blank)] − [Y(ED,Nva) − Y(ED,blank)]. There is no A-D pairing assumption across fractions. Positive values represent a larger Nva-related *relative* insoluble enrichment in the deletion background.

This normalization discards fraction-wide scale: multiplying every aggregate intensity in any sample by any positive factor leaves its shares unchanged. Thus an increase of every protein's aggregate mass could be invisible. Conversely, a few abundant proteins changing can move other proteins' shares. iBAQ shares are not calibrated aggregate mass fractions, and pellet contents can include constitutively insoluble material.

## Primary and predefined specificity results

'''
report+=table(s[:5],['analysis','n','q25','median','q75','positive_fraction','minimum','maximum'])
report+='''
Val uses its own complete-positive set. The common set requires positivity in all 48 measurements, allowing a like-for-like Nva-minus-Val comparison. Its median difference is −0.0252, with 48.57% positive; absence of a large descriptive shift is not an equivalence test. The full-primary minimum (−14.4044) disappears in the more restrictive common set, illustrating sensitivity of individual extremes to coverage; no trimming or outlier exclusion was performed.

## Coverage and QC

'''+table(cov,['universe','n','eligible_shared','fraction_shared','aggregate_complete','proteome_complete','required_measurements_per_protein'])+'\n'+table(qc,['fraction','stage','n'])
report+='''
Exclusion categories overlap; excluded_union is the non-overlapping total removed. Eligible sets contain 2,691 aggregate and 2,361 proteome accessions, with 2,221 shared (470 aggregate-only and 140 proteome-only). Requiring complete positivity removes 732 shared accessions from the primary universe. Primary coverage is 55.33% of eligible aggregate accessions and 63.07% of eligible proteome accessions. This selects reliably detected proteins and cannot represent proteins lost or gained below detection.

Retained positive iBAQ signal, relative to all eligible positive measurements in that fraction/sample:

| Universe | Fraction | Minimum | Maximum |
| --- | --- | --- | --- |
'''
for u in ['primary','secondary','common']:
    for f in ['aggregates','proteomes']:
        v=[float(r['retained_signal_fraction']) for r in den if r['universe']==u and r['fraction']==f]
        report+=f'| {u} | {f} | {min(v):.2%} | {max(v):.2%} |\n'
report+='''
Across all 24 ED columns per fraction, no eligible value was missing or nonfinite. Finite nonpositive values occurred in 215–433 eligible aggregate rows and 63–245 eligible proteome rows per sample. SAMPLE_QC.tsv gives each column; DENOMINATOR_QC.tsv gives the actual denominator and retained signal per sample. All normalized sums were one within 1e-12. Both actual archive SHA256 values match historical metadata. GATES.tsv records passing computational gates and the interpretation limitations separately. No post-outcome coverage cutoff was introduced.

## Label-split negative controls

For each treatment, the control is [Y(deletion,AB) − Y(deletion,CD)] − [Y(ED,AB) − Y(ED,CD)]. Both fraction averages are formed separately. These are diagnostics of label structure, not randomized null samples or biological replicate tests.

'''+table([r for r in s if r['analysis'].startswith('negative_control') or r['analysis'].startswith('primary_A') or r['analysis'].startswith('primary_C')],['analysis','n','q25','median','q75','positive_fraction'])
report+='''
The blank and Nva split controls are not exactly centered at zero (medians +0.1172 and +0.1008); their broad distributions caution against treating protein signs as independent evidence. AB-only and CD-only primary medians remain negative. The difference of these primary split estimates has median +0.0086, but this does not certify absence of batch effects. No arbitrary pass/fail tolerance was imposed on the controls.

## Leave-one-label-out robustness

The primary protein universe stays fixed. Eight checks omit a label in only one fraction; four additional checks omit it in both fractions as a sensitivity operation, without asserting pairing.

'''+table([r for r in s if r['analysis'].startswith('LOO_')],['analysis','n','median','q25','q75','positive_fraction','sign_agreement_with_primary'])
report+='''
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
'''
(P/'FINAL_PROXY_REPORT.md').write_text(report)
old=P/'README.md';backup=P/'README_PREEXISTING.md'
if old.exists() and not backup.exists():backup.write_bytes(old.read_bytes())
(P/'README.md').write_text('''# Bounded exploratory reanalysis — 2026-10-12

Start with [FINAL_PROXY_REPORT.md](FINAL_PROXY_REPORT.md). This is the new archival reanalysis, not preregistration confirmation. The historical VOID gate was not changed.

Primary result: 1,489 proteins; median interaction −0.4145, Q1 −1.2125, Q3 +0.2294 log2 units; 32.50% positive. Median direction survives all specified omissions. This is relative redistribution, not bulk aggregate mass or quantitative headroom.

Reproduce from any directory (Python, NumPy and pandas; exact versions in RUN_ENVIRONMENT.json):

```bash
OPENBLAS_NUM_THREADS=1 python /storage/kiran-stuff/proteostasis_law/proteostasis_rebuild/proxy_validation_20261012/PROXY_ANALYSIS.py
OPENBLAS_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 python /storage/kiran-stuff/proteostasis_law/proteostasis_rebuild/proxy_validation_20261012/TEST_PROXY_ANALYSIS.py
python /storage/kiran-stuff/proteostasis_law/proteostasis_rebuild/proxy_validation_20261012/BUILD_PROXY_REPORT.py
```

Scripts read fixed archival paths and write only next to themselves. No download is needed. Rerunning changes RUN_ENVIRONMENT.json's timestamp, so the delivered manifest then requires regeneration. A hash/schema/empty-set failure stops analysis and writes ANALYSIS_BLOCKER.txt; do not treat old result files as valid after a failed rerun. GATES.tsv is overwritten each run and records the current gate state.

- PROTOCOL.md, PROTOCOL_FREEZE.json, OUTCOME_READ_RECEIPT.json: frozen choices and run sequence.
- PROXY_ANALYSIS.py, ANALYSIS_RUN.txt: executed analysis and console summary.
- primary/secondary/common_PROTEIN_RESULTS.tsv: Y values and interactions.
- *_shares.tsv: per-sample compositions for each fixed universe.
- DISTRIBUTION_SUMMARY.tsv, ROBUSTNESS_PROTEINS.tsv, NEGATIVE_CONTROL_PROTEINS.tsv: descriptive results and controls.
- FILTER_QC.tsv, EXCLUSIONS.tsv, SAMPLE_QC.tsv, COVERAGE.tsv, DENOMINATOR_QC.tsv: exclusions, zeros, coverage and denominators.
- GATES.tsv: computational passes and inferential limitations.
- SOURCE_HASHES.tsv, SHA256SUMS: input and delivery SHA256 hashes.
- TEST_PROXY_ANALYSIS.py, PROXY_TEST_RESULTS.txt: eight passing tests, including independent raw reconstruction.
- PROXY_LITERATURE.tsv: sources and access limits; PROXY_LEDGER_PENDING.md: append-ready ledger entries, not centrally appended.

Verify delivery from this directory with `sha256sum -c SHA256SUMS`. The manifest covers regular files in this directory except itself and its verification log; it includes preserved earlier dossier files for inventory, not as new results. Earlier files such as PREREGISTERED_PROXY_RELATION.md and ONE_INTEGRATED_REPORT.md predate this run and do not establish preregistration for it. README_PREEXISTING.md preserves the earlier landing page. The new report governs interpretation of this reanalysis.
''')
(P/'PROXY_LEDGER_PENDING.md').write_text('''# Pending append — outside paths were read-only

The central ledger was not modified. These entries are prepared in its required format:

2026-10-12 | codex | RESULT | Bounded Pranjic archival compositional reanalysis gives Nva by ClpB-deletion interaction median -0.4145 log2 across 1489 complete-positive proteins with 32.50% positive and negative median in all label omissions; common-set Nva-minus-Val median -0.0252; not bulk aggregate mass or quantitative headroom | proteostasis_law/proteostasis_rebuild/proxy_validation_20261012/FINAL_PROXY_REPORT.md | Proteostasis envelope and fold condition
2026-10-12 | codex | DECISION | Frozen exploratory proxy protocol and executed raw-archive analysis do not restore the historical VOID gate or constitute genuine preregistration confirmation; FoldEco prior art and opposite-signed Fauvet model point estimate remain explicit | proteostasis_law/proteostasis_rebuild/proxy_validation_20261012/PROTOCOL.md | Proteostasis envelope and fold condition
''')
