# Bounded exploratory reanalysis — 2026-10-12

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
