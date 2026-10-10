# Evidence artifacts

Start with [EVIDENCE_REPORT.md](EVIDENCE_REPORT.md). `evidence.tsv` records source facts separately from interpretations. Empty values mean unavailable, never zero. Confidence concerns support for the bounded claim, not universal biological truth.

- `source_manifest.tsv`: DOIs, article/data identifiers, exact URLs, table/figure locations and retrieval state.
- `retrieval_attempts.tsv/json`: actual curl outcomes; guessed Elsevier filenames are explicitly candidates.
- `machine_deletion_essentiality.tsv`: all 15 machines, alleles, condition limits and unresolved single deletions.
- `machine_abundance.tsv`: existing Schmidt workbook-derived abundance, not throughput.
- `published_plating_table.tsv`: manual transcription of primary Table2 as returned by web search; one ND; superscript b means mucoid colonies. Incubation was 24 h at30–45°C and72 h at20°C.
- `plating_comparison.tsv`: arithmetic contrasts with no invented replicates or p-values.
- `model_vs_experiment.tsv`: existing model output beside unmatched experimental evidence; no empirical accuracy claim.
- `local_data_inventory.tsv`, `original_files_baseline.json`: provenance and preservation baseline.
- `sources/*.json`: saved web-tool responses, including errors/indexed text. They are not original PDF/XLS downloads.
- `checksums.sha256`: checksums of saved artifacts, excluding itself and the integrity report.

Rebuild derived files from this directory:

```bash
PYTHONDONTWRITEBYTECODE=1 python scripts/build_evidence.py
PYTHONDONTWRITEBYTECODE=1 python scripts/check_integrity.py
```

`scripts/retrieve.py` retries the recorded public targets and refuses to overwrite existing source files. Original source authentication, supplementary row reconciliation, clean single-deletion verification for three genes, and an independent Babu test remain outstanding. Do not run a client-overlap analysis until the original tables are available and reconciled.
