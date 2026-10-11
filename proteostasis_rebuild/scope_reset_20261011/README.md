# START HERE — completed proteostasis reconstruction

**Finished scientific reconstruction; new wet-lab validation unavailable. Manuscript verdict: theoretical reassessment and experimental constraints / not submission-ready.** Everything for this delivery is together here. Earlier branches were not modified.

Read [COMPLETE_INTEGRATED_REPORT.md](COMPLETE_INTEGRATED_REPORT.md) first for the original question, physical reconstruction, evidence, rejected claims, competing predictions, blockers and publishing decision. It ends with the explicit STOP / FINISHED scope boundary.

| Material | Purpose |
|---|---|
| [REVISED_MANUSCRIPT.md](REVISED_MANUSCRIPT.md) | Complete ~3,200-word manuscript, canonical authors, negative results, qualified evidence and reference IDs |
| [MATHEMATICAL_CHECKS.py](MATHEMATICAL_CHECKS.py) | Executable conservation, incremental workload, units, finite-pool positivity, occupancy and nonautomatic-response checks |
| [CHECK_RESULTS.txt](CHECK_RESULTS.txt) | **20 run, 20 passed, 0 failed, 0 errors, 0 skipped**; Python 3.12.11, SymPy 1.14.0 |
| [CLAIM_AND_EVIDENCE_MATRIX.tsv](CLAIM_AND_EVIDENCE_MATRIX.tsv) | 22 manuscript claim IDs, 34 evidence rows, exact source paths/lines, named literature and statuses |
| [ADVERSARIAL_AUDIT.md](ADVERSARIAL_AUDIT.md) | Actual limitations, original-scope checks, publication boundaries and repaired discrepancies |
| [ORIGINAL_FROZEN_MANUSCRIPT.md](ORIGINAL_FROZEN_MANUSCRIPT.md) | Byte-identical historical original; **contains withdrawn claims and is not the current scientific conclusion** |
| [SOURCE_MANIFEST.tsv](SOURCE_MANIFEST.tsv) | 27 source files with SHA-256 hashes; local sources remain in their original locations |
| [EVIDENCE_RECONCILIATION.md](EVIDENCE_RECONCILIATION.md) | Later Fauvet supplement inspection, precise uncertainty blocker and FoldEcoSlim prior-art correction |
| [ORIGIN_AND_SCOPE.md](ORIGIN_AND_SCOPE.md) | Earlier reset plus dated completion clarification |
| [MINIMAL_PHYSICAL_MODEL.md](MINIMAL_PHYSICAL_MODEL.md) | Detailed minimal accounting on which the manuscript is based |
| [DECISION_GATES.md](DECISION_GATES.md) | Prospective measurement gates, not completed experiments; source-status clarification included |
| [LEDGER_PENDING.md](LEDGER_PENDING.md) | Append-ready decision record; external ledger untouched under the write-only-here instruction |
| [VERIFY_DELIVERY.py](VERIFY_DELIVERY.py) | Checks source hashes, matrix anchors, references, authors, links and package consistency |
| [DELIVERY_VERIFICATION.txt](DELIVERY_VERIFICATION.txt) | Recorded package verification; distinct from mathematical and biological tests |

The canonical source remains [the original manuscript](../../envelope-paper/manuscript/MANUSCRIPT.md), read-only. Its root [historical README](../../README.md) and older reports may contain superseded statements. The report and revised manuscript here give the current verdict.

## Exact reproduction commands

From any directory:

```bash
cd /storage/kiran-stuff/proteostasis_law/proteostasis_rebuild/scope_reset_20261011
python3 MATHEMATICAL_CHECKS.py > CHECK_RESULTS.txt 2>&1
python3 VERIFY_DELIVERY.py > DELIVERY_VERIFICATION.txt 2>&1
```

The mathematical script uses the standard library and uses SymPy if installed. Without SymPy the two symbolic tests are explicitly skipped; exact rational checks remain. It does not fit biology, rerun old ensembles or write to earlier branches. The verification script reads original evidence sources; moving the whole package without those sources preserves the narrative and mathematical checks but not live provenance verification.

The remaining empirical gaps are the matched error-to-workload measurement, in-vivo arrivals/residence and functional throughput, and Fauvet total-fraction uncertainty/covariance with matched 30 °C rates. Some original client/interaction tables and the V04 duplicate-publication provenance also remain unresolved. Mathematical checks do not close these gaps. **STOP: the requested reconstruction is finished; no new theory branch or unperformed validation is implied.**
