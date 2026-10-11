# Start here: closed experimental decision, 2026-10-11

**design complete; no new wet-lab performed**

**no confirmatory new test from currently identified public processed data**

Read [EXPERIMENTAL_DECISION_REPORT.md](EXPERIMENTAL_DECISION_REPORT.md), then the single [ONE_EXPERIMENT_PROTOCOL.md](ONE_EXPERIMENT_PROTOCOL.md). The selected experiment is MG1655 IleRS editing-deficient versus repaired, NONE/Nva/Val, crossed with verified inducible additional ClpB and matched controls. True occupied site-time is not identifiable; a validated pulse-cohort DnaK-association integral is the explicitly bounded primary observable.

- [SYSTEM_SELECTION.md](SYSTEM_SELECTION.md), [IDENTIFIABILITY_MATRIX.tsv](IDENTIFIABILITY_MATRIX.tsv) and [STAGE1_STATUS.md](STAGE1_STATUS.md): public-data decision and prior contact.
- [ANALYSIS_PLAN.md](ANALYSIS_PLAN.md), [SAMPLE_SHEET_TEMPLATE.tsv](SAMPLE_SHEET_TEMPLATE.tsv), [FREEZE_RECORD_TEMPLATE.tsv](FREEZE_RECORD_TEMPLATE.tsv), [GATE_SPEC.tsv](GATE_SPEC.tsv): pilot-to-lock plan, metadata and exact decisions. Templates are not measurements; calibration and n are deliberately unknown.
- [SOURCE_ANCHORS.tsv](SOURCE_ANCHORS.tsv), [METHODS_SOURCE_AUDIT.md](METHODS_SOURCE_AUDIT.md): DOI/URL, local lines and source-access limits.
- [SELF_AUDIT.md](SELF_AUDIT.md), [CHECK_PROTOCOL_RESULTS.txt](CHECK_PROTOCOL_RESULTS.txt): scope, limitations and offline verification.

Run from this directory: `python CHECK_PROTOCOL.py > CHECK_PROTOCOL_RESULTS.txt`. The checker validates structure, factor balance, references, gates, preservation and absence of invented completion fields. It performs no model fitting or simulation and does not certify experimental feasibility. Confirmation is correctly **NOT READY** until pilot qualification and freeze.

Preexisting Stage1 inventories and retrieval logs remain unchanged. PREEXISTING_HASHES.tsv snapshots existing workspace files and previously inventoried source files before this pass's writes; ORIGINAL_SOURCE_HASHES.tsv remains the earlier source snapshot. All task additions are inside this folder. Central ledger writing is outside the permitted scope; LEDGER_PENDING.md contains an append-ready entry, not a claimed update.

Completion clarification: total activity and pulse-labeled abundance are separate measurements; pulse-cohort functional yield is UNIDENTIFIABLE without independent cohort activity or validated conditional synthesis-stop preserving physiology and specific-activity calibration. G6 cannot claim new-native yield without that validation. **no wet lab, confirmation not ready.**
