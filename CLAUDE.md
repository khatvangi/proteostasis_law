# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## active claim

Proteostasis law is an operational constraint on extant translation systems.
Translation output, decoding fidelity, folding success, aggregation burden, and
quality-control demand are coupled. Viability requires that combined burden stay
within finite buffering and clearance capacity (`B_total <= C_buffer`).

This is NOT a code-origin claim. It does not explain why codons are triplets
or the evolutionary origin of redundancy. Those claims are explicitly rejected.

## commands

### run all tests (15 consistency checks)
```
cd proteostasis-paper && python -m unittest discover -s tests -p 'test_*.py'
```

### run a single test
```
cd proteostasis-paper && python -m unittest tests.test_active_project.ActiveProjectConsistencyTests.test_name
```

No build step, linter, or dependencies beyond Python stdlib (csv, re, unittest, pathlib).

## repo layout

```
proteostasis-paper/          # <-- ALL ACTIVE WORK LIVES HERE
  manuscript/
    MANUSCRIPT.md            # canonical manuscript (PNAS-style, markdown)
    PROTEOSTASIS_LAW.md      # compact formal statement with burden terms
    PHYSICS_FRAMEWORK.md     # reduced dynamical model (dP/dt equation)
    ACTIVE_DRAFT.md          # working draft referencing figures
  evidence/
    CLAIM_MATRIX.md          # 10 claims (C1-C10) with support status
    EVIDENCE_LEDGER.md       # per-source extraction tracker (10 sources)
    SOURCE_EXTRACTION_TEMPLATE.md
    extractions/             # one .md per source paper (structured notes)
  data/
    processed/source_index.tsv  # one-row-per-dataset canonical index
    raw/landerer_2024/       # only staged raw data so far
    DATA_MANIFEST.md         # rules for what counts as active data
    SOURCE_INVENTORY.md      # source-level inventory by burden term
  figures/
    stability_regions.svg    # Fig. 1 — burden-capacity regime schematic
    codon_load_allocation.svg # Fig. 2 — synonymous codons as load choices
  tests/
    test_active_project.py   # consistency suite (see below)
  legacy/                    # archived triplet-era analyses
  CLAUDE.md                  # paper-level project instructions
  CRITICAL_EVALUATION.md     # honest assessment of claim strength

legacy_root_triplet_attempt/ # old triplet-origin figures, scripts, data
                             # NOT active evidence — treat as archive only
```

## evidence-first workflow

The project follows a strict order: data → evidence → paper.
The manuscript is downstream of evidence, not the driver.

Key files to check before modifying claims:
1. `evidence/CLAIM_MATRIX.md` — which claims are supported/partial/rejected
2. `evidence/EVIDENCE_LEDGER.md` — which sources have been extracted
3. `data/processed/source_index.tsv` — canonical dataset index

A claim can enter the manuscript only if it meets the promotion rule in
CLAIM_MATRIX.md (direct evidence, explicit organism, no triplet-origin logic).

## what the tests enforce

`test_active_project.py` checks internal consistency across the evidence system:
- source_index.tsv has unique IDs and all required fields populated
- ledger IDs match source_index IDs exactly
- "done" ledger entries have matching extraction files in evidence/extractions/
- extraction files contain all 7 required sections
- claim matrix has exactly C1-C10
- supported claims appear in the ledger
- PROTEOSTASIS_LAW.md includes all burden terms and the core inequality
- manuscript explicitly rejects triplet-origin claims
- active figures and raw data files exist
- ACTIVE_DRAFT.md and MANUSCRIPT.md reference figures correctly

If you add a new source, you must update source_index.tsv, EVIDENCE_LEDGER.md,
and create an extraction file — the tests will catch mismatches.

## burden terms (the formal vocabulary)

| Term | Meaning |
|------|---------|
| `B_error` | decoding ambiguity / mistranslation burden |
| `B_fold` | folding-failure burden |
| `B_agg` | aggregation / toxic load |
| `B_qc` | rescue and clearance demand |
| `C_buffer` | effective proteostasis network capacity |
| `B_total` | sum of burden terms; viability requires `B_total <= C_buffer` |

## legacy material

`legacy_root_triplet_attempt/` contains the old triplet-origin project:
Python figure scripts, metal-site analyses, Pareto tests, cross-species data.
These are archived computational assets. Do not treat them as evidence for the
current operational-constraint claim.
