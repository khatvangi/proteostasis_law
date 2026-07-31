# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## (!) CANONICAL LOCATION — read this first

**All active work lives in `envelope-paper/`** (canonical since 2026-07-29).

- manuscript: `envelope-paper/manuscript/MANUSCRIPT.md`
- how to reproduce: `envelope-paper/README.md`

Two earlier directories are superseded and must not be cited or submitted:

| Directory | Status |
|---|---|
| `envelope-paper/` | **ACTIVE** |
| `combined-paper/` | superseded — see `combined-paper/SUPERSEDED.md` |
| `proteostasis-paper/` | superseded — earlier evidence-ledger scaffold |
| `legacy_root_triplet_attempt/` | archive only |
| `manscuript-final.md` (root) | stale; its N=17,166 and several odds ratios are not reproducible |

## active claim

Extant translation operates inside a finite proteostasis envelope. Decoding
fidelity, folding success, aggregation burden, and quality-control demand are
coupled, and viability requires that combined burden stay within finite buffering
and clearance capacity (`B_total <= C_buffer`). Two bounding arguments place
E. coli roughly ONE order of magnitude *inside* that envelope — x25 at the
usage-weighted evaluation point, x1.9 to x25 across chaperone availability theta.
"Two orders" holds only at the bottom of the observed error window, which is the
evaluation point the paper rejects. At the x25 margin the distinguishing
prediction is a synthetic lethality reachable in wild type (+7.4% above additive
for 3x error against 3x rescue knockdown; 12 of 36 pairs survivable alone and
lethal together). Separately, measured per-codon mistranslation is organized at
the amino-acid level (mu), and the supply axis (nu) shows no comparable structure
-- including when mu's OWN per-amino-acid pattern is imposed on nu as the
benchmark (z = -2.75, p = 0.0034). The detectability confound was tested, not
conceded: regressing log mu on log sampling depth removes the 18.8% of variance
depth explains and the clustering persists on the residuals (z = -2.97,
p = 0.0026, 80% of the between-aa eta^2 retained). Residualizing is licensed by
the NEGATIVE depth-mu correlation, which is thin-sampling inflation rather than
error-driven detection. Do not reinstate "cannot be separated from
detectability" -- that claim was retracted on evidence.

This is NOT a code-origin claim. It does not explain why codons are triplets
or the evolutionary origin of redundancy. Those claims are explicitly rejected.

Two claims from the earlier combined draft are **rejected on the evidence** and
must not be reintroduced without new data:

- metal-binding-site codon deployment — a gene-level expression confound; dies
  against a within-gene background (0 of 4 significant, His OR 1.07)
- cross-species conservation of operational geometry — artifact of sharing
  E. coli mu across all three species plus non-independent tAI vectors

## (!) the rule

**A number may appear in the manuscript only if a script in
`envelope-paper/scripts/` can recompute it from `envelope-paper/data/raw/`.**

The previous draft's numbers traced to markdown summaries rather than to code.
Its verification pass checked draft-against-TSV but never TSV-against-code, and
so certified a corrupt supply axis, a superseded Monte Carlo run, and a retracted
statistic while reporting "17/20 pass". Never verify against a serialized summary.

## commands

### reproduce everything (from `envelope-paper/`)
```
python scripts/run_all.py            # 16 steps, assembles the paper, then the tests
python scripts/run_all.py --fast     # skip the 10,000-permutation runs
```
Individual steps are listed in `envelope-paper/README.md`. Two gates fail closed:
`01_validate_tai.py` (invalid supply axis) and `15_build_paper.py` (unresolved
table placeholder or missing figure).

### the paper is assembled, not hand-maintained
`manuscript/MANUSCRIPT.md` holds the prose, the figure embeds, and one
`<!-- TABLE:Table N -->` placeholder per main table. `scripts/15_build_paper.py`
fills those from `tables/TABLES.md` and emits `manuscript/PAPER.{md,html,pdf,docx}`.
Never paste a table body into the manuscript — a test fails if a markdown table
row appears in the prose. Table numbering changed in the v2 restructuring: the
two excluded analyses are no longer Tables 4/S4, and the retired anchoring grid is
no longer Table 6/S9.

### run the numeric test suite (103 tests, asserts the manuscript against data/computed/)
```
cd envelope-paper && python -m unittest discover -s tests -v
```

### legacy suite for the superseded scaffold (15 structural checks, no numbers)
```
cd proteostasis-paper && python -m unittest discover -s tests -p 'test_*.py'
```

Dependencies: numpy, pandas, scipy, seaborn, openpyxl. No build step or linter.

## (!) corrupt upstream data

`codon-deployment/data/computed/{ecoli,bsubtilis,scerevisiae}_tai_ws.tsv` are
corrupt — 44/60 E. coli values are bit-identical to yeast, and CTG is assigned a
near-floor weight. Do not use them as nu. See
`codon-deployment/data/computed/README_TAI_IS_CORRUPT.md`.
`structural-criticality/` symlinks all three and hardcodes them in `config.py`,
so every tAI-axis result there needs recomputation.

Validated E. coli nu: `envelope-paper/data/computed/nu_tai_ecoli_validated.tsv`.
No validated B. subtilis or S. cerevisiae vector exists yet.

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
