# Deposition contents — a finite proteostasis envelope for translation

Rebbeck, H. E., Paudyal, J. & Boggavarapu, K.
Department of Chemistry and Physics, McNeese State University, Lake Charles, LA 70609, USA

This archive is the complete computational record for the manuscript: every raw
input, every script, every computed output, the figures, the tables, the assembled
paper, and the test suite that asserts the manuscript against the computed data.

## Provenance

| | |
|---|---|
| Source repository | https://github.com/khatvangi/proteostasis_law |
| Path within it | `envelope-paper/` |
| Commit archived | see `COMMIT` in this directory |
| Built | see `BUILT` in this directory |
| Integrity | `SHA256SUMS` covers every file below |

## Reproducing everything from this archive

No network access is needed; all inputs are in `data/raw/`.

```bash
python scripts/run_all.py          # 16 steps (~12 min), then the test suite
python scripts/run_all.py --fast   # skips the 10,000-permutation runs (~2 min)
python -m unittest discover -s tests -v
```

Requires Python 3 with numpy, pandas, scipy, seaborn, openpyxl. `PAPER.html` and
`PAPER.docx` additionally need pandoc, and `PAPER.pdf` is printed from the HTML by
headless Chrome; all three are skipped with a message if the tool is absent, and
the markdown is always written.

Two gates fail closed. `01_validate_tai.py` exits non-zero if the translational
supply axis fails its validation test, and nothing downstream runs.
`15_build_paper.py` exits non-zero on an unresolved table placeholder or a missing
figure.

## Contents

```
manuscript/
  MANUSCRIPT.md            the paper: prose, figure embeds, table placeholders
  PAPER.md / .html / .pdf / .docx   assembled paper, tables inlined
  manuscript_v2_draft.md   the drafting source, kept verbatim
data/raw/                  staged inputs, self-contained (see below)
data/computed/             every number the manuscript states
scripts/                   01-15 plus run_all.py, figstyle.py, vendor/
tables/                    Tables 1-5 and S1-S7 as TSV, plus TABLES.md
figures/                   Fig 1-4 as PNG, PDF and SVG
tests/                     103 assertions of the manuscript against data/computed/
README.md                  full documentation, including what was excluded and why
```

### Raw inputs and their sources

| File | Source |
|---|---|
| `Data_S2_error_detection_rate.xlsx` | Landerer, Poehls & Toth-Petroczy, *Mol. Biol. Evol.* 41 (2024), PMC10939442 — supplementary Data S2 |
| `codon_error_rates_ecoli.tsv` | per-codon mistranslation rates extracted from the above |
| `ecoli_tai_reference.tsv` | tAI for *E. coli* K-12 (dos Reis et al. 2004 method; GtRNAdb gene copy numbers) |
| `ecoli_tai_CORRUPT_for_regression_test.tsv` | a known-bad supply vector, retained so the test suite can assert that it still fails validation |
| `ecoli_k12_cds.fna` | *E. coli* K-12 MG1655 coding sequences |
| `global_codon_usage_ecoli.tsv` | genome-wide codon usage |
| `metal_sites_ecoli_with_codons.csv` | annotated metal-binding sites (used only by the excluded analysis) |
| `arithmetic_results.json`, `two_pool_results.json`, `paired_mc_results.json` | staged outputs of the upstream `proteostasis-P1` bound calculations |
| `bsubtilis_tai_AS_USED_previously.tsv`, `scerevisiae_tai_AS_USED_previously.tsv`, `delta_a_by_species_AS_PUBLISHED.tsv`, `prior_tai_script_for_tRNA_counts.py.txt` | inputs to the excluded cross-species analysis, retained so its failure is reproducible |

## What this archive deliberately contains that the paper does not

`scripts/07_removed_results.py` reproduces two analyses that were **withdrawn on
verification** and are not results of the paper: metal-binding-site codon
enrichment (a gene-level expression confound; 0 of 4 tests survive an
expression-matched background) and cross-species conservation of operational
geometry (two artifacts, including 44 of 60 supply values being bit-identical
between two species). Their outputs are named `tables/Excluded_*.tsv` so they
cannot be mistaken for tables of the paper. `README.md` documents both in full,
under *Archived analyses*. They are kept because deleting them would hide a
correction.

## Verification properties of this archive

- **Every number in the manuscript is recomputable** from `data/raw/` by a script
  in `scripts/`. The test suite asserts the manuscript text against the computed
  artifacts, not against any intermediate summary.
- **A rebuild is byte-reproducible.** Re-running the pipeline reproduces every
  computed output, figure and table exactly. Figure jitter is seeded and build
  timestamps are stripped, so a rebuild-and-compare check has no exceptions to
  remember. `PAPER.pdf` and `PAPER.docx` are the exception in kind, not in
  content: Chrome and pandoc stamp a creation date into the container.
- **Consistency between sections is asserted, not just provenance per number.**
  A document whose sections are each individually sourced can still contradict
  itself; the suite compares the artifacts that must agree.

## Known limitations, stated plainly

- The framework's distinguishing prediction — supraadditive interaction between
  burden stages — is **untested experimentally**. It is confirmed in the model.
- Chaperone availability θ is **not measured**, so the headroom is a range
  (×1.9 to ×25) rather than a value.
- All per-codon mistranslation data are from *E. coli*; that half of the analysis
  is single-species.
- The reduced ODE is illustrative and not fitted; it establishes that a collapse
  threshold exists, not where any organism sits relative to it.

## Licence

**A licence must be selected at upload and is not asserted here.** Nothing in this
archive states one, and the authors have not chosen one yet. The conventional pair
for a deposition of this kind is CC-BY-4.0 for the manuscript, figures and data,
and a permissive code licence (MIT or BSD-3-Clause) for `scripts/` and `tests/`.
`scripts/vendor/two_pool_ode.py` is a verbatim copy of upstream project code and
carries whatever terms that project carries.
