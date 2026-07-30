# A finite proteostasis envelope for translation

Active manuscript: [`manuscript/MANUSCRIPT.md`](manuscript/MANUSCRIPT.md)

This is the canonical home of the proteostasis-law project as of 2026-07-29. It
supersedes `../combined-paper/` and `../proteostasis-paper/`.

## The claim

Extant translation operates inside a finite proteostasis envelope: viability
requires `B_total ≤ C_buffer`, where burden decomposes into decoding error,
folding failure, aggregation, and quality-control demand. A reduced dynamical
model with saturable rescue and superlinear aggregation feedback necessarily has
a saddle-node threshold. Two bounding arguments place E. coli roughly two orders
of magnitude *inside* its envelope. Separately, the genetic code is strongly
organized on the error axis and shows no structure on the supply axis.

This is **not** a code-origin claim. It says nothing about why codons are
triplets or why degeneracy exists.

## Reproducing everything

Run the scripts in numeric order from this directory. Nothing depends on
network access; all inputs are in `data/raw/`.

```bash
python scripts/01_validate_tai.py                      # ~10 s   validates the nu axis, aborts if it fails
python scripts/02_axis_structure.py                    # ~2 min  permutation tests (10,000 x 3 axes x 2 nulls)
python scripts/02_axis_structure.py --mu-stat median    # ~2 min  sensitivity to the mu summary statistic
python scripts/06_translation_burden.py                # ~5 s    burden magnitude from mu + codon usage
python scripts/07_removed_results.py                   # ~30 s   diagnostics for the two dropped results
cd scripts && python 03_fig1_envelope.py && python 04_fig2_axis.py && python 05_fig3_bounds.py
```

Then verify:

```bash
python -m unittest discover -s tests -v     # 31 tests
```

Requires Python 3 with numpy, pandas, scipy, seaborn, openpyxl. No build step.

## Layout

```
manuscript/MANUSCRIPT.md   the paper
scripts/                   every number and figure is generated here
  01_validate_tai.py       validates the tAI (supply) axis; ABORTS if it fails
  02_axis_structure.py     Delta_A, permutation nulls, variance decomposition
  03_fig1_envelope.py      Fig 1  ODE regimes and the saddle-node
  04_fig2_axis.py          Fig 2  mu clustering, nu null, variance partition
  05_fig3_bounds.py        Fig 3  paired bounds and headroom
  06_translation_burden.py burden magnitude implied by the mu data
  07_removed_results.py    reproduces the two results this paper drops
  figstyle.py              shared seaborn style
data/raw/                  staged inputs, self-contained
data/computed/             script outputs; the manuscript's only source of numbers
figures/                   PNG + PDF
tests/                     asserts the manuscript against data/computed/
```

## The rule this project now enforces

**A number may appear in the manuscript only if a script in `scripts/` can
recompute it from `data/raw/`.**

This exists because of how the previous version failed. Its numbers traced to a
markdown summary, which traced to another markdown summary. Its verification pass
compared the draft against precomputed TSVs but never compared those TSVs against
generating code — so it certified a corrupt supply axis, a superseded Monte Carlo
run, and a statistic the upstream project had explicitly retracted, and still
reported "17/20 pass". `tests/test_manuscript_numbers.py` closes that gap, and
the suite has been mutation-tested: altering a manuscript number, dropping a
citation, softening a negative result, reinstating the retracted statistic, or
substituting the corrupt tAI vector each make it fail.

## What was removed, and why

Two empirical results from the earlier combined draft are **excluded**, with the
corrected analyses reported in the Discussion and reproduced by
`scripts/07_removed_results.py`:

| Removed result | Why |
|---|---|
| Metal-binding-site codon enrichment (OR 1.28–1.40, p < 0.03) | Gene-level expression confound. Against non-metal positions of the *same* genes, 0 of 4 tests remain significant and His falls from OR 1.33 to 1.07. Separately, ~40% of annotated sites are mispositioned in sequence space. |
| Cross-species conservation of operational geometry (ρ = 0.63–0.92) | Two artifacts. Δ_A used E. coli μ for all three species (reproduced exactly), and 44 of 60 E. coli tAI values were bit-identical to the yeast ones. |

## Known gaps before submission

- Public deposition DOI is not yet assigned (`Data and code availability`).
- The framework's distinguishing prediction (supraadditivity) is untested. This
  is stated in the manuscript, not hidden.
- ν is validated for E. coli only. No validated B. subtilis or S. cerevisiae
  supply vector exists; see `../../codon-deployment/data/computed/README_TAI_IS_CORRUPT.md`.
