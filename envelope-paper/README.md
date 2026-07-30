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

One command rebuilds every number, figure, and table, then runs the tests.
Nothing needs network access; all inputs are in `data/raw/`.

```bash
python scripts/run_all.py            # ~3 min end to end, then 60 tests
python scripts/run_all.py --fast     # skip the two 10,000-permutation runs
python scripts/run_all.py --no-tests # rebuild only
```

The pipeline fails closed: if `01_validate_tai.py` rejects the supply axis it
exits non-zero and nothing downstream runs.

Individual steps, if you want them:

```bash
python scripts/01_validate_tai.py                     # ~3 s    validates nu, aborts if it fails
python scripts/02_axis_structure.py                   # ~65 s   permutation tests
python scripts/02_axis_structure.py --mu-stat median  # ~65 s   mu summary-statistic sensitivity
python scripts/06_translation_burden.py               # ~1 s    burden magnitude
python scripts/07_removed_results.py                  # ~4 s    the two dropped results
cd scripts && python 03_fig1_envelope.py && python 04_fig2_axis.py && python 05_fig3_bounds.py
python scripts/09_supraadditivity.py                  # ~28 s   the distinguishing prediction
cd scripts && python 10_fig4_supraadditivity.py       # Fig 4
python scripts/08_make_tables.py                      # ~1 s    tables
python scripts/11_headroom_sensitivity.py             # ~1 s    headroom across two axes
python scripts/12_chaperone_availability.py           # ~1 s    chaperone availability (theta)
python -m unittest discover -s tests -v               # 60 tests
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
  08_make_tables.py        Tables 1-5 and S1-S6
  09_supraadditivity.py    2x2 factorial: the distinguishing prediction
  10_fig4_supraadditivity.py  Fig 4  supraadditivity
  11_headroom_sensitivity.py  headroom vs evaluation point x chaperone anchoring
  12_chaperone_availability.py  headroom vs chaperone availability theta
  vendor/two_pool_ode.py   upstream model, verbatim (do not edit)
  run_all.py               rebuild everything, then run the tests
  figstyle.py              shared seaborn style
data/raw/                  staged inputs, self-contained
data/computed/             script outputs; the manuscript's only source of numbers
figures/                   Fig 1-4, PNG + PDF
tables/                    Tables 1-7 and S1-S6 as TSV, plus TABLES.md
tests/                     asserts the manuscript against data/computed/ and tables/
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
citation, softening a negative result, reinstating the retracted statistic,
substituting the corrupt tAI vector, drifting a generated table away from the
computed data, relabelling a null result as a signal, removing the "excluded"
marking from a dropped result, or reinstating the inflated ×158 headroom in the
abstract each make it fail. Headline claims are asserted against the *section*
that makes them, because a value that appears somewhere in a manuscript can still
be contradicted three sections later — which is how the earlier draft carried a
retracted statistic and its replacement at the same time.

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
- The framework's distinguishing prediction is untested **experimentally**.
  Result 5 confirms it holds in the model but is only ~0.2% above additive at
  wild-type margin, so the experiment must first compress the viability margin;
  the manuscript says so plainly rather than implying it is easy to test.
- Headroom is a range, not a point: ×1.9–×25 across the documented `C_tot`/`K_d`
  ranges crossed with chaperone availability θ. Any single figure quoted from it
  must say which evaluation point and which θ it came from.
- **θ is unmeasured.** The model omits nascent-chain folding, so it hands the whole
  chaperone pool to the damaged-protein pool. Table 7 exposes that as θ and gives
  the decision threshold: at θ ≥ 0.90 E. coli already sits in the regime where
  burden and capacity perturbations compound. One measurement — chaperone occupancy
  by nascent-chain folding in exponentially growing E. coli — collapses the range
  to a point. A model that represents nascent-chain folding explicitly, rather
  than absorbing it into θ, is the right next refinement.
- ν is validated for E. coli only. No validated B. subtilis or S. cerevisiae
  supply vector exists; see `../../codon-deployment/data/computed/README_TAI_IS_CORRUPT.md`.

## Figures and tables

| Artifact | Content |
|---|---|
| Fig 1 | envelope regimes and the saddle-node threshold |
| Fig 2 | mu clustering, nu null, between/within-amino-acid variance |
| Fig 3 | paired bounds, the paired ratio, and headroom |
| Fig 4 | the distinguishing prediction tested in the model |
| Table 1 | burden terms and their flux operationalization |
| Table 2 (+2b) | bounds on the tolerable per-codon error rate, and derived statistics |
| Table 3 | permutation tests: both axes, both nulls, both mu statistics |
| Table 4 | **removed result** — metal sites under both backgrounds |
| Table 5 | supraadditivity: interaction vs remaining margin |
| Table 6 | headroom across evaluation point and chaperone anchoring |
| Table 7 | headroom vs chaperone availability θ, over documented C_tot/K_d |
| Table S1 | per-codon (mu, nu) coordinates, all 59 analysed codons |
| Table S2 | per-amino-acid operational spread |
| Table S3 | validation of the nu axis, both candidate vectors |
| Table S4 | **removed result** — cross-species conservation of Delta |
| Table S5 | supraadditivity effect-size grid at the observed rate |
| Table S6 | the two capacity knobs are not equivalent |

Tables 4 and S4 report analyses this paper *excludes*; both are labelled as such
in `tables/TABLES.md`, and a test enforces that labelling.
