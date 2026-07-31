# A finite proteostasis envelope for translation

Active manuscript: [`manuscript/MANUSCRIPT.md`](manuscript/MANUSCRIPT.md)
Assembled paper (tables inlined, figures embedded): `manuscript/PAPER.md`,
`PAPER.html`, `PAPER.pdf`, `PAPER.docx` — all generated, do not edit.

This is the canonical home of the proteostasis-law project as of 2026-07-29. It
supersedes `../combined-paper/` and `../proteostasis-paper/`.

## The claim

Extant translation operates inside a finite proteostasis envelope: viability
requires `B_total ≤ C_buffer`, where burden decomposes into decoding error,
folding failure, aggregation, and quality-control demand. A reduced dynamical
model with saturable rescue and superlinear aggregation feedback necessarily has
a saddle-node threshold. Two bounding arguments place E. coli roughly **one**
order of magnitude inside its envelope (×25 at the usage-weighted evaluation
point, ×1.9 to ×25 across chaperone availability θ). At that margin the
framework's distinguishing prediction is a synthetic lethality reachable in wild
type. Separately, measured per-codon mistranslation is organized at the
amino-acid level, while translational supply shows no comparable structure within
a stated power bound.

This is **not** a code-origin claim. It says nothing about why codons are
triplets or why degeneracy exists.

## Reproducing everything

One command rebuilds every number, figure, and table, assembles the paper, and
runs the tests. Nothing needs network access; all inputs are in `data/raw/`.

```bash
python scripts/run_all.py            # ~10 min end to end, then 89 tests
python scripts/run_all.py --fast     # skip the 10,000-permutation runs
python scripts/run_all.py --no-tests # rebuild only
```

The pipeline fails closed twice: `01_validate_tai.py` exits non-zero if the
supply axis is invalid, and `15_build_paper.py` exits non-zero on an unresolved
table placeholder or a missing figure.

**A rebuild is byte-reproducible, with two exceptions.** Every computed output,
figure (PNG/PDF/SVG), table, `PAPER.md` and `PAPER.html` reproduces exactly, so
`git status` after a verification rebuild is empty and any diff is a real change.
This took fixing two things a rebuild-and-diff pass exposed: matplotlib writes a
build timestamp and address-derived SVG element ids (now `svg.hashsalt` plus
`metadata={"Date": None}`), and seaborn's stripplot drew Fig 2c's jitter from the
unseeded global RNG, so the same data gave a visibly different panel every build
(now `figstyle.JITTER_SEED`). The exceptions are `PAPER.pdf` and `PAPER.docx`:
Chrome and pandoc each stamp a creation date into the container, so those two
files differ on every build even when their content is identical. Expect them in
`git status`; ignore them.

Individual steps, if you want them:

```bash
python scripts/01_validate_tai.py                     # ~3 s    validates nu, aborts if it fails
python scripts/02_axis_structure.py                   # ~65 s   permutation tests
python scripts/02_axis_structure.py --mu-stat median  # ~65 s   mu summary-statistic sensitivity
python scripts/06_translation_burden.py               # ~1 s    burden magnitude
python scripts/07_removed_results.py                  # ~4 s    the two excluded analyses
cd scripts && python 03_fig1_envelope.py && python 04_fig2_axis.py && python 05_fig3_bounds.py
python scripts/09_supraadditivity.py                  # ~28 s   the distinguishing prediction
cd scripts && python 10_fig4_supraadditivity.py       # Fig 4
python scripts/11_headroom_sensitivity.py             # ~1 s    headroom across two axes
python scripts/12_chaperone_availability.py           # ~1 s    chaperone availability (theta)
python scripts/13_nu_power.py                         # ~3 min  axis power / minimum detectable effect
python scripts/14_mu_jackknife.py                     # ~1 min  mu jackknife + detectability
python scripts/08_make_tables.py                      # ~1 s    tables
python scripts/15_build_paper.py                      # ~5 s    assemble PAPER.md/.html/.pdf/.docx
python -m unittest discover -s tests -v               # 89 tests
```

Requires Python 3 with numpy, pandas, scipy, seaborn, openpyxl. `PAPER.html` and
`PAPER.docx` need pandoc; `PAPER.pdf` is printed from the HTML by headless
Chrome. All three are skipped with a message if the tool is absent — the markdown
is always written. No unicode-capable TeX engine is installed, which is why the
PDF route is HTML rather than LaTeX.

## Layout

```
manuscript/
  MANUSCRIPT.md            the paper: prose, figures, <!-- TABLE:... --> placeholders
  manuscript_v2_draft.md   the v2 restructuring, kept verbatim as the drafting source
  PAPER.{md,html,pdf,docx} generated: prose + inlined tables + embedded figures
scripts/                   every number and figure is generated here
  01_validate_tai.py       validates the tAI (supply) axis; ABORTS if it fails
  02_axis_structure.py     Delta_A, permutation nulls, variance decomposition
  03_fig1_envelope.py      Fig 1  ODE regimes and the saddle-node
  04_fig2_axis.py          Fig 2  mu clustering, nu null, variance partition
  05_fig3_bounds.py        Fig 3  paired bounds and headroom
  06_translation_burden.py burden magnitude implied by the mu data
  07_removed_results.py    reproduces the two analyses this paper excludes
  08_make_tables.py        Tables 1-5, S1-S7, and the three not-paper tables
  09_supraadditivity.py    2x2 factorial: the distinguishing prediction
  10_fig4_supraadditivity.py  Fig 4  supraadditivity
  11_headroom_sensitivity.py  headroom vs evaluation point x chaperone anchoring
  12_chaperone_availability.py  headroom vs chaperone availability theta
  13_nu_power.py           minimum detectable effect and power on each axis
  14_mu_jackknife.py       leave-one-codon-out jackknife + detectability exposure
  15_build_paper.py        inline the tables, embed the figures, render the paper
  vendor/two_pool_ode.py   upstream model, verbatim (do not edit)
  run_all.py               rebuild everything, then run the tests
  figstyle.py              shared seaborn style
data/raw/                  staged inputs, self-contained
data/computed/             script outputs; the manuscript's only source of numbers
figures/                   Fig 1-4, PNG + PDF
tables/                    Tables 1-5, S1-S7 as TSV, plus TABLES.md
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
figure, hand-typing a table row into the prose, citing a table that does not
exist, dropping a citation, mixing citation styles, softening a negative result,
reinstating the retracted statistic or the window-bottom saturation triple,
substituting the corrupt tAI vector, drifting a generated table or its caption
away from the computed data, relabelling a null result as a signal, removing the
"excluded" marking from a dropped analysis, or reinstating the inflated ×158
headroom in the abstract each make it fail.

Two structural rules follow from the same failure:

- **Assert relationships, not just values.** Per-number provenance is satisfied by
  a document whose sections are each sourced and mutually contradictory. That is
  literally what happened: `09_supraadditivity.py` hardcoded the evaluation point
  the headroom section spends a page rejecting, every number traced to code, the
  suite passed, and the paper understated its own central effect by 38×. Tests now
  compare the two artifacts that must agree.
- **Headline claims are scoped to their section.** A value that appears somewhere
  in a manuscript can still be contradicted three sections later.

Tables are never typed into the manuscript. `MANUSCRIPT.md` carries a
`<!-- TABLE:Table 3 -->` placeholder and `15_build_paper.py` fills it from
`tables/TABLES.md`; a test fails if a markdown table row appears in the prose.
Table captions interpolate their numbers from `data/computed/` for the same
reason — two captions went stale for a full commit after the manuscript had
corrected them.

## Archived analyses: what `scripts/07_removed_results.py` is

**Read this if you came to the code and found analyses the manuscript does not
discuss.** Two empirical results from the earlier combined draft failed
verification and were removed from the paper. The scripts and tables stay here
because the negatives are informative and because deleting them would hide a
correction, but they are **not results of this paper** and nothing in the
manuscript rests on them. The paper does not mention them: a reviewer should not
have to read about analyses that were withdrawn before submission.

| Archived analysis | Recomputed verdict |
|---|---|
| **Metal-binding-site codon enrichment.** Claimed all four two-codon liganding amino acids prefer the higher-supply synonym (OR 1.28–1.40, all p < 0.03). | Gene-level expression confound. Metalloproteins are enriched for abundant enzymes and abundant genes have stronger codon bias, so a genome-wide background confounds site-level selection with expression. Against non-metal positions of the *same* genes: Asp OR 1.19 (p = 0.17), Cys 1.22 (0.14), Glu 1.26 (0.14), His — the strongest published effect — 1.07 (0.62). 0 of 4 survive. A separate study extending the test to nine other definitions of structural criticality found all nine null. The annotation is also unreliable: only ~60% of annotated metal-binding residues sit at a position whose codon encodes the annotated amino acid, so roughly 40% are mispositioned in sequence space. |
| **Cross-species conservation of operational geometry.** Claimed per-amino-acid Δ_A is conserved across E. coli, B. subtilis, and S. cerevisiae (ρ = 0.63–0.92). | Two artifacts. Δ_A was computed in combined (μ, ν) space using **E. coli μ for all three species**, so the correlated vectors shared an identical coordinate by construction — reproduced exactly, to 8.9e-16. And the supposedly species-specific supply vectors were not independent: 44 of 60 E. coli tAI values are bit-identical to the yeast values. Recomputing on the supply axis alone with independently derived weights gives ρ = 0.46–0.73, one of three comparisons non-significant, in a different rank order. Testing the claim needs genuinely species-specific tRNA gene counts. |

Both failures were denominator failures rather than measurement failures: in the
metal-site case the site data and the test statistic were identical in the
surviving and the collapsing version, and only the comparison set differed.

Reproduce with `python scripts/07_removed_results.py`. Outputs:
`tables/Excluded_metal_site_backgrounds.tsv`,
`tables/Excluded_crossspecies.tsv` — named `Excluded_*` rather than `Table*` so
they cannot be mistaken for tables of the paper, and asserted by the test suite
to stay that way.

## Where the caveats live

The paper has **no Limitations section, by design.** Every caveat is stated at the
claim it bounds: mass-spectrometry detectability inside the codon result, θ inside
the margin result, the ν power narrowing inside the supply paragraph, the
shared-parameter concession inside the bounds result, the single-species scope
where the coordinates are defined, and the absence of archaeal and mammalian
evidence where capacity is established. A standalone section on top of that would
state each concession twice and hand a reviewer two slightly different phrasings
of the same admission to compare. A test asserts each caveat against the section
that carries it, so deleting one inline does not silently pass.

## Known gaps before submission

- Public deposition DOI is not yet assigned (`Data and code availability`).
- The framework's distinguishing prediction is untested **experimentally**. It is
  confirmed *in the model*, and at the paper's own operating point it is sharp:
  +7.4% above additive for a 3× error increase crossed with a 3× rescue
  knockdown, and 12 of 36 tested pairs survivable alone but lethal together. No
  sensitized background is required. (An earlier draft evaluated this factorial at
  the window bottom, reported +0.2%, and concluded the margin had to be compressed
  first. That was an artifact of the evaluation point.)
- Headroom is a range, not a point: ×1.9–×25 across the documented `C_tot`/`K_d`
  ranges crossed with chaperone availability θ. Any single figure quoted from it
  must say which evaluation point and which θ it came from.
- **θ is unmeasured.** The model omits nascent-chain folding, so it hands the whole
  chaperone pool to the damaged-protein pool. Table 3 exposes that as θ and gives
  the decision threshold: at θ ≥ 0.90 E. coli already sits in the regime where
  burden and capacity perturbations compound. One measurement — chaperone occupancy
  by nascent-chain folding in exponentially growing E. coli — collapses the range
  to a point. A model that represents nascent-chain folding explicitly, rather
  than absorbing it into θ, is the right next refinement.
- **The μ axis result cannot be attributed to decoding.** Mass-spectrometry
  sampling depth carries as much amino-acid-level structure as μ itself
  (η² = 0.560 against 0.556). The paper says so inside the codon result itself.
  Separating them needs per-codon rates from a method whose detection probability
  does not vary with amino acid identity.
- ν is validated for E. coli only. No validated B. subtilis or S. cerevisiae
  supply vector exists; see `../codon-deployment/data/computed/README_TAI_IS_CORRUPT.md`.

## Figures and tables

| Artifact | Content |
|---|---|
| Fig 1 | envelope regimes and the saddle-node threshold |
| Fig 2 | mu clustering, nu null, between/within-amino-acid variance |
| Fig 3 | paired bounds, the paired ratio, and headroom |
| Fig 4 | the distinguishing prediction tested in the model |
| Table 1 | burden terms and their flux operationalization |
| Table 2 (+2b) | bounds on the tolerable per-codon error rate, and derived statistics |
| Table 3 | headroom vs chaperone availability θ, over documented C_tot/K_d |
| Table 4 | supraadditivity: interaction vs remaining margin |
| Table 5 | permutation tests: both axes, both nulls, both mu statistics |
| Table S1 | per-codon (mu, nu) coordinates, all 59 analysed codons |
| Table S2 | per-amino-acid operational spread |
| Table S3 | validation of the nu axis, both candidate vectors |
| Table S4 | supraadditivity effect-size grid at the operating point |
| Table S5 | the two capacity knobs are not equivalent |
| Table S6 | what each axis test could have detected (minimum detectable effect) |
| Table S7 | leave-one-codon-out jackknife on the μ axis |
| *not paper tables* | `Excluded_metal_site_backgrounds`, `Excluded_crossspecies`, `Retired_anchoring_grid` |

The numbering changed with the v2 restructuring: the two excluded analyses are no
longer Tables 4 and S4, and the retired anchoring grid is no longer Table 6/S9.
`08_make_tables.py` deletes the files from the old numbering on every run, so a
stale `Table7_*.tsv` cannot sit in `tables/` looking current.
