# Proteostasis Law Theory

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21794565.svg)](https://doi.org/10.5281/zenodo.21794565)

Archived at [10.5281/zenodo.21794565](https://doi.org/10.5281/zenodo.21794565) (concept DOI, resolves to latest). Canonical manuscript: [`manuscript/bmb_v4.md`](manuscript/bmb_v4.md).

This repository rebuilds the Proteostasis Law from first principles. It is a clean, theory-first project: translation is treated as a multidimensional allocation problem, and viability is determined by the stability of the proteostasis system that receives translation's products.

## Canonical claim

Translation strategies occupy a Pareto trade-off surface involving productive throughput, decoding fidelity, folding yield, functional and chemical diversity, energetic cost, and quality-control demand. A strategy is biologically feasible only if its site-resolved damage influx, coupled to finite chaperone and degradation resources, admits a bounded stable state below specified damage thresholds. The Proteostasis Law is this feasibility restriction on the Pareto surface—not a claim that one scalar error rate determines fitness.

## Directory map

- `manuscript/`: **`bmb_v4.md` is canonical** and is the only manuscript in use. `MANUSCRIPT.md` and `COLLAPSE_BOUNDARY.md` are superseded earlier drafts, retained for provenance and banner-marked as such; see `manuscript/README.md`.
- `theory/`: law statement, mathematics, Pareto geometry, dynamics, predictions, scope, and nonclaims. `FOLD_THEOREM.md` derives where collapse occurs — the saddle-node is a constrained critical point of total removal on the aggregate nullcline (*critical point*, not maximum; see D011) — decomposes the margin into its saturation and sequestration deficits, extends the result to dividing cells, calibrates the growth-burden coupling to a measurement, and computes the Pareto surface the framework had only asserted.
- `empirical/`: experimental program, measurement-to-symbol map, and `GATE4_PROPOSAL.md` — what a first empirical test of the fold theorem would require, specified without reading any outcome.
- `notes/`: historical recovery and rejection ledger; provenance stays here.
- `scripts/`: the analysis code. `phase3/` derives and checks the theorem; `figures/` builds every figure in the manuscript.
- `tests/`: the check suite, asserting the manuscript's numbers against the code that produces them rather than against stored values.
- `data/figures/`: reduced arrays backing each figure, with per-file provenance (population, count, whether it is a subsample) and a sha256 manifest, so the figures rebuild on a clean checkout without the run root.

## Reproducing the figures

```
python scripts/figures/build_figure_data.py     # only script permitted to read the run root
python scripts/figures/fig_theorem.py           # Fig. 1  the theorem, two panels
python scripts/figures/fig_dilution.py          # Fig. 2  phi_enz flat, delta carries division
python scripts/figures/fig_saturation.py        # Fig. 3  saturation at the boundary
python scripts/figures/fig_front.py             # Fig. 4  the strategy front
python scripts/figures/fig_beta.py              # Fig. 5  the beta-indexed prediction
python scripts/figures/fig_identity.py          # Fig. S1 parallelism is bracket tolerance
```

Script names are deliberately **semantic rather than positional**. Manuscript figures
are numbered by order of first mention, so inserting one renumbers those after it;
each script declares its number in a single `FIGURE` constant and writes
`figures/figN.*` from that, so a reorder touches one line per script instead of
renaming files. `tests/figures/test_figure_wiring.py` asserts that the numbering
still follows first mention, that every figure is cited in the prose before it is
embedded, and that no caption number disagrees with the text number for the same
quantity.

Each writes PDF, SVG and PNG to `figures/` with timestamps stripped, so an unchanged
input gives byte-identical output and a rebuild produces a clean diff. No seaborn and
no dependence on a local matplotlib configuration: every rcParam that affects the
output is set explicitly, because a figure that renders differently on the
typesetter's machine than on ours is a defect.

## Building the submission PDF

```
python scripts/manuscript/to_latex.py
```

Writes `manuscript/bmb_v4.{tex,pdf}` and `manuscript/bmb_v4_supplementary.{tex,pdf}`.
The markdown is the single source; **the `.tex` files are build artefacts and must
never be edited.** The converter treats every backtick span as mathematics unless
it is declared in `CODE_SPANS`, so a new code-shaped span fails the build rather
than silently rendering as math, and it maps unicode from the document's own
inventory, so an unmapped character stops the build rather than vanishing.

The generated preamble uses `article` because no Springer class is installed here;
a banner at the top of the `.tex` gives the one-line change to `sn-jnl`.

Run the checks with `python -m pytest tests/ -q`.

## Project rules

1. Keep theorem-level results, model consequences, empirical hypotheses, and evolutionary speculation explicitly separated.
2. Preserve site and substitution identity in the damage influx; scalar mistranslation rates are coarse summaries only.
3. Conserve finite rescue resources and count ordinary nascent-chain folding occupancy.
4. Check dimensional consistency before combining quantities.
5. Treat the one-variable model only as an illustrative reduction; it does not establish bistability or hysteresis.
6. Add numerical claims only with traceable data, uncertainty, units, and executable analysis.
7. Do not reintroduce any item recorded in `notes/REJECTED_COMPONENTS.md` without new evidence and a recorded decision.

## Immediate next steps

1. Analyze existence, local stability, and fold boundaries in the conserved-resource dynamical system.
2. Choose measurable damage thresholds and preregister perturbation contrasts.
3. Assemble site-resolved decoding, expression, folding, and capacity measurements with uncertainty.
4. Fit and compare multi-proxy and single-variable predictors on held-out conditions.
5. Add literature citations at the marked placeholders and only then calibrate organism-specific models.

