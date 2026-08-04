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
- `scripts/`: future reproducible analysis code and conventions.
- `tests/`: future mathematical, numerical, and data-validation checks.

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

