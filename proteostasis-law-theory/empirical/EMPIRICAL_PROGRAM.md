# Empirical Program

## Objective

Test whether site-resolved translation damage and conserved rescue capacity predict bounded sub-threshold proteostasis better than scalar error or expression measures.

## Stage 1 — Operationalization

Choose one tractable organism and two to three conditions. Define `U`, `A`, viability thresholds, observation time, and the translation strategy perturbations before fitting. Calibrate units and measure ordinary nascent-chain chaperone occupancy.

## Stage 2 — Factorial perturbations

Construct perturbations spanning: (a) decoding spectrum or codon context, (b) expression flux, (c) site/substitution identity, (d) chaperone allocation, and (e) degradation capacity. Include matched-error-frequency pairs with different substitution spectra, and matched site substitutions at different expression levels. Randomize batches and include burden-neutral expression controls.

## Stage 3 — Dynamic measurement

Measure synthesis flux, decoding outcomes, soluble burden, aggregate burden, free and occupied resource pools, recovery time, viability, and function through time. Use pulse-chase or acute induction/withdrawal to distinguish influx, repair, and clearance.

## Stage 4 — Model fitting and tests

Fit the conserved-resource `U,A` model with hierarchical measurement error. Check identifiability by simulation before interpreting parameters. Estimate equilibria and Jacobian eigenvalues with uncertainty. Test the six predictions in `theory/PREDICTIONS.md` using prespecified contrasts.

Compare held-out predictive performance for:

1. mean mistranslation rate alone;
2. expression plus mean rate;
3. site- and identity-resolved influx;
4. influx plus folding load and conserved capacity.

## Stage 5 — Boundary challenge

Select strategies inferred at increasing distance from the feasibility boundary. Apply small burden or capacity perturbations and test predicted recovery time and sensitivity. A rescue experiment should alter the identified limiting allocation, not merely total protein expression.

## Decision criteria

Support requires reproducible out-of-sample improvement, correct qualitative interaction signs, and dynamic trajectories consistent with conserved resource occupancy. Failure modes, alternative models, and negative results will be retained. Organism-level generalization begins only after within-system validation.

**Literature placeholders:** [Citations needed: quantitative mistranslation spectra]; [Citations needed: cotranslational chaperone occupancy]; [Citations needed: aggregate clearance kinetics]; [Citations needed: early-warning behavior near cellular stability transitions].

