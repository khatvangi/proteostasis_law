# SUPERSEDED by ../envelope-paper/ on 2026-07-29

`P1_combined_capacity_codon_DRAFT.md` should not be submitted or cited. The
active manuscript is `../envelope-paper/manuscript/MANUSCRIPT.md`.

## Why

Verification against the upstream sources found four blocking problems.

1. **The ν (tAI) axis was corrupt.** `codon-deployment/data/computed/ecoli_tai_ws.tsv`
   had 44 of 60 values bit-identical to the yeast file and assigned CTG — E. coli's
   dominant Leu codon — a near-floor supply weight. Its generating script neither
   reproduced those files nor implemented dos Reis 2004 correctly. See
   `../../codon-deployment/data/computed/README_TAI_IS_CORRUPT.md`.

2. **Result 3 (metal sites) was a gene-level expression confound.** Against an
   expression-matched within-gene background, all four effects vanish (p = 0.14–0.62)
   and His falls from OR 1.33 to 1.07. The repo's own `metal_enrichment_legacy.tsv`
   already showed this with a narrower background. A follow-up project
   (`structural-criticality/`) independently failed to replicate it across nine
   other criticality definitions. Additionally ~40% of the annotated metal-binding
   residues sit at a sequence position whose codon does not encode the annotated
   amino acid.

3. **Result 4 (cross-species conservation) was doubly artifactual.** Δ_A was the
   2D (μ, ν) distance computed with E. coli μ for all three species, so the
   correlated vectors shared an identical coordinate by construction — reproduced
   exactly by `../envelope-paper/scripts/07_removed_results.py`.

4. **Result 6 quoted superseded and retracted numbers.** The ODE MC median
   1.7 × 10⁻² came from the `backup_uniformN` run (current: 2.59 × 10⁻² paired),
   and "paired P(arith < ODE) = 0.65" is the statistic
   `proteostasis-P1/PAIRED_MC_TASK.md` explicitly calls "not a valid paired
   statistic". The valid value is 0.768.

Also: 20 of the draft's 36 references were never cited, and its
error-minimization citation `[6,7]` pointed at a chaperone protocol and a
synonymous-edit paper.

## What was kept

The theory, the `dP/dτ` model and its saddle-node threshold, the two-bound
convergence, the headroom result, and the μ-axis clustering all reproduce and
carry over into `../envelope-paper/`. `B1_ODE_arithmetic_bound_resolution.md`
reached the right qualitative conclusion (the ODE bound sits *above* arithmetic,
and the closure mechanism is aggregation-death) but quoted the retracted paired
statistic and the backup-run median.

## Status of the codon-deployment material

Not abandoned, but not publishable as drafted. A corrected version would be a
negative-result paper: μ clustering holds, metal-site deployment does not, and
cross-species geometry cannot be tested without real GtRNAdb tRNA counts for
B. subtilis and S. cerevisiae. `proteostasis_family_inventory.md` and the
outline are retained for that purpose.
