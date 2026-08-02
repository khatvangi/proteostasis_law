# Falsifiable Predictions

These are empirical hypotheses, not theorem-level claims.

> **Phase 3 note.** `FOLD_THEOREM.md` sharpens prediction 1 into a contrast
> against the standard capacity-exhaustion alternative: collapse should occur
> while clearance runs at a small fraction of V_max (`s_a` of order 0.05-0.4),
> not at saturation. The observable is a ratio, so it does not require absolute
> capacity calibration. That file also imposes a design constraint absent below —
> growth rate must be held fixed, because growth rate sets the nascent-chain
> occupancy `nu`, whose drift compounds with the allocation knob.

1. **Burden-capacity synthetic interactions.** A translation perturbation that modestly increases site-weighted damage will have a disproportionately large effect when chaperone or degradation capacity is reduced. Interaction contrasts should exceed an additive null on a prespecified burden scale.
2. **Condition-dependent codon effects.** The same synonymous change will differ across nutrient, stress, and growth conditions because `E_p`, `Pr(r|c,u,e)`, nascent occupancy, and rescue allocation change. A context-invariant codon ranking would count against the model's condition-dependent form.
3. **Substitution-identity dependence.** At matched mean error frequency, perturbations producing different substituted residues or targeting sites of different criticality will yield different folding and aggregation burdens. A frequency-only model should lose predictive accuracy.
4. **Rescue-allocation compensation.** Increasing the specifically limiting rescue allocation should shift the feasibility boundary and preferentially rescue high-burden strategies; indiscriminate overexpression need not help if it imposes occupancy or energetic costs.
5. **Proximity-to-boundary scaling.** Recovery time, variance, and perturbation sensitivity should rise as a stable state approaches a generic fold boundary, with the precise scaling tested against model-specific alternatives.
6. **Multi-proxy superiority.** A preregistered predictor combining expression flux, codon- and condition-specific substitution spectra, site damage/criticality, folding load, chaperone occupancy, and degradation capacity should outperform single variables such as mean mistranslation rate or expression alone on held-out conditions.

### Discriminating failures

The framework would require revision if independently measured loads and conserved capacities systematically predict no feasible state in thriving cells; if substitution/site information adds no reproducible predictive value across well-powered tests; or if inferred boundary proximity fails to predict any dynamic sensitivity after measurement error and confounding are addressed.

