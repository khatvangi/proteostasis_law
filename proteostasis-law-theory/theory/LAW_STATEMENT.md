# Canonical Law Statement

Let a translation strategy `u` specify condition-dependent choices affecting initiation and elongation, codon-specific decoding, proofreading, amino-acid supply, and quality-control allocation. Let `P` be the attainable Pareto surface in objective space, whose coordinates include productive throughput, decoding fidelity, folding yield, functional/chemical diversity, energetic cost, and quality-control demand.

For environment `e`, strategy `u` produces a site-resolved damage influx `J_U(u,e)` and other ordinary folding loads. These drive a finite-resource proteostasis dynamical system

`dx/dt = F(x; u, e, theta)`,

with burden state `x`, parameters `theta`, conserved chaperone and degradation pools, and admissible damage region `D = {x: h_j(x) < H_j for all j}`.

## Proteostasis Law

**A translation strategy is viable in environment `e` only if the coupled proteostasis dynamics possess an attracting bounded invariant state (an equilibrium or bounded attractor) reached from the biologically relevant initial set and contained within `D`.**

The feasible subset of the Pareto surface is therefore

`P_feas(e) = {u in P: exists an attracting bounded invariant state B_u subset D with the relevant basin containing X_0}`.

This is a necessary feasibility condition, not a complete theory of fitness. Within `P_feas`, selection, drift, regulation, ecological interactions, and history may distinguish strategies.

## Claim status

- **Definition:** `P_feas` and the viability criterion above.
- **Theorem-level statement:** if no bounded attracting state below threshold exists for the specified model, parameters, and initial set, that strategy is not viable under those assumptions. This follows directly from the definition.
- **Model consequence:** folds or other stability-boundary crossings can remove feasible low-burden states in particular nonlinear models.
- **Empirical hypothesis:** observed viable strategies lie inside independently estimated feasibility boundaries.
- **Speculative implication:** evolution may concentrate strategies near, but generally not exactly on, boundaries when performance gains oppose robustness.

