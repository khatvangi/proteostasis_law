> **SUPERSEDED.** The canonical manuscript is [`bmb_v4.md`](bmb_v4.md). This file is kept for provenance and is not the version any published DOI refers to.

# The Proteostasis Law: A Stability Boundary on the Pareto Surface of Translation

## Abstract

Translation must produce proteins quickly and accurately while supporting diverse chemical functions, enabling productive folding, limiting energetic expenditure, and avoiding excessive quality-control demand. These objectives cannot generally be optimized independently. We propose that attainable translation strategies occupy a condition-dependent Pareto surface and that proteostasis imposes a feasibility boundary on this surface. The molecular input to that boundary is not a single mistranslation frequency. We define a proteome-integrated, site-resolved damage influx that weights protein synthesis flux, codon- and condition-specific substitution probabilities, site criticality, and the damage caused by a particular residue substitution. This influx drives a minimal dynamical system for soluble misfolded protein and aggregates, coupled to explicitly conserved chaperone and degradation pools. Ordinary nascent-chain folding occupies the same finite resources and therefore contributes to capacity demand even when it is not damage. The Proteostasis Law states that a translation strategy is viable only if the coupled dynamics admit a bounded attracting state below specified damage thresholds for biologically relevant initial conditions. This formulation separates a definition-level feasibility statement from conditional model behavior, empirical hypotheses, and evolutionary speculation. It predicts burden-by-capacity synthetic interactions, condition-specific codon effects, substitution-identity dependence, compensatory rescue allocation, increased sensitivity near stability boundaries, and superior prediction from multiple mechanistic proxies compared with any single scalar. The framework supplies a rigorous starting point for quantitative tests without claiming a universal geometry, organism-independent parameters, or a complete account of fitness.

## Introduction

Translation connects genetic information to the physical proteome. Its performance is often summarized by speed or error frequency, but a cell must solve a broader allocation problem. Ribosomes must deliver productive throughput; decoding and proofreading must control substitution spectra; nascent chains must reach functional conformations; the proteome must realize chemically diverse tasks; and all of this consumes energy and quality-control capacity. Improvement in one coordinate can worsen another. A slower, highly proofread strategy may reduce some errors while sacrificing throughput and energy. Faster synthesis can increase ordinary cotranslational chaperone occupancy even without changing the probability of mistranslation. The chemical identity and position of an error can matter more than its contribution to a proteome-wide mean.

The receiving system is also finite. Chaperones, proteases, disaggregases, and associated energy supplies are occupied by normal and abnormal clients. Damaged monomers can be refolded or degraded, but they can also nucleate aggregates; aggregates can grow, be disaggregated, or be cleared. A viable translation program must therefore be considered jointly with the state and allocation of proteostasis resources.

We develop a theory with two linked objects. First, translation strategies map to a multidimensional Pareto surface. Second, a dynamical feasibility filter maps their site-resolved damage and ordinary folding loads into bounded or unbounded proteostasis outcomes. The result is a feasible subset of the Pareto surface. This construction is deliberately a necessary condition for viability rather than a full fitness function. It is intended to organize measurements and perturbations around a precise question: does the translation–proteostasis system possess a stable, sub-threshold operating state?

Four logical categories are kept distinct throughout. Definitions and theorem-level statements follow from explicit mathematics. Model consequences depend on selected kinetic forms and parameter regimes. Empirical hypotheses can be falsified by measurement. Evolutionary implications are speculative unless separately supported.

## Translation as a Pareto problem

Let `u` denote a translation strategy available in environment `e`. It includes condition-dependent choices or states affecting initiation, elongation, codon-specific decoding, proofreading, amino-acid supply, and quality-control allocation. Define the objective map

`Phi(u,e) = (T, F, Y, V, -G, -Q)`,

where `T` is productive throughput, `F` decoding fidelity, `Y` folding yield, `V` functional and chemical diversity, `G` energetic cost, and `Q` quality-control demand. The negative signs orient all coordinates so that larger is preferable. Each coordinate requires an operational measurement and may itself be vector valued in a detailed implementation.

A strategy is Pareto efficient when no attainable alternative improves every coordinate weakly and at least one strictly. The resulting surface does not assume that all pairwise trade-offs occur, that the surface is smooth, or that a single axis captures its geometry. It is conditional on environment because substrate supply, expression, decoding, stress regulation, and resource allocation change.

Pareto efficiency alone does not imply viability. A strategy can be nondominated yet deliver loads that overwhelm its available proteostasis machinery. Conversely, several viable strategies may remain incomparable. Proteostasis therefore restricts rather than uniquely optimizes the translation surface.

## Formal law

Let `P_e` be the Pareto set in environment `e`. Translation strategy `u` generates site-resolved abnormal-protein influx and ordinary nascent-chain load. Let the coupled proteostasis state be `x`, with dynamics

`dx/dt = F(x;u,e,theta)`,

where `theta` contains kinetic and resource parameters. Let `D_H = {x: h_j(x)<H_j for every monitored damage coordinate j}` be the admissible region, and let `X_0` denote biologically relevant initial conditions.

**Proteostasis Law.** A translation strategy is viable in environment `e` only if the coupled proteostasis dynamics possess a bounded attracting invariant state `B_u` contained in `D_H`, and the relevant initial set lies in its basin of attraction.

The feasible subset is

`P_feas(e) = {u in P_e: there exists B_u subset D_H, B_u is bounded and attracting, and X_0 subset basin(B_u)}`.

This is a definition and necessary feasibility condition. The immediate theorem-level statement is correspondingly bounded: for a specified model, parameters, thresholds, and initial set, absence of such a state excludes viability under those assumptions. The law does not say that every member of `P_feas` has equal fitness or that every excluded strategy fails on the same time scale.

## Site-resolved damage influx

Consider protein `p`, site `i`, cognate residue `a(p,i)`, and encoded codon `c(p,i)`. Let `E_p(u,e)` be the translation flux of protein `p`; `Pr(r|c(p,i),u,e)` the probability that residue `r` is delivered; `kappa_{p,i}` a dimensionless site-criticality modifier; and `d(r,p,i,e)` the expected contribution of that substitution to entry into the modeled damaged pool. Define

`J_U(u,e) = sum_p E_p(u,e) sum_i kappa_{p,i} sum_{r != a(p,i)} Pr(r|c(p,i),u,e) d(r,p,i,e)`.

With `E_p` measured in molecules per volume per time and `d` expressed as damaged-pool units per translated molecule, `J_U` has influx units. The definitions of `kappa` and `d` must avoid double counting: for example, `kappa` may encode a preregistered dimensionless importance class while `d` provides the calibrated physical burden. Alternatively, they can be combined into one substitution-and-site damage kernel.

Non-error folding failure can enter as a distinct `J_fold`, and other damage sources as `J_other`, giving `J=J_U+J_fold+J_other`. Ordinary nascent chains that fold successfully are not added to `J`, but their chaperone occupancy affects available capacity.

A proteome-wide mean mistranslation rate is sometimes useful as a coarse approximation:

`epsilon_bar = [sum_p E_p sum_i sum_{r != a} Pr(r|c,u,e)]/[sum_p E_p L_p]`.

It is insufficient in general. Two strategies with the same `epsilon_bar` can differ in the proteins expressed, the sites affected, the substituted amino acids, and their folding or functional consequences. Reduction to `epsilon_bar` is justified only if these omitted weights are effectively constant or their variation is irrelevant at the resolution of the question—an assumption to test, not a property of translation.

## Finite-resource proteostasis dynamics

The smallest model that explicitly represents aggregation contains soluble damaged or misfolded monomer `U` and aggregate burden `A`:

`dU/dt = J - v_ref(U,C_f) - v_degU(U,D_f) - n(U) - g(U,A) + v_dis(A,C_f)`,

`dA/dt = n(U) + g(U,A) - v_dis(A,C_f) - v_degA(A,D_f)`.

Here refolding and soluble degradation remove `U`; nucleation and aggregate growth transfer it into `A`; disaggregation returns aggregate material to `U` in the displayed stoichiometry; and aggregate degradation clears `A`. One admissible family is

`v_ref=k_ref C_f U/(K_ref+U)`, `v_degU=k_U D_f U/(K_U+U)`,

`n=k_n U^m` with `m>1`, `g=k_gUA`,

`v_dis=k_dis C_f A/(K_dis+A)`, and `v_degA=k_A D_f A/(K_A+A)`.

These kinetics are a model choice, not the law. Alternative rate forms are acceptable if mass flow, units, and resource accounting remain explicit.

Let total chaperone satisfy

`C_tot=C_f+C_N+C_U+C_A+C_O`,

where `C_N` is occupancy by ordinary nascent folding, `C_U` and `C_A` are damaged-substrate complexes, and `C_O` represents other clients. Under rapid equilibrium, a closed expression such as

`C_f=C_tot/(1+N_f/K_N+U_f/K_CU+A_f/K_CA+O_f/K_CO)`

is valid only when the substrates in the denominator are free binding-competent concentrations. If `U` and `A` are total pools, free substrate and complexes must be solved together from their conservation equations. Total substrate cannot be inserted as if it were free substrate. The degradation pool obeys the analogous balance

`D_tot=D_f+D_U+D_A+D_O`,

or a more detailed catalytic queue model. Thus increased ordinary translation can reduce free folding capacity through `C_N` even when error influx is unchanged.

At an equilibrium `(U*,A*)`, local asymptotic stability requires the eigenvalues of the two-dimensional Jacobian to have negative real parts. Viability further requires that burden thresholds are not crossed and relevant initial conditions converge to the state. Numerical continuation can map the resulting boundary in load-allocation space. [Placeholder: formal propositions will be added after kinetic assumptions and invariant domains are fixed.]

## An illustrative one-variable reduction

For graphical intuition, aggregation and resource allocation can be collapsed into

`dM/dt = j - aM/(K+M) + bM^2`, with positive parameters.

The vector field is strictly convex. It can have two positive roots, a double root, or no positive root. When two roots exist, the low-burden root is stable and the upper root is unstable. The upper root is a separatrix-like threshold in this one-dimensional setting; trajectories above it grow rather than approach a second stable state. At the double root, stable and unstable equilibria meet and disappear in a fold. Therefore this reduction illustrates a stable operating state, an unstable threshold, and fold loss. It neither contains a stable high-burden attractor nor demonstrates hysteresis. More complex dynamics may produce those behaviors, but only with additional state structure and analysis.

## Consequences of the framework

The first consequence is non-equivalence of error spectra. Mean error frequency can remain unchanged while `J_U` changes because translation flux, site criticality, or substitution damage changes. The second is load coupling: productive nascent chains and damaged proteins can compete for a conserved chaperone pool. The third is conditional nonlinearity. In models with folds or sharp capacity boundaries, small changes in damage influx or resource allocation can have large effects near the boundary while producing weak effects far from it. These are model consequences when the relevant kinetics and parameter regimes are established; they are not universal theorems about every saturable rescue system.

The geometric consequence is that proteostasis carves a condition-specific feasible subset from the Pareto surface. Compensation can move a strategy back into this subset by lowering site-weighted influx, increasing effective refolding or degradation, changing occupancy, or improving clearance. Because these interventions have their own energetic and allocation costs, the rescued strategy can occupy a different Pareto position rather than simply reversing the original perturbation.

## Falsifiable predictions

The theory motivates six linked hypotheses.

First, burden and capacity perturbations should interact synthetically. A modest decoding or folding burden should be more harmful when the specifically relevant chaperone or degradation pool is constrained than predicted by a prespecified additive null.

Second, codon effects should be condition dependent. The same synonymous change can alter outcomes across nutrient or stress states because expression flux, decoding probabilities, nascent occupancy, and rescue allocation differ.

Third, substitution identity and site should matter at matched error frequency. Perturbations that generate different residues or affect sites with different damage kernels should yield distinguishable soluble and aggregate burdens.

Fourth, compensation should be allocation specific. Increasing the limiting rescue activity should shift the inferred feasibility boundary, whereas increasing a nonlimiting component may fail or impose a new cost.

Fifth, systems approaching a generic fold boundary should show increased recovery time and perturbation sensitivity; variance may also rise under controlled noise. Quantitative scaling is conditional on the local normal form and must be compared with alternatives.

Sixth, multi-proxy prediction should outperform single-variable prediction out of sample. A model containing expression, substitution spectra, site damage, folding load, occupied and free chaperone, and degradation capacity should predict burden dynamics and threshold crossings better than mean mistranslation rate, growth rate, or expression alone.

These hypotheses risk falsification. The framework would require revision if well-calibrated conserved capacities and loads consistently exclude healthy states; if site and substitution information provides no reproducible predictive gain; or if estimated proximity to a boundary has no dynamic consequence across adequately powered perturbations.

## Empirical program

An initial test should focus on one experimentally tractable system rather than seek immediate universality. Two or more environments would establish whether condition dependence is measurable. The design should independently vary expression or decoding burden, substitution/site identity, chaperone allocation, and degradation capacity. Matched-frequency error spectra are especially informative because they directly compare the site-resolved influx with a scalar error model.

Required measurements include protein-specific synthesis flux; codon- and condition-specific substitution probabilities; substitution- and site-sensitive folding or damage effects; non-error folding yield; soluble misfolded and aggregate burden through time; total, free, and occupied chaperone and degradation pools; and functional or viability thresholds. Orthogonal burden proxies are needed because no single assay cleanly identifies all damaged states. [Citations needed: quantitative decoding-error methods, cotranslational folding occupancy, and time-resolved aggregate clearance.]

The two-state model should be fit with measurement error and hierarchical variation, after synthetic-data identifiability analysis. Candidate models should be compared on held-out perturbations: mean error alone; expression plus mean error; site-resolved influx; and site-resolved influx plus conserved resource occupancy. Dynamic perturbation and washout data can separate influx from refolding, aggregation, and clearance more effectively than steady state alone.

Finally, strategies inferred to lie at different distances from the boundary should receive small burden and capacity challenges. Recovery time, burden amplification, and rescue specificity provide sharper tests than correlation with baseline growth. Threshold definitions, exclusion rules, contrasts, and model comparisons should be preregistered. [Placeholder: select organism, perturbation library, assay panel, and power calculation after feasibility experiments.]

## Discussion

The Proteostasis Law reframes translation as a constrained source coupled to a finite-capacity dynamical receiver. Its central contribution is not a claim that translation errors are universally catastrophic, nor a numerical ceiling. It is a formal viability criterion: Pareto-efficient translation strategies remain biologically available only where their complete loads and allocations support a bounded, stable, sub-threshold proteostasis state.

Site resolution is essential to this formulation. Expression determines how often a site is translated; codon and condition shape the substitution spectrum; residue identity and local structure shape the physical consequence. The scalar error rate discards precisely the dimensions through which translation couples to folding and aggregation. It remains useful as a baseline and may prove adequate in restricted regimes, but adequacy must be demonstrated empirically.

Explicit conservation is equally important. Chaperones and degradation machinery are pools with occupancy, not abstract rate constants that remain available at arbitrary load. Productive nascent folding competes for these resources. This creates a direct mechanistic link between throughput and quality-control capacity without declaring that faster translation must always be damaging.

The law is intentionally bounded. It does not identify the evolutionary optimum within the feasible set, establish a universal Pareto shape, or imply that every system operates close to instability. A possible evolutionary implication is that gains in throughput or economy may sometimes bring populations near a proteostasis boundary, with regulatory reserves providing robustness. That proposal is speculative and requires comparative evidence after within-system tests succeed.

Several extensions are natural: compartments, stress-response feedback, energy limitation, cell division, phase-separated states, and population heterogeneity. Each can alter invariant sets and thresholds. The immediate priority, however, is not maximal complexity. It is to determine whether a measured site-resolved influx and conserved two-state receiver correctly predict synthetic interactions, compensation, and boundary responses in a controlled system.

## Conclusion

Translation strategies are viable not merely when they are fast or accurate, but when their multidimensional performance and induced loads are dynamically supportable. Defining that support as a stable, bounded, sub-threshold state converts a broad intuition into a testable law. The decisive experiments will measure damage composition and finite capacity together, perturb both sides of the coupling, and ask whether the resulting feasibility boundary predicts outcomes beyond any single proxy.

