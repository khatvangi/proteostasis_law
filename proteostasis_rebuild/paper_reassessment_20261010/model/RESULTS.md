# Results: finite-capacity client routes

**Empirical predictive claim: FAIL / not established.** This is a validated computational discrimination of model assumptions, not a fitted E. coli failure probability. GroEL obligation and TF/DnaK route coverage use the observations being reproduced. Neither is held-out validation.

200 independent seeded parameter draws (seed 20261010); the same draws, demand and total scalar capacities are used across all three nested models. Primary load is 0.9 times the full architecture’s all-client productive radial capacity. This guarantees a feasible baseline by design, not by observation.

## Matched single-machine losses

Each target loses δ = 0.8 min(C_DnaK,C_GroEL,C_ClpB). Thus remaining scalar capacity is exactly equal and no capacity is created. Values below are fractions of total client demand; ranges are empirical 5th/50th/95th percentiles of the specified generator, not biological confidence intervals.

| Model | Baseline folded fraction, median | GroEL loss: folding drop q05/q50/q95 | DnaK loss: folding drop q05/q50/q95 | ClpB loss: folding drop q05/q50/q95 | Nonzero target contrast |
|---|---:|---|---|---|---:|
| scalar | 1 | 0/0/0 | 0/0/0 | 0/0/0 | 0/200 |
| independent | 1 | 0/0/0 | 0/0/0 | 0/0.0081159/0.0658733 | 183/200 |
| architecture | 1 | 0/0/0 | 0/0/0 | 0/0.0081159/0.0658733 | 184/200 |

Matched absolute work loss q05/q50/q95: 1.91633e-05/0.000202693/0.00244289 (baseline total capacity 1.1). The small ClpB pool often limits common perturbation size; a weak contrast here does not identify interchangeable machines.

## Pair interactions

For folded fraction F, excess joint loss is I = F(a)+F(b)−F(ab)−F(0), positive for supra-additive loss. Each endpoint loses 0.4 min(Ca,Cb). The ab, aa and bb arms all lose the same total capacity. aa/bb concentrate twice the endpoint dose on one endpoint. These are finite-dose interactions, not genetic-interaction scores.

| Pair | Model | I q05/q50/q95 | I > 1e-7 | ab/aa/bb contrast > 1e-7 |
|---|---|---|---:|---:|
| DnaK_GroEL | scalar | 0/0/0 | 0/200 | 0/200 |
| DnaK_GroEL | independent | 0/0/0 | 0/200 | 36/200 |
| DnaK_GroEL | architecture | 0/0/0 | 0/200 | 45/200 |
| DnaK_ClpB | scalar | 0/0/0 | 0/200 | 0/200 |
| DnaK_ClpB | independent | 0/0/0 | 0/200 | 182/200 |
| DnaK_ClpB | architecture | 0/0/0 | 0/200 | 182/200 |
| GroEL_ClpB | scalar | 0/0/0 | 0/200 | 0/200 |
| GroEL_ClpB | independent | 0/0/0 | 0/200 | 183/200 |
| GroEL_ClpB | architecture | 0/0/0 | 0/200 | 184/200 |
| TF_DnaK | scalar | 0/0/0 | 0/200 | 0/200 |
| TF_DnaK | independent | 0/0/0 | 0/200 | 23/200 |
| TF_DnaK | architecture | 0/0/0 | 0/200 | 38/200 |

Positive interactions can arise in a scalar model from an ordinary capacity threshold, even when none occur in this primary grid. They cannot establish a biological serial mechanism. Matched arm contrasts can reject the scalar reduction within a specified network, but do not validate that network.

## Obligation and upstream-overlap ablations

Complete knockouts below are **not** matched-total perturbations. Demand remains fixed to the original architecture baseline, including in ablations. All bypasses have finite machine costs.

| Route assumption | Perturbation | Essential fraction q05/q50/q95 | Folding fraction q05/q50/q95 | Essential quota shortfall |
|---|---|---|---|---:|
| full | KO_GroEL | 0/0/0 | 0.696503/0.965669/0.994059 | 200/200 |
| full | KO_TF | 0.722101/1/1 | 0.0808446/0.763094/1 | 17/200 |
| full | KO_DnaK | 1/1/1 | 0.912504/0.987923/0.998655 | 0/200 |
| full | KO_TF_DnaK | 0/0/0 | 0.0502947/0.433369/0.846554 | 200/200 |
| remove_obligate | KO_GroEL | 0.722101/1/1 | 0.809197/1/1 | 17/200 |
| remove_upstream | KO_TF | 0.722101/1/1 | 0.114008/0.925891/1 | 17/200 |
| remove_upstream | KO_DnaK | 1/1/1 | 0.912504/0.987923/0.998655 | 0/200 |
| remove_upstream | KO_TF_DnaK | 0/0/0 | 0.0591039/0.544394/0.969424 | 200/200 |
| no_serial | KO_TF_DnaK | 1/1/1 | 0.0524567/0.505111/0.897238 | 0/200 |

**Tautology audit.** Zero essential productive flux after GroEL knockout is a theorem of excluding every GroEL-free productive route for the essential-obligate class. Likewise, loss of overlap-class rescue after deleting both TF and DnaK follows from specifying those as its only upstream alternatives. Their reproduction is encoded. Single-to-double magnitudes, competition, and which finite-dose perturbation is worst are conditional computations; these are not yet independent empirical predictions. Calloni also reports impaired GroEL function, ribosome biogenesis and aggregation; these effects are not represented by this static LP. Removing upstream dependence only for overlap clients does not remove the separate upstream requirement imposed on obligate clients; no_serial isolates that latter assumption.

## Load and catalytic-range sensitivity

All rows use the same 200 trial identities when run with the default configuration; load is recalibrated to the width-specific unperturbed radial capacity. Width w means independent log10 multipliers Uniform[−w,w], not inferred kinetic distributions.

| Width | Load | Baseline essential shortfalls | Matched target contrast | GroEL worst (strict) | DnaK worst (strict) | ClpB worst (strict) |
|---|---:|---:|---:|---:|---:|
| 1 | 0.6 | 0/200 | 175/200 | 0/200 | 0/200 | 175/200 |
| 1 | 0.9 | 0/200 | 184/200 | 4/200 | 0/200 | 180/200 |
| 1 | 1.2 | 0/200 | 191/200 | 26/200 | 0/200 | 165/200 |
| 2 | 0.6 | 0/200 | 134/200 | 15/200 | 1/200 | 118/200 |
| 2 | 0.9 | 0/200 | 148/200 | 22/200 | 2/200 | 124/200 |
| 2 | 1.2 | 0/200 | 157/200 | 33/200 | 3/200 | 121/200 |

## Distinct outcomes and verification

Essential shortfall means failure to deliver the specified productive quota for the modeled essential class. It is not cell death. Degradation and finite-horizon sequestration can reduce unresolved influx but cannot satisfy that quota. Unserved flux is residual proteotoxic influx, not a measured proteotoxic-failure threshold; no aggregation dynamics or growth law is fitted. ATP is a synthetic route cost, not a measured cellular energy budget. Its minimization is fourth in the lexicographic objective, so it cannot replace productive rescue with cheaper disposal.

| Architecture arm | Degradation q05/q50/q95 | Sequestration q05/q50/q95 | Unserved q05/q50/q95 | ATP q05/q50/q95 |
|---|---|---|---|---|
| baseline | 0/0/0 | 0/0/0 | 0/0/0 | 0.000300501/0.00652253/0.148092 |
| matched_GroEL | 0/0/0 | 0/0/0 | 0/0/0 | 0.000300501/0.00652253/0.148092 |
| matched_DnaK | 0/0/0 | 0/0/0 | 0/0/0 | 0.000300501/0.00652253/0.148092 |
| matched_ClpB | 0/0/0 | 0/9.6209e-05/0.00123462 | 0/0/0 | 0.000199323/0.00560209/0.147875 |
| KO_TF_DnaK | 0/0/0.0497351 | 0.000225038/0.00989025/0.0686751 | 0/0/0.0563861 | 0.00340632/0.0365358/0.568437 |

Every route-allocation optimization checks primal feasibility, dual feasibility, nonnegative prices and strong duality. Maximum absolute certificate residual across all saved evaluations: **9.9e-12**. Ensemble assertions pass: scalar target invariance, equal remaining capacities and nested essential optima (secondary optima checked when primary quotas coincide).

The supported mathematical result is an explicit limitation of scalar reductions; the biological numerical ranking is not established. No pass is awarded for matching the observations used to specify routes. Required next evidence: client-resolved synthesis/turnover demands, productive route qualification, machine occupancy/service times, and held-out perturbation outcomes. Kerner client counts and capacity occupancy, Calloni interactions, and Babu genetic interactions are not full rate constraints.

Two direct confrontations with existing experiments — the imposed serial TF/DnaK requirement versus GroEL/ES overproduction rescue of the double knockout, and the effect of all-client radial normalization — are in [CONFRONTATION.md](CONFRONTATION.md), produced by `confront.py`.

See [THEORY.md](THEORY.md) for proofs, [README.md](README.md) for exact reproduction and scope, and [SOURCES.md](SOURCES.md) for provenance. Raw outputs and input draws are in `outputs/`.
