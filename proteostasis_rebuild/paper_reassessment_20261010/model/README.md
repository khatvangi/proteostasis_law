# Route-capacity reassessment

This directory is self-contained apart from a read-only input, `../../WO-07/bundles.tsv`. It preserves previous runs and `WORK_ORDERS.md`.

Run from the repository root:

```bash
python -m unittest discover -s paper_reassessment_20261010/model -v
python paper_reassessment_20261010/model/run.py --trials 200 --seed 20261010
python paper_reassessment_20261010/model/report.py
python paper_reassessment_20261010/model/validate.py
```

Dependencies: Python 3, NumPy and SciPy (`scipy.optimize.linprog`, HiGHS). Exact versions and input/code SHA256 values are recorded in `outputs/manifest.json`. A small diagnostic run writes elsewhere:

```bash
python paper_reassessment_20261010/model/run.py --trials 3 --out /tmp/route-smoke
```

The default run uses 200 independent trial identities, 28 capacity scenarios, three nested models, eight architecture-ablation evaluations per trial, and a 2×3×4 catalytic-width/load/matched-perturbation sensitivity grid. Each evaluation uses four LPs. Runs are single process; on this environment the entire suite is manageable on a workstation. `inputs.json` stores every primary demand-share vector, capacity, client-specific work matrix, multiplier, and radial margin. CSV files retain every trial, not only summaries. The report is generated from these CSVs and runs ensemble invariants before publishing. An alternative `--out` writes its report there rather than replacing the primary report.

## Explicit uncertainty design

- Obligate subset flux g is log-uniform on [0.005,0.35]. Within it, the essential-obligate flux fraction is uniform [0.02,0.4]. Neither range is inferred from client counts. Only this subset is GroEL-obligate; the rest is never made hard-obligate by a random route generator.
- Aggregate fraction is log-uniform [0.001,0.1]; overlap fraction is uniform [0.05,0.45]; remaining demand goes to other clients. All shares are positive and sum to one. This range ensures a restricted subset, not a claim that 35% is a measured upper bound.
- Effective catalytic multipliers are independent log-uniform [0.1,10]; sensitivity uses [0.01,100]. Multipliers include unknown cycle/occupancy conversion; TF throughput is a residence-time rate, not ATPase catalysis.
- Seven measured pool×multiplier values are normalized to total machine work capacity 1. Generic protease and finite-horizon storage each have capacity 0.05. They remain separate physical resources in the full LP and explicit optimistic pooled resources in its relaxations. No measured cellular units or ATP budget is inferred after normalization.
- Each machine/client service cost is independently log-uniform [0.5,2], held fixed across that client's alternative routes and all counterfactual models. This permits client-specific demand for different resources while keeping routes deterministic.
- Primary total demand is 0.9 times the architecture's all-client folding margin. Sensitivity uses 0.6, 0.9 and 1.2. Recalibration is only to the **unperturbed** architecture at each catalytic width. It removes an unknown global load scale and therefore prevents baseline feasibility from being counted as validation.
- Illustrative stage ATP coefficients: TF 0, DnaK 1, DnaJ/GrpE 0 additional, GroEL 7, GroES 0 additional, ClpB 2, protease 4, storage 0. These are accounting examples, not biological estimates. ATP route cost is the coefficient-weighted work sum. Changing them can select among equally productive/handled allocations but cannot change either preceding lexicographic optimum.

`THEORY.md` gives the LP, exact fractional dual, ordinary min-cut special case, nesting proof and scalar non-identification theorem. `SOURCES.md` separates measured inputs from assumptions. `RESULTS.md` is the numerical assessment and failed empirical-prediction gate. No distribution here is a posterior or a sample of observed cells.

## Validation scope

Tests include the small-demand tolerance regression caught by the first ensemble run (rejected outputs retained in `audit_default_tolerance/`), and cover serial AND versus OR, a hand-computable Hall cut, dual shadow prices, finite sinks, productive essential qualification, zero-resource failure, matched doses, nested feasible sets, invariant machine relabeling, seeded draws, and a constructive non-identification witness. Every LP additionally checks its certificate. Byte-identical small reruns check pipeline reproducibility on the recorded environment; other SciPy/HiGHS versions may choose different allocations at degenerate optima.

The generic solver accepts essential constraints through extra lower-bound rows; `evaluate` specifically uses the first client class as the essential aggregate. Do not silently apply that convenience routine to multiple distinct essential classes. Splitting the 13 essential proteins requires their separate demands and productive quotas.

The shared append-only programme-map ledger is outside the writable session paths. `LEDGER_PENDING.md` supplies the exact pending entry; it has not been appended to the central ledger.

## Reviewable route inventory

KJE means DnaK with its DnaJ and GrpE cofactors; GRO means GroEL with GroES. Each member consumes its own capacity. Members of a cofactor module are not asserted to be consecutive folding steps. Internal route names list consumed resources, while arrows below indicate the hypothesized stage order.

| Client class | Qualified productive alternatives | Other terminal routes |
|---|---|---|
| Essential obligate subset | TF → GRO **or** KJE → GRO | Generic degradation or finite storage; neither satisfies essential quota |
| Other obligate subset | TF → GRO **or** KJE → GRO | Generic degradation or finite storage |
| Upstream overlap | TF **or** KJE | Generic degradation or finite storage |
| Aggregate | ClpB with KJE recovery | Generic degradation or finite storage |
| Other clients | TF **or** KJE **or** GRO | Generic degradation or finite storage |

The mandatory assistance of the `other` class is a benchmark boundary: it denotes modeled machine-requiring clients, not every protein synthesized by the cell. Spontaneously folding, secreted and otherwise unmodeled proteins are outside the normalized demand. Thus the sampled obligate share is a fraction of modeled client influx, not a measured proteome-wide fraction. Distinct state routes can be added by extending this inventory, but every new qualified itinerary requires evidence or an explicit sensitivity assumption.
