# network_failure_test — results and adversarial audit

**AUDIT 2026-10-09.** Everything below is new work. The original code (`model.py`,
`pools.py`, `feedback.py`, `run_all.py`, `test_network.py`) and the original outputs in
`results/` are unchanged. The audit code and outputs are in `audit/`.

This is a structural model. Machine pools are measured (Schmidt 2016, glucose exponential,
via WO-07). Every per-unit throughput, every client–machine compatibility and the demand
mix are assumptions. Nothing here is biological evidence. The output says what follows
from the assumptions and which assumptions carry each result.

---

## Bottom line

1. **The 98.7% single-cut statistic comes from the route rules alone.** Probabilities
   chosen by hand fix it before any capacity, abundance or load enters. The exact
   probability from the rules is 97.99%. A 20,000-draw Monte Carlo gives 97.97%. The
   original 98.7% is the same quantity on 1,000 draws. The main driver is not ClpB but
   the GroEL/GroES class, which gets an alternative route only 20% of the time. The next
   drivers are DnaK/J/E and HtpG. Guaranteeing aggregates a ClpB-free route changes the
   rate by 0.5 points. Add one disjoint backup route per class and it falls to 0%.
   **The statistic says nothing about biology.**
2. **Most synthetic-lethal pairs are structural.** Roughly two-thirds come from route
   topology: knocking out both machines leaves some class with no route at all, at any
   load. For example, nascent chains get only TF and KJE by rule. The load-dependent
   remainder is real in the model. It survives removing every forced single cut and
   **increases** as routes are added. As a generic property of capacity-sharing networks
   it is robust. As a property of E. coli it is untested.
3. **Synthetic lethality does not tell a network apart from a pool.** The pooled-bucket
   strawman (A) shows synthetic-lethal pairs in 0–100% of instances by setting (100% at W=0, 33% at W=2). An unstructured
   random network shows more than the biology-ruled one (79% vs 57%).
4. **Multi-stage saturation is mostly an artefact of the route-preference assumption.**
   The saturation events never move the viability threshold, which is a theorem of the
   linear program. Per network, the number of events matches the original in only 34–49%
   of cases when route preference is changed. The first machine to saturate matches only
   30–47% of the time. Machines saturate in sequence in most networks with alternative
   routes. Which machines, and how many, is set by unmeasured quantities.
5. **The feedback experiment did not produce multiple thresholds.** None of the 772
   converged feedback runs has two or more jumps; the pre-set criterion was 25%. Single
   jumps appear in 2.0% (titration) and 0.5% (utilization) of converged runs. They need
   loop gain γβ ≫ 1. Jumps are 20–36× more common among runs that **failed to converge**,
   so those should not be counted.
6. **The Pareto trade-off holds within the model, but the published example
   understates it, and the ATP axis is invented.** Across 300 networks, 75% show an ATP
   versus reserve trade-off. The median ATP span is 24%, and a median of 34 of 40 front
   points are non-dominated. The saved example (seed 1) has only 2 non-dominated points
   and a 5.9% span. ATP weights are illustrative orderings, not measurements, and the
   "toxicity" axis is only reported, never optimized. So this is two objectives plus one
   observed quantity, not three objectives.
7. **The numbers reproduce exactly.** I reran 120 instances across 4 settings, 24
   feedback runs and the Pareto front. Every column matches the stored results to within
   6e-13. The 6 failed instances, all at α=0.3, are a solver bug, not a model outcome.
   HiGHS presolve wrongly reports "infeasible" when a demand share is below about 1e-7.
   With presolve off, all 6 solve. The original fallback loop could not catch this
   because all three HiGHS methods share the same presolve.

---

## 1. Code audit

| Item | Finding |
|---|---|
| Solver | 6/9,500 instances dropped (`bio_measured_W2_alpha0.3`, seeds 9, 65, 286, 370, 379, 416). The max-λ LP always has the feasible point x=0, λ=0, so "infeasible" is a solver error. All 6 solve with `presolve=False` (`audit/results/failed_instances_resolved.csv`). The bias from the drop is small (≤1.2% of that setting) but not random: it removes networks where one class has near-zero demand. The fallback counter reported 0 because failed rows are excluded before counting. |
| Silent errors | None found in the 9,494 successful rows. Zero-capacity knockouts equal the structural-cut count in 100% of rows. No knockdown exceeds the intact network. No 50% knockdown falls below 0.5. Intact λ*_C > 0 always. Ordering A ≥ C ≥ B holds in every row. |
| Replicates | All W and abundance settings reuse the **same 1,000 route sets**: seed s gives the same routes because routes are drawn before capacities. "0.987 in every setting" is one measurement repeated, not robustness. |
| A theorem presented as a result | "No synthetic lethality under 50% knockdown at ρ=0.5" (0.0 in every setting) is guaranteed. Halving any set of machines cannot push λ* below half, because λ*(c′) ≥ λ*(κc) = κλ*. It carries no information. |
| Feedback docstring | `feedback.py:25` says u_src = usage/c_eff. The code at line 66 divides by intact capacity. This affects interpretation, not the reported numbers. |
| Feedback couplings | Utilization feedback couples random machine pairs (e.g. HslV→TF). There is no mechanism behind it; it is a stress test of the solver, not a model of anything. |
| Pareto selection | `findParetoInstance` returns the first seed with a sequestration route and any ATP span. That makes it an existence example, not a representative one. |
| Tests | Original suite: 19/19 pass. The tests check LP correctness (Hall min-cut, duality, A ≥ C ≥ B, theorems). None tests whether a result depends on the generator's probabilities. |

## 2. What drives the single-cut statistic

A machine is a structural single cut when every route of some client class passes through
it. Losing it makes λ*_C = 0 at any load and any capacity. The statistic therefore depends
only on the route list.

Analytically, a network has no single cut only if all of the following hold. The
probabilities multiply because each rule uses its own random draw:

| Condition (original rules) | Probability |
|---|---|
| GRO class gets its KJE alternative | 0.20 |
| HTP primary is HtpG alone **and** a KJE alternative exists | 0.25 |
| DK class gets a GRO or TF alternative | 0.65 |
| DMG class gets ≥ 2 proteases | 0.784 |
| AGG class gets SEQ or a Lon-only route (ClpB-free) | 0.79 |
| **all five** | **0.0201 → 97.99% have a cut** |

Per-machine frequency of being a single cut (20,000 draws, `audit/results/struct_cut_machines.csv`):
GroEL/GroES 79.7%, DnaK/DnaJ/GrpE 74.3%, HtpG 50.7%, ClpB 20.8%, each protease 7%,
TF/IbpA/IbpB 0%.

Counterfactual architectures, all capacity-free (`struct_single_cuts.csv`, n=20,000 each):

| Architecture | Networks with a single cut |
|---|---|
| original rules | 97.97% |
| original + guaranteed ClpB-free aggregate route | 97.41% |
| original + one disjoint backup route per class ("nocut") | 0% (by construction) |
| unstructured random null (1–4 routes, 1–3 machines) | 90.5% |
| uniform optional-edge density q = 0 / 0.2 / 0.3 / 0.5 / 0.7 / 0.9 / 1.0 | 100 / 99.8 / 98.5 / 85.0 / 52.6 / 13.5 / 0% |

Reading: any sparse compatibility matrix produces single cuts. The cut rate is a dial set
by route density. Only above roughly 70% edge density does the rate fall below half.

## 3. Synthetic lethality: structural versus load-dependent

Architecture C, full knockouts, n=300 per setting, W=2, measured abundance unless stated
(`audit/results/cf_summary.csv`). The kinetic multipliers come from a separate random
stream, so capacity can be resampled with the structure held fixed.

| Setting | single cut | structural SL pair (ρ=0.5) | load-dependent SL pair, ρ = 0.25 / 0.5 / 0.9 |
|---|---|---|---|
| original rules | 99.7% | 45% | 21 / 21 / 17% |
| same structures, capacities resampled | 99.7% | 45% | 19 / 17 / 17% |
| original, W=0 | 99.7% | 41% | 21 / 18 / 18% |
| original, W=4 | 99.7% | 49% | 19 / 23 / 19% |
| original, abundance blind | 99.7% | 69% | 12 / 20 / 19% |
| original, abundance shuffled, W=0 | 99.7% | 60% | 15 / 21 / 19% |
| ClpB-free aggregate route | 98.7% | 52% | 22 / 21 / 19% |
| aggregate share fixed at 1% | 99.7% | 62% | 19 / 22 / 20% |
| **nocut** (no single cuts) | 0% | 84% | 26 / 29 / 19% |
| nocut, aggregate share 1% | 0% | 80% | 33 / 30 / 23% |
| nocut, abundance blind | 0% | 87% | 31 / 34 / 30% |
| density q=0.2 | 99.3% | 35% | 6 / 7 / 3% |
| density q=0.5 | 89.0% | 64% | 20 / 21 / 18% |
| density q=0.8 | 31.0% | 62% | 47 / 44 / 37% |
| density q=1.0 | 0% | 55% | 79 / 59 / 47% |
| unstructured null | 92.3% | 67% | 35 / 44 / 43% |

Notes: the original ρ=0.5 "any SL pair" figure of 56.7% splits into about 45% structural
plus 21% load-dependent, with overlap. A structural pair is lethal at any load, so its
presence at a given ρ is not load-dependent. Removing single cuts raises the structural
pair rate, because cuts move from size 1 to size 2.

Verdict: **load-dependent synthetic lethality is robust in this model class.** It does not
need forced edges and it grows with redundancy, which is what one expects when a dense
network of shared capacity loses two machines that each absorb the other's overflow. Its
frequency is set by route density, not by anything measured.

## 4. Bottlenecks: abundance versus kinetics

ClpB as the dominant bottleneck (54% of networks at W=2, 73% at W=0) comes from three
things together. Its measured pool is the smallest (0.0117 µM ClpB6). It is on every
default aggregate route. And the uniform Dirichlet demand gives aggregates about 1/6 of
all client flux. Take any one away and it collapses:

- flat abundance: 4%
- shuffled abundance: 7%
- aggregates at 1% of demand: 13%
- W=4 (kinetics dominate abundance): 30%

With structure fixed and only capacities resampled, the bottleneck set is unchanged in
only 43% of networks. **The bottleneck identity is a product of measured abundance and an
assumed demand mix; it is not a structural result.**

## 5. Multi-stage saturation

Below λ*, machines saturate one after another as overflow routes take up the load.
Viability is lost only at λ*: without feedback the feasible loads form one interval
[0, λ*] (original theorem, re-tested in `audit/test_audit.py`). With the same structures:

- route preference reversed: event count unchanged in 34% of networks, first saturating
  machine unchanged in 30%
- flat preference: 49% and 47%
- capacities resampled: 49% and 56%

The share of networks with ≥ 2 events moves 45 → 48 → 57% across the three preference
orders. It rises with route density: 20% at q=0.2, 100% at q=1.0.

Verdict: sequential saturation is generic to redundant networks. The specific sequence is
an artefact of the unmeasured preference order and throughputs. It is not a set of
failure thresholds.

## 6. Feedback experiment (re-analysis of `results/feedback.csv`)

| Feedback kind | runs | non-converged runs | ≥2 jumps (converged) | ≥1 jump (converged) | ≥1 jump among non-converged |
|---|---|---|---|---|---|
| none | 50 | 0 | 0% | 0% | — |
| titration | 400 | 43 | 0% | 2.0% | 39.5% |
| utilization | 400 | 35 | 0% | 0.5% | 20.0% |

Jumping runs have median log10(γβ) of 1.74 (titration) and 0.56 (utilization). Non-jumping
runs have 0.31 and −0.42. The jumps are runaway positive feedback with loop gain well above
1, the same mechanism as the hand-built positive control. Titration feedback cannot move
onset below λ*, because it acts only once load goes unserved. Utilization feedback moves it
earlier in 5% of runs, by construction.

The pre-set criterion was ≥2 jumps in ≥25% of converged runs. **It fails at 0%.**

## 7. Pareto front

The trade-off is computed by ε-constraint at load 0.5 λ*: minimize ATP subject to maximum
utilization ≤ t. Ensemble of 300 W=2 networks (`audit/results/pareto_ensemble.csv`):

- 25% degenerate: no trade-off at all
- non-degenerate: median ATP span 24%; 95th percentile 184%
- median 34/40 points non-dominated; the ATP floor is typically reached only at t=1, so
  the trade-off runs across the whole range
- sequestered flux changes along the front in 53% of networks

The stored seed-1 front: 2 non-dominated points, 5.9% span, and ATP flat from t=0.513 to 1.

Status: the model has a genuine two-objective trade-off between ATP and reserve. The ATP
costs are invented orderings (`ATP_ILLUSTRATIVE`), so its size means nothing
quantitatively. Toxicity is a passenger, so "multi-objective" in a three-axis sense is not
supported. Reserve margin at fixed load is 1 − λ/λ*, so it is not independent of λ*
(original test `test_reserve_margin_is_lambda_star`).

## 8. Verdict on the structural claim

| Claim | Status | Why |
|---|---|---|
| 98.7% of networks have a structural single-machine cut | **FAIL as a finding** (PASS as arithmetic) | Reproduces exactly, but it is a property of the hand-set edge probabilities (97.99% analytic). It falls to 0% with one backup route per class, and ClpB is not the driver. |
| Load-dependent synthetic lethal pairs exist | **PASS (model-generic)** | Survives removing forced cuts, ClpB exclusivity, abundance and kinetic changes. It grows with redundancy. Frequency set by assumed density. |
| Synthetic lethality distinguishes network from pool | **FAIL** | The pooled bucket also shows synthetic-lethal pairs (0–100% by setting; 100% at W=0). |
| ClpB/IbpB are the binding bottlenecks | **FAIL as structural** | Abundance × assumed aggregate demand share; 54% → 4–13% under counterfactuals. |
| Multiple saturation stages before failure | **PASS (soft) / FAIL as thresholds** | Generic, but the order and count depend on the unmeasured preference order, and they never move viability. |
| Multiple cascade thresholds under feedback (N3) | **FAIL** | 0% of converged runs; single jumps need loop gain ≫ 1. |
| Pooled capacity is a poor summary (N1) | **PASS, uninformative** | Guaranteed whenever pools span 3 decades and compatibility is sparse. A is a strawman. |
| Pareto trade-off | **PASS within model, illustrative only** | Real in 75% of networks; the ATP axis is unmeasured. |
| Reproducibility of reported numbers | **PASS** | 120 instances, 24 feedback runs and the Pareto front match exactly. |
| Completeness | **INCOMPLETE** | 6 instances dropped by a solver bug (fix: presolve off). Not rerun into the original results, which stay intact. The counterfactual ensembles are 300 networks each (binomial SE about 3 points), not 1,000. |

**Overall: the structural claim, that E. coli-like proteostasis fails through
distributed single cuts and synthetic-lethal pairs rather than one pooled threshold, is
BLOCKED on data.** The model can produce either picture. Which one it produces depends on
compatibility density, which is unmeasured.

## 9. Data that would discriminate

In order of leverage:

1. **The client–machine compatibility matrix.** For each client class, which machines can
   independently resolve it in vivo. This one input sets the single-cut rate (0–100%) and
   the synthetic-lethal frequency. Client overlap between DnaK, GroEL and TF, measured in
   cells lacking each machine, would turn the density knob q into a measured number.
2. **Pairwise genetic interactions among the 15 machines at graded expression.** A
   CRISPRi titration matrix of single and double knockdowns, with growth as readout, at
   two or more load levels (e.g. temperature, or a mistranslating background). The model
   predicts load-dependent synthetic lethality: pairs viable at low load and lethal at
   high load. A plain pool predicts the weakest link only. Full deletions cannot separate
   these, because structural pairs are lethal at any load.
3. **Single-gene essentiality of the 15 machines under the matched condition.** This
   is a cheap first check of the single-cut prediction. Under the original rules the
   model calls HtpG a single cut in 51% of networks and ClpB in 21%. Deletion-collection
   essentiality calls should be retrieved and compared before anything else. They were
   not retrieved in this audit.
4. **In vivo per-unit throughput (s_m) for at least ClpB, DnaJ and IbpB,** the
   machines the model most often finds binding. Without these, the bottleneck ranking is
   abundance only.
5. **Measured aggregate share of client flux in exponential growth.** Changing it from
   about 17% (uniform Dirichlet) to 1% moves ClpB from the main bottleneck to a minor one.

## Reproduce

```
cd network_failure_test
PYTHONDONTWRITEBYTECODE=1 python -m unittest test_network          # original, 19 tests
cd audit
PYTHONDONTWRITEBYTECODE=1 python -m unittest test_audit -v           # audit, 10 tests, ~20 s
python run_audit.py struct feedback failed                           # ~1.5 min, local
sbatch run_audit.sbatch                                               # pareto, cf, repro; ~15 min on 16 cores
```

Every part is deterministic in its seeds. Outputs are in `audit/results/`, with the
headline numbers in `audit_summary.json`. Software: scipy 1.14.0, numpy 2.2.6, HiGHS via
`scipy.optimize.linprog`.
