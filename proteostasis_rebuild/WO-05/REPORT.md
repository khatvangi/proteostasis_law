# WO-05 — error-to-burden semantics (reopened 2026-09-29)

## Verdict

**VERDICT: PASS** — G5.1–G5.5 are met. The earlier BLOCKED verdict (2026-09-28)
is withdrawn. Its *inference* was right: Data_S2 means are conditional on
detection. Its *reason* was wrong: it said the zero cells and denominators
were unavailable, but they are published in Landerer Data_S4, which the
earlier run never opened. With them, G5.5 can be computed, so it has been.

PASS certifies only that the gates were met. It does **not** establish an
E. coli headroom (see "What the number means" and the self-review).

## Sources (accessed 2026-09-29, all read-only)

| Source | Identifier | What was used |
|---|---|---|
| Landerer, Poehls, Toth-Petroczy 2024, *Mol Biol Evol* 41:msae048 | doi:10.1093/molbev/msae048, PMC10939442 | Methods definition, dataset count (80 E. coli), Fig. 1b legend; local full-text XML |
| Landerer supplement `msae048_supplementary_data.zip` | sha256 `c2b7cef8…c9c416`; local copy at `triplet-proof/reviewer_response/position_rates/raw/landerer2024/` | Data_S2 (byte-identical to the legacy `envelope-paper/data/raw` copy, sha256 `08ef9455…266b`), Data_S4 per-dataset `codon_counts.csv` / `substitution_errors.csv` |
| deTEL code, `git.mpi-cbg.de/tothpetroczylab/detel`, main @ `c3593d46` (2026-03-05) | cited in the paper's Code Availability | `eTEL/workflow/global_report.py:get_all_codon_count` (sha256 `989e0e0c…2b9080c5a`) |
| Stikeleather, Ali, Ho, Licknack, Lynch 2026, *Nucleic Acids Res* 54(13) | doi:10.1093/nar/gkag674, PMID 42406629, PMC13335486 | WT rate, definition, culture conditions; PMC full text via PubMed |

## G5.1 — three quantities, units, and what n means

All rates are dimensionless probabilities per codon translated. Burden flux is
a fraction of the modeled pool per second.

* **Raw decoding error** `e_raw`: the probability of a non-cognate outcome,
  synonymous outcomes included. MS never measures it. Under the legacy
  uniform `S`, `f_sub = e_raw (1-S)`.
* **Amino-acid substitution rate** `f_sub`: the probability per codon that
  the amino acid changes, summed over destinations.
* **MS detection rate** `r_MS`: per codon and per dataset, substitution-carrying
  PSMs covering the codon divided by all PSMs covering it (Landerer Methods).
  It is at amino-acid level, so there is no synonymous component. I/L and
  PTM-masked pairs are invisible and rare events are missed, which makes it a
  **lower bound on f_sub, not raw error**.

**How Data_S2 aggregates the 80 datasets** was settled from the raw
per-dataset counts, not from `(sd/se)^2`:

* Data_S4 has, for every one of the 80 E. coli datasets, `base_count`
  (covering PSMs) and `error_count` per codon. `detection_rate = error/base`
  holds to 1e-16. `error_count` equals the number of substitution-PSM rows for
  every (dataset, codon), with 0 mismatches. The total is 58,990 substitution
  PSMs, matching the paper's text.
* Cell census (61 sense codons × 80 datasets = 4,880 cells):
  - **19 no coverage** (codon row absent).
  - **3,015 covered with zero detected substitutions.**
  - 1,846 covered with ≥1 detected substitution.
* Data_S2 is reproduced **only** by the detected-only rule
  (`error_count > 0`). Under that rule the maximum deviations are 2.6e-9 (mean),
  4.2e-9 (sd) and 3.9e-10 (median), and `n = (sd/se)^2` equals the count for
  all 61 codons. The zero-inclusive rule misses by up to 1.8e-2 in the mean and
  its n never matches.
* **n is the number of datasets with ≥1 detected substitution at that codon.**
  It is not the number of covering datasets: almost every dataset covers
  almost every codon. For example, CTG, E. coli's most-used codon, has n = 52.
* Upstream mechanism: deTEL's `get_all_codon_count` maps `detection_rate`
  through `log10`, replaces `-inf` (rate 0) with NaN, and calls `dropna`, so
  covered zero-detection cells leave the per-dataset distribution (the
  Fig. 1b box plots are on a log axis). The script that wrote Data_S2 itself is
  not in the repository. The exact numerical reproduction above is what
  establishes the rule; the code shows how it arises.

So the prompt's alternative reading does not hold for Data_S2: that n < 80
reflects coverage gaps and covered zero-error datasets enter as 0. Covered
zeros are **excluded**. The paper's note that "not all substitutions are
present in all datasets" describes detection, not coverage.

## G5.2 — legacy mapping inventory (unchanged, 21/21 anchors found)

| File:line | Class | Reason |
|---|---|---|
| `envelope-paper/scripts/11_headroom_sensitivity.py:81` | double-discount | MS `mu` multiplied by `(1-S)` |
| `envelope-paper/scripts/12_chaperone_availability.py:72` | double-discount | MS `mu` multiplied by `(1-S)` |
| `envelope-paper/scripts/09_supraadditivity.py:99` | double-discount | MS-derived `f` is multiplied by `(1-S)` |
| `envelope-paper/scripts/vendor/two_pool_ode.py:262` | correct | Converts critical flux to a raw-error threshold |
| `envelope-paper/scripts/vendor/two_pool_ode.py:513` | ambiguous | `1e-4` input's raw-versus-substitution semantics are unstated |
| `proteostasis-P1/two_pool_ode.py:262` | correct | Converts critical flux to a raw-error threshold |
| `proteostasis-P1/two_pool_ode.py:513` | ambiguous | `1e-4` input's semantics are unstated |
| `proteostasis-P1/two_pool_ode.backup_uniformN.py:236` | correct | Converts critical flux to a raw-error threshold |
| `proteostasis-P1/two_pool_ode.backup_uniformN.py:483` | ambiguous | Crosscheck `1e-4` input's semantics are unstated |
| `proteostasis-P1/paired_mc.py:146` | ambiguous | `f_obs=1e-4` is mapped as raw without definition |
| `proteostasis-P1/paired_mc.py:148` | ambiguous | Same undefined `f_obs` is converted to effective damage |
| `proteostasis-P1/arithmetic_stress_test.py:72` | correct | Exact threshold on raw error |
| `proteostasis-P1/arithmetic_stress_test.py:80` | correct | Large-N threshold on raw error |
| `proteostasis-P1/arithmetic_stress_test.py:113` | correct | Threshold on raw error |
| `proteostasis-P1/arithmetic_stress_test.py:202` | ambiguous | Sets `(1-S)*p_m=1`; no longer the stated raw or substitution threshold |
| `proteostasis-P1/arithmetic_stress_test.py:242` | correct | Converts effective-damage threshold back to raw error |
| `proteostasis-P1/essential_bound.py:209` | correct | Threshold on raw error |
| `proteostasis-P1/essential_bound.py:216` | correct | Raw error to effective damage |
| `proteostasis-P1/essential_bound.py:241` | correct | Threshold on raw error |
| `proteostasis-P1/figures/fig2_arithmetic.py:114` | correct | Raw-error threshold; same expression at lines 109, 121, 123 |
| `envelope-paper/scripts/06_translation_burden.py:37` | correct | Uses MS `mu` directly as a substitution-level lower bound, no extra `(1-S)` |

## G5.3 — corrected mapping and guard

`error_semantics.py` types every rate (`RAW`, `SUB`, `MS`) and applies
`(1-S)` exactly once, to `RAW` only. The guarded legacy conversion raises
`TypeError` for `SUB` and `MS` inputs, and a unit test asserts this.

`flux()` now also takes `balanced_growth=True`, which uses the WO-01
per-protein synthesis rate `ln2/T_gen`. The default is `False` (legacy
`1/T_gen`), so each correction can be isolated.

## G5.4 — arithmetic threshold (unchanged)

The exact raw threshold at N=300, P=0.70, S=0.30, p_m=0.30 is
`0.005658142850513878`. The substitution-level threshold is
`0.003960699995359715`. The legacy 1.19e-3 is `-ln(0.70)/300`, which comes
from forcing `(1-S)p_m = 1` (`arithmetic_stress_test.py:202`).

## G5.5 — headroom recomputed

Model: the vendored legacy `two_pool_ode.py` (sha256 unchanged), with the
as_published anchoring (C_tot 50 uM, K_d 1 uM). `headroom_P = P_dagger/P*(J)`
is nonlinear in J, so every value comes from the steady-state solver, not
from rescaling. The old ×24.817 is reproduced exactly.

**Effect of each correction on the old ×24.817** (same input, the Data_S2
conditional mean usage-weighted to 6.334e-4):

| Mapping | J factor vs legacy | headroom_P |
|---|---|---|
| legacy: `(1-S)`, `1/T_gen` | 1 | **×24.817** |
| remove the `(1-S)` double discount only | 1/0.7 = 1.429 | **×17.320** |
| `ln2/T_gen` synthesis only | 0.693 | **×35.876** |
| both corrections | ln2/0.7 = 0.990 | **×25.064** |

The two mapping corrections cancel to within 1%. The ×25 "survives" the
mapping fixes only by that coincidence.

**Effect of correcting the empirical input** (both mapping corrections
applied). All inputs are E. coli eTEL aggregate statistics: MS detection rates
at substitution level, usage-weighted with the legacy genome codon counts.

| Empirical input | f (/codon) | headroom_P | headroom_A | linear f-margin |
|---|---|---|---|---|
| Data_S2 conditional mean (legacy mu; selected on detection) | 6.334e-4 | ×25.06 | ×280 | ×15.9 |
| **zero-inclusive mean over covering datasets (primary)** | 2.526e-4 | **×63.10** | ×1769 | ×40.0 |
| PSM-pooled Σerr/Σbase per codon | 1.414e-4 | ×112.87 | ×5657 | ×71.5 |
| median of the 80 per-dataset usage-weighted rates | 5.859e-5 | ×272.5 | ×32974 | ×172 |

The corrected critical substitution rate is f_crit = 1.0102e-2 /codon. It is
set by the imposed A_max gate (`aggregation_death`), not a fold (WO-04).
The global PSM-pooled rate over all codons is 1.193e-4.

**Per-dataset spread.** Each of the 80 datasets has its own usage-weighted
rate, from 6.7e-6 to 3.29e-3 (IQR 1.99e-5 to 2.00e-4). With both mapping
corrections, the headroom ranges **×2370 (min) / ×805 (q25) / ×273 (median) /
×80 (q75) / ×4.6 (max)**. Rate and PSM depth are negatively correlated across
datasets (Spearman −0.37).

### What the number means

It is an MS error-*detection* rate at substitution level, a lower bound on
`f_sub`. It is not raw decoding error. It pools 80 PRIDE datasets that differ
in strain, medium, growth phase and pipeline depth. It is therefore an **eTEL
aggregate statistic, not an E. coli physiological operating point**. Every
headroom above is conditional on the legacy model, which has a non-conservative
Phi inflow (WO-02), no dilution (WO-01) and an imposed A_max gate (WO-04).
Adding the missing dilution sink would further raise the headroom at any
input. That was not done here, because it belongs to the conservative model,
not to G5.5.

### Stikeleather 2026 — a separate, condition-specific estimate

* Wild-type E. coli mean translation-error rate **1.82e-3 /codon**, defined as
  total detected substitutions / total sites sampled.
  - The method is MS-detected, with I/L merged and chemical-artefact
    substitutions removed, so it is the same *kind* of quantity (`MS`,
    substitution level).
  - The abstract's "2 per 1000 amino acids" confirms the exponent.
  - SE 5.92e-5: the PMC text extraction lost the exponent. The mantissa
    matches, but the exponent is taken from the task statement and is
    **UNVERIFIED**.
* Condition: Xac-derived strain, LB (Miller), 37 °C, overnight, harvested in
  **stationary phase**, n = 3 biological replicates.
* The rate is 7.2× the zero-inclusive eTEL aggregate and 15.3× the global
  PSM-pooled eTEL rate. The authors attribute the difference to aggregation
  across 80 heterogeneous datasets.
* **No headroom is computed from it.** The legacy parameters (T_gen = 3600 s
  synthesis, chaperone pool) are exponential-growth parameters, and
  combining the two conditions is exactly what WO-07 forbids without a
  MISMATCHED flag. For orientation only, flagged MISMATCHED:
  f_crit/1.82e-3 = 5.55 on the linear f axis. This is not a stationary-phase
  margin.

## Adversarial self-review

1. *Did the reopening just move the goalposts to get a PASS?* No. G5.5 asks
   only that the ×25 be recomputed with the double discount removed. That is
   ×17.32 alone, or ×25.06 with the WO-01 ln2 correction. The earlier block
   rested on a factual error about data availability, now shown false:
   Data_S4 is in the same supplement as Data_S2. The earlier *semantic*
   inference, that Data_S2 is conditional on detection, is confirmed, not
   reversed.
2. *Does the old ×25 survive?* Numerically yes, at the same input: ×25.06.
   But this rests on a 1% cancellation between two independent errors, and
   the input is selection-biased upward by about 2.5×. As a statement about
   E. coli it is **not supported**. On the same eTEL corpus the corrected
   headroom runs from ×4.6 to ×2370 depending on the dataset, and from ×25 to
   ×273 depending on the aggregation rule. A single "×25", or "one order
   inside", is not a measured property of E. coli.
3. *Which direction is the truth likely to lie?* It is unresolved and cuts
   both ways.
   - Correcting selection bias lowers f and raises the headroom.
   - MS undercounts (I/L, PTM-masked, missed rare events), which raises the
     true `f_sub` and lowers the headroom.
   - The best-controlled single study (Stikeleather) gives a rate 7–15×
     higher than the eTEL aggregates, but in stationary phase.
   - The negative rate–depth correlation (−0.37) is unexplained. It could be
     heterogeneous biology or thin-sampling inflation, and the data here
     cannot separate them.
4. *Is the zero-inclusive mean the right estimand?* It gives equal weight to
   each covering dataset. The pooled ratio gives weight by PSM and is
   dominated by a few deep datasets: PXD002140 alone has 50M PSMs at about
   1e-5. Both are reported. Neither is a physiological operating point.
5. *Could the reconstruction match by accident?* The detected-only rule
   reproduces mean, sd, median and n for all 61 codons to about 1e-9 at once.
   The alternative misses by orders of magnitude more than that. The same
   holds for the yeast sheet (72 datasets; checked during analysis, not
   gated).
6. *Read-only?* Yes. All 30 legacy hashes from WO-00 re-verify unchanged, the
   WO-00 provenance check passes (35 claims, 0 failures), and every run used
   `PYTHONDONTWRITEBYTECODE=1`.

## Test evidence

`PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s WO-05 -p
'test_wo05.py' -v` runs 16 tests, all passing. The new tests cover:
- Data_S2 equals the detected-only mean, and the zero-inclusive rule fails.
- The cell census counts 3,015 covered zeros and 19 no-coverage cells.
- The estimators order as conditional > zero-inclusive > pooled.
- ×24.817 is reproduced.
- The J ratios hold exactly for each correction, and each correction moves
  the headroom in the stated direction.
- The stationary estimate never enters the headroom inputs.
- The quantity label forbids "baseline" and "raw" readings.

`run_wo05.py` exits 0 with G5.1–G5.5 true. The WO-01 to WO-04 suites still
pass.

VERDICT: PASS
