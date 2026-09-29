# WO-05 — error-to-burden semantics

## Verdict

**VERDICT: BLOCKED** — G5.1, G5.2, and G5.4 pass; G5.3's required
double-discount rejection tests pass. G5.5 cannot be evaluated from the
available empirical input semantics, so no ×25 headroom was recomputed and
the sequential loop stops here.

## G5.1 — three quantities and units

All rates are dimensionless probabilities per codon translated; the resulting
burden flux is in fraction of the modeled protein pool per second.

* **Raw decoding error** `e_raw`: probability that a decoding event inserts
  a non-cognate amino acid or tRNA outcome, including synonymous outcomes.
  Under the legacy uniform synonymous probability `S`, the model conversion
  to amino-acid substitution probability is `f_sub=e_raw*(1-S)`.
* **Amino-acid substitution rate** `f_sub`: probability per codon that the
  encoded amino acid changes, summed over destination amino acids. It is
  already post-synonymous filtering.
* **MS-detected rate** `r_MS`: Landerer et al.'s per-codon rate of detected
  PSMs carrying a substitution at the covered position divided by PSMs
  covering that codon position. This is amino-acid-level and destination
  summed, not raw decoding error; synonymous events are invisible, and
  I/L substitutions are mass-indistinguishable. The measured value is thus
  an incomplete observation of `f_sub`, not a raw-error rate.

The legacy codon table's `mu` values equal Data_S2's `mean` exactly. In
Data_S2, `(sd/se)^2` implies 3–68 contributing datasets per codon (maximum
deviation from an integer `5.13e-5`), below the 80 total datasets. That is
consistent with a mean across datasets that reported at least one detected
substitution for that codon; it is not a mean over all codon-covering PSMs
pooled across all datasets. The positive-detection selection is central to
the G5.5 block below.

## G5.2 — legacy mapping inventory

The automated inventory in `run_wo05.py:MAPPINGS` checked that every cited
anchor occurs at its stated line; all 21 anchors were found. Classification
is by the quantity actually supplied to the mapping:

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

`error_semantics.py` gives rates explicit kinds (`RAW`, `SUB`, `MS`). A raw
rate is converted with `(1-S)` once; substitution input passes through; MS
input is treated as a substitution-level observation (or an explicitly
detectability-adjusted value) without another synonymous discount. The
guarded legacy conversion raises `TypeError` for both `SUB` and `MS` inputs.
Tests confirm a raw input agrees with the legacy formula, while passing
substitution-level inputs through the legacy discount is rejected.

## G5.4 — arithmetic threshold

At `N=300`, `P_correct=0.70`, `S=0.30`, and `p_misfold=0.30`, the exact raw
threshold is
`(1-0.70**(1/300))/((1-0.30)*0.30) = 0.005658142850513878`.
The exact substitution-level threshold is `0.003960699995359715`. The
legacy-quoted `1.19e-3` is `-ln(0.70)/300 = 0.0011889164797957749`; it comes
from the separate Part A calculation that sets `(1-S)*p_m=1` (legacy
`arithmetic_stress_test.py:202`), omitting both stated factors. It is not the
threshold at the stated parameter values.

## G5.5 — blocked empirical headroom

The ×25 headline uses a usage-weighted average of Data_S2 codon means. Those
means are reported only for datasets with detected substitutions at that
codon, while the number of codon-covering datasets with zero detected
substitutions is not present in the supplied table. The `sd/se` inference
identifies fewer than 80 contributing datasets but does not identify the
missing zero cells or their coverage denominators. Consequently, weighting
the reported means cannot recover an unconditional empirical substitution
rate. Using these conditional means directly would select on positive
detections; multiplying by `(1-S)` would then add the separate double
discount. Neither yields a defensible corrected ×25 headroom.

`run_wo05.py:g55` therefore returns `BLOCKED` without computing headroom;
`test_g55_blocks_unresolved_empirical_denominator` asserts this behavior.
The runner exits with status 2 for the blocked gate. The old ×24.817 value is
not reproduced as a corrected value and is not replaced by an assumption.

## Adversarial self-review

The classification of MS as a substitution-level observation is supported by
the PSM numerator and destination-summed Data_S2 rates; MS invisibility means
it is not an unbiased estimate of all substitutions. The inference that the
means exclude zero-detection datasets is based on reported `sd/se` and is
not a substitute for the missing per-dataset coverage table. A zero-inclusive
mean could be computed only by assuming a common denominator of 80, but the
available data do not establish that every dataset covers every codon. That
assumption would invent the empirical input needed for the headline, so G5.5
remains blocked. G5.1, G5.2, and G5.4 are independently testable and passed;
they do not cure that missing denominator.

## Test evidence

`PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s WO-05 -p
'test_wo05.py' -v` — 9 tests passed. The runner reported G5.1, G5.2, and
G5.4 true, G5.5 `BLOCKED`, and exited 2 as designed.
