# STATUS — proteostasis rebuild (2026-09-29)

| WO | Verdict | One line |
|---|---|---|
| WO-00 | PASS | claim map, legacy hashes, legacy numbers reproduced |
| WO-01 | PASS | units; 1/T_gen overstates synthesis by 1/ln2; no dilution |
| WO-02 | PASS | conservative model; legacy Phi inflow has no donor |
| WO-03 | PASS | finite-pool / driven-cycle / competition chaperone models |
| WO-04 | PASS | legacy "collapse" is the imposed A_max gate, not a fold |
| WO-05 | **PASS** (reopened; was BLOCKED) | semantics fixed; x25 recomputed |
| WO-06 | PENDING, next eligible | literature parameter audit, not started |
| WO-07..10 | PENDING | |

## WO-05 in brief

- Landerer Data_S2 "mean" is a mean over datasets with **≥1 detected
  substitution**. This is verified from the per-dataset Data_S4 counts
  (3,015 covered zero cells excluded, 19 cells with no coverage). n = (sd/se)^2
  counts detection-positive datasets, not covering ones.
- Legacy ×24.817 → ×17.32 with the `(1-S)` double discount removed, ×35.88
  with `ln2/T_gen` alone, and **×25.06 with both**. The two corrections cancel
  to within 1%.
- With the selection-bias-free input (zero-inclusive mean, 2.53e-4/codon) the
  headroom is **×63.1**. With the PSM-pooled input it is ×112.9. Across the 80
  individual datasets it ranges from ×4.6 to ×2370.
- These are legacy-model-conditional eTEL aggregate statistics
  (MS-detected, substitution level, lower bound). **They are not an E. coli
  physiological headroom.** "×25 / one order inside" is not supported as a
  property of E. coli.
- Stikeleather 2026 (1.82e-3/codon) is kept as a separate **stationary-phase**
  estimate, and no headroom is computed from it.

## Open items carried forward

- WO-06: every parameter's source, organism and condition, with bibliographic
  checks.
- WO-07: an exponential vs stationary bundle is needed before any error rate
  (including Stikeleather's) is combined with growth-phase parameters.
- The legacy-model theta sweep (C23) and supraadditivity (C24) values have
  not been recomputed under the corrected mapping.
- The Stikeleather SE exponent is unverified (lost in the PMC text).
