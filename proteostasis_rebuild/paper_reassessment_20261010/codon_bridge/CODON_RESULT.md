# Codon result — 2026-10-10

Every number below is recomputed by `reanalyse.py` and `detection_controls.py` in this folder, from the
Landerer 2024 S4 per-dataset files and the Stikeleather 2026 `218_stats.csv` (SHA256 values in
`source_hashes.json`). Tests: `python -m unittest test_reanalyse -v`, 6/6 pass (run 2026-10-10).

## 1. What the data are

- **Landerer (exponential, mass spectrometry):** 80 E. coli datasets, 4,880 possible sense-codon cells.
  4,861 are covered; 3,015 of these have zero detected substitutions; 1,846 have at least one. In total,
  58,990 error PSMs over 4.94 × 10⁸ codon observations. The published Data_S2 per-codon "mean" is
  reproduced (max error 2.6 × 10⁻⁹) **only** by dropping the 3,015 covered zero cells. It is a mean
  conditional on detection, not a rate.
- Usage-weighted E. coli rate: **6.33 × 10⁻⁴ conditional** (the manuscript's value), **2.53 × 10⁻⁴
  zero-inclusive**, **1.41 × 10⁻⁴ pooled** (errors/observations). The manuscript's evaluation point is
  therefore inflated 2.5–4.5× by a zero-dropping convention in the source table.
- **Stikeleather (stationary phase, 3 biological replicates):** 6,692 errors over 3.68 × 10⁶ codons;
  replicate rates 1.71, 1.92, 1.83 × 10⁻³ (mean 1.82 × 10⁻³, SE 0.06 × 10⁻³). The per-replicate counts are
  exact integers (residual 7 × 10⁻¹²). This is a different growth state and a different detection pipeline.
  It is not an exponential-phase baseline and must not be pooled with Landerer. Its expression-level table
  has 11 rows where a zero mutation count carries a non-zero rate, so that table was not used.

## 2. Same-amino-acid codon contrasts (87 synonymous pairs)

The per-dataset rate difference was tested with an exact sign test, Holm-corrected over 174 tests (all
data, plus a depth-balanced subset with both depths ≥ 10,000 and within 2× of each other).

- Landerer: 17 all-data hits, 2 balanced hits. Stikeleather: 0 and 0.
- **Calibration (the important part).** Under a null with equal rates for every synonym in a
  dataset/amino acid, the observed unequal depths, and Poisson sampling (500 simulations,
  seed 20261010):
  - The **all-data sign test is badly miscalibrated.** At least one false Holm hit occurred in 100% of
    simulations (mean 6.5 false hits per simulation). The false hits concentrate in rare-codon pairs:
    ATA–ATC 84%, CTA–CTG 78%, ATA–ATT 64%, CCC–CCG 18%, CGA–CGC 12%. A shallow codon shows zero more often,
    so the sign of the difference is biased even when the rates are equal. **11 of the 17 all-data hits
    are not interpretable.**
  - The **depth-balanced test is calibrated:** family-wise false-hit rate 1%.
- **What survives.** Serine AGC has a higher detected rate than the other serine codons. Both balanced hits
  (AGC–TCC, 65 datasets, Holm p = 1 × 10⁻⁵; AGC–TCT, 67 datasets, Holm p = 1.1 × 10⁻⁴) have a 0% null hit
  rate, and so do the all-data AGC contrasts. The pooled rates are AGC 3.1 × 10⁻⁴ against TCC 2.6 × 10⁻⁵
  and TCT 1.8 × 10⁻⁵, a 12–17× difference. **96% of AGC's 2,208 detected errors are Ser→Asn**, which is
  consistent with misreading of AGC as AAC (a single middle-position near-cognate channel). The detected
  spectra of TCC/TCT are mostly Ala/Pro instead (total-variation distance 0.94–0.96).
- **Stikeleather does not reproduce the AGC excess.** There AGC is 1.2 × 10⁻⁴, lower than TCC
  (3.4 × 10⁻⁴), but its error counts are tiny (6 for AGC). With 3 replicates a sign test cannot go below
  p = 0.25, so this is **no test either way**, not a refutation.
- **Proline CCA vs CCG** (raised in the earlier draft): Landerer has 13 positive and 32 negative differences.
  The bootstrap CI on the mean (−2.1 × 10⁻⁴, +1.3 × 10⁻⁴) spans zero, the sign is not stable under
  leave-one-dataset-out, Holm p = 0.96, and there is no balanced pair. Stikeleather is 3/3 negative,
  p = 0.25. **No reliable CCA–CCG contrast exists.**

## 3. Remaining detection caveats (not removable with these files)

PSM quality: 15,705 of the error PSMs have PeptideProphet < 0.99 (minimum 0.05). PeptideProphet is not
a decoy-controlled substitution-localization FDR. A residual false-discovery contribution (isobaric or
near-isobaric modifications, mis-localization) is therefore **not excluded** for any codon, including
AGC→Asn (+27.01 Da). Undetected peptides cannot be recovered from the retained rows.

## 4. Claim triage for the manuscript

## 4. Claim triage for the manuscript

**Withdraw as stated: usage-weighted mu = 6.3e-4 as "the" E. coli per-codon error rate, and the x25
margin computed at it.** 6.3e-4 is a mean conditional on detection. The zero-inclusive mean is 2.5e-4
and the pooled rate 1.4e-4. The earlier WO-05 recalculation puts the margin at about x63 at the
zero-inclusive input and x113 at the pooled input, and those are still eTEL-aggregate, legacy-model
numbers, not an E. coli baseline. The paper cannot keep "one order of magnitude inside the envelope"
while using a rate defined by dropping zeros.

**Withdraw: any specific synonymous contrast other than serine AGC.** The all-data sign test used to
generate those contrasts has a 100% family-wise false-positive rate under equal rates with the observed
depths. The proline CCA-CCG contrast in particular has no support.

**Survives, narrowed: detected mistranslation is organized at the amino-acid level.** The one
depth-calibrated, exposure-balanced, per-dataset result in these data is that serine AGC is decoded less
accurately than TCC/TCT by about 12-17x, through a single Ser->Asn channel. That is a within-amino-acid,
synonymous-codon difference in measured error rate, which is the weaker form of the manuscript's claim.
It rests on one amino acid in one growth state and one detection pipeline.

**Not established, and the paper should say so: that a same-synonym substitution error changes the
folding outcome of the native protein.** I found no evidence for it, in these data or in the literature
checked here. The Landerer and Stikeleather files measure amino-acid misincorporation; they contain no
folding, aggregation or degradation readout at all. A synonymous codon change alters *which* wrong amino
acid is likely and *how often*, not the identity of the correctly made protein. The closest published
result runs the other way and is worth stating in the paper because it is the strongest independent
support for a codon-to-folding link: Reingewertz et al. 2024 (Protein Sci 33:e5087,
https://doi.org/10.1002/pro.5087; retrieved from PubMed, PMID 39074255) selected mutations in mouse DHFR
that reduce its dependence on GroEL for folding, and found that **synonymous** substitutions did so,
apparently by changing local translation rate. That is a codon-dependent folding/chaperone demand in an
extant system, which is exactly the mechanism the envelope argument needs - but it is co-translational
kinetics, not a mistranslation effect, and it is a heterologous client in E. coli, not a native protein.

**Stationary phase must be reported separately.** Stikeleather's 1.8e-3 is 3-13x higher than any
Landerer exponential estimate. It is a different physiological state. Pooling the two, or quoting the
higher value as "the" measured error rate, is not defensible.

## 5. What would settle the open points

1. A decoy-controlled, localization-scored reanalysis of the Landerer PSMs, which would say whether the
   AGC->Asn excess is decoding or a near-isobaric artifact. Needs the raw spectra, not these tables.
2. More than three replicates in stationary phase, so the sign test there can go below p = 0.25.
3. Any experiment pairing a synonymous substitution with a folding or chaperone-demand readout in a
   native E. coli protein. Reingewertz did this for heterologous mDHFR; nothing equivalent exists here.
