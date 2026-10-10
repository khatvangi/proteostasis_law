# Where the proteostasis paper stands after the 10 October 2026 reassessment

## The short version

Two of the manuscript's quantitative anchors do not survive: the error rate it evaluates at is
conditional on detection, and its synonymous-codon contrasts come from a test that fires on noise. One
result survives, narrower. The route model's central prediction was never tested by its own ensemble,
and one of its structural assumptions is contradicted by a 2004 experiment.

## What the data now say

The Landerer per-dataset files (80 E. coli datasets, 4,861 covered codon cells) reproduce the published
per-codon table only if the 3,015 cells that were covered but had zero detected substitutions are
discarded. Keeping them gives a usage-weighted rate of 2.5 × 10⁻⁴; pooling errors over observations gives
1.4 × 10⁻⁴. The manuscript's 6.3 × 10⁻⁴ is a mean over detections only, 2.5 to 4.5 times too high. At the
honest inputs the margin is roughly ×63 to ×113, not ×25, so "one order of magnitude inside the
envelope" has to go.

The synonymous comparisons are worse. Running the same per-dataset sign test on simulated data where
every synonym has an identical rate, at the real depths, gave at least one Holm-significant false hit in
100% of 500 simulations, averaging 6.5 each. They sit exactly where the paper's hits sit: rare codons
with lopsided depth (ATA–ATC fires 84% of the time under the null). Eleven of seventeen hits are
uninterpretable. The proline CCA–CCG contrast has no support — its bootstrap interval spans zero and its
sign flips when datasets are dropped.

One amino acid survives, and survives properly. Restricting to pairs seen at least 10,000 times each and
within twofold of each other gives a calibrated test (1% false-hit rate), and serine AGC is decoded 12–17
fold less accurately than TCC and TCT, Holm p of 1 × 10⁻⁵ and 1.1 × 10⁻⁴, with a zero false-hit rate under
the null. Ninety-six percent of AGC's 2,208 errors are one substitution, Ser→Asn, consistent with a
single near-cognate channel. That is a real within-amino-acid difference in decoding accuracy — in one
amino acid, one growth phase, one pipeline, with PSM scores that are not a decoy-controlled localization
FDR, so an isobaric artifact is not excluded. Stationary-phase data (1.8 × 10⁻³, three replicates) neither confirms nor
refutes it, and must be reported separately rather than pooled.

## What must be withdrawn

The 6.3 × 10⁻⁴ evaluation point and the ×25 margin computed at it. Every synonymous contrast except
serine AGC. And any suggestion that a same-synonym substitution error changes the folding outcome of the
native protein: there is no evidence for that here or in the literature I checked, and these datasets
contain no folding readout at all.

## The principle that is left, stated sharply

Client influx must be routed through a capacitated network, and the network is a hypergraph whose edges
are *qualified routes*, not machines. Feasibility is a Hall-type condition: every subset of client
classes must have demand no larger than the capacity of the routes that can serve it. Lethality of a set
of machines means that set is a cut in the qualified-route hypergraph. A single-pool capacity sum provably cannot say which machine is binding; the cut structure can.

Two experimental facts fix the content. Route qualification is *conditional*, not fixed: Δtig ΔdnaK
E. coli is dead at 30 °C, viable at 23 °C, and rescued at 30 °C by GroEL/GroES overproduction
(Vorderwülbecke 2004, FEBS Lett 559:181–187; Genevaux 2004, EMBO Rep 5:195–200). And demand per client
is codon-dependent: synonymous substitutions in mouse DHFR reduce its dependence on GroEL in E. coli,
apparently by changing local translation rate (Reingewertz 2024, Protein Sci 33:e5087). Conditional
cooperation plus codon-set demand weights is the whole claim.

The shipped route model fails that test twice. With TF and DnaK capacity at zero it delivers zero
essential flux in 200 of 200 draws even at 100× GroEL — the rescue is structurally inexpressible. And because demand is normalized to a margin insisting every aggregate be productively recovered through
the tiny ClpB pool, the folding machines end up in large excess: matched GroEL and DnaK losses change the
folded fraction by exactly zero in all 200 draws, and all four pair-interaction terms are zero. The
supra-additivity prediction was not refuted here; it was never tested. Dropping the aggregate quota
raises the demand scale by a median 5.7-fold.

## A prediction that would discriminate, and its blocker

Serial and parallel route structure make opposite sign predictions for growth interactions. Two machines
on the same serial path share a cut, so removing both should be no worse than multiplicative; two
machines offering parallel qualified routes to the same clients should be aggravating. On Babu's 2014
E. coli screen this is a sign test with no free parameters: tig × dnaK should be aggravating, while
dnaK × a groL hypomorph should not. The transcribed Kanemori protease table is consistent with the
parallel half (combined ΔhslVU Δ(clpPX-lon) plates at 3.3 × 10⁻⁵ against 1.3 and 0.22 alone), but that is
retrospective arithmetic on rounded published values. **The real test is blocked:** the Babu
supplementary XLS could not be retrieved — all 14 download attempts failed DNS resolution — and it lists
significant hits only, so absent pairs cannot be scored as neutral.

## Status, exactly

Recomputed and checked this run: codon reanalysis (6/6 tests pass), route model (15/15 tests pass, 200
trials, maximum duality gap 1 × 10⁻¹⁶, byte-identical rerun). Not done: original Kerner, Calloni and
Babu supplements, any decoy-controlled reanalysis of the mass spectra, and any client-resolved demand
flux. The canonical manuscript, legacy directories and WO files were not modified.
