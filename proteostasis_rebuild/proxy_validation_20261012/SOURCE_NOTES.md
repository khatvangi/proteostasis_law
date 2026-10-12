# Source notes and citation anchors
Retrieved/inspected 2026-10-12. Web captures are tool-rendered text, not downloaded origin bytes. Their hashes authenticate the capture only.
Direct HTTP failures and available-byte digests are recorded separately in HTTP_RETRIEVAL.tsv.

## Pranjic (published journal article)
Pranjic M, Spät P, Semanjski Curkovic M, Macek B, Gruic-Sovulj I, Mocibob M. Resilience and proteome response of Escherichia coli to high levels of isoleucine mistranslation. Int J Biol Macromol. 2024;262(Pt 1):130068. DOI 10.1016/j.ijbiomac.2024.130068. PMID 38340920; PXD041341 / JPST002117.
https://pubmed.ncbi.nlm.nih.gov/38340920/ ; https://www.sciencedirect.com/science/article/abs/pii/S0141813024008717
Publisher/PubMed indexed abstract establishes editing-defective IleRS, Val/Nva exposure, up to 20% Ile substitution by high-resolution MS, aggregation and proteostasis response. Not a measured arrival rate or machine occupancy. Captured in WEB_EVIDENCE_CAPTURE.txt.

## Martens readable citation (preprint)
Martens, Andrew T., and Vincent J. Hilser. Chaperone saturation mediates translation and protein folding efficiency. bioRxiv [Preprint], 2025 Jun 28;2025.06.25.661590, version 1. DOI 10.1101/2025.06.25.661590. PMID 40667059. PMCID PMC12262517.
https://pmc.ncbi.nlm.nih.gov/articles/PMC12262517/ ; https://pubmed.ncbi.nlm.nih.gov/40667059/
PMC indexed record explicitly labels it a preprint that has not been peer reviewed. Synonymous Luc changes preserve amino-acid sequence; reported approximate extrema: up to 70-fold Luc protein and at most about twofold chaperone response. These are author summaries, not matched new numerical cells. WEB_CITATION_CAPTURE.txt holds the retrieved metadata.
Figure 2 caption on PubMed: Luc/beta-galactosidase protein measured by immunoblot in arbitrary units; activity n=3 with 95% CIs; regulon panels use relative ribosome-protected expression and adjusted p<0.10; activity–sigma32 fit uses effective-variance standard errors. No exact plotted coordinate extraction performed.
GEO GSE104303 lists GSE104303_mRNA_mapped_read_counts.tsv.gz (47.4 KB) and GSE104303_mapped_read_counts.tsv.gz (177.6 KB); PMC lists NIHPP2025.06.25.661590V1-supplement-1.pdf (557.6 KB). Read-count filenames alone do not authenticate protein measurements. WEB_EVIDENCE_CAPTURE.txt preserves the GEO listing.

## Local originals and line coordinates
Reingewertz original full text is one flattened line, not a missing source: experiment_grounded_rebuild_20261010/novelty_audit/sources/D_htpg_codon_series/reingewertz2024.txt:1. reingewertz_readable.txt changes whitespace only to make citations usable. Table 1: lines 119–121; codons: 133–162; selection: 45–79.
Fauvet primary full text: experiment_grounded_rebuild_20261010/kinetics/fauvet2021_PMC8082187_fulltext.txt:73. DOI 10.3389/fmolb.2021.653073. Total-contrast covariance unresolved; do not compute p from per-protein SDs.
FoldEco DOI 10.1016/j.celrep.2012.02.011; FoldEcoSlim DOI 10.1002/pro.3639. FoldEcoSlim source lines 38–46 explicitly treat growth, free sigma32 and feedback; these are prior modeling concepts, not empirical tests of this dossier's inequalities.
