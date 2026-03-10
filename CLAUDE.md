# Proteostasis Law — Triplet Code Necessity

## status: ACTIVE — CRITICAL AUDIT NEEDED

proves triplet codons are mathematically necessary under proteostasis constraint.
K=47 mode-capacity bound (demand) exceeds doublet ceiling of 16.

## core result (SOLID)
- E. coli: K=47, B. subtilis: K=48, S. cerevisiae: K=48
- all 9 two-codon AAs have k_A=2, contributing 18 > 16 (doublet impossible)
- mode assignments from per-AA median split in (μ, tAI) space

## CRITICAL AUDIT (must fix before submission)
1. **fabricated μ values** — section 3.2 table has μ=0.135, 0.120 etc. that exist nowhere in the pipeline. actual Landerer values are 1000x-3000x smaller.
2. **wrong enriched codons** — Cys: pipeline says TGC, manuscript says TGT. Glu: pipeline says GAG, manuscript says GAA.
3. **PDB deduplication bug** — ligand counts in metal_codon_bias_summary.csv inflated ~9x by counting each PDB separately. must deduplicate by (gene, uniprot_pos).
4. **no manuscript file** — prose exists only in conversation context.

## fix checklist
- [ ] deduplicate ligand counts by (gene, uniprot_pos) in summarize_codon_bias_metals.py
- [ ] rerun Fisher's exact with deduplicated counts
- [ ] replace fabricated μ values with actual Landerer values
- [ ] fix enriched codons: Cys→TGC, Glu→GAG
- [ ] rewrite Pareto interpretation with corrected data
- [ ] regenerate Figure 3
- [ ] save manuscript as .md or .tex in the repo

## key data files
- `errors/codon_modes_ecoli.tsv` — 61 sense codons with μ, tAI, modes
- `errors/aa_mode_summary.tsv` — k_A per amino acid
- `errors/codon_error_rates.tsv` — Landerer μ values
- `metals/metal_codon_bias_summary.csv` — enrichment counts (NEEDS DEDUP)
- `cross_species/results/cross_species_comparison.tsv` — K across organisms

## related projects
- `triplet-proof/` — SGC optimality (Monte Carlo), COMPLETE
- `proteostasis-paper/` — manuscript repo (empty manuscript/ dir)
- `codon_project/` — circular permutation analysis
