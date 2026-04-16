# Context-Conditional Codon Deployment Under Proteostasis Constraint

---

## Abstract

Synonymous codons differ systematically in mistranslation rate and translational supply, and cells deploy them non-uniformly across sites. Whether this operational diversity is a structural requirement of the genetic code, or an incidental consequence of degeneracy, has not been addressed. We derive the synonymy necessity constraint: if translation viability under finite quality-control capacity requires context-conditional codon deployment, then the code must provide at least two codons per amino acid across the active alphabet. For a doublet code with one stop signal, the mean synonymy per amino acid is (15 - A)/A, which falls to 0.5 at A = 10 and to zero at A = 15. Doublet codes therefore cannot support within-amino-acid deployment for any alphabet approaching the canonical 20. This result is parameter-free. It is logically independent of, but converges with, the decoder-capacity constraint derived from block-partition analysis in the companion paper, which finds A_eff = 9 for doublet codes under G:U wobble.

The empirical test is whether extant triplet-code organisms actually deploy synonymy in structured ways. At metal-binding sites in E. coli (N = 17,166 residues), Asp and His favor the higher-mu synonym (OR = 1.27-1.30), accepting greater misincorporation risk for higher translational supply, while Cys favors the lower-mu synonym (OR = 1.29). Opposite preferences across amino acids at the same class of sites rule out single-axis optimization. Across E. coli, B. subtilis, and S. cerevisiae, slow leucine synonyms comprise 15-39% of Leu positions and the rank order of operational spread per amino acid is conserved (Spearman rho = 0.54-0.77) despite ~2 billion years of divergence.

Triplets are the shortest architecture compatible with both naming 20 amino acids and providing the within-amino-acid operational diversity that cells demonstrably use. We do not claim to reconstruct the historical path by which triplets arose; we claim that doublet architectures are excluded by two independent formal arguments, and that the synonymy triplets supply is a functionally deployed operational resource rather than idle redundancy.

---

## Introduction

The triplet genetic code is one of biology's deepest universals. All free-living organisms translate messenger RNA in three-nucleotide steps, yielding 64 codons that map to 20 canonical amino acids plus termination signals. The standard explanation for triplet length is arithmetic: a doublet code provides only 4^2 = 16 symbols, insufficient to uniquely specify 20+ distinct outputs, whereas triplets supply 4^3 = 64 with room to spare (1-3). That reasoning, while correct, answers only a naming question. It does not explain why the code exhibits structured degeneracy, why synonymous codons are distributed unevenly across amino acids, why the third ("wobble") position absorbs most redundancy, or why codon choice at individual sites correlates with expression level, protein structure, and cellular stress (4-6).

Classic frameworks address different facets of code structure. Stereochemical theories posit direct chemical affinity between codons and their cognate amino acids (3, 10). Coevolutionary accounts propose that the code expanded in concert with biosynthetic pathways (4, 5). Error-minimization analyses demonstrate that the standard code is unusually robust to mistranslation damage (6, 7). Information-theoretic treatments model the code as a noisy channel and derive bounds on alphabet size (8, 9). Frozen-accident views emphasize that once established, global reassignment is evolutionarily prohibitive (2, 16).

These frameworks address assignments, robustness, or reassignment dynamics: the mapping from codons to amino acids, and whether that mapping is good. They do not address a prior question: what amount of within-amino-acid operational multiplicity must a code possess before local context can be controlled at all? Error-minimization asks whether a given assignment cushions mistakes. We ask whether the code provides enough distinct synonyms per amino acid for cells to deploy different ones at different sites according to local demands.

Recent work on translation fidelity identifies a candidate use. Drummond and Wilke established that mistranslation-induced protein misfolding is a dominant selective pressure on coding sequences: highly expressed genes avoid codons prone to costly misincorporation, and codon usage bias correlates with predicted aggregation propensity (11). Balch and colleagues formalized proteostasis as a network-level property (a balance between synthesis, folding, quality control, and degradation) that can collapse under excess burden (12). If translation noise must remain below a finite quality-control capacity, and synonymous codons offer distinct noise profiles, then the cell's ability to choose among synonyms becomes functionally consequential.

Here we develop this idea in two parts. First, we formalize the proteostasis viability condition (J_in < J_crit) and derive the synonymy necessity constraint: if cells must deploy different synonyms at different sites to keep misfolding inflow below the collapse threshold, doublet codes lack sufficient within-amino-acid diversity to support this deployment, for any amino acid alphabet of 10 or more residues. Second, we test empirically whether extant triplet-code organisms deploy synonymy in context-dependent patterns. The formal constraint is parameter-free, depending only on alphabet size and code length. The empirical evidence is necessarily organism-specific: mu and tAI are measured in E. coli. Together they show that doublet codes are excluded by arithmetic and that the synonymy triplets provide is functionally used.

---

## Theory: A viability constraint on translation-induced misfolding

### The collapse threshold

A viable cell must keep the burden of misfold-prone species below the capacity of its quality-control machinery. We formalize this as the existence of a stable steady state for the misfold pool. Let P(t) denote the concentration of misfold-prone material. A minimal dynamical model takes the form:

dP/dt = J_in - R(P) - A(P)

where J_in is the influx of problematic translation products, R(P) is net rescue/removal by refolding, degradation, and sequestration, and A(P) is irreversible conversion into aggregation-like states. Rescue pathways are capacity-limited and saturable; aggregation introduces nonlinearity and positive feedback as load rises.

With saturating R(P) and superlinear A(P), the function F(P) = R(P) + A(P) generically develops a shape supporting multiple fixed points. As J_in increases, the stable low-P fixed point and an intermediate unstable fixed point collide and annihilate at a critical inflow J_crit. For J_in > J_crit, no low-P stable solution exists and trajectories undergo runaway collapse. The saddle-node bifurcation is generic to any system with saturating rescue and superlinear aggregation; we invoke it to justify a threshold structure for the misfolding burden, not to parameterize proteostasis dynamics in any organism. Scope and falsifiability are discussed below.

### Connecting J_in to code architecture

To connect J_in to genetics, we write the inflow as a proteome-weighted sum of site hazards. For proteins p synthesized at flux E_p, with sites i weighted by criticality kappa_{p,i}:

J_in ~ sum_p E_p sum_{i in p} kappa_{p,i} * q_{p,i}(c_{p,i})

where c_{p,i} is the codon used at site i and q_{p,i}(c) is the probability that decoding produces a proteostasis-relevant failure:

q_{p,i}(c) = mu(c) * [1 - S(c)] * p_let(p,i)

with mu(c) representing codon-specific misincorporation risk (14, 15), S(c) representing codon-level synonymy shielding (the probability that a single-nucleotide perturbation at codon c produces a synonymous outcome rather than an amino acid change), and p_let(p,i) capturing how likely a nonsynonymous event at that site generates misfolding. S(c) is a property of each individual codon's neighborhood in the code table; it differs from the aggregate third-position shielding fraction S_3, which averages over all sense codons at position 3 (S_3 = 0.69 for the standard code; Methods). This accounting identity links code architecture to the dynamical constraint. Triplet codes enable a dedicated wobble position where S(c) is systematically high, lowering effective hazard relative to doublet alternatives (Fig. 1C) (1).

### Operational axes and context-conditional selection

We characterize each codon by two measured quantities: mistranslation hazard mu(c) from mass-spectrometry misincorporation data in E. coli (15), and translational supply nu(c) proxied by the tRNA Adaptation Index tAI, which reflects tRNA gene copy number adaptation rather than directly measured ribosome transit times (13). Both quantities are defined relative to the extant E. coli translation apparatus, a triplet-code organism. The mu and nu values therefore describe the operational landscape of the existing triplet code, not properties that would hold for a hypothetical doublet or quadruplet system. We use them to characterize how triplet-code organisms deploy their available synonyms, not to derive the necessity of triplets from first principles.

The simplest model of codon selection assumes cells should prefer codons with lower mu and higher nu. Under such monotone optimization, any codon with both higher mu and lower nu than a synonym would be strictly inferior and never used. The strongest evidence against monotone selection comes from two-codon amino acids at high-criticality sites (Result II): Asp and His at metal-binding positions favor the higher-mu synonym, accepting greater misincorporation risk for higher translational supply, while Cys favors the lower-mu option. No single-axis rule explains both patterns. The genome-wide pattern confirms this: leucine, with six synonymous codons, deploys slow synonyms (CUA, UUA) at ~16% of positions in E. coli despite their scoring poorly on both mu and nu axes relative to CUG (Table S1).

The resolution is that translational supply requirements are context-targeted: at sites requiring co-translational folding pauses, low-supply codons are advantageous (19, 20). We formalize this as a site-specific cost:

L_i(c | A) = lambda_mu(x_i) * mu(c) + lambda_nu(x_i) * (nu(c) - nu*(x_i))^2

The quadratic form captures a key feature: oversupply relative to the site target is costly (missed folding windows), and undersupply is costly (reduced throughput). The interior optimum varies by site.

### Operational diversity

For each amino acid A, define the operational spread Delta_A as the mean pairwise distance among its synonyms in standardized (mu, nu) space:

Delta_A = [1 / C(n_A, 2)] sum_{j<k} sqrt((mu_j_tilde - mu_k_tilde)^2 + (nu_j_tilde - nu_k_tilde)^2)

This metric quantifies how much operational range the synonymous codon set provides. It is computed from E. coli data and characterizes the existing triplet code's operational landscape.

### The synonymy necessity constraint

The proteostasis viability condition (J_in < J_crit) combined with heterogeneous site demands implies that cells require within-amino-acid operational choice: the ability to deploy different synonyms at different sites according to local accuracy-supply tradeoffs. This yields a formal constraint on code architecture:

**Synonymy necessity constraint.** If translation viability under finite quality-control capacity requires context-conditional codon deployment, then the code must provide at least two codons per amino acid for every amino acid in the active alphabet. For a 4-letter nucleotide alphabet at code length L with at least one stop signal, the mean synonymy per amino acid is (4^L - 1 - A)/A. At L = 2 and A >= 8, this quantity falls to <= 0.9. At L = 2 and A >= 10, it is <= 0.5, insufficient for any distributed deployment. At L = 3, mean synonymy is >= 2.0 for all A <= 21, covering the entire canonical amino acid alphabet.

This constraint is parameter-free: it depends on the alphabet size and code length, not on the values of mu, nu, or any proteostasis parameter. It is logically independent of the decoder-capacity constraint (A_eff = 9 at L = 2 under G:U wobble) derived in the companion paper (Author et al.) from block-partition analysis. The two results exclude doublets from different premises: the decoder-capacity constraint says doublets cannot *name* enough amino acids; the synonymy necessity constraint says doublets cannot provide operational diversity even for the amino acids they *can* name. Their convergence on triplets as the minimal sufficient code length is independent confirmation from orthogonal arguments.

The synonymy necessity constraint and the proteostasis framework address different questions. The constraint says doublets lack synonymy. The framework says synonymy is functionally consequential. The empirical results (below) test the second claim within the existing triplet code.

---

## Results

### Result I: The synonymy necessity constraint excludes doublet architectures

The synonymy necessity constraint (Theory) is structural and parameter-free (Fig. 2A):

| Architecture | Sense codons | Synonymy for 15 AAs | Synonymy for 10 AAs |
|---|---|---|---|
| Doublet | <=15 | 0 | 0.5 per AA |
| Triplet | 61 | 3.1 per AA | 5.1 per AA |
| Quadruplet | ~255 | ~16 per AA (excess) | ~25 per AA (excess) |

The companion paper (Author et al.) derives the independent decoder-capacity constraint from block-partition analysis, showing A_eff = 9 for doublets under G:U wobble and proving beta invariance between triplet and quadruplet codes at matched adaptor count. We do not reproduce that analysis here. For the present paper, the relevant point is different: even if doublet codes could name 9-10 amino acids (as the decoder-capacity analysis shows they can), the resulting code would have zero or negligible synonymy per amino acid. Without synonymy, no within-amino-acid deployment strategy is possible, regardless of what that strategy might be. The synonymy necessity constraint and the decoder-capacity constraint exclude doublets from orthogonal premises.

---

### Result II: High-criticality sites exhibit context-dependent codon deployment

We tested whether synonymous codons are deployed non-uniformly at structurally critical sites using metal-coordinating residues as a proxy for high site criticality. Metal-binding positions are structurally brittle and often functionally essential. Mapping curated metal-binding residues onto E. coli coding sequences yielded N = 17,166 metal-binding residues across 1,847 proteins.

Across four canonical liganding amino acids, codon usage at metal-binding sites deviates from genome-wide background (Fig. 3A, 3B):

| Amino acid | Enriched codon | Odds ratio | 95% CI | p-value |
|---|---|---|---|---|
| Asp | GAC | 1.30 | [1.18, 1.43] | < 10^-3 |
| Cys | UGU | 1.29 | [1.15, 1.44] | < 10^-3 |
| Glu | GAA | 1.10 | [1.02, 1.18] | 0.01 |
| His | CAC | 1.27 | [1.14, 1.41] | < 10^-3 |

These odds ratios are statistically significant but modest in magnitude. For reference, codon usage bias in highly expressed E. coli genes produces odds ratios exceeding 3 for preferred codons (20). The metal-site signal is one force among several acting on codon choice at these positions; its diagnostic value lies in direction rather than strength.

The enriched codon is not uniformly the lower-mu synonym:

| Amino acid | Enriched codon | mu (enriched) | mu (alternative) | Consistent with mu-only? |
|---|---|---|---|---|
| Asp | GAC | 0.135 | GAU: 0.120 | No |
| Cys | UGU | 0.203 | UGC: 0.262 | Yes |
| Glu | GAA | 0.197 | GAG: 0.209 | Yes |
| His | CAC | 0.145 | CAU: 0.090 | No |

For Asp and His, the enriched codon has higher mu than its alternative. In (mu, nu) space, these codons have higher translational supply (tAI) at the cost of higher misincorporation risk (Fig. 3C). For Cys and Glu, enrichment is consistent with accuracy-first selection. Different amino acids at the same class of sites thus show opposite codon preferences. No single-axis optimization rule explains both patterns; each amino acid's enrichment reflects its own tradeoff geometry in (mu, nu) space.

---

### Result III: Non-monotone codon usage is conserved across species

We tested whether context-conditional selection is species-specific or universal using leucine, whose six codons span the full (mu, nu) range (Table S1, Fig. 2B).

| Species | Most common | Slow codon usage (tAI < 0.25) |
|---|---|---|
| *E. coli* K-12 | CUG (51%) | 39% |
| *B. subtilis* 168 | UUA (32%) | 15% |
| *S. cerevisiae* S288C | UUG (29%) | 31% |

In each species, slow Leu synonyms occur at non-trivial frequency (15-39% of positions). The pattern is qualitatively similar (slow codons always used) but quantitatively distinct: E. coli strongly favors CUG, S. cerevisiae favors UUG and UUA, and B. subtilis shows intermediate preference. Species-specific patterns reflect differing tRNA pools; the universal presence of slow synonyms supports non-monotone selection across the bacterial-eukaryotic divide.

Leu is the most permissive test case because its six codons span the full (mu, nu) range. To test whether the broader geometry of operational diversity is also conserved, we computed Delta_A for each amino acid in each species and asked whether amino acids with high diversity in one species also have high diversity in another. The rank correlations are positive across all three comparisons:

| Comparison | Spearman rho | p-value |
|---|---|---|
| *E. coli* vs *B. subtilis* | 0.68 | 0.001 |
| *E. coli* vs *S. cerevisiae* | 0.54 | 0.01 |
| *B. subtilis* vs *S. cerevisiae* | 0.77 | < 0.001 |

Amino acids with high Delta_A in one species tend to have high Delta_A in others despite ~2 billion years of divergence. This conservation is consistent with maintained selection on operational diversity rather than drift.

---

### Result IV: Observed operational diversity exceeds matched-null expectation

We asked whether the observed diversity structure could be a byproduct of degeneracy alone. We generated 10,000 matched null codes preserving degeneracy-class structure but randomizing codon-to-coordinate assignments (Fig. 4A):

z = 2.99, p = 0.0014

The observed diversity exceeds the null expectation. The z-score reflects moderate, not overwhelming, deviation. Two-codon amino acids universally achieve maximal separation in (mu, nu) space, while multi-codon amino acids show variable spread (Fig. 4B).

This null comparison provides supporting evidence but is not the primary logic. The core empirical finding is that codon deployment at high-criticality sites is non-random and non-monotone (Result II). The null test shows that within the triplet regime, the realized diversity structure is unlikely to be incidental.

---

## Discussion

### What the results show

Codon selection cannot be described by a single-axis optimum. At metal-binding sites, Asp and His favor the higher-mu synonym while Cys favors the lower-mu synonym. At six-codon families genome-wide, codons that score poorly on both measured axes are used at 15-39% of positions across species. Both patterns contradict the prediction that cells should always prefer the codon with lowest mu and highest nu. The resolution is context-conditional selection: each site incurs a cost with site-specific nu*(x), and the optimal codon varies with context. Sites requiring co-translational pauses favor low-supply codons; sites requiring rapid elongation favor high-supply codons; sites with extreme accuracy demands favor low-mu codons regardless of supply. Apparent contradictions in the codon usage literature (fast codons selected vs. slow codons selected) are both correct at different sites.

The functional consequence for code architecture is direct. Context-conditional deployment requires more than one synonym per amino acid, for every amino acid the cell must deploy this way. The synonymy necessity constraint formalizes this requirement and shows that doublet codes cannot meet it for any alphabet of 10 or more residues. The constraint is parameter-free: it depends only on alphabet size and code length, not on the specific values of mu, nu, or any proteostasis quantity. The empirical results in E. coli and the cross-species conservation of operational spread (Spearman rho = 0.54-0.77, ~2 Gya divergence) show that the premise of the constraint is satisfied in extant biology: cells do deploy synonymy context-dependently, and they do so in ways that track amino-acid identity rather than species identity. If this deployment is removed by collapsing the code to doublets, a dimension of cellular control present in existing organisms is eliminated.

The effect sizes at metal sites are modest (OR = 1.10-1.30), below the OR > 3 seen for preferred codons in highly expressed genes. This is expected: codon choice at any individual position reflects translational accuracy, speed, mRNA secondary structure, tRNA availability, mutational bias, and genetic drift, and the proteostasis-relevant component is one among several. The diagnostic value of the metal-site result is not the magnitude of any single enrichment but the directional heterogeneity across amino acids, which cannot be produced by any universal optimization rule.

### Relation to existing theories

Existing frameworks address code assignments (stereochemical and coevolutionary theories: 3-5, 10), robustness of assignments to error (error-minimization: 6, 7), or the difficulty of reassignment (frozen-accident: 2, 16). Each addresses a property of the mapping from codons to amino acids. None addresses what we ask: how much within-amino-acid operational multiplicity must the code possess for context-conditional deployment to be possible at all? Drummond and Wilke (11) showed that mistranslation costs shape codon usage within genes. The present results concern a prior question, what code architecture is required for such shaping to be possible, and find that the threshold separating doublets from triplets is where context-dependent selection becomes representable.

The companion paper (Author et al.) derives the decoder-capacity constraint from block-partition analysis: doublets have A_eff = 9 under G:U wobble, and triplets versus quadruplets tie on beta invariance at matched adaptor count. That result and the synonymy necessity constraint are logically independent. One derives from decoder physics (how many outputs can the translation apparatus distinguish), the other from operational demand (how many synonyms per amino acid must the code provide). They reach the same exclusion of doublets from orthogonal premises. The arithmetic excluding doublets does not exclude quadruplets, and the synonymy necessity constraint does not by itself prefer triplets over quadruplets either; quadruplets would provide more synonymy per amino acid. The companion paper's economy argument (25% shorter genomes, 25% lower frame damage) is what breaks the tie. The present paper shows that the synonymy triplets provide is not slack capacity, and asks an empirical question we cannot fully answer with extant data: whether quadruplet systems would utilize additional operational diversity or leave it unused.

### Scope, limitations, and what would falsify this

This paper addresses operational deployment of synonymous codons in extant translation systems. It does not reconstruct the evolutionary origin of the triplet code or the historical pathway by which the code was adopted. The pre-triplet world may have had different proteostasis demands, different amino acid alphabets, and different selective pressures; we make no claims about it.

The ODE model in Theory is not a quantitative proteostasis model. The saddle-node bifurcation at J_crit is generic to any system with saturating rescue and superlinear aggregation, and we do not estimate R(P) or A(P) from data or claim that any particular organism operates near the bifurcation. The ODE establishes the existence of a threshold structure; the falsifiable content of the paper is the empirical prediction that cells deploy synonymy context-dependently at structurally critical sites. If synonymous codons at high-criticality positions showed no departure from genome-wide frequencies, or if slow synonyms were absent from all species, the deployment hypothesis would be refuted. The synonymy necessity constraint would remain arithmetically valid but would lose its biological motivation: a code could formally satisfy it while cells made no functional use of the diversity.

All mu and nu values are measured in E. coli, limiting direct quantitative generalization to other organisms. Qualitative conclusions require only that synonyms differ along multiple axes, not precise values. The cost function L_i is phenomenological; we infer nu*(x) variation from usage patterns rather than measuring it directly. Metal-binding sites are a conservative proxy for high site criticality; extending to catalytic residues, buried cores, and domain boundaries would strengthen generality without changing the logical structure of the argument.

### Predictions

*Stress-dependent selection.* Lower effective proteostasis capacity (heat shock, chaperone depletion, aging) should intensify selective pressure on codon choice at high-kappa sites, because the margin between J_in and J_crit narrows.

*Vocabulary starvation.* Recoding experiments that eliminate synonymy (collapsing two-codon families to one) should produce fitness defects specifically under proteostasis stress, even when amino acid identity is preserved.

*Shielding-error alignment.* Position-specific decoding error rates should correlate with synonymy shielding: the noisiest position should be the one most enriched for synonymous absorption.

### Conclusion

Degeneracy is the diagnostic feature of the triplet code: 64 codons encoding 20 amino acids plus stops, with the ratio unevenly distributed across amino acids. Framed as slack capacity, this pattern is hard to explain. Framed as an operational resource that cells deploy context-dependently, the pattern is what the code must look like. Two-codon amino acids achieve maximal separation in (mu, nu) space. High-criticality sites use codon choice as a tunable accuracy-supply tradeoff. Amino acids retain the same operational ordering across ~2 billion years of divergence. The ratio that makes these patterns possible is excluded by doublet architectures, both from decoder physics (companion paper) and from the synonymy requirement derived here. The constraint that selects for triplets over doublets is not capacity alone. It is capacity plus operational multiplicity.

---

## Methods

### Overview and design logic

We implemented five linked analyses: (1) define codon operational coordinates and diversity metrics in E. coli; (2) state and evaluate the synonymy necessity constraint on doublet architectures; (3) test codon deployment at high-criticality metal-binding sites in E. coli; (4) evaluate cross-species conservation of slow-codon usage and diversity geometry; (5) test for non-accidental structure using matched-null ensembles.

---

### Codon operational coordinates

Genomic codons were converted T->U for presentation throughout figures and main text.

*Mistranslation risk mu(c).* Codon-specific misincorporation rates were compiled from mass-spectrometry measurements in E. coli (15). Each sense codon c was assigned a scalar mu(c) representing the probability of non-cognate amino acid incorporation per decoding event under standard growth conditions. These values describe the E. coli translation apparatus and may not generalize quantitatively to other organisms.

*Translational supply nu(c).* The tRNA Adaptation Index (tAI) served as a proxy for translational supply, reflecting the availability of cognate tRNAs for each codon (13). tAI reflects tRNA gene copy number adaptation and wobble pairing efficiency; it does not directly measure ribosome transit time at individual sites. For each codon c, tAI was computed from tRNA gene copy numbers and wobble pairing rules using the formulation of dos Reis et al. Species-specific tAI values were computed separately for each organism analyzed.

*Synonymy shielding.* Two related quantities are used. Codon-level synonymy shielding S(c) is the probability that a random single-nucleotide substitution at codon c produces a synonymous outcome (same amino acid). S(c) varies across codons; codons in the interior of large degeneracy blocks have higher S(c) than codons at block boundaries. The aggregate third-position shielding fraction S_3 averages over all 61 sense codons, considering only substitutions at the wobble position: for each codon, the three possible single-nucleotide changes at position 3 were enumerated and classified as synonymous or nonsynonymous. S_3 = 0.69 for the standard code. S(c) appears in the site-hazard factorization; S_3 characterizes the code-table-level architecture.

---

### Operational diversity metric

*Standardization.* Coordinates were standardized as z-scores across all 61 sense codons:

mu_tilde_j = (log mu_j - <log mu>) / sigma_{log mu}

nu_tilde_j = (nu_j - <nu>) / sigma_nu

*Operational spread Delta_A.* For each amino acid A with synonymous codon set C_A of size n_A:

Delta_A = [1 / C(n_A, 2)] sum_{j<k} sqrt((mu_tilde_j - mu_tilde_k)^2 + (nu_tilde_j - nu_tilde_k)^2)

---

### High-criticality site analysis

*Metal-binding residue identification.* Metal-coordinating residues were extracted from MetalPDB annotations (17) and mapped to E. coli K-12 MG1655 coding sequences via UniProt cross-references, yielding N = 17,166 metal-binding residues across 1,847 proteins.

*Enrichment testing.* For each two-codon amino acid with metal-liganding function (Cys, Asp, Glu, His), codon frequencies at metal-binding sites were compared to genome-wide background using Fisher's exact test. Enrichment was quantified by odds ratios with 95% confidence intervals.

---

### Cross-species analysis

Coding sequences were obtained from NCBI RefSeq (18): E. coli K-12 MG1655 (NC_000913), B. subtilis 168 (NC_000964), and S. cerevisiae S288C (GCF_000146045). Species-specific tAI values were computed from tRNA gene copy numbers (13).

*Slow-codon usage.* For leucine (6 codons), we computed the fraction of positions using slow synonyms (tAI < 0.25 in species-specific units).

*Diversity geometry conservation.* Delta_A was computed per amino acid in each species using species-specific tAI values. Spearman rank correlations were computed across species.

---

### Matched-null ensemble

10,000 null codes were constructed by: (1) preserving degeneracy-class structure (same number of 1-, 2-, 4-, and 6-codon amino acids); (2) randomly reassigning (mu, nu) values to codons within each degeneracy class; (3) recomputing Delta_A and total diversity for each null code. The observed diversity was compared to the null distribution.

---

## References

1. Crick, F. H. C. Codon-anticodon pairing: The wobble hypothesis. *J. Mol. Biol.* **19**, 548-555 (1966).
2. Crick, F. H. C. The origin of the genetic code. *J. Mol. Biol.* **38**, 367-379 (1968).
3. Woese, C. R. On the evolution of the genetic code. *Proc. Natl Acad. Sci. USA* **54**, 1546-1552 (1965).
4. Wong, J. T. F. A co-evolution theory of the genetic code. *Proc. Natl Acad. Sci. USA* **72**, 1909-1912 (1975).
5. Di Giulio, M. An extension of the coevolution theory of the origin of the genetic code. *Biol. Direct* **3**, 37 (2008).
6. Haig, D. & Hurst, L. D. A quantitative measure of error minimization in the genetic code. *J. Mol. Evol.* **33**, 412-417 (1991).
7. Freeland, S. J. & Hurst, L. D. The genetic code is one in a million. *J. Mol. Evol.* **47**, 238-248 (1998).
8. Tlusty, T. A model for the emergence of the genetic code as a transition in a noisy information channel. *J. Theor. Biol.* **249**, 331-342 (2007).
9. Tlusty, T. A colorful origin for the genetic code. *Phys. Life Rev.* **7**, 362-376 (2010).
10. Yarus, M., Widmann, J. J. & Knight, R. RNA-amino acid binding: A stereochemical era for the genetic code. *J. Mol. Evol.* **69**, 406-429 (2009).
11. Drummond, D. A. & Wilke, C. O. Mistranslation-induced protein misfolding as a dominant constraint on coding-sequence evolution. *Cell* **134**, 341-352 (2008).
12. Balch, W. E., Morimoto, R. I., Dillin, A. & Kelly, J. W. Adapting proteostasis for disease intervention. *Science* **319**, 916-919 (2008).
13. dos Reis, M., Savva, R. & Wernisch, L. Solving the riddle of codon usage preferences: a test for translational selection. *Nucleic Acids Res.* **32**, 5036-5044 (2004).
14. Mordret, E. *et al.* Systematic detection of amino acid substitutions in proteomes reveals mechanistic basis of ribosome errors and selection for translation fidelity. *Mol. Cell* **75**, 427-441.e5 (2019).
15. Landerer, M. *et al.* Quantitative analysis of codon-specific translation errors reveals a landscape of mistranslation hotspots and constraints. *Mol. Biol. Evol.* **41**, msae048 (2024).
16. Koonin, E. V. & Novozhilov, A. S. Origin and evolution of the universal genetic code. *Annu. Rev. Genet.* **51**, 45-62 (2017).
17. Andreini, C., Cavallaro, G., Lorber, C. & Rosato, A. MetalPDB: a database of metal sites in biological macromolecular structures. *Nucleic Acids Res.* **41**, D312-D319 (2013).
18. O'Leary, N. A. *et al.* Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation. *Nucleic Acids Res.* **44**, D733-D745 (2016).
19. Zhang, G., Hubalewska, M. & Bhardwaj, N. Transient ribosomal attenuation coordinates protein synthesis and co-translational folding. *Nat. Struct. Mol. Biol.* **16**, 274-280 (2009).
20. Plotkin, J. B. & Kudla, G. Synonymous but not the same: the causes and consequences of codon bias. *Nat. Rev. Genet.* **12**, 32-42 (2011).

---

## Figure Legends

*Figure 1. Proteostasis constraint and codon operational space.*
(A) Threshold structure in proteostasis dynamics. The misfold pool P admits a stable low-burden fixed point only when J_in < J_crit. The ODE motivates the threshold structure; parameters are not fit to data. The bifurcation is generic to systems with saturating rescue and superlinear aggregation.
(B) Synonymous codons in (mu, nu) space for E. coli. Each point represents a codon; codons for the same amino acid are connected. Operational spread Delta_A quantifies the dispersion of synonyms. mu: misincorporation risk; nu: translational supply proxy (tAI). Both axes are measured in the extant E. coli translation system.
(C) Synonymy shielding at the wobble position. Aggregate third-position shielding S_3 = 0.69: the fraction of single-nucleotide substitutions at third positions that are synonymous.

*Figure 2. Synonymy necessity constraint and cross-species conservation.*
(A) Synonymy capacity by code length. Doublet codes (16 codons) cannot encode 10+ amino acids with any synonymy. Triplet codes (64 codons) provide ~3 synonyms per amino acid. This is a parameter-free constraint independent of the proteostasis framework.
(B) Leucine codon usage across species. All species deploy slow synonyms (low tAI) at substantial frequencies (15-39%), inconsistent with monotone optimization. Species-specific preferences reflect different tRNA pools.

*Figure 3. Codon deployment at metal-binding sites.*
(A) Codon enrichment at metal-binding sites for two-codon amino acids. The enriched codon is not uniformly the lowest-mu option.
(B) Odds ratios with 95% confidence intervals. Effect sizes are modest (OR = 1.10-1.30); statistical significance reflects large sample size (N = 17,166).
(C) Enrichment patterns in (mu, nu) space. Asp and His favor the codon with higher translational supply (tAI) despite higher misincorporation risk; Cys and Glu favor the lower-mu option.

*Figure 4. Operational diversity structure.*
(A) Observed diversity exceeds matched-null expectation (z = 2.99, p = 0.0014). Deviation is moderate.
(B) Operational spread Delta_A per amino acid. Two-codon families achieve separation in (mu, nu) space; multi-codon families show variable spread.

---

## Supplementary Tables

*Table S1.* Leucine codon usage across species.

*Table S2.* Operational spread Delta_A per amino acid per species.
