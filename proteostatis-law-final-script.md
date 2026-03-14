# Proteostasis law: translation operates within a finite burden-capacity envelope

[Author names]

[Affiliations]

*To whom correspondence should be addressed. Email: [corresponding author]*

**Classification:** Biological Sciences / Systems Biology

**Keywords:** proteostasis | translation fidelity | synonymous codons | protein folding | quality control

---

## Significance Statement

Translation produces burden as well as product. As protein output rises, cells must simultaneously manage decoding errors, folding failures, aggregation, and quality-control demand. We assemble cross-organism evidence that extant translation systems operate within a finite proteostasis envelope. A reduced dynamical model with saturable rescue and positive-feedback damage predicts buffered, vulnerable, and overload regimes without requiring a universal numeric threshold. This framework explains why synonymous effects are often conditional on cellular state, why translational perturbations become more damaging when buffering is weakened, and why synonymous codons function as state-dependent load-allocation choices rather than interchangeable labels.

---

## Abstract

Translation is constrained not only by output but also by the burden it imposes on proteome integrity. Changes in translation throughput alter decoding error rates, cotranslational folding outcomes, aggregation propensity, and demand on rescue and clearance systems. We define proteostasis law as an operational constraint on translation: viability requires that the combined burden of decoding errors, folding failure, aggregation, and quality-control demand remain below the buffering and clearance capacity of the cell. Four lines of evidence support this framework. First, proteome-scale mass spectrometry shows that roughly one-fifth of *Escherichia coli* proteins contain at least one amino acid misincorporation, and codon-specific error burden spans three orders of magnitude. Second, synonymous recoding and chaperone perturbations alter folding yield, degradation susceptibility, and fitness without changing amino acid sequence. Third, weakening the Hsp70-associated chaperone network, depleting cotranslational quality-control factors, or overloading degradation machinery produces measurable fitness costs tied to protein production load. Fourth, thousands of synonymous edits in *E. coli* show fitness effects that depend strongly on growth condition rather than on fixed codon properties. A reduced burden model with finite-capacity rescue and positive damage feedback predicts buffered, vulnerable, and overload regimes, providing a physics-style interpretation of threshold-like deterioration observed experimentally. This framework supports a bounded operating-envelope view of translation in extant systems. It does not support claims about the evolutionary origin of triplet codons or a universal cross-organism proteostasis constant.

---

## Introduction

Every translated protein imposes cost beyond the resources consumed in its synthesis. Decoding errors substitute wrong amino acids. Nascent chains misfold. Misfolded species engage chaperones, proteasomes, and sequestration machinery. Aggregated material recruits still more quality-control capacity and can seed further aggregation. Because these processes share finite cellular components, translation throughput, decoding fidelity, folding success, aggregation burden, and quality-control demand are coupled variables rather than independently tunable traits (1–3).

That coupling creates a systems-level problem. High-throughput translation is beneficial for growth rate, but it also generates a proportionally larger flow of damaged products. The cell cannot maximize one without managing the other. This tradeoff has been studied gene by gene — Drummond and Wilke showed that mistranslation-induced misfolding shapes codon usage and coding-sequence evolution at individual loci (1) — but the constraint it imposes on translation as a whole has not been formalized as an operating principle.

We propose a specific formalization. **Proteostasis law** states that in extant cells, viability requires the total translation-derived burden to remain inside a finite buffering and clearance envelope. The claim is narrower than code-origin theories and wider than single-gene optimization models. It does not address why the genetic code has its present form or how many amino acids a code should encode. It addresses how a living translation apparatus avoids crossing from buffered operation into disproportionate quality deterioration.

The evidence that such a constraint operates is not new, but it is scattered across subdisciplines. Proteome-scale mass spectrometry has quantified mistranslation burden across hundreds of thousands of coding sites (4, 5). Synonymous recoding experiments have demonstrated that identical amino acid sequences can impose different folding and degradation burdens depending on codon choice and translational context (6). Chaperone genetics has shown that protein production costs are amplified when the Hsp70 network is weakened (7). Condition-dependent fitness studies have revealed that thousands of synonymous mutations in *E. coli* are neutral in one growth medium and costly in another (8). Each finding supports the general principle; none, alone, constitutes it.

The theoretical component of our proposal is also modest. A reduced dynamical model with saturable rescue and cooperative damage feedback generically produces three regimes — buffered, vulnerable, and overload — separated by fold bifurcations rather than sharp constants. This model provides a physical interpretation of threshold-like deterioration without requiring a single universal number for all organisms and conditions.

Two clarifications bound the scope from the outset. First, we use "law" in a systems sense: a proposed general constraint structure whose specific parameters are organism-dependent, not a claim that a single equation has been quantitatively validated across all species. Second, the evidence assembled here is cross-organism. Key studies span *E. coli* (4–6, 8–11), *Saccharomyces cerevisiae* (7, 12), and engineered expression systems (13). The manuscript argues for a general operating principle of extant translation, not a single-species parameterization.

---

## Theory

### The burden-capacity balance

At coarse resolution, the total translation-derived burden can be written as

> B_total = B_error + B_fold + B_agg + B_qc

where B_error denotes the decoding-error burden (amino acid misincorporation products), B_fold the folding-failure burden (kinetically trapped or misfolded intermediates), B_agg the aggregation or toxic-species load, and B_qc the demand on rescue and clearance systems. Viability requires

> B_total ≤ C_buffer

where C_buffer is the effective capacity of the proteostasis network. This inequality already implies that translation quality cannot be reduced to a single optimum on one axis; cells must operate inside a feasible region defined by load and buffering.

### Dynamical formulation

The corresponding dynamical statement tracks the effective burden pool P(t):

> dP/dt = J_err + J_fold + J_agg − R(P) − D(P) + A(P)

where J terms are inflow rates from different burden sources, R(P) and D(P) represent rescue (chaperone refolding) and clearance (proteolysis, sequestration), and A(P) captures positive feedback from aggregation or chaperone sequestration. Two generic biochemical properties govern the qualitative behavior. Rescue is saturable: chaperone capacity is finite, and proteasome throughput has a ceiling (2, 14). Aggregation introduces cooperative positive feedback: above a concentration threshold, misfolded species cross-seed, template, and sequester the quality-control machinery they depend on for clearance (3, 15).

A reduced one-variable model captures the essential structure:

> dP/dt = J − δP − V·P/(K + P) + α·P²

where J is effective inflow, δP is constitutive clearance, V·P/(K + P) is saturable rescue, and α·P² represents cooperative damage feedback. After nondimensionalization, the system is governed by three effective parameters: λ (effective inflow), ν (rescue capacity), and χ (cooperative damage strength). Such a system generically admits:

- **Buffered regime:** A single low-P attractor. Perturbations decay.
- **Vulnerable regime:** Two stable states separated by an unstable fixed point. The system is bistable; large perturbations can push it past a tipping point.
- **Overload regime:** The low-P branch is lost. Only the high-P (collapse) state remains.

The transition between buffered and vulnerable regimes occurs at a fold bifurcation (Fig. 1). The present evidence does not justify organism-specific numerical phase maps, but it does justify this qualitative stability structure.

### What the model does not do

We do not estimate R(P) or A(P) from data, nor do we fit the reduced model to any measured time course. The ODE serves as scaffolding: it justifies *why* a threshold structure (not just a slope) is expected from the known biochemistry of proteostasis. The empirical work that follows is about the components — whether the individual burden terms are real and measurable, whether buffering capacity is finite and consequential, and whether translational load and cellular state interact in the way the model predicts.

---

## Results

### Translation error burden is measurable and nontrivial

The first requirement of the framework is that B_error is substantial enough to matter for proteome integrity. Landerer et al. addressed this using a mechanistic model fitted to more than 100 mass spectrometry datasets from *E. coli* and *S. cerevisiae* (4). The central estimate is that 20–23% of synthesized proteins contain at least one amino acid misincorporation. The distribution of error burden across the codon table is far from uniform: codon-specific misincorporation rates span roughly 600-fold, from ~3 × 10⁻⁵ to ~2 × 10⁻² per decoding event (4, 5). Translation errors are not rare events that can be treated as a perturbation; they are a quantitatively significant and heterogeneous flow of damaged products.

This flow connects directly to downstream proteotoxic material. Ling et al. showed that aminoglycoside-induced mistranslation drives transient protein aggregation in *E. coli* (10). The experiment is a perturbation — drug-treated cells do not recapitulate normal physiology — but it demonstrates the causal chain: elevated decoding error → increased misfolding → observable aggregation. Evans et al. extended this by showing that increased mistranslation in *E. coli* engages a RpoS-dependent heat-shock response, confirming that decoding burden is coupled to stress-management machinery rather than remaining an isolated molecular defect (11).

### Folding burden depends on translational context

If translation burden were fully determined by the amino acid sequence, synonymous codon choice would be neutral for protein fate. It is not. Walsh et al. demonstrated that synonymous recoding of the *Neurospora crassa* clock protein FRQ perturbs cotranslational folding *in vivo*, increases degradation susceptibility, and impairs fitness in *E. coli* (6). The recoded protein has an identical amino acid sequence; the fitness cost flows through folding and downstream handling, not through a change in protein function.

The speed–yield tradeoff is explicit in the chaperone literature. Agashe et al. showed that trigger factor and DnaK increase the folding yield of multidomain substrates at the expense of folding speed (9). Chaperone-assisted folding is not free: it consumes time and quality-control capacity. A translation system that generates more misfolded intermediates per unit time draws more heavily on these limited resources, even if the final folded product is identical.

These results establish B_fold as a genuine component of translational burden. The cost of translation cannot be inferred from protein output alone, because the same amino acid sequence imposes different folding and quality-control demands depending on how it is synthesized.

### Buffering and clearance capacity are finite and consequential

The strongest evidence in this package concerns the right-hand side of the inequality: C_buffer is not an abstraction but an experimentally accessible quantity with measurable consequences when it is exceeded.

Farkas et al. provided the most direct test (7). Working in *S. cerevisiae*, they showed that weakening the Hsp70-associated chaperone network amplifies the cost of protein production. Strains with reduced Ssb (ribosome-associated Hsp70) function pay a disproportionate fitness price for high expression of aggregation-prone proteins. The cost is specific to burden: it does not scale with protein output per se but with the product of output and aggregation propensity. In the language of the theory, reducing C_buffer makes the same J_in more costly.

Chuang et al. demonstrated that damaged nascent proteins generate a real cotranslational clearance demand handled by the proteasome and translation-elongation-associated quality-control factors (12). This is not a hypothetical pathway. Cotranslational degradation is a measurable flux that competes for proteasome capacity with other substrates.

De Marco et al. showed the reciprocal: overexpressing chaperones (GroEL/GroES, DnaK/DnaJ/GrpE) in *E. coli* rescues the solubility of recombinant proteins that are otherwise insoluble at high expression levels (13). Expanding C_buffer allows the system to tolerate a higher J_in before quality collapse. These three studies are not numerically commensurate — they use different organisms, readouts, and perturbation strategies — but they converge on a single structural conclusion. Rescue and clearance are finite, biologically consequential, and load-bearing.

### Synonymous effects are conditional on cellular state

The reduced model predicts that the consequence of a given perturbation to J_in depends on the operating point. A synonymous change that is neutral in a well-buffered state can become costly when the system is near the vulnerable–overload boundary. Yang et al. tested this prediction at scale, introducing thousands of synonymous edits into the *E. coli* genome and measuring fitness under two growth conditions: glucose (well-buffered) and acetate (resource-constrained) (8). The results are striking: synonymous effects are strongly condition-dependent, with substantially larger fitness consequences in acetate than in glucose. Many mutations are neutral in one condition and significantly deleterious in the other.

Walsh et al. provide mechanistic support (6). Synonymous changes that alter cotranslational folding increase degradation susceptibility, but the magnitude of the fitness consequence depends on the cellular context in which the protein is expressed. A codon that is "safe" under favorable chaperone availability can become "risky" when folding support is limiting.

The implication is that synonymous codons are best described as load-allocation choices inside the burden-capacity envelope rather than as fixed-property labels. A given codon may favor lower error burden at the cost of reduced throughput, or higher throughput at the cost of larger folding and quality-control demand. The same codon-level perturbation can be mild in a buffered state and costly when the system operates near capacity (Fig. 2).

### Evidence for nonlinearity

The framework predicts that degradation of proteostasis quality should not be linear in burden: as the system approaches the fold bifurcation, small increases in load produce disproportionately large effects. The evidence is consistent with this prediction but does not prove a universal bifurcation law.

Martens and Hilser reported that synonymous codon changes can alter reporter output by up to ~70-fold while the chaperone heat-shock response increases only ~2-fold (16). The disproportion is expected from saturable buffering: once chaperones approach capacity, even a modest additional misfolding load escapes rescue.

Farkas et al. observed that the fitness cost of protein production increases sharply when combined with translational or folding stress — a synergistic interaction consistent with operating near a capacity boundary (7). The cost of a fixed expression load is not additive with background stress; it amplifies.

The manuscript should be read as claiming that threshold-like worsening is *expected* from finite-capacity rescue plus positive damage feedback, and that current data *partly* support this view. It should not be read as establishing a single cross-organism threshold constant or a measured phase diagram.

---

## Testable Predictions

The value of the framework extends beyond retrospective synthesis. Four predictions follow directly from the burden-capacity structure and can be tested with existing experimental tools.

**Prediction 1.** For a fixed coding sequence, synonymous variants that are tolerated in a well-buffered state should show larger fitness or quality defects when chaperone or clearance capacity is experimentally reduced. The specific prediction is an interaction: the fitness cost of a synonymous change × the degree of proteostasis impairment should be superadditive. Farkas et al. (7) have the experimental system to test this; extending their approach to defined synonymous variants would be direct.

**Prediction 2.** Increasing translation output without proportionally increasing rescue capacity should move cells toward a vulnerable regime in which small additional perturbations produce disproportionate losses in folding yield or fitness. Titrating inducer concentration while measuring both expression level and aggregation fraction in a fixed genetic background would trace this trajectory.

**Prediction 3.** Highly expressed or intrinsically difficult-to-fold proteins should show steeper condition-dependence of synonymous effects than low-expression or easy-to-fold proteins. The dataset of Yang et al. (8), combined with published aggregation propensity and expression-level data for *E. coli*, is already sufficient to test this prediction by stratifying their fitness measurements.

**Prediction 4.** Direct measurement of burden proxies (aggregation fraction, proteasome occupancy) and buffering proxies (free chaperone pool, clearance rate) in the same system should reveal a regime boundary more effectively than codon-only or expression-only analyses. Dual-reporter systems that simultaneously track protein output and aggregation load under titratable stress would provide the two-dimensional data the model requires.

---

## Discussion

The claim advanced here is conceptual but evidence-backed: translation output and proteome integrity are jointly constrained by finite buffering capacity. This constraint operates as an envelope, not as a single threshold constant. Its parameters — chaperone abundance, proteasome throughput, aggregation cooperativity — vary across organisms and conditions. What is general is the structure of the constraint, not a universal number.

### What the framework explains

Three features of translation biology acquire a unified interpretation under proteostasis law.

First, synonymous effects are conditional. A codon change that is neutral in one growth condition can be costly in another (8), not because the codon's molecular properties change, but because the cell's distance from its capacity boundary changes. The framework predicts the interaction term directly: codon effect × state vulnerability.

Second, translational perturbations become more damaging in weakly buffered states. Farkas et al. showed this for chaperone-depleted yeast (7); Ling et al. showed this for drug-induced mistranslation in *E. coli* (10). These are not different phenomena. They are different perturbations to the same coupled system.

Third, cells do not exploit all nominal coding freedom. Synonymous codons are not interchangeable labels precisely because different codons impose different loads along error, folding, and quality-control axes. The operating envelope constrains which synonymous choices are viable, and this constraint tightens as the cell's buffering margin shrinks.

### Relation to existing frameworks

Drummond and Wilke established that mistranslation-induced misfolding is a dominant selective pressure on coding sequences at the single-gene level (1). We extend their logic to the system level: the same costs that shape intra-gene codon usage also constrain how the entire translation apparatus must operate. The distinction is between optimization of individual genes (their contribution) and a feasibility constraint on the translation system as a whole (ours).

Balch et al. formalized proteostasis as a network property — a balance among synthesis, folding, trafficking, and degradation (2). Our contribution is to connect that network concept to the specific problem of translation-derived burden and to show that the resulting structure admits distinct operating regimes rather than a single optimum.

Error-minimization theories of the genetic code show that the standard codon table is unusually robust to point mutations and mistranslation (17, 18). This is compatible with our framework: a code that minimizes the phenotypic cost of errors also reduces J_in per unit of throughput, which helps keep the system inside the buffered regime. The two perspectives are complementary, not competing.

### What the framework does not explain

Proteostasis law does not address the evolutionary origin of the genetic code — neither why triplet codons are universal, nor why the standard code has its particular mapping of codons to amino acids. Those questions require historical and structural arguments beyond the scope of an operating-envelope framework for extant systems.

The framework also does not predict a universal threshold constant. We expect C_buffer to differ between organisms, growth phases, temperatures, and stress conditions. The prediction is structural (a finite envelope exists and has consequences) rather than parametric (the envelope has a specific numeric boundary).

### Cross-organism scope and its limits

The evidence assembled here is intentionally cross-organism. The error-burden measurements are from *E. coli* (4, 5); the chaperone-buffering genetics are from *S. cerevisiae* (7); the condition-dependent synonymous effects are from *E. coli* (8); the cotranslational clearance demand is from *S. cerevisiae* (12). This breadth supports the generality of the principle but limits the precision of any organism-specific quantitative claims. A fully parameterized burden-capacity model for any single organism will require matched measurements of error burden, chaperone capacity, and proteostasis output in the same system under the same conditions.

### Closing

Translation is a burden-generating process. That burden is managed by a finite network of rescue and clearance systems. When burden exceeds capacity, quality deteriorates — not linearly, but through the threshold-like dynamics that finite-capacity rescue and cooperative damage feedback produce. The evidence assembled here supports this constraint structure as a general operating principle. What it does not support is a universal numeric threshold or an origin story for the genetic code. The operating envelope is real and consequential, and it should be treated as such.

---

## Methods

### Scope and evidence selection

This manuscript assembles published evidence rather than presenting new primary data. Studies were selected to support or constrain the four core claims: that translation error burden is measurable, that folding burden depends on translational context, that rescue/clearance capacity is finite and consequential, and that synonymous effects are state-dependent. Studies were drawn from proteomics (4, 5), synonymous recoding experiments (6, 8), chaperone genetics (7, 9, 13), quality-control biochemistry (12), and aminoglycoside perturbation experiments (10, 11).

### Reduced dynamical model

The reduced burden model is

> dP/dt = J − δP − V·P/(K + P) + α·P²

where P is the effective burden-pool concentration, J the effective inflow, δ the constitutive clearance rate, V and K the maximal velocity and half-saturation constant for saturable rescue (chaperone/proteasome), and α the cooperative aggregation rate. Nondimensionalization yields dimensionless parameters λ = J/(δK), ν = V/(δK), and χ = αK/δ. Steady states satisfy

> λ = p + ν·p/(1 + p) − χ·p²

where p = P/K. Fold bifurcations (transitions between buffered and vulnerable regimes) occur where dλ/dp = 0, corresponding to simultaneous tangency of the inflow and total-sink curves. The model is analyzed qualitatively; no parameters are fitted to experimental data.

### Codon-specific error data

The statement that codon-specific error burden spans ~600-fold is from Landerer et al. (4), who compiled codon-specific misincorporation rates from reanalysis of mass spectrometry datasets. The estimate that 20–23% of proteins contain at least one misincorporation is from the same study's mechanistic model fits. We did not recompute these values; they are cited as published.

### Data and code availability

This manuscript presents no new primary data. All cited data are available through the original publications. The reduced model equations are presented in full in the Theory section; no code beyond standard numerical root-finding is required to reproduce the qualitative stability analysis.

---

## Figure Legends

**Figure 1. Stability structure of the burden-capacity model.** (**A**) The reduced model dP/dt = J − δP − VP/(K + P) + αP² produces buffered, vulnerable, and overload regimes. In the buffered regime (low J), a single stable low-P attractor exists. As J increases, a fold bifurcation creates a second (high-P) stable state; the system becomes vulnerable to large perturbations. Beyond a critical inflow, the low-P attractor is lost and only the high-P (collapse) state remains. (**B**) Qualitative phase portrait in λ–ν space at fixed χ, illustrating the dependence of regime boundaries on rescue capacity (ν). (**C**) Summary of the reduced equations and fold conditions. This is a conceptual figure, not a fitted phase diagram.

**Figure 2. Synonymous codons as state-dependent load-allocation choices.** Synonymous codons are not fixed-purpose labels. Codon choice determines the distribution of translational burden across error (misincorporation risk), folding (cotranslational misfolding propensity), and quality-control (rescue and clearance demand) axes. The same codon family can have different fitness consequences in buffered versus stressed states, making synonymous effects conditional on the burden-capacity context. Schematic illustrates how the operating point shifts relative to the capacity boundary under different cellular states.

---

## References

1. D. A. Drummond, C. O. Wilke, Mistranslation-induced protein misfolding as a dominant constraint on coding-sequence evolution. *Cell* **134**, 341–352 (2008).
2. W. E. Balch, R. I. Morimoto, A. Dillin, J. W. Kelly, Adapting proteostasis for disease intervention. *Science* **319**, 916–919 (2008).
3. M. S. Hipp, P. Kasturi, F. U. Hartl, The proteostasis network and its decline in ageing. *Nat. Rev. Mol. Cell Biol.* **20**, 421–435 (2019).
4. C. Landerer, J. Poehls, A. Toth-Petroczy, Fitness effects of phenotypic mutations at proteome-scale reveal optimality of translation machinery. *Mol. Biol. Evol.* **41**, msae048 (2024).
5. E. Mordret *et al.*, Systematic detection of amino acid substitutions in proteomes reveals mechanistic basis of ribosome errors and selection for translation fidelity. *Mol. Cell* **75**, 427–441.e5 (2019).
6. I. M. Walsh, M. A. Bowman, I. F. Soto Santarriaga, A. Rodriguez, P. L. Clark, Synonymous codon substitutions perturb cotranslational protein folding in vivo and impair cell fitness. *Proc. Natl. Acad. Sci. U.S.A.* **117**, 3528–3534 (2020).
7. Z. Farkas *et al.*, Hsp70-associated chaperones have a critical role in buffering protein production costs. *eLife* **7**, e29845 (2018).
8. D. D. Yang, L. M. Rusch, K. A. Widney, A. B. Morgenthaler, S. D. Copley, Synonymous edits in the *Escherichia coli* genome have substantial and condition-dependent effects on fitness. *Proc. Natl. Acad. Sci. U.S.A.* **121**, e2316834121 (2024).
9. V. R. Agashe *et al.*, Function of trigger factor and DnaK in multidomain protein folding: Increase in yield at the expense of folding speed. *Cell* **117**, 199–209 (2004).
10. J. Ling, C. Cho, L. T. Guo, H. R. Aerni, J. Rinehart, D. Söll, Protein aggregation caused by aminoglycoside action is prevented by a hydrogen peroxide scavenger. *Mol. Cell* **48**, 713–722 (2012).
11. C. R. Evans, J. Fan, X. Ling, Increased mistranslation protects *E. coli* from protein misfolding stress due to activation of a RpoS-dependent heat shock response. *FEBS Lett.* **593**, 3220–3227 (2019).
12. S. M. Chuang, L. Chen, D. Lambertson, M. Anand, T. G. Kinzy, K. Madura, Proteasome-mediated degradation of cotranslationally damaged proteins involves translation elongation factor 1A. *Mol. Cell. Biol.* **25**, 403–413 (2005).
13. A. de Marco, E. Deuerling, A. Mogk, T. Tomoyasu, B. Bukau, Chaperone-based procedure to increase yields of soluble recombinant proteins produced in *E. coli*. *BMC Biotechnol.* **7**, 32 (2007).
14. F. U. Hartl, A. Bracher, M. Hayer-Hartl, Molecular chaperones in protein folding and proteostasis. *Nature* **475**, 324–332 (2011).
15. J. Tyedmers, A. Mogk, B. Bukau, Cellular strategies for controlling protein aggregation. *Nat. Rev. Mol. Cell Biol.* **11**, 777–788 (2010).
16. A. T. Martens, V. J. Hilser, Chaperone saturation mediates translation and protein folding efficiency. *bioRxiv* 2025.06.25.661590 (2025). Preprint.
17. D. Haig, L. D. Hurst, A quantitative measure of error minimization in the genetic code. *J. Mol. Evol.* **33**, 412–417 (1991).
18. S. J. Freeland, L. D. Hurst, The genetic code is one in a million. *J. Mol. Evol.* **47**, 238–248 (1998).
