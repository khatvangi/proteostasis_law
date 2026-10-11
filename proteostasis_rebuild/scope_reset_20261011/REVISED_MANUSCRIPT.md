# Translation-associated proteostasis burden: theoretical reassessment and experimental constraints

Hannah E. Rebbeck, Janak Paudyal, and Boggavarapu Kiran

Department of Chemistry and Physics, McNeese State University, Lake Charles, LA 70609, USA

Corresponding author: kiran@mcneese.edu

**Status: theoretical reassessment and experimental constraints / not submission-ready.** Reconstruction completed 11 October 2026. No new biological experiment or fitted physiological model is reported. Claim identifiers C01–C22 refer to `CLAIM_AND_EVIDENCE_MATRIX.tsv`; literature identifiers are resolved in References.

## Abstract

Translation generates both ordinary folding demand and potentially additional processing associated with errors. The original proteostasis-envelope proposal asked whether finite processing capacity constrains this demand and whether experimentally attainable perturbations expose that constraint. Reassessment shows that its numerical physiological margin and generic collapse interpretation are unsupported. We reconstruct the question using counterfactual error-attributable outcome flux, machine-specific occupied site-time, and a conservative finite-pool realization. Conservation and steady occupancy impose necessary conditions but establish neither successful routing nor adequate function. Regulation can alter machinery abundance, yet neither a stress reporter nor abundance alone measures productive throughput. Published experiments support conditional chaperone cooperation, client-dependent effects and regulated buffering. They do not supply a matched measurement of translation-error-attributable workload. A previously constructed kinetic model also fails the point-estimate direction of an external insolubility comparison: deleting htpG in a ΔdnaKJ background reduces the reported insoluble fraction from 20.5% to 18.9%, whereas the sampled model usually predicts an increase. Uncertainty in that total contrast and matched growth remain unresolved. The surviving contribution is a transparent accounting and experimental constraint framework, with competing predictions for whether a calibrated error perturbation increases machine workload. It is not a validated capacity law or a submission-ready empirical paper.

## Significance

A protein substitution, a misfolded chain, an occupied chaperone and a loss of activity are different observations. Conflating them can manufacture a quantitative capacity margin from measurements that do not identify one. Separating these levels preserves the original experimental question while making its failure modes explicit. The framework specifies what must be measured to connect translation errors to processing and function, and records a negative model confrontation instead of treating agreement with selected observations as validation. Finite pools, kinetic competition and stress feedback have substantial prior art; the present contribution is their disciplined use in reassessing this particular translation-burden claim [FE12, SD17, FM17, FS19].

## Introduction

The canonical manuscript, *A finite proteostasis envelope for translation*, posed two questions: “First, is the envelope real and finite, and where does extant translation sit inside it? Second, does that position have consequences an experiment can reach?” [C01]. This is the manuscript's wording, not a recovered quotation from an earlier conversation. Its physical motivation remains legitimate: newly synthesized proteins enter folding and quality-control processes whose molecular components are finite. Translation errors can change the proteins entering those processes. However, a finite molecular pool does not by itself locate a cell relative to failure.

The original argument combined an aggregate burden, a reduced collapse model and literature-anchored numerical bounds. Subsequent audits removed the error-rate anchor, the physiological headroom estimate and the claim that generic aggregation feedback eliminates all bounded states. They also rejected a synonymous “shielding” factor for ribosomal mistranslation and the assumption that every amino-acid substitution is a misfolding event [C13, C16, C17]. These are changes in scientific content, not adjustments to confidence intervals.

We retain a bounded question: at measured protein-specific synthesis rates in a specified strain, medium, temperature and growth phase, what extra processing is attributable to a specified translation-error perturbation, which machinery does it occupy, and what functional change follows when machinery is altered? Ordinary folding and degradation form the baseline. Mature-protein damage is a separate source of demand. Synonymous variants are secondary interventions that can affect timing or error probability; neither mediator is presumed.

Published models already address extensive parts of this problem. FoldEco describes kinetic folding, aggregation, chaperoning and degradation [FE12]. Santra and colleagues incorporate growth and client-dependent folding properties [SD17]. FoldME constrains folding fluxes within a growth-optimization framework [FM17]. FoldEcoSlim includes a background proteome, growth dilution and a dynamic σ32-linked response [FS19]. We therefore claim neither a new finite-capacity principle nor novelty from adding feedback. Our purpose is to identify which statements survive the local evidence audit and which measurements could discriminate the remaining interpretations.

## Model derivation and equations

### Outcome flux and incremental workload

Let (s_i) be synthesis of client i, in chains per volume per time. Fix the condition and observation interval. Let (\pi_i(e)) be the probability that a newly synthesized chain has error pattern e, with e=0 denoting no error. For a specified outcome Y, let (p_{iY}(e)) be its probability under an intervention imposing e in a common processing environment. This is a causal estimand, not the frequency of errors among surviving proteins. Then

\[
J_Y=\sum_i s_i\sum_e\pi_i(e)p_{iY}(e),\qquad
J_Y^0=\sum_i s_i p_{iY}(0),
\tag{1}
\]
\[
\Delta J_{\mathrm{error},Y}=\sum_i s_i\sum_{e\ne0}\pi_i(e)[p_{iY}(e)-p_{iY}(0)].
\tag{2}
\]

Normalization of π gives equation (2) by subtraction [C02]. Error-free chains already fold, bind machinery and undergo degradation. Consequently neither a zero baseline nor unit probability of error-induced misfolding is justified. An error mixture could even reduce a particular processing outcome. Comparing induced and control mixtures by replacing π with its difference is legitimate only when synthesis and conditional outcomes remain unchanged. If competition or regulation changes those outcomes, the full expressions must be compared separately.

Outcome flux has units of chains per volume per time. To express machine demand, define (w_{ij}(e)) as expected cumulative occupied site-time on machine j per synthesized chain, including repeated visits. At stationarity with finite mean residence and a consistent treatment of exits,

\[
B_j=\sum_i s_i\sum_e\pi_i(e)w_{ij}(e),\qquad
\Delta B_{\mathrm{error},j}=\sum_i s_i\sum_{e\ne0}\pi_i(e)[w_{ij}(e)-w_{ij}(0)].
\tag{3}
\]

The units are occupied sites per volume [C03]. A vector over machines replaces a sum of unlike quantities. Error frequency, folding flux, aggregate concentration and occupancy cannot be added into one physical budget. Occupancy from background clients must be included when comparing B with an experimentally measured total pool. Equation (3) does not identify its own inputs, and a fixed-environment increment is not a universal additive decomposition when the environment responds.

### A minimal conservative realization

For one machine and multiple clients, let U, C, A and N denote free non-native, machine-bound, aggregated and native chains, respectively, all in chain-equivalent concentration. Let F be free site concentration, M total sites, and (F+\sum_iC_i=M). One bound chain occupies one site; conversion from protein subunits or oligomers to sites requires an explicit experimental convention. Let q be the fraction of synthesis entering U. Set (b_i=k_{\mathrm{on},i}FU_i) and (a_i=k_{a,i}U_i^2). Aggregation a is a transfer of chain equivalents, not a rate of creation of aggregate particles.

\[
\begin{aligned}
\dot U_i&=q_is_i+k_{u,i}N_i+r_iC_i-b_i-k_{f,i}U_i-d_{U,i}U_i-a_i-\mu U_i,\\
\dot C_i&=b_i-(r_i+f_i+d_{C,i}+\mu)C_i,\\
\dot A_i&=a_i-(d_{A,i}+\mu)A_i,\\
\dot N_i&=(1-q_i)s_i+k_{f,i}U_i+f_iC_i-(k_{u,i}+d_{N,i}+\mu)N_i,\\
\dot F&=\alpha-\mu F-\sum_i b_i+\sum_i(r_i+f_i+d_{C,i})C_i,\\
\dot M&=\alpha-\mu M.
\end{aligned}
\tag{4}
\]

Here μ is measured dilution, α supplies new free sites, and all rates are nonnegative [C04]. Capture and quadratic aggregation constants have inverse-concentration/time units; other rate constants have inverse-time units. Fixed M requires replacement (\alpha=\mu M). Machine degradation beyond dilution is omitted, and client degradation releases its site. The d terms are conditional effective loss rates: they do not establish unlimited protease capacity. Resolving protease competition would require additional measured pools, which are outside this reconstruction.

Let (P_i=U_i+C_i+A_i+N_i), and introduce bookkeeping counters with (\dot D_i=\sum_Xd_{X,i}X_i) and (\dot L_i=\mu P_i). Summation gives

\[
\dot P_i=s_i-\mu P_i-\sum_Xd_{X,i}X_i,\qquad
\frac{d}{dt}(P_i+D_i+L_i)=s_i.
\tag{5}
\]

All transfers cancel. The counters are not available client pools. At each zero-state boundary the outward loss terms vanish, leaving nonnegative inflow; F=0 stops capture. Moreover, equation (4) preserves (F+\sum C=M). These properties establish positivity and finite occupancy for admissible initial states. For constant positive μ and bounded synthesis, equation (5) bounds total client concentration. Aggregation alone cannot imply universal runaway in this closure. At zero dilution, boundedness needs an independently established removal mechanism [C05].

### Stationarity, function and response

At a steady state with finite residence,

\[
C_i=b_i\tau_i,\quad \tau_i=(r_i+f_i+d_{C,i}+\mu)^{-1},\quad
\sum_i b_i\tau_i\le M.
\tag{6}
\]

This is the occupancy identity associated with Little's law; binding flux counts visits [C06]. Productive formation satisfies

\[
J_{N,i}=(1-q_i)s_i+k_{f,i}U_i+f_iC_i
=(\mu+k_{u,i}+d_{N,i})N_i.
\tag{7}
\]

Only if an independent assay requires (N_i\ge N_{\mathrm{req},i}) does the corresponding lower bound on (J_{N,i}) follow [C07]. Native conformation need not confer native activity after mistranslation. An activity assay is therefore required; a death cutoff cannot be supplied by the model. Occupancy feasibility does not guarantee that kinetic routing achieves sufficient activity, and stationary equations cannot be imposed on an unrelaxed transient.

A limited regulatory closure is (\dot m=b(S)-\gamma_m m), (\alpha=k_t m), with measured active σ32 signal S, machinery transcript m, transcript loss γm and translation coefficient kt [C08]. New machinery enters F. This is an unfitted response to measured input, not a prediction of S. Increasing S does not necessarily increase α unless b responds and transcription, translation and assembly remain available. Nor does increased turnover imply increased native output. The closure applies only to regulated machinery; trigger factor is outside this regulon [T98, G04, L13, Z14]. We introduce no delay theorem or physiological response ceiling.

## Results

### Derived mathematics

The mathematical results are conservation, positivity, dimensional consistency, counterfactual subtraction and the stationary identities in equations (1)–(7). They are conditional deductions, not biological observations. `MATHEMATICAL_CHECKS.py` verifies them using exact rational arithmetic and symbolic algebra where available; its complete run is in `CHECK_RESULTS.txt`. The script also constructs **illustrative model consequences** in which error mixtures increase, leave unchanged or decrease occupied site-time, and in which capture by a nonproductive route lowers native yield. These witnesses disprove automatic sign claims; their parameter values are not estimates for cells [C09].

The earlier corrected sequestration reduction admits multiple finite states under particular assumptions. It supplies no universal collapse theorem and is not used to locate a physiological boundary here [C17]. Similarly, the earlier constant-demand argument distinguishes attainable supply from realized supply, but its whole-function demand construction and example numbers are not results of this translation-focused manuscript. A generic Hall-type subset count is not adopted: serial routes can consume more than one machine, invalidating that simplified criterion [C19].

### Verified retrospective experimental patterns

The source records support conditional cooperation. GroEL/GroES overproduction can rescue the documented trigger-factor/DnaK-deficient background under specified conditions [V04]. Thus obligatory serial passage through TF or DnaK cannot be imposed on every productive GroEL route. This observation informed model construction and is not independent validation. The V04 record carries a duplicate-publication flag whose bibliographic resolution remains a submission blocker. Calloni and colleagues report altered DnaK association after TF loss and extensive aggregation after combined chaperone loss [C2012]. Protein-identity counts and binding enrichment establish neither mass flux nor productive obligation [C10].

The σ32 literature supports regulation of chaperone activity and abundance [T98, G04, L13]. Engineered activation through feedback-deficient σ32 increases regulated machinery and improves selected recombinant-protein yields [Z14]. It establishes an attainable intervention in that experimental system, not the native feedback gain, functional assembly delay, or a maximum endogenous capacity [C11]. Abundance measurements likewise constrain available molecular amounts, not in-vivo service times [S2016].

For synonymous sequence changes, Reingewertz and colleagues report altered GroEL dependence of mouse DHFR variants expressed in E. coli [R24]. The translation-rate interpretation is computational, rather than a direct rate measurement. This is a specific heterologous-client result, not a calibrated proteome-wide demand weight. Separately, the archived reanalysis of Landerer data retains a depth-calibrated AGC substitution contrast under its specified exponential-phase pipeline [LA24]. We retain its qualitative existence, not a portable fold-change amplitude. The pipeline lacks a decoy-controlled localization error assessment and supplies no folding or occupancy readout. Neither result demonstrates a synonymous-error-to-workload causal chain [C12].

### Failed tests and unsupported extrapolations

The original error-rate estimate averaged detected substitutions while excluding covered zero-detection cells. Its use as an unconditional synthesis-error rate was invalid. Broad synonymous comparisons also failed a depth-matched null calibration. Corrected standing-frequency summaries do not repair the missing birth probabilities or provide a new headroom estimate [C13]. The prior all-client recovery normalization forced every aggregate through a small ClpB pool, leaving the intended interaction untested; it is excluded [C14].

The most direct negative confrontation concerns Fauvet's 30 °C data [H2021]. Reported insoluble protein decreases from 20.5% in ΔdnaKJ to 18.9% in ΔdnaKJ ΔhtpG, a 1.6-percentage-point decline. The earlier kinetic sweep instead increased aggregate after htpG deletion in 200/200 single-client and 199/200 two-client configurations. These counts describe sampled models, not biological replication or a universal theorem. The test fails at the point estimate [C15].

Its statistical interpretation remains open. The insoluble fraction includes membrane proteins; modeled aggregate is a different observable. Later-retrieved tables contain per-protein means and standard deviations across five biological replicates, correcting the earlier blanket statement that the supplement was unavailable. They do not provide the replicate-level total contrast or the cross-protein covariance needed to propagate uncertainty in a sum. Growth rescue at 37/39 °C does not establish dilution at 30 °C. These limitations prohibit a claim of statistically established refutation, but equally prohibit presenting the model as validated. Adjusting it to this endpoint would turn the endpoint into calibration data, not preserve a held-out test.

## Discussion

The original question remains experimentally meaningful after its numerical answer is withdrawn. Finite molecular pools limit occupancy; whether a particular translation perturbation materially uses those pools depends on the clients produced, their kinetic histories and the measured environment. A stress response is compatible with altered demand but does not identify its cause. Improved solubility is compatible with degradation as well as productive folding, and an increase in native abundance can coexist with altered specific activity. Joint readouts are required to distinguish these interpretations.

The strongest surviving contribution is consequently methodological: distinguish error events, client outcomes, site-time and function; formulate an explicit baseline; and retain contradictions when interpreting a model. These are established physical and causal principles applied to a failed quantitative synthesis. FoldEcoSlim makes especially clear that background-proteome competition and stress induction are already published [FS19]. Neither their combination nor a rewritten occupancy inequality establishes novelty. A publishable contribution would need a discriminating empirical result or a separately justified methodological advance [C18].

The evidence is also narrower than “chaperones are protective.” Routes can divert clients to disposal, cooperate conditionally, or compete for shared pools. More machinery need not improve every client. Conversely, loss of a route does not by itself disprove finite capacity when growth improves. Optimizing models and realized kinetic models answer different questions, and published optimization models can include costs that alter deletion effects. We make no blanket rejection of FoldME based on a simple route-removal counterexample.

## Methods and reproducibility

This is a records-based synthesis of the canonical manuscript, scope-reset documents, reassessment, work-order audit, corrected theory and kinetics records, and locally archived primary literature. No new literature search, raw-spectrum reanalysis, kinetic fit or experiment is represented as completed. The original manuscript is preserved byte-for-byte as `ORIGINAL_FROZEN_MANUSCRIPT.md`; the source manifest records its hash and the hashes of evidence anchors. The canonical author order and affiliation are retained without implying author approval of this revision.

`CLAIM_AND_EVIDENCE_MATRIX.tsv` links manuscript claim IDs to source paths, exact line intervals and literature IDs, distinguishing observed, derived, conditional, unsupported and withdrawn claims. “Observed” denotes published or archived retrospective evidence, not new observation by these authors. `SOURCE_MANIFEST.tsv` identifies the files used. Reproduce the mathematical checks from this directory with `python3 MATHEMATICAL_CHECKS.py > CHECK_RESULTS.txt 2>&1`. The checks do not rerun older simulation ensembles; those counts are quoted from their saved results. Mathematical success and empirical validation are separate verdicts [C20].

## Prospective discriminating tests

One bounded experiment addresses the primary gap [C21]. In matched cultures, apply a calibrated translation-error perturbation and an orthogonal fidelity rescue, measure nascent error spectra and protein-specific synthesis, and track client-specific occupied site-time, aggregation, turnover and activity. Cross the error intervention with independently controlled functional machinery expression. Measure early occupancy before induction and later active σ32 and assembled machinery. Prespecify an error increment, a minimum meaningful workload increment δB, an activity effect, the time window, and replicate-based uncertainty. Those values require assay information unavailable here; they are not invented numerical thresholds.

**Competing prediction A: additional processing is the mediator.** In the chosen condition, verified extra errors produce a reproducible positive ΔB above δB on the prespecified machine; machinery manipulation changes the associated activity loss with errors and synthesis held comparable. Fidelity rescue reverses the error-linked component. A sufficiently precise upper confidence bound below δB falsifies this material-workload prediction for the specified perturbation. An inseparable drug or synthesis effect makes the causal test uninterpretable, not positive.

**Competing prediction B: the perturbation has no material incremental workload on that machine.** Verified extra errors leave ΔB below δB because altered products have similar residence demand or leave by other routes; a reporter or activity change may still occur. A reproducible lower confidence bound above δB, reversed by fidelity rescue with synthesis and other injury controlled, falsifies this prediction. An interval crossing δB decides neither. This alternative does not assert that errors are harmless or that finite pools are absent.

The factorial interaction, (I=y_{11}-y_{10}-y_{01}+y_{00}), must use a prespecified readout and scale. Its sign is secondary: nonzero I alone does not demonstrate saturation, a physiological threshold or synthetic lethality. The separate Fauvet contradiction remains visible and unresolved; resolving its total-fraction uncertainty would constrain that kinetic model but would not replace the primary error-workload experiment [C22]. Scientific reconstruction is complete. New wet-lab validation is unavailable, and this manuscript remains not submission-ready.

## References

Identifiers below come from the existing verified ledgers or archived primary-source manifests; source-access limitations remain in the evidence matrix.

- **FE12** — Powers DL, Powers ET, Gierasch LM. FoldEco model. *Cell Reports* 1:265–276 (2012). DOI: 10.1016/j.celrep.2012.02.011; PMID: 22509487.
- **SD17** — Santra M, Farrell DW, Dill KA. Chaperone-mediated folding and growth model. *PNAS* 114:E2654 (2017). DOI: 10.1073/pnas.1620646114; PMID: 28292901.
- **FM17** — Chen et al. FoldME. *PNAS* 114:11548 (2017). DOI: 10.1073/pnas.1705524114; PMCID: PMC5664499.
- **FS19** — Pobre KFR, Powers DL, Ghosh K, Gierasch LM, Powers ET. Kinetic versus thermodynamic control of mutational effects on protein homeostasis: A perspective from computational modeling and experiment. *Protein Science* 28:1324–1339 (2019). DOI: 10.1002/pro.3639; PMID: 31074892; PMCID: PMC6566562. Erratum PMCID: PMC6863698 (acknowledgment correction).
- **V04** — Vorderwülbecke et al. GroEL/GroES rescue of TF/DnaK deficiency. *FEBS Letters* 559:181–187 (2004). DOI: 10.1016/S0014-5793(04)00052-3; PMID: 14960329. Duplicate-publication flag unresolved; not counted as independent replicate evidence.
- **C2012** — Calloni et al. DnaK client network. *Cell Reports* 1:251 (2012). DOI: 10.1016/j.celrep.2011.12.007; PMID: 22832197.
- **H2021** — Fauvet et al. Bacterial Hsp90 Facilitates the Degradation of Aggregation-Prone Hsp70–Hsp40 Substrates. *Frontiers in Molecular Biosciences* 8:653073 (2021). DOI: 10.3389/fmolb.2021.653073; PMCID: PMC8082187; proteomics accession: PXD010014.
- **R24** — Reingewertz et al. Mouse DHFR substitutions and GroEL dependence. *Protein Science* 33:e5087 (2024). DOI: 10.1002/pro.5087; PMID: 39074255.
- **LA24** — Landerer C, Poehls J, Toth-Petroczy A. Fitness effects of phenotypic mutations at proteome-scale reveal optimality of translation machinery. *Molecular Biology and Evolution* 41 (2024). PMCID: PMC10939442. Interpretation here uses the project's corrected codon analysis, not its withdrawn detection-conditional estimate.
- **S2016** — Schmidt et al. E. coli protein abundance measurements. *Nature Biotechnology* 34:104 (2016). DOI: 10.1038/nbt.3418; PMID: 26641532; PMCID: PMC4888949.
- **T98** — Tomoyasu T, Ogura T, Tatsuta T, Bukau B. DnaK/J and heat-shock control. *Molecular Microbiology* 30:567 (1998). DOI: 10.1046/j.1365-2958.1998.01090.x.
- **G04** — Guisbert E, Herman C, Lu CZ, Gross CA. GroEL/S and σ32 control. *Genes & Development* 18:2812 (2004). DOI: 10.1101/gad.1219204. Abstract-level source in the feedback ledger.
- **L13** — Lim et al. σ32 membrane targeting and regulation. *PLoS Biology* 11:e1001735 (2013). DOI: 10.1371/journal.pbio.1001735.
- **Z14** — Zhang et al. Engineered σ32 regulation and recombinant-protein folding. *ACS Chemical Biology* 9:1945 (2014). DOI: 10.1021/cb5004477.
