# Prior art: what is already done, and what is left to test

All records below were retrieved from PubMed on 2026-10-10 (title, authors, journal, DOI checked).
Scope statements about the three models come from reading the main texts: FoldEco from the PMC
XML; Santra and FoldME from the PMC HTML, plus FoldME's SI methods text. The Santra SI PDF and the
FoldEco Extended Procedures were **not** read. Anything marked "not stated" may be in those
unread parts.

## The capacity-and-routing network is not new

| model | what it already does | what it does not do (in text read) |
|---|---|---|
| **FoldEco**, Powers, Powers & Gierasch 2012, Cell Rep 1:265–276, [10.1016/j.celrep.2012.02.011](https://doi.org/10.1016/j.celrep.2012.02.011), PMID 22509487 | Mass-action ODE model of the E. coli cytosol: synthesis, U/M/N/aggregate states, nucleated aggregation, TF, KJE, GroEL/ES, ClpB+KJE disaggregation, Lon, and ClpAP-type degradation of tagged native protein. Finite chaperone concentrations. Degradation competes with folding. Fits Siller 2010 luciferase translation-rate data. Parameters largely in vitro. | Growth dilution explicitly excluded. No temperature. No HtpG, ClpXP, HslUV or Ibp. No bistability analysis, although KJE binds aggregates 1:1 (sequestration exists in the mechanism but is not analysed as feedback). Co-translational folding not modeled. |
| **Santra, Farrell & Dill 2017**, PNAS 114:E2654, [10.1073/pnas.1620646114](https://doi.org/10.1073/pnas.1620646114), PMID 28292901 | Builds on FoldEco. Kinetic model with growth rate λ, synthesis set by growth, and client sorting by misfolded-state stability and kinetic accessibility. Finite chaperone pools; GroEL near 100% occupancy at 40-min doubling. Concludes chaperone levels are "just high enough". Needs more chaperone at both fast and very slow growth. | Fixed 37 °C. No HtpG. No bistability. "Tipping point" means aggregation onset in abundance, not a bifurcation. "Economical" is a model conclusion, not a measurement. In-vitro-fitted rates; no in-vivo kinetic calibration. |
| **FoldME**, Chen et al. 2017, PNAS 114:11548, [10.1073/pnas.1705524114](https://doi.org/10.1073/pnas.1705524114), PMC5664499 | Genome-scale ME model (optimization, flux balance) with temperature-dependent folding (24–46 °C), DnaK and GroEL folding fluxes bounded by chaperone turnover, and proteome reallocation from ribosomes to chaperones. Reproduces the growth–temperature curve and the sharp collapse above 42 °C. | **Assumes the cell optimizes growth.** No explicit aggregate state, degradation or disaggregation. No HtpG. No kinetic routing competition. Folding parameters are theoretical predictions. |

Taken together: a finite-capacity chaperone network, kinetic competition of folding with
degradation and aggregation, protein-class sorting, growth-coupled demand, occupancy near
saturation, and temperature-driven collapse are **all published**. This project must not present
any of them as its contribution.

## Experimental anchors (also prior art, used as constraints)

| source | fact used |
|---|---|
| Kerner et al. 2005 Cell 122:209, [10.1016/j.cell.2005.05.028](https://doi.org/10.1016/j.cell.2005.05.028), PMID 16051146 | ~250 GroEL interactors; ~85 classed obligate, 13 essential; obligate class reported to occupy >75% of GroEL capacity |
| Fujiwara et al. 2010 EMBO J 29:1552, [10.1038/emboj.2010.52](https://doi.org/10.1038/emboj.2010.52), PMC3212837 | only ~60% of Class III are obligate in vivo under GroE depletion; obligate set limited by cavity size |
| Calloni et al. 2012 Cell Rep 1:251, [10.1016/j.celrep.2011.12.007](https://doi.org/10.1016/j.celrep.2011.12.007), PMID 22832197 | DnaK ~700 interactors; TF loss reroutes clients to DnaK; on GroEL depletion clients accumulate on DnaK and are then degraded; ΔdnaK Δtig gives extensive aggregation |
| Vorderwülbecke et al. 2004 FEBS Lett 559:181, [10.1016/S0014-5793(04)00052-3](https://doi.org/10.1016/S0014-5793(04)00052-3), PMID 14960329 | Δtig ΔdnaK52 lethal at 30 °C, viable at 23 °C, rescued at 30 °C by GroEL/ES overproduction but not by other chaperones. **PubMed tags this record "Duplicate Publication"; the primary version should be identified before citation.** |
| Genevaux et al. 2004 EMBO Rep 5:195, [10.1038/sj.embor.7400067](https://doi.org/10.1038/sj.embor.7400067) | Δtig ΔdnaK viable under specific conditions, in a very narrow temperature range |
| Kramer et al. 2004 J Bacteriol 186:3777, [10.1128/JB.186.12.3777-3784.2004](https://doi.org/10.1128/JB.186.12.3777-3784.2004) | TF N-domain fragments complement Δtig ΔdnaK |
| Fauvet et al. 2021 Front Mol Biosci 8:653073, [10.3389/fmolb.2021.653073](https://doi.org/10.3389/fmolb.2021.653073), PMC8082187 | ΔKJG has higher metabolic/respiratory enzyme levels than ΔKJ; HtpG proposed to mediate degradation of aggregation-prone DnaK substrates, possibly via HslUV. Growth improvement of ΔdnaKJ by ΔhtpG at 37/39 °C is per the earlier evidence audit (Fig 1A), not re-read here. |
| Reingewertz et al. 2024 Protein Sci 33:e5087, [10.1002/pro.5087](https://doi.org/10.1002/pro.5087), PMID 39074255 | synonymous and non-synonymous substitutions in mouse DHFR reduce GroEL dependence in E. coli; the synonymous effect is attributed to translation rate **by computational analysis only, not measured** |
| Kanemori et al. 1999 J Bacteriol 181:3674, [10.1128/jb.181.12.3674-3680.1999](https://doi.org/10.1128/jb.181.12.3674-3680.1999) (retrieved in the 2026-10-10 evidence run) | protease double-deletion plating is temperature-conditional, non-monotone |
| Schmidt et al. 2016 Nat Biotechnol 34:104, [10.1038/nbt.3418](https://doi.org/10.1038/nbt.3418) | machine abundances; concentrations, not capacities |
| Babu et al. 2014 PLoS Genet 10:e1004120, [10.1371/journal.pgen.1004120](https://doi.org/10.1371/journal.pgen.1004120) | genome-wide GI screen at 32 °C; supplementary S2 unretrieved, reserved as held-out |

## What this work could newly test (claimed novelty is provisional; a targeted citation search for each item is still owed)

1. **Necessary versus sufficient, stated as theorems.** The occupancy identity ΣJτ ≤ M^tot is
   Little's law. It is exact for any steady state, so it is a model-free necessary condition. A
   counterexample shows feasibility is not sufficient (EQUATIONS A4–A5). Santra reports occupancy
   but does not separate a necessary bound from the model's realized state. **Likely presentation
   novelty, not new physics.**
2. **A disposal route that can lower yield, with an analytic sign rule.** None of the three models
   includes HtpG. The rule sign(∂Y/∂H) = sign[(1−q)(k_a^eff+μ) − q(k_f+p_G k_G G)] (A8) turns
   Fauvet's observation into a conditional, signed prediction: ΔhtpG should help less, or hurt,
   when disaggregation is lost (ClpB) and when the disposal protease is gone. **Testable.**
3. **Optimizing versus kinetic models, separated by one class of data.** Any model where the cell
   attains the best point of its feasible region (FoldME-type) cannot show improvement on deleting
   a route without an added cost term. A kinetic-routing model can. Fauvet-type suppression is
   therefore discriminating evidence. This must be stated carefully, because FoldME has a
   proteome-cost term that could absorb it.
4. **Runaway only with a named feedback.** FoldEco contains the mechanism (KJE:aggregate 1:1) but
   no bifurcation analysis. A saddle-node at L_crit = (μ+δ)a²/(4bk_a) appears only when the
   sequestration coefficient b > 0. The test is hysteresis under graded load: present with
   feedback, absent without.
5. **Codon timing separated from substitution frequency.** σ (co-translational partition) and ε
   (misincorporation frequency) are different parameters. Proteomic error data constrain only ε.
   FoldEco mimics translation rate by a ~2x k_f change. None of the three treats synonymous codons.

Not claimed as new: a capacity network, kinetic routing, client classes, growth coupling,
temperature dependence, or the existence of collapse.
