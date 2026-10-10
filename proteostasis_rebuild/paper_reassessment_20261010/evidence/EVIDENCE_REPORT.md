# Evidence reassessment — 10 October 2026, Boron

**Experimental constraints exist; a validated client-demand/throughput network does not yet exist here.** The earlier audit's blanket claim that compatibility is unmeasured is too broad. Binding, perturbation-dependent folding and growth interactions constrain different parts of the problem. None licenses converting client counts into demand flux.

This is a source audit and limited empirical reanalysis, **not a completed replication of the three supplements**. All 14 recorded shell downloads failed DNS resolution (curl exit 6). The separate web tool recovered article text/indexed tables, but no original Kerner/Calloni client workbook or Babu interaction workbook. No missing experimental rows were generated. `retrieval_attempts.tsv` records exact attempted URLs and failures; `source_manifest.tsv` distinguishes confirmed identifiers from candidate filenames. Checksums of tool-response snapshots authenticate the saved responses, not unavailable publisher files.

## Existing project evidence

The project contains a checksum-verified Schmidt experimental workbook, bibliographic/full-text records, derived parameter tables, and synthetic network outputs. `local_data_inventory.tsv` classifies 77 data/record files; the baseline covers 137 original files. A filename scan of the wider `proteostasis_law` tree found no named Kerner, Calloni, Babu or Keio datasets. The existing Calloni PubMed XML is an abstract, not its client table.

`machine_abundance.tsv` re-extracts all 15 machine abundances using the existing checked loader. Its source is Schmidt Tables S6/S23, BW25113 glucose exponential growth; S23 volume is calculated. Functional-unit stoichiometries are model conventions, particularly for IbpA/B. Concentrations are **not capacities**. For example, DnaK is 10.8803 μM and ClpB 0.0116582 μM under those conventions. No matched in-vivo turnover conversion was recovered.

## Supported facts and limits

- **Kerner 2005:** approximately 250 interactors, approximately 85 classified obligate clients, 13 essential-client annotations and >75% capacity allocation are the original reported claims. The last quantity is a capacity-share/occupancy interpretation, not molecules per time. Supplemental-methods indexing identifies MC4100 derivatives; assay-level growth conditions and the original Table S3 rows remain unverified. [Original paper](https://doi.org/10.1016/j.cell.2005.05.028)
- **Calloni 2012:** MC4100 KHis cells, M63 SILAC, 37°C exponential growth yielded 674 interactors (672 excluding DnaJ/GrpE). S6/Fig4 document TF-loss rerouting; S7/Fig5 document GroE-depletion accumulation on DnaK and protein loss. S8/Fig6 compare aggregation at 30°C: 474 proteins for ΔdnaKJ versus 1,087 for ΔdnaKJ Δtig. These are protein-identity counts, not fluxes, biomass fractions or clean dnaK-only deletions. Pulse-chase binding kinetics also exist. No row-level client intersection was computed; Fig5's published overlap cannot substitute for an authenticated table join. [Deposited manuscript](https://pure.mpg.de/rest/items/item_1602419_3/component/file_1836221/content)
- **Babu 2014:** 3,816 high-confidence GIs involve at least one of 22 stability factors, not just pairs among the 15 machines. Screening used Hfr-C donors, Keio deletion/hypomorphic recipients, LB with Kan/Cm, 32°C, 36 h. Negative S-scores indicate aggravating growth interactions; they do not necessarily mean lethality. Positive interactions can reflect shared pathway membership. Table S10 overlap concerns GI neighborhoods, not physical clients. [Original article, Figs1/6](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004120)
- **Dependence correction:** Fujiwara's conditional MGM100 depletion study found 49 Class III proteins dependent and 34 independent, with Ypt1 untested (Fig2C; S3/S4). Thus the original ~85 estimate cannot be treated as 85 individually established obligate dependencies. [Experimental follow-up](https://pmc.ncbi.nlm.nih.gov/articles/PMC3212837/)

`evidence.tsv` separates physical interactome, inferred versus tested dependence, deletion/growth GI, abundance, occupancy and kinetic association. Co-purification alone does not establish direct binary contact or an independently sufficient rescue route. Degradation and aggregation are alternative outcomes; reduced total abundance alone does not measure degradation flux. There is no recovered absolute client-demand-flux row.

## Deletions and a concrete comparison

`machine_deletion_essentiality.tsv` covers all 15 genes. GroL/GroS have experimental inactivation/complementation evidence; GrpE has deletion-complementation evidence with suppressor exceptions. Nine other machines have condition-specific viability evidence, including a ClpB null whose exact boundaries were not verified. **Clean single-deletion conclusions remain unresolved for dnaJ, clpX and hslU**: accessible experiments delete multiple genes. This is not a 15-gene, condition-matched Keio validation. [GroE experiment](https://pmc.ncbi.nlm.nih.gov/articles/PMC209756/?page=2), [GrpE experiment](https://pubmed.ncbi.nlm.nih.gov/2651417/)

`model_vs_experiment.tsv` places these calls beside the existing `bio_measured_W2` outputs. The model generates an HtpG structural zero in 50.3% of networks; viable ΔhtpG W3110 cells constrain universal essentiality claims, but do not provide a matched accuracy estimate. Moreover, removing htpG improves ΔdnaKJ growth at 37/39°C (Fauvet Fig1A), competing with a model in which deleting capacity can only worsen growth. [Experimental study](https://www.frontiersin.org/journals/molecular-biosciences/articles/10.3389/fmolb.2021.653073/full)

`published_plating_table.tsv` transcribes 24 cells (23 numbers, one ND) from Kanemori Table2, explicitly as indexed published values, not raw colonies. At 30°C the relative plating efficiencies for WT, ΔhslVU, Δ(clpPX-lon) and their combination are 1.5, 1.3, 0.22 and 3.3×10⁻⁵. The calculated contrast `AB×WT/(A×B)` is 1.73×10⁻⁴. Because every strain is normalized to its own 42°C value, this measures the change in multiplicative interaction relative to 42°C, not absolute epistasis. Multi-gene deletions, rounded values and absent uncertainty preclude a p-value. SulA accumulation and suppressors provide a specific competing mechanism to generic folding overload. [Primary Table2 and experiment](https://pmc.ncbi.nlm.nih.gov/articles/PMC93843/)

## Three falsifiable predictions using existing experimental data

1. **Protease redundancy is conditional:** the combined deletion has lower relative plating efficiency than either deletion at 30°C and 45°C. Reject if either comparison reverses. The transcribed table supports both; this is a retrospective arithmetic check, not a blind prediction. It also rules out treating temperature as a simple monotone demand multiplier across these conditions.
2. **DnaK binding enrichment predicts dependence, imperfectly:** among proteins measured in both Calloni S1 and S3, greater baseline DnaK enrichment predicts >5% aggregation after ΔdnaKJ at 30°C, controlling for abundance and size. Prespecify a positive enrichment coefficient; reject if its confidence interval excludes positive values, otherwise report uncertainty. Fit on measured proteins only; use intact-protein IDs, exclude cofactors and do not turn missing detections into nonclients. **Executable after supplement retrieval; not tested here.** This is within-study validation, not independent replication.
3. **Proteomic dependence predicts aggravating GI:** freeze a gene-level dependence ranking from Calloni's perturbation tables, then test whether higher-ranked measured clients have more negative Babu S-scores with the corresponding chaperone perturbation than abundance-matched measured controls. Prespecify a negative association and a query-stratified permutation test. No negative association would challenge the proposed mapping from biochemical dependence to growth. **Babu can supply independent experimental outcomes, but no test was run.**

## Babu extraction and held-out status

[Table S2, supplement s023](https://journals.plos.org/plosgenetics/article/file?id=10.1371/journal.pgen.1004120.s023&type=supplementary) is explicitly XLS and should be structurally extractable; **its actual schema/rows were not accessed**, so successful extraction is unverified. Retrieve S1/S16 for alleles and Protocol S2/S16 for scoring; preserve donor/recipient direction and distinguish essential-gene hypomorphs from deletions. Of the model's 15 machines, 12 occur in the article's 22-factor list; GroL, GroS and GrpE are absent from that list, not proven unmeasured in every recipient position.

Freeze the predictor and exclusions before opening S2. Pair-level outcomes are uninspected, although the article's aggregate findings were read, so this is a **prospective held-out analysis plan, not pristine literature blindness**. S2 contains significant hits only: non-listed pairs must not be coded as neutral. Testing hit enrichment requires the tested-pair universe/full scores; without it, restrict any analysis to sign/rank among listed hits and acknowledge selection bias. Do not train on S9/S10 then call S2 independent. No graded-load test is possible from this one screening condition alone.

`SOURCE_INTEGRITY.md` reports the checks actually run. All outputs stay in this evidence folder; old runs and WORK_ORDERS.md remain unmodified. Ledger entries are staged in `LEDGER_PENDING.md`, because the global ledger lies outside the permitted write roots and this task restricts writes to this folder.

## Evidentiary levels (added 2026-10-10, this run)

`evidence_levels.tsv`, built by `scripts/levels.py`, assigns every row exactly one level. The level is
derived from `analysis_status`, not set by hand.

| Level | Meaning | Rows |
|---|---|---:|
| L1 direct data | recomputed in this project from a local original file whose hash we hold | 22 |
| L2 literature | a published value read from the article and not recomputed here | 28 |
| L3 assumed | a modelling convention or illustrative choice with no measurement behind it | 12 |

What this says plainly:

- **The only L1 evidence is three things.** The 15 Schmidt machine abundances (BW25113, glucose,
  exponential), the Landerer/Stikeleather codon error reanalysis (`codon_bridge/`), and two structural
  facts about the model itself (`model/confront.py`). Nothing else here was recomputed from an original
  file.
- **No L1 row measures a flux.** Abundances are concentrations; the codon tables are misincorporation
  counts with no folding, aggregation or growth readout. There is still no client-resolved demand flux
  anywhere in this project, so the model's demand vector remains L3 throughout.
- **Every quantitative route and share in the model is L3.** Twelve assumed inputs are enumerated by
  name, including the two the reassessment found to be wrong or load-bearing: the serial TF/DnaK
  requirement for obligate clients (contradicted by GroEL/ES overproduction rescue of Δtig ΔdnaK, see
  `model/CONFRONTATION.md`) and the all-client radial demand normalization.
- **Deletion essentiality is L2 and condition-bound.** Clean single-deletion calls are still unresolved
  for *dnaJ*, *clpX* and *hslU*; every accessible experiment deletes more than one gene.
- Nothing was invented to fill a gap. The Kerner, Calloni and Babu supplements were not retrieved, and
  no row stands in for them.
