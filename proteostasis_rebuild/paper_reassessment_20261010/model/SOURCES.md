# Evidence and assumption boundary

Primary sources checked on 2026-10-10; the benchmark reads only the local WO-07 TSV and does not download or alter legacy inputs.

| Evidence | Use here | Not justified by that evidence |
|---|---|---|
| [Schmidt et al., Nature Biotechnology 2016](https://pubmed.ncbi.nlm.nih.gov/26641532/), DOI 10.1038/nbt.3418; local `WO-07/bundles.tsv`, EXPONENTIAL/MATCHED | Seven separate glucose-condition pools, in functional-unit concentrations | Throughput, client-specific service times, ATP demand, or headroom |
| [Kerner et al., Cell 2005](https://pubmed.ncbi.nlm.nih.gov/16051146/), DOI 10.1016/j.cell.2005.05.028 | About 85 obligate clients among about 250 associated clients; 13 essential; motivates a restricted obligate subset | 85/250 is not the proteome-wide demand share, and 13/85 is not its essential flux fraction |
| [Calloni et al., Cell Reports 2012](https://pubmed.ncbi.nlm.nih.gov/22832197/), DOI 10.1016/j.celrep.2011.12.007 | DnaK cooperation with upstream TF and downstream GroEL; combined TF/DnaK loss motivates a cooperation hypothesis | Exact routes of every client, demand fractions, rate coefficients, or a quantitative out-of-sample prediction |
| [Babu et al., PLOS Genetics 2014](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1004120) | 3,816 high-confidence interactions involving chaperones/proteases; not used as training data | GI scores as machine-work coefficients, or LP flux as fitness |

Kerner's reported >75% occupation of GroEL capacity by obligate clients is **not** a 75% share of all incoming client flux. Converting it to a flux constraint needs client-specific service times and the experimental condition. The benchmark does not impose it, so it is architecture-informed rather than a fitted representation of the complete Kerner dataset. A future fit must account for this occupancy observation rather than silently replacing it with count fractions.

Calloni reports roughly 700 DnaK interactors, including aggregation-prone clients, and a combined-loss phenotype involving downstream GroEL function and ribosome biogenesis. The model's TF-or-KJE upstream module followed by GroEL is a coarse hypothesis, not a literal complete measured pathway map. In particular, requiring upstream assistance for every modeled obligate client is an extra assumption, tested by `no_serial`. The overlap class's exclusive TF/KJE alternatives are also an assumption, tested by adding a finite GroEL alternative. No gene-level client list or complete perturbation dataset is available in this implementation.

The ClpB/KJE recovery route is an illustrative architecture choice, not a quantitative constraint extracted from these datasets. Other recovery mechanisms could change its sensitivity. Generic protease and storage resources have explicitly synthetic capacities; no Lon, ClpXP or Ibp measurement is claimed. HtpG is omitted to keep the scaffold minimal, not because its abundance or routes are zero. Generalization to a full proteostasis network needs additional qualified routes and capacities.

The older `network_failure_test` results are not evidence used by this model. Its arbitrary route probabilities and cut frequencies are neither imported nor re-estimated.
