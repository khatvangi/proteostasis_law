# Measurement Map

| Symbol or object | Meaning | Candidate measurement | Key uncertainty/control |
|---|---|---|---|
| `E_p(u,e)` | protein-specific translation flux | pulse labeling, ribosome profiling calibrated to synthesis | elongation bias, absolute scale, turnover |
| `Pr(r|c,u,e)` | residue delivered for codon and condition | targeted proteomics, reporter libraries, decoding assays | detection bias, rare-event calibration |
| `kappa_{p,i}` | site criticality modifier | deep mutational scanning, conservation-informed assays | circularity, context transfer |
| `d(r,p,i,e)` | substitution-specific damage contribution | variant folding, solubility, function, degradation assays | scale definition, epistasis, double counting |
| `J_fold` | non-error folding failure influx | pulse-chase folding yield and insolubility | distinguish delay from failure |
| `N_f`, `C_N` | free/occupied nascent folding load | ribosome-associated chaperone occupancy | binding stoichiometry and dwell time |
| `C_tot`, `C_f` | total/free chaperone pool | quantitative proteomics plus occupancy sensors | compartmentation, paralogs |
| `D_tot`, `D_f` | total/free degradation pool | abundance plus activity/queue assays | substrate competition, ATP dependence |
| `U` | soluble damaged/misfolded burden | conformation sensors, limited proteolysis, solubility fraction | proxy specificity |
| `A` | aggregate burden | microscopy, sedimentation, aggregate proteomics | reversible condensates vs aggregates |
| rate constants | rescue, aggregation, clearance kinetics | time-resolved perturbation and pulse-chase fits | identifiability, model mismatch |
| `H_j` | damage thresholds | growth, viability, function, lineage persistence | endpoint choice and adaptation |

All measurements must retain condition labels, biological replication, uncertainty, and units. Orthogonal proxies should be used where no direct measurement exists.

