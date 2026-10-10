# Confronting the route architecture with two existing experiments

Produced by `confront.py` (200 draws, seed 20261010); numbers in `confront.json`.
Read with [RESULTS.md](RESULTS.md).

## 1. The imposed serial upstream requirement is falsified as a hard constraint

In the shipped architecture both obligate client classes can only be rescued through
`TF -> GroEL/ES` or `DnaK/DnaJ/GrpE -> GroEL/ES`. Setting TF and DnaK capacity to zero therefore
leaves no productive route, whatever GroEL can do.

Measured here: with TF and DnaK capacity at zero and GroEL/GroES capacity multiplied by 1, 2, 5, 10
and 100, essential productive flux is exactly zero in **200 of 200 draws at every multiplier**.
The model cannot express GroEL-mediated rescue of a TF/DnaK double loss at all.

That is contradicted by published experiment. GroEL/GroES overproduction restores growth of
*E. coli* Δtig ΔdnaK52 at 30 °C and suppresses aggregation of proteins ≥60 kDa that normally
require TF and DnaK (Vorderwülbecke et al. 2004, FEBS Lett 559:181–187,
https://doi.org/10.1016/S0014-5793(04)00052-3). Growth below 30 °C alone also abrogates the
synthetic lethality. The double deletion is itself viable under permissive conditions
(Genevaux et al. 2004, EMBO Rep 5:195–200, https://doi.org/10.1038/sj.embor.7400067), and
TF fragments containing the N domain complement it (Kramer et al. 2004, J Bacteriol
186:3777–3784, https://doi.org/10.1128/JB.186.12.3777-3784.2004). References retrieved from
PubMed. So TF/DnaK upstream service is a *conditional*, not obligatory, precondition for
GroEL-mediated folding of these clients.

The shipped `no_serial` control goes to the opposite extreme: with the obligate classes routed
directly to GroEL/ES at normal cost, the double knockout is fully served in **200 of 200 draws at
1× GroEL**, so it cannot express the synthetic lethality either.

### Recommended compatible alternative, and how far it gets

Add a *finite and expensive* direct `GroEL/ES` route for obligate clients, costing κ times the
normal GroEL service, instead of deleting the route. Sweeping κ (same draws, double knockout at
baseline versus 10× GroEL/GroES):

| κ | double KO still short of the essential quota | rescued by 10× GroEL/GroES |
|---:|---:|---:|
| 1 | 0/200 | 200/200 |
| 2 | 1/200 | 200/200 |
| 5 | 14/200 | 200/200 |
| 10 | 31/200 | 200/200 |
| 20 | 46/200 | 199/200 |
| 50 | 71/200 | 186/200 |

Overproduction rescue is recovered immediately, at any κ. Baseline lethality is not: even at κ = 50
only 71 of 200 draws fail at baseline. **No single κ reproduces both halves of the experiment in a
majority of draws.** The honest reading is that a static cost penalty is the wrong device: the
experimental lethality is explicitly temperature-conditional (lethal at 30 °C, viable at 23 °C,
rescued by GroEL overproduction at 30 °C). A model that reproduces it needs the direct route's
throughput, or the client influx, to depend on condition — not a structurally absent route. Until
that is added, the `KO_TF_DnaK` rows in RESULTS.md should be read as encoded assumptions, not
predictions, and the double knockout should not be described as lethal without naming a temperature.

## 2. All-client radial normalization is set by the smallest pool, not by folding

Primary demand is 0.9 × λ*, where λ* is the largest scale at which **every** class, including the
aggregate class, is served productively. The aggregate class's only productive route consumes ClpB,
whose measured pool is 0.0117 µM against DnaK's 10.88 µM (Schmidt BW25113 glucose exponential).

Measured here: dropping the aggregate class from the radial requirement raises λ* by a median
**5.74×** (q05–q95: 1.0–141×); it is more than 2× higher in **138 of 200 draws**. The aggregate
share of demand in those draws has median 0.010.

The consequence is visible throughout RESULTS.md: at the primary load the baseline folded fraction
is 1 in every model, matched GroEL and DnaK losses produce a folding drop of exactly 0 in all 200
draws, and the pair-interaction term I is 0 for all four pairs in all three models. Only ClpB loss
does anything. The ensemble is normalized into a regime where the folding machines are in large
excess, so **the framework's central supra-additivity prediction is not tested by this run at all** —
not rejected, untested. Survival does not require productive ClpB-mediated recovery of every
aggregate; aggregation and degradation are legitimate outcomes. Normalizing demand to a margin that
insists otherwise is the wrong reference point. Recommended: calibrate λ* on the folding classes
only, and treat aggregate handling as an outcome rather than a quota.
