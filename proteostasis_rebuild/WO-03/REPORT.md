# WO-03 — chaperone allocation and kinetics alternatives

## gates (fixed before analysis)

G3.1 exact finite-pool equilibrium, residual < 1e-10, audit value C_f = 6.58872
at M_T = 50 · G3.2 driven cycle governed by K_M = (k_off+k_cat+mu)/k_on, shown
to differ from K_d · G3.3 nascent competition conserves chaperone and reduces
exactly at zero load · G3.4 separate benchmark / kinetics / competition and map
each to DnaK, GroEL, ClpB · G3.5 QSS reductions checked against full ODE.

## deliverables

`chaperone.py`, `run_wo03.py` → `wo03_results.json`, `test_wo03.py` (9 OK).

## results

**G3.1 equilibrium benchmark.** Exact finite-pool root over M_T = 0–300 uM:
residuals ≤ 7.7e-14. At C_T = 50, K_d = 1, M_T = 50 uM: C_f = 6.588723 (exact)
vs 0.980392 (legacy closure); legacy-shaped folding rate ratio 1.753816. Both
match the prior audit to 6 figures, recomputed independently. Negative control:
at M_T = 10 uM the legacy closure implies 45.5 uM of bound chaperone against
10 uM of client — it binds chaperone to clients that do not exist.

**G3.2 driven cycle.** Closed cycle C + U ⇌ B → C + U at C_T = M_T = 50 uM,
K_d = 1 uM, full ODE to steady state:

| k_cat/k_off | K_M (uM) | rel. err, quadratic with K_M | rel. err, with K_d |
|---:|---:|---:|---:|
| 0.001 | 1.001 | 1e-15 | 7e-5 |
| 0.1 | 1.1 | 0 | 0.007 |
| 1 | 2 | 5e-16 | 0.060 |
| 10 | 11 | 6e-16 | 0.382 |

The steady state of a driven cycle is a non-equilibrium steady state; its
occupancy is set by K_M, and K_d is recovered only when k_cat ≪ k_off.

Four-state DnaK-like cycle (ATP state fast/weak K_dT = 10 uM, ADP state
slow/tight K_dD = 1 uM; hydrolysis stimulated on the complex; illustrative
rates): driven K_eff = 0.21 uM — **tighter than either nucleotide state**. With
hydrolysis tuned to satisfy detailed balance (Kolmogorov criterion), K_eff =
5.5 uM, between the two states as equilibrium requires. So the direction of the
error in using a peptide K_d is not even fixed: a driven Hsp70 can bind tighter
than any equilibrium K_d (this is the "ultra-affinity" argument, attributed from
memory to De Los Rios & Barducci 2014 eLife — citation UNVERIFIED until WO-06),
while a fast-completing cycle binds weaker (K_M > K_d).

**G3.3 nascent competition.** Extended model with nascent clients X and
complexes BX. Symbolic: P_T and C_T identities hold; with nu_c = 0 and
X = BX = 0 the first five equations are term-for-term the WO-02 equations and
dX/dt = dBX/dt = 0 (so the reduction is exact, not just numerically close;
numerical trajectories agree to 8.7e-10). The legacy's theta becomes an output,
theta = BX/C_T; it increases monotonically with nascent load (0 → 0.088 for
nu_c 0 → 0.4 at illustrative nascent rates). These theta values say nothing
about E. coli; they show the quantity is computable once nascent-chain
engagement rates are measured.

**G3.5 QSS.** Finite-pool tQSSA (B eliminated via K_M, client total W = U+B as
slow variable) against the full WO-02 ODE over 10 doublings: max rel. err
6.9e-5 in U, 4.2e-3 in A after the first 0.1 doubling. At steady state the QSS
relation for B is exact, not approximate.

## G3.4 — what each layer can and cannot stand for

| layer | can represent | cannot represent |
|---|---|---|
| equilibrium K_d benchmark | an in-vitro chaperone–peptide titration (e.g. the Pierpaoli-type K_d the legacy cites) | any ATP-consuming chaperone in vivo; capacity under flux |
| simple cycle, K_M | a generic ATP-driven foldase as an enzyme: DnaK/DnaJ/GrpE or GroEL/ES treated as one effective cycle with partition phi (productive vs released unfolded) | nucleotide-state asymmetry; GroEL size limit and client specificity; co-chaperone limitation |
| four-state cycle | DnaK-like nucleotide cycle with hydrolysis-driven locking; shows affinity is a kinetic, not thermodynamic, property | GrpE/DnaJ concentrations explicitly; multiple binding sites per client |
| competition (X, BX) | nascent-chain engagement of the same pool (DnaK; trigger factor would be a separate ribosome-bound pool) | trigger factor's ribosome tethering; cotranslational vs post-translational timing |
| ClpB | **not represented by any of the above.** ClpB acts on aggregates, with DnaK, so it belongs in k_dis. In WO-02 k_dis is a constant; biologically it should fall when free DnaK is scarce. That dependence is a candidate positive-feedback route (aggregates → DnaK sequestration → less disaggregation) and is tested in WO-04 | — |

## adversarial self-review

- *Are the cycle results just algebra dressed as simulation?* The K_M result is
  algebraically forced; the ODE confirms the code, not the biology. What is not
  forced is its consequence: the legacy's use of an equilibrium K_d as the
  saturation constant of an ATP-driven rescue arm has no justification, and the
  error can go either way (K_M > K_d or ultra-affinity K_eff < K_d).
- *Could the ultra-affinity demonstration be an artefact of my chosen rates?*
  Its existence is rate-dependent, yes — that is why the detailed-balance
  control is included: the same code, with drive removed, obeys the
  equilibrium bound. The magnitude (0.21 vs 1 uM) is illustrative only.
- *Is the competition model adequate?* It adds one client class sharing one
  pool. Real cells have trigger factor, DnaK, and GroEL with different client
  sets; a single X is the minimal structure that turns theta from an input into
  an output. It is not a quantitative model of allocation.
- *Anything weakened?* G3.3 said "reduces exactly". The first run showed only
  numerical agreement (8.7e-10). I added the symbolic identity so the gate is
  met as written rather than approximately.

VERDICT: PASS — G3.1–G3.5 all met.
