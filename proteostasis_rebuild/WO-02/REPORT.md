# WO-02 — source–sink conservative equations

## gates (fixed before analysis)

G2.1 symbolic conservation of P_T and C_T · G2.2 dimension check of every RHS
term · G2.3 symbolic nonnegative-orthant invariance · G2.4 ≥200 random
integrations conserve totals (rel. < 1e-6), no state < −1e-9 · G2.5 legacy Phi
donor test and phantom-mass quantification.

## deliverables

`model.py` (single sympy definition, lambdified for numerics), `EQUATIONS.md`,
`proofs.py` → `proofs.json`, `test_wo02.py`.

## results

| gate | result |
|---|---|
| G2.1 | `sympy.simplify(dP_T/dt − (s_P − k_d U − k_dA A − mu P_T)) = 0`; same for C_T |
| G2.2 | 0 of all RHS terms fail the WO-01 checker; every term is uM/s |
| G2.3 | all five faces pass. eps, phi ∈ [0,1] handled by writing eps = a/(a+b), phi = c/(c+d); the cleared face polynomials have only nonnegative coefficients |
| G2.4 | 200 random parameter sets (log-uniform over ≥2 decades each), LSODA rtol 1e-10: 0 failures, max rel. err P_T 3.0e-13, C_T 9.5e-10, most negative state 0.0 |
| G2.5 | summing the legacy pool equations and removing the named source and sinks leaves exactly `J_bare(Phi−1)` (symbolic identity). No pool is debited |

**Size of the phantom inflow in the legacy.** At the legacy's own operating point
(mu = 6.33e-4, baseline): Phi = 1.033, so 3.2% of inflow is created from nothing.
At the legacy threshold P_dagger = 0.0274: Phi = 1.94, so **48% of inflow at the
threshold has no source**. The location of the legacy threshold therefore rests
substantially on non-conserved mass. (Its mechanism label at that point is
`aggregation_death`, the A_max gate; WO-04 separates the two.)

**What the conservative form excludes by construction.** Because every source is
a constant (`eps s_P`, `(1−eps) s_P`, `s_C`), a state-dependent amplification of
inflow cannot be written in this model. Any positive feedback must come from a
named mechanism acting on removal (e.g. chaperone sequestration, growth
coupling). Whether one exists is WO-03/WO-04's question, not an assumption here.

## tests

`python -m unittest test_wo02 -v` → 9 OK (includes two negative controls: a
Phi-like donorless source breaks the P_T identity; a sink not vanishing at U=0
is caught by the face test). `python proofs.py` → exit 0.

## adversarial self-review

- *Hard-coded pass.* My first G2.5 returned `donor_pool_found: False` as a
  literal. That is an assertion, not a test. Replaced with a symbolic residual
  on the WO-01 transcription of the legacy equations; it is now computed and
  checked against `J_bare(Phi−1)`.
- *Is G2.4 circular?* Partly. The augmented Q variable integrates the claimed
  total balance written independently of the species equations, so a coding
  error in a species term would show. It cannot catch a mistake shared by both
  (e.g. a wrong conservation law); G2.1 covers that symbolically, and the C_T
  check is against a closed form.
- *Does conservation make the model right?* No. It makes it admissible. The
  choice of species, the `k_a U²` aggregation law, and treating the chaperone
  pool as one effective species are modelling decisions with no validation here.
  Dropping the legacy `A_sat` saturation is a decision too: the legacy gave no
  mechanism for it, and keeping it would require one.
- *Is `s_P = mu P_T` right?* It assumes balanced exponential growth with
  negligible turnover of native protein. It is wrong for stationary phase,
  where mu → 0 and synthesis is not proportional to dilution. WO-07 must use a
  separate stationary bundle rather than setting mu small in this form.
- *Scenario values.* `SCENARIO` in model.py is illustrative and every entry is
  UNVERIFIED; P_T = 3000 uM is an order-of-magnitude placeholder. No WO-02 result
  depends on it (all gates are symbolic or use random parameters).

VERDICT: PASS — G2.1–G2.5 all met.
