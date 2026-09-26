# Independent evidence and falsification audit

Date: 2026-09-26

This is an independent second-wave audit of the proteostasis-only claims. It is
not a manuscript draft and it does not modify the first-wave redraw or audit.
The requested relative path `../proteostasis-P1` was absent in this checkout;
the same upstream files were found and read at
`/storage/kiran-stuff/proteostasis-P1/`. That path discrepancy is recorded
rather than silently normalized.

## Executive finding

The arithmetic reference value is wrong as quoted. With `N=300`,
`P_correct=0.7`, `p_misfold=0.3`, and `S=0.3`, the exact raw-error threshold is

```text
f_raw = [1 - P_correct^(1/N)] / [(1-S) p_misfold]
      = 0.0056581428505 per codon.
```

The log approximation is `0.0056615070466`. The quoted `0.00119` is instead
`-ln(0.7)/300`, which drops both the misfold probability and the synonymous
filter. The upstream code contains the correct denominator at
`arithmetic_stress_test.py:14-19,70-81`; the manuscript and summary preserve
the inconsistent value at `MANUSCRIPT.md:104-116` and
`arithmetic_summary.md:14-19`.

This is a parameter/identity failure, not evidence that the model is viable or
invalid. The number must be relabeled as a model-implied arithmetic threshold
after the error definition is fixed.

## 1. Measured amino-acid substitutions versus raw decoding errors

The upstream ODE mapping is
`J_bare = f_codon N_prot (1-S_avg) p_baseline / T_gen`
(`LITERATURE_ANCHORS.md:76-85`; also `two_pool_ode.py:261-263`). That mapping is
valid only if `f_codon` is a raw decoding-error probability before the
synonymous outcomes are removed.

There are two non-equivalent quantities:

```text
e_raw          = all decoding errors per codon
f_substitution = nonsynonymous amino-acid substitutions per codon
f_substitution ≈ e_raw (1-S)
J_bare         ∝ f_substitution p_misfold
```

Therefore:

* If the empirical input is raw decoding error, multiplying by `1-S=0.7` once
  is appropriate.
* If the empirical input is an empirically measured amino-acid substitution
  rate, multiplying it by `0.7` again is a double discount and must be removed.
  The measured value should enter as `f_substitution`.

Landerer et al. 2024 describe proteome-scale amino-acid misincorporation
observations and a fitted model over more than 100 mass-spectrometry datasets,
with 20–23% of proteins expected to contain at least one substitution; that is
not automatically a raw synonymous-inclusive decoding rate. Their primary
article is [PMC10939442](https://pmc.ncbi.nlm.nih.gov/articles/PMC10939442/);
the manuscript’s extraction/definition is also stated at
`MANUSCRIPT.md:285-286`.

The Stikeleather et al. 2026 primary record is especially explicit: it reports
the wild-type rate as **1.82e-3 per codon** (SE 5.92e-5) and gives the
substitution spectra (`PMC13335486`, lines 241-249). It also states that
cultures were grown overnight to stationary phase (lines 132-135), and warns
that pooled prior datasets combine strains and growth conditions (lines
306-314). That measured substitution rate cannot be put through a second
`(1-S)` filter. It also cannot be treated as an exponential-growth calibration
without changing the condition bundle.

The standalone executable checks are in
[bound_identity_checks.py](bound_identity_checks.py) and
[test_bound_identity_checks.py](test_bound_identity_checks.py).

## 2. Direct arithmetic reproduction

The upstream point formula is `arithmetic_stress_test.py:14-19` and the exact
implementation is `:70-75`. Independent evaluation gives:

| quantity | value per codon | interpretation |
|---|---:|---|
| exact raw-error threshold | `5.6581428505e-3` | includes one `(1-S)` and `p_misfold` |
| `-ln(P)/(N(1-S)p_misfold)` | `5.6615070466e-3` | large-N approximation |
| substitution threshold, if input is already substitution-level | `3.962...e-3` | no synonymous filter; `f_sub` is the input |
| quoted `-ln(P)/N` | `1.1889164798e-3` | neither a raw-error nor substitution threshold under the stated model |
| exact unfiltered `1-P^(1/N)` | `1.1882099986e-3` | close to the quoted number, but still not the stated model |

The arithmetic result at `arithmetic_summary.md:16` therefore cannot be
described as reproducing the 0.00119 value while also claiming the stated
`p_misfold` and `S` parameters. The proteome-integrated result at
`arithmetic_summary.md:32-40` is a different calculation and should not be
used to conceal this point-estimate identity error.

## 3. `Prot_tot=300 uM`: not total E. coli protein

The ODE sets `Prot_tot_uM=300` at `two_pool_ode.py:95-117`, and the literature
anchor calls that “total soluble cytoplasmic protein” at
`LITERATURE_ANCHORS.md:45-50`. The independent BioNumbers record [BNID
104726](https://bionumbers.hms.harvard.edu/bionumber.aspx?id=104726&s=n&v=18)
reports **4 mM** for total protein in exponentially growing *E. coli* on
glucose medium, with its calculation and assumptions shown on the record.
Thus 300 uM is approximately 13.3-fold below the cited total concentration.

The defensible classification, if 300 uM remains in the reduced model, is
“effective accessible misfolded/client pool,” not total cellular protein. That
interpretation must be stated and independently measured or calibrated.

This difference materially rescales the model’s algebraic feedback without
proving a viability outcome. At fixed `P`, `C_tot=50 uM`, `K_d=1 uM`,
`k_obs_max=0.01/s`, `k_agg=1000 M^-1 s^-1`, `k_deg=3e-4/s`:

| fixed `P` | `Prot_tot` | `M=P Prot_tot` | `C_free` | `Phi` |
|---:|---:|---:|---:|---:|
| 0.001 | 300 uM | 0.3 uM | 38.46 uM | 1.030 |
| 0.001 | 4,000 uM | 4.0 uM | 10.00 uM | 1.426 |
| 0.01 | 300 uM | 3.0 uM | 12.50 uM | 1.314 |
| 0.01 | 4,000 uM | 40.0 uM | 1.22 uM | 7.903 |
| 0.10 | 300 uM | 30 uM | 1.61 uM | 5.635 |
| 0.10 | 4,000 uM | 400 uM | 0.125 uM | 284.96 |
| 0.25 | 300 uM | 75 uM | 0.658 uM | 18.57 |
| 0.25 | 4,000 uM | 1,000 uM | 0.050 uM | 1290.10 |

These are scaling diagnostics only. They are not a recomputation of the
paper’s viability claim. The ODE’s `M`, `C_free`, `v_agg`, and `Phi` definitions
are at `two_pool_ode.py:8-16,127-150`.

## 4. `A_max=0.25` is not an established measured death threshold

The anchor file presents `A_max=0.15-0.35`, baseline 0.25, as a fraction of
aggregated proteome at which viability collapses (`LITERATURE_ANCHORS.md:104-116`).
The ODE then uses that number to restrict the saddle-node domain
(`two_pool_ode.py:27-36,209-258`). This makes `A_max` a gate on the headline
threshold, not a harmless descriptive parameter.

Primary-record audit:

* Ciryam et al. 2013 is [Cell Reports 5:781-790, PMID
  24183671](https://pubmed.ncbi.nlm.nih.gov/24183671/), as the primary record
  confirms. Its result concerns supersaturation and a metastable
  neurodegeneration-relevant subproteome; it does not establish a universal
  *E. coli* “20% aggregate mass causes death” threshold. The anchor’s citation
  is therefore materially overextended.
* The accessible Bednarska record is [“Protein aggregation in bacteria: the
  thin boundary between functionality and toxicity,” Microbiology 159,
  1795–1806, PMID 23894132](https://pubmed.ncbi.nlm.nih.gov/23894132/).
  It supports the qualitative statement that bacterial aggregation can reduce
  fitness/viability, but it is not the cited “Mol Cell 52:617” primary record
  and does not validate a universal 0.20–0.25 whole-proteome mass cutoff.

The exact citation bundle should therefore be labeled **uncertain / not
substantiated for quantitative use**. Until a matched *E. coli* experiment
measures aggregate fraction and viability under the same perturbation, use
`A_max` only as an explicit sensitivity parameter. Do not call it a measured
cell-death threshold, and do not report aggregation-death bounds as empirical
maximum tolerable rates.

## 5. Model-structure audit and falsification limits

The two-pool equations are explicitly stated at
`two_pool_ode.py:4-25`:

```text
dP/dt = J_bare Phi(P) - R(P) - drain(P)(1-A_sat)
dA/dt = drain(P)(1-A_sat) - k_clear A
```

Important assumptions and failures of formal closure are:

1. **Positive feedback is imposed.** `Phi(P)=1+v_agg/(v_fold+k_deg)` is a
   closure at `two_pool_ode.py:13-15,140-141`; it assumes aggregation-linked
   client production/recirculation increases inflow. It is not derived from a
   measured translation or folding flux. The manuscript’s reduced model makes
   the same status clear at `MANUSCRIPT.md:70-80`: `+chi x^2` is a minimal
   approximation and is not fitted.
2. **`P+A<=1` is not enforced.** `P` and `A` are described as fractions, but the
   ODE has no projection, complement pool, or invariant proof. The cap is only
   an operational search boundary for `A`; `A_qs` itself is uncapped at
   `two_pool_ode.py:152-159`, and the anchor explicitly acknowledges possible
   `A>1` at `LITERATURE_ANCHORS.md:116`. A formally bounded model must include
   an explicit conservation/occupancy state or prove the domain invariant.
3. **Bounded operational model versus formal unbounded model.** Calling the
   threshold “death” adds an external absorbing boundary at `A_max`; the
   uncapped ODE does not itself generate that biological event. Results must
   distinguish mathematical runaway, operational cap crossing, and measured
   loss of viability.
4. **No nascent-chain competition.** The model treats `C_tot` as available to
   the damaged pool, while the manuscript acknowledges that nascent-chain
   folding is not represented (`MANUSCRIPT.md:234-236`). That omits the main
   route by which translation flux can compete with rescue capacity.
5. **No adaptation.** There is no stress-induced chaperone synthesis,
   translational slowdown, degradation reprogramming, growth-rate feedback, or
   dilution adaptation in the state equations. The manuscript’s factorial is
   therefore a stock-model prediction, not a cell-level confirmation.
6. **No condition matching.** The Stikeleather measurement is stationary-phase
   (`PMC13335486`, lines 132-135), whereas the plan itself correctly warns that
   exponential-phase chaperone and dilution parameters cannot simply be held
   fixed (`PLAN_P1_REPAIR.md:68-90,224-228`). Error rate, chaperone pool, and
   dilution must be varied as a condition triple.

The ODE can still be useful as a falsifiable reduced model, but only for
predictions conditional on these assumptions. It does not provide experimental
confirmation. Farkas et al. support the qualitative buffering premise in yeast
[eLife/PMCID PMC5788500](https://pmc.ncbi.nlm.nih.gov/articles/PMC5788500/),
and McDonald et al. show that different mistranslations can stress distinct
proteostasis branches [NAR 2025/PMCID
PMC12082455](https://pmc.ncbi.nlm.nih.gov/articles/PMC12082455/). Those studies
support mechanism plausibility, not the numerical *E. coli* threshold.

## 6. What would be decisive: matched-condition experiment

The minimal decisive design is a matched-condition *E. coli* experiment with
two orthogonal perturbations and multiple levels, not a single error-rate
comparison.

### Factors and gradients

Use a translational-fidelity perturbation that changes the amino-acid error
spectrum (for example, calibrated mistranslating tRNA/ribosome alleles or a
validated fidelity perturbation), crossed with proteostasis capacity reduction
(for example, titratable DnaK/GroEL/ClpB capacity or a matched chemical/genetic
perturbation). The core factorial is:

| | control capacity | reduced capacity |
|---|---|---|
| control fidelity | 0, 1 | 2 |
| elevated misincorporation | 3 | 4 |

Run at least 4–6 graded levels of each factor, with the perturbations crossed
in a response-surface design. Use the same strain background, medium,
temperature, carbon source, growth phase, dilution rate, sampling density, and
induction history in all cells. Include an unperturbed wild-type and matched
vector/allele controls. Analyze exponential and stationary phase separately;
do not pool them.

### Required simultaneous measurements

1. **Error spectrum:** targeted/untargeted high-resolution proteomics reporting
   source→destination amino-acid substitutions and codon-resolved rates. Report
   raw decoding proxies and substitution rates separately, with detection
   limits and missed-substitution correction.
2. **Translation flux:** growth rate, ribosome occupancy/polysomes or nascent
   chain labeling, total protein synthesis, and dilution rate. This tests the
   `J_bare` mapping rather than inferring flux from OD alone.
3. **Free and bound chaperone:** quantitative DnaK/GroEL/ClpB abundance plus
   client-bound versus free fractions (native pull-down/crosslinking or an
   orthogonal calibrated occupancy assay). Measure ATP state if feasible.
4. **Aggregation compartment:** soluble/insoluble fractionation plus imaging
   and aggregate proteomics, distinguishing diffuse misfolded monomer,
   inclusion bodies, membrane-associated material, and polar aggregates. The
   response must report aggregate mass/fraction, not only a fluorescent proxy.
5. **Viability:** CFU, single-cell membrane integrity, division/growth arrest,
   and recovery after stress removal. Predefine whether “death” means loss of
   CFU, irreversible arrest, or another endpoint.
6. **Stress response:** heat-shock/proteostasis transcript and protein panels
   (including DnaK/GroEL/ClpB, sigma-factor response, degradation pathways),
   measured on the same time course. This tests the adaptation omission.

### Decisive tests

Fit the observed time courses jointly to (i) the stock two-pool model, (ii) a
bounded-conservation variant with `P+A<=1`, and (iii) an adaptive extension
with measured capacity changes. Pre-register the following falsifiers:

* no monotone relation between measured substitution flux and misfolded/aggregate
  flux after controlling for translation rate;
* no capacity-dependent interaction in the matched 2x2/gradient design;
* an apparent threshold that moves with aggregate compartment or stress
  response rather than with total `A`;
* systematic `P+A>1` or failure to predict the measured aggregate fraction;
* a fitted `Phi` that is unnecessary, has the wrong sign, or requires
  condition-specific free parameters without independent measurements;
* stationary and exponential conditions requiring incompatible “universal”
  `A_max` or `Prot_tot` values.

Farkas 2018 and McDonald 2025 make the design biologically motivated, but the
pre-existing experiments are not this matched quantitative test. The novelty
needed for a defensible claim is a mechanistic, conditional *E. coli*
relationship between measured error spectrum, translation flux, free/bound
capacity, aggregate compartment, stress adaptation, and a predeclared
viability threshold—if such a threshold exists.

## 7. Go/no-go for full ODE reimplementation

### No-go now if

* the project cannot define whether `f` is raw decoding error or measured
  amino-acid substitution rate;
* `Prot_tot=300 uM` remains labeled as total protein, or no effective-pool
  calibration is available;
* `A_max=0.25` remains a claimed measured death threshold without a matched
  primary calibration;
* the reimplementation still allows `P+A>1` while calling states fractions;
* the experiment cannot measure error spectrum, translation flux, chaperone
  occupancy, aggregate compartment, viability, and stress response under the
  same conditions.

### Go only if

1. the arithmetic identity and error semantics are corrected and tested;
2. the ODE is explicitly bounded or the unbounded formalism is renamed and
   treated as a mathematical sensitivity model;
3. `Prot_tot`, `A_max`, `Phi`, `k_agg`, `k_clear`, and the nascent-chain omission
   are assigned independent measurements or honest uncertainty classes;
4. the matched 2x2/gradient data show a reproducible, condition-specific
   capacity–mistranslation interaction and the stock model predicts it without
   post hoc tuning; and
5. the result is reported as a conditional mechanistic threshold, not as
   experimental confirmation or a universal E. coli death boundary.

**Current decision: NO-GO for a full ODE reimplementation as a claim-producing
exercise. GO for a small, independent calibration/identity harness and the
matched experiment above.**
