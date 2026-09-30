# WO-07 — matched-condition bundles (2026-09-30)

## gates (fixed before analysis, WORK_ORDERS.md lines 143-153; not edited)

- G7.1 separate exponential and stationary bundles; every entry carries a
  condition tag and a match flag (MATCHED / MISMATCHED / UNMEASURED).
- G7.2 no bundle mixes phases silently: an automated check fails if a bundle
  draws a value tagged with the other phase without a MISMATCHED flag.
- G7.3 the effect of each mismatch on a model output is quantified as a range,
  not hidden in a point value.

## deliverables

| file | what it is |
|---|---|
| `bundles.tsv` | 80 rows, 20 columns, one row per (bundle, parameter, candidate); 27 parameters per bundle |
| `build_bundles.py` | writes the TSV. Every value is read from a saved file. The only literals are the in-vitro DnaK constants and the Stikeleather rate, and the checker confirms each literal appears in its WO-06 audit row |
| `derived.py` → `derived.json` | the copy-weighted length distribution and translated codon usage (Schmidt 2016 copies × K-12 CDS), plus the eTEL aggregates re-weighted |
| `mismatch_effects.py` → `mismatch_effects.json` | G7.3: one effect id per mismatch with output ranges; reference outputs with honest kinds; counterfactuals for prohibited borrowings |
| `check_bundles.py` | validator, exit 1 on any error. It also fails if the TSV differs from a fresh `build()` |
| `test_wo07.py` | 38 tests: 23 negative controls, 7 of which reproduce bypasses found by the independent reviewer |

- Nothing outside `WO-07/` was written, and no network was used.
- Inputs are files saved by WO-05/WO-06, plus two read-only legacy data files: the K-12 CDS FASTA (sha256 `9bd477a9…`, recorded in `derived.json`) and the UniProt length table (hashed in WO-00).
- A rerun of the pipeline is byte-identical (sha256 of all outputs compared).

## method

### anchors: what each bundle is matched TO

| bundle | anchor | why |
|---|---|---|
| EXPONENTIAL | *E. coli* BW25113, M9 glucose, 37 °C, balanced exponential growth | Schmidt 2016 glucose is the only condition in that source whose total protein is an independent measurement |
| STATIONARY | *E. coli* Xac (the paper's wild type), LB (Miller), 37 °C, overnight to stationary | the Stikeleather 2026 error input defines this bundle |

### two error quantities, kept apart

This split was added after review.

- **`substitution_freq_standing_f`** is what mass spectrometry on a harvested culture measures: substitutions per codon in the **standing** proteome. Both Landerer and Stikeleather measure this.
- **`synthesis_error_rate_per_codon`** is what a burden flux needs. It is **UNMEASURED in both bundles**.
- The two coincide only in balanced growth with no selective degradation of erroneous chains.
- In a stationary harvest most chains were made before growth stopped. So the Stikeleather value is MATCHED as a standing frequency and is *not* a stationary per-synthesis rate. Outputs that need a synthesis rate (O2, O2s, O3) list that substitution among their held assumptions.

### flags

- **MATCHED**: the value is *E. coli*, its source-recorded phase is the bundle's, and its **source-recorded** strain, medium and temperature equal the anchor. WO-06 must record it as VERIFIED with org_cond_match MATCHED.
- **MISMATCHED**: any other row that carries a number. The row lists its axes (PHASE_OPPOSITE, PHASE_MIXED, IN_VITRO, ORGANISM, STRAIN, MEDIUM, TEMPERATURE, TIME_IN_PHASE, GROWTH_RATE, VOLUME_ASSUMED, QUANTITY, SELECTION, NORMALIZATION_FROM_EXPONENTIAL) and effect ids.
- **UNMEASURED**: the value is `NA`, with no source.
- Every (bundle, parameter) pair has **exactly one** MATCHED or UNMEASURED row. MISMATCHED candidates sit next to an UNMEASURED row and never replace it.

### how G7.2 is enforced

The checker never trusts what a row says about itself. From the row's `source_key` it re-derives:

- the phase: the WO-06 audit phase field (a whole-word parser; hyphenated or two-phase text is AMBIGUOUS), Schmidt's per-condition phase, or the Schmidt condition behind a `derived.json` entry;
- for MATCHED rows, the organism, strain, medium, temperature and WO-06 status.

It then enforces these rules:

- An opposite-phase value must be MISMATCHED with PHASE_OPPOSITE. MIXED, IN_VITRO and AMBIGUOUS values must be MISMATCHED. A sourceless row must be UNMEASURED.
- Stikeleather (PMID 42406629) is rejected anywhere in EXPONENTIAL, whether found by source string, audit id, or the audit row's record.
- A pool row names exactly one machine. Any `+`, `sum`, `combined` or `C_tot` is rejected.
- Each MISMATCHED row's effect must vary *that* parameter, in *that* bundle.

## results

### the two bundles

| | EXPONENTIAL | STATIONARY |
|---|---:|---:|
| parameters accounted for | 27 | 27 |
| MATCHED | **13** | **1** |
| MISMATCHED candidate rows | 12 | 16 |
| UNMEASURED | 13 | 25 |

**EXPONENTIAL, MATCHED (13)**

- growth rate 0.58 h⁻¹
- synthesis rate s_P = µ·P_T = 0.478 µM chains s⁻¹ (a balanced-growth closure)
- P_T = 2965 µM
- eight separate machine pools
- copy-weighted length (mean 263.4 codons)
- translated codon usage

Every one is an abundance, a growth rate, or derived from them. **No kinetic parameter is matched.**

**STATIONARY, MATCHED (1)**

- the Stikeleather standing substitution frequency, 1.82 × 10⁻³ per codon (SE 5.92 × 10⁻⁵)

**UNMEASURED in both bundles**

- the synthesis error rate
- GroEL cycle rate and ClpB disaggregation rate
- k_d, k_a, p_misfold, φ, k_mis, k_dA, A_max
- the in-vivo DnaK cycle rate and K_M

**UNMEASURED in EXPONENTIAL only**

- the standing substitution frequency. None of the 80 Landerer datasets carries a phase annotation in the saved records, so all four eTEL aggregates are MISMATCHED (PHASE_MIXED) candidates.

**UNMEASURED in STATIONARY only**

- growth rate, s_P, P_T, every pool, and the length weighting of the anchor culture

### findings made while building the bundles

1. **Hidden cross-phase borrowing inside the stationary source.**
   - Schmidt 2016 measured total protein mass per cell for glucose only. It scaled every other condition "assuming that the volumetric protein concentration is condition independent".
   - So every absolute stationary concentration (P_T 3527 µM, DnaK 35.5 µM, …) carries the **glucose-exponential** volumetric protein concentration. The reviewer checked this independently: the calculated stationary volume cancels out.
   - These rows are MISMATCHED with NORMALIZATION_FROM_EXPONENTIAL.
   - Ratios to P_T cancel the normalization exactly: in `S_NORM`, O4 stays invariant while O2s scales linearly.
   - Schmidt's stationary cultures are also a different strain and medium (glucose M9 then 1–3 days starved) from the anchor (Xac, LB, overnight).
2. **The Stikeleather value is a standing-proteome frequency.** This was a reviewer ERROR, now fixed by the two-quantity split above.
3. **Two WO-06 quantity mismatches could be fixed from files on disk.**
   - Weighting lengths by glucose copies (2327 proteins, 99.8% of quantified copies) gives a mean of **263.4 codons**, against 307.6 genome-unweighted (legacy N = 300).
   - Translated codon weights shift the eTEL aggregates by 1.01–1.24× and O1 by 1.01–1.22×. As a control, the genomic-weight aggregates reproduce WO-05 to machine precision.
4. **Stationary net growth is negative** (−0.01 ± 0.003 h⁻¹). There is therefore **no dilution sink**, so O5 = 0 over the whole span. My first version took |µ| as a dilution rate; the reviewer caught the sign error, and it is fixed.

### chaperone pools, machine by machine

These are Schmidt 2016 BW25113 values per whole-cell volume, in functional units, never summed.

| machine | exp. glucose (µM) | exp. 20-condition span (µM) | stationary 1 d (µM) | stat./glucose | chains per unit, exp. (DERIVED) | chains per unit, stat. (SENSITIVITY) | cycle/turnover rate in saved records |
|---|---:|---:|---:|---:|---:|---:|---|
| DnaK (monomer) | 10.88 | 8.29–25.1 | 35.45 | 3.3× | 273 | 99 | in vitro only: T→R 0.04 s⁻¹, R→T 1.0 s⁻¹ |
| GroEL14 | 0.950 | 0.84–2.60 | 1.970 | 2.1× | 3122 | 1790 | UNMEASURED |
| GroES7 | 1.236 | — | 3.362 | 2.7× | 2398 | 1049 | UNMEASURED |
| ClpB6 | 0.01166 | 0.0090–0.049 | 0.1062 | 9.1× | 254,343 | 33,203 | UNMEASURED |
| DnaJ2 | 0.201 | — | 0.339 | 1.7× | 14,731 | 10,414 | UNMEASURED |
| GrpE2 | 1.727 | — | 3.207 | 1.9× | 1717 | 1100 | UNMEASURED |
| trigger factor | 14.58 | 6.41–21.8 | 7.90 | 0.54× | 203 | 446 | UNMEASURED |
| HtpG2 | 0.613 | — | 0.984 | 1.6× | 4839 | 3584 | UNMEASURED |

The machines differ in abundance by four orders of magnitude. Their responses to phase run in opposite directions: ClpB rises 9×, while trigger factor falls 2×.

`CF_SUMMED_POOL` measures the cost of summing:

| pool used | O2 σ_DnaK | legacy O3 headroom |
|---|---:|---:|
| legacy 50 µM | 0.015 | 108 |
| DnaK alone | 0.069 | 63.6 |

The legacy summed pool overstates DnaK capacity 4.6×. The ClpB6 value assumes full hexamer assembly, which at about 0.07 µM protomer is not guaranteed.

### G7.3: effect of each mismatch

The outputs are defined in the `mismatch_effects.py` docstring.

- **O1**: fraction of standing chains with at least one MS-detectable substitution.
- **O4**: P_T / pool of one machine.
- **O2**: σ_DnaK = s_P·φ_sub / (k_cat·DnaK_T).
- **O2s**: the break-even synthesis rate k_cat·DnaK_T / φ_sub.
- **O3**: legacy headroom.
- **O5**: dilution share.
- **O6**: DnaK bound fraction.

Every effect varies a MISMATCHED input, so **every effect range is a SENSITIVITY** (or LEGACY_CONDITIONAL for O3), never a prediction.

**Reference outputs** carry their own kind. Only exponential O4 and s_P are DERIVED_OUTPUT, because every input they use is MATCHED. Exponential O1 is a SENSITIVITY, because f is MIXED-phase.

**EXPONENTIAL** (reference: MATCHED values; f = eTEL zero-inclusive, k_cat = 0.04 s⁻¹ and K = 1 µM as declared reference points)

| effect | mismatched quantity, span | output: reference → range |
|---|---|---|
| X_F_ETEL | f, 5.66e-5 … 5.09e-4 (four MIXED-phase aggregates; MS lower bounds, **open above**) | O1: 0.063 → 0.015–0.122 (8.3×)<br>O2: 0.069 → 0.016–0.134<br>O3: 63.6 → 30.6–286 |
| X_F_ETEL (secondary) | 80 per-dataset rates, 6.7e-6 … 3.3e-3 | O1: 0.0018–0.52 (291×) |
| X_KCAT | k_cat, 0.003 … 1.0 s⁻¹ (every in-vitro DnaK rate on file) | O2: 0.0028–0.92 (333×)<br>O3: 3.6–2508 (689×) |
| X_K | K, 0.06 … 107 µM (in vitro) | O6: 0.51–0.999<br>O3: 10.6–67.5 |
| X_PT_MILO | P_T, 2965 … 6642 µM (generic estimate) | O4 and O2: ×2.24<br>O3: 32.3–63.6 |
| X_PT_BNID | P_T, 2965 → 4000 µM (B/r, 40 min doubling, assumed 1 µm³) | O4 and O2: ×1.35<br>O3: 49.4–63.6 |
| X_N_GENOME | genome-unweighted vs copy-weighted N | O1: 0.063–0.073<br>O3: 54.3–63.6 |
| X_CODON_GENOMIC | genomic vs translated codon weights (conditional-mean aggregate, the largest mover) | O1: 0.122–0.148 (1.22×); 1.01–1.22× across aggregates |

**STATIONARY** (reference: Stikeleather f; everything else MISMATCHED or illustrative)

| effect | mismatched quantity, span | output: reference → range |
|---|---|---|
| S_N | copy-weighted 1 d vs genome-unweighted N | O1: 0.310–0.393<br>O2s: 3.6–4.6 µM s⁻¹ |
| S_KCAT | k_cat, 0.003 … 1.0 s⁻¹ | O2s: 0.34–114 µM s⁻¹ |
| S_MEDIUM_DnaK | DnaK 16.5 … 75.5 µM, using the exponential LB:glucose scale (**borrowed, and confounded with growth rate**) | O4: 47–214<br>O2s: 2.1–9.7 µM s⁻¹ |
| S_MEDIUM_* (7 other machines) | same construction | e.g. ClpB O4: 14,127–87,317 |
| S_MEDIUM_P_T | P_T 3118 … 3852 µM (scale is a normalization artefact) | O4 DnaK: 88–109 |
| S_NORM | normalization factor 0.5–2 (ILLUSTRATIVE) | O4: **invariant**<br>O2s: 2.3–9.1 µM s⁻¹ |
| S_K | K, 0.06 … 107 µM | O6: 0.53–0.999 |
| S_MU | net growth −0.013 … −0.007 h⁻¹ | O5 = 0 (no dilution); exponential glucose gives 0.35 at the same illustrative k_d |

The sampling error on the matched stationary input is small by comparison: at f ± 2 SE, O1 = 0.295–0.325.

### Stikeleather kept separate from exponential capacity parameters

The Stikeleather value is the STATIONARY bundle's only MATCHED entry. It is rejected in EXPONENTIAL under any flag or route, as the tests show.

`CF_STIKELEATHER_INTO_EXP` quantifies the prohibited borrowing. It is reported and never used:

- f would be 3.6–32× the eTEL aggregates;
- O1 would be 0.35, against 0.015–0.12;
- the legacy O3 headroom would be 7.5, against 30.6–286, a 4–38× cut.

`CF_EXP_POOLS_INTO_STAT` quantifies the reverse borrowing. It would move O2s from 4.57 to 1.40 µM s⁻¹ and O4 from 99 to 324.

### the legacy headroom against the bundles (LEGACY_CONDITIONAL)

The control reproduces WO-05: as-published parameters with the zero-inclusive input give **63.10** (test `test_legacy_control_reproduces_wo05`).

| legacy run | headroom |
|---|---:|
| as-published parameters (control) | 63.10 |
| correct P_T only (300 → 2965 µM) | 19.65 |
| every matched exponential value substituted (P_T, DnaK alone, T_gen = ln2/µ, copy-weighted N), with the in-vitro k_cat 0.04 s⁻¹ | 63.63 |

The last two rows agree with the control only by **compensation, not confirmation**. Across the in-vitro k_cat span alone the headroom moves 689×. The legacy headroom is not constrained by any matched data.

## gate evaluation

| gate | evidence | met |
|---|---|---|
| G7.1 | The bundles are separate. Every row has a `condition_tag` (phase;strain;medium;T) and a three-word `match_flag`. Each of the 27 parameters appears in both bundles with exactly one MATCHED or UNMEASURED row. Tests: `test_g71_both_bundles_and_flags`, `test_every_param_accounted_in_both` | yes |
| G7.2 | Phase is re-derived from the source, and an opposite-phase value without MISMATCHED fails. The negative controls are stationary→exponential marked MATCHED, the same with a false `value_phase`, exponential→stationary marked MATCHED or UNMEASURED, a flag without PHASE_OPPOSITE, MIXED or in-vitro values marked MATCHED, Stikeleather in EXPONENTIAL, and a hyphenated phase string. All are caught. The positive control (properly flagged) raises no phase error. The real table has 0 errors | yes |
| G7.3 | All 28 MISMATCHED rows link to effects that vary their own parameter in their own bundle. Each effect has at least one output with a finite lo ≤ hi. No effect output may be labelled DERIVED_OUTPUT (checker plus test) | yes |

### findings on the gates (the gates themselves were not changed)

- **G7.3, in vitro.** The in-vitro ranges span the *saved in-vitro constants*. They are not a bound on the in-vitro → in-vivo error, which no saved record measures. The true size of the IN_VITRO mismatch is unbounded by data.
- **G7.3, stationary medium and normalization.** These spans use a borrowed exponential scale or an illustrative factor, and are labelled accordingly.
- **G7.2, trust boundary.** Phase and condition come from WO-06's fields, so the check cannot detect a source whose own phase statement is wrong.
- **Value check.** The audit-row value check still accepts any number that appears in the audit row's value, range, measured or supports fields. For L02, "22" or "50" (peptide lengths and nM) would pass. The committed TSV is additionally pinned by the regenerate-and-compare check, so a hand edit fails. A bad literal in `build_bundles.py` would pass only if its number also appears in that audit row's text.

## independent review

A read-only reviewer (separate agent) independently recomputed from the xlsx and the saved XML:

- all 8 exponential pools, P_T, µ and s_P;
- the copy-weighted N;
- the stationary rows;
- the Stikeleather conditions;
- the Schmidt stationary medium;
- O1, σ and the S_NORM invariance.

It found no fabricated value. It reported 5 ERRORs and 4 validator bypasses. All were verified against the saved text and fixed:

1. **Stikeleather was MATCHED as if it were a stationary synthesis error rate.** It is a standing-proteome frequency. Fixed by the two-quantity split. The anchor strain was also corrected from "Xac-derived" to Xac.
2. **Sign error in O5.** |−0.01 h⁻¹| was used as a dilution rate. Fixed: negative growth gives no dilution sink.
3. **O1 in EXPONENTIAL was labelled DERIVED_OUTPUT** although f is MIXED-phase. Fixed: reference kinds are now computed from input dependencies.
4. **S_N O1 was labelled DERIVED_OUTPUT with a DATA span** although both endpoints are mismatched. Fixed: all effect outputs are SENSITIVITY or LEGACY_CONDITIONAL.
5. **Stationary O1 and O2s used abundance weights under a "newly made chains" label.** Fixed: O1 is now defined on the standing proteome, and O2s lists the synthesis-rate substitution among its held assumptions.
6. **Validator bypasses** (each is now a regression test):
   - MATCHED conditions and WO-06 status were self-reported. A yeast k_d, a B/r P_T, or a retyped E05 could pass. Now they come from the source.
   - `vector:` values were unchecked.
   - The phase parser misread "post-exponential", and the Stikeleather guard was string-only.
   - The effect was not tied to its parameter.
   - The ×1e3 rule applied to all units; it is now restricted to concentrations.

Reviewer weaknesses also adopted:

- BNID axes are completed (medium, growth rate, assumed volume).
- The codon effect is shown for the aggregate that moves most.
- The eTEL span is marked open above.
- The S_MEDIUM scales are labelled as confounded with growth rate, and as a normalization artefact for P_T.

## adversarial section: is either bundle complete enough for quantitative physiological predictions?

**No. Neither bundle is.** WO-07 PASS certifies bookkeeping and mismatch transparency, not biological completeness.

### EXPONENTIAL

- 13 of 27 parameters are MATCHED, all abundances or growth.
- The WO-02 model needs K_M, k_cat, φ, k_d, k_a, k_dis, k_dA, k_mis and ε. **None is matched.**
  - Seven have no value in any saved record.
  - DnaK k_cat and K exist only in vitro.
  - ε needs p_misfold and a per-synthesis error rate, and both are UNMEASURED. The standing frequency itself is MIXED-phase only.
- So no steady state, stability classification, fold location or headroom of the rebuild model can be computed for *E. coli* in exponential growth. Any such number is a sensitivity.
- **The only DERIVED_OUTPUTs are O4 (chains per machine unit) and s_P.**
- O1 spans 8× across eTEL aggregates and 291× across datasets.
- O2, the simplest capacity-vs-load index, spans 333× across the in-vitro k_cat span alone.

### STATIONARY

- 1 of 27 parameters is MATCHED, and it is a standing frequency, not the per-synthesis rate a burden flux needs.
- s_P is UNMEASURED, and s_P = µ·P_T fails at negative net growth. **No burden flux can be formed from the stationary data.**
- Pools come from another strain, medium and starvation time, and their absolute values inherit an exponential normalization.
- **No stationary output is DERIVED.**

### weaknesses that survive even in MATCHED rows

1. **Volume.** Concentrations are per *calculated* whole-cell volume, including the periplasm. Absolute µM values are lower bounds on cytoplasmic ones.
   - O4 uses whole-cell P_T. Periplasmic mass is 15% of total in stationary phase against 6% in LB (Schmidt), so part of the exponential-vs-stationary O4 difference is compartment, not load.
2. **Coverage.** About 55% of genes are quantified, so P_T is a lower bound on chains.
3. **Untested closures.** s_P = µ·P_T and "copies = synthesis weights" both neglect turnover, and their error is not quantified here.
   - A saved record with direct *E. coli* synthesis rates exists: Li 2014, PMC4006352 in WO-06/records, ribosome profiling. It could test these closures and is the first thing WO-08 should use.
4. **Scalar f.** O1 applies one f to all chains. Per-codon and per-protein heterogeneity would shift it, in an undetermined direction.
5. **One anchor.** The exponential anchor is one strain in one medium. Across 20 exponential conditions DnaK varies 3×, so "exponential *E. coli*" is not one parameter set either.
6. **O2 assumptions.** O2 assumes one DnaK cycle per client. The in-vivo DnaJ:DnaK ratio of 1:27 makes the in-vitro, DnaJ-triggered T→R rate a questionable in-vivo turnover.

**For WO-08.** MISMATCHED spans must be sampled as spans with their status carried through (G8.3), not as calibrated distributions. In stationary phase, anything that needs s_P or a per-synthesis error rate cannot be sampled from data at all.

## limitations

- Only files saved by WO-05/WO-06 were used. PRIDE metadata for the 80 Landerer datasets was not fetched. That is the concrete route to a phase-annotated exponential standing frequency.
- One source (Schmidt 2016) supplies every pool in both phases, so lab and strain effects cannot be separated from phase effects.
- Schmidt accessions are mapped to the MG1655 CDS. BW25113 is a K-12 derivative; its sequence differences are not checked.

## reproduce

```
cd proteostasis_rebuild
PYTHONDONTWRITEBYTECODE=1 python WO-07/derived.py
PYTHONDONTWRITEBYTECODE=1 python WO-07/build_bundles.py
PYTHONDONTWRITEBYTECODE=1 python WO-07/mismatch_effects.py
PYTHONDONTWRITEBYTECODE=1 python WO-07/check_bundles.py
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s WO-07 -p 'test_wo07.py' -v
```

VERDICT: PASS — G7.1, G7.2 and G7.3 are met. PASS certifies that the bundles are phase-separated, fully flagged, machine-resolved, and that every mismatch is quantified as a labelled range. It does **not** certify either bundle as complete: EXPONENTIAL matches no kinetic parameter, and STATIONARY matches one standing error frequency and nothing else.
