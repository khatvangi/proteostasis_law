# WO-06 — literature parameter audit (2026-09-29)

## gates (fixed before analysis, WORK_ORDERS.md)

- G6.1 every parameter used by the legacy model and the rebuild has: value,
  unit, organism, condition (phase/medium/temperature), source, what was
  actually measured, verification status.
- G6.2 every citation is checked against a bibliographic record (PubMed or
  equivalent); unmatched → UNVERIFIED, mismatched → MISCITED. Nothing is marked
  VERIFIED from memory.
- G6.3 an automated check fails if any row lacks organism/condition/status.

## deliverables

| file | what it is |
|---|---|
| `parameter_audit.tsv` | 59 rows, 28 columns, one row per (parameter, source) pair |
| `build_audit.py` | writes the TSV; every computed number is read from the two JSON files below |
| `fetch_records.py` | the only networked step. Saves PubMed efetch records, PMC full text and NCBI ecitmatch results (with controls) to `records/` |
| `records/` | 27 PubMed records, 17 PMC full texts, `citmatch.json`, BioNumbers BNID 104726 page, UniProt UP000000625 record, Schmidt 2016 supplementary xlsx (sha256 `3280a13f…`) |
| `schmidt_pools.py` → `schmidt_pools.json` | measured chaperone pools and total protein in 22 conditions |
| `audit_checks.py` → `audit_derived.json` | derived checks: codon-usage recount, code-implied S, k_deg half-lives, unit conversions, pool sums |
| `check_audit.py` | the G6.3 validator (exit 1 on any error) |
| `test_wo06.py` | 22 tests, including negative controls |
| `claim_updates.tsv` | proposed CLAIM_REGISTER status changes. The register itself is not edited, because preexisting files are preserved |

## method

1. **Inventory.** Every numeric parameter was read from code, not from summaries:
   - legacy `two_pool_ode.py` (`Params` and `sample_params`), `arithmetic_stress_test.py`, `paired_mc.py`
   - envelope scripts 06/09/11/12
   - rebuild WO-02 `SCENARIO`, WO-03 cycle and competition rates, the WO-04 sampling domain, and the WO-05 error inputs

   `test_wo06.G61Coverage` parses those sources and fails if any parameter name is missing from the table.
2. **Bibliographic check.** Each citation *as the legacy wrote it* was sent to
   NCBI ecitmatch (journal|year|volume|page|author). Each probe expected to
   fail has a same-form **control** (same author at the real journal/year, or
   same journal/volume at a real page) that must resolve. All 6 controls fire,
   so the NOT_FOUNDs are informative. The PubMed record of each matched paper is
   saved, and the validator compares its journal|year|volume with the row.
3. **Value check.** A value counts as "located" only if the quoted text was found in
   the saved abstract, PMC full text, or data file. Where only the abstract was
   retrievable (Lorimer 1996, Mogk 1999, Hoffmann 2001, Belle 2006) and the
   value is not in it, the row says `value_located = NO`.
4. **Status vocabulary.** It separates bibliographic error from biological mismatch:
   - `citation_check`: MATCHED / MISCITED / AMBIGUOUS / UNMATCHED / NO_CITATION / LOCAL_DATA. This is bibliographic only.
   - `org_cond_match`: MATCHED / ORGANISM_MISMATCH / CONDITION_MISMATCH / QUANTITY_MISMATCH. This compares the source with the model's use.
   - `verification_status`:
     - VERIFIED
     - MISMATCHED_CONDITION: right paper, value present, but another organism, condition or quantity
     - NOT_SUPPORTED: right paper, but the value is not in it or is contradicted by it
     - MISCITED: wrong bibliographic details, or the paper is credited with a property it does not have
     - UNVERIFIED
     - ASSUMED / ILLUSTRATIVE
   - `measurement_type` UNKNOWN marks rows whose citation cannot be matched at all.

## results

### status counts (59 rows)

| status | all rows | legacy-model rows (L*, A*, excluding WO-06 check rows) |
|---|---:|---:|
| VERIFIED | 10 | **1** (the UniProt length distribution) |
| MISMATCHED_CONDITION | 5 | 0 |
| NOT_SUPPORTED | 1 | 1 |
| MISCITED | 6 | 6 |
| UNVERIFIED | 11 | 9 |
| ASSUMED | 12 | 8 |
| ILLUSTRATIVE | 14 | 0 (all rebuild) |

**None of the legacy two-pool model's 13 scalar parameter values is verified
for E. coli in balanced growth.** The 10 VERIFIED rows split as follows:

- 1 legacy data file (UniProt lengths)
- 1 rebuild attribution (De Los Rios & Barducci)
- 1 stationary-phase error rate (Stikeleather)
- 7 reference values added by this audit for WO-07: BNID 104726, plus Schmidt 2016 total protein, DnaK, GroEL, ClpB, co-chaperones and growth rates

Milo 2013 is MISMATCHED_CONDITION: it is a generic, not E. coli-specific, estimate.

### outright miscitations vs organism/condition mismatches

**Miscited** (the reference as written does not exist, resolves to another
paper, or credits the paper with something it does not say):

| row | legacy citation | what the record shows |
|---|---|---|
| L02, L04 | "Pierpaoli et al. 1997 EMBO J" for k_obs_max 3e-3–8.4e-2 /s and K_d 0.06–2 uM | No such record. Both ranges appear **verbatim** in Pierpaoli, Gisler & Christen **1998 Biochemistry 37:16741**. The actual 1997 Pierpaoli paper is J Mol Biol 269:757 |
| L09 | "Drummond & Wilke 2009 Cell" for p_misfold 0.3 | No such record. The papers are Cell 2008;134:341 and Nat Rev Genet 2009;10:715. The legacy *range* 0.1–0.5 matches Cell 2008's "~10–50% of random substitutions disrupt protein function". That is a secondary citation about loss of function from mutations, not misfolding per mistranslation. Neither paper gives a misfolding probability, and the 2009 review says the failure rate "remains essentially unknown" |
| L13a | "Ciryam 2013 PNAS 110:E3453" | No such record. Ciryam's only 2013 PNAS paper is 110:E132 (cotranslational folding). The supersaturation paper is Cell Rep 5:781, human, with no threshold |
| L13b | "Bednarska 2013 Mol Cell 52:617" | Mol Cell 52:617 is Aakre et al., a *Caulobacter* toxin that blocks replication. Bednarska 2013 is a Microbiology review with no 20% threshold in its abstract |
| L01b | "Christiano 2014 for bacteria" | The record is the yeasts *S. cerevisiae* and *S. pombe* |

**Unmatched** (no record found; this is not proof the paper doesn't exist):
- Yamanaka 2017 Curr Biol (L13c)
- Stirling 2018 Cell Rep 25:2242 (L13d)
- "BioNumbers / Milo & Phillips" without a BNID, used for Prot_tot and N (L06a, L08a)
- "Cohen / Meisl / Knowles" (L05c)

**Right paper, wrong conditions or quantity** (the source is described correctly but measured something else; the MISMATCHED_CONDITION label is used only where the value was actually seen):
- Belle 2006: yeast bulk native-protein half-lives (L01a). The value is not in the abstract, so the status is UNVERIFIED.
- Pierpaoli 1998: in vitro, 25 °C, nucleotide-free DnaK, short peptides (L02, L04). This is *also* miscited.
- Landerer Data_S2: pooled over 80 datasets and selected on detection (E01, E02).
- Genomic rather than translation-weighted codon usage (E05).
- Mogk 1999: heat-shock recovery (L11).

**Right paper, but the value is not in it or is contradicted:**
- Upadhyay 2012: the full text was searched and contains no rate constant (L05a, NOT_SUPPORTED).
- Hoffmann, Posten & Rinas 2001 (L05b): the abstract says the aggregation step "was found to follow **first order** kinetics" under one production assumption. It is UNVERIFIED rather than NOT_SUPPORTED because only the abstract was retrieved, and fitted constants may be in the full text.

### the anchors the task asked to revisit

1. **Total protein, legacy 300 µM.** Wrong by about 10×.
   - BNID 104726 gives **4 mM** (E. coli B/r, balanced growth, glucose minimal, 37 °C; calculated from Neidhardt 1996).
   - Schmidt 2016 glucose gives **2.97 mM** (BW25113, whole-cell volume). Across the 22 exponential conditions it gives 2.72–3.24 mM, but only the glucose value is independent: the others were adjusted assuming constant volumetric concentration.
   - Milo 2013 gives 2–4 million proteins per µm³, i.e. **3.3–6.6 mM**. This is a generic cross-organism estimate.
   - The legacy value is 9.9× (Schmidt), 13.3× (BNID) and 11–22× (Milo) low.
   - Where "300" came from is not established. Two leads only: the BNID 104726 page points to BNID 104678, a "300 mM" metabolite pool, and Zimmerman & Trach 1991 give macromolecules at 0.3–0.4 g/ml. The rebuild placeholder P_T = 3000 µM is consistent with these but stays ILLUSTRATIVE.
   - Leverage: Prot_tot multiplies k_agg in the aggregation drain and sets M for chaperone competition.
2. **A_max = 0.25.** No verified source: all four citations fail (2 MISCITED, 2 UNMATCHED). It is ASSUMED. Because the legacy threshold *is* this gate (WO-04), the legacy headline threshold rests on an unsourced number.
3. **DnaK/GroEL/ClpB pools, K_d and kinetic constants.**
   - Measured pools (Schmidt 2016, BW25113, µM per whole-cell volume):

     | protein | glucose | LB | exponential range | stationary |
     |---|---:|---:|---:|---:|
     | DnaK | 10.9 | 23.2 | 8.3–25.1 | 35 |
     | GroEL protomer | 13.3 | 25.5 | 11.7–36.5 | — |
     | GroEL 14-mer | 0.95 | 1.8 | — | — |
     | ClpB hexamer | 0.012 | — | 0.009–0.049 | 0.09–0.11 |

   - DnaK:DnaJ is about 27:1 in glucose.
   - The legacy "combined pool ≈ 50 µM" is near the top of the DnaK + GroEL *protomer* sum (24 µM glucose, 49 µM LB). It is 2–4× the functional-unit sum (12 µM glucose, 25 µM LB), and it treats machines with different clients as one pool.
   - The "DnaK peaks at 30–50 µM" figure is uncited. It matches stationary phase only.
   - Lorimer 1996's record says GroEL/ES suffice for "no more than 5%" of proteins. The "GroEL ≈ 30 µM" figure is not located. Among exponential conditions, only 42 °C glucose (36.5 µM protomer) exceeds it.
   - K_d 0.06–2 µM is the nucleotide-free R-state peptide K_d in vitro. The ATP state in the same abstract is 2.2–107 µM.
   - k_obs_max is an in-vitro peptide *binding* rate, not a folding or turnover rate. The DnaK-ATP kobs 0.001–7.9 s⁻¹ in the same abstract is quoted from Gisler 1998, not measured there. The rate-limiting DnaK cycle step (T→R) is 0.04 s⁻¹ in vitro (Pierpaoli 1997 JMB).
4. **Aggregation, disaggregation, degradation, clearance rates.**
   - k_agg has no verified second-order constant (see above).
   - k_clear: the "95% in 2 h" figure is not located, and the setting is heat-shock recovery with induced chaperones. The ClpB hexamer is about 0.01 µM in exponential growth.
   - k_deg: both citations are yeast, and one is miscited as bacteria. The legacy's own range (t½ 11.6–116 min) mostly lies outside the "1–10 h" consensus it cites, and the baseline (38.5 min) lies fully outside it. Goldberg 1972 (E. coli) supports faster degradation of abnormal proteins, but no rate constant was located.
   - The rebuild's k_d, k_a and k_dis inherit these values and remain ILLUSTRATIVE.
5. **p_misfold = 0.3.** Miscited (there is no D&W Cell 2009).
   - The 0.1–0.5 range matches a Cell 2008 statement about random *mutations* disrupting *function*.
   - The "30%" in NRG 2009 is the share of new proteins rapidly degraded in one early study; a later study found "at most a few percent".
   - Neither is a misfolding probability per mistranslation event, and the same authors state the failure rate is unknown.
6. **N = 300.** Labelled "median", but the legacy's own UniProt file gives median 271 and mean 307.6. D&W 2009 says 335 codons on average. None of these is synthesis-weighted, which is the burden-relevant weighting.
7. **Generation time and dilution.** T_gen is uncited (ASSUMED). The 30–180 min envelope lies within measured batch doubling times (0.36–2.7 h, Schmidt Table S23). The legacy has no dilution term (WO-01).
8. **Synonymous factor S = 0.30.** Uncited. Under uniform single-nucleotide misreading, the standard code gives 0.245–0.255 (all positions) or 0.69–0.72 (third position only). S depends on the misreading spectrum, which is not uniform, and 0.30 matches neither limit.
9. **Landerer and Stikeleather error inputs.**
   - Landerer: record matched, values reproduced (WO-05). The condition is MIXED, and the Data_S2 mean is selected on detection, so it is MISMATCHED_CONDITION as an E. coli exponential rate.
   - Stikeleather: 1.82e-3 /codon, **SE 5.92 × 10⁻⁵ now read in the saved PMC XML**. This closes the WO-05 open item about the SE exponent. Conditions: Xac strain, LB (Miller), 37 °C, stationary. It is VERIFIED for its stated stationary role, and MISMATCHED if combined with exponential parameters.

### consequences for earlier WOs (recorded here; their files are untouched)

- The WO-03 attribution to De Los Rios & Barducci 2014 eLife 3:e02218 is now
  VERIFIED against the record (R16).
- The WO-03 four-state rates are ILLUSTRATIVE, as WO-03 said. Their hydrolysis
  rate (10 s⁻¹) is 250× the measured rate-limiting T→R step (0.04 s⁻¹), so the
  0.21 µM ultra-affinity *magnitude* should not be quoted. Only its existence
  was claimed.
- The WO-05 Stikeleather SE exponent is verified (see item 9).

## gate evaluation

| gate | evidence | met |
|---|---|---|
| G6.1 | 59 rows. Every row has value, unit, organism, strain, phase, medium, temperature, source, measured_in_source, measurement_type and verification_status. Coverage test: every parameter name in legacy `Params`/`sample_params`, the arithmetic `Baseline`, WO-02 `PARAM_NAMES`/`SCENARIO`, the WO-03 cycle and competition dicts, and the WO-04 `sample` keys appears in the table | yes |
| G6.2 | Every citation attached to a parameter (LITERATURE_ANCHORS, WO-03, WO-05) has a saved PubMed record or an ecitmatch result with a fired control. Unmatched → UNVERIFIED (5 rows), bibliographically mismatched → MISCITED (5 rows, plus Christiano misattributed as bacterial). Rule R9 enforces exactly this mapping. The validator also rejects VERIFIED without a matched saved record, a located value, and matched organism/condition | yes |
| G6.3 | `check_audit.py` exits 0 on the table. `test_missing_organism_fails` blanks each of organism, strain, phase, medium, temperature, unit, source, measurement type and status in turn, and each is caught | yes |

**Finding on G6.1 (not a change to the gate).** For a parameter with no source,
"organism/condition" cannot be a real value. The table writes the explicit
token `NO_SOURCE`, and the validator allows it only on ASSUMED or ILLUSTRATIVE
rows. I read this as meeting the gate's intent (no silent gaps), not its most
literal reading. A stricter reader could call those rows ungated. There are 26
such rows, and all are labelled as placeholders.

## adversarial self-review

An independent read-only reviewer checked 25 rows against the saved records.
It ran its own positive-control NCBI queries and recomputed the Schmidt
numbers by hand: DnaK glucose 10.880 µM and GroEL 13.298 µM, both matching
`schmidt_pools.json`. It found **no fabricated value**. It found five
mislabelled rows and three wording errors, all verified against the saved text
and corrected:

1. **L09 was wrong in my favour of harshness.** I wrote "no 0.3 in either
   paper" because my search looked for the point value and missed the range.
   Cell 2008 does contain "~10–50%" (the legacy range), for a different
   quantity. value_located is now YES and the interpretation is rewritten. The
   MISCITED status stands.
2. **L01a** was labelled MISMATCHED_CONDITION without the value having been
   seen, which contradicts the vocabulary. It is now UNVERIFIED, and new rule
   R5b (with a negative-control test) makes that label require a located value.
3. **L05b** was NOT_SUPPORTED on an abstract only. It is now UNVERIFIED; the
   first-order statement is kept as a note.
4. **L06c (Milo)** was VERIFIED with MATCHED conditions. A generic
   cross-organism estimate is MISMATCHED_CONDITION.
5. **L03a note** said 30 µM GroEL is matched "in LB/42C". Only 42 °C glucose
   exceeds it (LB is 25.5 µM).
6. The Schmidt total was labelled "dataset 2" but sums the combined table
   (dataset-1-only rows are about 0.5% of copies). The label was fixed; the
   value is unchanged.
7. E01: Landerer's 20–23% is not stated as E. coli-specific. The note was fixed.
8. L02: the DnaK-ATP kobs is quoted from Gisler 1998, not measured in the
   cited record. The evidence now names both candidate "release" sources.

My own falsification attempts on the ten highest-leverage or weakest-provenance
parameters (A_max, Prot_tot, k_agg, k_clear, C_tot, K_d, k_obs_max, p_misfold,
S, error inputs):

- *Could "MISCITED" be a matcher artefact?* The first ecitmatch run put a dummy
  author in every probe. For journal+year-only probes (Pierpaoli EMBO J, D&W
  Cell 2009, Yamanaka), NOT_FOUND was therefore meaningless. I caught this
  before the verdict, reran with real authors, and added a same-form control
  for each failing probe. All six controls fire, so the NOT_FOUNDs are
  informative. They still do not prove that a paper is absent from every
  database.
- *Is Prot_tot really wrong, or a definition difference?* The legacy calls it
  "total soluble cytoplasmic protein". Even a soluble fraction cannot be 10×
  below total protein; none of the three independent records is below 2.7 mM.
- *Is A_max perhaps right even though uncited?* This audit cannot say. It shows
  only that no source supports it. The one nearby verified number (Mogk 1999:
  15–25% of detected *species* aggregate under heat stress, and cells recover)
  is a species count, not a mass fraction, so it neither supports nor refutes
  0.25.
- *Are the Schmidt pools the right comparison for C_tot?* Only partly. They are
  MS estimates over calculated whole-cell volumes, in one strain. The "free,
  folding-competent" pool the model needs is not measured anywhere here. They
  bound total abundance; they do not calibrate C_T.
- *Did I weaken a gate?* G6.1's organism/condition field for unsourced
  parameters is filled with an explicit NO_SOURCE token, not a value. That is
  recorded above as a finding. The validator restricts it to placeholder rows,
  so it cannot hide a sourced parameter.

## limitations

- Several papers were checked only at abstract level because no
  machine-readable full text was retrievable: Lorimer 1996, Mogk 1999,
  Hoffmann 2001, Belle 2006, Bednarska 2013 and Goldberg 1972. For these,
  "value not located" means *not in the retrievable text*, not *absent from
  the paper*.
- ecitmatch NOT_FOUND with a fired control is strong evidence that the citation
  as written is wrong. It does not prove that no such paper exists anywhere
  (Yamanaka, Stirling).
- Schmidt 2016 concentrations use calculated whole-cell volumes (Volkmer &
  Heinemann 2011, record saved) and cover about 55% of genes. They are MS
  estimates, not ground truth, and only strain BW25113 is used here.
- The `leverage` column is a qualitative judgement from where each parameter
  enters the legacy headline (WO-04/WO-05). It is not a sensitivity computation;
  that is WO-08.
- The manuscript's 27-item reference list was not audited except where an entry
  sources a parameter.

## reproduce

```
cd proteostasis_rebuild
PYTHONDONTWRITEBYTECODE=1 python WO-06/fetch_records.py      # network; records/ already saved
PYTHONDONTWRITEBYTECODE=1 python WO-06/schmidt_pools.py
PYTHONDONTWRITEBYTECODE=1 python WO-06/audit_checks.py
PYTHONDONTWRITEBYTECODE=1 python WO-06/build_audit.py
PYTHONDONTWRITEBYTECODE=1 python WO-06/check_audit.py
PYTHONDONTWRITEBYTECODE=1 python -m unittest discover -s WO-06 -p 'test_wo06.py' -v
```

VERDICT: PASS — G6.1, G6.2 and G6.3 are met. PASS certifies that the audit was done to the gates. It does not certify any parameter: only 1 of the legacy model's parameter rows is VERIFIED, and that is the length distribution.
