# STATE — running log

Append-only log of loop transitions. Machine-readable mirror: `loop_state.json`.

## 2026-09-27

- scaffold written: WORK_ORDERS.md (gates fixed before analysis), STATE.md,
  STATUS.md, CLAIM_REGISTER.tsv, loop_state.json.
- pre-existing repo state noted: `git status` in proteostasis_law showed the
  `proteostasis-paper` submodule already modified BEFORE this run began. Not
  touched by this run.
- WO-00 started.
- WO-00 PASS (33 claims, anchors verified; 30 legacy files hashed; 7/7 legacy
  numbers reproduced at rel. err 0). One self-caught gate loosening fixed before
  verdict (see WO-00/REPORT.md). WO-01 started.
- WO-01 PASS (legacy dimensionally consistent; non-conservative Phi inflow,
  missing dilution, and a 1/ln2 synthesis-rate overstatement found; C34, C35
  added). WO-02 started.
- WO-02 PASS (conservative model; totals proved symbolically; orthant invariant;
  200 integrations conserve totals; legacy Phi inflow has no donor and is 48%
  of inflow at the legacy threshold). Self-caught hard-coded G2.5 replaced by a
  computed test. WO-03 started.
- WO-03 PASS (exact finite-pool binding; driven cycle set by K_M not K_d;
  four-state cycle shows drive-dependent ultra-affinity with detailed-balance
  control; nascent competition makes theta an output; tQSSA err < 0.5%).
  WO-04 started.

## 2026-09-28 (backfilled 2026-09-29 from WO reports; not logged at the time)

- WO-04 PASS (see WO-04/REPORT.md). WO-05 started.
- WO-05 BLOCKED on G5.5: Data_S2 means judged conditional on detection from
  (sd/se)^2 < 80; zero cells and denominators judged unavailable. Loop stopped.

## 2026-09-29

- WO-05 reopened adversarially on user instruction.
- Found: Landerer Data_S4 (same supplement as Data_S2) has per-dataset
  base_count/error_count for all 80 E. coli datasets. Data_S2 reproduced to
  ~1e-9 ONLY when covered zero-error cells are dropped (3015 such cells; 19
  cells without coverage). n = (sd/se)^2 = datasets with >=1 detection.
  Upstream deTEL global_report.py drops rate-0 rows via log10 -> NaN -> dropna.
  Prior inference confirmed; prior block reason (data unavailable) was false.
- G5.5 computed in the legacy model: x24.817 reproduced; remove (1-S) only
  x17.320; ln2/T_gen only x35.876; both x25.064; with zero-inclusive input
  x63.10; PSM-pooled x112.87; per-dataset spread x4.6 to x2370. All labelled
  eTEL-aggregate, legacy-model-conditional, not an E. coli baseline.
- Stikeleather 2026 NAR (doi:10.1093/nar/gkag674) recorded as a separate
  STATIONARY-phase estimate (1.82e-3 /codon); no headroom computed from it.
- CLAIM_REGISTER: C05 CORRECTED, C06 REJECTED, C10 CONDITIONAL, C21-C24
  CORRECTED. Provenance check 35/0; legacy hashes 30/30 unchanged.
- WO-05 PASS. WO-06 eligible, NOT started in this session (a citation audit
  run short on time would risk VERIFIED marks without real record checks).
  Loop paused.
- WO-06 started (same day, later session).
- WO-06 PASS. 59-row parameter_audit.tsv; every citation checked against a
  saved PubMed/PMC/BioNumbers/UniProt record or an ecitmatch probe with a fired
  control; check_audit.py 0 errors; 22 tests OK; WO-01..05 suites OK;
  provenance 35/0; legacy hashes unchanged.
- Key findings: none of the 13 legacy two-pool scalar parameters is verified
  for E. coli balanced growth. Prot_tot 300 uM is ~10x low (BNID 104726 4 mM;
  Schmidt 2016 glucose 2.97 mM). A_max has no source (2 MISCITED, 2
  UNMATCHED). Pierpaoli "1997 EMBO J" K_d/k_obs ranges are verbatim in
  Pierpaoli 1998 Biochemistry (in-vitro R-state peptide binding, 25 C).
  D&W "2009 Cell", Ciryam "PNAS 110:E3453", Bednarska "Mol Cell 52:617"
  miscited; Christiano 2014 is yeast. Stikeleather SE 5.92e-5 now verified.
- One self-caught probe flaw fixed before verdict (dummy author in ecitmatch);
  independent reviewer corrected 5 labels (no fabricated values).
- Proposed register changes in WO-06/claim_updates.tsv (CLAIM_REGISTER.tsv and
  STATUS.md not edited: preexisting files preserved per instruction).
- WO-07 eligible, NOT started. Loop paused.

## 2026-09-30

- WO-07 started.
- WO-07 PASS. bundles.tsv: 80 rows, 27 params per bundle. EXPONENTIAL 13
  MATCHED (all abundances/growth; no kinetic parameter), 12 MISMATCHED, 13
  UNMEASURED. STATIONARY 1 MATCHED (Stikeleather standing substitution
  frequency), 16 MISMATCHED, 25 UNMEASURED. check_bundles.py 0 errors; 38
  tests OK (23 negative controls); WO-01..06 suites OK; provenance 35/0;
  legacy hashes 30/30 unchanged; pipeline rerun byte-identical.
- Key findings: Schmidt stationary absolute concentrations carry the
  glucose-exponential volumetric normalization (hidden cross-phase borrowing;
  ratios to P_T cancel it). MS error rates are standing-proteome frequencies;
  the per-synthesis error rate is UNMEASURED in both phases. Legacy summed
  C_tot overstates DnaK 4.6x. Legacy headroom with matched values (63.6) equals
  as-published (63.1) only by compensation; P_T alone gives 19.7; in-vitro
  k_cat span moves it 689x. Stikeleather-into-exponential would cut it 4-38x
  (counterfactual, prohibited in bundles).
- Independent reviewer: no fabricated value; 5 ERRORs + 4 validator bypasses
  verified and fixed before verdict (Stikeleather semantics, |mu| sign error,
  DERIVED_OUTPUT mislabels, self-reported MATCHED conditions); each bypass is a
  regression test.
- PASS certifies bookkeeping/mismatch transparency, NOT biological
  completeness: neither bundle supports quantitative physiological prediction.
- WO-08 eligible, NOT started. Loop paused.
