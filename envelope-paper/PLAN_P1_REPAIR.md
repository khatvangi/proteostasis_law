# P1 repair plan — envelope paper to submission

Working document, created 2026-07-31 after the external review and Kiran's triage.
Delete at submission. Rule for this run: **no new analyses.** Anything that
surfaces outside this list goes to `HISTORY.md` and the run continues.

---

## Part 0 — what I verified before planning

### Item 8 is confirmed, and it is worse than the review states

`bounds_summary.json` — the source of every number in Result 2 and Fig. 3a,b — is
written by `scripts/05_fig3_bounds.py`, a **figure script**, which reads three
staged files from `data/raw/`:

```
scripts/05_fig3_bounds.py:47   paired_mc_results.json
scripts/05_fig3_bounds.py:48   arithmetic_results.json
scripts/05_fig3_bounds.py:49   two_pool_results.json
```

None of the 16 pipeline steps computes a bound. So the clean-rebuild verification
proves the assembly is deterministic; it does not show the numbers derive from
anything. The tests assert transcription, not derivation, for the paper's spine.

**Beyond the review's account:** `scripts/06_translation_burden.py:32` also reads
`arithmetic_results.json`, for the proteome length distribution. Result 1 — the
usage-weighted 6.3 × 10⁻⁴, the 613-fold span, the 16–18% figure — therefore rests
on staged output too, not only Result 2. The absolute path at
`scripts/vendor/two_pool_ode.py:69–83` loads the same proteome table from
`/storage/kiran-stuff/proteostasis-P1/`. One fix serves all three: vendor
`ecoli_proteome_lengths.tsv` into `data/raw/` as a genuine primary input.

**Good news the review missed:** the generators exist upstream and are importable.

| Generator | Upstream file | Entry point |
|---|---|---|
| paired Monte Carlo | `proteostasis-P1/paired_mc.py` | `run_mc(n=5000, seed=17)` |
| combinatorial bound + its MC | `proteostasis-P1/arithmetic_stress_test.py` | `part_D(lengths, n=5000, seed=17)` |
| two-pool threshold ensemble | `proteostasis-P1/two_pool_ode.py` | already vendored |
| proteome lengths | `proteostasis-P1/ecoli_proteome_lengths.tsv` | 206 KB, vendorable |

### The July 2026 citation resolves — correction to the triage

`10.1093/nar/gkag674` **does** resolve. It is Stikeleather, Ali, Ho, Licknack &
Lynch, *Translation accuracy in E. coli*, **Nucleic Acids Research 54(13)**,
online 6 July 2026. Verified against the OUP article page and the preprint
(bioRxiv 2025.04.18.649569 = PMC12258882 = PubMed 40661473).

Three consequences:

1. **"Xac" is an *E. coli* strain**, not *Xanthomonas*. The paper studies three
   ribosomal variants (wild-type, restrictive, error-prone) derived from Xac
   *E. coli*. The review's description was accurate and the internal-incoherence
   objection does not hold.
2. **The wild-type rate is 1.82 × 10⁻³ per codon (SE 5.92 × 10⁻⁵)** — confirmed
   from both sources. The operating-point problem is real and precise.
3. **The preprint and the review's citation are the same study**, so there is one
   source, not two independent ones. It is peer-reviewed and from the Lynch lab,
   which is stronger than an unresolvable DOI — but it is not corroboration.

The same paper reports that codons enriched in highly expressed genes are **not**
more accurately translated, contradicting Sun & Zhang 2022
(`10.1126/sciadv.abl9812`), which found the opposite in *E. coli* proteomic data.
Both must be cited; the μ result sits between them.

### New, and it changes what step 5 can conclude

**The 1.82 × 10⁻³ figure is a stationary-phase measurement.** Cultures were grown
overnight to stationary phase. The two-pool model is calibrated for
*exponentially growing* E. coli: GroEL ≈ 30 µM in exponential growth, DnaK
30–50 µM, and the model's clearance term includes growth dilution, which is near
zero in stationary phase. Dropping a stationary-phase error rate into an
exponential-phase parameterization while holding capacity fixed is not a harder
test of the model, it is a mismatched one, and it would manufacture part of the
collapse the triage anticipates.

So step 5 must vary **condition as a triple** — error rate, chaperone pool,
dilution — not the error rate alone. Two defensible evaluations:

- exponential K-12: μ from Landerer (condition to be confirmed), current capacity
  parameters, 6.3 × 10⁻⁴;
- stationary Xac: 1.82 × 10⁻³ with stationary-phase capacity and dilution.

If wild type still lands past its own collapse threshold under a *self-consistent*
stationary-phase parameterization, the reductio stands and the framing changes to
relative structure. If it does not, the model survives and the higher rate becomes
a second operating point rather than a refutation. Either way the answer is
earned rather than assumed.

### The remaining blockers, as measured

| Blocker | Status |
|---|---|
| Literal `nan` in a table | Confirmed. `tables/Table4_supraadditivity.tsv` 2 rows, `TableS4_supraadditivity_grid.tsv` 12 rows, and `PAPER.md:206` renders `\| nan \| nan \| nan \|`. Root cause: the `d4` formatter in `08_make_tables.py` has no finite guard, unlike the `D_both` formatter beside it. |
| Figure order | Confirmed 1 → 3 → 4 → 2 (`MANUSCRIPT.md` lines 82, 126, 144, 168). Resolves itself when Fig. 2 leaves. |
| Draft-history language | Confirmed, 8 instances at lines 120, 128 (×3), 132, 136, 146, 160. |
| Absolute path in vendored ODE | Confirmed, `vendor/two_pool_ode.py:69–83`. |
| Deposition DOI | Still `TODO`, `MANUSCRIPT.md:246`. |
| Licence | Still none anywhere. |

### A collision to handle deliberately

Fixing the binding equation means editing `two_pool_ode.py:127–129`, and
`tests/test_manuscript_numbers.py` asserts the vendored copy is **byte-identical**
to upstream. Do not delete that test to make room. Either fix upstream and
re-vendor (preferred — it is your own project, and the byte-identity guarantee
stays meaningful), or keep upstream frozen and convert the test to assert that the
only deviation is the documented binding patch.

---

## The steps

Owner column: **CC** = mechanical, hand to Claude Code. **K** = needs Kiran.

| # | Step | Gates | Owner | Estimate |
|---|---|---|---|---|
| 1 | Codon section → P2 | — | CC | 1 session |
| 2 | Import the bound generators; purge computed inputs from `data/raw/` | blocks 3,4,5 | CC | 1 session |
| 3 | Binding equation → quadratic root | needs 2 | CC | half a day |
| 4 | `C_N` explicit; retire θ | needs 2,3 | **K** | ~1 week |
| 5 | Rerun at both operating points; decide framing | needs 2,3,4 | K + CC | 1 session |
| 6 | Formalism, terminology, blockers | independent of 3–5 | CC | 1 session |
| 7 | Prose pass, then PLOS Computational Biology | needs all | K | — |

### Step 1 — codon section to P2

**Leaves:** scripts 01, 02, 04, 13, 14; Fig. 2; Table 5; Tables S1, S2, S3, S6,
S7; the detectability, power, jackknife and residualization analyses; the ν
validation gate; test classes `AxisPowerAndRobustness` (14 tests),
`TaiAxisIsValidated` (4), and the codon-scoped tests inside
`ManuscriptContainsClaim`, `ArtifactsExist` and `TablesAgreeWithData`
(~30–35 tests total migrate).

**Also delete:** step 07 and the removed-results diagnostics. The archived-analyses
disclosure in `README.md` stays — it is the record — but its tables can go to P2
with the codon work, since the metal-site analysis is a codon-deployment result.

**Stays in P1, and this is the coupling the triage does not spell out:** the μ
*data* is not the μ *analysis*. `06_translation_burden.py` needs
`codon_error_rates_ecoli.tsv` and `global_codon_usage_ecoli.tsv` for Result 1's
6.3 × 10⁻⁴, 613-fold span and 16–18%. Both papers will share those inputs. Copy
them into P2 rather than symlinking — this project has already been burned once by
`structural-criticality/` symlinking a corrupt tAI vector. Sharing an input
dataset raises no duplicate-publication issue; sharing an analysis would.

**Result:** pipeline 16 steps → ~9, both slow permutation steps gone, run time
~12 min → ~2 min, and the three open codon questions leave P1's critical path.

### Step 2 — the generators (the gate)

1. Import `paired_mc.py` and the combinatorial solve as numbered steps, reading
   from `data/raw/` and writing `data/computed/`.
2. Vendor `ecoli_proteome_lengths.tsv` into `data/raw/`; delete the absolute path.
3. Move `arithmetic_results.json`, `two_pool_results.json`,
   `paired_mc_results.json` out of `data/raw/`. They are computed output wearing
   the name of an input. Keep them under `data/reference_upstream/` only as a
   regression target: the reimported generators must reproduce them, and a test
   should assert that.
4. **New test: no file in `data/raw/` is written by any script.** Implement by
   scanning every script for writes whose path resolves under `data/raw/`, plus a
   manifest of raw files with checksums that the suite verifies unchanged.
5. Point `06` at the vendored lengths table, not at `arithmetic_results.json`.

Until this lands, nothing downstream is trustworthy, including anything step 3
or 5 computes.

### Step 3 — binding equation

Replace `c_free()` at `vendor/two_pool_ode.py:127–129` with the quadratic root.
Expected direction, from the triage's own arithmetic: saturation 97.4% → 98.0% at
baseline, 71.4% → 80.7% at θ = 0.9, `C_tot` **more** inert, the throughput-versus-
concentration asymmetry stronger, and the bottom of the headroom range above ×1.9.
The correction runs one way everywhere, so no current conclusion depends on the
error — fix it because it is unphysical and quoted to three figures.

Downstream to regenerate: 09, 11, 12, Fig. 3c, Fig. 4, Tables 3, 4, S4, S5, and
every saturation number in the prose.

### Step 4 — `C_N` (the only new science)

`C_T = C_f + C_U + C_A + C_N`, with `C_N` the commitment to ordinary nascent-chain
folding. This retires θ as a swept parameter and replaces the paper's softest
number with a computed one. Kiran's step.

### Step 5 — both operating points

Per Part 0: vary condition as a triple, not the error rate alone. Confirm the
Landerer measurement's growth condition first. **Do not tune parameters** until
the model survives or fails on its own calibration.

### Step 6 — formalism, terminology, blockers

- Vector capacity: `D(x) ≤ C`, `H = min_i C_i/D_i`, replacing the additive scalar
  sum. One subsection; no quantitative result depends on it, since the numbers come
  from the two-pool model.
- "any reduced model" → "this minimal model, when saturable rescue is coupled to
  sufficiently strong superlinear amplification".
- Derive or label `Φ(P)`: state what the amplification physically represents, or
  mark it as a posited closure.
- Terminology: mechanistically distinct (not independent) bounds; simulation
  intervals (not 95% CIs); model-implied viability thresholds (not maximum
  tolerable rates).
- `nan` fix + guard. **Implementation note, deviating slightly from the literal
  instruction:** forbid `nan`/`inf`/empty in *rendered* tables (TABLES.md, PAPER.md)
  absolutely; in the full-precision TSVs, permit a non-finite value only in a row
  flagged `collapsed_both`, and require the flag. Blanket-forbidding it in the TSVs
  would force a number into a genuinely undefined cell — which is exactly the error
  that produced spurious negative interactions earlier in this project.
- Figure renumbering (free once Fig. 2 leaves).
- `PAPER.docx` / `PAPER.pdf` untracked — **already done** (commit `4fad6dc`).
- Strip the 8 draft-history phrases, including "The defensible statement is" and
  the falsifiability aside.
- References: add Stikeleather 2026 (NAR 54:13, `10.1093/nar/gkag674`), Sun &
  Zhang 2022 (`10.1126/sciadv.abl9812`), the mistranslating-tRNA branch study
  (NAR 2025, `gkaf428`), and co-translational chaperone binding
  (Nat Commun 2025, `s41467-025-59067-9`).
- Check the author name against your other publications; keep `Classification`
  and `Significance statement` only if PLOS wants them (it does not).
- Deposition DOI; licence.

---

## P2 already exists, and the migration unblocks it

Decision resolved by the record, not by my inference. `/storage/kiran-stuff/PORTFOLIO_PLAN.md`
and `/storage/kiran-stuff/dashboard/papers/P2_error_matched.md` define P2 as
**"Synonymous codons are error-matched and supply-diverse"**, workspace
`/storage/kiran-stuff/codon-deployment/`, draft at `codon-deployment/JME/MANUSCRIPT.md`
(2026-04-16), single-author, target MBE or GBE. **That is where the codon section
goes.** Not `codon_project/`, and not a new sibling.

P2 has been halted since April 2026: `codon-deployment/DATA_AUDIT.md` reads
**"STOP — Issues 1 and 2 require resolution before figures"**. Both are things the
migrating material already fixes:

| P2's blocking issue | What arrives with the codon section |
|---|---|
| **Issue 1** — μ values wrong by 2–4 orders of magnitude; two of four within-pair orderings reversed. Audit says they "need replacement from Landerer 2024 Data_S2". | Exactly that: `data/raw/codon_error_rates_ecoli.tsv` and `data/computed/codon_axes.tsv`, extracted from Data_S2 and asserted by tests, plus the median-statistic sensitivity. |
| **Issue 2** — metal-enrichment p-values not significant; `metal_sites_ecoli_with_codons.csv` has a "codon mapping bug (~45% wrong codons)". | `07_removed_results.py` settles it: 0 of 4 survive an expression-matched within-gene background, His 1.33 → 1.07. And its independent diagnostic found ~40% of annotated sites mispositioned — corroborating the audit's ~45% from a different direction. |
| **Issue 3** — cross-species tAI files fake (resolved by recomputing from GtRNAdb). | The validation gate itself (`01_validate_tai.py`), which any future species vector must pass, plus the regression guard on the corrupt vector. |

**Two of P2's three headline results in the portfolio plan are dead.**
`PORTFOLIO_PLAN.md:13` lists "μ-clustering z = −3.59; metal-site OR 1.28–1.40;
cross-species ρ 0.63–0.92". The second and third were refuted in this repository —
the metal-site effect is an expression confound and the cross-species correlation
was a shared-μ artifact with 44 of 60 supply values bit-identical between species.
Anyone resuming P2 from the dashboard would rebuild a paper on two refuted results.
The archived-analyses material must therefore migrate **with** the codon section:
those are P2's negative results, not P1's, and P2 is where they carry weight.

Only the μ-clustering result survives, and it is now much stronger than z = −3.59:
the depth residualization, the non-uniform ν benchmark and the 400-replicate power
curve all arrive with it. The three questions that were blocking P1 become P2's
principal robustness section.

**Authorship consequence to settle before the move.** P1 is Rebbeck, Paudyal &
Boggavarapu; P2's draft is single-author. Moving analysis out of a three-author
paper into a single-author one is an authorship decision, not a filing decision.

## Decisions needed

1. **Upstream ODE:** fix `proteostasis-P1/two_pool_ode.py` and re-vendor, or fork
   with a documented single deviation?
2. **P1's venue** — see the correction below.
3. **P2 authorship**, given that material moves between papers.
4. Licence for the deposition, and reserve the Zenodo DOI (blocks submission, not
   the repair).

## Correction: P1's venue was already decided, and I said it was not

`PORTFOLIO_PLAN.md:12` assigns **P1 to "MBE or J Mol Biol"**, and line 33 sequences
it as "P1 to MBE or J Mol Biol (independent timeline)". Line 23 rules PNAS out
explicitly: "no NAS sponsor is available". I reported that no venue was on record
and recommended PLOS Computational Biology, having searched only inside
`proteostasis_law/`; the portfolio plan is one directory above, in
`/storage/kiran-stuff/`.

The substantive point survives the correction, and the migration sharpens it. MBE
is *Molecular Biology and Evolution*, and it was a reasonable assignment while P1
still contained the μ/ν code-organization material. After step 1 there is no
evolutionary argument left in P1 — it becomes a dynamical-systems constraint paper
about extant translation, with an explicit disclaimer that it makes no code-origin
claim. J Mol Biol fits that; PLOS Computational Biology fits it better still,
because the distinguishing prediction is untested and that journal takes
model-plus-reanalysis as complete. MBE now fits **P2**, which is where the
evolutionary content is going, and which the portfolio also aims at MBE — two
submissions to the same journal from the same program is worth avoiding on
reviewer-pool grounds, by the plan's own reasoning at line 23.

Recommendation: P1 → J Mol Biol or PLOS Comp Biol, P2 → MBE. Kiran's call.

---

## Parked — surfaced during planning, out of scope for this run

- Landerer's growth condition is unconfirmed; needed before step 5 can call
  6.3 × 10⁻⁴ an exponential-phase number.
- The scalar `f_codon` coarseness (review §6) is a real limitation but a modelling
  extension, not a repair; a sentence in Limitations, not a rewrite.
- `structural-criticality/` still has tAI-axis results resting on the corrupt
  vector.
