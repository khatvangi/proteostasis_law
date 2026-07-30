# HISTORY

## 2026-07-29

**Q:** Check the combined-paper work — science, code, results. It is almost ready.

**A:** Reviewed `combined-paper/P1_combined_capacity_codon_DRAFT.md` against the
upstream sources (`codon-deployment/`, `proteostasis-P1/`). Not submission-ready.
One blocking data bug and several serious issues:

1. **BLOCKER — the ν (tAI) axis is corrupt.** `codon-deployment/data/computed/ecoli_tai_ws.tsv`
   is 45/61 codons bit-identical to the yeast file; E. coli CTG gets tAI 0.059
   (should be high — dominant Leu codon), AAG gets 1.00 while AAA gets 0.43.
   `scripts/01_compute_species_tai.py` neither reproduces those files nor
   implements dos Reis 2004 correctly (its wobble indices are permuted;
   r = −0.19 vs the canonical formula). So there is no valid provenance for ν.
2. **Result 4 is confounded twice.** Δ_A in `delta_a_by_species.tsv` is the 2D
   (μ, ν) distance and uses *E. coli μ for all three species* (verified to
   0.0000 max diff), so the cross-species correlation shares an identical
   coordinate. Combined with the 45/61 shared ν values, ρ(E. coli, yeast) = 0.92
   is largely artifactual. ν-only + independent tAI gives ρ = 0.58 / 0.46 (ns) / 0.73.
3. **Result 3 (metal sites) may be an expression confound.** The repo's own
   `metal_enrichment_legacy.tsv` runs the same test against a narrower
   (metalloprotein-gene-sized) background and gets OR 1.21–1.34, all p > 0.08,
   with His at OR 0.99. The published version uses a genome-wide background;
   only then is it significant. Also: ~40% of proteins in
   `metal_sites_ecoli_with_codons.csv` have completely mis-mapped codons
   (bimodal per-protein concordance, 46% overall), and N does not reconcile
   (paper 1,276 vs 1,109 concordant vs 2,397 all vs 1,003 in metal_raw_counts.csv).
4. **Result 6 quotes superseded numbers.** ODE MC median 1.7e-2 is the
   `backup_uniformN` run; current is 1.99e-2 (paired 2.59e-2). "Paired
   P(arith < ODE) = 0.65" is the statistic `PAIRED_MC_TASK.md` explicitly calls
   "not a valid paired statistic"; the corrected value is 0.768.
5. **20 of 36 references are never cited**, and the Discussion's error-minimization
   citation "[6,7]" points to a chaperone protocol and a synonymous-edit paper
   (should be 32/33) — in-text numbering was not remapped when the two source
   reference lists were merged.
6. **No code.** `combined-paper/` is prose + figure images only. No script
   anywhere generates Fig1_envelope, Fig2_axis_asymmetry, or Fig5_bound_convergence,
   and the whole statistics layer in `codon-deployment/` is precomputed TSVs with
   no generating script. The 15 passing tests cover the *old* `proteostasis-paper/`
   and check zero numbers.

What holds up: the μ axis (μ clustering z = −3.59 reproduced as −3.54; robust to
mean→median, z = −2.57; 613-fold span confirmed), the metal-site enrichment
*direction* (robust to all three tAI variants), and the ODE headroom result
(×158 / ×11,091 confirmed in `two_pool_results.json`).

**Q:** Make this work count — the paper and the differential equation got derailed
and demoted; rescue it.

**A:** Rebuilt the paper around its spine as `envelope-paper/` (now canonical).
Chose "envelope + μ-axis" scope after establishing that Result 3 dies under an
expression-matched background (0 of 4 significant, His OR 1.33 → 1.07) and that
`structural-criticality/` had already failed to replicate it across nine other
criticality definitions.

Delivered:
- `scripts/01_validate_tai.py` — validates the ν axis against a translational-
  selection test and **aborts** if it fails. The corrupt vector fails (12/18
  amino acids, p = 0.12); the reference passes (15/18, p = 0.0038).
- `scripts/02_axis_structure.py` — regenerates Δ_A, both permutation nulls, and a
  variance decomposition, from raw inputs. With validated ν: μ z = −3.56
  (p = 0.0010), ν z = −0.06 (p = 0.45), η² = 0.556. μ survives median
  summarization (z = −2.59).
- `scripts/03/04/05` — Fig 1 (saddle-node regimes), Fig 2 (axis structure),
  Fig 3 (paired bounds + headroom). Fig 1 and the bounds figure previously had no
  generating script anywhere.
- `scripts/06`, `scripts/07` — burden magnitude (16–18%, not the drafted 20–25%)
  and reproduction of both removed results.
- `manuscript/MANUSCRIPT.md` — Results 3 and 4 removed and reported as negatives;
  Result 6 repointed at current JSONs and the valid paired statistic (0.768, not
  the retracted 0.654); 27 references, all cited, all with identifiers.
- `tests/test_manuscript_numbers.py` — 31 tests, mutation-tested: altering a
  number, dropping a citation, softening a negative result, reinstating the
  retracted statistic, or substituting the corrupt tAI each make it fail.

Corrections to my own earlier review: the tAI/yeast identity is **44 of 60**, not
45 of 61 (my quick check counted absent codons as matches), and the metal-site
`codon` column is faithful to the CDS 98.3% of the time — the broken step is
`uniprot_pos`, not the codon lookup.

Also fixed a science error carried in the old Fig 1 legend: the transition is
loss of the low-burden state, not bistability — `+χx²` dominates at large x, so
there is no second stable attractor.

Housekeeping: `envelope-paper/` named canonical in CLAUDE.md and README;
`combined-paper/SUPERSEDED.md` records the four blockers;
`README_TAI_IS_CORRUPT.md` written in codon-deployment and a warning appended to
`structural-criticality/CLAUDE.md`, whose `in_disulfide` tAI result needs
recomputation.

**Q:** Generate the figures, tables etc.

**A:** Figures already existed (Fig 1-3, PNG+PDF). Added the table set and a driver:
- `scripts/08_make_tables.py` → `tables/`: Tables 1-4 and S1-S4 as full-precision
  TSV plus `TABLES.md` with typeset values formatted identically to the
  manuscript, so drift between table and prose is detectable.
- Table 3 fixes two presentation errors found while building it: the ν rows were
  duplicated across the mean/median μ blocks (ν doesn't depend on the μ statistic,
  so one test was displaying as two), and non-significant rows were labelled
  "clustered"/"spread" on the strength of their sign alone — now "no signal".
- `scripts/run_all.py` — one-command rebuild of everything, fails closed if the ν
  axis validation rejects (nothing downstream runs).
- Manuscript now cites Table 1-4 and S1-S4 at the right places.
- Test suite 31 → 41. New tests assert tables match `data/computed/`, the
  manuscript matches the tables, ν rows aren't duplicated, null results aren't
  labelled as signals, and Tables 4/S4 stay marked "excluded from the paper".
  Mutation-tested: table drift, relabelling a null, making a within-gene p
  significant, and dropping the "excluded" marking each fail.

**Q:** What next?

**A:** Recommended and ran the supraadditivity test — the paper's one
distinguishing prediction, previously asserted only verbally. 2x2 factorial on the
vendored two-pool ODE (error rate up vs rescue throughput down), readout
log10(min(P_dagger/P*, A_max/A*)).

Result: the interaction is supraadditive everywhere it is defined, never
subadditive — so the framework is coherent on its own terms. But at wild-type
margin it is **+0.2% of additive** (undetectable). It grows to +4.2% at margin
1.50, and at margin 1.19 (x16 headroom) 3-fold perturbations that are each
survivable become **jointly lethal** (338 of 676 grid combinations). This merges
predictions 1 and 2: the experiment must compress the margin first, and then looks
for synthetic lethality, not a small additivity deviation. Predictions rewritten
accordingly; the paper now says the prediction is harder to test than earlier
drafts implied.

Two self-corrections during the analysis:
- My first capacity knob was `C_tot`, which turned out nearly inert — the folding
  arm is 97.9% saturated at the operating point (47.5 uM free chaperone vs
  0.052 uM misfolded protein), so a 3-fold pool cut leaves v_fold at 96%.
  Switched to `k_obs_max` (throughput) and reported both. The saturation is now a
  named Limitation: it contradicts the capacity evidence cited in Result 1
  (refs 4, 21), so the model likely understates how tightly C_buffer binds.
- My first collapse handling substituted D_both = baseline margin as a "lower
  bound", which produced spurious NEGATIVE interactions (-17%) for cases that are
  in fact the strongest supraadditivity. Collapsed cells now carry no numeric
  interaction at all, and a test enforces that.

Added: Result 5, Fig 4, Table 5, Tables S5/S6, `scripts/09`, `scripts/10`,
`scripts/vendor/two_pool_ode.py` (verbatim upstream, with a test asserting it has
not drifted), and `scripts/run_all.py` extended. Tests 41 -> 50; mutation-tested
(overstating the wild-type effect, dropping the model-only caveat, removing the
saturation limitation, numbering a collapsed cell, or introducing subadditivity
each fail). Verified by deleting every generated artifact and rebuilding from raw.

**Q:** Commit, and do things to strengthen the work.

**A:** Committed as `be5fa1d` on master (a cron auto-syncs every 30 min, so
branching would have fought the automation; left the unrelated proteostasis-paper
submodule change alone).

Strengthening found a real problem with the paper's headline number. The ×158
headroom was obtained by evaluating the two-pool model at f = 1e-4 — the BOTTOM
of the quoted observed window — while the paper derives a usage-weighted mean of
6.33e-4 from its own mu data. Same class of error as quoting the arithmetic bound
at its deterministic point rather than its MC median.

Corrected (`scripts/11_headroom_sensitivity.py`, Table 6, rebuilt Fig 3c):

| evaluated at | as published | across six chaperone anchorings |
|---|---|---|
| 1e-4 (window bottom) | ×158 | ×40–×158 |
| **6.33e-4 (usage-weighted mu)** | **×25** | **×4.6–×25** |
| 1.1e-3 (unweighted mean) | ×14 | ×1.4–×14 |
| 1e-3 (window top) | ×16 | ×2.1–×16 |

So "roughly two orders of magnitude inside" holds only at the most favourable
corner. The defensible claim is roughly ONE order — ×25 at the internally
consistent evaluation point.

This weakens the margin claim but substantially strengthens the paper: at margin
1.39 log10, E. coli sits only +0.20 above the 1.19 margin at which 3-fold error
and 3-fold rescue perturbations become jointly lethal (Result 5). The
distinguishing prediction is reachable with a modest sensitization, not remote.
Abstract, significance statement, Result 3, Result 5, prediction 2, Fig 3c legend
and Methods all rewritten.

Also strengthened the test suite after mutation testing showed two escapes:
asserting a value appears SOMEWHERE passes even when the abstract and results
disagree. Added `section()` scoping so headline claims are checked against the
section that makes them, plus a test that Result 3's inline table matches
`headroom_sensitivity.tsv` cell by cell. That last test caught a bug in itself —
"10⁻⁴" is a substring of "6.3 × 10⁻⁴", so it had been matching the wrong row.

Tests 50 -> 53. Verified by deleting every generated artifact and rebuilding.
