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

**Q:** Go (re-anchor the chaperone parameters).

**A:** Read `proteostasis-P1/LITERATURE_ANCHORS.md` first, and it overturned my own
diagnosis. C_tot (30-80 uM, baseline 50; Lorimer 1996 GroEL, DnaK levels), K_d
(0.06-2 uM, baseline 1; Pierpaoli 1997 EMBO J) and k_obs_max are all properly
sourced and used in range. So the 97.9%-saturated folding arm is NOT a parameter
error, and my earlier claim that the parameterization "contradicts the paper's own
capacity evidence" was wrong. Retracted it in the manuscript.

The real diagnosis is structural: the rescue term hands the ENTIRE chaperone pool
to the damaged-protein pool because the model does not represent nascent-chain
folding at all. So C_tot in this model is chaperone *available to the damaged
pool*, and using the total silently assumes availability is 100%.

Made that assumption an explicit parameter instead of re-anchoring by invention:
C_avail = C_tot(1 - theta), swept with C_tot and K_d over their documented ranges
(`scripts/12`, Table 7). At the usage-weighted mu:

- theta = 0 (the implicit assumption): saturation 0.974, headroom x24.8, margin 1.39
- theta >= 0.90: margin reaches the supraadditivity onset (1.19), headroom x15.7
- theta >= 0.98: folding arm unsaturates (c_free 0.53 uM < K_d), headroom x4.6
- full documented grid: x1.9 to x25.3

theta is NOT measured and is not estimated. What this buys is a falsifiable
threshold: measure chaperone occupancy by nascent-chain folding in exponentially
growing E. coli, and if theta >= 0.90 then E. coli already operates in the regime
where burden and capacity perturbations compound. That collapses the headroom range
to a point.

Result 3, Limitations, Fig 3 legend and Methods rewritten; Methods now carries the
parameter provenance explicitly. Tests 53 -> 60, including one asserting the
retracted overstatement does not return. One of my own earlier tests failed
correctly during this pass -- it asserted the Limitations section still contained
the phrase I had just retracted -- and was updated to check the corrected framing.
Mutation-tested (reverting to the parameter criticism, presenting theta as
measured, dropping the range from Limitations, or dropping the pinning measurement
each fail). Verified by wiping all artifacts and rebuilding.

**Q:** Review: one blocking defect (Result 5 evaluated at the point Result 3
rejects), three consistency errors, two unaddressed objections.

**A:** All confirmed and fixed. The blocking defect was mine and the review's
numbers reproduced exactly.

**Blocking.** `scripts/09` hardcoded `OBSERVED_F = 1e-4` — the window bottom, the
evaluation point Result 3 spends a page rejecting. So Result 5 was anchored at
margin 2.20 while the paper's own margin is 1.39. Re-anchored (the point is now
READ from translation_burden.json; 1e-4 kept as a labelled comparison row):

- interaction, 3x error / /3 rescue: +0.196% -> **+7.44%** (a 38x understatement)
- max non-collapsed: +13.0% at 5x/20x -> **+9.59% at 2x/5x**
- single-viable, jointly lethal: 5/36 -> **12/36**, mildest 1.5x error with /10
- saturation triple: 97.9%/47.5/0.052 uM -> **97.4%/37.6/0.331 uM** (6 places),
  now matching Table 7's theta=0 row; "~900-fold excess" -> ~113-fold

This inverted the experimental conclusion in the paper's favour. Prediction 1 no
longer requires a sensitized background: synthetic lethality is present at
wild-type margin. Worst part: the test suite was ENFORCING the inconsistency by
asserting the strings "+0.2% of the additive expectation" and "97.9% saturated".
Replaced with a test that the evaluation point matches Result 3's, plus one that
the anchor is not hardcoded again.

**Objection 1 (nu power).** Correct that nu is coarse: 21/59 distinct values,
largest tied group 13, vs 59/59 for log mu. But the test still has the power:
nu rejects at >=35% tightening (19.6% below null), and mu's observed clustering
sits 37.4% below its null — a 17.8 pp margin. So the contrast is not a resolution
artifact. Boundary conceded: under a subset model (effect in half the amino acids)
power reaches only 0.75 and never 0.80, so the claim is narrowed to "does not
operate on nu uniformly across amino acids at comparable magnitude".

**Objection 2 (mu leverage / detectability).** Leverage fails: leave-one-codon-out
jackknife is significant in 59/59 deletions, and dropping CCC gives z = -3.67
(marginally stronger); span reported both ways, 613-fold and 286-fold without CCC.
Detectability is real and damaging: sampling depth carries **eta^2 = 0.560**
between amino acids against 0.556 for log mu, and Spearman(depth, log mu) = -0.37
(p = 0.004). Survives dropping the thinnest quartile (z = -2.85) but cannot be
separated from the biology with these data. Now stated in Result 4 and Limitations.

**Consistency.** Abstract x4.6-x25 -> x1.9-x25 (the theta grid); Table 6 retired to
S9 (its anchoring axis duplicated theta and reached its low end at C_tot = 1-2 uM,
outside the documented 30-80; verified the two cells are identical to 8e-8);
README "two orders" fixed; Table 3's 2D rows declared uninterpreted by
construction; nu p-values labelled descriptive with the ~0.68 hypothesised-direction
value stated; "median factor of 3.1" distinguished from the median ratio 3.00.

Tests 60 -> 73. Verified by wiping every artifact and rebuilding (15 steps).

## 2026-07-30

**Q:** Save the pasted v2 restructuring as `manuscript_v2_draft.md`, then add
figures and make a complete manuscript with figures, tables, and references.

**A:** Saved the draft verbatim at `envelope-paper/manuscript/manuscript_v2_draft.md`
and rebuilt the paper around it.

**Table renumbering.** The v2 inventory drops the two excluded analyses as paper
tables and retires the anchoring grid, so: Table 3 = chaperone availability (was 7),
Table 4 = supraadditivity (was 5), Table 5 = axis tests (was 3), S4-S7 shift down
one. `08_make_tables.py` was rewritten with CONTENT-named functions (a number now
appears once, in the output map) so a renumbering cannot leave a `table7()`
emitting Table 3; it also deletes the ten files from the old numbering on every
run, since a stale `Table7_*.tsv` in `tables/` is indistinguishable from a live one.
The two excluded analyses are now `Excluded_*.tsv`.

**Two stale captions found and fixed.** `TABLES.md` still described the knob table
with the retracted window-bottom triple (97.9% saturated, 47.5 uM, 0.052 uM) and
the effect grid as "f = 1e-4" — both corrected in the manuscript a commit earlier.
Every caption number is now interpolated from `data/computed/`, and a test asserts
the retracted strings never come back. Same failure mode as the blocking defect:
correcting a premise and leaving a consumer stale.

**Assembly instead of hand-maintenance.** `MANUSCRIPT.md` carries the prose,
embedded figures and one `<!-- TABLE:Table N -->` placeholder per main table;
`15_build_paper.py` fills them from `TABLES.md` and emits `PAPER.md`, `PAPER.html`
(figures base64-embedded), `PAPER.pdf` (17 pp, printed by headless Chrome — no
unicode-capable TeX engine is installed) and `PAPER.docx`. A test fails if a
markdown table row ever appears in the prose.

**Restored what v2 had dropped.** The v2 rewrite lost the excluded-analyses
disclosure, the Limitations section, and with them the only citations of refs 25
and 26 (which would have shipped uncited). Both sections are back: *Two analyses
excluded during verification* (the within-gene ORs, ~60% site concordance, 44 of 60
identical tAI values) and *Limitations* (E. coli-only mu, tAI as a proxy, the
theta range, the detectability confound).

**Also fixed two more stale claims** in prose I had not swept: `envelope-paper/README.md`
said "two orders of magnitude inside" and, under Known gaps, still carried the
retracted "only ~0.2% above additive ... the experiment must first compress the
viability margin". Root `CLAUDE.md` said "two orders" as well.

Tests 73 -> 89, and mutation-tested: reverting the interaction to the window-bottom
value, dropping a figure, hand-typing a table row, citing a nonexistent table,
uncrediting refs 25/26, reinstating the saturation triple, dropping the theta
concession, stating the headroom as a point in the abstract, and reverting a
caption to a typed value all fail. One escape found and closed (a bracket-style
citation slipping past the new parenthetical parser). Verified by deleting every
generated artifact and rebuilding all 16 steps.

**Q (same day, on review):** The reference argument was backwards — an uncited
reference is repointed or deleted, never propped up by restoring the text it sat
in. Cut both restored sections. Check first whether Limitations duplicates the
inline treatments. Put the excluded-analyses disclosure in the repository README,
not the paper. And confirm the Chrome-printed PDF is print resolution and that no
metal-site table came back.

**A:** Audited the restored Limitations section item by item: seven of its nine
items duplicated an inline treatment word for word (detectability in the codon
result, theta and the headroom range in the margin result, the nu power narrowing
in the supply paragraph, the unfitted ODE in Theory, the order-of-magnitude bound
in the bounds result, the untested prediction in two places). Two carried
something the inline text did not — the single-species scope of the codon result
and the absence of archaeal/mammalian evidence — and those moved to the claims
they bound. Then deleted the section. Deleted the excluded-analyses section too
and rewrote the disclosure as *Archived analyses* in `envelope-paper/README.md`,
where a reader who finds `scripts/07_removed_results.py` will be, with the
recomputed ORs, the ~60% concordance and the 44-of-60 tAI identity.

Refs 25 and 26 repointed to the folding-burden paragraph beside Walsh and Buhr,
carrying the claim that B_fold is a burden component output does not predict —
a stronger placement than the generic tAI caveat they had before.

Build checks. No collision: the excluded work has no numbered table, only
`Excluded_*.tsv`. Chrome had embedded the full 300 dpi rasters (2132x783 for a
7.09 in panel, ~285 dpi at printed width), so resolution was adequate — but
figstyle now also writes SVG and 15_build_paper.py substitutes it for the typeset
output, so the figures are vector. One 26x536 raster remains, the Fig 4 colorbar
ramp, which carries no text. PDF 1.1 MB -> 0.61 MB, 16 pp.

Tests 89 -> 94. The Limitations-scoped assertions were retargeted to the sections
that now carry each claim rather than deleted, plus NoStandaloneCaveats (neither
section may return, no trace of the withdrawn work in the paper) and a test that
every inline caveat survived the deletion — so removing a caveat inline cannot
pass just because the section that duplicated it is gone. Seven mutations
confirmed failing, including moving 25/26 back to a Methods caveat.

**Verification rebuild (same day).** The final 16-step rebuild from `data/raw/`
passed — all steps ok, 94 tests ok — but `git status` afterwards showed twelve
modified artifacts, which is exactly the noise that hides a real change. Two
causes, one substantive:

- matplotlib writes `<dc:date>` into svg/pdf and derives svg element ids from
  object addresses, so every rebuild rewrote every figure. Fixed with
  `svg.hashsalt` and `metadata={"Date": None}` / `{"CreationDate": None}`.
- **Fig 2c's stripplot jitter was unseeded.** seaborn draws jitter from the numpy
  global RNG, so the same data produced a visibly different panel on every build.
  Fixed with `figstyle.JITTER_SEED = 42`.

Everything except `PAPER.pdf` and `PAPER.docx` now reproduces byte-for-byte
(Chrome and pandoc stamp a creation date into the container). Verified by building
twice and diffing md5s. Tests 94 -> 96.

**Q (review round 3):** Untrack the two nondeterministic renderer outputs rather
than documenting an exception. Then the three open items that can still move a
claim: run the depth-residualized permutation test instead of conceding the
detectability confound; impose mu's observed per-amino-acid pattern on nu rather
than a uniform shrinkage; raise the 40-replicate power estimate. Also: does the
three-z-values discrepancy carry a caption note, and is the parent repo still
auto-committing?

**A:** All three analyses run; all three moved claims.

**Detectability — the concession was not earned, and is now withdrawn.**
Regressing log mu on log sampling depth (slope -0.93, p = 6.0e-4) removes the
18.8% of variance depth explains; on the residuals the clustering persists,
Delta = 0.777 vs null 1.130, z = -2.97, p = 0.0026, with 80% of the
between-amino-acid eta^2 retained (0.446 vs 0.556). Residualizing is licensed by
the sign of the correlation, exactly as argued: negative means thin-sampling
inflation, not error-driven detection, so depth is a nuisance regressor rather
than a mediator. "Cannot be separated from detectability" is retracted; the
remaining caveat is one sentence about a detectability component depth does not
track.

**nu power benchmark.** mu's per-amino-acid tightness runs 0.06 to 2.03
(median 0.52) with 2 of 18 amino acids wider than chance -- strongly non-uniform.
Imposing that exact pattern on nu: z = -2.75, p = 0.0034, DETECTED. Calibration
(same pattern on a structureless axis) gives z = -3.39 +- 0.63 against mu's own
-3.54, so the transfer is faithful and nu's coarse resolution attenuates without
erasing. The claim strengthens from "a uniform effect of mu's magnitude would have
been seen" to "mu's effect in the form it actually takes would have been seen".

**Power at 400 replicates (parallel, Wilson intervals).** 0.51 [0.46, 0.56] at a
60% tightening and 0.83 [0.79, 0.86] at 80%. The old wording -- "does not reach
0.80 anywhere on our grid" -- was not merely imprecise at 40 reps, it was WRONG:
the point estimate at the extreme grid point exceeds 0.80. Rewritten.

**Three z values.** -3.56 / -3.55 / -3.54 are the same statistic at 10,000 /
4,000 / 2,000 permutations. Now stated in Methods, with the identity that makes it
checkable: the two analyses that use 10,000 agree to -3.5643673490070. A test
asserts that agreement across scripts.

**PAPER.pdf/docx untracked and gitignored.** Documenting an exception installs a
standing carve-out into the one check whose value is that it needs no judgment.
Rebuild-and-diff is now exception-free.

**Auto-commit cron: still live.** `*/30 * * * * /storage/kiran-stuff/git-auto-sync.sh`,
which runs `git add -A`, commits and pushes. It fired twice during this session;
the 22:30 commit's entire content was the two nondeterministic renderer outputs.
Flagged for a decision -- untracking those removes the empty-commit class but not
the risk of publishing a half-edited state.

Tests 96 -> 103.
