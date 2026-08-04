# History

## 2026-07-31

**Q:** Host-level operational cleanup and a single clean Phase 1 submission on boron:
find and gracefully kill any live `run_experiment_[a-d].py` processes from the two
duplicate launch attempts (recorded PIDs 743-746 and 749-752), quarantine the existing
`results/phase1/{A,B,C,D}`, confirm the repaired test suite passed, launch exactly one
clean run of A-D, verify at 20 s, write a manifest plus JOBS.md/STATUS.md, and commit.
Do not touch scientific code, tests, configs, theory, or manuscript.

**A:** No live experiment processes existed on the host — nothing was killed. The
recorded PIDs came from a sandboxed session with its own PID namespace; on boron,
743/749/750 have since been reused by kernel worker threads (`kworker/*-kblockd`,
PPid 2), so killing by recorded PID would have signalled unrelated kernel objects.
Identity was verified via `/proc/<pid>/status` before any action. The empty A-D
directories were moved (not deleted) to
`results/phase1/quarantine_preclean_20260731T162755-0500/`.

The test gate passed (exit code 0, `69 passed`), so one instance each of A-D was
launched as `systemd-run --user` transient units into
`results/phase1/run_20260731T162946-0500/`, with thread limits verified inside each
running process. At +20 s all four were active with empty stderr and exactly 36 live
processes (4 parents + 32 workers), matching the configs. Experiment A has since
finished successfully (11/11 internal-consistency checks); B, C, D still running.

Committed as `a17dfaf` on `master`. Not pushed — the repository has no remote
configured at all. Results and ops logs are excluded by the `results/` gitignore rule,
so the manifest was written but not committed.

**Note for future sessions:** a second agent session was concurrently active in this
repo during this work. Re-check for duplicate launches before starting any new run.

**Q:** Complete the unfinished operational submission for Phase 2A: run the Phase 2
and full pytest suites and fix only real defects, run a smoke matched benchmark,
write `theory/MATCHED_IMPLEMENTATION_PROTOCOL.md` / `PHASE2_JOBS.md` / STATUS
updates, create a timestamped matched run root, launch the full boron n=2000
matched benchmark as a persistent systemd unit, submit the free-limit counterpart
on nitrogen via Slurm, write a submission manifest, and commit. Do not touch the
c02 audit or Phase 1 D.

**A:** Suites passed with no pre-existing failures: `tests/phase2` 76 passed,
full `tests` 145 passed. T0 was rerun and reproduced its headline exactly
(RHS 4.800103e-06, Jacobian 8.263311e-06 at epsilon 1e-6; slopes 1.0023 / 1.0028),
10/10 checks, exit 0.

The smoke run exposed one genuine defect. Every cell summary recorded only
`tsv_sha256`, a byte hash of a file whose rows each carry a wall-clock `seconds`
column, so that hash can never match between boron's free arm and nitrogen's free
arm even when every computed value is identical — which is exactly the cross-host
comparison the benchmark exists to make. Diagnosed by column-wise diff across
three replicates: `seconds` was the *only* differing column, including between
`workers=1` and `workers=4`. Fix was additive — a `payload_sha256` over the result
columns plus `payload_excludes` — with four new tests, one of which asserts that
`seconds` is the only non-deterministic column so the exclusion set cannot drift
silently. No scientific criterion, T0 threshold, or mapping was touched. Suites
after the fix: 80 and 149 passed.

Smoke then passed on all criteria: schema equal to `ROW_FIELDS`, `payload_sha256`
identical across a repeat run and across worker counts, zero numerical failures,
and 8/8 exact label *and* admissibility agreement between the boron and free arms
at epsilon 1e-6 under both protocols.

`theory/MATCHED_IMPLEMENTATION_PROTOCOL.md` did not exist even though `mapping.py`
and `t0_equivalence.py` cite its sections 3 and 5 by number — those were dangling
references to the derivation licensing the tolerances. It was written with the
term-by-term flux identity verified against `model.py`, the a priori tolerance
derivations, and explicit statements that nitrogen is the epsilon = 0 limit and
that percentages are not comparable outside the matched benchmark.

Boron runs the full 28-cell factorial as unit
`phase2-matched-boron-20260731T180846.service` (MainPID 1446897, 16 workers, one
thread each, stderr empty, exactly one instance). Nitrogen runs the 14-cell
free-limit counterpart as Slurm job 4, RUNNING at 1 CPU / 4 GB. Both hosts were
proven to regenerate a bit-identical n=2000 seed-20260731 sample matrix under
different numpy builds (2.2.6 vs 1.26.4) — without that, no cross-host cell could
be compared. Run root `results/phase2/matched_20260731T175912-0500/`, manifest
`SUBMISSION_MANIFEST.json`. No result was interpreted; both jobs were still
running at commit time.

**Note:** the reported audit parent PID 1361326 with 24 workers had no matching
process on boron; its `c02` stage had already written `done in 1132.4s` before this
session began. Its directory was left untouched and is gitignored. Phase 1 D
(PID 1269455) was still running and was not signalled.

**Q:** Operational continuation of the "already-running" Phase 2B multistability audit:
verify the live c02 process (parent PID 1361326, `--workers 24 --n-controls 200`, started
17:35:30), create a launcher + manifest under the gitignored audit root only, run it as a
persistent `systemd-run --user` unit that waits for c02 to exit, gates on completeness,
then runs c03/c04/c05 and d01. No tracked edits, no commits, no interference with the
matched-benchmark task.

**A:** The stated premise was wrong on two points, both verified before acting. PID 1361326
does not exist on boron (`/proc/1361326` absent), and a full `/proc` sweep matching on the
string `c02_dense_multistability` in the command line found no live process at all — the
same PID-namespace mismatch recorded earlier today. More importantly, **c02 had already
finished**: its log reads `304 draws (68 candidate / 36 zero_stable / 200 control), 24
workers` and `done in 1132.4s`, and 17:35:32 + 1132.4 s = 17:54:24, exactly the c02 output
mtime. So the existing JSON/TSV are not a smoke run — they are the full run's own output
(304 records, 0 errors, 608 points).

That inverted the `smoke_pre_full/` step: moving those files would have destroyed the
completed result and broken c03/c04, which both open `c02_dense_multistability.json` from
the audit root (`c03_basins.py:238`, `c04_taxonomy.py:164`). The continuation therefore
took the documented "unsafe -> leave in place, validate mtime and record count" branch,
with the decision and its reasoning recorded in the manifest. The move logic is still
implemented and still checks live open FDs; it simply did not fire.

The completeness gate was not weakened. All 15 checks pass, with record and error counts
recomputed from the `records` array rather than read out of c02's own summary line.

Launched as `proteostasis-phase2-cont-20260731T180145-0500.service`, main PID 1424672,
running `ops/run_continuation.sh -> ops/continuation.py` under
`results/phase2/audit_20260731T172711-0500/ops/`. At +105 s: active/running, 25 c03
processes (1 parent + 24 workers), 67 (draw, evaluation point) pairs in flight, empty
stderr, exactly one such unit. `git status` byte-identical to session start; nothing
tracked was touched, no commit, no signal sent to any process, Phase 1 D and the
matched-benchmark run untouched.

**Q:** Recover Phase 1 experiment D without disturbing the live run: add a checkpointed,
timeout-bounded runner that is scientifically equivalent to `run_experiment_d.py`, prove
the equivalence with tests, smoke it, commit source/tests/docs only, then launch the full
60-background recovery as a persistent systemd unit with 16 concurrent subprocesses, one
BLAS thread each, and a 3600 s per-background wall limit.

**A:** The diagnosis held up under measurement. `Pool.map` returns only when every task
returns and writes nothing before that, so the live run had spent 5 h 40 min wall and
9 h 26 min CPU with a completely empty output directory. Timing every background
individually showed the median finishes in ~25 s and only four draws — 14, 19, 37, 44 —
exceed 45 s. Four pathological backgrounds out of sixty were blocking the other fifty-six.

The equivalence was made structural rather than asserted: the new runner imports
`run_experiment_d` and calls that module's own `_backgroundTask`, `_pairSummary` and
`PAIRS`, so there is no second copy of the model to drift. A test greps the runner's
source and fails if any of those is ever redefined there. Checkpoint rows are stored as
JSON, not TSV, because `writeTable`'s `%.12g` truncates float64 — formatting happens once
at merge, so the merged `interactions.tsv` is byte-identical to a direct write.

One real defect was found and fixed while writing the tests: the controller did not
forward `--n-backgrounds` to its children, so a reduced-size run would have handed each
child a different LHS matrix. The sample-hash gate would have caught it as a rejected
checkpoint rather than producing wrong numbers, but the flag was silently broken.

38 new tests, 187 total (149 pre-existing, none modified, none regressed). The strongest
are byte-for-byte equality of the merged table against a direct `writeTable`, and
float-for-float equality of the checkpointed path against a direct `_backgroundTask` call
on four measured-fast backgrounds. Timeout is tested through the real kill path, not a
mock. Smoke on backgrounds 5/48/56 reproduced identical artifact hashes across two
independent runs and resumed in 0.1 s versus 18.1 s with every checkpoint file untouched.

Committed as `ee64a3f` (source + tests) and a follow-up docs commit for
`PHASE1_D_CHECKPOINTED.md`. Launched as
`proteostasis-phase1-D-checkpointed-20260731T223225-0500.service`, MainPID 1622784, into
`results/phase1/D_checkpointed_20260731T223225-0500/`. At +40 s: active, exactly one unit,
1 controller + 16 workers all inside the unit's cgroup, 19 backgrounds already
checkpointed, stderr empty everywhere. The original unit
`proteostasis-phase1-20260731T162946-0500-D.service` (MainPID 1269455) was never
signalled and is still running with its output directory still empty.

**Note on timeout semantics:** a background that exceeds the wall limit is recorded as
`unresolved_timeout` — it contributes no rows, enters no interaction summary, and is not
counted in `n_errors`. It is a statement about the budget, not about the model, and must
not be reported as a numerical failure.

**Note on the duplicate-launch check:** the first pre-launch scan matched its own shell,
because the `-c` string of the scanning bash contains the pattern being searched for.
This is the exact pitfall documented in `results/phase1/ops/scan_live_experiments.sh`.
The check was redone as an argv-element test rather than a substring match.

## 2026-08-01

**Q:** Finish only the tracked documentation, tests and commit for the already-computed
experiment D closure. Do not recompute the 10,000-bootstrap analysis unless a validation
file fails, and do not alter raw Phase 1 results. Write `EXPERIMENT_D_FINAL.md` and
`PHASE2_CLOSURE_FINAL.md`, add a compact tracked validator, update STATUS.md, run the
suites, and commit as `Close Phase 2 interaction analysis`.

**A:** The bootstrap was not recomputed — `validation.json` reports 36/36 with
`all_passed: true`, so the existing output stood. Committed on `master` as
`Close Phase 2 interaction analysis`; full suite **249 passed**,
`scripts/phase2/check_d_closure.py` **177/177**.

A previous timed-out pass had already drafted both documents. Rather than trust them,
every descriptive claim was re-derived from the float64 checkpoint payloads, and that
found three real defects.

The largest was in the sensitivity section of *both* documents: "the chaperone-only
additive claim fails to reach a majority of backgrounds under every column" is false.
`conditional` is 0.6087 and the `usable` bounds are [0.5833, 0.6250] — both above one
half. Only the `requested` column fails to reach a majority (≤ 0.5000). The claim
actually fails because the grouped bootstrap CI for the majority fraction, [0.457,
0.739], straddles one half, and the three sensitivity columns are arithmetic
reweightings of the same 28 backgrounds rather than a second test. Conflating the two
would have implied the sensitivity analysis was doing inferential work it cannot do.

Second, `9.5e-10` as the "worst absolute discrepancy" of the `%.12g` TSV round-trip does
not reproduce from any column. The actual figures are 4.9e-12 *relative* (which is what
the validation check asserts) and, because burdens run to the censor at 1000, 4.8e-09
absolute in `excess_multiplicative` and 5.0e-10 in `excess_additive`.

Third, the draft overstated the influx-vs-nascent contrast. `nascent × total_capacity`
*is* supported on all three nulls, so "the synergy requires damage influx" and "the
interaction is driven by damage influx, not folding-machinery occupancy" were too strong.
The contrast that survives is the categorical one — synthetic collapse 13/46 against
43/46 — plus an order-of-magnitude smaller graded excess on the escape-free subset
(0.0028 against 0.0348). Scoped accordingly in the refined law, in §5.5, and in P4.

Everything else in the drafts verified exactly: the 470 lethal-single cells with Bliss
excess never negative, `mult_pred ≥ add_pred` in 100 % of both-singles-damaging cells,
`censored_12` identical to the blowup set, 40 of the 80 anti-additive cells carrying a
censored single, the 31/15 protective split with 25 and 3 additive majorities, and
Spearman −0.60. A numeric diff over the carried-forward §1–§3 confirmed zero drift from
the pending note, so every negative result survived intact; two caveats the draft had
trimmed (the weak markers' CIs touching 0.5, and the c05 design note) were restored.

The new tracked validator `scripts/phase2/check_d_closure.py` exists because the
evidence is gitignored while the prose is not, so nothing would otherwise notice the two
drifting apart. It re-checks the pinned source and artifact hashes, the counts, the nine
majority counts and verdicts, the collapse rates, the arithmetic of all 13 sensitivity
bounds, and that both documents still state those numbers. It carries no copy of the
estimator. Three negative controls were run to prove it is not vacuous, and the
absent-results branch returns `None` and exits 0 with an explicit `SKIP` rather than
reporting a success it did not verify — that branch is itself tested.

One test-gate defect was fixed: `TestDFinalOutputsAgreeWithRecomputation` skipped on
`D_final/d_final_results.json` but its `setUpClass` globs checkpoints out of the Phase 1
run root, so it would have errored rather than skipped had `results/phase1/` been cleaned
independently.

**Note:** no claim ledger exists in this repository (only `theory/SCOPE_AND_NONCLAIMS.md`),
so none was created. `OPS_SUBMISSION.md` was left uncommitted — it is an untracked
artefact of an earlier session and outside this task's scope.

## 2026-08-03

**Q:** Repair D013's caveat, run antecedent check A1 (does mistranslation damage
removal capacity, and does the identity survive), then preregister and run
spatial sequestration.

**A:** D013's bistability numbers were all computed under constant dilution — the
regime that predicts zero growth loss by construction. Caveat moved to the head of
the entry and pinned by test in both DECISIONS.md and STATUS.md.

A1 (D027): the identity survives capacity self-damage at machine precision, in
both a parametric mode (capacity falls with j) and a state-dependent one
(capacity falls with burden), down to 1.8% of nominal capacity. There is no
corrected form — the row operation needs only state-independent influx and exact
mass balance, so how parameters depend on j is irrelevant. What dies is the
SHORTCUT: {G=0} moves with the load, so fold-finding goes from two equations to
three. j_crit falls to 0.32x/0.13x but the slowing exponent does not move
(paired Wilcoxon p = 0.312) — collapse happens sooner, not steeper. Also recorded
that the numerical check carries no truncation term and so cannot fail unless mass
balance does (slope in h is -0.97, not +2).

**Q:** Re-anchor the band to Lindner's data before running sequestration; drop
the two-commit ordering rule.

**A:** The band was wrong by a factor of thirty. ">30% OF the loss of reproductive
ability" is a SHARE of the aging effect, not the effect. The full text gives the
effect as -3.95 ± 0.5% and the aggregate share as ~30-40%, so the measured
quantity is 1.2-1.8%. The misreading had survived two analysis cycles and been
written into POSTDICTION_PROTOCOL as a worked example. Band is now derived by
`aggregateAttributableLoss()`: [0.0104, 0.0178].

**Q:** Run sequestration.

**A:** FAILED, and it is the sharpest failure so far (D029). D024 verified on the
extended system first. 0 of 384 qualified cells in band, in both growth-cost
regimes; every bistable cell predicts 0.482-1.000, i.e. 27x-92x too severe.
Sequestration is not inert — with a_s exempt from the growth cost it gives 12
bistable cells against the control's 3 — but the high state it creates is always
complete arrest. Inverting the growth law: the model's second attractor sits
7.5x-254x more aggregate-laden than the cell that was measured. Next candidate,
from the observation: a size-limited deposit. Collision count now three, all
informative failures; the fold theorem is untouched by any of them.

Also found and recorded a process defect: under the ORIGINAL band the literal
criterion would have PASSED, carried entirely by k_seq=0 control cells, because
D028 never required the mechanism under test to be switched on. Now protocol
rule 6.

## 2026-08-04 (continued)

**Q:** Zenodo, then figures for bmb_v4.md.

**A:** Deposit minted — concept DOI 10.5281/zenodo.21794565, v1.0.0
10.5281/zenodo.21794566, uploaded file MD5-verified. Byline corrected to Kiran
Boggavarapu from ORCID 0000-0003-0751-6459 (given Kiran, family Boggavarapu); the
repo metadata and deposit were already right and the manuscript was the outlier.

**Q:** Figure 1, and is panel (b)'s horizontal tangent computable?

**A:** Initially no — `G_u > 0` structurally (all four terms of `G` push the same
way), so the aggregate coordinate has no horizontal tangent, confirmed across 40
draws. But that stopped one derivative short: Cramer on the other column gives
`du*/dj = -G_a/det J`, so the SOLUBLE coordinate turns where `G_a = 0`. That point
exists, is unique, and is generic (`det J = R_a.G_u = 2.027e-03`), and it is the
same locus where root-finding loses the nullcline. Panel (b) plots both
coordinates, all computed, with a x182 zoom inset. NOT a prediction:
`j_turn/j_crit = 0.9990`, and on the stable branch in only 6 of 30 draws.
Section 3.1 rewritten as three cases in Wang et al.'s vocabulary (Haldane type,
excluded by sign structure; output = u outside their framework; det J = 0 the
fold).

**Q:** Check whether section 5's headline used the degrading metric.

**A:** It did, and worse. TWO distinct defects (D035, D036): the metric divides by
`max(|det J|, |cross|)`, both of which vanish at a saddle-node (corr with |eig| =
-0.262, it degrades as the bracket tightens); AND three headline numbers came from
an undeclared 20-of-325 subsample, all three optimistic. Forensics: the draw was
uniform and seeded, so nothing selected them — two deviations are coin flips
(P = 0.501, 0.438) and the third, a maximum, understated with probability 0.939 =
1 - 20/325, deterministically. Rule: never subsample an extremum.

**Q:** Sweep every extremum; the rule is under-applied.

**A:** D037/D038. Found a fifth affected number (solver max, 6.652e-07 -> 7.56e-07)
and then a sixth (`s_u` width 0.890 is the regulation experiment's, not the kinetic
box's 0.876). Section 5 now reports median/p99/max per row with populations named,
since a maximum grows with population size and is not comparable across
denominators. FOLD_THEOREM.md decided LIVE and corrected in place rather than
bannered.

**Q:** Figures.

**A:** Fig 1 (theorem, two panels), Fig 2 (saturation), Fig 4 (the beta-indexed
prediction) built; data/figures/ builder extracts reduced arrays with provenance
per file. D039: Figure 2's specified `s_a` screen is NOT defensible — the
distribution has no gap and the median slides continuously with any floor
(0.090 to 0.355), and the first build's screen created a 5x figure-text
divergence. No screen, no exclusion; the complete 2884 reproduce section 6 exactly.

**Remaining:** Figures 3, S1, 5; then Zenodo v1.1.0 so the archive captures the
corrected state in one version.

**Q:** The screen finding is not really about a figure. Check whether any other
standing limitation in FOLD_THEOREM.md or STATUS.md is phrased as licensing an
action; record that the sixth number was a different failure type; state the 13.6x
position rather than leaving it to be reconstructed under review. Then figures 3,
S1, 5.

**A:** D040 and D041.

Swept every limitation in `theory/`, `notes/` and STATUS.md, classified as
observation / restriction / direction. **One** in the licensing form — the `s_a`
screen already found. Near misses left alone and listed: "should not be quoted
without X" is a restriction (adds a requirement, opposite direction); the `phi` ->
`phi_enz` substitution is an algebraic identity; PREDICTIONS.md's "should"s are
hypotheses about the world. Recorded as an audit form and pinned by
`tests/phase3/test_limitation_forms.py`, which permits the banned phrasing where
it is QUOTED and forbids it where asserted — an absence test would have failed on
exactly the text that fixes the defect.

The sweep also found that **D038's relabelling had reached the manuscript only**.
STATUS.md and FOLD_THEOREM.md still read "up to 13.6x", "as much as 13.6x" and
"max 0.623" against a corrected manuscript — the lineage-split failure, inside the
correction that introduced the warning against it. Corrected at all sites, and the
position is now written down: the design was not sized to bound the per-network
spread and we cannot bound it.

**Q:** (during Figure S1) reconcile its numbers against the text first.

**A:** D041. Section 5's normalisation contrast **does not reproduce**. -0.262 and
+0.060 appear in none of 48 correlation definitions over the complete 325, nor on
the 20-state draw, recomputed from the run root through `determinantIdentity`
itself. Measured: **-0.835 / +0.041** log-log, -0.221 / +0.073 raw. The medians and
p90s in the same table reproduce exactly, which is what hid it — a wrong population
moves every entry, an ad-hoc computation moves only some. The claim STRENGTHENS:
the degrading metric degrades more decisively than reported, and every correction
moves against the earlier reporting. `figS1.py:captionNumbers` now owns every
number in that paragraph and the test asserts the manuscript against the function.

A **seventh** affected number: "the single state bracketed to eigenvalue -4.2e-9"
is a MINIMUM written as a definite description, missed by two exhaustive extremum
sweeps on grammar alone. Full population reaches 8.10e-10 with three states
tighter, so the uniqueness claim was false as well. Also: the corrections paragraph
said "Three of these values" and listed four.

**Figures 3 and S1 built.** Fig 3 — front spans 0.227-0.965, optimum solved on the
boundary at `j/j_crit = 1.000000` against the grid's 0.8975; connector drawn grey
rather than colour-mapped because nothing was solved between the 13 points. Fig S1
— the correlation panel, and it supplies a fact the text lacked: the log-log slope
is **1.00**, so parallelism error is *proportional* to bracket looseness, a sharper
statement of bracket tolerance than the correlation alone.

**Remaining:** the optional sixth figure (spec not recoverable from the repo
record), then Zenodo v1.1.0.

**Q:** Build Figure 5 (phi_enz and delta against mu, supporting Corollary 3), then
wire all the figures into the prose — callouts that read the panel, no duplication,
four specific missing statements, caption/text cross-check, placement and size.
Then cut v1.1.0.

**A:** D042.

**Numbering.** The dilution figure supports Corollary 3 in section 4.2, which
precedes sections 6 and 7, so by order of first mention it is **Figure 2** and the
saturation, front and beta figures shift to 3, 4, 5. It was specified as "Figure 5
(optional sixth)" on BUILD order, which is not mention order — the same confusion
that had the beta figure built before the two preceding it. Filenames are now
semantic (`fig_dilution.py`, `saturation.tsv`) with the number in one `FIGURE`
constant per script, so a reorder — including dropping this figure on page count —
touches one line instead of renaming files. Manifest hashes unchanged by the
rename, confirming only names moved.

**Four missing statements added.** Section 3.1 now gives Fig. 1b's two ordered
singularities with values and, in one sentence, that the root-finder's failure and
the horizontal tangent are one locus. Section 6 states the no-screen decision and
the sensitivity ladder in the prose. Section 5's contrast was already corrected.

**Duplication resolved per case:** table keeps the medians, caption keeps the widths
(section 9 now points at the panel instead of restating 0.876); the regulation
widths stay in prose because no figure carries that experiment.

**Two divergences found by wiring, neither by the test suite.**
1. `fig_beta` used a single damping of 0.35 against a table built from the per-beta
   measurement — every figure row ~1.5% off section 8.4. Damping is beta-dependent
   (0.3462 at beta=1, ~0.355 below), now stored at full precision and interpolated.
   Storing it at 3 decimals was tried first and lost agreement with D033 in the
   fourth digit, which is how the precision requirement was found rather than
   assumed. A slow test recomputes the cache so it is not an orphan.
2. A Fig. 5 callout claimed the requirement never comes within fifteenfold of the
   measured load. True at the three marked beta, **false at the left edge of the
   figure's own x-axis**, where it closes to 3.19x. Corrected to "3x to 214x".
   A smaller one at the same site: 15.94x had been promoted to a claim of "16x or
   more" by borrowing the table's rounding.

`tests/figures/test_figure_wiring.py` pins numbering-by-first-mention, prose
citation before embed, caption-vs-text agreement against the generator, and the
84/174 mm width and 234 mm height limits.

**Q:** Prepare the LaTeX version and PDF; create the supplementary and finish it.

**A:** D043. `scripts/manuscript/to_latex.py` generates `bmb_v4.{tex,pdf}` and
`bmb_v4_supplementary.{tex,pdf}` from the markdown. The markdown stays the single
source; the `.tex` files carry a do-not-edit banner. Supplementary is its own
document with `\thefigure` = S1, S2, ... so the sequences cannot collide.

276 backtick spans, exactly TWO of them literal code, so the rule is inverted:
math by default, and a code-shaped span not in `CODE_SPANS` fails the build rather
than rendering as math.

**Four silent corruptions, every one found by rendering pages and reading them:**
1. `j_crit` unbraced in math is `j_c` followed by "rit" — a different symbol, no
   error, dozens of occurrences.
2. An unescaped `%` in Figure 5's caption (the only one quoting percentages)
   opened a LaTeX comment, swallowed the closing brace of `\caption{}`, made
   pandoc escape the whole raw block, and the orphaned `\centering` **centred
   every paragraph from section 10 to the end**. The build reported "figures 6".
3. `width=\linewidth` scaled the 84 mm panels to the 160 mm text block, enlarging
   their 7 pt labels by 1.9x and discarding the whole width discipline.
4. pandoc refuses `$...$` with a space before the closing delimiter and escapes
   both dollars; the Greek substitution left a trailing space.

Plus double-numbered headings ("0.1 1. Introduction") and a byline printed twice.
The heading fix asserts that LaTeX's counter reproduces the markdown's own numbers
before stripping them, because the prose cross-references them.

None of these raised an error. **Four of the last seven real defects in this
project were found by looking at output, none by a passing test.**

Main text 21 pages, supplementary 1. `tests/manuscript/test_to_latex.py`, 15
checks.
