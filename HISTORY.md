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
