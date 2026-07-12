# Proteostasis paper family — collected state (scan of /storage/kiran-stuff)

Scanned 2026-07-11. This is the "collect" pass: current on-disk state of the
proteostasis manuscript and everything complementary to it. Facts only; open
decisions are flagged, not resolved.

## The family (six-paper program + workspaces)

| ID | Short title | Workspace | Stage (on disk) | Canonical manuscript file |
|----|-------------|-----------|-----------------|---------------------------|
| **P1** | Proteostasis Law | `proteostasis-P1/` (compute) + `proteostasis_law/` (drafts) | compute done; drafting not started in-workspace | **two candidates — see below** |
| P2 | Error-matched / supply-diverse codon deployment | `codon-deployment/` | draft, open post-audit issues | `codon-deployment/JME/MANUSCRIPT.md` (2026-04-16) |
| **P3** | Triplet architecture / MC error-minimization | `triplet-proof/` | **near-submit (JME)** — the "complementary" paper | `triplet-proof/manuscript_JME.md` (2026-03-13) |
| P4 | aaRS promiscuity | `aaRS/` | near-submit | `aaRS/paper-jan-26/MANUSCRIPT_COMBINED.md` |
| P5 | Decoder capacity / β-invariance | — | PARKED | none |
| P6 | Intersection of constraints (why triplets) | `triplet-intersection/` | compute; skeleton approved | `triplet-intersection/manuscript_outline.md` (2026-04-16) |

Additional related dirs found in scan (not in the six-paper grid):
`why-codons-are-triplets/` (earlier P6 iteration, 2026-03-13),
`codon_metal_analysis/`, `codon_project/`, `CodonFM/`, `foldon_project/`,
`exon_to_folding/`. (P2's workspace `codon-deployment/` is in the grid above,
not here.)

## P1 — the proteostasis paper itself

P1 is spread across **two workspaces and two divergent manuscripts**:

- **Compute** lives in `proteostasis-P1/` — two-pool ODE robustness extension.
  All sub-tasks report done (arithmetic, essential-gene bound, paired MC, v2
  mass-balance cap, empirical-N). Latest run = empirical-N bootstrap of N_prot.
- **Manuscript A** — `proteostasis_law/proteostasis-paper/manuscript/MANUSCRIPT.md`
  (2026-03-15, ~3,700 words). Title (from the file itself): *"The proteostasis
  principle: translation operates within a finite burden-capacity envelope."*
  Burden-capacity framework
  (B_total ≤ C_buffer). Passes all 15 consistency tests. Has 2 SVG figures.
  Carries a detailed `CRITICAL_EVALUATION.md` with a concrete fix list.
- **Manuscript B** — `proteostasis_law/manscuript-final.md` (2026-04-16, ~5,400
  words, newest, named "final", filename has a typo). Title: *"Context-Conditional
  Codon Deployment Under Proteostasis Constraint."* Codon-deployment content
  (μ/tAI axes, metal-site enrichment, cross-species Δ_A) framed under a
  proteostasis viability condition, PLUS a "synonymy necessity constraint"
  excluding doublet codes.

These are **different papers** (burden-capacity vs codon-deployment), not two
drafts of one. Deciding which is the P1 spine — and where Manuscript B's
content belongs relative to P2/P3/P6 — is the first open decision.

## Scientific-readiness issues found (blockers)

1. **P1 directional contradiction (stale README).**
   `proteostasis-P1/README.md` states the ODE bound "sits ~30–100× **below**"
   the arithmetic bound. But the latest `two_pool_summary.md` (post-empirical-N)
   shows the opposite: arithmetic window `[1e-4, 1e-3]`, ODE operational bound
   `1.00e-2`, with "inverse gap (ODE / arithmetic lower): ×100.03" and
   "×10.00" for the upper — i.e. the ODE bound is **10–100× above (looser)**,
   not below (tighter). The README predates the empirical-N rerun and is stale.
   Must be reconciled before P1 drafts, and the "two bounds converge" narrative
   re-examined (they now differ by 1–2 orders of magnitude in the loosening
   direction).

2. **manscuript-final.md figures do not exist.** Fig 1–4 are referenced in text
   and legends, but `proteostasis-paper/figures/` holds only the 2 SVGs for
   Manuscript A. All 4 codon-deployment figures would need to be generated.

3. **Triplet-origin scope tension.** The repo's own CLAUDE.md/claim-matrix
   *rejects* C10 (proteostasis explains triplet-code origin). Manuscript B
   derives a doublet-exclusion / synonymy-necessity argument that brushes this
   guardrail. The triplet-necessity claim is owned by **P3 (triplet-proof)** and
   **P6 (intersection)**. Overlap must be resolved so the papers don't
   self-collide or double-claim.

4. **MANUSCRIPT.md critical-eval fix list (Manuscript A).** From
   `CRITICAL_EVALUATION.md`: downgrade "principle" → "framework/constraint";
   add a measurement-operationalization paragraph (how to measure B_total /
   avoid double-counting); strengthen alternative explanations for
   state-dependent synonymous effects; add rough parameter bounding for the ODE;
   split a formal Limitations subsection; add the neutral-forces (mutation/drift)
   caveat for codon usage.

5. **Authorship / target journal unknown** for P1 (sole-author vs co-authors;
   MBE vs JMB) per the dashboard session-boot notes.

## triplet-proof (P3) — the complementary paper

`triplet-proof/` — *"Triplet architecture enables deep error-minimization in
the genetic code"* (JME submission). Has: `manuscript_JME.md`,
`supplementary_materials_JME.pdf` (2026-05-12, most recent activity),
`cover_letter_JME.md`, `submission_checklist_JME.md`, figures/, full MC code
(1M-sample runs). Headline results: SGC beats 1M random codes (p<1e-6);
architecture (doublet→triplet) shifts z by 11.0 vs 1.8 for alphabet 10→20;
wobble absorbs 69% of single-nt errors as synonymous. This is the furthest-along
of the triplet-family papers and the natural companion cite for P1's
architecture claims.

## The stale dashboard
`/storage/kiran-stuff/dashboard/` (2026-04-17) still lists P1's next blocker as
"finish paired MC + essential bound" — but those are now done. Confirms the
dashboard is a stale cross-reference, consistent with `protein-brain/SPEC.md`.
