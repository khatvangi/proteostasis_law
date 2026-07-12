# P1 combined outline — "Capacity + Codon" proteostasis paper

**Draft outline for review. Nothing written yet — this is the map so you can
react before prose.** Conventions borrowed from `protein-brain/SPEC.md`: real
file paths cited, unknowns marked `VERIFY`, no fabricated numbers (every figure
here is traced to a source file on disk).

---

## The merge in one sentence

Take the **burden-capacity envelope theory** from Manuscript A
(`proteostasis-paper/manuscript/MANUSCRIPT.md`) and the **context-conditional
codon-deployment evidence** from Manuscript B (`manscuript-final.md`), and fuse
them into a single paper whose claim is: *extant translation operates inside a
finite proteostasis envelope, and synonymous codon choice is one of the levers
cells use to stay inside it.* The doublet-exclusion / "synonymy necessity"
argument is **removed** and cited out to P3 (`triplet-proof`, submitting now)
and P6 (`triplet-intersection`).

## Working title (pick one; all avoid "principle" and "triplet")
1. *Translation operates within a finite proteostasis envelope, and codon choice
   is a load-allocation lever within it*
2. *A burden-capacity envelope for translation: synonymous codons as
   state-dependent load allocation*
3. *Proteostasis as an operational constraint on translation: theory and
   codon-deployment evidence*

Target journal: **MBE or JMB** — `VERIFY` (dashboard lists both; P1 session-boot
says "ask Kiran"). Authorship: `VERIFY` (sole-author vs co-authors — dashboard
flags unknown).

---

## What comes from where

| Section | Source | Fate in combined paper |
|---|---|---|
| Burden-capacity inequality `B_total ≤ C_buffer`, 4 burden terms | A (Theory) | **Core theory, kept** |
| Reduced ODE, saddle-node → buffered/vulnerable/overload regimes | A (Theory) | **Kept**, reframed as threshold-existence only |
| Error burden measurable (Landerer) | A + B | **Kept, merged** |
| Buffering/clearance finite (Farkas, Chuang, de Marco) | A | **Kept** |
| Synonymous effects state-dependent (Yang, Walsh) | A | **Kept** |
| μ/ν (tAI) codon coordinates in E. coli | B (Methods) | **Kept** — the operational axes |
| Metal-site codon enrichment (Asp/His high-μ, Cys/Glu low-μ) | B (Result II) | **Kept — headline empirical result** |
| Cross-species Δ_A conservation (ρ=0.54–0.77) | B (Result III) | **Kept** |
| Matched-null (z=2.99, p=0.0014) | B (Result IV) | **Kept** as support, not primary |
| Synonymy necessity constraint / doublet exclusion | B (Theory, Result I) | **REMOVED → cite P3/P6** |
| Companion "decoder-capacity A_eff=9" references | B | **Removed** (that's P5/P6 territory) |
| supraadditivity prediction | A (Predictions) | **Kept**, labeled untested |

---

## Section-by-section outline

### Abstract (~220 words)
- Translation imposes burden as well as product; four coupled burden terms.
- State the envelope: viability requires `B_total ≤ C_buffer`; near saturation
  the burden becomes threshold-like (reduced ODE, saddle-node).
- Codon choice enters as a **load-allocation lever**: synonyms differ in
  mistranslation hazard μ and translational supply ν, and cells deploy them
  non-uniformly and **non-monotonically** at high-criticality sites.
- Empirical spine: metal-site enrichment shows opposite μ-preferences across
  amino acids (rules out single-axis optimization); operational spread is
  conserved across ~2 Gyr.
- Explicitly **not** claimed: triplet-code origin; universal cross-organism
  threshold. (Keeps the C9/C10 guardrail from `CLAIM_MATRIX.md`.)

### Significance statement (A already has one — adapt)

### 1. Introduction
- Para 1: translation = output + burden (from A intro).
- Para 2: the four burden terms and the coupling claim (cascade), stated as a
  framework/hypothesis, **not** a "principle" (critical-eval fix C-1).
- Para 3: synonymous codons are not interchangeable — they carry distinct μ/ν;
  Drummond–Wilke mistranslation-cost framing (from B intro).
- Para 4: the question — how does codon deployment sit inside the envelope? Two
  parts: (i) the envelope exists and is finite; (ii) codon choice is used to
  allocate load within it.
- Para 5 (**neutral-forces caveat, critical-eval fix I-5**): much codon usage is
  mutation/drift; the load-allocation reading applies to the *selected*
  component. State this up front.

### 2. Theory
- 2.1 Burden-capacity envelope (`B_total ≤ C_buffer`); define B_error, B_fold,
  B_agg, B_qc, C_buffer (table from repo CLAUDE.md).
- 2.2 **Operationalization paragraph (critical-eval fix C-2 + rec #2):** how each
  term maps to an observable, how double-counting is avoided (a misfolded
  protein is simultaneously B_fold/B_agg/B_qc). This is the fix that makes the
  inequality testable rather than verbal.
- 2.3 Collapse threshold: reduced ODE `dx/dτ = λ − x − νx/(1+x) + χx²`,
  saddle-node bifurcation → J_crit. **Frame as threshold-existence only**, not a
  fitted organism model (critical-eval fix I-4; honest per A's own scope note).
- 2.4 Bridge to codon architecture: `J_in = Σ_p E_p Σ_i κ_{p,i} q_{p,i}`,
  `q(c) = μ(c)[1−S(c)] p_let` — the accounting identity from B that links codon
  choice to the envelope.
- 2.5 Load-allocation cost `L_i(c|A) = λ_μ μ(c) + λ_ν (ν(c)−ν*(x_i))²` — the
  site-specific tradeoff that predicts non-monotone deployment.

### 3. Results
- **R1. Burden terms are measurable and capacity is finite.** Merge A's
  error-burden (Landerer msae048) + finite-capacity (Farkas 2018) evidence.
- **R2. Codon operational axes and non-monotone deployment at metal sites.**
  B's Result II — N=17,166 metal-binding residues; Asp GAC / His CAC enriched at
  *higher* μ, Cys UGU / Glu GAA at *lower* μ. Opposite directions = no
  single-axis rule. (`VERIFY` all μ/OR numbers against
  `codon-deployment/` sources — see blocker B4, P2 audit found μ errors.)
- **R3. Operational diversity is conserved across species.** B's Result III —
  slow-Leu usage 15–39%; Spearman ρ = 0.54–0.77 across E. coli/B. subtilis/
  S. cerevisiae.
- **R4. Deployment structure exceeds matched-null.** B's Result IV — z=2.99,
  p=0.0014, framed as supporting not primary.
- **R5 (optional). Threshold behaviour is consistent with the two-pool bound.**
  Pull the P1 two-pool ODE result *only if the ODE/arithmetic reconciliation
  (blocker B1) is resolved.* See "Open blockers."

### 4. Discussion
- What the results show: codon choice is not single-axis; envelope is finite.
- **Alternative explanations (critical-eval fix I-3 + rec #3):** tRNA pool
  remodeling, mRNA stability, metabolic coupling as non-proteostatic accounts of
  state-dependent synonymous effects. Argue parsimony, don't assert proof.
- Relation to existing theories (error-minimization, robustness) — and a
  one-paragraph pointer to the companion triplet-architecture paper (P3) for the
  code-length question, which this paper does **not** address.
- 4.x **Formal Limitations subsection (critical-eval fix rec #5):** E. coli-
  weighted evidence; ODE not fit to data; supraadditivity untested; μ measured
  in one organism.

### 5. Predictions
- Stress-dependent selection; vocabulary starvation; shielding-error alignment
  (from B); supraadditivity (from A, labeled untested — critical-eval fix C-1).

### 6. Methods
- Merge B's Methods (μ from Landerer/Mordret MS data; ν=tAI via dos Reis;
  MetalPDB→UniProt mapping; matched-null construction; Δ_A metric) with A's
  ODE/nondimensionalization. Add reproducibility statement (critical-eval fix
  M-5): evidence ledger + consistency tests available.

### Figures (combined set — **most need to be built**, see blocker B2)
1. **Envelope + threshold** (adapt A's `stability_regions.svg`).
2. **Codon operational space (μ,ν)** — new; B Fig 1B.
3. **Metal-site deployment** — new; B Fig 3 (enrichment + OR + μ,ν positions).
4. **Cross-species conservation** — new; B Fig 2B + Δ_A ranks.
5. **Matched-null** — new; B Fig 4.

---

## Open blockers to resolve before/while drafting

- **B1 — ODE vs arithmetic bound direction (SCIENCE).** `proteostasis-P1/README.md`
  says the ODE bound sits ~30–100× *below* the arithmetic bound; the latest
  `two_pool_summary.md` (post-empirical-N) shows it ~10–100× *above* (looser:
  arithmetic `[1e-4,1e-3]` vs ODE `1.00e-2`). These contradict. Decide the real
  relationship and rewrite the "two bounds converge" narrative before any
  two-pool result enters this paper (R5). Until then, leave R5 out.
- **B2 — Figures.** 4 of 5 combined figures do not exist as files; only A's 2
  SVGs are on disk. Must be generated from `codon-deployment/` + `proteostasis-P1/`
  data.
- **B3 — Triplet content excision.** Remove B's synonymy-necessity constraint and
  companion-paper A_eff=9 references; replace with a citation to P3. Verify no
  orphaned cross-references remain.
- **B4 — μ-value audit.** The P2 data audit (`codon-deployment` DATA_AUDIT) found
  μ values wrong by 2–4 orders of magnitude and Cys/Glu orderings reversed. R2's
  numbers must be re-pulled from Landerer 2024 Data_S2 before they can be trusted.
- **B5 — "principle" → "framework/constraint"** throughout (critical-eval C-1).
- **B6 — Authorship + target journal** (`VERIFY` with Kiran).

## Sibling papers (context, not merged)
- **P3 `triplet-proof`** — triplet architecture / error-minimization; submitting
  now; owns the code-length claim this paper cites out to.
- **Foldon paper `foldon_project`** — 4-layer contact architecture (PNAS/PLOS);
  separate. Note its codon layer-selection result is **negative** (p=0.43), which
  is *consistent with* — and worth a one-line acknowledgement in — this paper's
  claim that codon load-allocation is site-criticality-driven, not
  folding-layer-driven.
