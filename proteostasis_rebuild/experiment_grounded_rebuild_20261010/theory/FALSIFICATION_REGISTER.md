# Falsification register

Pre-registered before any manuscript. Nothing here has been run. Each entry states what the theory
predicts, the **minimum** readouts that discriminate the models, what result kills the prediction,
and what a negative result would and would not mean.

A prediction is non-tautological only if a competing model makes a different call. Each entry names
the competitor. Entries P1–P3 are prospective. Section R lists the retrospective checks that must
**not** be presented as predictions, because the theory was built knowing them.

---

## P1. The HtpG disposal sign flips with genetic background

**Prediction.** From the sign rule (EQUATIONS A8), deleting *htpG* improves growth when the
disposal odds exceed the client's productive alternatives, and stops helping — or hurts — when
the chain's alternative fate becomes loss anyway. Two specified moves change the inequality in
opposite directions:

| background | what it changes | predicted ΔhtpG effect vs ΔdnaKJ alone |
|---|---|---|
| ΔdnaKJ (reference) | — | improves growth (Fauvet; retrospective, see R1) |
| ΔdnaKJ ΔclpB | lowers recovery of aggregates, raises k_a^eff | improvement **shrinks or reverses** |
| ΔdnaKJ ΔhslUV | removes the disposal protease, lowers effective q | improvement **shrinks** |
| ΔdnaKJ + GroEL/ES overproduced | raises productive alternative p_G k_G G | improvement **grows** |

**Minimum readouts.** Growth rate (or plating efficiency) for these four backgrounds ± *htpG*, at
one permissive and one restrictive temperature. Eight to ten strains. No proteomics needed for the
sign test.

**Competing calls.** A model where chaperones are beneficial by construction predicts ΔhtpG hurts
in every background. An optimizing-feasibility model predicts ΔhtpG never helps (route removal
cannot raise the optimum) unless a proteome-cost term carries it, in which case the cost term
predicts a background-independent effect proportional to HtpG's proteome share, not the ordering
above. The ordering is the discriminator.

**Killed if.** The ΔclpB and GroEL-overproduction arms move the same way, or ΔhtpG helps equally in
all four backgrounds.

**A negative result means** the sign rule's ingredients are wrong or dominated by something else —
most likely σ32 induction, which changes every pool at once and is not in the model. It would not
by itself refute the necessary conditions in §2.

---

## P2. Genetic-interaction sign follows route structure, not abundance

**Prediction.** Two machines on the same serial path share a cut, so losing both is not worse than
multiplicative. Two machines offering parallel qualified routes to the same clients aggravate.
Concretely: *tig* × *dnaK* aggravating; *dnaK* × a *groL* hypomorph not aggravating; *hslVU* ×
*clpPX-lon* aggravating.

**Minimum readouts.** S-scores (or plating efficiencies) for these three pairs plus their singles,
at one temperature. This is a sign test with no free parameters. Babu 2014 Table S2 would supply it
if retrieved — **currently blocked**, and it lists significant hits only, so absent pairs cannot be
scored neutral. Alternatively ~12 plates made directly.

**Competing calls.** An abundance-ranked single-pool model predicts the most aggravating pairs
involve the scarcest machines (ClpB) regardless of client overlap. That is a different ranking.

**Killed if.** The signs come out opposite, or ranking by abundance predicts the data as well as
ranking by route overlap.

**A negative result means** the serial/parallel abstraction is too coarse for growth data — likely,
since growth integrates everything. It would not refute the quota or occupancy conditions.

---

## P3. Hysteresis under a graded load sweep, only if sequestration is real

**Prediction.** Ramp an aggregation-prone load up and then back down through the same levels. If
aggregates hold chaperone (β > 0), the down-ramp retraces a different curve and a high-aggregate
state persists at loads that were previously fine — a saddle-node with history dependence
(L_crit = (μ+δ)a²/(4bk_a)). If β = 0, up and down coincide; overload shows smooth decline or a
single kink, no hysteresis.

**Minimum readouts.** One inducible aggregation-prone substrate with ~6 graded induction levels,
up then down, measuring growth rate and insoluble fraction at each step. Plus **one measurement
that is currently missing entirely**: how much chaperone is held in the insoluble fraction per
aggregated chain (β) — a co-sedimentation stoichiometry.

**Competing calls.** FoldEco contains the sequestration mechanism (KJE binds aggregates 1:1) but
reports no bistability; Santra and FoldME report none. A no-feedback model predicts no hysteresis
at any load.

**Killed if.** No hysteresis appears at any load while β is measured to be substantially non-zero;
or hysteresis appears where β ≈ 0 (which would mean a different feedback and the named mechanism
is wrong).

**A negative result means** the specific sequestration feedback is not operative, which is a
genuine result and must be reported as such. It is **not** evidence that proteostasis has no
thresholds — and the theory explicitly declines to claim universal collapse thresholds either way.

---

## R. Retrospective checks — not predictions

| id | fact | why it is not a prediction |
|---|---|---|
| R1 | ΔhtpG improves ΔdnaKJ growth (Fauvet 2021) | the disposal route was put in the model to accommodate it |
| R2 | Δtig ΔdnaK lethal 30 °C, viable 23 °C, rescued by GroEL/ES (Vorderwülbecke 2004; Genevaux 2004) | used to reject the hard serial requirement; it is a structural constraint |
| R3 | protease double deletion aggravating at 30 °C, non-monotone in T (Kanemori 1999) | used to reject monotone temperature-as-load |
| R4 | synonymous mDHFR variants reduce GroEL dependence (Reingewertz 2024) | motivates the codon-timing variable; also the translation-rate link there is computational, not measured |

---

## Blockers before any manuscript

1. **No in-vivo arrival flux or service time** for any machine–client pair. Until then the
   occupancy condition has an unknown left side and no margin may be quoted.
2. **β unmeasured**, so no collapse threshold.
3. **Babu S2 unretrieved**, so P2 is blocked on existing data.
4. **Calloni and Kerner supplements unretrieved**, so client-resolved demand remains unavailable.
5. **No phase-matched parameter set exists.** Abundances are exponential-phase; kinetic rates are
   in vitro; the error tables are standing frequencies. **No E. coli parameter set may be called
   calibrated.**
6. **σ32 regulation is not in the model**, so every overload statement is open-loop.
7. **Multi-client sign transfer under shared saturated pools is unproven** (EQUATIONS A8, Result 4).

Drafting rule: the manuscript waits until at least one of P1–P3 has been run, and the novelty
claims in PRIOR_ART §"what this work could newly test" have had a targeted citation search.
