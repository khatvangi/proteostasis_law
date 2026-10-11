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

## P3. A bounded hysteresis loop under a graded load sweep, only if sequestration is real and capture is fast enough

*Revised 2026-10-10 after the independent audit (CORRECTIONS.md). The earlier P3 predicted a
saddle-node "runaway". That came from letting free chaperone go negative, and it is withdrawn.*

**Prediction (reduced model, EQUATIONS A6, Theorem A6).** Ramp an aggregation-prone load up, then
back down through the same levels. Hysteresis appears only if **both** conditions hold:
(i) aggregates hold chaperone (β > 0), and (ii) at full pool, chaperone capture beats spontaneous
aggregation plus dilution (k_c M^tot > k_a + μ). Where it appears, it has a definite shape:
- on the way up, the insoluble fraction jumps at L_crit;
- on the way down, it stays high until the lower load L_sat;
- the upper branch is **finite and linear in load**, A_hi = k_a L/((k_a + μ)(μ + δ)), with no
  divergence while growth or aggregate clearance continues;
- on the upper branch the chaperone held in the insoluble fraction **plateaus at the total pool**;
- the loop's relative width is ((x − 1)/(x + 1))², with x = k_c M^tot/(k_a + μ).

If (i) fails, or (ii) fails, the up and down sweeps coincide, with a smooth response or a single
kink.

**Minimum readouts.** One inducible aggregation-prone substrate, up then down. Use a fine load
ladder near the transition, not ~6 levels: the loop can be a few percent of L_crit wide (2.8 % at
x = 1.4). Six coarse levels can step straight over it, and a null from them would be
uninformative. At each step measure the following.
1. **Growth rate μ**, because dilution enters both L_crit and A_hi.
2. **Insoluble fraction of the substrate**, for A.
3. **Chaperone in the insoluble fraction per aggregated chain** (co-sedimentation stoichiometry,
   for β), and **whether it saturates at the total pool**. This tests the cap directly. Unbounded
   growth of held chaperone past M^tot would contradict mass balance and would mean the assay is
   wrong.
4. **Aggregate clearance after induction is switched off** (a chase), for μ + δ.
5. **The capture-versus-aggregation ratio at low load**: the fraction of the substrate aggregating
   with the chaperone present versus depleted or deleted, which bounds k_c M^tot against
   k_a + μ. Condition (ii) is checked against this measurement before the sweep is read.

**Competing calls.** FoldEco contains the sequestration mechanism (KJE binds aggregates 1:1) but
reports no bistability. Santra and FoldME report none. That is consistent with condition (ii)
failing at their parameters, but this has not been checked and is not claimed. A no-feedback model
predicts no hysteresis at any load. A "runaway" model predicts unbounded growth of the insoluble
fraction above threshold while μ > 0. The corrected theory rules that out unless μ + δ → 0.

**Killed if.**
- β is measured clearly non-zero, condition (ii) is measured to hold, the ladder is fine enough to
  resolve the predicted width, and no loop appears.
- A loop appears where β ≈ 0, or where condition (ii) clearly fails. That would mean a different
  feedback is operating, and the named mechanism would be wrong.
- The upper branch keeps growing at fixed load while μ + δ is measurably positive. That would
  falsify the capped-pool model itself.

**A negative result means** the sequestration feedback, or the capture condition, is not operative
at that load. That is a real result and must be reported as such. It is **not** evidence that
proteostasis has no thresholds, and the theory declines to claim universal thresholds either way.

**Does the hysteresis prediction survive the correction?** Yes, mathematically, in narrower form.
The runaway and "two states or none" are gone. A loop between two finite states survives, under
one more unmeasured condition than before, and with an upper branch whose level is predicted in
closed form. It is now a sharper and more falsifiable claim, but a more fragile one: the window may
be narrow, and the fold condition (ii) holds only in the reduced model. With the client-bound term,
recycling of disaggregated chains, or several pools, it has not been derived. **No experiment has
been run.** The verdict is "survives as a conditional prediction", not "supported".

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
2. **β and k_c unmeasured**, so neither the existence of the A6 loop (k_c M^tot > k_a + μ) nor its
   loads (L_sat, L_crit) can be stated, and no collapse threshold is claimed. The fold condition is
   derived only for the reduced one-variable model: client-bound term omitted, δ lumped.
3. **Babu S2 unretrieved**, so P2 is blocked on existing data.
4. **Calloni and Kerner supplements unretrieved**, so client-resolved demand remains unavailable.
5. **No phase-matched parameter set exists.** Abundances are exponential-phase; kinetic rates are
   in vitro; the error tables are standing frequencies. **No E. coli parameter set may be called
   calibrated.**
6. **σ32 regulation is not in the model**, so every overload statement is open-loop.
7. **Multi-client sign transfer under shared saturated pools is unproven** (EQUATIONS A8, Result 4).
8. **Growth closure μ(x) unspecified.** Under strictly growth-proportional synthesis, slow growth is
   never viable (EQUATIONS A3.2), so "graceful decline to μ → 0" is not a prediction. Basal
   synthesis must be measured, or bounded, before any slow-growth claim.

Drafting rule: the manuscript waits until at least one of P1–P3 has been run, and the novelty
claims in PRIOR_ART §"what this work could newly test" have had a targeted citation search.
