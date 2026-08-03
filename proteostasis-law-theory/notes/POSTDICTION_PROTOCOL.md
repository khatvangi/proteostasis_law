# Post-diction protocol

A post-diction is an attempt to explain an already-published observation with a
model that was not built to explain it. It is the cheapest evidence this project
can generate, and the easiest to fake without noticing.

These rules were not designed in advance. They are what the first attempt
(D026, Lindner et al. 2008) discovered by getting them wrong, and they are
written down so the discovery does not stay tacit.

---

## 1. Find the measured NUMBER in the paper, not the phenomenon

"Bistable, therefore rejuvenation" would have passed. ">30% loss of reproductive
ability" is what caught the wrong regime: the regime that produced the
bistability was constant dilution, which predicts a loss of exactly zero.

A phenomenon can be matched by any model loose enough. A number cannot.

**Corollary A — quote the number in the form the paper states it.** The Lindner
abstract says ">30% of the loss of reproductive ability (aging) in these cells
**relative to the new-pole progeny**". That is a RATIO BETWEEN TWO LINEAGES, not
a loss against an unburdened reference. A post-diction that computes the absolute
quantity is answering a different question, however close the number looks.

**Corollary B — establish what the percentage is a percentage OF, from the full
text, before using it.** This is the rule that cost the most. ">30% **of** the
loss of reproductive ability" was read through two entire analysis cycles as
">30% loss of reproductive ability". It is not. The full text gives the aging
effect as `[Δ(GR_old − GR_new)]mean/GR_mean = −3.95 ± 0.5%` and the
aggregate-attributable share as `Agg/(Agg + Pole) ≈ 30–40%`. The measured
quantity is therefore **≈1.2–1.8%**, not 30%. The abstract alone could not settle
it and was not enough. **An abstract is not a source for a number.**

## 2. Record the predicted BAND and the parameter regime BEFORE the run

Written down before the run. Where it sits in git history does not matter — the
band being fixed in advance is the substance; commit ordering was ceremony and is
not required.

**State the upper edge, not only the lower one**, and derive both from data
wherever the source permits. A bound reported as ">30%" does not by itself
exclude 95%; if the analysis is going to treat 95% as too severe, that has to be
settled in advance. Prefer a data-derived edge to a judgment about how authors
phrase bounds — in the Lindner case the full text bounded the share at ≈40%,
which the abstract did not.

## 3. Record every regime tried, including abandoned ones

The count of regimes tried is part of the result. A match found in the fourth
growth law examined is weaker evidence than a match found in the first, and the
reader can only apply that discount if the denominator is reported.

## 4. A match obtained under a regime the project has already rejected on other
   grounds counts as a FAILURE

Not a partial success, not a promising direction. If the only setting that
reproduces the observation is one the project has independently ruled out, the
model does not explain the observation — it explains it away.

## 6. The pass criterion must require the MECHANISM to be on

Include the mechanism-off control in the grid, and require the in-band cell to
have the mechanism switched on. A control cell cannot demonstrate a mechanism no
matter which band it lands in.

This is not pedantry. D028 fixed the band, the growth law and the falsifier, and
still admitted a "pass" that was carried entirely by `k_seq = 0` cells — the
two-state model, reproducing an earlier result, with sequestration doing nothing
(D029). Report both criteria: the literal preregistered one as the audit trail,
and the mechanism-on one as the scientific verdict.

## 5. Retro-applied to the Lindner entry (D026)

Two corrections to the first pass, both recorded rather than quietly fixed:

- **The first pass used `D.Growth(0.04)`** — constant dilution, `k_mu = inf`.
  That produced a clean-looking quantitative match and it was the wrong regime,
  because constant dilution means growth cannot respond to burden and therefore
  predicts zero reproductive loss by construction. It was caught mid-analysis,
  not by a check. Rule 4 exists because of it.

- **D026's magnitude verdict was CORRECT, and an earlier version of this section
  wrongly withdrew it.** That earlier version argued D026 had no licence to call
  the hyperbolic-feedback loss of "48-95% (median 51%)" too severe, since ">30%"
  is a lower bound with no declared ceiling. That argument rested on reading the
  abstract as ">30% loss". The full text says the measured quantity is ≈30-40% OF
  a 3.95 ± 0.5% aging effect, i.e. **≈1.2-1.8%**. A predicted 48-95% is therefore
  too severe by a factor of roughly **30 to 60**, and D026's instinct was right
  even though its stated reasoning was imprecise.

  The correction of the correction is the lesson: a wrong reading of a number
  survived being questioned, being defended, and being written into a protocol as
  a worked example. Only the full text stopped it. See Corollary B.

  D026's other ground stands independently and was never in doubt: the
  physiological law was MONOSTABLE in four of six settings, so there was no
  second attractor for a daughter to be rejuvenated into.

The relative quantity D026 computed, `1 - mu_high/mu_low`, was correct — it is
the old-pole lineage against the new-pole lineage, which is what the paper
reports. That part needs no correction.

---

## Checklist

Before the run:

- [ ] the paper's measured number, quoted verbatim, with DOI
- [ ] **what that number is a percentage OF**, established from the full text
- [ ] the quantity in the model that corresponds to it, and why it corresponds
- [ ] the match band, **both edges**, derived from data where the source permits
- [ ] the mechanism-OFF control is in the grid, and the pass criterion requires
      the mechanism ON
- [ ] the parameter regime required, and which regimes are disqualified in advance
- [ ] what result would falsify, stated as a specific outcome
- [ ] committed

After the run:

- [ ] every regime tried, including abandoned ones, with counts
- [ ] the verdict against the pre-stated band, not against a band chosen after
- [ ] if it failed: what the failure names as the next missing mechanism
- [ ] tests pinning the result, so a negative cannot drift into a positive
