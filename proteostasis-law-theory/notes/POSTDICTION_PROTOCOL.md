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

**Corollary — quote the number in the form the paper states it.** The Lindner
abstract says ">30% of the loss of reproductive ability (aging) in these cells
**relative to the new-pole progeny**". That is a RATIO BETWEEN TWO LINEAGES, not
a loss against an unburdened reference. A post-diction that computes the absolute
quantity is answering a different question, however close the number looks.

## 2. Record the predicted BAND and the parameter regime BEFORE the run

In a commit that precedes the results commit. Two separate commits, in that
order, so the ordering is in the history rather than in a claim about intent.

**State the upper edge, not only the lower one.** A bound reported as ">30%" does
not by itself exclude 95%. If the analysis is going to treat 95% as too severe,
that judgment has to be declared in advance, because after the run it is
indistinguishable from moving the goalposts. This rule exists because D026 broke
it in the direction of excess severity — see §5.

## 3. Record every regime tried, including abandoned ones

The count of regimes tried is part of the result. A match found in the fourth
growth law examined is weaker evidence than a match found in the first, and the
reader can only apply that discount if the denominator is reported.

## 4. A match obtained under a regime the project has already rejected on other
   grounds counts as a FAILURE

Not a partial success, not a promising direction. If the only setting that
reproduces the observation is one the project has independently ruled out, the
model does not explain the observation — it explains it away.

## 5. Retro-applied to the Lindner entry (D026)

Two corrections to the first pass, both recorded rather than quietly fixed:

- **The first pass used `D.Growth(0.04)`** — constant dilution, `k_mu = inf`.
  That produced a clean-looking quantitative match and it was the wrong regime,
  because constant dilution means growth cannot respond to burden and therefore
  predicts zero reproductive loss by construction. It was caught mid-analysis,
  not by a check. Rule 4 exists because of it.

- **D026's magnitude verdict was not licensed by the number.** It reported the
  hyperbolic-feedback loss as "48-95% (median 51%), more severe than reported",
  treating ">30%" as though larger values disconfirmed it. ">30%" is a LOWER
  bound; 51% is consistent with it unless an upper edge is declared, and D026
  declared none. The real failure of that post-diction is the one that stands:
  the physiological law was MONOSTABLE in four of six settings, so there was no
  second attractor for a daughter to be rejuvenated into. Rule 2's upper-edge
  requirement exists because of this.

The relative quantity D026 computed, `1 - mu_high/mu_low`, was correct — it is
the old-pole lineage against the new-pole lineage, which is what the paper
reports. That part needs no correction.

---

## Checklist

Before the run:

- [ ] the paper's measured number, quoted verbatim, with DOI
- [ ] the quantity in the model that corresponds to it, and why it corresponds
- [ ] the match band, **both edges**, with the reasoning for each
- [ ] the parameter regime required, and which regimes are disqualified in advance
- [ ] what result would falsify, stated as a specific outcome
- [ ] committed

After the run:

- [ ] every regime tried, including abandoned ones, with counts
- [ ] the verdict against the pre-stated band, not against a band chosen after
- [ ] if it failed: what the failure names as the next missing mechanism
- [ ] tests pinning the result, so a negative cannot drift into a positive
