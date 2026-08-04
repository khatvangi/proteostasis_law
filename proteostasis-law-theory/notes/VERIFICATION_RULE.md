# Before pinning a number, state what would make it the wrong number

A test that asserts a value is reproduced says nothing about whether the value
was worth computing. This repository has more than 340 tests and they do their
job, which is to prevent drift. They cannot detect that a quantity is the wrong
quantity, and they did not.

## The evidence

Three of the last four real defects were found by RENDERING or READING, not by
running the suite:

1. **The fold sat off its own nullcline** (D034). Root-finding in `a` at fixed
   `u` loses the curve at its turning point and returns two disconnected pieces.
   Every numerical test passed; the figure showed the fold floating in the gap.
2. **The caption metric read exactly 1.0** (D035). `rel_err` divides by
   `max(|det J|, |cross|)`, and both vanish at a saddle-node, so it is `0/0`
   precisely where a caption would quote it. Nothing asserted the metric was
   meaningful, only that it was reproduced.
3. **Three headline numbers came from a 6 % subsample** (D036), all three
   flattering. No test asks how many states a statistic was computed on. The
   forensics matter: the selection was `sample(n=20, random_state=1)`, a uniform
   draw, so nothing chose those states. Two of the three deviations are coin
   flips -- measured `P(subsample flatters)` of 0.501 and 0.438 over 4000 redraws.
   The third, a maximum, understated with probability 0.939, which is exactly
   `1 - 20/325`. One deterministic bias, two chance deviations, and the
   deterministic one is the general rule below.

The fourth, the Lindner abstract misreading, was found by reading the full text.

## The rule

**Before pinning a number with a test, write down what would make it the wrong
number.** Then check those things, and pin those too where they are checkable.

For a verification statistic the standard questions are:

- [ ] **Does the normalisation vanish where the quantity is evaluated?** Any
      residual divided by something that goes to zero at the point of interest is
      `0/0` exactly where it is quoted.
- [ ] **Does the statistic improve or degrade as it approaches the object being
      verified?** Measure it: correlate the residual against a proximity measure.
      A negative correlation means the number flatters itself.
- [ ] **What population is it computed on, and is that the whole of it?** Record
      the count and whether it is a subsample, at the point of use. Failure mode:
      an INCOMPLETE population. "2884 folds" headed a table computed over 325.
- [ ] **Is it the population the surrounding sentence is about?** A separate
      failure mode, and the one the completeness check cannot catch, because
      nothing is missing. The `s_u` p5–p95 width of 0.890 is correct and complete
      over the regulation experiment's 30 networks; §9 quoted it as a general
      property of the kinetic box, whose 2884 folds give 0.876. Right statistic,
      wrong experiment. The two values agreed, so nothing moved — which is
      precisely why it survived four earlier corrections. Ask both questions: is
      the population whole, and is it the population this sentence claims.
- [ ] **Is it an extremum?** A maximum or minimum over a subsample is optimistic
      BY CONSTRUCTION, not by chance: it understates with probability `1 - k/N`.
      Never subsample an extremum. Medians and correlations on a subsample are
      noisy but unbiased; extrema are biased.
- [ ] **Is it an extremum in disguise?** Two sweeps of every max, min, worst-case
      and "throughout" both missed "the single state bracketed to eigenvalue
      -4.2e-9", which is a MINIMUM written as a definite description. "The single
      state where...", "the one draw that...", "the tightest / closest / last /
      only..." are extrema in prose and carry the same bias. The full population
      was 5x tighter and had three such states, so the uniqueness claim was false
      as well. Grep for the grammar, not only for the word.
- [ ] **Can a committed script recompute it right now?** Not "was it computed
      correctly" but "does the code that produces it exist in the repository". §5's
      normalisation contrast (-0.262 and +0.060) was computed in-session, never
      committed, and does not reproduce under any of 48 definitions (D041). The
      medians beside it in the same table reproduced exactly, which is what hid it:
      a wrong population moves every entry, an ad-hoc computation moves only some.
- [ ] **Does this number appear elsewhere in the repository, and does that copy
      agree?** Grep it. A corrected value in one document and the original in
      another is a lineage split inside a single repository, which is the failure
      that has already cost this project once. Correct or mark superseded at every
      site, and prefer marking in the decision log, which is a record rather than
      a statement of current fact.
- [ ] **Would this check fail if the implementation were wrong?** The self-damage
      identity check carries no truncation term (D027), so it cannot fail unless
      mass balance fails. That is worth knowing and worth saying; it does not make
      the check useless, it makes it a check of something narrower than it looks.

## Pin the property, never the token

A test that asserts a string is ABSENT fails the moment the text becomes honest
about its own history. Section 5 states its superseded values explicitly, as
corrections, which is the right presentation -- and two tests asserting those
strings were absent failed on correct text. String-absence tests actively
penalise the correction discipline this project runs on.

Assert the property the string stood for. Not "1.436e-07 does not appear" but
"1.436e-07 does not appear as a current claim": absent from the results table,
present only inside the sentence that supersedes it.

## Auditing extrema across a document

"Never subsample an extremum" is an audit criterion, not a note about one row.
Enumerate every maximum, minimum, worst-case and "throughout" in the text, and for
each record its population, whether that population is complete, and a p99 beside
it. This paper had eight such quantities and the audit found a fifth affected
number after four were already corrected -- the solver maximum, understated at
6.652e-07 by a 20-state draw against 7.56e-07 over the full 325.

The deeper problem is that an extremum is a weak statistic even when computed
correctly, because it grows with population size. A maximum over 2884 draws will
exceed a maximum over 325 from the same distribution, so two maxima with different
denominators are not comparable, and neither is comparable to a reader's rerun at a
third size. **Where a maximum is doing bounding work, report a high quantile beside
it or instead of it.** A p99 is stable under resampling and carries the same claim.

## A limitation is an observation, not a licence

The worst defect in the figure work was not introduced by the figure work. A
standing limitation in `theory/FOLD_THEOREM.md` read:

> Some draws collapse at `s_a` near 0.003, where aggregation is so fast the
> low-burden branch barely exists. **These should be screened rather than
> averaged into a median.**

The first sentence is an observation and may well be right. The second is an
instruction, and it authorises an exclusion that no evidence supports. Acting on
it produced a 5x divergence between Figure 2 and §6 before the screen was tested
and abandoned (D039). The limitation had sat in the repository for months as a
**pre-authorised free parameter**, waiting for the first person to act on it.

**A suspicion about data quality does not license a threshold.** The two claims
must be stated separately, because they have different evidential status and the
first can survive the second's failure — as it did here.

The audit form: enumerate every limitation, caveat and open-question note in the
repository, and classify each as

- **an observation** — "X is true of the data" — safe;
- **a restriction** — "do not quote X without Y" — safe, and the opposite
  direction, since it adds a requirement rather than removing data;
- **a direction** — "X should be screened / dropped / excluded / substituted" —
  a free parameter with prior authorisation. Rewrite as an observation, and if
  the operation is genuinely wanted, make it a decision entry with a test.

Sweeping `theory/`, `notes/` and `STATUS.md` on 2026-08-04 found exactly one in
the third form: the bullet above. Near misses that are safe and were left alone:
"should not be quoted as a result without a larger sample" and "should not be
quoted without stating which law produced them" are restrictions; "every `phi` in
Consequences 2 and 3 should be read as `phi_enz` at `delta = 0`" is a direction in
grammar but an algebraic identity in content, so it changes no number; and the
`should`s in `theory/PREDICTIONS.md` are hypotheses about the world, not
instructions about the data.

## Where this applies

Everywhere a number enters the manuscript, not only in post-diction. The
post-diction protocol governs how an empirical target is chosen; this governs
whether a computed quantity means what it is being asked to mean. They are
different failures and both have occurred here.
