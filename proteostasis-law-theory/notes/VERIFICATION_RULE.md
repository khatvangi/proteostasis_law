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
      the count and whether it is a subsample, at the point of use.
- [ ] **Is the population the one the surrounding text describes?** Two
      populations under one sentence is how "2884 folds" came to head a table
      computed over 325.
- [ ] **Is it an extremum?** A maximum or minimum over a subsample is optimistic
      BY CONSTRUCTION, not by chance: it understates with probability `1 - k/N`.
      Never subsample an extremum. Medians and correlations on a subsample are
      noisy but unbiased; extrema are biased.
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

## Where this applies

Everywhere a number enters the manuscript, not only in post-diction. The
post-diction protocol governs how an empirical target is chosen; this governs
whether a computed quantity means what it is being asked to mean. They are
different failures and both have occurred here.
