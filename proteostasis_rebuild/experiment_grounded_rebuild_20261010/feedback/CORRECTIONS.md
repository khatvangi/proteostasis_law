# Corrections, 2026-10-10

## Result

**Tests: 21 run, 21 pass, 0 fail.** Command: `cd feedback && python -m unittest -v test_feedback`.
Run on Python 3.12.11 in 126.5 s; the output is in `test_output.txt`. Before these corrections the
suite had 12 tests, all passing; 9 were added.

All parameters are dimensionless placeholders. Nothing here is an E. coli estimate.

## Where the corrections came from

The source is the independent review `../feedback_review/REVISIONS_REQUIRED.md`, with its script
`check_feedback_claims.py` and its output `check_output.txt`. It was written by a separate agent
that did not edit these files. Nothing in `feedback_review/` or `../independent_review/` was
modified.

I re-ran the review script against the corrected `feedback.py` (exit 0, about 8 min; output in
`review_recheck_output.txt`). It is identical to the original `check_output.txt` except section 3,
where the code's h′(0) now equals the true value for n = 1, 2 and 0.5.

The review's claims were not taken on trust. Each one used here was reproduced by a test in this
folder, or by the re-run above.

## What was wrong, and what replaces it

| # | Was | Now | Checked by |
|---|---|---|---|
| 1 | "G < 1 is necessary but not sufficient; a crossing can exist if e² > a²+μ²" (MINIMAL_DYNAMICS §3) | **Wrong.** The linearised loop is stable at every delay **iff G ≤ 1**. G < 1 forces e < μ, so both crossing coefficients are positive. At G = 1 the only candidate is ω = 0, and f(0) = aμ(1+G) ≠ 0. G > 1 gives exactly one crossing frequency, every crossing destabilises, so it is stable iff τ < τ*. Full proof in §3. Stated as linear and local. | `TestDelayTheorem` (5,000 random loops; equality case; nonconstant-history integration at G = 0.5, 1, 2.5) |
| 2 | "e < μ" or "e² < a²+μ²" read as safe | Neither suffices. Witness: a = b = 1, c = 2, e = 0.5, μ = 1 gives G = 2.5, τ* = 1.267 | `test_strong_gain_is_unstable_at_every_delay_beyond_tau_star` (decays at 0.9τ*, grows at 1.1, 2 and 3τ*) |
| 3 | "The pool can never exceed Cmax"; "at any state served ≤ k·Cmax" (§4); `findCapMax` "largest pool feedback can ever hold" | The band is **forward-invariant and attracting**: C(t) − Cmax ≤ (C(0) − Cmax)e^(−μt). The bound holds only if C(0) ≤ Cmax and growth is constant. Two counterexamples follow. | `TestCeiling` (3 tests) |
| 3a | | *Initial condition:* C(0) = 3Cmax stays above Cmax until t ≈ 3.3, and served flux reaches 6.3 > k·Cmax = 4.2 > L | `TestCeiling` |
| 3b | | *Nonconstant growth:* a μ upshift (0.25 → 0.5) starts the pool above its new ceiling | `TestCeiling` |
| 3c | | In fixed mode Cmax → ∞ as μ → 0 | `TestCeiling` |
| 4 | POSITIVE_INVARIANCE Consequence 1: trajectories outside "enter the C-band" | Kept, but now stated with a finite-time proof (h bounded below 1 while U is bounded). Also says plainly that the box is not an absolute bound. | `TestCeiling` |
| 5 | `findSigmaActivity` returned h′ = 0 at F = 0 for every n | It now returns −1/Ks for n = 1, 0 for n > 1, and −inf for n < 1. Harmless for earlier results, since F > 0 whenever C > 0. | `test_sigma_derivative_at_zero_free_chaperone`. Restoring the old line makes this test fail (mutation check). |
| 6 | Collapsed state under growth arrest presented without its condition (§5) | Now stated as **conditional on synthesis ∝ growth**. In the placeholder case, synthesis in the collapsed state is ~1000× below the served state. Fixed mode has no ceiling instead. Not measured for E. coli. | `test_collapse_rests_on_synthesis_proportional_to_growth` |
| 7 | "Tomoyasu 1998 shows high gain" (§3, EXPERIMENTAL_CONTROL) | Tomoyasu shows sensitivity, not G > 1. So whether any delay limit exists is unknown. | text only |
| 8 | Served flux read as folded or native output (§4, §7, PRIOR_ART_MAP) | "Served" is chaperone turnover. The model has no native-product or aggregate state. | text only |
| 9 | "Raising δ lowers served flux (tested)" | Now a proof (constant μ, fixed mode), plus the test | existing test; review re-run 0/100 exceptions |
| 10 | Assumptions mixed into the text | An assumption-vs-LIT table was added to MINIMAL_DYNAMICS. Spare capacity = k·F is marked as true by construction. `findLoopTerms` now refuses growth mode. | `test_collapse_rests_on_synthesis_proportional_to_growth` (raises) |

## One place where I did not follow the review

The review asked for the wording "the U-diagonal of the growth-mode Jacobian turns negative". At the
saddle of the placeholder case that entry is **+0.0095**, so it turns positive. The original "a turns
negative" was right in content, since the entry is −a. The text now says "∂(dU/dt)/∂U becomes
positive".

Also, the review suggested forbidding n < 1 in `findSigmaActivity`. I return the true value (−inf)
instead. Forbidding it would make the review's own script crash, and the value is correct either way.

## Scope

None of this is a statement about E. coli:

- **Delay result.** Linear and local. The Hopf direction (whether the oscillation beyond τ* is small
  and stable or jumps to a large one) was not checked. Distributed delays were not checked, and
  neither was the growth-mode model with delay.
- **Global convergence.** Shown only for τ = 0, fixed production and constant μ.
- **Parameters.** No parameter value is calibrated. G, τ, k, Cmax and the β's are unmeasured.

## Files changed

- `feedback.py`
- `test_feedback.py`
- `MINIMAL_DYNAMICS.md` (rewritten)
- `POSITIVE_INVARIANCE.md`
- `EXPERIMENTAL_CONTROL.md` (one sentence)
- `PRIOR_ART_MAP.md` (two cells)

New files: `CORRECTIONS.md`, `test_output.txt`, `review_recheck_output.txt`.

Base commit: 675cbfa. These changes are not committed.
