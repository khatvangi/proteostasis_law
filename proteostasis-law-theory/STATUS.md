# Status

The paper is complete and is the artifact. This file says where things are and
carries the caveats that have to stay written down. It keeps no second copy of
the paper's results.

## Canonical

| what | where |
|---|---|
| manuscript source | `manuscript/MANUSCRIPT_BMB_v5.md` |
| built paper | `manuscript/bmb_v5.pdf`, `bmb_v5_supplementary.pdf` |
| build | `python scripts/manuscript/to_latex.py` |
| decisions, in order | `DECISIONS.md` (D001–D071) |
| provenance of the ensembles | `data/PROVENANCE.md` |
| target | *Bulletin of Mathematical Biology* |

## Why the results are not restated here

Every quantity in the paper is recomputed by a script in `scripts/` from
`data/computed/` or the phase-1 run root, and asserted by the test suite. This
file used to restate those quantities, and drifted: it carried saturation
fractions, a shortfall split and a ceiling factor that no longer matched the
manuscript. D047 identified the general form of that failure — **two ranges for
one quantity will always let prose land in the wrong place** — after it put the
paper's only falsifiable prediction out by fivefold. A second copy of the
numbers is that fault with more steps, so this file no longer keeps one. The
caveats below are the exception, and the reason is given there: they are
answers, not results, and they have no generator.

The Phase 0 to Phase 2 record, including its stale figures, is preserved
verbatim in `STATUS_PHASE0_TO_2.md`. Read it as history, not as a source.

## Tests

```
python -m unittest discover -s tests/theory       # the theorem layer
python -m unittest discover -s tests/manuscript   # build and canonical-file
python -m unittest discover -s tests/phase3       # the numeric body
python -m unittest discover -s tests/phase2
python -m unittest discover -s tests/figures
python -m unittest discover -s tests -p 'test_*.py'
```

The subdirectories are not packages, so `-t .` fails on them; run each with
`-s` as above.

## Pinned caveats

These are kept here, and not moved to the history file, because they are answers
a later session must meet in writing rather than re-derive. Tests assert them.

**Bistability under division.** A bistable window was found at one parameter point,
not surveyed.

**That parameter point is unphysical, in a specific and measured way.** It uses
constant dilution (`k_mu = inf`), under which growth rate cannot respond to
burden, so the same regime predicts **exactly zero** growth-rate loss at any
aggregate load. That contradicts the one dosage-resolved measurement the project
holds (D015: 3.2 % loss at <0.1 % misfolded) and the observation the first
post-diction tried to explain (D026/D028: a 1.2–1.8 % aggregate-attributable
growth deficit — not the ">30 %" an earlier reading of that abstract took it to
be). The regime that produces the bistability is the same regime that gets the
measured quantity wrong. Under the physiological laws it does not survive as
stated — linear arrest gives no bounded high-burden state, hyperbolic feedback
is monostable in four of six settings.
Pinned by test so it cannot be edited away.

**The per-network spread of `phi` is not bounded.** The nested design
was not sized to bound the per-network spread, and the stated position is that
we cannot bound it. Any figure quoted for it is a largest-observed value over
ten draws and understates its population (D038). A referee pushing on it should
meet this sentence, not a number invented on the spot.

## Supersession

Phases A through C (D045–D071) rebuilt the theorem layer and re-derived every
quantity the manuscript reports. Where an earlier decision, note or status
entry disagrees with the manuscript, the manuscript and the generating script
are correct and the earlier text is superseded. In particular:

- `STATUS_PHASE0_TO_2.md` — superseded in full by the manuscript.
- `OPS_SUBMISSION.md` — a record of a failed launch that produced no output;
  it is not the provenance of any number. See `data/PROVENANCE.md`.
- `PHASE2_CLOSURE_FINAL.md` §8 — its open list is closed.
- Any claim that the identity needs (H1b), that the fold state lies in a centre
  manifold, that all four contributions to `G_u` share a sign, or that the
  Section 8.3 requirement sits *above* the measured Δ*rpoH* load — all four were
  wrong and are corrected in D062, D067 and the residual sweep.
