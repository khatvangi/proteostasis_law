# Status

The paper is complete and is the artifact. This file says where things are; it
deliberately carries no numbers.

## Canonical

| what | where |
|---|---|
| manuscript source | `manuscript/MANUSCRIPT_BMB_v5.md` |
| built paper | `manuscript/bmb_v5.pdf`, `bmb_v5_supplementary.pdf` |
| build | `python scripts/manuscript/to_latex.py` |
| decisions, in order | `DECISIONS.md` (D001–D068) |
| provenance of the ensembles | `data/PROVENANCE.md` |
| target | *Bulletin of Mathematical Biology* |

## Why there are no numbers here

Every quantity in the paper is recomputed by a script in `scripts/` from
`data/computed/` or the phase-1 run root, and asserted by the test suite. This
file used to restate those quantities, and drifted: it carried saturation
fractions, a shortfall split and a ceiling factor that no longer matched the
manuscript. D047 identified the general form of that failure — **two ranges for
one quantity will always let prose land in the wrong place** — after it put the
paper's only falsifiable prediction out by fivefold. A second copy of the
numbers is that fault with more steps, so this file no longer keeps one.

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

## Supersession

Phases A through C (D045–D068) rebuilt the theorem layer and re-derived every
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
