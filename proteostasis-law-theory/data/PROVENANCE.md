# Provenance of the phase-1 run

Every number in the manuscript traces to one run root. This file records
its identity, because the run root itself is gitignored and the deposit is
not reproducible without knowing which run produced it.

## The run in use

    results/phase1/run_20260731T162946-0500

| experiment | git commit | tree dirty | config sha256 | host | runtime (s) | outputs hashed |
|---|---|---|---|---|---|---|
| A | `73e1c0ab3415` | yes | `58380520700d` | boron | 170.5 | 6 |
| B | `a17dfafdd966` | yes | `86e14f20b7dd` | boron | 489.7 | 5 |
| C | `a17dfafdd966` | yes | `5a116aa37614` | boron | 2126.5 | 2 |
| D | `850726c76afe` | yes | `18280a039eea` | boron | 32538.5 | 3 |

Full per-experiment records, including every output file's SHA-256, are in
`<run>/{A,B,C,D}/provenance.json`.

**The four experiments ran at three different commits, not one.** A at
`73e1c0ab`, B and C at `a17dfafd`, D at `850726c7`. The two ensembles the paper
leans on most -- the load grid (B) and the kinetic box (C) -- share a commit,
which is the case that matters for comparing them; but this is not a
single-commit run and must not be described as one. `OPS_SUBMISSION.md` records
a single "canonical repository HEAD", and that is true only of experiment A.

**The working tree was dirty at launch in all four.** The commit identifies the
baseline; it does not fully determine the code that ran. Both facts are
limitations of this deposit and are stated rather than glossed.

## What OPS_SUBMISSION.md is, and is not

`OPS_SUBMISSION.md` documents an EARLIER submission the same afternoon
(16:12:20) whose outputs went to `results/phase1/{A,B,C,D}`. Those four
directories were **empty** -- the launch was reaped by the execution sandbox and
produced no scientific output -- and were quarantined at 16:28:56 into
`results/phase1/quarantine_preclean_20260731T162755-0500`, moved rather than
deleted. The run that produced the data was launched at 16:29:46.

So `OPS_SUBMISSION.md` is a record of scheduler probing, environment capture and
a failed launch. It is kept for that history. **It is not the provenance of any
number in the manuscript**, and reading it as such would attribute the results to
a run that generated nothing.
