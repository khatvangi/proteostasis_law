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

## Source snapshots (D2)

`scripts/analysis/snapshot_phase1_sources.py` writes the tracked source of each
experiment's driver and the `proteostasis` package to
`data/provenance_snapshots/{A,B,C,D}/`, with a SHA-256 per file in `index.json`.
This replaces a twelve-hex commit prefix the reader would have to resolve with a
fixed artifact in the deposit.

**The launch-time diff is not among them and cannot be.** `provenance.json`
records `dirty: true` for all four experiments and stores no patch. A diff
between the recorded commit and the next one would be a plausible-looking guess
rather than the code that ran, and is deliberately not produced.

**Experiment A's recorded commit contains no code.** `73e1c0ab` holds 17 files,
every one of them markdown or `.gitignore`; no Python existed in the repository
at that point, so A ran entirely from an uncommitted tree and its recorded
commit is not a baseline in any useful sense. Its snapshot is taken from
`a17dfafd`, the earliest commit containing `run_experiment_a.py`, and
`index.json` marks this with `recorded_commit_held_no_source`. The table above
lists A's commit as `73e1c0ab` with `tree dirty: yes`, which is true and
understates the problem; this entry is the correction.

**Where the gap sits.** Phase 1 produced the raw ensembles. Every derived number
in the manuscript is recomputed from those raw outputs by tracked code under the
test suite, so the unrecorded state affects the generation of the ensembles and
not the analysis built on them. Re-running the 5000-draw box would not close it
either, since a new run would carry a new provenance rather than recover the old.
