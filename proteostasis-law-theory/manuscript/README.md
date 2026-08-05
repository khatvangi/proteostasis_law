# manuscript/

**Canonical: [`MANUSCRIPT_BMB_v5.md`](MANUSCRIPT_BMB_v5.md).** It is the only
manuscript in use. `scripts/manuscript/to_latex.py` reads it and writes
`bmb_v5.{tex,pdf}` and `bmb_v5_supplementary.{tex,pdf}`; those are build
artefacts and must never be edited.

| file | status |
|---|---|
| `MANUSCRIPT_BMB_v5.md` | **canonical** — target: Bulletin of Mathematical Biology |
| `bmb_v4.md` | superseded; the v4 draft, retained for provenance |
| `SECTION_8_4_LINEAGE_PREDICTION.md` | superseded; folded into §8.3 |
| `COLLAPSE_BOUNDARY.md` | superseded; earlier phase 3 write-up |
| `MANUSCRIPT.md` | superseded; phase 0 framework paper |

Numbers in the canonical manuscript are recomputable from `scripts/phase3/` and
`scripts/analysis/`; the decision log in `DECISIONS.md` records what each claim
rests on and what was withdrawn.

This file has named the wrong canonical manuscript three times. A test now
asserts that both READMEs name the file the converter actually reads, so a
lineage split fails the suite instead of sitting in the documentation.
