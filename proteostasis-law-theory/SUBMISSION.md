# Submission

Target: *Bulletin of Mathematical Biology*, research article.
Preprint: bioRxiv, posted before or at submission.

## Files

| file | purpose |
|---|---|
| `manuscript/bmb_v5.pdf` | the paper, 30 pages, Springer `sn-jnl` class |
| `manuscript/bmb_v5_supplementary.pdf` | supplementary, 2 pages |
| `manuscript/bmb_v5.tex` | LaTeX source, for production |
| `manuscript/springer/sn-jnl.cls`, `sn-mathphys-ay.bst` | class and bibliography style |
| `figures/fig{1..5}.pdf`, `figS{1,2}.pdf` | figures, vector |
| `COVER_LETTER.md` | cover letter |

Rebuild with `python scripts/manuscript/to_latex.py`. The build fails closed on an
undeclared display, an unbalanced figure float, a table left as a `longtable`,
and any count that departs from `EXPECTED`.

## Metadata

- Title: An Exact Fold Condition for Mass-Balanced Models of Protein Quality Control
- Author: Kiran Boggavarapu, ORCID 0000-0003-0751-6459
- Affiliation: Department of Chemistry and Physics, McNeese State University, Lake Charles, LA 70609, USA
- Corresponding: kiran@mcneese.edu
- MSC: 92C40 (primary); 37G10, 34C23, 92C42 (secondary)
- Keywords: proteostasis; saddle-node bifurcation; conserved resource; clearance network; growth dilution; Hopf bifurcation
- Code: https://github.com/khatvangi/proteostasis-law-theory, MIT
- Archive: https://doi.org/10.5281/zenodo.21794565 (concept DOI)
- Competing interests: none
- Funding: to be stated by the author
- Data: no new experimental data generated

## Gate record

The seven-point read-not-run pass is recorded in `DECISIONS.md` D070, including
what was read as rendered and what was checked mechanically. Four defects were
found and fixed by it; none was a number.

## (!) The availability statement does not yet point at this paper

Checked 2026-08-05 by resolving the DOI and the URL in the rendered PDF. Both
are real and neither is a placeholder — and both currently resolve to the
**previous** version of this work:

| what the PDF says | what it resolves to today |
|---|---|
| `doi.org/10.5281/zenodo.21794565` | a record titled *"An Exact **Collapse Threshold** for **Conserved-Resource** Models…"*, described as the "**collapse boundary**", archiving GitHub tag `v1.0.0` |
| `github.com/khatvangi/proteostasis-law-theory` | public, MIT, at `aa97b09` — **5 commits behind local**, README still naming `bmb_v4.md` canonical |

Both are titles and claims this paper explicitly supersedes, and the archived
code predates Phases A–C, so it cannot recompute a single number in the current
manuscript. `.zenodo.json` and `CITATION.cff` locally carry the correct title
and description; they have simply never been deposited. **Do not submit before
the steps below.**

## Done

1. ~~Push `master`~~ — pushed, remote at `da0a74e`.
2. ~~Bump the version~~ — `2.0.0` in `.zenodo.json` and `CITATION.cff`, agreeing.
   `1.2.0` was not available: that tag points at `0550f6ff`, which predates
   Phase A, so a Release cut from it would archive a tree older than the one
   Zenodo already holds.
3. ~~Tag~~ — `v2.0.0` annotated on `da0a74e` and pushed. Not on `1ccbc64`: that
   tree still declares `1.2.0`, so a deposit from it would carry the wrong
   version.

## (!) Stop: the Zenodo integration is off

`GET /repos/khatvangi/proteostasis-law-theory/hooks` returns `[]`. There is no
Zenodo webhook on this repository, and the only GitHub Release that exists is
`v1.1.0` — the deposit `10.5281/zenodo.21798343` that `CITATION.cff` records.

**Zenodo archives releases published after a repository is enabled; it does not
reach back for earlier ones.** Publishing `v2.0.0` now would therefore mint
nothing and consume the version number, leaving the concept DOI still resolving
to the superseded record. The Release is deliberately not created.

## What only the author can do, in this order

1. **Enable the repository on zenodo.org** (Account → GitHub → toggle
   `khatvangi/proteostasis-law-theory` on). This is an authenticated web action.
   Confirm afterwards that a webhook appears:
   `gh api repos/khatvangi/proteostasis-law-theory/hooks` must be non-empty.
2. **Publish the GitHub Release for the existing `v2.0.0` tag** — a tag alone
   never triggers a deposit. The Release picks up `.zenodo.json`, so the
   corrected title, description and version land automatically.
3. **Re-resolve the concept DOI and READ the record**, before anything else.
   It must show a title containing "Fold Condition", a description with no
   "collapse boundary", and version `2.0.0`. Two things to watch:
   - a published Zenodo version cannot be withdrawn, only superseded, so if the
     record comes through wrong the remedy is another version, not an edit;
   - if re-enabling mints a **new concept DOI** rather than continuing
     `10.5281/zenodo.21794565`, the DOI printed in the manuscript is then the
     wrong one and the availability sentence must be corrected and the PDF
     rebuilt before submission.
   This check cannot live in the suite. Every test here verifies local state;
   nothing verifies the deposit, the remote, or any other external artifact,
   and nothing can, because those are network acts. It is a permanent manual
   gate item.
4. Add the minted version DOI to `CITATION.cff` under `identifiers`.
5. Post to bioRxiv, record the DOI, and put it in the cover letter's preprint
   sentence — which currently asserts a deposit that does not exist yet.
6. Submit through Editorial Manager with the files above.
7. Decide whether to suggest reviewers. None is proposed here; naming real
   people is the author's call, not a thing to generate.
