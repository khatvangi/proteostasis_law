# Submission

Target: *Bulletin of Mathematical Biology*, research article.
Preprint: bioRxiv, posted before or at submission.

## Files

| file | purpose |
|---|---|
| `manuscript/bmb_v5.pdf` | the paper, 31 pages, Springer `sn-jnl` class |
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
- Funding: none received (declared in the manuscript's Statements and Declarations)
- Data: no new experimental data generated

## Gate record

The seven-point read-not-run pass is recorded in `DECISIONS.md` D070, including
what was read as rendered and what was checked mechanically. Four defects were
found and fixed by it; none was a number.

## The deposit — done, with one field outstanding

Deposited 2026-08-06 as **10.5281/zenodo.21830300**, a new version of concept
`10.5281/zenodo.21794565`. Two things verified from the Zenodo API:

- the **concept DOI is preserved**, so the string printed in the manuscript
  resolves to this deposit and will resolve to every later one;
- the archive is `proteostasis-law-theory-2.0.1.zip`, md5
  `5d27fa3496890ddeb64a1fa9488fb434`, byte-for-byte identical to the local build
  from tag `v2.0.1`.

**(!) The record's metadata is still the previous version's.** Zenodo's "New
version" form pre-fills from the prior deposit, and it was published unedited,
so the landing page reads *"An Exact Collapse Threshold for Conserved-Resource
Models…"*, describes the "collapse boundary", and has an empty version field. A
referee following the DOI lands on a record named after the superseded paper.

Files are immutable once published; **metadata is not**. Open
`https://zenodo.org/records/21830300`, click Edit, and set:

- Title: `An Exact Fold Condition for Mass-Balanced Models of Protein Quality Control: analysis code`
- Version: `2.0.1`
- Description: the fold-condition text from `.zenodo.json`

No new version is needed for this. Do it before submitting.

**Do not enable the Zenodo GitHub integration.** This concept was created by
manual upload; enabling the integration would mint a SECOND record with its own
concept DOI, competing with the one printed in the paper.

## Done

1. ~~Push `master`~~ — pushed; the remote carries the reviewed HEAD.
2. ~~Bump the version~~ — see item 4. `1.2.0` was not available: that tag points
   at `0550f6ff`, which predates Phase A, so a Release cut from it would archive
   a tree older than the one Zenodo already holds.
3. ~~Tag~~ — the release tag is **`v2.0.1`**, annotated on the reviewed HEAD and
   pushed. `v2.0.0` remains at `da0a74e` and **must not be released**: it
   predates the three external review passes, so a deposit from it would archive
   the 356-word abstract, no Declarations section, `(G1)–(G4)`, the withdrawn
   ceiling-factor contrast, and four citations since corrected. It was never
   deposited, so it is inert; it is left in place rather than moved because a
   pushed ref that changes meaning is the lineage split this repo keeps finding.
4. ~~Version~~ — `2.0.1` in `.zenodo.json` and `CITATION.cff`, agreeing.

## (!) Stop: the Zenodo integration is off

`GET /repos/khatvangi/proteostasis-law-theory/hooks` returns `[]`. There is no
Zenodo webhook on this repository, and the only GitHub Release that exists is
`v1.1.0` — the deposit `10.5281/zenodo.21798343` that `CITATION.cff` records.

**Zenodo archives releases published after a repository is enabled; it does not
reach back for earlier ones.** Publishing against a dead webhook mints nothing
and consumes the version number, leaving the concept DOI still resolving to the
superseded record. No Release has been created; enable the integration first.

## What only the author can do, in this order

1. **Enable the repository on zenodo.org** (Account → GitHub → toggle
   `khatvangi/proteostasis-law-theory` on). This is an authenticated web action.
   Confirm afterwards that a webhook appears:
   `gh api repos/khatvangi/proteostasis-law-theory/hooks` must be non-empty.
2. **Publish the GitHub Release for the `v2.0.1` tag** — a tag alone
   never triggers a deposit. The Release picks up `.zenodo.json`, so the
   corrected title, description and version land automatically.
3. **Re-resolve the concept DOI and READ the record**, before anything else.
   It must show a title containing "Fold Condition", a description with no
   "collapse boundary", and version `2.0.1`. Two things to watch:
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
