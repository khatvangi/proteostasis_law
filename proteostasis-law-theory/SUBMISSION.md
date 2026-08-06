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

## What only the author can do, in this order

1. **Push `master`.** The remote is 5 commits behind: `cb2044e`, `3f6fa95`,
   `4993870`, `085f097`, `540f2a3` — every phase the paper depends on.
2. **Decide the version.** `.zenodo.json` and `CITATION.cff` both say `1.2.0`,
   which is the tag the remote already has. Phases A–C changed theorem
   statements and retracted claims, so the next number is a judgement call;
   bump both files together or the canonical-file test fails.
3. **Publish a GitHub *Release*, not just a tag.** Tags `v1.0.0`, `v1.1.0` and
   `v1.2.0` exist on the remote, but no GitHub Release does, which is why Zenodo
   never advanced past the record it minted from `v1.0.0`. Only a Release
   triggers the deposit. The Release picks up `.zenodo.json`, so the corrected
   title and description land automatically.
4. **Re-resolve the concept DOI and read it.** The record must show the title
   containing "Fold Condition", a description with no "collapse boundary", and
   the version from step 2. This check cannot be automated in the suite —
   `test_canonical_file` asserts the *local* metadata and is blind to the
   deposit, which is how this survived.
5. Post to bioRxiv, record the DOI, and put it in the cover letter's preprint
   sentence — which currently asserts a deposit that does not exist yet.
6. Submit through Editorial Manager with the files above.
7. Decide whether to suggest reviewers. None is proposed here; naming real
   people is the author's call, not a thing to generate.
