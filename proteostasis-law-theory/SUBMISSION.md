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

## What only the author can do

Submission is an outward-facing act and is not automated here:

1. Post to bioRxiv and record the DOI.
2. Add the bioRxiv DOI to the cover letter's preprint sentence.
3. Cut a Zenodo release so the version DOI matches the submitted state, and
   confirm `.zenodo.json` and `CITATION.cff` still agree on the version.
4. Submit through Editorial Manager with the files above.
5. Decide whether to suggest reviewers. None is proposed here; naming real
   people is the author's call, not a thing to generate.
