# Submission package

*An Exact Fold Condition for Mass-Balanced Models of Protein Quality Control*
Kiran Boggavarapu — *Bulletin of Mathematical Biology*

| file | what it is |
|---|---|
| `MANUSCRIPT.pdf` | the paper, 31 pages, figures and tables embedded |
| `SUPPLEMENTARY.pdf` | supplementary, 2 pages, Figures S1–S2 |
| `COVER_LETTER.md` | cover letter |
| `figures/Figure1–5.{pdf,png}` | main figures, for separate upload |
| `figures/FigureS1–S2.{pdf,png}` | supplementary figures |
| `latex_source/` | source, for production |

## latex_source

It compiles as it stands. The layout mirrors the repository because the `.tex`
refers to figures as `../figures/`:

```
latex_source/
  manuscript/   bmb_v5.tex, bmb_v5_supplementary.tex, sn-jnl.cls, sn-mathphys-ay.bst
  figures/      fig1.pdf … figS2.pdf
```

```
cd latex_source/manuscript && pdflatex bmb_v5 && pdflatex bmb_v5
```

Verified by compiling a copy of this directory in isolation, with no access to
the repository: 31 and 2 pages, matching the shipped PDFs, no missing files.

`manuscript_source.md` is the markdown the `.tex` is generated from. It is the
single source; the `.tex` is a build artefact and should not be hand-edited.

## Two things to check before submitting

- The concept DOI in the availability statement currently resolves to the
  **previous** version of this work, under its earlier title. It is a real DOI,
  not a placeholder; the deposit simply has not been updated. See `SUBMISSION.md`
  in the repository root for the sequence that fixes it.
- The cover letter states that a preprint is deposited on bioRxiv. Post it and
  insert the DOI, or cut the sentence.