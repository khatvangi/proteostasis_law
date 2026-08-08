"""Assemble `FINAL_SUBMISSION/` from the canonical build, and assert it.

The package deliberately carries each figure TWICE: once as `figures/FigureN.pdf`
for upload to a portal that wants them separately, and once as
`latex_source/figures/figN.pdf` where the `.tex` expects to find it. That is the
right shape for a submission bundle and it is a divergence surface -- regenerate
one figure and one copy goes stale in silence, which is the lineage split this
project keeps catching, relocated into the delivery folder.

So the copies are not made and forgotten. Every pair is checksummed after the
copy, and every shipped artifact is checksummed against the canonical file it
came from. A mismatch fails the build rather than shipping.

The first version of this package was assembled by hand and did not compile: the
`.tex` asks for `../figures/fig1.pdf` and the folder held `Figure1.pdf`, with
the class file in a subdirectory rather than beside the `.tex`. The layout below
mirrors the repository for exactly that reason, and `--compile` proves it by
building a copy in a scratch directory with no access to the repo.
"""

from __future__ import annotations

import hashlib
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
OUT = REPO_ROOT / "FINAL_SUBMISSION"

MAIN_FIGS = ["fig1", "fig2", "fig3", "fig4", "fig5"]
SUPP_FIGS = ["figS1", "figS2"]


def _sha(p: Path) -> str:
    return hashlib.sha256(p.read_bytes()).hexdigest()


def _copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def build() -> dict:
    if OUT.exists():
        shutil.rmtree(OUT)
    (OUT / "figures").mkdir(parents=True)
    (OUT / "latex_source" / "manuscript").mkdir(parents=True)
    (OUT / "latex_source" / "figures").mkdir(parents=True)

    m, f = REPO_ROOT / "manuscript", REPO_ROOT / "figures"
    shipped = {
        OUT / "MANUSCRIPT.pdf": m / "bmb_v5.pdf",
        OUT / "SUPPLEMENTARY.pdf": m / "bmb_v5_supplementary.pdf",
        OUT / "COVER_LETTER.md": REPO_ROOT / "COVER_LETTER.md",
        OUT / "latex_source" / "manuscript" / "bmb_v5.tex": m / "bmb_v5.tex",
        OUT / "latex_source" / "manuscript" / "bmb_v5_supplementary.tex":
            m / "bmb_v5_supplementary.tex",
        OUT / "latex_source" / "manuscript" / "sn-jnl.cls":
            m / "springer" / "sn-jnl.cls",
        OUT / "latex_source" / "manuscript" / "sn-mathphys-ay.bst":
            m / "springer" / "sn-mathphys-ay.bst",
        OUT / "latex_source" / "manuscript_source.md":
            m / "MANUSCRIPT_BMB_v5.md",
    }
    for dst, src in shipped.items():
        _copy(src, dst)

    # Editorial Manager takes the cover letter as a document, not markdown.
    # Generated rather than kept by hand so it cannot drift from the .md, which
    # is the source; a hand-maintained .docx is a second lineage for the one
    # part of the package a human actually reads first.
    docx = OUT / "COVER_LETTER.docx"
    r = subprocess.run(["pandoc", str(REPO_ROOT / "COVER_LETTER.md"),
                        "-o", str(docx), "--from=markdown", "--to=docx"],
                       capture_output=True, text=True)
    if r.returncode != 0 or not docx.exists():
        raise SystemExit(f"cover letter did not convert to docx:\n{r.stderr}")

    # The abstract as a standalone page, for portals that ask for it as a file.
    # Page 1 of the built PDF IS the abstract page -- title, byline, abstract,
    # keywords, MSC -- so it is extracted rather than re-rendered. A separately
    # typeset abstract would be a second lineage for the paper's most-read
    # paragraph, and would drift the first time the abstract was edited.
    r = subprocess.run(["pdftotext", "-layout", "-f", "2", "-l", "2",
                        str(m / "bmb_v5.pdf"), "-"],
                       capture_output=True, text=True)
    if "Introduction" not in r.stdout:
        raise SystemExit("page 2 is not the Introduction: the abstract no longer "
                         "ends on page 1, so ABSTRACT.pdf would be truncated")
    tmp = OUT / "ABSTRACT_%d.pdf"
    r = subprocess.run(["pdfseparate", "-f", "1", "-l", "1",
                        str(m / "bmb_v5.pdf"), str(tmp)],
                       capture_output=True, text=True)
    got = OUT / "ABSTRACT_1.pdf"
    if r.returncode != 0 or not got.exists():
        raise SystemExit(f"abstract page did not extract:\n{r.stderr}")
    got.rename(OUT / "ABSTRACT.pdf")

    # the two copies of each figure, and the record that they are one file
    pairs = []
    for stem in MAIN_FIGS + SUPP_FIGS:
        label = "Figure" + stem.replace("fig", "")
        upload = OUT / "figures" / f"{label}.pdf"
        compile_ = OUT / "latex_source" / "figures" / f"{stem}.pdf"
        _copy(f / f"{stem}.pdf", upload)
        _copy(f / f"{stem}.pdf", compile_)
        _copy(f / f"{stem}.png", OUT / "figures" / f"{label}.png")
        pairs.append((stem, upload, compile_, f / f"{stem}.pdf"))

    # --- the assertions -------------------------------------------------
    bad = []
    for dst, src in shipped.items():
        if _sha(dst) != _sha(src):
            bad.append(f"{dst.name} differs from its canonical source {src}")
    for stem, upload, compile_, canon in pairs:
        h = {_sha(upload), _sha(compile_), _sha(canon)}
        if len(h) != 1:
            bad.append(f"{stem}: the upload copy, the compile copy and "
                       f"figures/{stem}.pdf are not the same file")
    if bad:
        raise SystemExit("submission package is inconsistent:\n  "
                         + "\n  ".join(bad))

    return {"files": len(shipped) + 2 + 3 * len(pairs), "figure_pairs": len(pairs)}


def compileInIsolation() -> tuple[int, int]:
    """build latex_source in a scratch tree with no access to the repository."""
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp) / "latex_source"
        shutil.copytree(OUT / "latex_source", work)
        here = work / "manuscript"
        pages = []
        for stem in ("bmb_v5", "bmb_v5_supplementary"):
            for _ in range(2):
                r = subprocess.run(
                    ["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                     stem],
                    cwd=here, capture_output=True, text=True)
            if r.returncode != 0:
                tail = "\n".join(r.stdout.splitlines()[-25:])
                raise SystemExit(f"latex_source does not compile ({stem}):\n{tail}")
            out = subprocess.run(["pdfinfo", str(here / f"{stem}.pdf")],
                                 capture_output=True, text=True).stdout
            pages.append(int([l for l in out.splitlines()
                              if l.startswith("Pages")][0].split()[-1]))
        return pages[0], pages[1]


def main() -> int:
    info = build()
    print(f"FINAL_SUBMISSION assembled: {info['files']} files, "
          f"{info['figure_pairs']} figure pairs checksum-matched")
    if "--compile" in sys.argv:
        a, b = compileInIsolation()
        print(f"isolated compile: {a} pages main, {b} pages supplementary")
        for pdf, n in ((OUT / "MANUSCRIPT.pdf", a),
                       (OUT / "SUPPLEMENTARY.pdf", b)):
            out = subprocess.run(["pdfinfo", str(pdf)], capture_output=True,
                                 text=True).stdout
            got = int([l for l in out.splitlines()
                       if l.startswith("Pages")][0].split()[-1])
            if got != n:
                raise SystemExit(
                    f"{pdf.name} ships {got} pages but latex_source builds {n}")
        print("shipped PDFs and the isolated build agree on page count")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
