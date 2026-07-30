#!/usr/bin/env python3
"""
assemble the complete paper: prose + inlined tables + embedded figures.

MANUSCRIPT.md carries the prose and a placeholder for each main table:

    <!-- TABLE:Table 3 -->

this script replaces each placeholder with the generated block from
tables/TABLES.md, so a table body is never typed into the manuscript by hand.
that is deliberate: the one thing this project has failed at repeatedly is a
number that was correct where it was computed and stale where it was quoted.

outputs, all in manuscript/:
    PAPER.md     assembled markdown, tables inline
    PAPER.html   self-contained (figures base64-embedded), readable in a browser
    PAPER.pdf    printed from the html by headless chrome, if available
    PAPER.docx   for submission systems that want a word file

fails closed: an unresolved placeholder, a missing table block, or a missing
figure is an error, not a warning. run scripts/08_make_tables.py and the figure
scripts first.
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
MSDIR = ROOT / "manuscript"
MS = MSDIR / "MANUSCRIPT.md"
TABLES = ROOT / "tables" / "TABLES.md"
FIGS = ROOT / "figures"

PLACEHOLDER = re.compile(r"^<!-- TABLE:(.+?) -->$", re.M)

# printed-page styling for the html and, through it, the pdf. single column,
# because the figures are already sized for a journal column and a two-column
# css layout would scale them down twice.
CSS = """
:root { --ink:#1a1a1a; --muted:#555; --rule:#d0d0d0; }
body { font-family: "Georgia","Times New Roman",serif; color:var(--ink);
       max-width: 44rem; margin: 0 auto; padding: 2.5rem 1.5rem 4rem;
       line-height: 1.55; font-size: 11.5pt; }
h1 { font-size: 1.7rem; line-height:1.25; margin-bottom:.3rem; }
h2 { font-size: 1.15rem; margin-top: 2.2rem; border-bottom:1px solid var(--rule);
     padding-bottom:.25rem; }
h3 { font-size: 1.02rem; margin-top: 1.8rem; font-style: italic; }
p { margin: .7rem 0; text-align: justify; hyphens: auto; }
blockquote { border-left:3px solid var(--rule); margin:1.2rem 0; padding:.3rem 1rem;
             color:var(--muted); font-size:.92em; }
code, pre { font-family:"DejaVu Sans Mono",monospace; font-size:.88em; }
pre { background:#f6f6f4; padding:.7rem .9rem; overflow-x:auto; border-radius:3px; }
img { max-width:100%; height:auto; display:block; margin:1.4rem auto .4rem; }
table { border-collapse: collapse; width:100%; font-size:.78em; margin:.8rem 0 1.4rem;
        font-family:"Helvetica Neue",Arial,sans-serif; display:block; overflow-x:auto; }
th, td { border-top:1px solid var(--rule); padding:.28rem .45rem; text-align:left;
         vertical-align:top; }
thead th { border-bottom:1.5px solid var(--ink); border-top:1.5px solid var(--ink);
           font-weight:600; }
em strong, strong em { font-variant: small-caps; }
@media print {
  body { max-width:none; font-size:10pt; padding:0; }
  h2, h3 { break-after: avoid; }
  img, table, pre { break-inside: avoid; }
}
"""


def table_blocks(text):
    """
    split TABLES.md into {"Table 3": block}. a block runs from its "## Table N."
    heading to the next level-2 heading, so "### Table 2b." stays inside Table 2.
    """
    out, label, buf = {}, None, []
    for line in text.splitlines():
        m = re.match(r"^## (Table S?\d+[a-z]?)\.", line)
        if line.startswith("## ") or line.startswith("# "):
            if label:
                out[label] = "\n".join(buf).strip()
            label, buf = (m.group(1) if m else None), []
        if label:
            buf.append(line)
    if label:
        out[label] = "\n".join(buf).strip()
    return out


def as_paper_block(block):
    """demote the table's own headings to bold captions inside the flowing text."""
    lines = []
    for line in block.splitlines():
        m = re.match(r"^#{2,3} (.+?)\s*$", line)
        lines.append(f"**{m.group(1)}**" if m else line)
    return "\n".join(lines)


def assemble():
    ms = MS.read_text()
    blocks = table_blocks(TABLES.read_text())

    missing = []

    def sub(m):
        label = m.group(1).strip()
        if label not in blocks:
            missing.append(label)
            return m.group(0)
        return as_paper_block(blocks[label])

    paper = PLACEHOLDER.sub(sub, ms)
    if missing:
        sys.exit(f"ERROR: no block in tables/TABLES.md for: {missing}. "
                 "run scripts/08_make_tables.py first.")
    if "<!-- TABLE:" in paper:
        sys.exit("ERROR: an unresolved table placeholder survived assembly")

    # every figure the manuscript embeds must exist
    for rel in re.findall(r"!\[[^\]]*\]\(([^)]+)\)", paper):
        p = (MSDIR / rel).resolve()
        if not p.exists():
            sys.exit(f"ERROR: manuscript embeds a missing figure: {rel}")

    # the status note is provenance for the repository, not part of the paper
    paper = re.sub(r"\n> \*\*Status\.\*\*.*?\n\n---\n", "\n", paper, flags=re.S)
    return paper


def run(cmd, what):
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return True
    except FileNotFoundError:
        print(f"  skipped {what}: {cmd[0]} not installed")
    except subprocess.CalledProcessError as e:
        print(f"  skipped {what}: {cmd[0]} failed\n{e.stderr[-600:]}")
    return False


def main():
    paper = assemble()
    md = MSDIR / "PAPER.md"
    md.write_text(paper)
    n_tab = paper.count("\n|")
    print(f"wrote {md.relative_to(ROOT)}  "
          f"({len(paper.splitlines())} lines, {n_tab} table rows)")

    css = MSDIR / ".paper.css"
    css.write_text(CSS)
    html = MSDIR / "PAPER.html"
    ok_html = run(["pandoc", str(md), "-o", str(html), "--standalone",
                   "--embed-resources", "--css", str(css),
                   "--metadata", "title=A finite proteostasis envelope for translation",
                   f"--resource-path={MSDIR}"], "PAPER.html")
    if ok_html:
        print(f"wrote {html.relative_to(ROOT)}  "
              f"({html.stat().st_size / 1e6:.1f} MB, figures embedded)")

    run(["pandoc", str(md), "-o", str(MSDIR / "PAPER.docx"),
         f"--resource-path={MSDIR}"], "PAPER.docx")
    if (MSDIR / "PAPER.docx").exists():
        print(f"wrote manuscript/PAPER.docx")

    # pdflatex cannot set this text (superscripts, greek, µ, minus sign) and no
    # unicode-capable tex engine is installed, so print the html instead
    chrome = shutil.which("google-chrome") or shutil.which("chromium")
    if ok_html and chrome:
        pdf = MSDIR / "PAPER.pdf"
        if run([chrome, "--headless", "--disable-gpu", "--no-sandbox",
                "--no-pdf-header-footer", f"--print-to-pdf={pdf}",
                html.as_uri()], "PAPER.pdf") and pdf.exists():
            print(f"wrote {pdf.relative_to(ROOT)}  "
                  f"({pdf.stat().st_size / 1e6:.1f} MB)")
    css.unlink(missing_ok=True)


if __name__ == "__main__":
    main()
