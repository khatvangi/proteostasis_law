"""Convert `manuscript/bmb_v4.md` to submission LaTeX, and build the PDF.

WHY A CONVERTER AND NOT A HAND-PORT
The markdown is the single source. A hand-made .tex would be a second lineage for
every number in the paper, which is the failure this project keeps catching. This
script is the only thing that writes `manuscript/bmb_v4.tex`; that file is a build
artefact and must never be edited.

THE HARD PART IS THE INLINE CODE SPANS
276 backtick spans, of which exactly TWO are literal code and the rest are
mathematics written in backticks because markdown has no math. Guessing which is
which with a regex would silently mangle the paper, so the rule is inverted: every
span is MATH unless it appears in `CODE_SPANS`, and any span that looks
code-shaped but is not listed raises. A new code span therefore fails the build
instead of quietly rendering as math.

UNICODE
pdflatex is the only engine available here, so no character may reach it
unmapped. Prose unicode is mapped by `\newunicodechar` declarations GENERATED FROM
THE DOCUMENT'S OWN INVENTORY, so a character that is added later and has no
mapping stops the build rather than vanishing.

CLASS
`article` with geometry, because no Springer class is installed on this machine.
The preamble is written so that swapping in `sn-jnl.cls` or `svjour3.cls` is a
one-line change; see the banner at the top of the generated file.
"""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
SRC = REPO_ROOT / "manuscript" / "bmb_v4.md"
OUT_TEX = REPO_ROOT / "manuscript" / "bmb_v4.tex"
OUT_PDF = REPO_ROOT / "manuscript" / "bmb_v4.pdf"

# the only backtick spans that are literal code. everything else is mathematics.
CODE_SPANS = frozenset({
    "scripts/figures/",
    "scripts/figures/fig_identity.py:captionNumbers",
})

# a span shaped like code must be declared above or the build stops
_CODE_SHAPED = re.compile(r"\.py\b|\.tsv\b|\.md\b|scripts/|tests/|data/|[a-z_]+\.[A-Z]")

_SUP = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5", "⁶": "6",
        "⁷": "7", "⁸": "8", "⁹": "9", "⁻": "-", "⁺": "+"}
_SUB = {"₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4", "₅": "5", "₆": "6",
        "₇": "7", "₈": "8", "₉": "9"}

_GREEK = {"α": r"\alpha", "β": r"\beta", "γ": r"\gamma", "δ": r"\delta",
          "ε": r"\varepsilon", "θ": r"\theta", "ι": r"\iota", "λ": r"\lambda",
          "μ": r"\mu", "ν": r"\nu", "ρ": r"\rho", "σ": r"\sigma", "τ": r"\tau",
          "φ": r"\varphi", "χ": r"\chi", "Δ": r"\Delta"}

# operators, inside math mode
_MATH_OPS = {"∇": r"\nabla ", "×": r"\times ", "·": r"\cdot ", "−": "-",
             "≤": r"\le ", "≥": r"\ge ", "→": r"\to ", "∥": r"\parallel ",
             "±": r"\pm ", "∂": r"\partial ", "…": r"\dots ", "∎": r"\square",
             "√": r"\surd ", "∈": r"\in ", "≈": r"\approx ", "–": "-", "—": "---"}

# prose (text mode) mappings, emitted as \newunicodechar declarations
_TEXT_MAP = {
    "×": r"$\times$", "−": r"$-$", "–": "--", "—": "---", "∇": r"$\nabla$",
    "·": r"$\cdot$", "≤": r"$\le$", "≥": r"$\ge$", "→": r"$\to$",
    "∥": r"$\parallel$", "±": r"$\pm$", "∂": r"$\partial$", "…": r"\dots",
    "∎": r"$\square$", "√": r"$\surd$", "§": r"\S", "°": r"$^\circ$",
    "∈": r"$\in$", "≈": r"$\approx$", "″": "''", "′": "'",
    "ö": r'\"{o}', "é": r"\'{e}", "ü": r'\"{u}', "ä": r'\"{a}', "å": r"\aa{}",
    "’": "'", "‘": "`", "“": "``", "”": "''",
}
_TEXT_MAP.update({k: f"${v}$" for k, v in _GREEK.items()})
_TEXT_MAP.update({k: f"$^{{{v}}}$" for k, v in _SUP.items()})
_TEXT_MAP.update({k: f"$_{{{v}}}$" for k, v in _SUB.items()})

# multi-letter operators that must not be italicised as a product of letters.
# LaTeX already defines the first group; the second needs \operatorname.
_BUILTIN_OPS = ("det", "sin", "cos", "max", "min", "ln", "log", "exp")
_OTHER_OPS = ("row", "ceiling", "cross")


def _runsToScript(s: str, table: dict, wrapper: str) -> str:
    """collapse a run of super/subscript characters into one ^{...} or _{...}."""
    out, i = [], 0
    while i < len(s):
        if s[i] in table:
            j = i
            while j < len(s) and s[j] in table:
                j += 1
            out.append(wrapper % "".join(table[c] for c in s[i:j]))
            i = j
        else:
            out.append(s[i])
            i += 1
    return "".join(out)


def mathify(span: str) -> str:
    """render one backtick span as LaTeX inline math.

    Order matters and each step exists because a real span broke without it:
      * `f_{iota,I}` already carries a LaTeX subscript group, so existing
        `_{...}`/`^{...}` groups are protected BEFORE set braces like `{G = 0}`
        are escaped, or the escaping destroys them.
      * `j_crit` in math means `j_c` followed by `rit`. Multi-character
        subscripts MUST be braced; this is silent corruption otherwise, and it
        appears dozens of times in this paper.
      * pandoc will not parse `$...$` as math when a space precedes the closing
        delimiter, and escapes both dollars instead, so the result is stripped.
    """
    s = span

    # 1. protect sub/superscript groups the source already wrote
    protected: list[str] = []

    def _protect(m):
        protected.append(m.group(0))
        return f"\x00{len(protected) - 1}\x00"

    s = re.sub(r"[_^]\{[^{}]*\}", _protect, s)
    # 2. remaining braces are set delimiters
    s = s.replace("{", r"\{").replace("}", r"\}")
    s = re.sub(r"\x00(\d+)\x00", lambda m: protected[int(m.group(1))], s)

    # 3. unicode scripts -> real scripts
    s = _runsToScript(s, _SUP, "^{%s}")
    s = _runsToScript(s, _SUB, "_{%s}")

    # 4. letters and operators
    for k, v in _GREEK.items():
        s = s.replace(k, v + " ")
    for k, v in _MATH_OPS.items():
        s = s.replace(k, v)
    s = s.replace("~", r"\sim ")

    # 5. starred equilibria, then multi-character scripts
    s = re.sub(r"(?<=[A-Za-z\)\]])\*", "^{*}", s)
    s = re.sub(r"_(?!\{)([A-Za-z][A-Za-z0-9]*)", r"_{\1}", s)
    s = re.sub(r"\^\(([^()]*)\)", r"^{(\1)}", s)

    # 6. named operators, upright
    for name in _BUILTIN_OPS:
        s = re.sub(rf"(?<![A-Za-z\\]){name}(?![A-Za-z}}])", rf"\\{name} ", s)
    for name in _OTHER_OPS:
        s = re.sub(rf"(?<![A-Za-z\\]){name}(?![A-Za-z}}])",
                   rf"\\operatorname{{{name}}}", s)

    s = s.replace("%", r"\%").replace("&", r"\&").replace("#", r"\#")
    s = re.sub(r"\s+", " ", s).strip()
    s = re.sub(r"\s+(?=[_^])", "", s)
    return f"${s}$"


def texttt(span: str) -> str:
    s = span.replace("_", r"\_").replace("/", "/\\allowbreak ")
    return rf"\texttt{{{s}}}"


def convertSpans(md: str) -> tuple[str, int, int]:
    """replace every backtick span; math by default, code only when declared."""
    n_math = n_code = 0
    unknown = []

    def repl(m):
        nonlocal n_math, n_code
        span = m.group(1)
        if span in CODE_SPANS:
            n_code += 1
            return texttt(span)
        if _CODE_SHAPED.search(span):
            unknown.append(span)
            return span
        n_math += 1
        return mathify(span)

    out = re.sub(r"`([^`\n]+)`", repl, md)
    if unknown:
        raise SystemExit(
            "code-shaped spans not declared in CODE_SPANS:\n  "
            + "\n  ".join(sorted(set(unknown)))
            + "\nadd them to CODE_SPANS or rewrite them; the build will not guess.")
    return out, n_math, n_code


def convertFigures(md: str) -> tuple[str, int]:
    """turn each image + bold caption pair into a real float.

    the manuscript's width discipline is 84 mm or 174 mm, so the float asks for
    the figure's own width rather than scaling it to the text block.
    """
    n = 0
    pattern = re.compile(
        r"!\[Figure ([0-9S]+)\]\(\.\./figures/(fig[0-9S]+)\.pdf\)\n\n"
        r"\*\*Fig\. \1\*\* (.*?)(?=\n\n)", re.S)

    def repl(m):
        nonlocal n
        n += 1
        num, stem, caption = m.group(1), m.group(2), m.group(3)
        caption = " ".join(caption.split())
        return (
            "\\begin{figure}[htbp]\n\\centering\n"
            f"\\includegraphics[width=\\linewidth,max width=\\linewidth]{{../figures/{stem}.pdf}}\n"
            f"\\caption{{{caption}}}\n\\label{{fig:{stem}}}\n\\end{{figure}}")

    out = pattern.sub(repl, md)
    return out, n


def residualUnicode(text: str) -> list[str]:
    return sorted({c for c in text if ord(c) > 127})


def preamble(chars: list[str]) -> str:
    lines = [rf"\newunicodechar{{{c}}}{{{_TEXT_MAP[c]}}}" for c in chars
             if c in _TEXT_MAP]
    missing = [c for c in chars if c not in _TEXT_MAP]
    if missing:
        raise SystemExit("no text-mode mapping for: "
                         + ", ".join(f"U+{ord(c):04X} {c!r}" for c in missing))
    return "\n".join(lines)


HEADER = r"""%% GENERATED FILE -- do not edit.
%% Written by scripts/manuscript/to_latex.py from manuscript/bmb_v4.md.
%% The markdown is the source; every number in it is recomputed by a script in
%% scripts/ and asserted by tests/.
%%
%% TO SUBMIT TO BULLETIN OF MATHEMATICAL BIOLOGY: replace the \documentclass
%% line below with Springer's, which is not installed on the build machine:
%%     \documentclass[sn-mathphys-num]{sn-jnl}
%% and drop the geometry line. Nothing else here is class-specific.
\documentclass[11pt,a4paper]{article}
\usepackage[margin=25mm]{geometry}

\usepackage[T1]{fontenc}
\usepackage[utf8]{inputenc}
\usepackage{lmodern}
\usepackage{microtype}
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage[export]{adjustbox}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}
\usepackage{textcomp}
\usepackage{newunicodechar}
\usepackage[hidelinks]{hyperref}
\usepackage[font=small,labelfont=bf]{caption}
\usepackage{parskip}

\setlength{\emergencystretch}{3em}
\providecommand{\tightlist}{%
  \setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}

%% unicode that survives into prose, mapped from this document's own inventory
"""


def build() -> dict:
    md = SRC.read_text()
    md, n_fig = convertFigures(md)
    md, n_math, n_code = convertSpans(md)

    body = subprocess.run(
        ["pandoc", "--from=markdown+raw_tex+pipe_tables", "--to=latex",
         "--wrap=preserve", "--top-level-division=section"],
        input=md, capture_output=True, text=True, check=True).stdout

    # pandoc escapes nothing inside our raw figure blocks, but it does turn the
    # setext title into a section; the title/byline block is rebuilt explicitly
    body = re.sub(r"^\\section\{An Exact Collapse Threshold.*?\}\n", "", body,
                  count=1, flags=re.S)

    chars = residualUnicode(body)
    tex = (HEADER + preamble(chars) + "\n\n"
           + r"\title{An Exact Collapse Threshold for Conserved-Resource Models"
             "\n" + r"       of Protein Quality Control}" + "\n"
           + r"\author{Kiran Boggavarapu\\[2pt]"
             "\n" + r"{\small Department of Chemistry and Physics, McNeese State "
             r"University, Lake Charles, LA 70609, USA}\\" + "\n"
           + r"{\small \texttt{kiran@mcneese.edu}}}" + "\n"
           + r"\date{}" + "\n\n"
           + r"\begin{document}" + "\n" + r"\maketitle" + "\n\n"
           + body + "\n" + r"\end{document}" + "\n")
    OUT_TEX.write_text(tex)
    return {"spans_math": n_math, "spans_code": n_code, "figures": n_fig,
            "unicode_mapped": len(chars), "tex_bytes": len(tex)}


def compile() -> int:
    """two pdflatex passes, for hyperref and any cross-references."""
    for i in range(2):
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode",
                            "-halt-on-error", "-output-directory",
                            str(OUT_TEX.parent), str(OUT_TEX)],
                           capture_output=True, text=True)
        if r.returncode != 0:
            tail = [l for l in r.stdout.splitlines() if l.startswith("!")
                    or "l." == l[:2]]
            print("\n".join(r.stdout.splitlines()[-40:]))
            print("\n--- errors ---\n" + "\n".join(tail[:20]))
            return r.returncode
    for ext in (".aux", ".log", ".out", ".toc"):
        p = OUT_TEX.with_suffix(ext)
        if p.exists():
            p.unlink()
    return 0


if __name__ == "__main__":
    info = build()
    print("bmb_v4.tex written")
    for k, v in info.items():
        print(f"  {k:16s} {v}")
    rc = compile()
    if rc == 0:
        print(f"  pdf              {OUT_PDF.name} "
              f"({OUT_PDF.stat().st_size // 1024} kB)")
    sys.exit(rc)
