"""Convert `manuscript/MANUSCRIPT_BMB_v5.md` to submission LaTeX, and build the PDF.

WHY A CONVERTER AND NOT A HAND-PORT
The markdown is the single source. A hand-made .tex would be a second lineage for
every number in the paper, which is the failure this project keeps catching. This
script is the only thing that writes `manuscript/bmb_v5.tex`; that file is a build
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
SRC = REPO_ROOT / "manuscript" / "MANUSCRIPT_BMB_v5.md"
OUT_TEX = REPO_ROOT / "manuscript" / "bmb_v5.tex"
OUT_PDF = REPO_ROOT / "manuscript" / "bmb_v5.pdf"
OUT_SUPP_TEX = REPO_ROOT / "manuscript" / "bmb_v5_supplementary.tex"

# the only backtick spans that are literal code. everything else is mathematics.
CODE_SPANS = frozenset({
    "scripts/figures/",
    "scripts/figures/fig_identity.py:captionNumbers",
})

# a span shaped like code must be declared above or the build stops
_CODE_SHAPED = re.compile(r"\.py\b|\.tsv\b|\.md\b|scripts/|tests/|data/|[a-z_]+\.[A-Z]")

_SUP = {"⁰": "0", "¹": "1", "²": "2", "³": "3", "⁴": "4", "⁵": "5", "⁶": "6",
        "⁷": "7", "⁸": "8", "⁹": "9", "⁻": "-", "⁺": "+", "ⁿ": "n"}
_SUB = {"₀": "0", "₁": "1", "₂": "2", "₃": "3", "₄": "4", "₅": "5", "₆": "6",
        "₇": "7", "₈": "8", "₉": "9"}

_GREEK = {"α": r"\alpha", "β": r"\beta", "γ": r"\gamma", "δ": r"\delta",
          "ε": r"\varepsilon", "θ": r"\theta", "ι": r"\iota", "λ": r"\lambda",
          "μ": r"\mu", "ν": r"\nu", "ρ": r"\rho", "σ": r"\sigma", "τ": r"\tau",
          "φ": r"\varphi", "χ": r"\chi", "Δ": r"\Delta",
          "κ": r"\kappa", "π": r"\pi", "ω": r"\omega"}

# operators, inside math mode
_MATH_OPS = {"∇": r"\nabla ", "×": r"\times ", "·": r"\cdot ", "−": "-",
             "≤": r"\le ", "≥": r"\ge ", "→": r"\to ", "∥": r"\parallel ",
             "±": r"\pm ", "∂": r"\partial ", "…": r"\dots ", "∎": r"\square",
             "∈": r"\in ", "≈": r"\approx ", "–": "-", "—": "---",
             "≠": r"\neq ",
             # U+2016 NORM. symmetric, so one mapping serves both delimiters;
             # \lVert would be wrong on the closing one. NOT U+2225 above:
             # `∇R ∥ ∇G` is a relation, `‖∇G‖` is a magnitude.
             "‖": r"\|"}

# prose (text mode) mappings, emitted as \newunicodechar declarations
_TEXT_MAP = {
    "×": r"$\times$", "−": r"$-$", "–": "--", "—": "---", "∇": r"$\nabla$",
    "·": r"$\cdot$", "≤": r"$\le$", "≥": r"$\ge$", "→": r"$\to$",
    "∥": r"$\parallel$", "±": r"$\pm$", "∂": r"$\partial$", "…": r"\dots",
    "‖": r"$\|$",
    "∎": r"$\square$", "√": r"$\surd$", "§": r"\S", "°": r"$^\circ$",
    "∈": r"$\in$", "≈": r"$\approx$", "≠": r"$\neq$", "″": "''", "′": "'",
    "ö": r'\"{o}', "é": r"\'{e}", "ü": r'\"{u}', "ä": r'\"{a}', "å": r"\aa{}",
    "’": "'", "‘": "`", "“": "``", "”": "''",
}
_TEXT_MAP.update({k: f"${v}$" for k, v in _GREEK.items()})
_TEXT_MAP.update({k: f"$^{{{v}}}$" for k, v in _SUP.items()})
_TEXT_MAP.update({k: f"$_{{{v}}}$" for k, v in _SUB.items()})

# multi-letter operators that must not be italicised as a product of letters.
# LaTeX already defines the first group; the second needs \operatorname.
_BUILTIN_OPS = ("det", "sin", "cos", "max", "min", "ln", "log", "exp", "sup")
# `tr J` set t, r and J as three italic variables. LaTeX has no \tr, so it goes
# here and not above -- putting it in the builtin list emitted \tr and the build
# stopped, which is the behaviour that list is supposed to have.
# Re and Im are the same case as tr: two italic variables otherwise.
_OTHER_OPS = ("row", "ceiling", "cross", "tr", "Re", "Im")


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


def _singleBars(s: str) -> str:
    r"""``\\|`` is markdown's escape for a literal pipe inside a table cell.

    Passing it through to LaTeX math sets a DOUBLE bar, so |tr J|, |d2R/ds2| and
    |w_1| -- all scalars -- rendered as norms. The author wrote absolute value.
    """
    return s.replace(r"\|", "|")


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
    s = _singleBars(span)

    # 1. protect sub/superscript groups the source already wrote
    protected: list[str] = []

    def _protect(m):
        protected.append(m.group(0))
        return f"\x00{len(protected) - 1}\x00"

    s = re.sub(r"[_^]\{[^{}]*\}", _protect, s)
    # 2. remaining braces are set delimiters
    s = s.replace("{", r"\{").replace("}", r"\}")
    s = re.sub(r"\x00(\d+)\x00", lambda m: protected[int(m.group(1))], s)

    # `\sqrt` needs its argument braced. This MUST run after the brace escaping
    # above: done before it, the braces this introduces were themselves escaped
    # into set delimiters, so `±i√det J` reached the PDF as `\sqrt\{\det J\}`
    # -- an empty radical followed by a literal set. Rendered, caught, moved.
    s = re.sub(r"√\s*\\?([A-Za-z]+)\s+([A-Za-z])", r"\\sqrt{\\\1 \2}", s)
    s = re.sub(r"√\s*\(([^()]*)\)", r"\\sqrt{\1}", s)
    s = s.replace("√", r"\surd ")

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
        # The caption bypasses pandoc, so its markdown emphasis is converted here
        # -- the first build printed a literal "**(a)**" in Figure 1. Backtick
        # spans are PROTECTED first: `(u*, a*)` is mathematics, and the emphasis
        # rule matched from the first star to the second and turned it into
        # \emph{, a}, silently deleting the starred equilibria. That survived a
        # full visual read of the article build.
        held: list[str] = []

        def _hold(mm):
            held.append(mm.group(0))
            return f"\x00{len(held) - 1}\x00"

        caption = re.sub(r"`[^`\n]+`", _hold, caption)
        caption = re.sub(r"\*\*(.+?)\*\*", r"\\textbf{\1}", caption)
        caption = re.sub(r"(?<!\*)\*([^*]+?)\*(?!\*)", r"\\emph{\1}", caption)
        caption = re.sub(r"\x00(\d+)\x00", lambda mm: held[int(mm.group(1))],
                         caption)
        # A literal % opens a LaTeX comment. Figure 5's caption quotes several
        # percentages, so the unescaped % swallowed the closing brace of
        # \caption{...}; pandoc could not find \end{figure}, escaped the WHOLE
        # raw block into \textbackslash begin\{figure\}, and the orphaned
        # \centering then centred every paragraph to the end of the document.
        # The build reported "figures 6" throughout. Hence the assertion below.
        caption = re.sub(r"(?<!\\)%", r"\\%", caption)
        caption = re.sub(r"(?<!\\)&", r"\\&", caption)
        return (
            "\\begin{figure}[htbp]\n\\centering\n"
            # NATURAL SIZE, clamped. The figures are built at exactly 84 mm or 174 mm
            # because that is what the journal specifies, and every font size in
            # them was chosen for that width. `width=\linewidth` rescaled an 84 mm
            # panel to the 160 mm text block, enlarging its 7 pt labels by 1.9x and
            # throwing away the entire width discipline. `max width` shrinks only
            # what would otherwise overflow.
            f"\\includegraphics[max width=\\linewidth]{{../figures/{stem}.pdf}}\n"
            f"\\caption{{{caption}}}\n\\label{{fig:{stem}}}\n\\end{{figure}}")

    out = pattern.sub(repl, md)
    return out, n


def naturalTableWidths(tex: str) -> tuple[str, int, int]:
    """let LaTeX size table columns instead of pandoc's equal fractions.

    Pandoc derives column widths from the markdown separator row, which is
    uniform here, so a 5-column table gets five 20% columns and the first one
    wraps "direct solver against continuation sweep, relative" over four lines
    while the numeric columns sit half empty. These tables are narrow, so `l`
    columns are better; anything that then overflows shows up as an overfull hbox
    rather than as silently bad typesetting.
    """
    n = 0

    def repl(m):
        nonlocal n
        n += 1
        cols = m.group(1).count(r"\raggedright")
        return r"\begin{longtable}[]{@{}" + "l" * cols + "@{}}"

    # the closing `@{}}` sits at the end of the LAST column line, not on its own,
    # so the span is matched non-greedily rather than line by line
    out = re.sub(r"\\begin\{longtable\}\[\]\{@\{\}\n(.*?)@\{\}\}",
                 repl, tex, flags=re.S)
    # `n` counts CONVERSIONS, not tables: pandoc emits simple `@{}ll@{}` columns
    # directly for narrow tables, and those never match the pattern above. In v4
    # the paper had exactly two tables and both needed conversion, so the two
    # numbers coincided and the mislabel was invisible; v5 has four tables of
    # which two needed conversion, and the assertion named "tables" passed at 2
    # while the document contained 4. Count the environments themselves.
    n_tables = out.count(r"\begin{longtable}")
    # pandoc wraps header cells in minipages sized to \linewidth; with `l`
    # columns that would stretch each header to the full text width
    out = re.sub(r"\\begin\{minipage\}\[b\]\{\\linewidth\}\\raggedright\n(.*?)\n\\end\{minipage\}",
                 lambda m: m.group(1).strip(), out, flags=re.S)
    return out, n_tables, n


def residualUnicode(text: str) -> list[str]:
    return sorted({c for c in text if ord(c) > 127})


def longtableToTabular(tex: str) -> tuple[str, int]:
    """turn pandoc's longtables into ordinary `table` floats.

    `longtable` runs its own output routine. On a page that ALSO carries a
    top-placed float, that interaction displaced the page folio into the middle
    of the body text -- page 6 printed its "6" over the "u" of "input-output",
    and the page carried no folio at its foot. Every count assertion passed and
    the LaTeX log reported no overfull box; only rendering the page showed it.

    These tables are short and none of them spans a page, so a float is both the
    correct environment and the one Springer expects. The ``\\LTcaptype`` wrapper
    that suppressed longtable's counter goes with it.
    """
    n = 0

    def repl(m):
        nonlocal n
        n += 1
        cols, body = m.group(1), m.group(2)
        body = re.sub(r"\\(endhead|endlastfoot|endfirsthead)\b", "", body)
        body = body.replace(r"\noalign{}", "")
        body = re.sub(r"\\(top|mid|bottom)rule", "", body)
        rows = [ln.rstrip() for ln in body.strip().split("\n") if ln.strip()]
        head, rest = rows[0], rows[1:]
        return ("\\begin{table}[htbp]\n\\centering\n"
                + "\\begin{tabular}{" + cols + "}\n\\toprule\n"
                + head + "\n\\midrule\n" + "\n".join(rest)
                + "\n\\bottomrule\n\\end{tabular}\n\\end{table}")

    # the column spec is `@{}lll@{}`, which CONTAINS braces, so it cannot be
    # matched with `[^}]*` -- that stops at the first `}` and the whole pattern
    # fails silently, converting nothing. Take everything to the last `}` on the
    # line instead.
    out = re.sub(
        r"\{\\def\\LTcaptype\{none\}[^\n]*\n"
        r"\\begin\{longtable\}\[\]\{([^\n]+)\}\n(.*?)\\end\{longtable\}\n\}",
        repl, tex, flags=re.S)
    return out, n


def preamble(chars: list[str]) -> str:
    lines = [rf"\newunicodechar{{{c}}}{{{_TEXT_MAP[c]}}}" for c in chars
             if c in _TEXT_MAP]
    missing = [c for c in chars if c not in _TEXT_MAP]
    if missing:
        raise SystemExit("no text-mode mapping for: "
                         + ", ".join(f"U+{ord(c):04X} {c!r}" for c in missing))
    return "\n".join(lines)


HEADER = r"""%% GENERATED FILE -- do not edit.
%% Written by scripts/manuscript/to_latex.py from manuscript/MANUSCRIPT_BMB_v5.md.
%% The markdown is the source; every number in it is recomputed by a script in
%% scripts/ and asserted by tests/.
%%
%% CLASS: Springer Nature sn-jnl, vendored at manuscript/springer/sn-jnl.cls
%% (LPPL). The option is sn-mathphys-ay -- AUTHOR-YEAR, which is what Bulletin of
%% Mathematical Biology uses and how the reference list below is written. The
%% Springer template ships with sn-mathphys-num uncommented; selecting it would
%% silently reformat the reference list into numbered style and break every
%% in-text citation, which are author-year throughout.
\documentclass[pdflatex,sn-mathphys-ay]{sn-jnl}

%% the class loads fontenc, geometry, hyperref, natbib, amsthm and others itself
\usepackage{amsmath,amssymb}
\usepackage{graphicx}
\usepackage[export]{adjustbox}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{array}
\usepackage{textcomp}

\usepackage{newunicodechar}
%% pandoc marks an uncaptioned longtable with \LTcaptype{none}; the caption
%% package then tries to step a counter of that name. Declaring it makes the
%% marker the no-op pandoc intends rather than a fatal error.
\newcounter{none}

\setlength{\emergencystretch}{3em}
\providecommand{\tightlist}{%
  \setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}

%% unicode that survives into prose, mapped from this document's own inventory
"""


# ---------------------------------------------------------------------------
# display equations
# ---------------------------------------------------------------------------
#
# The markdown writes displayed equations as fenced code blocks, because markdown
# has no display math. Pandoc renders those as VERBATIM, so every equation in the
# paper was typeset in typewriter face with ASCII operators, and the widest one
# ran past the right margin.
#
# Layout of a displayed equation is a judgement -- where to break, what to align
# on, when a matrix is a matrix -- so these are written out rather than derived.
# The dictionary is keyed by the EXACT source block, so editing an equation in the
# markdown invalidates its key and stops the build until the layout is revisited.
# That is the intended coupling: one lineage for the numbers, an explicit and
# self-invalidating mapping for the typesetting.

DISPLAY_MATH = {
"""du/dt = j − v_ref − v_degU − n − g + v_dis
da/dt =              n + g − v_dis − v_degA""":
r"""\begin{align*}
\frac{du}{dt} &= j - v_{ref} - v_{degU} - n - g + v_{dis},\\
\frac{da}{dt} &= n + g - v_{dis} - v_{degA}.
\end{align*}""",

# the closure, printed in free concentrations. the four equations are solved
# JOINTLY, so they are set as one aligned block rather than four displays:
# splitting them once invited the reading that each is evaluated in turn, which
# is the substitution D004 forbids.
"""u_f = u /(1 + C_f/K_CU + D_f/K_DU)     C_f = C_tot /(1 + ν + u_f/K_CU + a_f/K_CA)
a_f = a /(1 + C_f/K_CA + D_f/K_DA)     D_f = D_tot /(1 + u_f/K_DU + a_f/K_DA)""":
r"""\begin{align*}
u_f &= \frac{u}{1 + C_f/K_{CU} + D_f/K_{DU}},
&\qquad
C_f &= \frac{C_{tot}}{1 + \nu + u_f/K_{CU} + a_f/K_{CA}},\\[2pt]
a_f &= \frac{a}{1 + C_f/K_{CA} + D_f/K_{DA}},
&\qquad
D_f &= \frac{D_{tot}}{1 + u_f/K_{DU} + a_f/K_{DA}}.
\end{align*}""",

"""v_ref  = k_ref C_f u_f/(K_ref + u_f)      refolding
v_degU = k_U D_f u_f/(K_U + u_f)          soluble degradation
n      = k_n u_f^m,  m > 1                primary nucleation
g      = k_g u_f a_f                      elongation
v_dis  = k_dis C_f a_f/(K_dis + a_f)      disaggregation
v_degA = k_A D_f a_f/(K_A + a_f)          aggregate clearance""":
r"""\begin{align*}
v_{ref}  &= \frac{k_{ref}\,C_f\,u_f}{K_{ref} + u_f}   &&\text{refolding}\\
v_{degU} &= \frac{k_U\,D_f\,u_f}{K_U + u_f}           &&\text{soluble degradation}\\
n        &= k_n u_f^{\,m},\quad m > 1                 &&\text{nucleation}\\
g        &= k_g\,u_f\,a_f                             &&\text{aggregate growth}\\
v_{dis}  &= \frac{k_{dis}\,C_f\,a_f}{K_{dis} + a_f}   &&\text{disaggregation}\\
v_{degA} &= \frac{k_A\,D_f\,a_f}{K_A + a_f}           &&\text{aggregate clearance}
\end{align*}""",

"""R(u,a) = v_ref + v_degU + v_degA        G(u,a) = da/dt""":
r"""\begin{equation*}
R(u,a) = v_{ref} + v_{degU} + v_{degA},
\qquad
G(u,a) = \frac{da}{dt}
\end{equation*}""",

"""j_crit = R(u*, a*).""":
r"""\begin{equation*}
j_{crit} = R(u^{*}, a^{*}).
\end{equation*}""",

"""det J = −(∇R × ∇G)""":
r"""\begin{equation*}
\det J = -(\nabla R \times \nabla G)
\end{equation*}""",

# the proof display now ends the proof, so the tombstone lives inside it
"""det J = det [ ∂(du/dt)/∂u   ∂(du/dt)/∂a ]  =  det [ −R_u  −R_a ]
            [ ∂(da/dt)/∂u   ∂(da/dt)/∂a ]        [  G_u   G_a ]

      = −(R_u G_a − R_a G_u) = −(∇R × ∇G). ∎""":
r"""\begin{align*}
\det J
&= \det\begin{bmatrix}
     \partial(du/dt)/\partial u & \partial(du/dt)/\partial a\\[2pt]
     \partial(da/dt)/\partial u & \partial(da/dt)/\partial a
   \end{bmatrix}
 = \det\begin{bmatrix} -R_u & -R_a\\[2pt] G_u & G_a \end{bmatrix}\\[4pt]
&= -(R_u G_a - R_a G_u) = -(\nabla R \times \nabla G). \qquad\qquad\square
\end{align*}""",

"""r(s) = j_crit + ½ r''(0) s² + O(s³),""":
r"""\begin{equation*}
r(s) = j_{crit} + \tfrac{1}{2} r''(0)\, s^2 + O(s^3),
\end{equation*}""",

"""w·D²F(v,v) = −D²R(v,v) + λ D²G(v,v),""":
r"""\begin{equation*}
w \cdot D^2F(v,v) = -D^2R(v,v) + \lambda\, D^2G(v,v),
\end{equation*}""",

"""|w₁|/‖w‖ = (1 + (1+λ)²)^(−1/2),""":
r"""\begin{equation*}
\frac{|w_1|}{\lVert w \rVert} = \bigl(1 + (1+\lambda)^2\bigr)^{-1/2},
\end{equation*}""",

"""w·F_j = 1 − R_j + λ G_j,""":
r"""\begin{equation*}
w \cdot F_j = 1 - R_j + \lambda\, G_j,
\end{equation*}""",

"""r''(0) = −w·D²F(v,v),""":
r"""\begin{equation*}
r''(0) = -\,w \cdot D^2F(v,v),
\end{equation*}""",

"""r'(s) = ∇R·γ' = det J /‖∇G‖""":
r"""\begin{equation*}
r'(s) = \nabla R \cdot \gamma' = \frac{\det J}{\lVert \nabla G \rVert}
\end{equation*}""",

"""det J = −det [ ∇R ; ∇G ; ∇C ]""":
r"""\begin{equation*}
\det J = -\det\left[\,\nabla R \,;\, \nabla G \,;\, \nabla C\,\right]
\end{equation*}""",

"""R_tot = R + μ(u,a)·(u + a)""":
r"""\begin{equation*}
R_{tot} = R + \mu(u,a)\cdot(u + a)
\end{equation*}""",

"""j_crit = C_enz · φ_enz /(1 − δ),     φ_enz = R_enz(u*,a*)/C_enz,     δ = R_dil(u*,a*)/j_crit""":
r"""\begin{equation*}
j_{crit} = \frac{C_{enz}\,\varphi_{enz}}{1 - \delta},
\qquad
\varphi_{enz} = \frac{R_{enz}(u^{*},a^{*})}{C_{enz}},
\qquad
\delta = \frac{R_{dil}(u^{*},a^{*})}{j_{crit}}
\end{equation*}""",

"""j ≤ ( √(1 + 4εC_0) − 1 )/(2ε)   →   √(C_0/ε)   for large ε.""":
r"""\begin{equation*}
j \le \frac{\sqrt{1 + 4\varepsilon C_0} - 1}{2\varepsilon}
\;\longrightarrow\;
\sqrt{C_0/\varepsilon}
\quad\text{for large }\varepsilon.
\end{equation*}""",
}


def convertDisplayMath(md: str) -> tuple[str, int]:
    """replace every fenced block with its declared LaTeX, or stop."""
    n = 0
    missing = []

    def repl(m):
        nonlocal n
        src = m.group(1).rstrip("\n")
        if src not in DISPLAY_MATH:
            missing.append(src)
            return m.group(0)
        n += 1
        return DISPLAY_MATH[src] + "\n"

    out = re.sub(r"^```\n(.*?)^```\n", repl, md, flags=re.S | re.M)
    if missing:
        raise SystemExit(
            "display block with no declared layout:\n\n"
            + "\n\n".join(missing)
            + "\n\nAdd it to DISPLAY_MATH. Equation layout is not derived, and a "
              "fenced block left alone is typeset as verbatim.")
    return out, n


def convertHeadings(md: str) -> tuple[str, int, int]:
    """turn markdown headings into raw LaTeX sectioning commands.

    The markdown numbers its own headings ("## 1. Introduction") and the prose
    cross-references those numbers ("Section 4.2", "section 8.4"). LaTeX also
    numbers sections, so leaving the manual numbers in printed "0.1  1.
    Introduction". The manual number is therefore stripped and LaTeX's counter
    takes over -- but only after ASSERTING that the counter will produce the same
    number, since a mismatch would silently break every cross-reference in the
    paper.

    Unnumbered headings (Abstract, References, ...) become starred so they do not
    advance the counter.
    """
    lines, sec, sub = [], 0, 0
    for line in md.split("\n"):
        m3 = re.match(r"^### (\d+)\.(\d+)\s+(.*)$", line)
        m2 = re.match(r"^## (\d+)\.\s+(.*)$", line)
        mu3 = re.match(r"^### (?!\d)(.*)$", line)
        mu2 = re.match(r"^## (?!\d)(.*)$", line)
        if m2:
            sec += 1
            sub = 0
            if int(m2.group(1)) != sec:
                raise SystemExit(
                    f"heading numbering would shift: markdown says section "
                    f"{m2.group(1)}, LaTeX would emit {sec} -- every 'Section N' "
                    f"cross-reference in the prose depends on these agreeing")
            lines.append(rf"\section{{{m2.group(2).strip()}}}")
        elif m3:
            sub += 1
            if (int(m3.group(1)), int(m3.group(2))) != (sec, sub):
                raise SystemExit(
                    f"heading numbering would shift at {m3.group(0)!r}: LaTeX "
                    f"would emit {sec}.{sub}")
            lines.append(rf"\subsection{{{m3.group(3).strip()}}}")
        elif mu2:
            lines.append(rf"\section*{{{mu2.group(1).strip()}}}")
        elif mu3:
            lines.append(rf"\subsection*{{{mu3.group(1).strip()}}}")
        else:
            lines.append(line)
    return "\n".join(lines), sec, sub


def splitFrontMatter(md: str) -> tuple[str, str, str]:
    """separate title block, abstract and body.

    The markdown front matter carries the title, byline, MSC and keywords, and
    LaTeX's \\maketitle supplies the first two. Letting pandoc convert the block
    as prose printed the byline TWICE in the first build. The abstract likewise
    has to become an `abstract` environment, or it renders as a numbered
    subsection -- it appeared as "0.1 Abstract".
    """
    m = re.search(r"^## Abstract\n+(.*?)\n+---\n", md, re.S | re.M)
    if not m:
        raise SystemExit("front matter not in the expected shape")
    front = md[:m.start()]
    abstract = m.group(1).strip()
    body = md[m.end():]
    msc = re.search(r"\*\*MSC\*\* (.+)", front)
    kw = re.search(r"\*\*Keywords\*\* (.+)", front)
    meta = ""
    if kw:
        meta += r"\keywords{" + kw.group(1).strip().replace("; ", ", ") + "}\n\n"
    if msc:
        meta += r"\pacs[MSC Classification]{" + msc.group(1).strip() + "}\n"
    return abstract, meta, body


def _pandoc(md: str) -> str:
    return subprocess.run(
        ["pandoc", "--from=markdown+raw_tex+pipe_tables", "--to=latex",
         "--wrap=preserve", "--top-level-division=section"],
        input=md, capture_output=True, text=True, check=True).stdout


def splitSupplementary(md: str) -> tuple[str, str]:
    """lift the supplementary section out of the body into its own document.

    Journals take supplementary material as a separate file, and a supplementary
    figure that sits inside the main PDF gets numbered and typeset as if it were
    part of the paper. The section is removed from the body here and rendered as a
    standalone document with `\thefigure` redefined to S1, S2, ... so its numbers
    cannot collide with the main sequence.
    """
    # match the SECTION, not one spelling of its title: v4 called it
    # "Supplementary figure", v5 calls it "Supplementary material", and pinning
    # the exact token made a rename look like a missing section.
    m = re.search(r"^## Supplementary\b[^\n]*\n(.*?)(?=^## |\Z)", md,
                  re.S | re.M)
    if not m:
        raise SystemExit("no '## Supplementary ...' section found")
    supp = m.group(1).strip()
    body = (md[:m.start()] + md[m.end():])
    # the horizontal rule that separated it goes too
    body = re.sub(r"\n---\n\n---\n", "\n---\n", body)
    return supp, body


# sections kept in the markdown as project provenance but not part of the paper
INTERNAL_SECTIONS = ("Reference verification",)


def stripInternalSections(md: str) -> tuple[str, int]:
    """drop working sections from the submission without losing them upstream.

    `### Reference verification` records which references were checked against
    PubMed and what was corrected. That is provenance the repository should keep
    and a referee should not receive after the reference list.
    """
    n = 0
    for name in INTERNAL_SECTIONS:
        pat = re.compile(rf"^#{{2,3}} {re.escape(name)}\n.*?(?=^#{{1,3}} |\Z)",
                         re.S | re.M)
        md, k = pat.subn("", md)
        n += k
    return md, n


def _renderBody(md: str, *, headings: bool) -> tuple[str, dict]:
    """markdown -> LaTeX body, with the checks that caught real corruption."""
    n_sec = n_disp = 0
    md, n_strip = stripInternalSections(md)
    if headings:
        md, n_sec, _ = convertHeadings(md)
    md, n_disp = convertDisplayMath(md)
    md, n_fig = convertFigures(md)
    md, n_math, n_code = convertSpans(md)
    tex = _pandoc(md)
    if r"\textbackslash begin\{figure\}" in tex:
        raise SystemExit(
            "pandoc escaped a raw figure block instead of passing it through. "
            "Something in a caption breaks its LaTeX scanner -- an unescaped %, "
            "or an unbalanced brace. The build must not continue: the orphaned "
            "\\centering silently centres the rest of the document.")
    n_graphics = tex.count(r"\includegraphics")
    if n_graphics != n_fig:
        raise SystemExit(f"{n_fig} figures converted but {n_graphics} survived "
                         "pandoc")
    tex, n_tab, n_conv = naturalTableWidths(tex)
    tex, n_float = longtableToTabular(tex)
    if n_float != n_tab:
        raise SystemExit(f'{n_tab} tables but {n_float} converted to floats; a longtable left in place can displace the page folio')
    return tex, {"sections": n_sec, "figures": n_fig, "tables": n_tab,
                 "tables_rewidthed": n_conv,
                 "displays": n_disp, "stripped": n_strip,
                 "spans_math": n_math, "spans_code": n_code}


TITLE = ("An Exact Fold Condition for Mass-Balanced Models\n"
         "       of Protein Quality Control")
# the ORCID is in the markdown front matter and was being dropped, because this
# byline is hardcoded rather than parsed.
#
# NOT via sn-jnl's \orcid: that macro is `\gdef\orcid#1{\href{#1}{\orcidlogo}}`,
# which expects a URL and typesets `Orcidlogo.eps` -- a file the class ships
# separately and we have not vendored. Passing it a bare identifier killed the
# build. Springer collects the ORCID through the submission system anyway, so
# the identifier goes in the author note as text, where a reader can see it and
# no missing graphic can break the compile.
# The ORCID is NOT emitted into the PDF, and that is a decision rather than an
# omission. Four placements were tried and all failed:
#   \orcid{}        -- sn-jnl defines it as \href{#1}{\orcidlogo} and typesets
#                      Orcidlogo.eps, which the class ships separately and we
#                      have not vendored; the build dies.
#   redefined \orcid -- still dies: the macro cannot sit in that position in the
#                      class's author block whatever it expands to.
#   \footnotetext    -- silently swallowed before \maketitle. The build SUCCEEDED
#                      and the identifier never reached the PDF.
#   \presentaddress  -- renders, under the label "Present Address:", which is
#                      false. A wrong label is worse than an absent one.
# The identifier is carried by CITATION.cff, .zenodo.json and the markdown front
# matter, and Springer collects it through the submission system.
ORCID = "0000-0003-0751-6459"   # recorded here for the metadata files only
BYLINE = (r"\author*[1]{\fnm{Kiran} \sur{Boggavarapu}}"
          + r"\email{kiran@mcneese.edu}" + "\n"
          + r"\affil*[1]{\orgdiv{Department of Chemistry and Physics}, "
            r"\orgname{McNeese State University}, \orgaddress{\city{Lake Charles}, "
            r"\state{LA}, \postcode{70609}, \country{USA}}}" + "\n"
          )


def _document(chars: list[str], title: str, front: str, body: str,
              extra_preamble: str = "") -> str:
    """sn-jnl wants the abstract and keywords as PREAMBLE commands, not
    environments, and the title block uses \\fnm/\\sur/\\affil markup."""
    return (HEADER + preamble(chars) + "\n" + extra_preamble + "\n"
            + rf"\title[Fold condition for quality-control models]{{{title}}}"
            + "\n" + BYLINE + "\n\n" + front + "\n"
            + r"\begin{document}" + "\n" + r"\maketitle" + "\n\n"
            + body + "\n" + r"\end{document}" + "\n")


# What the finished paper contains. The build asserts these rather than reporting
# them, because "figures 6" appeared in the log through every one of the runs in
# which one figure was silently mangled -- the count was wrong, persistently, and
# nothing compared it against anything. A count assertion is far cheaper than
# reading 21 pages and would have caught most of what reading caught.
EXPECTED = {
    "sections": 10,          # numbered ## 1 .. ## 10
    "figures_main": 5,       # fig1 theorem, fig2 dilution, fig3 saturation,
                             # fig4 hopf, fig5 beta
    "figures_supp": 2,       # figS1 identity, figS2 pareto front
    "tables": 4,             # genericity, verification, section 6, beta
    # of the four tables, those pandoc emitted as longtables and this script
    # converted to `table` floats. Table 1 joined the set when task B7 widened
    # its caption and added a row; the conversion is what SHOULD happen, and
    # the count is asserted so that a table silently STAYING a longtable --
    # which is how a page folio once landed in the body text -- still fails.
    "tables_rewidthed": 3,
    "displays": 16,
    "spans_code": 0,         # v5 quotes no file paths in the body
    "stripped": 0,           # v5 carries no internal-only section
}
PAGES_MAIN = (18, 26)      # tolerance, not a target
PAGES_SUPP = (1, 2)


def checkCounts(info: dict) -> None:
    wrong = {k: (info[k], v) for k, v in EXPECTED.items() if info.get(k) != v}
    if wrong:
        raise SystemExit("the document does not contain what it should:\n  "
                         + "\n  ".join(f"{k}: got {g}, expected {e}"
                                       for k, (g, e) in wrong.items()))


def checkBalanced(tex: str, name: str) -> None:
    """an unclosed \\caption{ is how a whole float turned into escaped text."""
    for env in ("figure", "longtable", "equation*", "align*"):
        b, e = tex.count(rf"\begin{{{env}}}"), tex.count(rf"\end{{{env}}}")
        if b != e:
            raise SystemExit(f"{name}: {env} unbalanced, {b} begin / {e} end")
    for m in re.finditer(r"\\caption\{", tex):
        depth, i = 0, m.end() - 1
        while i < len(tex):
            if tex[i] == "{" and tex[i - 1] != "\\":
                depth += 1
            elif tex[i] == "}" and tex[i - 1] != "\\":
                depth -= 1
                if depth == 0:
                    break
            i += 1
        else:
            raise SystemExit(f"{name}: unclosed \\caption{{ at offset {m.start()}")


def build() -> dict:
    md = SRC.read_text()
    abstract_md, meta_md, body_md = splitFrontMatter(md)
    supp_md, body_md = splitSupplementary(body_md)

    body, info = _renderBody(body_md, headings=True)
    supp, sinfo = _renderBody(supp_md, headings=False)
    abstract_md, n_am, n_ac = convertSpans(abstract_md)
    meta_md, _, _ = convertSpans(meta_md)
    abstract, meta = _pandoc(abstract_md), _pandoc(meta_md)

    chars = residualUnicode(body + abstract + meta)
    front = (r"\abstract{" + abstract.strip() + "}" + "\n\n" + meta.strip() + "\n")
    OUT_TEX.write_text(_document(chars, TITLE, front, body))

    # the supplementary numbers its figures S1, S2, ... so they cannot collide
    schars = residualUnicode(supp)
    supp_front = (r"\abstract{Supplementary material accompanying the paper above. "
                  r"Every number in the caption is recomputed by the script named in "
                  r"it and asserted by the accompanying test suite.}" + "\n")
    OUT_SUPP_TEX.write_text(_document(
        schars, TITLE, supp_front, supp,
        extra_preamble=(r"\renewcommand{\thefigure}{S\arabic{figure}}" + "\n"
                        # one figure on a two-page document should not be pushed
                        # to a page of its own by the default float fractions
                        + r"\renewcommand{\topfraction}{0.92}" + "\n"
                        + r"\renewcommand{\textfraction}{0.06}" + "\n"
                        + r"\renewcommand{\floatpagefraction}{0.85}" + "\n")))

    checkBalanced(OUT_TEX.read_text(), "main")
    checkBalanced(OUT_SUPP_TEX.read_text(), "supplementary")
    return {"spans_math": info["spans_math"] + sinfo["spans_math"] + n_am,
            "spans_code": info["spans_code"] + sinfo["spans_code"] + n_ac,
            "sections": info["sections"],
            "figures_main": info["figures"], "figures_supp": sinfo["figures"],
            "tables": info["tables"], "displays": info["displays"],
            "tables_rewidthed": info["tables_rewidthed"],
            "stripped": info["stripped"], "unicode_mapped": len(chars)}


# pdflatex stamps /CreationDate, /ModDate and a random /ID into every PDF, so two
# builds of identical input differ byte-wise and every rebuild is a binary diff.
# The figures already strip their timestamps; SOURCE_DATE_EPOCH does the same here.
# The value is arbitrary and fixed -- it only has to be stable.
_EPOCH = "1785801600"     # 2026-08-04T00:00:00Z


def compile(tex: Path) -> int:
    """two pdflatex passes, for hyperref and any cross-references."""
    import os
    texinputs = str(REPO_ROOT / "manuscript" / "springer") + ":"
    env = dict(os.environ, SOURCE_DATE_EPOCH=_EPOCH, FORCE_SOURCE_DATE="1",
               TEXINPUTS=texinputs + os.environ.get("TEXINPUTS", ""))
    for _ in range(2):
        r = subprocess.run(["pdflatex", "-interaction=nonstopmode",
                            "-halt-on-error", "-output-directory",
                            str(tex.parent), str(tex)],
                           capture_output=True, text=True, env=env)
        if r.returncode != 0:
            print("\n".join(r.stdout.splitlines()[-40:]))
            print("\n--- errors ---")
            print("\n".join(l for l in r.stdout.splitlines() if l.startswith("!")))
            return r.returncode
    for ext in (".aux", ".log", ".out", ".toc"):
        p = tex.with_suffix(ext)
        if p.exists():
            p.unlink()
    return 0


def checkPages() -> None:
    import shutil
    if not shutil.which("pdfinfo"):
        return
    for pdf, (lo, hi) in ((OUT_TEX.with_suffix(".pdf"), PAGES_MAIN),
                          (OUT_SUPP_TEX.with_suffix(".pdf"), PAGES_SUPP)):
        out = subprocess.run(["pdfinfo", str(pdf)], capture_output=True,
                             text=True).stdout
        n = int(re.search(r"Pages:\s+(\d+)", out).group(1))
        if not lo <= n <= hi:
            raise SystemExit(f"{pdf.name}: {n} pages, expected {lo}-{hi}")
        print(f"  {pdf.name:32s} {n} pages, "
              f"{pdf.stat().st_size // 1024} kB")


if __name__ == "__main__":
    info = build()
    checkCounts(info)
    for k, v in info.items():
        print(f"  {k:16s} {v}")
    rc = 0
    for tex in (OUT_TEX, OUT_SUPP_TEX):
        rc |= compile(tex)
    if rc == 0:
        checkPages()
    sys.exit(rc)
