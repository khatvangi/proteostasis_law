"""the markdown -> LaTeX conversion, and the corruptions it must not repeat.

Each test here corresponds to something that actually went wrong while building
the PDF, and every one of them was SILENT: the build reported success and the
damage was only visible in the rendered pages.
"""

from __future__ import annotations

import subprocess
import sys
import unittest
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "scripts" / "manuscript"))

import to_latex as T  # noqa: E402


class TestMathConversion(unittest.TestCase):
    def testMultiCharacterSubscriptsAreBraced(self):
        """`j_crit` unbraced means j_c followed by 'rit'. Silent corruption."""
        self.assertIn("j_{crit}", T.mathify("j_crit"))
        self.assertIn("s_{ref}", T.mathify("s_ref"))
        self.assertIn("C_{enz}", T.mathify("C_enz"))
        # a single-character subscript needs no braces and must not gain any
        self.assertEqual(T.mathify("C_0"), "$C_0$")

    def testNoSpaceBeforeTheClosingDollar(self):
        """pandoc silently escapes both dollars instead of parsing math."""
        for span in ("μ", "δ", "φ_enz", "det J = −(∇R × ∇G)", "μ₀ = 0.3"):
            out = T.mathify(span)
            self.assertTrue(out.startswith("$") and out.endswith("$"))
            self.assertFalse(out[-2].isspace(), f"space before $ in {out!r}")
            self.assertFalse(out[1].isspace(), f"space after $ in {out!r}")

    def testExistingScriptGroupsSurviveBraceEscaping(self):
        """`f_{ι,I}` carries its own group; escaping set braces destroyed it."""
        self.assertIn(r"f_{\iota", T.mathify("x'_o = ±f_{ι,I} det(H)/det(J)"))
        self.assertIn(r"\{G = 0\}", T.mathify("{G = 0}"))

    def testNamedOperatorsAreUpright(self):
        self.assertIn(r"\det", T.mathify("det J"))
        self.assertIn(r"\sin", T.mathify("sin θ"))
        self.assertIn(r"\max", T.mathify("max(|det J|, |cross|)"))

    def testStarredEquilibria(self):
        self.assertEqual(T.mathify("u*"), "$u^{*}$")
        self.assertIn("u^{*}", T.mathify("j_crit = R(u*, a*)"))


class TestSpanClassification(unittest.TestCase):
    def testTheManuscriptClassifiesWithoutGuessing(self):
        """a new code-shaped span must fail the build, not render as math."""
        md = (_REPO_ROOT / "manuscript" / "bmb_v4.md").read_text()
        _, n_math, n_code = T.convertSpans(md)
        self.assertEqual(n_code, len(T.CODE_SPANS))
        self.assertGreater(n_math, 250)

    def testAnUndeclaredCodeSpanRaises(self):
        with self.assertRaises(SystemExit):
            T.convertSpans("see `scripts/phase3/newthing.py` for details")


class TestWholeBuild(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.info = T.build()
        cls.main = T.OUT_TEX.read_text()
        cls.supp = T.OUT_SUPP_TEX.read_text()

    def testFiguresAreSplitBetweenTheTwoDocuments(self):
        self.assertEqual(self.info["figures_main"], 5)
        self.assertEqual(self.info["figures_supp"], 1)

    def testNoRawBlockWasEscapedByPandoc(self):
        """an unescaped % in a caption orphaned \\centering and centred §10 on."""
        for doc, name in ((self.main, "main"), (self.supp, "supplementary")):
            self.assertNotIn(r"\textbackslash begin", doc, f"{name} has an "
                             "escaped raw block")
            self.assertEqual(doc.count(r"\begin{figure}"),
                             doc.count(r"\end{figure}"), f"{name} floats unbalanced")
            self.assertEqual(doc.count(r"\centering"), doc.count(r"\begin{figure}"),
                             f"{name} has a \\centering outside a float")

    def testFiguresAreNotRescaled(self):
        """84 and 174 mm are the journal's widths; \\linewidth threw them away."""
        self.assertIn(r"\includegraphics[max width=\linewidth]", self.main)
        self.assertNotIn(r"width=\linewidth,", self.main)

    def testSupplementaryNumbersItsFiguresSeparately(self):
        self.assertIn(r"\renewcommand{\thefigure}{S\arabic{figure}}", self.supp)
        self.assertNotIn(r"\thefigure", self.main)

    def testEveryHeadingLevelAndNumberSurvived(self):
        self.assertEqual(self.info["sections"], 10)
        self.assertEqual(self.main.count(r"\section{"), 10)

    def testNoUnmappedUnicodeReachesPdflatex(self):
        for doc in (self.main, self.supp):
            declared = {c for c in doc if ord(c) > 127}
            body = doc[doc.index(r"\begin{document}"):]
            used = {c for c in body if ord(c) > 127}
            self.assertTrue(used <= declared)
            for c in used:
                self.assertIn(f"\\newunicodechar{{{c}}}", doc,
                              f"U+{ord(c):04X} {c!r} has no mapping")

    def testBothPdfsExistAndAreNonTrivial(self):
        for pdf in (T.OUT_TEX.with_suffix(".pdf"),
                    T.OUT_SUPP_TEX.with_suffix(".pdf")):
            self.assertTrue(pdf.exists(), f"{pdf.name} not built")
            self.assertGreater(pdf.stat().st_size, 50_000)

    def testThePdfBuildIsReproducible(self):
        """pdflatex stamps a date and a random /ID; both are pinned."""
        src = (_REPO_ROOT / "scripts" / "manuscript" / "to_latex.py").read_text()
        self.assertIn("SOURCE_DATE_EPOCH", src)
        self.assertIn("FORCE_SOURCE_DATE", src)

    def testTheGeneratedTexSaysItIsGenerated(self):
        for doc in (self.main, self.supp):
            self.assertIn("GENERATED FILE -- do not edit", doc)


if __name__ == "__main__":
    unittest.main()
