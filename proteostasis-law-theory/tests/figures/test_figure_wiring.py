"""how the figures attach to the prose, not what they contain.

Three things are pinned here and nothing else is:

1. NUMBERING FOLLOWS FIRST MENTION. The dilution figure was built last and is
   cited in section 4.2, so it is Figure 2 and everything after it shifted. Script
   filenames are deliberately semantic (`fig_dilution.py`, not `fig5.py`) so that
   a future reorder touches one FIGURE constant per script instead of six files.
2. EVERY FIGURE IS CITED IN PROSE before it is embedded. A referee should never
   meet a panel they were not sent to.
3. NO CAPTION NUMBER DISAGREES WITH THE TEXT NUMBER for the same quantity, and
   both are asserted against the generator rather than against each other. Caption
   drifting from text by one digit is the specific failure that got through on P3.
"""

from __future__ import annotations

import re
import sys
import unittest
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (_REPO_ROOT / "scripts", _REPO_ROOT / "scripts" / "figures",
           _REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import fig_theorem, fig_dilution, fig_saturation  # noqa: E402
import fig_front, fig_beta, fig_identity  # noqa: E402

_DOC = (_REPO_ROOT / "manuscript" / "bmb_v4.md").read_text()

# script module -> the number it claims in the manuscript
_SCRIPTS = {
    "fig_theorem": fig_theorem, "fig_dilution": fig_dilution,
    "fig_saturation": fig_saturation, "fig_front": fig_front,
    "fig_beta": fig_beta, "fig_identity": fig_identity,
}
_MAIN = ("fig1", "fig2", "fig3", "fig4", "fig5")


def _proseMentions(n: str) -> list:
    """offsets of every citation of figure `n` OUTSIDE its own caption.

    `Fig. 1a` counts as citing Figure 1, so the lookahead excludes only a further
    digit; and the caption is `**Fig. N**`, whose two asterisks sit immediately
    before the match, with a newline before them.
    """
    return [m.start() for m in re.finditer(rf"Fig\. {n}(?![0-9])", _DOC)
            if _DOC[max(0, m.start() - 2): m.start()] != "**"]


def _captionOf(stem: str) -> str:
    """the caption paragraph following the embed of `stem`."""
    label = "S1" if stem == "figS1" else stem[3:]
    marker = f"**Fig. {label}**"
    start = _DOC.index(marker)
    return _DOC[start: _DOC.index("\n\n", start)]


class TestNumberingFollowsFirstMention(unittest.TestCase):
    def testEachScriptDeclaresExactlyOneFigureNumber(self):
        claimed = sorted(m.FIGURE for m in _SCRIPTS.values())
        self.assertEqual(claimed, sorted(list(_MAIN) + ["figS1"]))

    def testEmbedsAppearInNumericOrder(self):
        embeds = re.findall(r"!\[Figure ([0-9S]+)\]\(\.\./figures/(fig[0-9S]+)\.pdf\)",
                            _DOC)
        labels = [e[0] for e in embeds]
        self.assertEqual(labels, ["1", "2", "3", "4", "5", "S1"])
        for label, stem in embeds:
            self.assertEqual(stem, f"fig{label}", "embed path disagrees with label")

    def testFirstProseMentionOrdersTheMainFigures(self):
        """the rule: numbered by order of first mention, captions excluded."""
        firsts = []
        for stem in _MAIN:
            n = stem[3:]
            hits = _proseMentions(n)
            self.assertTrue(hits, f"Fig. {n} is never cited in prose")
            firsts.append(min(hits))
        self.assertEqual(firsts, sorted(firsts),
                         "first prose mentions are out of numeric order")

    def testEveryFigureIsCitedBeforeItIsEmbedded(self):
        for stem in _MAIN:
            n = stem[3:]
            embed = _DOC.index(f"![Figure {n}](")
            cite = min(_proseMentions(n))
            self.assertLess(cite, embed,
                            f"Fig. {n} is embedded before the text sends anyone to it")

    def testSupplementaryIsCitedInSectionFive(self):
        cite = _DOC.index("Fig. S1")
        self.assertLess(cite, _DOC.index("## Supplementary figure"))


class TestCaptionsAgreeWithTextAndWithTheGenerator(unittest.TestCase):
    """every shared number checked against the code, never caption-against-text."""

    def testDilutionBandMatchesBothPlaces(self):
        o = fig_dilution.build()
        band = f"{o['phi_lo']:.4f}–{o['phi_hi']:.4f}"
        self.assertEqual(band, "0.1245–0.1343")
        self.assertIn(band, _DOC, "section 4.2 does not carry the computed band")
        self.assertAlmostEqual(o["delta_hi"], 0.3915, places=4)
        self.assertIn("0 → 0.39", _DOC)

    def testSaturationWidthsAndFloorLadderMatch(self):
        o = fig_saturation.build()
        cap = _captionOf("fig3")
        for key in ("$s_{\\mathrm{ref}}$", "$s_u$", "$s_a$"):
            w = o["stats"][key]["width"]
            self.assertIn(f"{w:.3f}", cap, f"width {w:.3f} missing from the caption")
        # the ladder is quoted in section 6's prose; endpoints must be the computed ones
        floors = dict(o["sensitivity"])
        self.assertIn(f"{floors[1e-4]:.3f}", _DOC)
        self.assertIn(f"{floors[2e-2]:.3f}", _DOC)
        # medians live in the table only -- the caption must not restate them
        self.assertNotIn("0.175", cap)
        self.assertNotIn("0.056", cap)

    def testFrontNumbersMatchBothPlaces(self):
        o = fig_front.build()
        cap = _captionOf("fig4")
        for tok in (f"{o['front_lo']:.3f}", f"{o['front_hi']:.3f}",
                    f"{o['exact_ratio']:.6f}", f"{o['grid_best_ratio']:.4f}"):
            self.assertIn(tok, cap, f"{tok} missing from the Fig. 4 caption")
            self.assertIn(tok, _DOC, f"{tok} missing from the text")

    def testBetaBandMatchesTableAndCaption(self):
        """caption rounds to 2 dp and the section 8.4 table to 3; both from one solve."""
        o = fig_beta.build()
        cap = _captionOf("fig5")
        # the generator returns proteome FRACTIONS; both table and caption are in %
        lo1, hi1 = (100 * v for v in o["at_beta_1"])
        lo25, hi25 = (100 * v for v in o["at_beta_025"])
        self.assertIn(f"{lo1:.2f}–{hi1:.2f} %", cap)
        self.assertIn(f"{lo25:.2f}–{hi25:.2f} %", cap)
        # the table carries the same solve at two more digits
        self.assertIn(f"{lo1:.4f}", _DOC)
        self.assertIn(f"{hi25:.4f}", _DOC)
        # the prose quotes a closest approach over the WHOLE plotted range, not at
        # the marked beta -- an earlier draft quoted the marked value and was wrong
        # by a factor of five at the left edge of its own figure.
        self.assertAlmostEqual(o["closest_ratio"], 3.19, places=2)
        self.assertAlmostEqual(o["widest_ratio"], 214.0, delta=1.0)
        self.assertIn("between 3× and 214× below", _DOC)
        # 15.94, which the table rounds to 16x -- the prose says "at least
        # fifteenfold" rather than sixteen, because 15.94 is not 16
        self.assertGreater(100 * o["rpoH"][0] / hi25, 15.0)
        self.assertLess(100 * o["rpoH"][0] / hi25, 16.0)
        self.assertIn("at least fifteenfold below", _DOC)

    def testIdentityCaptionMatchesSectionFive(self):
        import pandas as pd
        df = pd.read_csv(_REPO_ROOT / "data/figures/identity.tsv", sep="\t")
        c = fig_identity.captionNumbers(df)
        cap = _captionOf("figS1")
        sec5 = _DOC[_DOC.index("## 5. Numerical verification"):
                    _DOC.index("## 6. Where the boundary sits")]
        for tok in ("+0.9960", "−0.835", "+0.041", "1.54×10⁻²", "1.29×10⁻⁹"):
            self.assertIn(tok, cap, f"{tok} missing from the S1 caption")
            self.assertIn(tok, sec5, f"{tok} missing from section 5")
        self.assertAlmostEqual(c["corr_parallelism"], 0.9960, places=4)
        self.assertAlmostEqual(c["max_corr_loglog"], -0.835, places=3)
        # the retracted pair must not survive in the caption at all
        for tok in ("−0.262", "+0.060"):
            self.assertNotIn(tok, cap)


class TestPhysicalSize(unittest.TestCase):
    """BMB: 84 mm double-column or 174 mm single-column, never taller than 234 mm."""

    def testEveryFigureIsAPermittedWidthAndWithinHeight(self):
        import matplotlib
        matplotlib.use("Agg")
        from matplotlib.backends.backend_pdf import PdfPages  # noqa: F401
        from pypdf import PdfReader
        mm = 25.4 / 72.0     # pdf points -> mm
        for stem in list(_MAIN) + ["figS1"]:
            pdf = _REPO_ROOT / "figures" / f"{stem}.pdf"
            self.assertTrue(pdf.exists(), f"{stem}.pdf not built")
            box = PdfReader(str(pdf)).pages[0].mediabox
            w, h = float(box.width) * mm, float(box.height) * mm
            self.assertLessEqual(h, 234.0 + 1.0, f"{stem} is {h:.1f} mm tall")
            self.assertTrue(w <= 84.0 + 1.5 or abs(w - 174.0) <= 1.5,
                            f"{stem} is {w:.1f} mm wide, neither 84 nor 174")


if __name__ == "__main__":
    unittest.main()
