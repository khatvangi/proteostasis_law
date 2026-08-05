"""figures 3 and S1, and the section-5 numbers figure S1 now owns (D041).

The point of the S1 tests is NOT that the figure renders. It is that every number
in section 5's normalisation paragraph is asserted against `captionNumbers`, which
recomputes it, rather than against a value stored anywhere. The pair -0.262 /
+0.060 survived in the manuscript precisely because no test could tell the
difference between a number and a number-shaped string.
"""

from __future__ import annotations

import re
import sys
import unittest
from pathlib import Path

import pandas as pd

_REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (_REPO_ROOT / "scripts", _REPO_ROOT / "scripts" / "figures",
           _REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import fig_front  # noqa: E402
import fig_identity  # noqa: E402

_MANUSCRIPT = (_REPO_ROOT / "manuscript" / "MANUSCRIPT_BMB_v5.md").read_text()


def _section(doc: str, header: str) -> str:
    start = doc.index(header)
    rest = doc[start + len(header):]
    cut = rest.find("\n## ")
    return rest if cut < 0 else rest[:cut]


class TestFigureS1OwnsItsNumbers(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        df = pd.read_csv(_REPO_ROOT / "data/figures/identity.tsv", sep="\t")
        cls.c = fig_identity.captionNumbers(df)
        cls.sec = _section(_MANUSCRIPT, "\n## 5. Numerical verification")

    def testThePopulationIsCompleteAndNamed(self):
        self.assertEqual(self.c["n"], 325)
        self.assertIn("load grid, 325", self.sec)

    def testTheParallelismCorrelationMatchesTheComputation(self):
        self.assertAlmostEqual(self.c["corr_parallelism"], 0.9960, places=4)
        self.assertIn("+0.9960", self.sec)

    def testTheNormalisationContrastMatchesTheComputation(self):
        """the exact numbers D041 corrected, asserted against the generator."""
        self.assertAlmostEqual(self.c["max_corr_loglog"], -0.835, places=3)
        self.assertAlmostEqual(self.c["grad_corr_loglog"], +0.041, places=3)
        self.assertAlmostEqual(self.c["max_corr_raw"], -0.221, places=3)
        self.assertAlmostEqual(self.c["grad_corr_raw"], +0.073, places=3)
        for token in ("−0.835", "+0.041", "−0.221", "+0.073"):
            self.assertIn(token, self.sec, msg=f"{token} missing from section 5")

    def testTheUnreproducibleValuesAppearOnlyAsHistory(self):
        """pin the property: banned as a current claim, required as a correction."""
        for tok in ("−0.262", "+0.060"):
            if tok not in self.sec:
                continue
            idx = self.sec.index(tok)
            window = self.sec[max(0, idx - 400): idx]
            self.assertRegex(
                window, r"earlier draft|previously reported|superseded",
                msg=f"{tok} appears without being marked as superseded",
            )

    def testTheMedianSplitStatesTheDefectWithoutACorrelation(self):
        self.assertAlmostEqual(self.c["max_tighter_half"], 3.133e-07, places=9)
        self.assertAlmostEqual(self.c["max_looser_half"], 1.455e-07, places=9)
        self.assertGreater(self.c["max_degradation"], 2.0)
        self.assertIn("3.13×10⁻⁷", self.sec)
        self.assertIn("1.46×10⁻⁷", self.sec)

    def testTheTightestBracketIsTheFullPopulationsAndNotTheDraws(self):
        """a minimum written as a definite description is still a minimum."""
        self.assertAlmostEqual(self.c["tightest_eig"], 8.0998e-10, places=13)
        self.assertAlmostEqual(self.c["tightest_sin"], 7.7469e-09, places=12)
        self.assertIn("8.10×10⁻¹⁰", self.sec)
        self.assertIn("7.75×10⁻⁹", self.sec)
        # the subsample's minimum must not stand as a current claim
        self.assertNotRegex(self.sec, r"single state bracketed")

    def testTheCorrectionsParagraphCountsWhatItLists(self):
        """the count said three and listed four for two revisions running."""
        start = self.sec.index("Five of these values")
        para = self.sec[start: self.sec.index("\n\n", start)]
        # each correction is written as "<quantity> from <old> to <new>"
        moved = re.findall(r"\bfrom\b.{0,60}?\bto\b", para)
        self.assertEqual(len(moved), 5,
                         f"paragraph says five and lists {len(moved)}: {moved}")

    def testTheMaximaAreTheFullPopulationsAndSlopeIsReported(self):
        self.assertAlmostEqual(self.c["grad_max"], 1.2924e-09, places=12)
        self.assertAlmostEqual(self.c["max_max"], 1.5448e-02, places=6)
        self.assertIn("1.54×10⁻²", self.sec)
        self.assertIn("slope is 1.00", self.sec)


class TestFigure3ReconcilesWithSectionSeven(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.o = fig_front.build()

    def testTheFrontRangeMatchesTheText(self):
        self.assertAlmostEqual(self.o["front_lo"], 0.227, places=3)
        self.assertAlmostEqual(self.o["front_hi"], 0.965, places=3)
        self.assertIn("0.227 to 0.965", _MANUSCRIPT)

    def testTheOptimumIsSolvedOnTheBoundaryNotReadOffTheGrid(self):
        self.assertAlmostEqual(self.o["exact_ratio"], 1.0, places=6)
        # the grid cannot supply this number, which is why the figure does not ask it to
        self.assertLess(self.o["grid_best_ratio"], 0.90)
        self.assertGreater(self.o["exact_throughput"], self.o["grid_best_throughput"])

    def testTheFrontIsTheWholeFeasibleSet(self):
        self.assertEqual(self.o["n_feasible"], 469)
        self.assertEqual(self.o["n_front"], 13)


class TestFiguresAreDeterministic(unittest.TestCase):
    def testRebuildingGivesIdenticalBytes(self):
        for mod in (fig_front, fig_identity):
            first = mod.build()["hashes"]
            second = mod.build()["hashes"]
            self.assertEqual(first, second, f"{mod.__name__} is not byte-stable")


if __name__ == "__main__":
    unittest.main()
