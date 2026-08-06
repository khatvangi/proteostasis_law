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
    # NOTE (2026-08-05): the prose-content assertions that lived here were
    # deleted rather than rewritten. They pinned v4 sentences in sections v5
    # restructured, and rewriting them would have meant writing tests against
    # text that has not settled. What survives in this file is the class of
    # assertion that is about correctness rather than wording: populations
    # named and sized, no maximum without a p99 beside it, no figure script
    # reading the gitignored run root, and caption numbers recomputed from
    # their generator. Re-pin the prose once the manuscript stops moving.

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


    def testTheTightestBracketIsTheFullPopulationsAndNotTheDraws(self):
        """a minimum written as a definite description is still a minimum."""
        self.assertAlmostEqual(self.c["tightest_eig"], 8.0998e-10, places=13)
        self.assertAlmostEqual(self.c["tightest_sin"], 7.7469e-09, places=12)
        self.assertIn("8.10×10⁻¹⁰", self.sec)
        self.assertIn("7.75×10⁻⁹", self.sec)
        # the subsample's minimum must not stand as a current claim
        self.assertNotRegex(self.sec, r"single state bracketed")



class TestFigure3ReconcilesWithSectionSeven(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.o = fig_front.build()

    def testTheFrontRangeMatchesTheText(self):
        self.assertAlmostEqual(self.o["front_lo"], 0.227, places=3)
        self.assertAlmostEqual(self.o["front_hi"], 0.965, places=3)
        # The caption stated this range twice, at three and at four decimals,
        # and this check pinned the three-decimal copy. Task R9 removed the
        # duplicate — two printings of one quantity are how prose lands in the
        # wrong place (D047) — so the surviving four-decimal statement is what
        # is asserted, and it is asserted once.
        self.assertIn(f"{self.o['front_lo']:.4f} to {self.o['front_hi']:.4f}",
                      _MANUSCRIPT)
        self.assertEqual(_MANUSCRIPT.count("0.2271 to 0.9652"), 1,
                         "the front range is printed more than once again")

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
