"""Figure 1: the numbers in the caption, and the honesty of panel (b).

Captions on an earlier project drifted from the figures because they were written
by hand against figures generated separately. Every number quoted in Figure 1's
caption is recomputed here from the same functions the figure calls, and matched
against the caption text in the manuscript.
"""

import re
import sys
import unittest
from pathlib import Path

import numpy as np

_REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (_REPO_ROOT / "scripts", _REPO_ROOT / "scripts" / "phase3",
           _REPO_ROOT / "scripts" / "figures"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import fold_theorem as FT  # noqa: E402
import boundary_structure as B  # noqa: E402
import dilution as D  # noqa: E402
import _figstyle as F  # noqa: E402

_MS = (_REPO_ROOT / "manuscript" / "bmb_v4.md").read_text()


class TestFigureOneNumbers(unittest.TestCase):
    def setUp(self):
        self.p = M.Params().validate()
        self.fold = FT.foldSolve(self.p)

    def testFoldMatchesTheCaption(self):
        jc, us, as_ = self.fold
        self.assertAlmostEqual(jc, 0.154239, places=6)
        self.assertAlmostEqual(us, 0.416573, places=6)
        self.assertAlmostEqual(as_, 0.264969, places=6)
        self.assertIn("0.1542", _MS)
        self.assertIn("0.4166", _MS)
        self.assertIn("0.2650", _MS)

    def testGradientsAreParallelToTheQuotedPrecision(self):
        _, us, as_ = self.fold
        s = FT.determinantIdentity(us, as_, self.p)["sin_angle"]
        self.assertLess(s, 1e-9)
        self.assertIn("3.5", _MS)          # 3.5e-10, quoted in the caption

    def testCaptionDoesNotQuoteTheIllConditionedResidual(self):
        """`rel_err` is 0/0 at a saddle-node and returns 1.0 regardless (D027)."""
        _, us, as_ = self.fold
        self.assertAlmostEqual(
            FT.determinantIdentity(us, as_, self.p)["rel_err"], 1.0, places=6)
        cap = _MS[_MS.index("**Fig. 1**"):]
        cap = cap[:cap.index("\n\n")]
        self.assertNotIn("relative error", cap)


class TestPanelBIsHonest(unittest.TestCase):
    """no horizontal tangent exists, so none may be plotted as computed."""

    def testGuIsStrictlyPositiveOnTheNullcline(self):
        p = M.Params().validate()
        g = D.Growth(0.0)
        vals = []
        for u in np.geomspace(1e-3, 3.0, 60):
            for a in B.allNullclineRootsAt(u, p, g):
                if a <= 0.0:
                    continue
                vals.append(FT._centralGradient(FT.aggregateG, u, a, p)[0])
        v = np.array(vals)
        self.assertGreater(len(v), 50)
        self.assertGreater(v.min(), 0.0)

    def testTheScriptReportsNoHomeostasisPoint(self):
        import fig1
        stats = np.array([[0.1, 0.1, 0.1, 0.5, 1.0], [0.2, 0.2, 0.2, 0.25, 1.0]])
        self.assertFalse(fig1.homeostasisPointExists(stats)["exists"])
        stats[0, 3] = -1.0
        self.assertTrue(fig1.homeostasisPointExists(stats)["exists"])

    def testTheTurnIsComputedAndGeneric(self):
        """du*/dj = -G_a/det J vanishes at the nullcline turn, where det J != 0."""
        import fig1
        p = M.Params().validate()
        UU, AA, GG, _ = fig1.fields(p, n=90)
        stats = fig1.branchStats(p, fig1.nullclineSegments(UU, AA, GG))
        turn = fig1.turningPoint(p, stats)
        self.assertIsNotNone(turn)
        self.assertEqual(turn["n_sign_changes"], 1)
        self.assertLess(abs(turn["G_a"]), 1e-9)
        self.assertGreater(turn["det_J"], 0.0)          # generic, and stable
        self.assertAlmostEqual(turn["det_J"], turn["Ra_times_Gu"],
                               delta=1e-3 * abs(turn["det_J"]) + 1e-9)
        jc, _, _ = FT.foldSolve(p)
        self.assertLess(turn["j"], jc)                  # the turn precedes the fold

    def testCaptionUsesComputedDataNotASchematic(self):
        cap = _MS[_MS.index("**Fig. 1**"):]
        cap = cap[:cap.index("\n\n")]
        self.assertIn("not a schematic", cap)
        self.assertIn("0.9990", cap)
        self.assertIn("2.027e-03", cap)

    def testSectionFiveUsesTheNonDegradingResidual(self):
        """the max-normalised metric worsens as the bracket tightens (D027).

        The superseded values are NOT required to be absent -- §5 states them
        explicitly as corrections, which is the honest presentation. What must
        hold is that they never appear as a CURRENT claim, i.e. only inside the
        sentence that supersedes them, and never in the results table.
        """
        s = _MS[_MS.index("## 5. Numerical verification"):_MS.index("## 6.")]
        table = s[s.index("| quantity |"):s.index("\n\n", s.index("| quantity |"))]
        for old in ("1.436", "+0.9987", "8.2×10⁻¹⁰"):
            self.assertNotIn(old, table, f"{old} still in the results table")
        correction = s[s.index("previously reported"):]
        correction = correction[:correction.index("\n\n")]
        for old in ("1.436×10⁻⁷", "+0.9987", "8.2×10⁻¹⁰"):
            self.assertIn(old, correction, f"{old} quoted outside the correction")
        self.assertIn("2.34×10⁻¹⁰", table)
        self.assertIn("load grid, 325", table)

    def testSectionThreeOneUsesTheirVocabulary(self):
        s = _MS[_MS.index("### 3.1"):_MS.index("## 4.")]
        self.assertIn("Haldane", s)
        self.assertIn("det H = G_u", s)
        # output = u is outside their framework and must be flagged as such
        self.assertIn("we do not call it infinitesimal homeostasis", s)


class TestFigureHygiene(unittest.TestCase):
    def testSeedIsPinned(self):
        self.assertEqual(F.SEED, 20260804)

    def testWidthsMatchTheSubmissionGuidelines(self):
        self.assertAlmostEqual(F.W_DOUBLE / F.MM, 84.0, places=6)
        self.assertAlmostEqual(F.W_SINGLE / F.MM, 174.0, places=6)
        self.assertAlmostEqual(F.H_MAX / F.MM, 234.0, places=6)

    def testNoFigureScriptReadsFromResults(self):
        """the run root is gitignored; a clean checkout must still reproduce.

        `build_figure_data.py` is the single exception and the reason the rule
        holds for everything else: it reads the run root WHEN PRESENT and writes
        the reduced arrays into data/figures/, which is tracked. Figures read
        those. If it ever stops being the only exception, this test fails.
        """
        allowed = {"build_figure_data.py"}
        for f in sorted((_REPO_ROOT / "scripts" / "figures").glob("*.py")):
            if f.name in allowed:
                continue
            src = f.read_text()
            self.assertNotIn('"results"', src, f"{f.name} reads results/")
            self.assertNotIn("'results'", src, f"{f.name} reads results/")
            self.assertNotIn("phase1RunDir", src, f"{f.name} reads the run root")

    def testStyleDoesNotDependOnALocalMatplotlibConfig(self):
        src = (_REPO_ROOT / "scripts" / "figures" / "_figstyle.py").read_text()
        self.assertIn("rcdefaults", src)
        self.assertIn("seaborn", src)      # only as the documented exclusion
        self.assertNotIn("import seaborn", src)


if __name__ == "__main__":
    unittest.main()
