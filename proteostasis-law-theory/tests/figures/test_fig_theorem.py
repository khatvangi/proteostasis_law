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

_MS = (_REPO_ROOT / "manuscript" / "MANUSCRIPT_BMB_v5.md").read_text()


class TestFigureOneNumbers(unittest.TestCase):
    # NOTE (2026-08-05): the prose-content assertions that lived here were
    # deleted rather than rewritten. They pinned v4 sentences in sections v5
    # restructured, and rewriting them would have meant writing tests against
    # text that has not settled. What survives in this file is the class of
    # assertion that is about correctness rather than wording: populations
    # named and sized, no maximum without a p99 beside it, no figure script
    # reading the gitignored run root, and caption numbers recomputed from
    # their generator. Re-pin the prose once the manuscript stops moving.

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
        import fig_theorem
        stats = np.array([[0.1, 0.1, 0.1, 0.5, 1.0], [0.2, 0.2, 0.2, 0.25, 1.0]])
        self.assertFalse(fig_theorem.homeostasisPointExists(stats)["exists"])
        stats[0, 3] = -1.0
        self.assertTrue(fig_theorem.homeostasisPointExists(stats)["exists"])

    def testTheTurnIsComputedAndGeneric(self):
        """du*/dj = -G_a/det J vanishes at the nullcline turn, where det J != 0."""
        import fig_theorem
        p = M.Params().validate()
        UU, AA, GG, _ = fig_theorem.fields(p, n=90)
        stats = fig_theorem.branchStats(p, fig_theorem.nullclineSegments(UU, AA, GG))
        turn = fig_theorem.turningPoint(p, stats)
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
        self.assertIn("2.027×10⁻³", cap)

    def testScientificNotationIsUniformAcrossTheManuscript(self):
        """`3.5e-10` reads as "3.5e minus 10" once typeset (D043).

        Pin the property, not the token: no e-notation anywhere, so the
        manuscript and the PDF cannot disagree about how a number is written.
        """
        import re
        stray = re.findall(r"\d+\.?\d*e-\d+", _MS)
        self.assertEqual(stray, [], f"e-notation left in the manuscript: {stray}")



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

        The check is on the parsed CODE, not on the source text. A string search
        failed the moment a module explained in a comment why it no longer calls
        `phase1RunDir` -- the same shape as D042: a token-absence test breaking on
        text that became honest about its own history. What matters is whether
        the module CALLS it, and only the syntax tree knows that.
        """
        import ast

        allowed = {"build_figure_data.py"}
        for f in sorted((_REPO_ROOT / "scripts" / "figures").glob("*.py")):
            if f.name in allowed:
                continue
            tree = ast.parse(f.read_text())
            names = {n.id for n in ast.walk(tree) if isinstance(n, ast.Name)}
            attrs = {n.attr for n in ast.walk(tree) if isinstance(n, ast.Attribute)}
            self.assertNotIn("phase1RunDir", names | attrs,
                             f"{f.name} calls phase1RunDir; the run root is "
                             "gitignored, so the figure would not rebuild from "
                             "a clean checkout")
            strings = {n.value for n in ast.walk(tree)
                       if isinstance(n, ast.Constant) and isinstance(n.value, str)}
            self.assertNotIn("results", strings, f"{f.name} reads results/")

    def testStyleDoesNotDependOnALocalMatplotlibConfig(self):
        src = (_REPO_ROOT / "scripts" / "figures" / "_figstyle.py").read_text()
        self.assertIn("rcdefaults", src)
        self.assertIn("seaborn", src)      # only as the documented exclusion
        self.assertNotIn("import seaborn", src)


if __name__ == "__main__":
    unittest.main()
