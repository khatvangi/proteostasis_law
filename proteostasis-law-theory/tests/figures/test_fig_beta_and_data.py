"""Figure 4's caption numbers, and the provenance of data/figures/."""

import hashlib
import json
import sys
import unittest
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (_REPO_ROOT / "scripts", _REPO_ROOT / "scripts" / "phase3",
           _REPO_ROOT / "scripts" / "figures"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import asymmetric_division as A  # noqa: E402
import fig_beta  # noqa: E402

_MS = (_REPO_ROOT / "manuscript" / "MANUSCRIPT_BMB_v5.md").read_text()
_DATA = _REPO_ROOT / "data" / "figures"
_MANIFEST = _DATA / "MANIFEST.json"


class TestBetaFigureCaption(unittest.TestCase):
    def testTheThreeQuotedIntervalsAreRecomputed(self):
        """damping is per-beta; a single value put every row ~1.5% off (D042)."""
        for beta, want in ((1.0, "0.05–0.08"), (0.5, "0.09–0.16"),
                           (0.25, "0.18–0.31")):
            lo, hi = A.requiredAggregateFractionBeta(beta, A.dampingAtBeta(beta))
            self.assertEqual(f"{100*lo:.2f}–{100*hi:.2f}", want)
            self.assertIn(want, _MS)

    def testEveryStoredDampingIsInTheMeasuredRange(self):
        for beta, d in A.DAMPING_BY_BETA:
            self.assertGreaterEqual(d, 0.346)
            self.assertLessEqual(d, 0.356)

    def testCaptionRefusesALowerBoundOnBeta(self):
        cap = _MS[_MS.index("**Fig. 5**"):]
        cap = cap[:cap.index("\n\n")]
        self.assertIn("No lower limit on `β` is drawn", cap)
        self.assertIn("46.5", cap)
        self.assertNotIn("0.145", cap)

    def testCaptionCallsTheWildTypeFigureABound(self):
        cap = _MS[_MS.index("**Fig. 5**"):]
        cap = cap[:cap.index("\n\n")]
        self.assertIn("a bound and not a value", cap)


@unittest.skipUnless(_MANIFEST.is_file(), "data/figures/ not built")
class TestFigureDataProvenance(unittest.TestCase):
    def setUp(self):
        self.m = json.loads(_MANIFEST.read_text())

    def testHashesMatch(self):
        for e in self.m["files"]:
            got = hashlib.sha256((_DATA / e["file"]).read_bytes()).hexdigest()
            self.assertEqual(got, e["sha256"], f"{e['file']} changed")

    def testEveryFileCarriesItsPopulationAndCount(self):
        for e in self.m["files"]:
            self.assertIn("population", e)
            self.assertIn("n_states", e)
            self.assertIn("is_subsample", e)
            self.assertGreater(e["n_states"], 0)

    def testNothingIsASubsample(self):
        """D036: three headline numbers came from an undeclared 20-state sample."""
        for e in self.m["files"]:
            self.assertFalse(e["is_subsample"], f"{e['file']} is a subsample")

    def testThePopulationsAreDistinguishedAndSizedAsTheManuscriptStates(self):
        by = {e["file"]: e for e in self.m["files"]}
        self.assertEqual(by["saturation.tsv"]["n_states"], 2884)
        self.assertEqual(by["saturation.tsv"]["population"], "kinetic_box")
        self.assertEqual(by["identity.tsv"]["n_states"], 325)
        self.assertEqual(by["identity.tsv"]["population"], "load_grid")
        self.assertIn("load grid", _MS)
        self.assertIn("kinetic box", _MS)


class TestSectionFiveNamesItsPopulations(unittest.TestCase):
    def testEveryTableRowNamesAPopulation(self):
        s = _MS[_MS.index("## 5. Numerical verification"):_MS.index("## 6.")]
        self.assertIn("| quantity | population | median | p99 | max |", s)
        self.assertIn("load grid, 325", s)
        self.assertIn("kinetic box, 2884", s)

    def testEveryMaximumIsReportedBesideAP99(self):
        """an extremum grows with population size; a p99 does not (D037)."""
        s = _MS[_MS.index("## 5. Numerical verification"):_MS.index("## 6.")]
        table = s[s.index("| quantity |"):s.index("\n\n", s.index("| quantity |"))]
        for row in table.split("\n"):
            if not row.startswith("|") or "---" in row or "quantity" in row:
                continue
            cells = [c.strip() for c in row.strip("|").split("|")]
            # population column must never be blank
            self.assertTrue(cells[1], f"row names no population: {row}")
            # if a max is given, a p99 must be given too
            if cells[-1] != "—":
                self.assertNotEqual(cells[-2], "—",
                                    f"max without a p99 beside it: {row}")
        self.assertIn("stable under resampling", s)
        self.assertIn("7.56×10⁻⁷", table)

    def testTheTableCarriesTheFullPopulationValues(self):
        s = _MS[_MS.index("## 5. Numerical verification"):_MS.index("## 6.")]
        table = s[s.index("| quantity |"):s.index("\n\n", s.index("| quantity |"))]
        self.assertIn("+0.9960", table)
        self.assertIn("1.63×10⁻⁹", table)
        self.assertNotIn("+0.9987", table)
        self.assertNotIn("8.2×10⁻¹⁰", table)

    def testTheCorrectionsAreStatedRatherThanSilentlySwapped(self):
        """a section about exactness should survive arithmetic done by a reader."""
        s = _MS[_MS.index("## 5. Numerical verification"):_MS.index("## 6.")]
        self.assertIn("20-state random subsample", s)
        self.assertIn("None of these changes weakens any claim", s)


if __name__ == "__main__":
    unittest.main()


class TestSaturationFigureMatchesSectionSix(unittest.TestCase):
    """the screen that would have split figure from text is not applied (D039)."""

    def testNoScreenAndNoExclusion(self):
        import fig_saturation
        self.assertEqual(fig_saturation.S_A_ZERO, 0.0)

    def testMediansReproduceTheTextExactly(self):
        import numpy as np
        import pandas as pd
        d = pd.read_csv(_REPO_ROOT / "data" / "figures" / "saturation.tsv",
                        sep="\t")
        for col, want in (("s_ref", 0.175), ("s_u", 0.155), ("s_a", 0.056)):
            self.assertAlmostEqual(float(np.median(d[col])), want, places=3)
            self.assertIn(f"{want:.3f}", _MS)

    def testTheWidthsAreQuotedWithTheirPopulation(self):
        """the kinetic-box widths moved into the caption; the regulation ones did not.

        Section 9 used to restate the 0.876 width that the panel shows better, so it
        now points at the figure instead. The REGULATION widths stay in the prose,
        because they belong to a different experiment that no figure carries -- and
        that distinction is exactly what the sixth affected number turned on (D039).
        """
        cap = _MS[_MS.index("**Fig. 3**"):]
        cap = cap[:cap.index("\n\n")]
        self.assertIn("all 2884 folds of the kinetic box", cap)
        self.assertIn("0.876", cap)
        self.assertNotIn("0.876", _MS[_MS.index("## 9. Predictions"):])
        self.assertIn("that experiment's own 30 networks", _MS)
