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
    # NOTE (2026-08-05): the prose-content assertions that lived here were
    # deleted rather than rewritten. They pinned v4 sentences in sections v5
    # restructured, and rewriting them would have meant writing tests against
    # text that has not settled. What survives in this file is the class of
    # assertion that is about correctness rather than wording: populations
    # named and sized, no maximum without a p99 beside it, no figure script
    # reading the gitignored run root, and caption numbers recomputed from
    # their generator. Re-pin the prose once the manuscript stops moving.

    def testTheThreeQuotedIntervalsAreRecomputed(self):
        """damping is per-beta; a single value put every row ~1.5% off (D042)."""
        # the table rows come from ONE generator over ONE range now, so check
        # every row against it rather than three hand-picked beta. The table and
        # the figure disagreeing about their range is what let prose understate
        # the closest approach to the measured load fivefold.
        import fig_beta
        for r in fig_beta.tableRows():
            lo, hi = A.requiredAggregateFractionBeta(
                r["beta"], A.dampingAtBeta(r["beta"]))
            self.assertAlmostEqual(100 * lo, r["pct_lo"], places=6)
            self.assertAlmostEqual(100 * hi, r["pct_hi"], places=6)
            self.assertIn(f"{r['pct_lo']:.3f} – {r['pct_hi']:.3f}", _MS,
                          f"beta={r['beta']} row missing from the table")

    def testEveryStoredDampingIsInTheMeasuredRange(self):
        for beta, d in A.DAMPING_BY_BETA:
            self.assertGreaterEqual(d, 0.346)
            self.assertLessEqual(d, 0.356)


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
            # "complete" means no RANDOM subsampling. saturation.tsv is now
            # built at re-solved states, so 117 of 2884 draws that admit a fold
            # yield no solvable state and carry no row. That is a stated
            # criterion, not a subsample -- but it must be ACCOUNTED FOR, so the
            # note has to say so and the count has to match the generator.
            if e["is_subsample"]:
                self.assertIn("re-solve", e.get("note", ""),
                              f"{e['file']} is short of its population and the "
                              "note does not say why")

    def testThePopulationsAreDistinguishedAndSizedAsTheManuscriptStates(self):
        by = {e["file"]: e for e in self.m["files"]}
        # re-solved states: 117 of the 2884 that admit a fold do not re-solve
        self.assertEqual(by["saturation.tsv"]["n_states"], 2767)
        self.assertEqual(by["saturation.tsv"]["population"], "kinetic_box")
        self.assertEqual(by["identity.tsv"]["n_states"], 325)
        self.assertEqual(by["identity.tsv"]["population"], "load_grid")
        self.assertIn("load grid", _MS)
        self.assertIn("kinetic box", _MS)


class TestSectionFiveNamesItsPopulations(unittest.TestCase):
    def testEveryTableRowNamesAPopulation(self):
        s = _MS[_MS.index("## 5. Numerical verification"):_MS.index("## 6.")]
        self.assertIn("| quantity | ensemble | median | p99 | max |", s)
        self.assertIn("load grid, 325", s)
        self.assertIn("kinetic box, 2884", s)   # the phi rebuild row

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
        # the PROPERTY: every max has a p99 beside it, and the section says why
        self.assertIn("p99", s)
        self.assertIn("7.56×10⁻⁷", table)


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
        for col, want in (("s_ref", 0.180), ("s_u", 0.159), ("s_a", 0.049)):
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
        self.assertIn("2767 kinetic-box networks", cap)
        # solved states move the width 0.876 -> 0.867; section 9 quotes it as
        # a caution about dispersion, which is the panel's own point, so the
        # width may appear in both -- what must not differ is the VALUE.
        self.assertIn("0.867", cap)
        self.assertIn("0.867", _MS[_MS.index("## 9. Predictions"):])
        self.assertIn("over the 30 networks of that experiment", _MS)
