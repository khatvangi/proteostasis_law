"""phase 3: caveats that must not be edited away.

D013's bistability numbers were computed under constant dilution, a regime that
by construction predicts zero growth-rate loss -- which is exactly the quantity
D015 and D026 measure as non-zero. the caveat therefore has to sit AT the claim,
ahead of the numbers, not in a trailing note a reader may never reach. these
tests assert both presence and position, following the K6/H3 pattern in
test_gate4_discriminating.py.
"""

import unittest
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]


def _sectionOf(doc, header, nextHeaderPrefix="## "):
    """return the text of one '## ...' section, exclusive of the next header."""
    start = doc.index(header)
    rest = doc[start + len(header):]
    cut = rest.find("\n" + nextHeaderPrefix)
    return header + (rest if cut < 0 else rest[:cut])


class TestD013CarriesItsCaveatInline(unittest.TestCase):
    def setUp(self):
        doc = (_REPO_ROOT / "DECISIONS.md").read_text()
        self.entry = _sectionOf(doc, "## D013 —")

    def testTheCaveatNamesConstantDilutionAndZeroGrowthLoss(self):
        low = self.entry.lower()
        self.assertIn("constant dilution", low)
        self.assertIn("k_mu = inf", low)
        self.assertRegex(self.entry, r"ZERO growth-rate loss")

    def testTheCaveatCrossReferencesD015AndD026(self):
        self.assertIn("D015", self.entry)
        self.assertIn("D026", self.entry)

    def testTheCaveatPrecedesTheNumbersItQualifies(self):
        """inline at the claim, not a trailing note -- position is the point."""
        caveat = self.entry.index("CAVEAT AT THE CLAIM")
        firstNumber = self.entry.index("0.194")
        self.assertLess(caveat, firstNumber)

    def testTheSharpFormIsStatedNotJustTheFactOfIdealisation(self):
        self.assertRegex(
            self.entry,
            r"regime that produces the bistability is the same regime\s+"
            r"that gets the measured quantity wrong",
        )


class TestStatusRestatesTheCaveat(unittest.TestCase):
    def setUp(self):
        self.doc = (_REPO_ROOT / "STATUS.md").read_text()

    def testStatusSaysThePointIsUnphysicalInAMeasuredWay(self):
        self.assertRegex(
            self.doc, r"That parameter point is unphysical, in a specific and\s+measured way"
        )
        self.assertIn("`k_mu = inf`", self.doc)

    def testStatusKeepsBothTheNotSurveyedAndTheUnphysicalQualifiers(self):
        """'not surveyed' was necessary but is no longer sufficient."""
        idx = self.doc.index("found at one parameter point,\nnot surveyed.")
        tail = self.doc[idx: idx + 1200]
        self.assertIn("exactly zero", tail)
        self.assertIn("D015", tail)
        self.assertIn("D026", tail)
        self.assertIn("Pinned by test so it cannot be edited away.", tail)


if __name__ == "__main__":
    unittest.main()
