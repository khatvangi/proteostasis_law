"""structural checks on the Gate 4 proposal.

the proposal is prose, but three of its properties are load-bearing and easy to
lose in editing: that it commits to reading no outcome, that its kill criteria
include VOID conditions rather than only reject conditions, and that the
substitution forced by the feasible instrument (s_u in place of s_a) is stated
rather than glossed. these assert those, nothing more.
"""

import re
import unittest
from pathlib import Path

DOC = Path(__file__).resolve().parents[2] / "empirical" / "GATE4_PROPOSAL.md"


class TestGate4Proposal(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.t = DOC.read_text()

    def testExists(self):
        self.assertTrue(DOC.is_file())

    def testDeclaresNoOutcomeWasRead(self):
        self.assertIn("no outcome value", self.t.lower())

    def testHasKillCriteriaIncludingVoidConditions(self):
        for k in ("K1.", "K2.", "K3."):
            self.assertIn(k, self.t)
        self.assertGreaterEqual(self.t.count("VOID"), 2,
                                "kill criteria must include VOID, not only reject")

    def testStatesTheArmSubstitutionExplicitly(self):
        """s_u is measurable, s_a is not; quoting s_a while measuring s_u is the
        failure mode this asserts against."""
        self.assertIn("s_u", self.t)
        self.assertIn("s_a", self.t)
        self.assertIn("H1'", self.t)

    def testRecordsBothFeasibilityCitations(self):
        for doi in ("10.1021/acssynbio.4c00612", "10.1088/1478-3975/13/2/025002"):
            self.assertIn(doi, self.t)

    def testKeepsTheChaperoneArmMarkedUntested(self):
        self.assertRegex(self.t, r"chaperone arm remains unmeasurable|folding side")
