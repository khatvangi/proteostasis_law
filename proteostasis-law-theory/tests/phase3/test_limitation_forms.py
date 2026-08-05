"""phase 3: no limitation may be phrased as a licence (D040).

a limitation that says "these should be screened" is a free parameter with prior
authorisation sitting in the repository waiting to be used. one did, for months,
and figure 2 acted on it and produced a 5x figure-text divergence (D039).

these tests pin the PROPERTY and not the token, per the rule in
notes/VERIFICATION_RULE.md -- the banned phrase is allowed to appear where the
text is recording its own history, and is banned only as a current instruction.
"""

import re
import unittest
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]

# every document that states limitations or caveats a later session could act on
_LIMITATION_DOCS = (
    "theory/FOLD_THEOREM.md",
    "STATUS.md",
    "notes/VERIFICATION_RULE.md",
)

# markers that turn a directive into a record of a superseded directive
_SUPERSESSION = ("earlier revision", "an earlier", "read \"", "> ", "D039", "D040")


class TestScreeningIsNotLicensed(unittest.TestCase):
    def setUp(self):
        self.fold = (_REPO_ROOT / "theory/FOLD_THEOREM.md").read_text()

    def testTheLimitsBulletSeparatesTheObservationFromTheExclusion(self):
        """the draws being marginal and the draws being removable are two claims."""
        bullet = self.fold[self.fold.index("Some draws collapse at `s_a` near 0.003"):][:1400]
        low = bullet.lower()
        self.assertIn("not a licence", low)
        self.assertIn("d039", low)
        # the tested reason a screen fails has to sit here, not only in the log
        self.assertIn("no gap", low)
        self.assertRegex(bullet, r"factor of four")

    def testTheDirectiveSurvivesOnlyAsItsOwnHistory(self):
        """pin the property: banned as a current instruction, allowed as a record."""
        for hit in re.finditer(r"should be screened", self.fold):
            before = self.fold[max(0, hit.start() - 300): hit.start()].lower()
            self.assertTrue(
                any(m.lower() in before for m in _SUPERSESSION),
                msg="'should be screened' appears without a supersession marker",
            )

    def testNoOtherDocumentReinstatesAScreen(self):
        """a directive form is permitted where it is QUOTED, never where asserted.

        the rule document has to name the banned form in order to ban it, and the
        FOLD_THEOREM bullet has to name it to record that it was withdrawn. the
        distinguishing property is quotation, not the absence of the string --
        an absence test here would fail on exactly the text that fixes the defect.
        """
        banned = re.compile(
            r"(should|must|need to) be (screened|excluded|dropped|discarded|filtered)"
        )
        for rel in _LIMITATION_DOCS:
            doc = (_REPO_ROOT / rel).read_text()
            for hit in banned.finditer(doc):
                line = doc[doc.rfind("\n", 0, hit.start()) + 1:
                           doc.find("\n", hit.end())]
                before = doc[max(0, hit.start() - 300): hit.start()].lower()
                quoted = line.lstrip().startswith(">") or '"' in line
                self.assertTrue(
                    quoted or any(m.lower() in before for m in _SUPERSESSION),
                    msg=f"{rel}: directive-form limitation asserted, not quoted:\n{line}",
                )


class TestExtremaCarryTheirDirectionEverywhere(unittest.TestCase):
    """D038 relabelled these; the relabelling had reached the manuscript only."""

    def testThirteenPointSixIsNeverStatedAsABound(self):
        for rel in ("STATUS.md", "theory/FOLD_THEOREM.md", "manuscript/MANUSCRIPT_BMB_v5.md"):
            doc = (_REPO_ROOT / rel).read_text()
            self.assertNotRegex(
                doc, r"(up to|as much as|at most)\s+\*{0,2}13\.6",
                msg=f"{rel}: 13.6x stated as a bound",
            )

    def testEverySiteQuotingThirteenPointSixNamesItsSample(self):
        for rel in ("STATUS.md", "theory/FOLD_THEOREM.md", "manuscript/MANUSCRIPT_BMB_v5.md"):
            doc = (_REPO_ROOT / rel).read_text()
            if "13.6" not in doc:
                continue
            low = doc.lower()
            self.assertTrue(
                "largest-observed" in low or "largest observed" in low
                or "observed" in low,
                msg=f"{rel}: 13.6x quoted without saying it is an observed value",
            )
            self.assertRegex(doc, r"10 (kinetic )?draws|ten draws")

    def testTheUnboundedPositionIsStatedNotImplied(self):
        """a referee pushing on 13.6x should meet a written answer, not a new one."""
        for rel in ("STATUS.md", "theory/FOLD_THEOREM.md"):
            low = (_REPO_ROOT / rel).read_text().lower()
            self.assertIn("not sized to bound", low)

    def testSixTwoThreeIsALowerBoundAtEverySite(self):
        for rel in ("STATUS.md", "manuscript/MANUSCRIPT_BMB_v5.md"):
            doc = (_REPO_ROOT / rel).read_text()
            if "0.623" not in doc:
                continue
            window = doc[max(0, doc.index("0.623") - 300): doc.index("0.623") + 400]
            self.assertRegex(window.lower(), r"lower bound|largest observed")
            self.assertNotRegex(window, r"max 0\.623")


class TestTheSubsampleDefaultIsDisarmed(unittest.TestCase):
    """correcting the numbers while leaving the mechanism armed fixes a symptom."""

    def testVerifyAgainstRunDefaultsToTheCompletePopulation(self):
        import inspect
        import sys
        sys.path.insert(0, str(_REPO_ROOT / "scripts"))
        sys.path.insert(0, str(_REPO_ROOT / "scripts" / "phase3"))
        import fold_theorem as FT
        sig = inspect.signature(FT.verifyAgainstRun)
        self.assertIsNone(sig.parameters["n_identity"].default,
                          "the 20-state default that caused D036 is back")

    def testTheRetractedMetricIsLabelledWhereItIsPrinted(self):
        src = (_REPO_ROOT / "scripts/phase3/fold_theorem.py").read_text()
        self.assertIn("RETRACTED metric", src)
        self.assertIn("identity_grad_err_median", src)


class TestBothPopulationQuestionsAreAsked(unittest.TestCase):
    """the sixth affected number was a complete population under the wrong sentence."""

    def setUp(self):
        self.rule = (_REPO_ROOT / "notes/VERIFICATION_RULE.md").read_text()

    def testCompletenessAndAttributionAreSeparateChecks(self):
        low = self.rule.lower()
        self.assertIn("is that the whole of it", low)
        self.assertIn("is it the population the surrounding sentence is about", low)

    def testTheAttributionCheckSaysWhyCompletenessCannotCatchIt(self):
        idx = self.rule.lower().index("surrounding sentence is about")
        tail = self.rule[idx: idx + 900].lower()
        self.assertIn("nothing is missing", tail)
        self.assertIn("0.890", tail)
        self.assertIn("0.876", tail)

    def testTheLicenceRuleIsRecordedAsAnAuditForm(self):
        low = self.rule.lower()
        self.assertIn("a limitation is an observation, not a licence", low)
        self.assertIn("does not license a threshold", low)
        # the three-way classification is the operable part
        for form in ("an observation", "a restriction", "a direction"):
            self.assertIn(form, low)


if __name__ == "__main__":
    unittest.main()
