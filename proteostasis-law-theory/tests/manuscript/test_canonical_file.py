"""The READMEs must name the file the build actually reads.

Task B8.7. Both READMEs declared `bmb_v4.md` canonical while the converter read
`MANUSCRIPT_BMB_v5.md` and every current number came from the latter. That is a
lineage split sitting in the documentation, and it is the THIRD recurrence of the
same failure -- which is the argument for asserting it rather than fixing it
again.

The assertion is deliberately against `to_latex.SRC` rather than against a
hardcoded name. A future v6 renames one constant and this test follows it; a
future v6 that renames the constant and forgets the README fails here.
"""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

_REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO_ROOT / "scripts" / "manuscript"))

import to_latex as T  # noqa: E402


class TestCanonicalFileIsNamedConsistently(unittest.TestCase):

    def setUp(self):
        self.canonical = T.SRC.name
        self.stem = T.OUT_TEX.stem

    def testTheConverterReadsAFileThatExists(self):
        self.assertTrue(T.SRC.exists(), f"{self.canonical} is missing")

    def testBothReadmesNameTheCanonicalManuscript(self):
        for rel in ("README.md", "manuscript/README.md"):
            text = (_REPO_ROOT / rel).read_text()
            self.assertIn(self.canonical, text,
                          f"{rel} does not name {self.canonical}")

    def testNeitherReadmeStillDeclaresASupersededFileCanonical(self):
        """the specific shape of the bug: 'X is canonical' for the wrong X."""
        import re
        for rel in ("README.md", "manuscript/README.md"):
            text = (_REPO_ROOT / rel).read_text()
            for m in re.finditer(r"[Cc]anonical:?\**\s*\[?`([^`]+)`", text):
                self.assertEqual(m.group(1), self.canonical,
                                 f"{rel} declares {m.group(1)} canonical")
            for m in re.finditer(r"`([^`]+)` is canonical", text):
                self.assertEqual(m.group(1), self.canonical,
                                 f"{rel} declares {m.group(1)} canonical")

    def testTheRootReadmeNamesTheBuildProducts(self):
        text = (_REPO_ROOT / "README.md").read_text()
        self.assertIn(f"manuscript/{self.stem}.", text,
                      f"README does not name the build output {self.stem}")

    def testTheZenodoMetadataMatchesTheCitationVersion(self):
        """.zenodo.json and CITATION.cff must not disagree about the version."""
        import json
        z = json.loads((_REPO_ROOT / ".zenodo.json").read_text())
        cff = (_REPO_ROOT / "CITATION.cff").read_text()
        self.assertIn(f"version: {z['version']}", cff,
                      "CITATION.cff and .zenodo.json disagree on the version")
        # the title and the description must describe the same claim
        self.assertIn("Fold Condition", z["title"])
        self.assertNotIn("collapse boundary", z["description"].lower())


if __name__ == "__main__":
    unittest.main()
