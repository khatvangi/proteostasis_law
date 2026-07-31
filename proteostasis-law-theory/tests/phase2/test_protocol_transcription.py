"""both root protocols must BE the originals, not approximations of them.

the 2x2 factorial attributes cell-to-cell differences to exactly two causes:
model form and root protocol.  that attribution is only sound if each protocol
in `scripts/phase2/protocols.py` is the thing it claims to be:

  P_BORON    == `proteostasis.equilibria.findEquilibria` (checked against the
                shipped function directly, on the shipped model);
  P_NITROGEN == nitrogen's `lhs_sweep.classify_one` (checked against a dump
                produced by EXECUTING that function on nitrogen).
"""

import json
import unittest
from pathlib import Path

import numpy as np

import _context  # noqa: F401
from phase2 import protocols
from phase2.lhs import parametersForEpsilon, sampleMatrix
from phase2.mapping import nitrogenToBoron
from phase2.models import BoronAdapter, FreeLimitAdapter
from proteostasis import equilibria as boron_equilibria
from proteostasis.model import Params

REFERENCE = (Path(__file__).resolve().parents[2] / "data" / "phase2"
             / "nitrogen_protocol_reference.json")
MATRIX = sampleMatrix(60, 20260731)


class TestBoronProtocolTranscription(unittest.TestCase):
    """the generic multistart must reproduce the shipped one exactly.

    the generic version exists because `equilibria.findEquilibria` is hard-wired
    to `model.Params` and cannot solve the free-limit arm.  if the two ever
    diverged, cell 3 vs its free counterpart would measure a transcription bug
    and call it a model-form effect.
    """

    def testReproducesFindEquilibria(self):
        for eps in (1e-6, 0.3, 2.0):
            for i in (0, 4, 17, 33, 51):
                q = parametersForEpsilon(MATRIX[i], eps)
                adapter = BoronAdapter(q, eps)
                p = Params(**nitrogenToBoron(q, eps)).validate()

                mine = protocols.findEquilibriaBoron(adapter)
                theirs = boron_equilibria.findEquilibria(p)

                self.assertEqual(len(mine), len(theirs),
                                 f"root count differs at eps={eps}, sample {i}")
                for a, b in zip(mine, theirs):
                    self.assertEqual(a["u"], b.u)
                    self.assertEqual(a["a"], b.a)
                    self.assertEqual(a["stable"], b.stable)
                    self.assertEqual(a["residual"], b.residual)
                    self.assertEqual(a["eig_real_max"], b.eig_real_max)

    def testConstantsMatchTheShippedOnes(self):
        self.assertEqual(protocols.BORON_GRID, 9)
        self.assertEqual((protocols.BORON_LO, protocols.BORON_HI), (1e-7, 1e4))
        self.assertEqual(protocols.BORON_RES_TOL, 1e-8)
        self.assertEqual(protocols.BORON_DEDUPE_RTOL, 1e-5)

    def testBoronProtocolNeverLeavesTheOrthant(self):
        """log coordinates make positivity structural; the continuation must
        never fire under this protocol."""
        for eps in (1e-6, 2.0):
            for i in (2, 22):
                adapter = BoronAdapter(parametersForEpsilon(MATRIX[i], eps), eps)
                protocols.classifyBoron(adapter, t_end=5.0e3)
                self.assertEqual(adapter.continuation_evaluations, 0)


class TestNitrogenProtocolTranscription(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ref = json.loads(REFERENCE.read_text())

    def testReferenceWasProducedOnNitrogen(self):
        self.assertEqual(self.ref["environment"]["node"], "nitrogen")
        self.assertEqual(self.ref["seed"], 20260731)

    def testLabelsAgreeWithRealNitrogenClassifyOne(self):
        """nitrogen's four labels collapse onto ours: both threshold labels are
        `stable_attractor`, since the threshold split is applied downstream."""
        collapse = {"stable_subthreshold": protocols.LABEL_STABLE,
                    "stable_overthreshold": protocols.LABEL_STABLE,
                    "no_bounded_attractor_operational": protocols.LABEL_NONE,
                    "numerical_failure": protocols.LABEL_FAIL}
        checked = 0
        for case in self.ref["cases"]:
            eps = case["epsilon"]
            for row in case["rows"]:
                i = row["sample"]
                adapter = FreeLimitAdapter(parametersForEpsilon(MATRIX[i], eps), eps)
                got = protocols.classifyNitrogen(adapter)
                self.assertEqual(got["label"], collapse[row["label"]],
                                 f"eps={eps} sample={i}")
                checked += 1
        self.assertEqual(checked, 240)

    def testRootsAgreeWithRealNitrogenToMachinePrecision(self):
        for case in self.ref["cases"]:
            eps = case["epsilon"]
            for row in case["rows"]:
                if row["u_eq"] is None:
                    continue
                i = row["sample"]
                adapter = FreeLimitAdapter(parametersForEpsilon(MATRIX[i], eps), eps)
                got = protocols.classifyNitrogen(adapter)
                self.assertIsNotNone(got["root"])
                self.assertAlmostEqual(got["root"]["u"], float(row["u_eq"]), places=12)
                self.assertAlmostEqual(got["root"]["a"], float(row["a_eq"]), places=12)

    def testThresholdSplitReproducesNitrogensOwn(self):
        """nitrogen's componentwise D on free coordinates == our adm_comp_free."""
        for case in self.ref["cases"]:
            eps = case["epsilon"]
            for row in case["rows"]:
                if row["label"] not in ("stable_subthreshold", "stable_overthreshold"):
                    continue
                i = row["sample"]
                adapter = FreeLimitAdapter(parametersForEpsilon(MATRIX[i], eps), eps)
                got = protocols.classifyNitrogen(adapter)
                uf, af = got["root"]["u_free"], got["root"]["a_free"]
                ours = bool(uf < 1.0 and af < 1.0 and uf >= 0 and af >= 0)
                self.assertEqual(ours, row["label"] == "stable_subthreshold",
                                 f"eps={eps} sample={i}")

    def testConstantsMatchTheShippedOnes(self):
        self.assertEqual(protocols.NITROGEN_GUESSES,
                         ((0.01, 0.01), (0.1, 0.05), (0.5, 0.2), (1.0, 1.0)))
        self.assertEqual(protocols.NITROGEN_RES_TOL, 1e-7)
        self.assertEqual(protocols.NITROGEN_RETURN_TOL, 1e-3)
        self.assertEqual(protocols.NITROGEN_T_END, 150.0)


class TestProtocolsAreModelBlind(unittest.TestCase):

    def testBothProtocolsAcceptBothModelForms(self):
        for eps in (1e-6, 2.0):
            q = parametersForEpsilon(MATRIX[1], eps)
            for adapter in (BoronAdapter(q, eps), FreeLimitAdapter(q, eps)):
                for fn in (protocols.classifyNitrogen,
                           lambda a: protocols.classifyBoron(a, t_end=5.0e3)):
                    res = fn(adapter)
                    self.assertIn(res["label"], (protocols.LABEL_STABLE,
                                                 protocols.LABEL_NONE,
                                                 protocols.LABEL_FAIL))


if __name__ == "__main__":
    unittest.main()
