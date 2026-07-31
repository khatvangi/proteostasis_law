"""the local free-limit transcription must BE nitrogen's model, not resemble it.

`scripts/phase2/nitrogen_limit.py` is a hand transcription of
nitrogen:.../src/proteostasis_model.py.  T0 and the whole free arm of the
benchmark are worthless if that transcription drifted, so it is checked against
a dump produced by EXECUTING the real module on nitrogen
(`data/phase2/nitrogen_reference.json`), not against a reading of it.

the reference was generated with the nitrogen source at sha256
9499a3a8f8642fa40a6360fb21039a6554c18411da7e82cf0d11f86f4df1c790, which is the
hash recorded in that workspace's own MANIFEST.md section 8.
"""

import json
import unittest
from pathlib import Path

import numpy as np

import _context  # noqa: F401
from phase2 import nitrogen_limit as nl
from phase2.mapping import NitrogenParams

REFERENCE = Path(__file__).resolve().parents[2] / "data" / "phase2" / "nitrogen_reference.json"
EXPECTED_SOURCE_SHA = "9499a3a8f8642fa40a6360fb21039a6554c18411da7e82cf0d11f86f4df1c790"


class TestTranscriptionFidelity(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.ref = json.loads(REFERENCE.read_text())

    def testReferenceWasProducedOnNitrogen(self):
        self.assertEqual(self.ref["environment"]["node"], "nitrogen")

    def testRecordedSourceHashIsStillPinned(self):
        pinned = (Path(__file__).resolve().parents[2]
                  / "scripts" / "phase2" / "NITROGEN_SOURCE.sha256").read_text()
        self.assertIn(EXPECTED_SOURCE_SHA, pinned)

    def testRightHandSideIsBitIdentical(self):
        """bit-identical, not merely close.

        the reference was computed under python 3.10 / numpy 1.26 and this runs
        under a different build, so exact equality is a real statement about the
        expression tree and its evaluation order, not a tautology.
        """
        checked = 0
        for case in self.ref["cases"]:
            q = NitrogenParams(**case["nitrogen_params"])
            for row in case["rows"]:
                got = nl.rhsVector([row["u"], row["a"]], q)
                want = np.array([float(v) for v in row["rhs"]])
                np.testing.assert_array_equal(got, want)
                checked += 1
        self.assertEqual(checked, 48)

    def testJacobianIsBitIdentical(self):
        for case in self.ref["cases"]:
            q = NitrogenParams(**case["nitrogen_params"])
            for row in case["rows"]:
                got = nl.jacobian(row["u"], row["a"], q).ravel()
                want = np.array([float(v) for v in row["jac"]])
                np.testing.assert_array_equal(got, want)

    def testResourceFractionsAreBitIdentical(self):
        for case in self.ref["cases"]:
            q = NitrogenParams(**case["nitrogen_params"])
            for row in case["rows"]:
                cf, df, _ = nl.resourceFractions(row["u"], row["a"], q)
                self.assertEqual(cf, float(row["cf"]))
                self.assertEqual(df, float(row["df"]))


class TestFreeLimitStructure(unittest.TestCase):
    """properties the epsilon = 0 face must have by construction."""

    def setUp(self):
        self.q = NitrogenParams(c_u=0.5, c_a=0.5, d_u=0.5, d_a=0.5, n_load=0.5)

    def testResourcePoolsReconstructTheirTotals(self):
        for u, a in ((0.0, 0.0), (0.3, 0.2), (5.0, 7.0)):
            res = nl.poolBalanceResiduals(u, a, self.q)
            for k, v in res.items():
                self.assertLess(abs(v), 1e-14, k)

    def testFreeCoordinatesAreTheStateItself(self):
        """there is no bound substrate in this convention; that IS the model."""
        for u, a in ((0.0, 0.0), (0.3, 0.2), (5.0, 7.0)):
            self.assertEqual(nl.freeCoordinates(u, a, self.q), (u, a))

    def testAnalyticJacobianMatchesComplexStep(self):
        """the check nitrogen's own suite runs, reproduced here."""
        rng = np.random.default_rng(11)
        for _ in range(20):
            u, a = rng.uniform(1e-3, 3.0, size=2)
            got = nl.jacobian(u, a, self.q)
            h = 1e-30
            cols = []
            for k in range(2):
                x = np.array([complex(u), complex(a)])
                x[k] += 1j * h
                cf = nl.rhs(x[0], x[1], self.q)
                cols.append(np.array([v.imag / h for v in cf]))
            want = np.column_stack(cols)
            np.testing.assert_allclose(got, want, rtol=2e-10, atol=1e-14)

    def testNascentOccupancyIsNotDamageInflux(self):
        """raising n_load must lower free chaperone and leave j untouched."""
        lo = NitrogenParams(n_load=0.1)
        hi = NitrogenParams(n_load=3.0)
        self.assertEqual(lo.j, hi.j)
        self.assertGreater(nl.resourceFractions(0.2, 0.1, lo)[0],
                           nl.resourceFractions(0.2, 0.1, hi)[0])


if __name__ == "__main__":
    unittest.main()
