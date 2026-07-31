"""the benchmark design: determinism, matching, and the four named cells.

these tests do not run the full benchmark.  they assert the properties that make
its output interpretable -- that both hosts draw the same sample, that every cell
sees the same draw, that admissibility is reported under all four criteria rather
than one, and that cell 1 versus cell 2 really is a code-equivalence test.
"""

import unittest

import numpy as np

import _context  # noqa: F401
from phase2 import protocols, run_matched_benchmark as bench
from phase2.lhs import (DEFAULT_N, DEFAULT_SEED, SAMPLE_HASHES, describeDesign,
                        parametersForEpsilon, sampleHash, sampleMatrix,
                        sampledCoordinates)
from phase2.mapping import EPSILON_LADDER, LHS_COORDINATES, PINNED_NITROGEN
from phase2.models import BoronAdapter, FreeLimitAdapter


class TestSamplingDeterminism(unittest.TestCase):

    def testPinnedHashesReproduce(self):
        """these hashes were verified bit-identical on boron (scipy 1.14.0) and
        nitrogen (scipy 1.15.2).  if this fails, the two hosts are no longer
        scoring the same parameter draws and no cross-host cell may be compared.
        """
        for n, want in SAMPLE_HASHES.items():
            self.assertEqual(sampleHash(sampleMatrix(n, DEFAULT_SEED)), want, f"n={n}")

    def testRepeatedCallsAreBitIdentical(self):
        a = sampleMatrix(64, DEFAULT_SEED)
        b = sampleMatrix(64, DEFAULT_SEED)
        np.testing.assert_array_equal(a, b)

    def testDesignRecordFlagsTheHashMatch(self):
        d = describeDesign(DEFAULT_N, DEFAULT_SEED)
        self.assertTrue(d["sample_matrix_hash_pinned_match"])
        self.assertEqual(d["n_samples"], 2000)
        self.assertEqual(d["seed"], 20260731)


class TestCoordinateTransforms(unittest.TestCase):

    def testRangesAreRespectedAtTheCorners(self):
        for name, kind, lo, hi in LHS_COORDINATES:
            k = [c[0] for c in LHS_COORDINATES].index(name)
            z0 = np.zeros(7)
            z1 = np.ones(7)
            self.assertAlmostEqual(sampledCoordinates(z0)[name], lo, places=12)
            self.assertAlmostEqual(sampledCoordinates(z1)[name], hi, places=12)
            self.assertIn(kind, ("log", "linear"))
            del k

    def testSevenCoordinatesExactly(self):
        self.assertEqual(len(LHS_COORDINATES), 7)

    def testPinnedCoordinatesAreActuallyPinnedAcrossSamples(self):
        M = sampleMatrix(32, DEFAULT_SEED)
        sampled = {c[0] for c in LHS_COORDINATES}
        first = parametersForEpsilon(M[0], 0.3).asDict()
        for i in range(1, 32):
            q = parametersForEpsilon(M[i], 0.3).asDict()
            for k, v in first.items():
                if k in sampled or k in ("c_u", "c_a"):
                    continue          # c_u, c_a depend on ref_capacity by design
                self.assertEqual(q[k], v, f"{k} moved between samples")

    def testPinnedValuesMatchNitrogensDefaults(self):
        q = parametersForEpsilon(sampleMatrix(4, DEFAULT_SEED)[0], 1.0).asDict()
        for k, v in PINNED_NITROGEN.items():
            self.assertEqual(q[k], v, k)


class TestFactorialStructure(unittest.TestCase):

    def testCellNamesAreUniqueAcrossTheWholeFactorial(self):
        names = [bench.cellName(f, e, p)
                 for f in ("boron", "free")
                 for e in EPSILON_LADDER
                 for p in ("boron", "nitrogen")]
        self.assertEqual(len(names), len(set(names)))
        self.assertEqual(len(names), 28)

    def testTheFourAuditCellsExist(self):
        self.assertEqual(bench.cellName("boron", 1e-6, "nitrogen"),
                         "boron_eps1e-06_nitrogen")
        self.assertEqual(bench.cellName("free", 1e-6, "nitrogen"),
                         "free_eps1e-06_nitrogen")
        self.assertEqual(bench.cellName("boron", 1e-6, "boron"),
                         "boron_eps1e-06_boron")
        self.assertEqual(bench.cellName("boron", 2.0, "boron"), "boron_eps2_boron")

    def testEpsilonLadderIsTheSpecifiedOne(self):
        self.assertEqual(EPSILON_LADDER, (1e-6, 1e-3, 1e-2, 1e-1, 0.3, 1.0, 2.0))


class TestAdmissibilityIsReportedNotMixed(unittest.TestCase):

    def testAllFourCriteriaArePresentInEveryRow(self):
        for k in ("adm_burden_total", "adm_comp_total",
                  "adm_burden_free", "adm_comp_free"):
            self.assertIn(k, bench.ROW_FIELDS)

    def testBurdenAndComponentwiseAreGenuinelyDifferentSets(self):
        """(0.7, 0.7) is componentwise admissible and burden inadmissible."""
        adm = bench._admissibility(0.7, 0.7, 0.7, 0.7)
        self.assertTrue(adm["adm_comp_total"])
        self.assertFalse(adm["adm_burden_total"])

    def testThresholdIsSharedByBothConfigs(self):
        """boron's burden_threshold and nitrogen's threshold_u/a are both 1.0."""
        self.assertEqual(bench.THRESHOLD_H, 1.0)
        self.assertEqual(PINNED_NITROGEN["threshold_u"], 1.0)
        self.assertEqual(PINNED_NITROGEN["threshold_a"], 1.0)

    def testLabelsDoNotEncodeAThreshold(self):
        for label in (protocols.LABEL_STABLE, protocols.LABEL_NONE, protocols.LABEL_FAIL):
            self.assertNotIn("threshold", label)


class TestCellOneEqualsCellTwo(unittest.TestCase):
    """the code-equivalence test, on a small deterministic subsample.

    the full n = 2000 version is the benchmark itself; this is the guard that
    keeps the claim true as the code changes.
    """

    def testLabelsAndRootsAgreeAtTheAnchor(self):
        M = sampleMatrix(40, DEFAULT_SEED)
        eps = 1e-6
        for i in range(40):
            q = parametersForEpsilon(M[i], eps)
            b = protocols.classifyNitrogen(BoronAdapter(q, eps))
            f = protocols.classifyNitrogen(FreeLimitAdapter(q, eps))
            self.assertEqual(b["label"], f["label"], f"sample {i}")
            if b["root"] is None:
                self.assertIsNone(f["root"])
                continue
            for k in ("u", "a"):
                rel = abs(b["root"][k] - f["root"][k]) / max(abs(f["root"][k]), 1e-12)
                self.assertLess(rel, 1e-4, f"sample {i} coordinate {k}")

    def testAdmissibilityAgreesAtTheAnchor(self):
        M = sampleMatrix(40, DEFAULT_SEED)
        eps = 1e-6
        for i in range(40):
            q = parametersForEpsilon(M[i], eps)
            b = protocols.classifyNitrogen(BoronAdapter(q, eps))
            f = protocols.classifyNitrogen(FreeLimitAdapter(q, eps))
            if b["root"] is None or f["root"] is None:
                continue
            ab = bench._admissibility(b["root"]["u"], b["root"]["a"],
                                      b["root"]["u_free"], b["root"]["a_free"])
            af = bench._admissibility(f["root"]["u"], f["root"]["a"],
                                      f["root"]["u_free"], f["root"]["a_free"])
            self.assertEqual(ab, af, f"sample {i}")


class TestFreeArmNeedsNoBoronCode(unittest.TestCase):
    """the nitrogen host runs the free arm only and must not need `proteostasis`.

    asserted structurally: the boron import happens inside `BoronAdapter`, so
    building and running a free adapter must not touch it.
    """

    def testFreeAdapterDoesNotImportProteostasis(self):
        import subprocess
        import sys
        from pathlib import Path

        scripts = Path(__file__).resolve().parents[2] / "scripts"
        code = (
            "import sys; sys.path.insert(0, %r)\n"
            "import phase2.models as m, phase2.protocols as pr\n"
            "from phase2.lhs import sampleMatrix, parametersForEpsilon\n"
            "q = parametersForEpsilon(sampleMatrix(4, 20260731)[0], 1e-6)\n"
            "pr.classifyNitrogen(m.FreeLimitAdapter(q, 1e-6))\n"
            "assert 'proteostasis' not in sys.modules, sorted(k for k in sys.modules "
            "if 'proteo' in k)\n"
            "print('clean')\n" % str(scripts)
        )
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
        self.assertEqual(out.returncode, 0, out.stderr)
        self.assertIn("clean", out.stdout)


if __name__ == "__main__":
    unittest.main()
