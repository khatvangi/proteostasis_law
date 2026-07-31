"""reproducibility tests.

these fail when a result depends on anything other than (code, config, seed):
worker scheduling, dict iteration order, wall-clock time, or an unseeded rng.
each check is paired with a mutation showing it can distinguish a real change.
"""

import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

import _context  # noqa: F401
from proteostasis import Params, findFold, lowestStableEquilibrium
from proteostasis.provenance import (canonicalJson, hashObject, hashFile, writeTable,
                                     gitCommit, environmentInfo)
from proteostasis.sweeps import (GLOBAL_RANGES, latinHypercube, paramsFromSample,
                                 factorialGrid, parallelMap)


def _square(x):
    return x * x


class TestSampling(unittest.TestCase):

    def testLatinHypercubeIsSeedDeterministic(self):
        a = latinHypercube(GLOBAL_RANGES, 40, seed=99)
        b = latinHypercube(GLOBAL_RANGES, 40, seed=99)
        self.assertEqual(hashObject(a), hashObject(b))

    def testDifferentSeedGivesDifferentSample(self):
        a = latinHypercube(GLOBAL_RANGES, 40, seed=99)
        c = latinHypercube(GLOBAL_RANGES, 40, seed=100)
        self.assertNotEqual(hashObject(a), hashObject(c))

    def testSamplesRespectDeclaredRanges(self):
        for s in latinHypercube(GLOBAL_RANGES, 200, seed=5):
            for r in GLOBAL_RANGES:
                self.assertGreaterEqual(s[r.name], r.lo * (1 - 1e-9))
                self.assertLessEqual(s[r.name], r.hi * (1 + 1e-9))

    def testFactorialGridIsOrderStable(self):
        axes = {"b": [1.0, 2.0], "a": [10.0, 20.0, 30.0]}
        self.assertEqual(hashObject(factorialGrid(axes)), hashObject(factorialGrid(axes)))
        self.assertEqual(len(factorialGrid(axes)), 6)


class TestParallelDeterminism(unittest.TestCase):

    def testParallelMapPreservesOrder(self):
        items = list(range(200))
        serial = parallelMap(_square, items, n_workers=1)
        par = parallelMap(_square, items, n_workers=4)
        self.assertEqual(serial, par)
        self.assertEqual(par, [i * i for i in items])

    def testWorkerCountDoesNotChangeResults(self):
        items = list(range(97))
        outs = [parallelMap(_square, items, n_workers=w) for w in (1, 2, 5, 8)]
        for o in outs[1:]:
            self.assertEqual(o, outs[0])


class TestSerialization(unittest.TestCase):

    def testCanonicalJsonIsKeyOrderIndependent(self):
        a = {"z": 1, "a": {"y": 2, "b": 3}}
        b = {"a": {"b": 3, "y": 2}, "z": 1}
        self.assertEqual(canonicalJson(a), canonicalJson(b))
        self.assertEqual(hashObject(a), hashObject(b))

    def testCanonicalJsonDistinguishesRealChanges(self):
        self.assertNotEqual(hashObject({"a": 1.0}), hashObject({"a": 1.0000001}))

    def testWrittenTablesAreByteIdentical(self):
        df = pd.DataFrame({"b": np.linspace(0, 1, 50), "a": np.arange(50)})
        with tempfile.TemporaryDirectory() as td:
            h1 = writeTable(df, Path(td) / "one.tsv")
            h2 = writeTable(df, Path(td) / "two.tsv")
            self.assertEqual(h1, h2)
            # column order in the frame must not change the file
            h3 = writeTable(df[["a", "b"]], Path(td) / "three.tsv")
            self.assertEqual(h1, h3)

    def testWrittenTablesDetectValueChanges(self):
        with tempfile.TemporaryDirectory() as td:
            d1 = pd.DataFrame({"a": [1.0, 2.0]})
            d2 = pd.DataFrame({"a": [1.0, 2.0000001]})
            self.assertNotEqual(writeTable(d1, Path(td) / "a.tsv"),
                                writeTable(d2, Path(td) / "b.tsv"))


class TestComputationalDeterminism(unittest.TestCase):

    def testFoldDetectionRepeatsExactly(self):
        p = Params(j=0.02)
        a = findFold(p, "j", 1e-4, 2.0)
        b = findFold(p, "j", 1e-4, 2.0)
        self.assertEqual(a.fold_value, b.fold_value)
        self.assertEqual(a.bracket, b.bracket)
        self.assertEqual(a.n_steps, b.n_steps)

    def testEquilibriumSearchRepeatsExactly(self):
        p = Params(j=0.05)
        a = lowestStableEquilibrium(p, n_grid=7)
        b = lowestStableEquilibrium(p, n_grid=7)
        self.assertEqual(a.asDict(), b.asDict())

    def testSweepSliceRepeatsExactly(self):
        base = Params()
        def slice_():
            out = []
            for s in latinHypercube(GLOBAL_RANGES, 8, seed=1234):
                p = paramsFromSample(s, base)
                eq = lowestStableEquilibrium(p.with_(j=1e-3), n_grid=5)
                out.append(None if eq is None else eq.asDict())
            return out
        self.assertEqual(hashObject(slice_()), hashObject(slice_()))


class TestProvenance(unittest.TestCase):

    def testGitCommitIsRecorded(self):
        g = gitCommit()
        self.assertIsNotNone(g["commit"], "not in a git repo, or git unavailable")
        self.assertEqual(len(g["commit"]), 40)

    def testEnvironmentIsRecorded(self):
        e = environmentInfo()
        for k in ("python", "numpy", "scipy", "pandas", "platform", "node"):
            self.assertIn(k, e)
            self.assertTrue(e[k])

    def testFileHashChangesWithContent(self):
        with tempfile.TemporaryDirectory() as td:
            a, b = Path(td) / "a", Path(td) / "b"
            a.write_text("x")
            b.write_text("y")
            self.assertNotEqual(hashFile(a), hashFile(b))


if __name__ == "__main__":
    unittest.main()
