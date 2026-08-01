"""scientific-equivalence and checkpoint-integrity tests for the D recovery runner.

the recovery runner exists to change WHEN and WHERE experiment D's work happens,
never WHAT it computes. these tests are the gate on that claim:

  * the LHS matrix it regenerates is the original's, element for element;
  * a background computed through the checkpoint path is numerically identical
    to calling `run_experiment_d._backgroundTask` directly -- every float, not
    just every summary statistic;
  * a valid checkpoint is reused, never silently recomputed;
  * a checkpoint that is corrupted, truncated, or belongs to a different config,
    source version, sample or index is REJECTED rather than trusted;
  * a timed-out background is classified `unresolved_timeout`, contributes no
    rows, and is excluded from every interaction summary;
  * the merged summary is what `run_experiment_d.main` would have written for
    the same set of backgrounds.

the two "fast usable" indices below were chosen by measurement, not assumption:
they are the quickest draws that still exercise the full 108-row factorial.
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

import pandas as pd

import _context  # noqa: F401
import run_experiment_d as original
import run_experiment_d_checkpointed as ckpt
from proteostasis.model import Params
from proteostasis.provenance import canonicalJson, hashObject
from proteostasis.sweeps import latinHypercube, rangesFromConfig

REPO_ROOT = Path(__file__).resolve().parents[1]
CONFIG = REPO_ROOT / "configs" / "phase1" / "experiment_d.json"
RUNNER = REPO_ROOT / "scripts" / "run_experiment_d_checkpointed.py"

#: measured ~17 s each and produce the full 3 pairs x 6 x 6 factorial
FAST_USABLE = (48, 56)
#: measured ~0.25 s each; exercise the early-return (unusable) path
FAST_UNUSABLE = (5, 35)
FAST = tuple(sorted(FAST_USABLE + FAST_UNUSABLE))


def _loadCfg():
    with open(CONFIG) as fh:
        return json.load(fh)


class TestSampleMatrix(unittest.TestCase):
    """the recovery run must evaluate the same parameter draws as the original."""

    def testSampleMatrixIsIdenticalToOriginalConstruction(self):
        cfg = _loadCfg()
        # rebuilt exactly as run_experiment_d.main builds it
        base_o = Params(**cfg.get("base_params", {})).validate()
        ranges_o = rangesFromConfig(cfg.get("param_ranges"))
        samples_o = latinHypercube(ranges_o, cfg["n_backgrounds"], cfg["seed"])

        base_c, samples_c = ckpt.sampleMatrix(cfg)
        self.assertEqual(base_o, base_c)
        self.assertEqual(len(samples_o), 60)
        self.assertEqual(samples_o, samples_c)
        self.assertEqual(hashObject(samples_o), hashObject(samples_c))

    def testEveryBackgroundSampleMatchesElementwise(self):
        cfg = _loadCfg()
        ranges = rangesFromConfig(cfg.get("param_ranges"))
        samples_o = latinHypercube(ranges, cfg["n_backgrounds"], cfg["seed"])
        _b, samples_c = ckpt.sampleMatrix(cfg)
        for i, (a, b) in enumerate(zip(samples_o, samples_c)):
            self.assertEqual(sorted(a), sorted(b), f"parameter names differ at {i}")
            for k in a:
                # exact float equality, not assertAlmostEqual: a resampled draw
                # would differ in the low bits and must fail here
                self.assertEqual(a[k], b[k], f"background {i} parameter {k}")

    def testDifferentNBackgroundsGivesDifferentMatrix(self):
        """guards the reason --n-backgrounds must be forwarded to children."""
        cfg = _loadCfg()
        _b, full = ckpt.sampleMatrix(cfg)
        _b, small = ckpt.sampleMatrix(dict(cfg, n_backgrounds=3))
        self.assertNotEqual(full[0], small[0])


class TestCheckpointedEquivalence(unittest.TestCase):
    """the expensive gate: checkpoint path vs. direct call, float for float."""

    @classmethod
    def setUpClass(cls):
        cls.cfg = _loadCfg()
        cls.base, cls.samples = ckpt.sampleMatrix(cls.cfg)
        cls.identity = ckpt.runIdentity(cls.cfg, cls.samples)
        cls.tmp = Path(tempfile.mkdtemp(prefix="d_ckpt_test_"))
        cls.ckdir = cls.tmp / "run"

        # (1) the original computation, called directly, with the original's own
        #     module-level context -- this is the reference
        original._CTX.update(cfg=cls.cfg, base=cls.base)
        cls.direct = {i: original._backgroundTask((i, cls.samples[i])) for i in FAST}

        # (2) the same backgrounds through the checkpointed worker path
        for i in FAST:
            ckpt.runOneBackground(cls.cfg, cls.ckdir, i)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def testCheckpointRowsAreNumericallyIdenticalToDirectCall(self):
        for i in FAST:
            with self.subTest(background=i):
                rows_direct, _meta = self.direct[i]
                on_disk = json.loads((ckpt.backgroundDir(self.ckdir, i) / "rows.json").read_text())
                # canonicalJson round-trips float64 through repr, so string
                # equality here IS bitwise equality of every value
                self.assertEqual(canonicalJson(rows_direct), canonicalJson(on_disk))

    def testCheckpointMetaIsNumericallyIdenticalToDirectCall(self):
        for i in FAST:
            with self.subTest(background=i):
                _rows, meta_direct = self.direct[i]
                on_disk = json.loads((ckpt.backgroundDir(self.ckdir, i) / "meta.json").read_text())
                self.assertEqual(canonicalJson(meta_direct), canonicalJson(on_disk))

    def testFastUsableBackgroundsReallyProduceTheFullFactorial(self):
        """if these degenerated to empty, the identity test above would be vacuous."""
        n_cells = len(original.PAIRS) * len(self.cfg["burden_levels"]) * len(self.cfg["capacity_levels"])
        for i in FAST_USABLE:
            rows, meta = self.direct[i]
            self.assertTrue(meta["usable"], f"background {i} expected usable")
            self.assertEqual(len(rows), n_cells)
        for i in FAST_UNUSABLE:
            rows, meta = self.direct[i]
            self.assertFalse(meta["usable"])
            self.assertEqual(rows, [])

    def testResumeDoesNotRecomputeAValidCheckpoint(self):
        stamps = {}
        for i in FAST:
            d = ckpt.backgroundDir(self.ckdir, i)
            stamps[i] = {p.name: (p.stat().st_mtime_ns, p.read_bytes())
                         for p in sorted(d.iterdir())}

        for i in FAST:
            ckpt.runOneBackground(self.cfg, self.ckdir, i)          # second call

        for i in FAST:
            d = ckpt.backgroundDir(self.ckdir, i)
            for p in sorted(d.iterdir()):
                self.assertIn(p.name, stamps[i])
                mtime, content = stamps[i][p.name]
                # an untouched file proves the work was skipped, not redone
                self.assertEqual(p.stat().st_mtime_ns, mtime, f"{i}/{p.name} rewritten")
                self.assertEqual(p.read_bytes(), content, f"{i}/{p.name} changed")

    def testForceRecomputesAndStillAgreesWithTheDirectCall(self):
        i = FAST_UNUSABLE[0]
        d = ckpt.backgroundDir(self.ckdir, i)
        before = (d / "checkpoint.json").stat().st_mtime_ns
        ckpt.runOneBackground(self.cfg, self.ckdir, i, force=True)
        self.assertNotEqual((d / "checkpoint.json").stat().st_mtime_ns, before)
        rows_direct, meta_direct = self.direct[i]
        self.assertEqual(canonicalJson(rows_direct),
                         canonicalJson(json.loads((d / "rows.json").read_text())))
        self.assertEqual(canonicalJson(meta_direct),
                         canonicalJson(json.loads((d / "meta.json").read_text())))

    # ------------------------------------------------------------------
    # checkpoint rejection
    # ------------------------------------------------------------------

    def _copyCheckpoint(self, i):
        dst = Path(tempfile.mkdtemp(prefix="d_ckpt_bad_", dir=self.tmp))
        shutil.copytree(ckpt.backgroundDir(self.ckdir, i),
                        ckpt.backgroundDir(dst, i))
        return dst

    def testUncorruptedCopyIsAccepted(self):
        """control for the rejection tests below."""
        i = FAST_UNUSABLE[0]
        dst = self._copyCheckpoint(i)
        ok, reason, payload = ckpt.checkpointStatus(dst, i, self.identity, self.samples[i])
        self.assertTrue(ok, reason)
        self.assertIsNotNone(payload)

    def testCorruptedRowsFileIsRejected(self):
        i = FAST_USABLE[0]
        dst = self._copyCheckpoint(i)
        p = ckpt.backgroundDir(dst, i) / "rows.json"
        blob = json.loads(p.read_text())
        blob[0]["burden_12"] = float(blob[0]["burden_12"]) + 1.0    # silent tamper
        p.write_text(canonicalJson(blob) + "\n")
        ok, reason, _ = ckpt.checkpointStatus(dst, i, self.identity, self.samples[i])
        self.assertFalse(ok)
        self.assertEqual(reason, "rows.json hash mismatch")

    def testTruncatedRowsFileIsRejected(self):
        i = FAST_USABLE[0]
        dst = self._copyCheckpoint(i)
        p = ckpt.backgroundDir(dst, i) / "rows.json"
        p.write_text(p.read_text()[: len(p.read_text()) // 2])
        ok, reason, _ = ckpt.checkpointStatus(dst, i, self.identity, self.samples[i])
        self.assertFalse(ok)
        self.assertEqual(reason, "rows.json hash mismatch")

    def testCorruptedMetaFileIsRejected(self):
        i = FAST_UNUSABLE[0]
        dst = self._copyCheckpoint(i)
        p = ckpt.backgroundDir(dst, i) / "meta.json"
        p.write_text(p.read_text().replace("false", "true", 1))
        ok, reason, _ = ckpt.checkpointStatus(dst, i, self.identity, self.samples[i])
        self.assertFalse(ok)
        self.assertEqual(reason, "meta.json hash mismatch")

    def testUnparseableCheckpointIsRejected(self):
        i = FAST_UNUSABLE[0]
        dst = self._copyCheckpoint(i)
        (ckpt.backgroundDir(dst, i) / "checkpoint.json").write_text("{not json")
        ok, reason, _ = ckpt.checkpointStatus(dst, i, self.identity, self.samples[i])
        self.assertFalse(ok)
        self.assertIn("unreadable checkpoint.json", reason)

    def testConfigHashMismatchIsRejected(self):
        i = FAST_UNUSABLE[0]
        dst = self._copyCheckpoint(i)
        other = ckpt.runIdentity(dict(self.cfg, burden_threshold=2.0), self.samples)
        ok, reason, _ = ckpt.checkpointStatus(dst, i, other, self.samples[i])
        self.assertFalse(ok)
        self.assertEqual(reason, "config hash mismatch")

    def testSourceHashMismatchIsRejected(self):
        i = FAST_UNUSABLE[0]
        dst = self._copyCheckpoint(i)
        other = dict(self.identity, source_hash="0" * 64)
        ok, reason, _ = ckpt.checkpointStatus(dst, i, other, self.samples[i])
        self.assertFalse(ok)
        self.assertEqual(reason, "source hash mismatch")

    def testSampleHashMismatchIsRejected(self):
        i = FAST_UNUSABLE[0]
        dst = self._copyCheckpoint(i)
        wrong = dict(self.samples[i])
        wrong["nu"] = wrong["nu"] * 1.000001
        ok, reason, _ = ckpt.checkpointStatus(dst, i, self.identity, wrong)
        self.assertFalse(ok)
        self.assertEqual(reason, "sample hash mismatch")

    def testBackgroundIndexMismatchIsRejected(self):
        """a checkpoint filed under the wrong index must not be adopted."""
        i, j = FAST_UNUSABLE
        dst = Path(tempfile.mkdtemp(prefix="d_ckpt_idx_", dir=self.tmp))
        shutil.copytree(ckpt.backgroundDir(self.ckdir, i), ckpt.backgroundDir(dst, j))
        ok, reason, _ = ckpt.checkpointStatus(dst, j, self.identity, self.samples[j])
        self.assertFalse(ok)
        self.assertEqual(reason, "background index mismatch")

    def testMissingCompletionMarkerIsRejected(self):
        i = FAST_UNUSABLE[0]
        dst = self._copyCheckpoint(i)
        (ckpt.backgroundDir(dst, i) / "DONE").unlink()
        ok, reason, _ = ckpt.checkpointStatus(dst, i, self.identity, self.samples[i])
        self.assertFalse(ok)
        self.assertEqual(reason, "no completion marker")

    def testStaleCompletionMarkerIsRejected(self):
        """a DONE that certifies some other checkpoint must not vouch for this one."""
        i = FAST_UNUSABLE[0]
        dst = self._copyCheckpoint(i)
        (ckpt.backgroundDir(dst, i) / "DONE").write_text("0" * 64 + "\n")
        ok, reason, _ = ckpt.checkpointStatus(dst, i, self.identity, self.samples[i])
        self.assertFalse(ok)
        self.assertEqual(reason, "completion marker does not match checkpoint.json")

    def testMissingPayloadIsRejected(self):
        i = FAST_UNUSABLE[0]
        dst = self._copyCheckpoint(i)
        (ckpt.backgroundDir(dst, i) / "rows.json").unlink()
        ok, reason, _ = ckpt.checkpointStatus(dst, i, self.identity, self.samples[i])
        self.assertFalse(ok)
        self.assertEqual(reason, "missing rows.json")

    def testRejectedCheckpointIsRecomputedNotTrusted(self):
        i = FAST_UNUSABLE[0]
        dst = self._copyCheckpoint(i)
        p = ckpt.backgroundDir(dst, i) / "meta.json"
        p.write_text(canonicalJson({"background": i, "error": "", "usable": True}) + "\n")
        ckpt.runOneBackground(self.cfg, dst, i)      # must not adopt the tampered file
        ok, reason, payload = ckpt.checkpointStatus(dst, i, self.identity, self.samples[i])
        self.assertTrue(ok, reason)
        self.assertEqual(canonicalJson(payload["meta"]), canonicalJson(self.direct[i][1]))

    # ------------------------------------------------------------------
    # merge
    # ------------------------------------------------------------------

    def testMergedSummaryMatchesOriginalSummaryLogic(self):
        """merge output vs. what run_experiment_d.main would write for these rows."""
        merged = ckpt.mergeCheckpoints(self.cfg, self.ckdir, self.identity, self.samples)

        # rebuild the original summary from the DIRECT results, in the same
        # ascending-background order Pool.map would have preserved
        rows, metas = [], []
        for i in FAST:
            r, m = self.direct[i]
            rows += r
            metas.append(m)
        df, dfm = pd.DataFrame(rows), pd.DataFrame(metas)
        expected = dict(
            experiment="D",
            n_backgrounds=int(len(dfm)),
            n_usable_backgrounds=int(dfm["usable"].sum()) if "usable" in dfm else 0,
            n_errors=int((dfm["error"] != "").sum()) if "error" in dfm else 0,
            unusable_reasons=(dfm.loc[~dfm["usable"], "reason"].value_counts().to_dict()
                              if "reason" in dfm and (~dfm["usable"]).any() else {}),
            n_cells=int(len(df)),
            burden_threshold=self.cfg["burden_threshold"],
            by_pair={label: original._pairSummary(df[df["pair"] == label])
                     for label, _, _ in original.PAIRS} if len(df) else {},
            overall=original._pairSummary(df) if len(df) else {},
        )

        got = {k: v for k, v in merged.items() if k != "recovery"}
        self.assertEqual(canonicalJson(expected), canonicalJson(got))

        # and the same for what actually landed on disk
        on_disk = json.loads((self.ckdir / "summary.json").read_text())
        on_disk.pop("recovery")
        self.assertEqual(canonicalJson(expected), canonicalJson(on_disk))

    def testMergedInteractionsTableMatchesADirectWrite(self):
        ckpt.mergeCheckpoints(self.cfg, self.ckdir, self.identity, self.samples)
        rows = []
        for i in FAST:
            rows += self.direct[i][0]
        ref = Path(tempfile.mkdtemp(prefix="d_ckpt_ref_", dir=self.tmp)) / "interactions.tsv"
        from proteostasis.provenance import writeTable
        writeTable(pd.DataFrame(rows), ref)
        self.assertEqual(ref.read_bytes(), (self.ckdir / "interactions.tsv").read_bytes())

    def testMergeRecordsRecoveryCountsSeparately(self):
        merged = ckpt.mergeCheckpoints(self.cfg, self.ckdir, self.identity, self.samples)
        rec = merged["recovery"]
        self.assertEqual(rec["n_backgrounds_requested"], 60)
        self.assertEqual(rec["n_backgrounds_completed"], len(FAST))
        self.assertEqual(rec["n_unresolved_timeout"], 0)
        self.assertEqual(rec["n_failed_process"], 0)
        self.assertEqual(rec["config_hash"], self.identity["config_hash"])
        self.assertEqual(rec["source_hash"], self.identity["source_hash"])
        self.assertEqual(rec["sample_matrix_hash"], self.identity["sample_matrix_hash"])

    def testCheckpointRecordsTheExactSampleParameters(self):
        for i in FAST:
            d = ckpt.backgroundDir(self.ckdir, i)
            ck = json.loads((d / "checkpoint.json").read_text())
            self.assertEqual(ck["background_index"], i)
            self.assertEqual(ck["sample"], self.samples[i])
            self.assertEqual(ck["sample_hash"], hashObject(self.samples[i]))
            meta = json.loads((d / "meta.json").read_text())
            # sample_index survives into the background row as `background`
            self.assertEqual(meta["background"], i)
            for k, v in self.samples[i].items():
                self.assertEqual(meta[f"p_{k}"], float(v))


class TestTimeoutClassification(unittest.TestCase):
    """a background that runs out of wall clock is unresolved, not failed."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="d_ckpt_timeout_"))
        # 0.5 s cannot even import numpy, so every background is guaranteed to
        # hit the wall limit -- this exercises the real kill path, not a mock
        cls.proc = subprocess.run(
            [sys.executable, str(RUNNER), "--config", str(CONFIG),
             "--outdir", str(cls.tmp), "--n-backgrounds", "3",
             "--concurrency", "3", "--background-timeout", "0.5"],
            cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=600)
        cls.summary = json.loads((cls.tmp / "summary.json").read_text())

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def testControllerExitsCleanly(self):
        self.assertEqual(self.proc.returncode, 0, self.proc.stderr[-2000:])

    def testAllThreeAreClassifiedUnresolvedTimeout(self):
        rec = self.summary["recovery"]
        self.assertEqual(rec["n_unresolved_timeout"], 3)
        self.assertEqual(rec["unresolved_timeout_backgrounds"], [0, 1, 2])
        self.assertEqual(rec["n_backgrounds_completed"], 0)

    def testTimeoutIsNotCountedAsANumericalFailure(self):
        self.assertEqual(self.summary["n_errors"], 0)
        self.assertEqual(self.summary["recovery"]["n_failed_process"], 0)
        self.assertEqual(self.summary["unusable_reasons"], {})

    def testTimedOutBackgroundsContributeNoRowsAndNoSummary(self):
        self.assertEqual(self.summary["n_cells"], 0)
        self.assertEqual(self.summary["n_backgrounds"], 0)
        self.assertEqual(self.summary["overall"], {})
        self.assertEqual(self.summary["by_pair"], {})

    def testUnresolvedMarkersAreWrittenAndSayWhy(self):
        for i in range(3):
            p = self.tmp / "unresolved" / f"background_{i:04d}.json"
            self.assertTrue(p.exists(), p)
            rec = json.loads(p.read_text())
            self.assertEqual(rec["status"], "unresolved_timeout")
            self.assertEqual(rec["timeout_s"], 0.5)
            self.assertIn("not a numerical failure", rec["note"])

    def testNoCheckpointIsLeftBehindForATimedOutBackground(self):
        for i in range(3):
            self.assertFalse((ckpt.backgroundDir(self.tmp, i) / "DONE").exists())

    def testProgressFileRecordsTheTimeouts(self):
        prog = json.loads((self.tmp / "progress.json").read_text())
        self.assertEqual(sorted(prog["unresolved_timeout"]), [0, 1, 2])
        self.assertEqual(prog["running"], [])
        self.assertEqual(prog["pending"], [])


class TestBackgroundSelection(unittest.TestCase):
    """--only restricts scheduling; it must not perturb the sample matrix."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="d_ckpt_only_"))
        cls.proc = subprocess.run(
            [sys.executable, str(RUNNER), "--config", str(CONFIG),
             "--outdir", str(cls.tmp), "--only", ",".join(str(i) for i in FAST_UNUSABLE),
             "--concurrency", "2", "--background-timeout", "300"],
            cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=900)

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def testOnlySelectedBackgroundsAreEvaluated(self):
        self.assertEqual(self.proc.returncode, 0, self.proc.stderr[-2000:])
        got = sorted(int(p.name.split("_")[1])
                     for p in (self.tmp / "backgrounds").iterdir())
        self.assertEqual(got, sorted(FAST_UNUSABLE))

    def testSelectionKeepsTheFullSizeSampleMatrix(self):
        """the draws must be the ones a 60-background run would use."""
        manifest = json.loads((self.tmp / "run_manifest.json").read_text())
        cfg = _loadCfg()
        _b, samples = ckpt.sampleMatrix(cfg)
        self.assertEqual(manifest["n_backgrounds"], 60)
        self.assertEqual(manifest["sample_matrix_hash"], hashObject(samples))
        self.assertEqual(manifest["selected_backgrounds"], sorted(FAST_UNUSABLE))
        for i in FAST_UNUSABLE:
            ck = json.loads((ckpt.backgroundDir(self.tmp, i) / "checkpoint.json").read_text())
            self.assertEqual(ck["sample"], samples[i])

    def testOutOfRangeSelectionIsRejected(self):
        proc = subprocess.run(
            [sys.executable, str(RUNNER), "--config", str(CONFIG),
             "--outdir", str(self.tmp), "--only", "999"],
            cwd=str(REPO_ROOT), capture_output=True, text=True, timeout=300)
        self.assertNotEqual(proc.returncode, 0)
        self.assertIn("outside 0..59", proc.stderr)


class TestSourceIdentity(unittest.TestCase):

    def testSourceHashCoversTheScientificModules(self):
        h = ckpt.sourceHashes()
        for rel in ("scripts/run_experiment_d.py", "scripts/proteostasis/model.py",
                    "scripts/proteostasis/simulate.py", "scripts/proteostasis/equilibria.py",
                    "scripts/proteostasis/sweeps.py"):
            self.assertIn(rel, h)
            self.assertEqual(len(h[rel]), 64)

    def testSourceHashExcludesTheOrchestratorItself(self):
        """orchestration edits must not invalidate finished science."""
        self.assertNotIn("scripts/run_experiment_d_checkpointed.py", ckpt.sourceHashes())

    def testRunnerReusesTheOriginalScientificFunctions(self):
        """equivalence is structural: there is no second copy of the model."""
        self.assertIs(ckpt.original, sys.modules["run_experiment_d"])
        self.assertIs(ckpt.original._backgroundTask, original._backgroundTask)
        self.assertIs(ckpt.original._pairSummary, original._pairSummary)
        self.assertEqual(ckpt.original.PAIRS, original.PAIRS)
        src = RUNNER.read_text()
        for forbidden in ("def _score(", "def _perturb(", "def _backgroundTask(",
                          "def _pairSummary("):
            self.assertNotIn(forbidden, src,
                             f"{forbidden} reimplemented in the recovery runner")


if __name__ == "__main__":
    unittest.main()
