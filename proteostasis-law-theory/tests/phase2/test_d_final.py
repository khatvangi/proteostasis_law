"""experiment D closure: sign conventions, background-level inference, hashes.

three kinds of test live here and they fail for different reasons.

MODEL-FREE tests pin the statistics themselves -- the sign convention of each
null, the fact that a grouped bootstrap really does resample groups, the exact
arithmetic of the sensitivity bounds, and the decision rule's insistence on
three coherent legs.  these run on a clean checkout and are the ones that pin
the method.

DESIGN tests pin the analysis constants against `configs/phase1/experiment_d.json`
and `run_experiment_d.PAIRS`, so the closure cannot silently drift from the
experiment it summarises.

ARTIFACT tests assert against the gitignored run root and its D_final outputs.
they skip when those are absent rather than fail, which is the same convention
the rest of `tests/phase2` uses.
"""

import json
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

import _context  # noqa: F401
from _context import REPO_ROOT
import check_d_closure as closure
import d_final as D
import run_experiment_d as original
from proteostasis.provenance import hashFile, hashObject, loadConfig

PINNED = json.loads((REPO_ROOT / "scripts" / "phase2" / "D_RUN_HASHES.json").read_text())
RUN_ROOT = REPO_ROOT / PINNED["run_root"]
D_FINAL = (REPO_ROOT / "results" / "phase2" / "closure_20260731T220024-0500" / "D_final")


def _frame(background, excess, column="excess_additive", **extra):
    """minimal frame with the columns the helpers actually read."""
    n = len(background)
    d = pd.DataFrame({"background": background, column: excess,
                      "synthetic_collapse": np.zeros(n, dtype=bool)})
    for k, v in extra.items():
        d[k] = v
    return d


# --------------------------------------------------------------------------
# sign conventions -- task 1
# --------------------------------------------------------------------------

class TestNullSignConventions(unittest.TestCase):

    def testThreeNullsWithTheDeclaredDirections(self):
        self.assertEqual([n.name for n in D.NULLS],
                         ["additive", "multiplicative", "bliss"])
        self.assertEqual(D.NULL_BY_NAME["additive"].worse, "greater")
        self.assertEqual(D.NULL_BY_NAME["multiplicative"].worse, "greater")
        self.assertEqual(D.NULL_BY_NAME["bliss"].worse, "less")

    def testWorseDirectionMatchesTheOriginalPairSummary(self):
        """`_pairSummary` is the definition; this module must not restate it.

        the original hardcodes (null, column, direction) triples.  if either
        side is ever edited, this test is what catches the divergence.
        """
        import inspect
        src = inspect.getsource(original._pairSummary)
        self.assertIn('("additive", "excess_additive", "greater")', src)
        self.assertIn('("multiplicative", "log_excess_multiplicative", "greater")', src)
        self.assertIn('("bliss", "excess_bliss", "less")', src)
        for n in D.NULLS:
            self.assertIn(f'"{n.excess_col}", "{n.worse}"', src)

    def testTiesAreNeverCountedAsWorse(self):
        z = np.array([0.0, 0.0, 0.0])
        for n in D.NULLS:
            self.assertEqual(n.worseMask(z).sum(), 0, n.name)

    def testBlissNegativeMedianAndMajorityWorseAreTheSameStatement(self):
        """the apparent paradox in the raw summary, reproduced and resolved.

        70 % of the values are negative.  on the survival scale negative IS
        worse, so `frac_worse > 0.5` and `median < 0` agree; they would only
        conflict if the additive convention were wrongly applied to bliss.
        """
        v = np.concatenate([np.full(70, -0.01), np.full(30, +0.01)])
        bliss = D.NULL_BY_NAME["bliss"]
        self.assertAlmostEqual(float(bliss.worseMask(v).mean()), 0.70)
        self.assertLess(float(np.median(v)), 0.0)
        # and the additive reading of the same numbers gives the opposite answer
        self.assertAlmostEqual(float(D.NULL_BY_NAME["additive"].worseMask(v).mean()), 0.30)

    def testWorseSignIsPlusOneForBurdenScalesAndMinusOneForSurvival(self):
        self.assertEqual(D.NULL_BY_NAME["additive"].worseSign, 1.0)
        self.assertEqual(D.NULL_BY_NAME["multiplicative"].worseSign, 1.0)
        self.assertEqual(D.NULL_BY_NAME["bliss"].worseSign, -1.0)


# --------------------------------------------------------------------------
# design constants -- task 2
# --------------------------------------------------------------------------

class TestDesignConstants(unittest.TestCase):

    def setUp(self):
        self.cfg = loadConfig(REPO_ROOT / PINNED["config_path"])

    def testLevelsMatchTheExperimentConfig(self):
        self.assertEqual(list(D.BURDEN_LEVELS), list(self.cfg["burden_levels"]))
        self.assertEqual(sorted(D.CAPACITY_LEVELS), sorted(self.cfg["capacity_levels"]))

    def testCellArithmetic(self):
        self.assertEqual(D.CELLS_PER_PAIR, 36)
        self.assertEqual(D.DOUBLE_CELLS_PER_PAIR, 25)
        # exactly one level per axis is the unperturbed 1.0, which is what makes
        # 36 - 25 = 11 the single-perturbation count
        self.assertEqual(D.BURDEN_LEVELS.count(1.0), 1)
        self.assertEqual(D.CAPACITY_LEVELS.count(1.0), 1)
        self.assertEqual(D.CELLS_PER_PAIR - D.DOUBLE_CELLS_PER_PAIR, 11)

    def testPairsMatchTheExperimentModule(self):
        self.assertEqual(tuple(p[0] for p in original.PAIRS), tuple(D.PAIRS))

    def testCountsInThePinnedManifestAreSelfConsistent(self):
        c = PINNED["counts"]
        self.assertEqual(c["n_usable_backgrounds"] + c["n_unusable_backgrounds"],
                         c["n_backgrounds_completed"])
        self.assertEqual(c["n_backgrounds_completed"] + c["n_unresolved_timeout"]
                         + c["n_failed_process"], c["n_backgrounds_requested"])
        self.assertEqual(c["n_cells"],
                         c["n_usable_backgrounds"] * len(D.PAIRS) * D.CELLS_PER_PAIR)
        self.assertEqual(c["n_double_cells"],
                         c["n_usable_backgrounds"] * len(D.PAIRS) * D.DOUBLE_CELLS_PER_PAIR)
        self.assertEqual(c["n_unresolved_timeout"],
                         len(c["unresolved_timeout_backgrounds"]))


# --------------------------------------------------------------------------
# background-level summaries and subsets -- task 3
# --------------------------------------------------------------------------

class TestPerBackground(unittest.TestCase):

    def testFractionMedianAndMajorityAreComputedWithinBackground(self):
        d = _frame([1, 1, 1, 2, 2, 2], [1.0, 2.0, -1.0, -1.0, -2.0, 3.0])
        pb = D.perBackground(d, D.NULL_BY_NAME["additive"])
        self.assertEqual(list(pb.index), [1, 2])
        self.assertAlmostEqual(pb.loc[1, "frac_worse"], 2 / 3)
        self.assertAlmostEqual(pb.loc[2, "frac_worse"], 1 / 3)
        self.assertAlmostEqual(pb.loc[1, "median_excess"], 1.0)
        self.assertAlmostEqual(pb.loc[2, "median_excess"], -1.0)
        self.assertTrue(pb.loc[1, "majority_worse"])
        self.assertFalse(pb.loc[2, "majority_worse"])

    def testExactlyHalfIsNotAMajority(self):
        d = _frame([1, 1], [1.0, -1.0])
        pb = D.perBackground(d, D.NULL_BY_NAME["additive"])
        self.assertAlmostEqual(pb.loc[1, "frac_worse"], 0.5)
        self.assertFalse(bool(pb.loc[1, "majority_worse"]))

    def testNonFiniteValuesAreDroppedNotZeroed(self):
        d = _frame([1, 1, 1], [np.inf, np.nan, 2.0])
        pb = D.perBackground(d, D.NULL_BY_NAME["additive"])
        self.assertEqual(int(pb.loc[1, "n_cells"]), 1)
        self.assertAlmostEqual(pb.loc[1, "median_excess"], 2.0)

    def testBackgroundWithNoCellsIsDroppedEntirely(self):
        d = _frame([1, 2], [np.nan, 2.0])
        pb = D.perBackground(d, D.NULL_BY_NAME["additive"])
        self.assertEqual(list(pb.index), [2])


class TestSubsets(unittest.TestCase):

    def _d(self):
        return pd.DataFrame({
            "background": [1, 1, 1, 1],
            "burden_0": [1.0, 1.0, 1.0, 1.0],
            "burden_1": [2.0, 0.5, 2.0, 2.0],
            "burden_2": [3.0, 3.0, 0.5, 3.0],
            "survival_1": [0.5, 0.0, 0.5, 0.5],
            "survival_2": [0.5, 0.5, 0.0, 0.5],
            "censored_12": [False, False, False, True],
            "status_12": ["converged", "converged", "timeout", "blowup"],
        })

    def testNamedSubsets(self):
        d = self._d()
        np.testing.assert_array_equal(D.subsetMask(d, "all"), [True] * 4)
        np.testing.assert_array_equal(D.subsetMask(d, "both_singles_damaging"),
                                      [True, False, False, True])
        np.testing.assert_array_equal(D.subsetMask(d, "bliss_informative"),
                                      [True, False, False, True])
        np.testing.assert_array_equal(D.subsetMask(d, "uncensored"),
                                      [True, True, True, False])
        np.testing.assert_array_equal(D.subsetMask(d, "clean"),
                                      [True, False, False, False])
        np.testing.assert_array_equal(D.subsetMask(d, "converged"),
                                      [True, True, False, False])

    def testUnknownSubsetRaises(self):
        with self.assertRaises(ValueError):
            D.subsetMask(self._d(), "whatever")

    def testWhenBothSinglesAreDamagingTheMultiplicativeNullIsTheStricterOne(self):
        """the arithmetic lemma behind the `both_singles_damaging` subset.

        with x = b1/b0 and y = b2/b0, mult_pred - add_pred = b0*(x-1)*(y-1).
        so the multiplicative prediction sits at or above the additive one
        exactly when both singles push the same way.  where a single is
        PROTECTIVE the multiplicative null predicts LESS burden and is therefore
        easier to exceed -- which is why a higher `frac_worse` under the
        multiplicative null is not automatically stronger evidence.
        """
        rng = np.random.default_rng(0)
        b0 = rng.uniform(0.01, 1.0, 5000)
        x = rng.uniform(0.5, 3.0, 5000)
        y = rng.uniform(0.5, 3.0, 5000)
        add = b0 * x + b0 * y - b0
        mult = b0 * x * y
        both_up = (x > 1) & (y > 1)
        self.assertTrue(bool((mult[both_up] >= add[both_up]).all()))
        one_down = ((x - 1) * (y - 1)) < 0
        self.assertTrue(bool((mult[one_down] < add[one_down]).all()))


# --------------------------------------------------------------------------
# grouped bootstrap
# --------------------------------------------------------------------------

class TestGroupedBootstrap(unittest.TestCase):

    def testFastPathTakesLiterallyTheSameResamplesAsTheReference(self):
        rng = np.random.default_rng(7)
        groups = [rng.normal(size=25) for _ in range(46)]
        ref = D.clusterBootstrap(groups, D.statMedianOfGroupMedians, n_boot=500, seed=11)
        fast = D.scalarBootstrap([float(np.median(g)) for g in groups], D.statMedian,
                                 n_boot=500, seed=11)
        for k in ("point", "lo", "hi"):
            self.assertAlmostEqual(ref[k], fast[k], places=12, msg=k)

    def testFastPathMatchesTheReferenceForTheMajorityStatistic(self):
        rng = np.random.default_rng(3)
        groups = [rng.normal(0.2, 1.0, 25) for _ in range(46)]
        ref = D.clusterBootstrap(groups, D.statFractionOfGroupsWithMajority(1.0),
                                 n_boot=500, seed=5)
        fast = D.scalarBootstrap([float((g > 0).mean()) for g in groups],
                                 D.statFractionAboveHalf, n_boot=500, seed=5)
        for k in ("point", "lo", "hi"):
            self.assertAlmostEqual(ref[k], fast[k], places=12, msg=k)

    def testPooledRatePointEstimateIsTheCellLevelFraction(self):
        k = np.array([10, 0, 25, 3], dtype=float)
        n = np.array([25, 25, 25, 25], dtype=float)
        r = D.pooledRateBootstrap(k, n, n_boot=200, seed=1)
        self.assertAlmostEqual(r["point"], k.sum() / n.sum())

    def testSameSeedReproducesTheIntervalExactly(self):
        v = np.linspace(-1, 1, 46)
        a = D.scalarBootstrap(v, D.statMedian, n_boot=300, seed=42)
        b = D.scalarBootstrap(v, D.statMedian, n_boot=300, seed=42)
        self.assertEqual(a, b)

    def testResamplingGroupsIsWiderThanPretendingCellsAreIndependent(self):
        """the whole reason this module exists.

        each background here is 25 copies of one value, i.e. perfectly
        correlated within background -- the extreme of what the real design has.
        a cell-level bootstrap would see 1150 observations and produce a much
        narrower interval than the 46 genuinely independent ones support.
        """
        rng = np.random.default_rng(1)
        per_bg = rng.normal(0.3, 1.0, 46)
        grouped = D.pooledRateBootstrap((per_bg > 0).astype(float) * 25,
                                        np.full(46, 25.0), n_boot=2000, seed=9)
        cells = np.repeat((per_bg > 0).astype(float), 25)
        naive = D.pooledRateBootstrap(cells, np.ones(len(cells)), n_boot=2000, seed=9)
        self.assertAlmostEqual(grouped["point"], naive["point"], places=12)
        self.assertGreater(grouped["hi"] - grouped["lo"],
                           3.0 * (naive["hi"] - naive["lo"]))

    def testEmptyInputIsNotAnError(self):
        r = D.scalarBootstrap([], D.statMedian, n_boot=10, seed=1)
        self.assertEqual(r["n_groups"], 0)
        self.assertTrue(np.isnan(r["point"]))


class TestSignTest(unittest.TestCase):

    def testTiesAreDiscardedAndTheCountsAreReported(self):
        v = np.array([1.0, 1.0, 1.0, -1.0, 0.0, 0.0])
        r = D.signTestBackgrounds(v, D.NULL_BY_NAME["additive"])
        self.assertEqual((r["n_worse"], r["n_better"], r["n_tied"]), (3, 1, 2))

    def testBlissCountsNegativesAsWorse(self):
        v = np.array([-1.0, -1.0, 1.0])
        r = D.signTestBackgrounds(v, D.NULL_BY_NAME["bliss"])
        self.assertEqual((r["n_worse"], r["n_better"]), (2, 1))

    def testBalancedEvidenceGivesPOne(self):
        v = np.array([1.0, -1.0, 1.0, -1.0])
        self.assertAlmostEqual(
            D.signTestBackgrounds(v, D.NULL_BY_NAME["additive"])["p_two_sided"], 1.0)

    def testAllTiedGivesNoPValueRatherThanAFalseOne(self):
        r = D.signTestBackgrounds(np.zeros(5), D.NULL_BY_NAME["additive"])
        self.assertTrue(np.isnan(r["p_two_sided"]))


# --------------------------------------------------------------------------
# sensitivity bounds -- task 4
# --------------------------------------------------------------------------

class TestUnresolvedBounds(unittest.TestCase):

    def testTheActualClosureNumbers(self):
        b = D.unresolvedBounds(43, n_usable=46, n_unresolved=2, n_unusable=12,
                               n_requested=60)
        self.assertAlmostEqual(b["conditional"], 43 / 46)
        self.assertAlmostEqual(b["usable_lo"], 43 / 48)
        self.assertAlmostEqual(b["usable_hi"], 45 / 48)
        self.assertAlmostEqual(b["requested_lo"], 43 / 60)
        self.assertAlmostEqual(b["requested_hi"], 45 / 60)

    def testTheConditionalEstimateAlwaysLiesInsideTheUsableBounds(self):
        for k in range(0, 47):
            b = D.unresolvedBounds(k, 46, 2, 12, 60)
            self.assertLessEqual(b["usable_lo"], b["conditional"] + 1e-12)
            self.assertGreaterEqual(b["usable_hi"], b["conditional"] - 1e-12)

    def testRequestedBoundsAreNeverAboveTheUsableOnes(self):
        for k in (0, 13, 28, 43, 46):
            b = D.unresolvedBounds(k, 46, 2, 12, 60)
            self.assertLessEqual(b["requested_lo"], b["usable_lo"] + 1e-12)
            self.assertLessEqual(b["requested_hi"], b["usable_hi"] + 1e-12)

    def testNothingIsImputed(self):
        """with zero unresolved backgrounds the interval collapses to a point."""
        b = D.unresolvedBounds(43, 46, 0, 12, 58)
        self.assertAlmostEqual(b["usable_lo"], b["usable_hi"])
        self.assertAlmostEqual(b["usable_lo"], 43 / 46)

    def testImpossibleNumeratorRaises(self):
        with self.assertRaises(ValueError):
            D.unresolvedBounds(47, 46, 2, 12, 60)


# --------------------------------------------------------------------------
# the decision rule -- task 5
# --------------------------------------------------------------------------

class TestSupportDecision(unittest.TestCase):

    ADD = D.NULL_BY_NAME["additive"]
    BLISS = D.NULL_BY_NAME["bliss"]

    def testAllThreeLegsRequired(self):
        good_med = dict(point=0.05, lo=0.01, hi=0.09)
        good_maj = dict(point=0.8, lo=0.65, hi=0.92)
        self.assertTrue(D.supportDecision(good_med, good_maj, self.ADD)["supported"])
        # CI touching zero
        r = D.supportDecision(dict(point=0.05, lo=-0.01, hi=0.09), good_maj, self.ADD)
        self.assertFalse(r["supported"])
        self.assertIn("median_ci_excludes_zero", r["failing"])
        # majority CI straddling one half
        r = D.supportDecision(good_med, dict(point=0.6, lo=0.45, hi=0.75), self.ADD)
        self.assertFalse(r["supported"])
        self.assertIn("majority_ci_above_half", r["failing"])

    def testThreeCriteriaAreDocumented(self):
        self.assertEqual(len(D.DECISION_CRITERIA), 3)

    def testBlissRequiresTheNegativeDirection(self):
        med = dict(point=-0.02, lo=-0.04, hi=-0.005)
        maj = dict(point=0.8, lo=0.65, hi=0.92)
        self.assertTrue(D.supportDecision(med, maj, self.BLISS)["supported"])
        flipped = dict(point=0.02, lo=0.005, hi=0.04)
        r = D.supportDecision(flipped, maj, self.BLISS)
        self.assertFalse(r["supported"])
        self.assertIn("median_direction", r["failing"])
        self.assertIn("median_ci_excludes_zero", r["failing"])

    def testAMajorityPointEstimateAboveHalfIsNotEnoughOnItsOwn(self):
        r = D.supportDecision(dict(point=0.05, lo=0.01, hi=0.09),
                              dict(point=0.61, lo=0.46, hi=0.76), self.ADD)
        self.assertFalse(r["supported"])


# --------------------------------------------------------------------------
# artifact identity -- source and result hashes
# --------------------------------------------------------------------------

class TestPinnedSourceHashes(unittest.TestCase):

    def testTrackedSourceFilesStillHashToThePinnedValues(self):
        """the D result is only attributable to this tree while these match."""
        for rel, want in PINNED["source_files"].items():
            self.assertEqual(hashFile(REPO_ROOT / rel), want, rel)

    def testCombinedSourceHash(self):
        self.assertEqual(hashObject(PINNED["source_files"]), PINNED["source_hash"])

    def testConfigStillHashesToThePinnedValue(self):
        self.assertEqual(hashObject(loadConfig(REPO_ROOT / PINNED["config_path"])),
                         PINNED["config_hash"])

    def testSampleMatrixStillReproduces(self):
        import run_experiment_d_checkpointed as ck
        _base, samples = ck.sampleMatrix(loadConfig(REPO_ROOT / PINNED["config_path"]))
        self.assertEqual(len(samples), PINNED["counts"]["n_backgrounds_requested"])
        self.assertEqual(hashObject(samples), PINNED["sample_matrix_hash"])


@unittest.skipUnless(RUN_ROOT.is_dir(), f"D run root absent: {RUN_ROOT}")
class TestRunArtifacts(unittest.TestCase):

    def testArtifactHashes(self):
        for name, want in PINNED["artifact_sha256"].items():
            self.assertEqual(hashFile(RUN_ROOT / name), want, name)

    def testSummaryCountsMatchThePinnedManifest(self):
        s = json.loads((RUN_ROOT / "summary.json").read_text())
        c = PINNED["counts"]
        self.assertEqual(s["n_usable_backgrounds"], c["n_usable_backgrounds"])
        self.assertEqual(s["n_backgrounds"], c["n_backgrounds_completed"])
        self.assertEqual(s["n_cells"], c["n_cells"])
        self.assertEqual(s["n_errors"], 0)
        self.assertEqual(s["recovery"]["n_unresolved_timeout"], c["n_unresolved_timeout"])
        self.assertEqual(s["recovery"]["unresolved_timeout_backgrounds"],
                         c["unresolved_timeout_backgrounds"])
        self.assertEqual(s["recovery"]["n_failed_process"], 0)
        self.assertEqual(s["unusable_reasons"], {c["unusable_reason"]:
                                                 c["n_unusable_backgrounds"]})

    def testTimeoutsAreNotCountedAsFailuresAnywhere(self):
        """the one semantic rule the whole recovery rests on."""
        s = json.loads((RUN_ROOT / "summary.json").read_text())
        self.assertEqual(s["n_errors"], 0)
        self.assertEqual(s["recovery"]["failed_process_backgrounds"], [])
        self.assertEqual(s["n_backgrounds"] + s["recovery"]["n_unresolved_timeout"],
                         s["recovery"]["n_backgrounds_requested"])
        for i in PINNED["counts"]["unresolved_timeout_backgrounds"]:
            rec = json.loads((RUN_ROOT / "unresolved" /
                              f"background_{i:04d}.json").read_text())
            self.assertEqual(rec["status"], "unresolved_timeout")


@unittest.skipUnless(RUN_ROOT.is_dir(), f"D run root absent: {RUN_ROOT}")
class TestBackgroundLevelReproducesFromCheckpoints(unittest.TestCase):
    """the D_final point estimates, recomputed here from the checkpoints.

    these do not read D_final's own output for the values they assert; they
    recompute them, then compare.  a stale or hand-edited D_final table fails.
    """

    @classmethod
    def setUpClass(cls):
        rows = []
        for p in sorted((RUN_ROOT / "backgrounds").glob("background_*/rows.json")):
            rows += json.loads(p.read_text())
        cls.df = pd.DataFrame(rows)
        cls.d = D.doubleCells(cls.df)

    def testCellCountsAreTheAdvertisedOnes(self):
        c = PINNED["counts"]
        self.assertEqual(len(self.df), c["n_cells"])
        self.assertEqual(len(self.d), c["n_double_cells"])
        self.assertEqual(self.df["background"].nunique(), c["n_usable_backgrounds"])
        self.assertEqual(set(self.d.groupby(["background", "pair"]).size().unique()),
                         {D.DOUBLE_CELLS_PER_PAIR})

    def testCellLevelFractionsReproduceTheShippedSummary(self):
        """the background-level analysis works on the same values the run did."""
        s = json.loads((RUN_ROOT / "summary.json").read_text())
        for null in D.NULLS:
            got = float(null.worseMask(self.d[null.excess_col].to_numpy()).mean())
            self.assertAlmostEqual(got, s["overall"][f"{null.name}_frac_worse_than_null"],
                                   places=12, msg=null.name)

    def testTheMultiplicativeSummaryFieldIsALogScaleMedian(self):
        """named `multiplicative_median_excess`, actually the median LOG excess.

        this is a labelling defect in `_pairSummary`, not a sign error: the
        column it summarises is `log_excess_multiplicative`.  pinned here so the
        closure keeps reporting the scale explicitly.
        """
        s = json.loads((RUN_ROOT / "summary.json").read_text())
        log_med = float(np.median(self.d["log_excess_multiplicative"].to_numpy()))
        lin_med = float(np.median(self.d["excess_multiplicative"].to_numpy()))
        self.assertAlmostEqual(s["overall"]["multiplicative_median_excess"], log_med,
                               places=12)
        self.assertNotAlmostEqual(s["overall"]["multiplicative_median_excess"], lin_med,
                                  places=6)

    def testCensoringIsExactlyTheEscapeSet(self):
        """so a censored cell understates, never overstates, supra-additivity."""
        self.assertTrue(bool((self.d["censored_12"].astype(bool)
                              == (self.d["status_12"] == "blowup")).all()))

    def testBlissTiesAreExactlyTheCellsWhereASingleIsAlreadyLethal(self):
        tie = self.d["excess_bliss"].to_numpy() == 0.0
        lethal = ((self.d["survival_1"].to_numpy() == 0.0)
                  | (self.d["survival_2"].to_numpy() == 0.0))
        self.assertTrue(bool(np.all(lethal[tie])))
        self.assertTrue(bool((self.d.loc[lethal, "excess_bliss"].to_numpy() >= 0).all()))

    def testSyntheticCollapseIsExactlyItsDefinition(self):
        want = (self.d["viable_1"].astype(bool) & self.d["viable_2"].astype(bool)
                & ~self.d["viable_12"].astype(bool))
        self.assertTrue(bool((self.d["synthetic_collapse"].astype(bool) == want).all()))

    def testHeadlineBackgroundLevelEstimates(self):
        """the numbers EXPERIMENT_D_FINAL.md reports, recomputed from scratch."""
        expect = {
            ("influx_x_total_capacity", "additive"): (46, 1.0),
            ("influx_x_total_capacity", "multiplicative"): (44, 0.9565217391304348),
            ("influx_x_total_capacity", "bliss"): (41, 0.8913043478260869),
            ("influx_x_chaperone_only", "additive"): (28, 0.6086956521739131),
            ("influx_x_chaperone_only", "multiplicative"): (35, 0.7608695652173914),
            ("influx_x_chaperone_only", "bliss"): (28, 0.6086956521739131),
            ("nascent_x_total_capacity", "additive"): (37, 0.8043478260869565),
            ("nascent_x_total_capacity", "multiplicative"): (38, 0.8260869565217391),
            ("nascent_x_total_capacity", "bliss"): (37, 0.8043478260869565),
        }
        for (pair, null_name), (n_major, frac) in expect.items():
            pb = D.perBackground(self.d[self.d["pair"] == pair], D.NULL_BY_NAME[null_name])
            self.assertEqual(len(pb), 46, f"{pair}/{null_name}")
            self.assertEqual(int(pb["majority_worse"].sum()), n_major, f"{pair}/{null_name}")
            self.assertAlmostEqual(float(pb["majority_worse"].mean()), frac, places=12)

    def testSyntheticCollapseByBackgroundAndByCell(self):
        pb = D.perBackgroundCollapse(self.d)
        self.assertEqual(int(pb["any_collapse"].sum()), 43)
        self.assertAlmostEqual(float(pb["any_collapse"].mean()), 43 / 46, places=12)
        self.assertAlmostEqual(float(pb["n_collapse"].sum() / pb["n_cells"].sum()),
                               0.22057971014492753, places=12)
        by_pair = {p: int(D.perBackgroundCollapse(self.d[self.d["pair"] == p])
                          ["any_collapse"].sum()) for p in D.PAIRS}
        self.assertEqual(by_pair, {"influx_x_total_capacity": 43,
                                   "influx_x_chaperone_only": 36,
                                   "nascent_x_total_capacity": 13})

    def testBootstrapIntervalsReproduceUnderThePinnedSeed(self):
        pb = D.perBackground(self.d[self.d["pair"] == "nascent_x_total_capacity"],
                             D.NULL_BY_NAME["additive"])
        ci = D.scalarBootstrap(pb["median_excess"].to_numpy(), D.statMedian,
                               D.N_BOOT, D.BOOTSTRAP_SEED)
        self.assertAlmostEqual(ci["point"], 0.001964866936136067, places=15)
        self.assertAlmostEqual(ci["lo"], 0.0006893231502037162, places=15)
        self.assertAlmostEqual(ci["hi"], 0.0040603597991905604, places=15)
        self.assertGreater(ci["lo"], 0.0)
        again = D.scalarBootstrap(pb["median_excess"].to_numpy(), D.statMedian,
                                  D.N_BOOT, D.BOOTSTRAP_SEED)
        self.assertEqual(ci, again)


@unittest.skipUnless(RUN_ROOT.is_dir() and (D_FINAL / "d_final_results.json").exists(),
                     f"D run root or D_final outputs absent: {RUN_ROOT}, {D_FINAL}")
class TestDFinalOutputsAgreeWithRecomputation(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.res = json.loads((D_FINAL / "d_final_results.json").read_text())
        rows = []
        for p in sorted((RUN_ROOT / "backgrounds").glob("background_*/rows.json")):
            rows += json.loads(p.read_text())
        cls.d = D.doubleCells(pd.DataFrame(rows))

    def testCountsBlock(self):
        c = self.res["counts"]
        self.assertEqual(c["n_usable"], 46)
        self.assertEqual(c["n_unusable"], 12)
        self.assertEqual(c["n_unresolved_timeout"], 2)
        self.assertEqual(c["n_failed_process"], 0)
        self.assertEqual(c["n_double_cells"], 3450)

    def testEveryPrimaryEstimateRecomputes(self):
        for r in self.res["estimates"]:
            if r["subset"] != "all":
                continue
            null = D.NULL_BY_NAME[r["null"]]
            pb = D.perBackground(self.d[self.d["pair"] == r["pair"]], null)
            self.assertEqual(int(len(pb)), r["n_backgrounds"])
            self.assertAlmostEqual(float(np.median(pb["median_excess"])),
                                   r["bg_median_excess"], places=12)
            self.assertAlmostEqual(float((pb["frac_worse"] > 0.5).mean()),
                                   r["frac_backgrounds_majority_worse"], places=12)
            self.assertAlmostEqual(float(pb["n_worse"].sum() / pb["n_cells"].sum()),
                                   r["cell_level_frac_worse"], places=12)

    def testVerdictsAreDerivedFromTheThreeLegs(self):
        by = {(v["pair"], v["null"]): v for v in self.res["verdicts"]}
        self.assertEqual(len(by), 9)
        for v in by.values():
            self.assertIn(v["verdict"], ("supported", "supported_only_conditionally",
                                         "not_supported"))
            if v["verdict"] == "supported":
                self.assertTrue(v["primary_supported"])
                self.assertEqual(v["primary_failing"], "")
            if v["verdict"] == "supported_only_conditionally":
                self.assertFalse(v["primary_supported"])
                self.assertTrue(v["neutral_supported"])

    def testTheTwoPairsThatFailTheirPrimaryTestAreRecordedAsSuch(self):
        by = {(v["pair"], v["null"]): v["verdict"] for v in self.res["verdicts"]}
        self.assertEqual(by[("influx_x_chaperone_only", "additive")],
                         "supported_only_conditionally")
        self.assertEqual(by[("influx_x_chaperone_only", "bliss")],
                         "supported_only_conditionally")
        self.assertEqual(by[("influx_x_total_capacity", "additive")], "supported")

    def testSensitivityBoundsNeverImputeTheUnresolvedBackgrounds(self):
        for r in self.res["sensitivity"]:
            self.assertLessEqual(r["usable_lo"], r["conditional"] + 1e-12)
            self.assertGreaterEqual(r["usable_hi"], r["conditional"] - 1e-12)
            self.assertEqual(r["n_unresolved"], 2)
            self.assertEqual(r["n_requested"], 60)

    def testNoCellLevelPValueIsReportedAnywhere(self):
        blob = json.dumps(self.res)
        self.assertNotIn("cell_level_p", blob)
        for r in self.res["estimates"]:
            self.assertLessEqual(r["sign_test_n_worse"] + r["sign_test_n_better"]
                                 + r["sign_test_n_tied"], 46)


# --------------------------------------------------------------------------
# the tracked bridge from the two closure documents to the gitignored run
# --------------------------------------------------------------------------

class TestClosureCheckSkipsRatherThanFails(unittest.TestCase):
    """model-free: this half must behave correctly on a clean checkout.

    `check_d_closure` is the only tracked thing that can notice the closure prose
    drifting from the run that produced it.  a validator that silently reported
    success when the evidence was missing would be worse than none at all, so the
    absent-results branch is tested explicitly rather than assumed.
    """

    def testAbsentRunRootYieldsNoneNotAnEmptyPass(self):
        original_root = closure.RUN_ROOT
        try:
            closure.RUN_ROOT = Path("/nonexistent/experiment/d/run/root")
            self.assertIsNone(closure.runChecks())
            self.assertFalse(closure.resultsPresent())
        finally:
            closure.RUN_ROOT = original_root

    def testTheHeadlineExpectationsAreTheOnesTheDocumentsClaim(self):
        """pinned here so a headline cannot be relaxed by editing the checker."""
        m = closure.HEADLINE_MAJORITY
        self.assertEqual(m[("influx_x_total_capacity", "additive")], (46, True))
        self.assertEqual(m[("influx_x_total_capacity", "multiplicative")], (44, True))
        self.assertEqual(m[("influx_x_total_capacity", "bliss")], (41, True))
        self.assertEqual(m[("nascent_x_total_capacity", "additive")], (37, True))
        self.assertEqual(m[("nascent_x_total_capacity", "multiplicative")], (38, True))
        self.assertEqual(m[("nascent_x_total_capacity", "bliss")], (37, True))
        # chaperone-only: multiplicative is headline, additive and bliss are not
        self.assertEqual(m[("influx_x_chaperone_only", "multiplicative")], (35, True))
        self.assertEqual(m[("influx_x_chaperone_only", "additive")], (28, False))
        self.assertEqual(m[("influx_x_chaperone_only", "bliss")], (28, False))
        self.assertEqual(closure.HEADLINE_COLLAPSE["overall"][0], 43)
        self.assertEqual(closure.HEADLINE_COLLAPSE["nascent_x_total_capacity"][0], 13)


@unittest.skipUnless(closure.resultsPresent(),
                     "experiment D closure output absent; the tracked documents "
                     "were not checked against a run")
class TestTrackedDocumentsMatchTheShippedOutput(unittest.TestCase):

    def testEveryClosureCheckPasses(self):
        checks = closure.runChecks()
        bad = [f"{n} [{d}]" for n, ok, d in checks if not ok]
        self.assertEqual(bad, [], f"{len(bad)} of {len(checks)} closure checks fail")

    def testTheCheckIsNotVacuous(self):
        """a validator asserting nothing would pass against anything at all."""
        self.assertGreater(len(closure.runChecks()), 150)

    def testAWrongHeadlineNumberWouldBeCaught(self):
        """negative control on the validator itself, not on the data.

        this is the test that makes the class above meaningful: it proves the
        headline comparison is live, so a silent drift between the documents and
        `D_final/` cannot pass as 177 green checks.
        """
        key = ("influx_x_chaperone_only", "additive")
        original_value = closure.HEADLINE_MAJORITY[key]
        try:
            closure.HEADLINE_MAJORITY[key] = (35, True)   # the multiplicative count
            bad = [n for n, ok, _ in closure.runChecks() if not ok]
            self.assertIn("majority influx_x_chaperone_only/additive", bad)
            self.assertIn("primary_supported influx_x_chaperone_only/additive", bad)
        finally:
            closure.HEADLINE_MAJORITY[key] = original_value
        self.assertEqual([n for n, ok, _ in closure.runChecks() if not ok], [])


if __name__ == "__main__":
    unittest.main()
