import math
import sys
import unittest
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import error_semantics as es  # noqa: E402
import run_wo05 as r  # noqa: E402
import landerer_s4 as L  # noqa: E402

KW = dict(N_prot=300.0, p_misfold=0.3, T_gen=1800.0, S=0.3)


class TestTypes(unittest.TestCase):
    def test_substitution_through_one_minus_S_fails(self):  # G5.3
        for kind in (es.SUB, es.MS):
            with self.assertRaises(TypeError):
                es.substitution_from_raw(es.ErrorRate(1e-3, kind), 0.3)
            with self.assertRaises(TypeError):
                es.flux_checked_legacy(es.ErrorRate(1e-3, kind), **KW)

    def test_raw_discounted_exactly_once(self):
        raw = es.ErrorRate(1e-3, es.RAW)
        self.assertAlmostEqual(es.flux(raw, **KW), es.flux_legacy(1e-3, **KW), places=18)

    def test_ms_not_discounted(self):
        ms = es.ErrorRate(1e-3, es.MS)
        self.assertAlmostEqual(es.flux(ms, **KW) * (1 - KW["S"]), es.flux_legacy(1e-3, **KW),
                               places=18)

    def test_detectability_only_raises_rate(self):
        ms = es.ErrorRate(1e-3, es.MS)
        self.assertGreater(es.substitution_from_ms(ms, 0.5).value, ms.value)
        with self.assertRaises(ValueError):
            es.substitution_from_ms(ms, 1.5)

    def test_bad_kind_rejected(self):
        with self.assertRaises(ValueError):
            es.ErrorRate(1e-3, "per_codon")


class TestNumbers(unittest.TestCase):
    def test_g54_threshold(self):
        g = r.g54()
        self.assertAlmostEqual(g["threshold_raw_exact"], 5.6581428505e-3, places=12)
        self.assertAlmostEqual(g["legacy_quoted_minus_lnP_over_N"], 1.1889164798e-3, places=12)
        self.assertAlmostEqual(g["threshold_sub_exact"], 0.7 * g["threshold_raw_exact"], places=15)

    def test_g52_anchors(self):
        rows = r.g52()
        self.assertTrue(all(x["anchor_found"] for x in rows))
        dd = {(Path(x["file"]).name, x["line"]) for x in rows if x["class"] == "double-discount"}
        self.assertEqual(dd, {("11_headroom_sensitivity.py", 81),
                              ("12_chaperone_availability.py", 72),
                              ("09_supraadditivity.py", 99)})

    def test_g51_conditional_mean(self):
        g = r.g51()
        self.assertTrue(g["legacy_mu_equals_DataS2_mean"])
        self.assertTrue(g["legacy_data_S2_identical_to_supplement"])
        self.assertLess(g["n_datasets_per_codon_max"], 80)
        self.assertEqual(g["total_substitution_psms"], 58990)   # Landerer text: 58,990
        self.assertEqual(g["numerator_mismatches_codon_counts_vs_substitution_rows"], 0)


class TestLandererS4(unittest.TestCase):
    """the aggregation question is answered from per-dataset counts, not (sd/se)^2."""

    @classmethod
    def setUpClass(cls):
        cls.cc, cls.se = L.load_s4("ecoli")
        cls.s2 = L.load_s2("ecoli")

    def test_s2_is_detected_only_mean(self):
        rec = L.reconstruct_s2(self.cc, self.s2)
        det = rec["detected_only_error_count_gt_0"]
        self.assertLess(det["max_abs_dev_mean"], 1e-8)
        self.assertLess(det["max_abs_dev_sd"], 1e-8)
        self.assertTrue(det["n_equals_sd_over_se_sq"])
        # the zero-inclusive rule, which the user-proposed reading implies, does NOT match
        zi = rec["zero_inclusive_covering"]
        self.assertGreater(zi["max_abs_dev_mean"], 1e-3)
        self.assertFalse(zi["n_equals_sd_over_se_sq"])

    def test_zero_cells_are_covered_not_missing(self):
        c = L.cell_census(self.cc, 80)
        self.assertEqual(c["datasets"], 80)
        self.assertEqual(c["covered_zero_cells"], 3015)
        self.assertEqual(c["no_coverage_cells"], 19)
        self.assertLess(c["max_abs_rate_minus_err_over_base"], 1e-15)

    def test_estimator_ordering(self):
        w = L.usage_weights()
        rt = L.per_codon_rates(self.cc)
        cond = L.usage_weighted(rt.conditional_mean, w)
        zero = L.usage_weighted(rt.zero_inclusive_mean, w)
        pool = L.usage_weighted(rt.pooled_ratio, w)
        self.assertAlmostEqual(cond, 6.334247974475959e-4, delta=1e-10)   # legacy mubar
        self.assertGreater(cond, zero)          # selection on detection inflates
        self.assertGreater(zero, pool)
        self.assertTrue((rt.conditional_mean >= rt.zero_inclusive_mean - 1e-15).all())


class TestG55(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.g = r.g55()

    def test_legacy_x25_reproduced(self):
        self.assertAlmostEqual(self.g["old_x"], 24.81728310207098, places=6)

    def test_corrected_flux_has_no_one_minus_S(self):
        f = self.g["inputs"]["conditional_mean_DataS2"]
        J_leg = f["headroom"]["legacy_(1-S)_1/T"]["J"]
        J_fix = f["headroom"]["no_(1-S)_1/T"]["J"]
        self.assertAlmostEqual(J_fix * 0.7, J_leg, delta=1e-18)
        J_both = f["headroom"]["no_(1-S)_ln2/T"]["J"]
        self.assertAlmostEqual(J_both, J_fix * math.log(2), delta=1e-18)

    def test_each_correction_direction(self):
        e = self.g["effect_on_old_x"]
        self.assertLess(e["remove_(1-S)_only"], self.g["old_x"])     # more inflow, less room
        self.assertGreater(e["ln2_only"], self.g["old_x"])           # less inflow, more room
        self.assertAlmostEqual(e["remove_(1-S)_only"], 17.3204, places=3)
        self.assertAlmostEqual(e["ln2_only"], 35.8765, places=3)
        self.assertAlmostEqual(e["both"], 25.0643, places=3)
        self.assertAlmostEqual(e["zero_inclusive_input_both"], 63.104, places=2)
        self.assertAlmostEqual(e["pooled_input_both"], 112.873, places=2)

    def test_stationary_not_mixed(self):
        st = self.g["stikeleather_2026"]
        self.assertIn("STATIONARY", st["condition"])
        self.assertTrue(st["headroom"].startswith("NOT COMPUTED"))
        self.assertNotIn("stikeleather", " ".join(self.g["inputs"]).lower())

    def test_label_is_not_a_baseline(self):
        self.assertIn("NOT an E. coli physiological", self.g["quantity_label"])
        self.assertIn("NOT raw", self.g["quantity_label"])

if __name__ == "__main__":
    unittest.main()
