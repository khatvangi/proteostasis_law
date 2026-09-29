import sys
import unittest
from pathlib import Path

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import error_semantics as es  # noqa: E402
import run_wo05 as r  # noqa: E402

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
        self.assertLess(g["n_datasets_per_codon_max"], 80)

    def test_g55_blocks_unresolved_empirical_denominator(self):
        g = r.g55()
        self.assertEqual(g["status"], "BLOCKED")
        self.assertFalse(g["computed_headroom"])
        self.assertIn("eligible-dataset denominator", g["reason"])


if __name__ == "__main__":
    unittest.main()
