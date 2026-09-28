import sys
import unittest
from pathlib import Path

import sympy as sp

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parent))
import units as u  # noqa: E402
import variables as v  # noqa: E402
import legacy_units_audit as la  # noqa: E402


class TestChecker(unittest.TestCase):
    def test_rejects_scale_mix(self):
        x, y = sp.symbols("x y")
        with self.assertRaises(u.DimensionError):
            u.unit_of(x + y, {"x": u.UM, "y": u.M})

    def test_rejects_dimension_mix(self):
        x, y = sp.symbols("x y")
        with self.assertRaises(u.DimensionError):
            u.unit_of(x + y, {"x": u.UM, "y": u.PER_S})

    def test_rejects_exp_of_dimensional(self):
        x = sp.symbols("x")
        with self.assertRaises(u.DimensionError):
            u.unit_of(sp.exp(x), {"x": u.S})

    def test_undeclared_symbol_raises(self):
        with self.assertRaises(u.DimensionError):
            u.unit_of(sp.Symbol("ghost"), {})

    def test_mass_action_rate(self):
        k, a, b = sp.symbols("k a b")
        got = u.unit_of(k * a * b, {"k": u.PER_UM_PER_S, "a": u.UM, "b": u.UM})
        self.assertTrue(got.same(u.UM_PER_S))


class TestRebuildDeclarations(unittest.TestCase):
    def test_every_quantity_has_a_unit(self):  # G1.1
        for name, (unit, meaning) in {**v.STATES, **v.PARAMS}.items():
            self.assertIsInstance(unit, u.Unit, name)
            self.assertTrue(meaning, name)

    def test_no_fraction_states(self):  # G1.1
        for name, (unit, _) in v.STATES.items():
            self.assertTrue(unit.same(u.UM), f"state {name} is not a concentration")

    def test_conservation_laws_listed(self):  # G1.3
        keys = " ".join(v.CONSERVATION)
        self.assertIn("P_T", keys)
        self.assertIn("C_T", keys)


class TestLegacyAudit(unittest.TestCase):
    def setUp(self):
        self.res = la.run()

    def test_all_legacy_pieces_dimensionally_consistent(self):  # G1.2
        self.assertTrue(self.res["all_dimension_ok"], self.res["pieces"])

    def test_negative_control(self):
        self.assertTrue(self.res["negative_control_rejected"])

    def test_every_flux_classified(self):  # G1.4
        for f in self.res["fluxes"]:
            self.assertIn(f["conservative"], (True, False))
            self.assertTrue(f["note"])

    def test_phi_excess_is_nonconservative(self):
        phi = [f for f in self.res["fluxes"] if f["flux"] == "J_bare*(Phi-1)"][0]
        self.assertFalse(phi["conservative"])
        self.assertIsNone(phi["donor"])

    def test_ln2_factor(self):
        r = self.res["synthesis_rate_semantics"]["legacy_over_balanced"]
        self.assertAlmostEqual(r, 1 / 0.6931471805599453, places=12)


if __name__ == "__main__":
    unittest.main()
