import unittest
import numpy as np
from model import Params, rhs, jacobian, finite_difference_jacobian, equilibria, classify_equilibrium, integrate


class ModelTests(unittest.TestCase):
    def test_analytic_numerical_jacobian_agree(self):
        rng = np.random.default_rng(8128)
        for _ in range(40):
            y = np.exp(rng.uniform(-3, 2, 2))
            np.testing.assert_allclose(jacobian(y, Params()), finite_difference_jacobian(y, Params()), rtol=2e-6, atol=2e-7)

    def test_positive_quadrant_is_invariant_on_boundaries(self):
        p = Params()
        for a in np.geomspace(1e-8, 1e3, 60):
            self.assertGreaterEqual(rhs(0, [0, a], p)[0], 0.0)
        for u in np.geomspace(1e-8, 1e3, 60):
            self.assertGreaterEqual(rhs(0, [u, 0], p)[1], 0.0)

    def test_numerical_positivity(self):
        p = Params()
        for y0 in ([0, 0], [1e-10, 4], [4, 1e-10], [2, 2]):
            sol, ok, _ = integrate(p, y0, 200, 500)
            self.assertTrue(ok, sol.message)
            self.assertGreaterEqual(np.min(sol.y), -1e-8)

    def test_equilibrium_residual_and_stability(self):
        roots = equilibria(Params())
        self.assertTrue(roots)
        for root in roots:
            self.assertLess(np.linalg.norm(rhs(0, root, Params()), np.inf), 1e-8)
            label, eig = classify_equilibrium(root, Params())
            self.assertIn(label, {"stable", "unstable", "nonhyperbolic"})
            self.assertTrue(np.all(np.isfinite(eig)))

    def test_nascent_occupancy_is_not_influx(self):
        low, high = Params(N=0.0), Params(N=3.0)
        self.assertAlmostEqual(rhs(0, [0, 0], low)[0], low.J)
        self.assertAlmostEqual(rhs(0, [0, 0], high)[0], high.J)
        self.assertLess(rhs(0, [1, 0], high)[0], 0.5 + rhs(0, [1, 0], low)[0])


if __name__ == "__main__":
    unittest.main()

