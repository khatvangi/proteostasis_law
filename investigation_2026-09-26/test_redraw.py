import numpy as np
from redraw import (
    approximate_free,
    equilibrium_polynomial,
    finite_pool,
    folding_rate,
    g,
    gprime,
    positive_equilibria,
    positive_zero,
    stationary_points,
)


def test_finite_pool_mass_balance_and_equilibrium():
    cf, cb = finite_pool(50.0)
    mf = 50.0 - float(cb)
    assert np.isclose(cf + cb, 50.0)
    assert np.isclose(mf + cb, 50.0)
    assert np.isclose(float(cf) * mf / float(cb), 1.0)


def test_approximation_is_not_exact_finite_pool():
    cf, _ = finite_pool(50.0)
    assert float(cf) > float(approximate_free(50.0))


def test_scalar_stationary_and_threshold():
    xs = stationary_points()
    assert len(xs) == 1
    assert np.isclose(xs[0], 3.890758215117349, rtol=1e-10)
    assert np.isclose(g(xs[0]), 4.80218919587308, rtol=1e-10)
    assert np.isclose(positive_zero(), 9.264593793858728, rtol=1e-10)
    assert np.isclose(g(positive_zero()), 0.0, atol=1e-10)
    assert gprime(1.0) > 0 and gprime(5.0) < 0


def test_equilibrium_polynomial_has_correct_lambda_coefficient():
    assert np.allclose(equilibrium_polynomial(2.0), [-0.15, 0.85, 3.0, -2.0])
    assert not np.allclose(equilibrium_polynomial(2.0), [-0.15, 0.85, 5.0, -2.0])


def test_lambda_two_roots_match_original_g_and_stability():
    roots = positive_equilibria(2.0)
    assert np.allclose(roots, [0.5808674541, 7.9669676044], rtol=1e-9, atol=1e-10)
    for root in roots:
        assert abs(g(root) - 2.0) < 1e-9
    assert gprime(roots[0]) > 0
    assert gprime(roots[1]) < 0


def test_all_real_lambda_two_polynomial_roots_match_original_g():
    roots = np.roots(equilibrium_polynomial(2.0))
    assert np.all(np.abs(roots.imag) < 1e-12)
    for root in roots.real:
        assert abs(g(root) - 2.0) < 1e-9
    assert np.isclose(sorted(roots.real), [-2.8811683919, 0.5808674541, 7.9669676044],
                      rtol=1e-9, atol=1e-10).all()


def test_finite_pool_rate_ratios_are_mass_balance_consistent():
    for M, expected_ratio in [(0.0, 1.0), (10.0, 1.19042), (50.0, 1.75382), (300.0, 1.16534)]:
        cf_exact, cb = finite_pool(M)
        ratio = float(folding_rate(cf_exact) / folding_rate(approximate_free(M)))
        assert np.isclose(ratio, expected_ratio, rtol=0.0, atol=5e-6)
        mf = M - float(cb)
        assert np.isclose(float(cf_exact) * mf, float(cb), rtol=0.0, atol=1e-10)
