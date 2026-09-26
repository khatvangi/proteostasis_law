import numpy as np
from redraw import finite_pool, approximate_free, g, gprime, stationary_points


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
    assert np.isclose(g(9.264593793858728), 0.0, atol=1e-10)
    assert gprime(1.0) > 0 and gprime(5.0) < 0
