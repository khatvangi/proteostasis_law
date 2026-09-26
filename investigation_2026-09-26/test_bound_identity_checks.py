import math

from bound_identity_checks import (
    f_crit_exact,
    f_crit_log_approx,
    pool_state,
    raw_error_from_substitution,
    substitution_error_from_raw,
)


def test_reference_bound_is_not_quoted_value():
    exact = f_crit_exact(300, 0.7, 0.3, 0.3)
    approx = f_crit_log_approx(300, 0.7, 0.3, 0.3)
    assert math.isclose(exact, 0.005658142850513878, rel_tol=1e-12)
    assert math.isclose(approx, 0.0056615070466465465, rel_tol=1e-12)
    assert math.isclose(-math.log(0.7) / 300, 0.0011889164797957749)
    assert exact > 4 * (-math.log(0.7) / 300)


def test_synonymous_filter_is_applied_once_and_only_to_raw_error():
    raw = 1e-3
    measured_substitution = substitution_error_from_raw(raw, 0.3)
    assert measured_substitution == 7e-4
    assert raw_error_from_substitution(measured_substitution, 0.3) == raw
    # If 7e-4 was measured as a substitution rate, multiplying by 0.7 again
    # creates a false 4.9e-4 effective rate.
    assert measured_substitution * 0.7 == 4.9e-4


def test_total_pool_scaling_is_material_to_feedback():
    low = pool_state(300, 0.1)
    high = pool_state(4000, 0.1)
    assert high.c_free_uM < low.c_free_uM
    assert high.phi > low.phi * 40
