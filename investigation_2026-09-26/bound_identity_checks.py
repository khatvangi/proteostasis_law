"""Independent numerical checks for the second-wave proteostasis audit.

This file is deliberately standalone: it does not import or modify either
upstream model.  It checks the arithmetic identity, the measured-vs-raw error
mapping, and the effect of treating 300 uM as an accessible pool versus the
~4 mM total E. coli protein concentration anchor.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


def f_crit_exact(N: float, p_correct: float, p_misfold: float,
                 s_syn: float) -> float:
    """Raw per-codon error threshold when s_syn is a separate filter."""
    return (1.0 - p_correct ** (1.0 / N)) / ((1.0 - s_syn) * p_misfold)


def f_crit_log_approx(N: float, p_correct: float, p_misfold: float,
                      s_syn: float) -> float:
    return -math.log(p_correct) / (N * (1.0 - s_syn) * p_misfold)


def substitution_error_from_raw(raw_error: float, s_syn: float) -> float:
    """Convert a raw decoding error into a nonsynonymous substitution rate."""
    return raw_error * (1.0 - s_syn)


def raw_error_from_substitution(substitution_error: float, s_syn: float) -> float:
    """Recover raw error only when the measured value is substitution-level."""
    return substitution_error / (1.0 - s_syn)


@dataclass(frozen=True)
class PoolState:
    prot_tot_uM: float
    P: float
    c_free_uM: float
    v_fold_s: float
    v_agg_s: float
    phi: float


def pool_state(prot_tot_uM: float, P: float, *, c_tot_uM: float = 50.0,
               kd_uM: float = 1.0, k_obs_max_s: float = 1e-2,
               k_agg_M_s: float = 1e3, k_deg_s: float = 3e-4) -> PoolState:
    """Evaluate only the upstream ODE's algebraic feedback terms."""
    M_uM = P * prot_tot_uM
    c_free = c_tot_uM / (1.0 + M_uM / kd_uM)
    v_fold = k_obs_max_s * c_free / (c_free + kd_uM)
    v_agg = k_agg_M_s * M_uM * 1e-6
    phi = 1.0 + v_agg / (v_fold + k_deg_s)
    return PoolState(prot_tot_uM, P, c_free, v_fold, v_agg, phi)


def main() -> None:
    N, pc, pm, S = 300.0, 0.70, 0.30, 0.30
    exact = f_crit_exact(N, pc, pm, S)
    approx = f_crit_log_approx(N, pc, pm, S)
    print(f"exact_raw_fcrit={exact:.15g}")
    print(f"log_approx_raw_fcrit={approx:.15g}")
    print(f"quoted_no_filters={-math.log(pc)/N:.15g}")
    print(f"exact_substitution_fcrit={f_crit_exact(N, pc, pm, 0.0):.15g}")
    for P in (0.001, 0.01, 0.1, 0.25):
        for total in (300.0, 4000.0):
            s = pool_state(total, P)
            print(f"pool total={total:g} P={P:g} M_uM={total*P:g} "
                  f"c_free_uM={s.c_free_uM:.9g} phi={s.phi:.9g}")


if __name__ == "__main__":
    main()
