"""
WO-02 conservative source-sink proteostasis model (concentrations in uM, s).

states  N native, U free non-native, B chaperone-client complex,
        A aggregate (monomer-equivalents), C free chaperone
totals  P_T = N + U + B + A,  C_T = C + B

  dN/dt = (1-eps) s_P + phi k_cat B - k_mis N                         - mu N
  dU/dt =  eps s_P + k_mis N + ((1-phi) k_cat + k_off) B - k_on C U
           - k_d U - k_a U^2 + k_dis A                                 - mu U
  dB/dt =  k_on C U - (k_off + k_cat) B                                - mu B
  dC/dt =  s_C + (k_off + k_cat) B - k_on C U                          - mu C
  dA/dt =  k_a U^2 - (k_dis + k_dA) A                                  - mu A

every internal flux appears with opposite signs in exactly two equations; the
only sources are s_P and s_C, the only sinks degradation (k_d U, k_dA A) and
dilution (mu x). there is no state-dependent source: nothing like the legacy
Phi factor can be written in this form.

the symbolic system is defined once (RHS_SYM) and lambdified for numerics, so
the object that is proved conservative is the object that is integrated.
"""
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "WO-01"))
import units as un  # noqa: E402
import variables as var  # noqa: E402

STATE_NAMES = ["N", "U", "B", "A", "C"]
PARAM_NAMES = ["s_P", "eps", "s_C", "mu", "k_on", "k_off", "k_cat", "phi",
               "k_d", "k_a", "k_dis", "k_dA", "k_mis"]

X = sp.symbols(STATE_NAMES, nonnegative=True)
N, U, B, A, C = X
PAR = sp.symbols(PARAM_NAMES, nonnegative=True)
(s_P, eps, s_C, mu, k_on, k_off, k_cat, phi, k_d, k_a, k_dis, k_dA, k_mis) = PAR

# named fluxes (uM/s). listing them separately makes the donor/receiver
# structure inspectable by tests.
FLUX = {
    "synth_native": (1 - eps) * s_P,
    "synth_nonnative": eps * s_P,
    "synth_chaperone": s_C,
    "bind": k_on * C * U,
    "release": k_off * B,
    "cycle_native": phi * k_cat * B,
    "cycle_nonnative": (1 - phi) * k_cat * B,
    "unfold": k_mis * N,
    "degrade_U": k_d * U,
    "aggregate": k_a * U**2,
    "disaggregate": k_dis * A,
    "degrade_A": k_dA * A,
}

# stoichiometry: flux -> {state: coefficient}. dilution is added separately.
STOICH = {
    "synth_native": {"N": 1},
    "synth_nonnative": {"U": 1},
    "synth_chaperone": {"C": 1},
    "bind": {"C": -1, "U": -1, "B": 1},
    "release": {"B": -1, "U": 1, "C": 1},
    "cycle_native": {"B": -1, "N": 1, "C": 1},
    "cycle_nonnative": {"B": -1, "U": 1, "C": 1},
    "unfold": {"N": -1, "U": 1},
    "degrade_U": {"U": -1},
    "aggregate": {"U": -1, "A": 1},
    "disaggregate": {"A": -1, "U": 1},
    "degrade_A": {"A": -1},
}

# a protein molecule's client moiety is counted in P_T; the chaperone moiety
# of B is counted in C_T. these weights define the two conserved totals.
PROTEIN_W = {"N": 1, "U": 1, "B": 1, "A": 1, "C": 0}
CHAP_W = {"N": 0, "U": 0, "B": 1, "A": 0, "C": 1}


def build_rhs():
    rhs = {s: sp.Integer(0) for s in STATE_NAMES}
    for f, st in STOICH.items():
        for s, c in st.items():
            rhs[s] += c * FLUX[f]
    for s, x in zip(STATE_NAMES, X):
        rhs[s] += -mu * x
    return [sp.simplify(rhs[s]) for s in STATE_NAMES]


RHS_SYM = build_rhs()
P_T = N + U + B + A
C_T = C + B

_rhs_num = sp.lambdify((X, PAR), RHS_SYM, "numpy")
_jac_num = sp.lambdify((X, PAR), sp.Matrix(RHS_SYM).jacobian(X), "numpy")


def rhs(t, x, p):
    """numeric RHS; p is a sequence in PARAM_NAMES order."""
    return np.asarray(_rhs_num(x, p), dtype=float)


def jac(t, x, p):
    return np.asarray(_jac_num(x, p), dtype=float)


def pvec(d):
    """dict -> parameter vector in PARAM_NAMES order (missing keys raise)."""
    return np.array([float(d[k]) for k in PARAM_NAMES])


# illustrative scenario for later WOs. NOT a calibrated parameter set: every
# value is UNVERIFIED until WO-06 audits it. legacy-derived values are marked.
LN2 = float(np.log(2))
SCENARIO = {
    "mu": LN2 / 3600.0,        # 60 min doubling (legacy T_gen)
    "P_T": 3000.0,             # uM, order of magnitude only; legacy used 300
    "C_T": 50.0,               # uM, legacy C_tot
    "k_on": 1.0,               # 1/(uM s); with k_off=1 gives legacy K_d = 1 uM
    "k_off": 1.0,              # 1/s
    "k_cat": 1e-2,             # 1/s, legacy k_obs_max
    "phi": 1.0,                # legacy has no unproductive-cycle partition
    "k_d": 3e-4,               # 1/s, legacy k_deg
    "k_a": 1e-3,               # 1/(uM s) = legacy k_agg 1e3 /(M s)
    "k_dis": 4e-4,             # 1/s, legacy k_clear read as disaggregation
    "k_dA": 0.0,
    "k_mis": 0.0,
    "eps": 0.04,               # placeholder; WO-05 derives it
}


def scenario_params(**over):
    d = {**SCENARIO, **over}
    d["s_P"] = d["mu"] * d["P_T"]     # balanced growth, turnover neglected
    d["s_C"] = d["mu"] * d["C_T"]     # chaperone pool at its balanced size
    return d
