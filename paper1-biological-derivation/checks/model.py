"""
paper 1 biological derivation, stages 0-4: the one flux table.

every equation of stages 1, 2 and 4 is ASSEMBLED from this table
(stoichiometry x flux, then dilution -mu*x on every pool of the stage). no
right-hand side is typed by hand, so the object that is checked for
conservation is the object whose steady states are analysed.

a flux is a transfer of monomer-equivalents (client) or of machine molecules
(uM/s) that either moves material between two named pools, or crosses the
boundary (synthesis in; proteolytic degradation out). dilution -mu*x is a
concentration sink: molecules go to daughter cells, they are not destroyed.

stage flux lists
  S1   stage 1: first-order (unlimited) machinery
  S2   stage 2: finite chaperone, protease and disaggregase cycles
  S2C  stage 2 with protease and disaggregase at their first-order limit.
       this is term-for-term the rebuild core WO-02 (= WO-04 V0)
  S4   stage 4: S2 + reversible chaperone binding to aggregate (C + A <-> CA)
  S4C  S2C + the same binding step. term-for-term WO-04 V1

no parameter value in this directory is a claim about E. coli. numerical
checks draw parameters from declared mathematical test domains.
"""
import sys

import numpy as np
import sympy as sp

sys.dont_write_bytecode = True

# ---------------------------------------------------------------- symbols
STATE_NAMES = ["N", "U", "B", "C", "DU", "D", "A", "EA", "E", "CA"]
XS = {n: sp.Symbol(n, nonnegative=True) for n in STATE_NAMES}
N, U, B, C, DU, D, A, EA, E, CA = (XS[n] for n in STATE_NAMES)

PARAM_NAMES = ["s_P", "eps", "mu", "k_mis", "k_a", "k_dA",
               "k_r", "k_d", "k_dis",
               "s_C", "k_on", "k_off", "k_cat", "phi",
               "s_D", "k_onD", "k_offD", "k_catD",
               "s_E", "k_onE", "k_offE", "k_catE",
               "k_onA", "k_offA"]
PS = {n: sp.Symbol(n, nonnegative=True) for n in PARAM_NAMES}
(s_P, eps, mu, k_mis, k_a, k_dA, k_r, k_d, k_dis, s_C, k_on, k_off, k_cat, phi,
 s_D, k_onD, k_offD, k_catD, s_E, k_onE, k_offE, k_catE, k_onA, k_offA) = (
    PS[n] for n in PARAM_NAMES)

# ---------------------------------------------------------------- flux table
# name -> (rate, stoichiometry {pool: coefficient}, stage introduced)
# an empty receiver side means the material leaves across the boundary.
FLUX = {
    # stage 1 -- present in every stage
    "syn_N":  ((1 - eps) * s_P,        {"N": 1},                    1),
    "syn_U":  (eps * s_P,              {"U": 1},                    1),
    "unfold": (k_mis * N,              {"N": -1, "U": 1},           1),
    "agg":    (k_a * U**2,             {"U": -1, "A": 1},           1),
    "degA":   (k_dA * A,               {"A": -1},                   1),
    # stage 1 first-order machinery (replaced at stage 2)
    "resc":   (k_r * U,                {"U": -1, "N": 1},           1),
    "degU":   (k_d * U,                {"U": -1},                   1),
    "dis":    (k_dis * A,              {"A": -1, "U": 1},           1),
    # stage 2 chaperone cycle (DnaK/DnaJ/GrpE or GroEL/ES as one effective cycle)
    "syn_C":  (s_C,                    {"C": 1},                    2),
    "bind":   (k_on * C * U,           {"C": -1, "U": -1, "B": 1},  2),
    "rel":    (k_off * B,              {"B": -1, "C": 1, "U": 1},   2),
    "cycN":   (phi * k_cat * B,        {"B": -1, "C": 1, "N": 1},   2),
    "cycU":   ((1 - phi) * k_cat * B,  {"B": -1, "C": 1, "U": 1},   2),
    # stage 2 protease cycle (Lon / ClpXP acting on non-native monomer)
    "syn_D":  (s_D,                    {"D": 1},                    2),
    "bindD":  (k_onD * D * U,          {"D": -1, "U": -1, "DU": 1}, 2),
    "relD":   (k_offD * DU,            {"DU": -1, "D": 1, "U": 1},  2),
    "catD":   (k_catD * DU,            {"DU": -1, "D": 1},          2),
    # stage 2 disaggregase cycle (ClpB-like, counterfactual independent pool)
    "syn_E":  (s_E,                    {"E": 1},                    2),
    "bindE":  (k_onE * E * A,          {"E": -1, "A": -1, "EA": 1}, 2),
    "relE":   (k_offE * EA,            {"EA": -1, "E": 1, "A": 1},  2),
    "catE":   (k_catE * EA,            {"EA": -1, "E": 1, "U": 1},  2),
    # stage 4 feedback A: chaperone binds aggregate surface
    "bindA":  (k_onA * C * A,          {"C": -1, "A": -1, "CA": 1}, 4),
    "relA":   (k_offA * CA,            {"CA": -1, "C": 1, "A": 1},  4),
}

CORE = ["syn_N", "syn_U", "unfold", "agg", "degA"]
CHAP = ["syn_C", "bind", "rel", "cycN", "cycU"]
PROT = ["syn_D", "bindD", "relD", "catD"]
DISG = ["syn_E", "bindE", "relE", "catE"]
FEEDA = ["bindA", "relA"]
STAGES = {
    "S1":  CORE + ["resc", "degU", "dis"],
    "S2":  CORE + CHAP + PROT + DISG,
    "S2C": CORE + ["degU", "dis"] + CHAP,
    "S4":  CORE + CHAP + PROT + DISG + FEEDA,
    "S4C": CORE + ["degU", "dis"] + CHAP + FEEDA,
}

# conserved-quantity weights. a complex carries one client monomer-equivalent
# AND one machine molecule, so it counts in both totals.
WEIGHTS = {
    "client":       {"N": 1, "U": 1, "B": 1, "DU": 1, "A": 1, "EA": 1, "CA": 1},
    "chaperone":    {"C": 1, "B": 1, "CA": 1},
    "protease":     {"D": 1, "DU": 1},
    "disaggregase": {"E": 1, "EA": 1},
}
# boundary role of each flux for each total: +1 synthesis, -1 degradation, 0 internal
BOUNDARY = {
    "client":       {"syn_N": 1, "syn_U": 1, "degA": -1, "degU": -1, "catD": -1},
    "chaperone":    {"syn_C": 1},
    "protease":     {"syn_D": 1},
    "disaggregase": {"syn_E": 1},
}

# ---------------------------------------------------------------- dimensions
# (concentration exponent, time exponent). every rhs term must be (1, -1).
CONC = {"s_P": (1, -1), "s_C": (1, -1), "s_D": (1, -1), "s_E": (1, -1),
        "eps": (0, 0), "phi": (0, 0),
        "k_a": (-1, -1), "k_on": (-1, -1), "k_onD": (-1, -1), "k_onE": (-1, -1),
        "k_onA": (-1, -1)}
DIMS = {n: (1, 0) for n in STATE_NAMES}
DIMS.update({n: CONC.get(n, (0, -1)) for n in PARAM_NAMES})


def stage_states(stage):
    """pools touched by the fluxes of a stage, in STATE_NAMES order."""
    used = {p for f in STAGES[stage] for p in FLUX[f][1]}
    return [n for n in STATE_NAMES if n in used]


def assemble(stage):
    """rhs dict {pool: expression}: sum of stoichiometry x rate, minus dilution."""
    names = stage_states(stage)
    rhs = {n: sp.Integer(0) for n in names}
    for f in STAGES[stage]:
        rate, st, _ = FLUX[f]
        for pool, c in st.items():
            rhs[pool] += c * rate
    for n in names:
        rhs[n] += -mu * XS[n]
    return rhs


def total(kind, stage):
    return sum(w * XS[p] for p, w in WEIGHTS[kind].items() if p in stage_states(stage))


def expected_total_rate(kind, stage):
    """source - boundary sinks - mu * total, written from the BOUNDARY table,
    independently of the stoichiometry."""
    out = -mu * total(kind, stage)
    for f, sign in BOUNDARY[kind].items():
        if f in STAGES[stage]:
            out += sign * FLUX[f][0]
    return out


# ---------------------------------------------------------------- numerics
def lambdas(stage):
    """numeric rhs and jacobian of a stage, state order stage_states(stage)."""
    names = stage_states(stage)
    xs = [XS[n] for n in names]
    ps = [PS[n] for n in PARAM_NAMES]
    rhs = assemble(stage)
    vec = sp.Matrix([rhs[n] for n in names])
    f = sp.lambdify((xs, ps), list(vec), "numpy")
    J = sp.lambdify((xs, ps), vec.jacobian(xs), "numpy")
    return names, f, J


def pvec(p):
    return [float(p.get(n, 0.0)) for n in PARAM_NAMES]


def lu(rng, lo, hi):
    return float(np.exp(rng.uniform(np.log(lo), np.log(hi))))


def sample_stage2(rng):
    """declared MATHEMATICAL test domain for the finite-machinery stages. the
    chaperone and client ranges are WO-04's (bifurcation.py:sample); protease
    and disaggregase ranges are declared here on the same pattern. machine
    totals are drawn and s_X = mu X_T is set, so X_T is the drawn value."""
    p = {"eps": lu(rng, 1e-4, 0.5), "mu": lu(rng, 1e-5, 1e-3),
         "k_mis": lu(rng, 1e-7, 1e-4), "k_a": lu(rng, 1e-6, 1e-1),
         "k_dA": lu(rng, 1e-7, 1e-3), "phi": rng.uniform(0.05, 1.0),
         "k_on": lu(rng, 1e-2, 10), "k_off": lu(rng, 1e-2, 10), "k_cat": lu(rng, 1e-3, 1),
         "k_onD": lu(rng, 1e-2, 10), "k_offD": lu(rng, 1e-2, 10), "k_catD": lu(rng, 1e-3, 1),
         "k_onE": lu(rng, 1e-2, 10), "k_offE": lu(rng, 1e-2, 10), "k_catE": lu(rng, 1e-3, 1),
         "P_T": lu(rng, 300, 5000), "C_T": lu(rng, 1, 100),
         "D_T": lu(rng, 0.1, 100), "E_T": lu(rng, 0.01, 10)}
    p["s_P"] = p["mu"] * p["P_T"]
    for m in "CDE":
        p["s_" + m] = p["mu"] * p[m + "_T"]
    return p
