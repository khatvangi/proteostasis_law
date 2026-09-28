"""
WO-03 chaperone allocation: three layers kept deliberately separate.

1. EQUILIBRIUM BENCHMARK. closed C + M <-> CM with K_d = k_off/k_on. exact
   finite-pool occupancy from the quadratic, versus the legacy closure
   C_free = C_T/(1 + M_T/K_d) (two_pool_ode.py:129). this is what a titration
   of chaperone against peptide measures. it is not in vivo capacity.

2. ATP-DRIVEN CYCLE. (a) the minimal cycle of WO-02, C + U <-> B -> C + (N|U),
   whose steady-state occupancy is set by K_M = (k_off + k_cat + mu)/k_on, not by
   K_d. (b) a four-state DnaK-like cycle (ATP/ADP free and bound, hydrolysis,
   nucleotide exchange) solved exactly as a Markov chain at fixed client
   concentration, to show that a driven cycle's effective affinity is not the
   equilibrium affinity of either nucleotide state.

3. NASCENT-CHAIN COMPETITION. WO-02 extended with nascent clients X that bind
   the same chaperone pool. the legacy's theta (fraction of chaperone committed
   elsewhere, 12_chaperone_availability.py:67) becomes an OUTPUT, BX/C_T.

all numbers used below are illustrative; none is a measured E. coli value.
"""
import sys
from pathlib import Path

import numpy as np
import sympy as sp

sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "WO-02"))
import model as m  # noqa: E402


# ---------------------------------------------------------------- layer 1
def bound_exact(C_T, M_T, K):
    """physical root of B^2 - (C_T + M_T + K) B + C_T M_T = 0, written in the
    cancellation-free form. valid for K = K_d (equilibrium) or K = K_M (tQSSA)."""
    s = C_T + M_T + K
    return 2.0 * C_T * M_T / (s + np.sqrt(s * s - 4.0 * C_T * M_T))


def free_exact(C_T, M_T, K):
    b = bound_exact(C_T, M_T, K)
    return C_T - b, M_T - b, b


def free_legacy(C_T, M_T, K_d):
    """two_pool_ode.py:129 -- treats total M as if it were free ligand."""
    return C_T / (1.0 + M_T / K_d)


def mass_balance_residuals(C_T, M_T, K):
    Cf, Mf, b = free_exact(C_T, M_T, K)
    return {"chaperone": np.abs(Cf + b - C_T),
            "client": np.abs(Mf + b - M_T),
            "mass_action": np.abs(Cf * Mf - K * b) / np.maximum(K * b, 1e-300)}


# ---------------------------------------------------------------- layer 2a
def K_M(k_on, k_off, k_cat, mu=0.0):
    return (k_off + k_cat + mu) / k_on


def cycle_closed_rhs(t, y, k_on, k_off, k_cat):
    """closed cycle, phi = 0 (client always returned to U): C + U <-> B -> C + U.
    conserved: C + B = C_T, U + B = M_T. steady state is a NON-equilibrium
    steady state whenever k_cat > 0."""
    C, U, B = y
    v = k_on * C * U - (k_off + k_cat) * B
    return [-v, -v, v]


# ---------------------------------------------------------------- layer 2b
def four_state_occupancy(U, k_onT, k_offT, k_onD, k_offD, k_h, k_h0, k_ex):
    """stationary distribution of one chaperone over (CT, CD, BT, BD) at fixed
    free client U. returns fraction bound (BT + BD).

    CT + U <-> BT  (k_onT, k_offT)     ATP state: fast, weak
    CD + U <-> BD  (k_onD, k_offD)     ADP state: slow, tight
    BT -> BD  k_h   (client/J-protein-stimulated hydrolysis)
    CT -> CD  k_h0  (basal hydrolysis)
    BD -> BT, CD -> CT  k_ex  (nucleotide exchange, ADP -> ATP)
    with k_h = k_h0 = k_ex = 0 the ATP and ADP populations decouple; the driven
    cycle is the case k_h > k_h0."""
    # generator matrix Q[i, j] = rate i -> j ; order CT, CD, BT, BD
    Q = np.zeros((4, 4))
    Q[0, 2] = k_onT * U
    Q[2, 0] = k_offT
    Q[1, 3] = k_onD * U
    Q[3, 1] = k_offD
    Q[2, 3] = k_h
    Q[0, 1] = k_h0
    Q[3, 2] = k_ex
    Q[1, 0] = k_ex
    np.fill_diagonal(Q, -Q.sum(1))
    # solve pi Q = 0, sum pi = 1
    Aeq = np.vstack([Q.T, np.ones(4)])
    rhs = np.append(np.zeros(4), 1.0)
    pi, *_ = np.linalg.lstsq(Aeq, rhs, rcond=None)
    return float(pi[2] + pi[3]), pi


def K_eff(U, occ):
    """effective dissociation constant that would give occupancy occ at U."""
    return U * (1.0 - occ) / occ


# ---------------------------------------------------------------- layer 3
# nascent-chain competition: extend the WO-02 system with
#   X  : nascent / newly released chains that are chaperone clients (uM)
#   BX : chaperone-nascent complex (uM)
# new chains: a fraction nu_c of the correctly made flux (1-eps) s_P enters X
# instead of N. X folds spontaneously (k_fX), binds chaperone (k_onX, k_offX),
# is completed by the cycle (k_catX) to N, or misfolds to U (k_xu).
X_, BX_ = sp.symbols("X BX", nonnegative=True)
nu_c, k_onX, k_offX, k_catX, k_fX, k_xu = sp.symbols(
    "nu_c k_onX k_offX k_catX k_fX k_xu", nonnegative=True)
EXT_STATES = list(m.X) + [X_, BX_]
EXT_PARAMS = list(m.PAR) + [nu_c, k_onX, k_offX, k_catX, k_fX, k_xu]
EXT_PARAM_NAMES = m.PARAM_NAMES + ["nu_c", "k_onX", "k_offX", "k_catX", "k_fX", "k_xu"]


def build_ext_rhs():
    N, U, B, A, C = m.X
    base = list(m.RHS_SYM)
    to_X = nu_c * (1 - m.eps) * m.s_P
    bindX = k_onX * C * X_
    relX = k_offX * BX_
    catX = k_catX * BX_
    foldX = k_fX * X_
    misX = k_xu * X_
    dN = base[0] - to_X + catX + foldX
    dU = base[1] + misX
    dB = base[2]
    dA = base[3]
    dC = base[4] - bindX + relX + catX
    dX = to_X - bindX + relX - foldX - misX - m.mu * X_
    dBX = bindX - relX - catX - m.mu * BX_
    return [dN, dU, dB, dA, dC, dX, dBX]


EXT_RHS = build_ext_rhs()
_ext_num = sp.lambdify((EXT_STATES, EXT_PARAMS), EXT_RHS, "numpy")
_ext_jac = sp.lambdify((EXT_STATES, EXT_PARAMS), sp.Matrix(EXT_RHS).jacobian(EXT_STATES), "numpy")


def ext_rhs(t, x, p):
    return np.asarray(_ext_num(x, p), dtype=float)


def ext_jac(t, x, p):
    return np.asarray(_ext_jac(x, p), dtype=float)


def ext_pvec(d):
    return np.array([float(d[k]) for k in EXT_PARAM_NAMES])


def ext_conservation():
    N, U, B, A, C = m.X
    dPT = sum(EXT_RHS[i] for i in (0, 1, 2, 3, 5, 6))
    dCT = EXT_RHS[4] + EXT_RHS[2] + EXT_RHS[6]
    PT = N + U + B + A + X_ + BX_
    CT = C + B + BX_
    okP = sp.simplify(dPT - (m.s_P - m.k_d * U - m.k_dA * A - m.mu * PT)) == 0
    okC = sp.simplify(dCT - (m.s_C - m.mu * CT)) == 0
    return okP, okC


def ext_reduces_symbolically():
    """with nu_c = 0 and X = BX = 0 the first five equations must be the WO-02
    equations term for term, and dX/dt = dBX/dt = 0 (so X, BX stay 0)."""
    sub = {nu_c: 0, X_: 0, BX_: 0}
    same = all(sp.simplify(EXT_RHS[i].subs(sub) - m.RHS_SYM[i]) == 0 for i in range(5))
    stays = all(sp.simplify(EXT_RHS[i].subs(sub)) == 0 for i in (5, 6))
    return same, stays


def steady_state(rhs, jac, p, x0, T=None):
    """integrate to steady state with a stiff solver, then polish with Newton."""
    from scipy.integrate import solve_ivp
    from scipy.optimize import root
    mu = p[m.PARAM_NAMES.index("mu")]
    T = T or 40.0 / mu
    sol = solve_ivp(rhs, (0, T), x0, method="LSODA", jac=jac, args=(p,),
                    rtol=1e-10, atol=1e-12)
    r = root(lambda x: rhs(0, x, p), sol.y[:, -1], jac=lambda x: jac(0, x, p),
             method="hybr", tol=1e-14)
    x = r.x if r.success else sol.y[:, -1]
    return x, float(np.max(np.abs(rhs(0, x, p))))
