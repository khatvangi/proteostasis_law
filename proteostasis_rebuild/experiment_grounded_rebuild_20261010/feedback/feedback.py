"""minimal two-variable sigma32-style negative feedback on chaperone capacity.

this module exists to check algebra (signs, conservation, bounds, delay threshold).
it carries NO physiological parameter values. every number used in the tests is a
dimensionless placeholder chosen to exercise a case, and must never be reported as
an estimate for E. coli.

state
  U  unfolded client that needs a chaperone (per volume)
  C  total sigma32-regulated chaperone pool (one lumped machine)

derived
  F = C*K/(K+U)   free chaperone (quasi-steady binding, saturating service)
  k*F             spare service capacity (unused turnover)

dynamics
  dU/dt = L - k*C*U/(K+U) - (mu(U)+delta)*U
  dC/dt = production(F(t-tau), U) - mu(U)*C

  production, mode "fixed":   b0 + bm*h(F)              (rate per volume)
  production, mode "growth":  mu(U)*P*(phi0 + phim*h(F)) (allocation x growth)
  h(F) = 1/(1+(F/Ks)^n)       sigma32 activity, repressed by free chaperone

  L      client influx needing service
  k      max turnover per chaperone
  K      client load at half-saturation of the chaperone
  delta  first-order degradation of unfolded client (protease route)
  tau    delay from sigma32 activity to new functional chaperone
"""
import numpy as np
from scipy.optimize import brentq


def findSigmaActivity(F, Ks, n):
    # hill repression of sigma32 activity by free chaperone; returns h and dh/dF
    r = (F / Ks) ** n
    h = 1.0 / (1.0 + r)
    dh = -n * r / (F * (1.0 + r) ** 2) if F > 0 else 0.0
    return h, dh


def findFree(C, U, K):
    # free chaperone: the fraction not occupied by client
    return C * K / (K + U)


def findGrowth(U, p):
    # growth (= dilution) rate; optional load-dependent growth inhibition
    if p.get("Ktox") is None:
        return p["mu"]
    return p["mu"] / (1.0 + (U / p["Ktox"]) ** p.get("ntox", 1))


def findProduction(U, F, p):
    # chaperone synthesis rate driven by sigma32 activity
    h, _ = findSigmaActivity(F, p["Ks"], p["n"])
    if p.get("mode", "fixed") == "fixed" and p.get("Ktox") is not None:
        # fixed synthesis with slowed growth would pile chaperone above Cmax; refuse it
        raise ValueError("growth inhibition (Ktox) requires mode='growth'")
    if p.get("mode", "fixed") == "growth":
        return findGrowth(U, p) * p["P"] * (p["phi0"] + p["phim"] * h)
    return p["b0"] + p["bm"] * h


def findRates(U, C, p, Fdelayed=None):
    # right-hand side; Fdelayed lets a caller pass F(t - tau)
    mu = findGrowth(U, p)
    F = findFree(C, U, p["K"]) if Fdelayed is None else Fdelayed
    dU = p["L"] - p["k"] * C * U / (p["K"] + U) - (mu + p.get("delta", 0.0)) * U
    dC = findProduction(U, F, p) - mu * C
    return dU, dC


def findCapMax(p):
    # largest chaperone pool feedback can ever hold: h = 1 everywhere
    if p.get("mode", "fixed") == "growth":
        return p["P"] * (p["phi0"] + p["phim"])
    return (p["b0"] + p["bm"]) / p["mu"]


def findCapMin(p):
    # basal pool with sigma32 fully repressed: h = 0
    if p.get("mode", "fixed") == "growth":
        return p["P"] * p["phi0"]
    return p["b0"] / p["mu"]


def findLoopTerms(U, C, p):
    """linearisation of the fixed-production, constant-growth model.

    returns a, b, c, e, mu such that the linear system is
      u' = -a u - b x
      x' =  c u(t-tau) - e x(t-tau) - mu x
    a: client removal slope   b: service gain per chaperone
    c: load -> synthesis gain e: chaperone self-repression gain via free pool
    """
    K, k, mu = p["K"], p["k"], p["mu"]
    F = findFree(C, U, K)
    _, dh = findSigmaActivity(F, p["Ks"], p["n"])
    a = k * C * K / (K + U) ** 2 + mu + p.get("delta", 0.0)
    b = k * U / (K + U)
    c = p["bm"] * dh * (-C * K / (K + U) ** 2)   # dh<0, dF/dU<0 -> c>0
    e = -p["bm"] * dh * K / (K + U)              # dh<0, dF/dC>0 -> e>0
    return a, b, c, e, mu


def findJacobian(U, C, p):
    # undelayed jacobian assembled from the loop terms
    a, b, c, e, mu = findLoopTerms(U, C, p)
    return np.array([[-a, -b], [c, -(mu + e)]])


def findJacobianNumeric(U, C, p, eps=1e-7):
    # central finite differences of the full right-hand side (any mode)
    J = np.zeros((2, 2))
    x = np.array([U, C], float)
    for j in range(2):
        d = np.zeros(2)
        d[j] = eps * max(1.0, abs(x[j]))
        fp = np.array(findRates(*(x + d), p))
        fm = np.array(findRates(*(x - d), p))
        J[:, j] = (fp - fm) / (2 * d[j])
    return J


def findCOnUNullcline(U, p):
    # the C that makes dU/dt = 0 at given U
    mu = findGrowth(U, p)
    return (p["L"] - (mu + p.get("delta", 0.0)) * U) * (p["K"] + U) / (p["k"] * U)


def findEquilibria(p, grid=4000):
    # all equilibria: scan dC/dt along the U-nullcline for sign changes
    # constant growth: U <= L/(mu+delta) at any equilibrium. with growth inhibition
    # the dilution flux can fall at high U, so the caller must give the scan range Uscan
    Uhi = p.get("Uscan") or p["L"] / (findGrowth(0.0, p) + p.get("delta", 0.0))
    Us = np.geomspace(Uhi * 1e-9, Uhi * (1 - 1e-12), grid)
    G = lambda U: findRates(U, findCOnUNullcline(U, p), p)[1]
    vals = np.array([G(U) for U in Us])
    roots = []
    for i in range(len(Us) - 1):
        # only states with a non-negative chaperone pool are physical
        if findCOnUNullcline(Us[i], p) < 0 or findCOnUNullcline(Us[i + 1], p) < 0:
            continue
        if np.sign(vals[i]) != np.sign(vals[i + 1]):
            U = brentq(G, Us[i], Us[i + 1], xtol=1e-14 * Uhi)
            roots.append((U, findCOnUNullcline(U, p)))
    return roots


def findFluxPartition(U, C, p):
    # fate of influx at a state: chaperone-served, protease-degraded, diluted unfolded
    served = p["k"] * C * U / (p["K"] + U)
    degraded = p.get("delta", 0.0) * U
    diluted = findGrowth(U, p) * U
    return served, degraded, diluted


def findLoopGain(a, b, c, e, mu):
    # zero-frequency loop gain G = (e*a + b*c)/(a*mu)
    return (e * a + b * c) / (a * mu)


def findCriticalDelay(a, b, c, e, mu):
    """smallest delay at which a root of
       (s+a)(s+mu) + (e s + e a + b c) exp(-s tau) = 0
    crosses the imaginary axis. returns (omega, tau) or None if stable for all tau."""
    A2 = a * a + mu * mu - e * e
    A0 = (a * mu) ** 2 - (e * a + b * c) ** 2
    zs = np.roots([1.0, A2, A0])
    zs = [z.real for z in zs if abs(z.imag) < 1e-12 and z.real > 0]
    best = None
    for z in zs:
        w = np.sqrt(z)
        P = (1j * w + a) * (1j * w + mu)
        Q = e * 1j * w + e * a + b * c
        # exp(-i w tau) = -P/Q  ->  tau = (-arg(-P/Q) mod 2pi)/w
        tau = (-np.angle(-P / Q)) % (2 * np.pi) / w
        if best is None or tau < best[1]:
            best = (w, tau)
    return best


def findLinearDdeGrowth(a, b, c, e, mu, tau, T, dt=1e-3):
    # euler integration of the linearised delayed loop; returns log of |x| growth
    nd = int(round(tau / dt))
    n = int(round(T / dt))
    u = np.zeros(n + nd + 1)
    x = np.zeros(n + nd + 1)
    u[: nd + 1] = 1e-3      # constant history
    for i in range(nd, n + nd):
        ud, xd = u[i - nd], x[i - nd]
        u[i + 1] = u[i] + dt * (-a * u[i] - b * x[i])
        x[i + 1] = x[i] + dt * (c * ud - e * xd - mu * x[i])
    amp = lambda s: np.max(np.abs(np.r_[u[s], x[s]]))
    m = n // 4
    return np.log(amp(slice(-m, None)) / amp(slice(nd + m, nd + 2 * m)))
