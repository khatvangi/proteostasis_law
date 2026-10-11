"""independent checks of feedback/MINIMAL_DYNAMICS.md. every number is a dimensionless
placeholder; nothing here is an E. coli estimate. run: python check_feedback_claims.py
seeds are fixed; output is deterministic up to floating point.

linearised delayed loop (feedback.findLoopTerms convention):
  u' = -a u - b x
  x' =  c u(t-tau) - e x(t-tau) - mu x
characteristic function f(s) = (s+a)(s+mu) + (e(s+a)+bc) exp(-s tau)
"""
import sys, os
import numpy as np
from scipy.integrate import solve_ivp

# read-only import of the audited module (only for cross-checks of its own functions)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "feedback"))
import feedback as fb


def findGain(a, b, c, e, mu):
    return (e * a + b * c) / (a * mu)


def findCrossingSquares(a, b, c, e, mu):
    # |P(iw)|^2 = |Q(iw)|^2  <=>  z^2 + A2 z + A0 = 0 with z = w^2
    A2 = a * a + mu * mu - e * e
    A0 = (a * mu) ** 2 - (e * a + b * c) ** 2
    return A2, A0, [z.real for z in np.roots([1.0, A2, A0]) if abs(z.imag) < 1e-12 and z.real > 0]


def findCheb(N):
    # trefethen chebyshev differentiation matrix on [-1,1], x[0] = 1
    x = np.cos(np.pi * np.arange(N + 1) / N)
    cc = np.r_[2.0, np.ones(N - 1), 2.0] * (-1.0) ** np.arange(N + 1)
    X = np.tile(x, (N + 1, 1)).T
    D = np.outer(cc, 1 / cc) / (X - X.T + np.eye(N + 1))
    return D - np.diag(D.sum(axis=1)), x


def findRightmostRoot(a, b, c, e, mu, tau, N=60):
    # spectral discretisation of the DDE generator (breda-style), then newton polish
    A0 = np.array([[-a, -b], [0.0, -mu]])
    A1 = np.array([[0.0, 0.0], [c, -e]])
    if tau == 0:
        lam = np.linalg.eigvals(A0 + A1)
        return lam[np.argmax(lam.real)]
    D, _ = findCheb(N)
    D = D * 2.0 / tau                      # theta = tau*(x-1)/2 on [-tau, 0]
    M = np.kron(D, np.eye(2))
    M[:2, :] = 0.0
    M[:2, :2] = A0                         # theta = 0
    M[:2, -2:] = A1                        # theta = -tau
    lam = np.linalg.eigvals(M)
    lam = lam[np.argsort(-lam.real)][:6]   # polish the few rightmost, keep the best
    f = lambda s: (s + a) * (s + mu) + (e * (s + a) + b * c) * np.exp(-s * tau)
    df = lambda s: 2 * s + a + mu + (e - tau * (e * (s + a) + b * c)) * np.exp(-s * tau)
    out = []
    for s in lam:
        for _ in range(50):
            s = s - f(s) / df(s)
        if abs(f(s)) < 1e-8 * (1 + abs(s) ** 2):
            out.append(s)
    return max(out, key=lambda s: s.real) if out else lam[0]


def drawLoop(rng, n):
    return 10 ** rng.uniform(-3, 3, size=(n, 5))


def checkThreshold():
    print("== 1. delay-independent stability: exact condition ==")
    rng = np.random.default_rng(1)
    P = drawLoop(rng, 200000)
    nlo = nlo_bad = nlo_e = nhi = nhi_ok = 0
    for a, b, c, e, mu in P:
        G = findGain(a, b, c, e, mu)
        _, _, zs = findCrossingSquares(a, b, c, e, mu)
        if G < 1:
            nlo += 1
            nlo_e += (e >= mu) or (e * e > a * a + mu * mu)
            nlo_bad += len(zs) > 0
        elif G > 1:
            nhi += 1
            nhi_ok += len(zs) == 1
    print(f"G<1 draws: {nlo}; with e>=mu or e^2>a^2+mu^2: {nlo_e}; with an imaginary-axis crossing: {nlo_bad}")
    print(f"G>1 draws: {nhi}; with exactly one crossing frequency: {nhi_ok}")

    # spectral check: G<1 stays stable at every delay tried, including very long ones
    rng = np.random.default_rng(2)
    worst = -np.inf
    tried = 0
    while tried < 40:
        a, b, c, e, mu = 10 ** rng.uniform(-1, 1, size=5)
        if findGain(a, b, c, e, mu) >= 1:
            continue
        tried += 1
        for tau in [0.0, 0.1, 1.0, 10.0, 50.0]:
            N = int(min(200, 40 + 3 * tau * max(a, mu, e, 1)))
            worst = max(worst, findRightmostRoot(a, b, c, e, mu, tau, N).real)
    print(f"40 G<1 loops x 5 delays (0..50): largest rightmost-root real part = {worst:.3e} (must be < 0)")

    # G>1: crossing at tau*, stable just below, unstable just above, and cross-check fb
    rng = np.random.default_rng(3)
    tried = agree = flips = 0
    while tried < 40:
        a, b, c, e, mu = 10 ** rng.uniform(-1, 1, size=5)
        if findGain(a, b, c, e, mu) <= 1.05:
            continue
        tried += 1
        w, tau = fb.findCriticalDelay(a, b, c, e, mu)
        s_at = findRightmostRoot(a, b, c, e, mu, tau, 120)
        agree += abs(s_at.real) < 1e-6 and abs(abs(s_at.imag) - w) < 1e-6 * max(1, w)
        lo = findRightmostRoot(a, b, c, e, mu, 0.97 * tau, 120).real
        hi = findRightmostRoot(a, b, c, e, mu, 1.03 * tau, 120).real
        flips += lo < 0 < hi
    print(f"40 G>1 loops: root on imaginary axis at fb tau* in {agree}/40; stable at 0.97tau*, unstable at 1.03tau* in {flips}/40")

    # boundary G = 1 exactly: no crossing, stable for every finite tau, margin -> 0
    a, mu, e = 1.0, 1.0, 0.5
    bc = a * mu - e * a                     # makes e*a + b*c = a*mu
    b, c = 1.0, bc
    print(f"G=1 boundary (a=mu=1, e=0.5, bc=0.5): G={findGain(a,b,c,e,mu):.12f}, crossings={findCrossingSquares(a,b,c,e,mu)[2]}")
    for tau in [1.0, 10.0, 100.0]:
        N = int(40 + 4 * tau)
        print(f"   tau={tau:6.1f}  rightmost root real part = {findRightmostRoot(a,b,c,e,mu,tau,N).real:.3e}")


def checkWrongInequalities():
    print("== 2. counterexamples to candidate inequalities ==")
    # e < mu and e^2 < a^2+mu^2 hold, but G > 1, so a finite tau* exists
    a, b, c, e, mu = 1.0, 1.0, 2.0, 0.5, 1.0
    w, tau = fb.findCriticalDelay(a, b, c, e, mu)
    print(f"a=1,b=1,c=2,e=0.5,mu=1: e<mu, e^2<a^2+mu^2, G={findGain(a,b,c,e,mu):.2f} -> tau*={tau:.4f}, w={w:.4f}")
    print(f"   rightmost Re at 0.9tau*: {findRightmostRoot(a,b,c,e,mu,0.9*tau).real:+.4e};  at 1.1tau*: {findRightmostRoot(a,b,c,e,mu,1.1*tau).real:+.4e}")
    # the doc's exception 'G<1 but e^2 > a^2+mu^2': shown impossible in section 1; try a direct search
    rng = np.random.default_rng(4)
    P = drawLoop(rng, 200000)
    hits = sum(1 for a, b, c, e, mu in P if findGain(a, b, c, e, mu) < 1 and e * e > a * a + mu * mu)
    print(f"random search, 2e5 draws: G<1 and e^2>a^2+mu^2 found {hits} times")
    # with e^2 > a^2 + mu^2 the gain is always > 1
    a, b, c, e, mu = 1.0, 0.01, 0.01, 2.0, 1.0
    print(f"a=1,mu=1,e=2,bc=1e-4: e^2>a^2+mu^2 forces G={findGain(a,b,c,e,mu):.4f} > 1; fb tau*={fb.findCriticalDelay(a,b,c,e,mu)[1]:.4f}")


def checkSigmaDerivative():
    print("== 3. h'(F) at F = 0 in feedback.findSigmaActivity ==")
    for n in [1, 2, 0.5]:
        _, dh = fb.findSigmaActivity(0.0, 1.0, n)
        true = {1: -1.0, 2: 0.0, 0.5: -np.inf}[n]
        print(f"   n={n}: code returns {dh}, true derivative at F=0 = {true}")


def findOde(t, y, p):
    # undelayed fixed-production model, self-contained (no fb dependency)
    U, C = y
    F = C * p["K"] / (p["K"] + U)
    h = 1.0 / (1.0 + (F / p["Ks"]) ** p["n"])
    mu = p["mu"](t) if callable(p["mu"]) else p["mu"]
    dU = p["L"] - p["k"] * C * U / (p["K"] + U) - (mu + p["delta"]) * U
    dC = p["b0"] + p["bm"] * h - mu * C
    return [dU, dC]


def checkCeiling():
    print("== 4. 'C never above Cmax' / 'served <= k*Cmax at any state' ==")
    p = dict(L=5.0, k=2.0, K=1.0, mu=0.5, delta=0.0, b0=0.05, bm=1.0, Ks=1.0, n=2)
    Cmax = (p["b0"] + p["bm"]) / p["mu"]
    print(f"   placeholder: Cmax={Cmax:.3f}, k*Cmax={p['k']*Cmax:.3f}, L={p['L']} (L > k*Cmax)")
    # (a) initial pool above the ceiling (e.g. a pre-induced cell moved to new conditions)
    sol = solve_ivp(findOde, [0, 20], [1.0, 3 * Cmax], args=(p,), dense_output=True, rtol=1e-9, atol=1e-12)
    t = np.linspace(0, 20, 4001)
    U, C = sol.sol(t)
    served = p["k"] * C * U / (p["K"] + U)
    above = t[C > Cmax]
    print(f"   (a) C(0)=3Cmax: C>Cmax until t={above.max():.2f}; max served={served.max():.3f} vs k*Cmax={p['k']*Cmax:.3f}; "
          f"max served/L={served.max()/p['L']:.3f}")
    # (b) growth-rate step: steady state at mu=0.25, then mu jumps to 0.5 (fixed production)
    p1 = dict(p, mu=0.25)
    ss = solve_ivp(findOde, [0, 400], [1.0, 1.0], args=(p1,), rtol=1e-10, atol=1e-12).y[:, -1]
    sol = solve_ivp(findOde, [0, 20], ss, args=(p,), dense_output=True, rtol=1e-9, atol=1e-12)
    U, C = sol.sol(t)
    print(f"   (b) mu step 0.25->0.5: C(0)={ss[1]:.3f} vs new Cmax={Cmax:.3f}; C>Cmax until t={t[C > Cmax].max() if (C > Cmax).any() else 0:.2f}")
    # (c) fixed production with mu -> 0: ceiling (b0+bm)/mu diverges
    for mu in [0.5, 0.05, 0.005]:
        print(f"   (c) fixed mode, mu={mu}: Cmax=(b0+bm)/mu={(p['b0']+p['bm'])/mu:.1f}")


def checkDegradationLowersServed():
    print("== 5. raising delta always lowers chaperone-served flux (constant mu, fixed mode) ==")
    rng = np.random.default_rng(5)
    lg = lambda lo, hi: 10 ** rng.uniform(lo, hi)
    n = bad = 0
    for _ in range(100):
        p = dict(L=lg(-1, 1), k=lg(-1, 1), K=lg(-1, 1), mu=lg(-2, 0), delta=lg(-3, 0),
                 b0=lg(-3, -1), bm=lg(-1, 1), Ks=lg(-1, 1), n=int(rng.choice([1, 2, 4])))
        U1, C1 = fb.findEquilibria(p)[0]
        U2, C2 = fb.findEquilibria(dict(p, delta=p["delta"] * 2))[0]
        s1 = fb.findFluxPartition(U1, C1, p)[0]
        s2 = fb.findFluxPartition(U2, C2, dict(p, delta=p["delta"] * 2))[0]
        n += 1
        bad += s2 > s1 * (1 + 1e-10)
    print(f"   delta doubled in {n} draws: served flux rose in {bad}")


if __name__ == "__main__":
    checkThreshold()
    checkWrongInequalities()
    checkSigmaDerivative()
    checkCeiling()
    checkDegradationLowersServed()
