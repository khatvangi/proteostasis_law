#!/usr/bin/env python3
"""First-pass, reproducible proteostasis redraws and numerical checks.

Only standard NumPy/Matplotlib are required.  The calculations are deliberately
kept independent of the legacy ODE implementation.
"""
from pathlib import Path
import math
import numpy as np

OUT = Path(__file__).resolve().parent


def g(x, rho=4.0, chi=0.15):
    return x + rho * x / (1.0 + x) - chi * x * x


def gprime(x, rho=4.0, chi=0.15):
    return 1.0 + rho / (1.0 + x) ** 2 - 2.0 * chi * x


def finite_pool(M_total, C_total=50.0, Kd=1.0):
    """Return (C_free, C_bound) in concentration units.

    Equations: C_T=C_f+C_b, M_T=M_f+C_b, Kd=C_f*M_f/C_b.
    The numerically stable physical root is obtained by solving the quadratic
    in C_b and selecting 0 <= C_b <= min(C_T,M_T).
    """
    M_total = np.asarray(M_total, dtype=float)
    disc = (C_total + M_total + Kd) ** 2 - 4.0 * C_total * M_total
    C_bound = 2.0 * C_total * M_total / (
        C_total + M_total + Kd + np.sqrt(disc)
    )
    return C_total - C_bound, C_bound


def approximate_free(M_total, C_total=50.0, Kd=1.0):
    return C_total / (1.0 + np.asarray(M_total, dtype=float) / Kd)


def folding_rate(C_free, k_obs_max=1.0, Kd=1.0):
    return k_obs_max * C_free / (C_free + Kd)


def stationary_points(rho=4.0, chi=0.15):
    # (1 + x)^2 g'(x) = 0 gives a cubic; retain real x > 0.
    # Expanded: -2 chi x^3 + (1 - 4 chi) x^2 + (2 - 2 chi) x + (1 + rho).
    coeff = [-2 * chi, 1 - 4 * chi, 2 - 2 * chi, 1 + rho]
    roots = np.roots(coeff)
    return sorted(float(r.real) for r in roots if abs(r.imag) < 1e-9 and r.real > 0)


def make_plots():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                         "axes.spines.right": False})

    # Panel 1: g(x) and the scalar flow F=lambda-g.
    x = np.linspace(0, 12, 1600)
    rho, chi = 4.0, 0.15
    xs = stationary_points(rho, chi)
    xmax = xs[0]
    lams = [0.0, 2.0, 4.80218919587308, 5.2]
    fig, ax = plt.subplots(figsize=(7.2, 4.5))
    ax.plot(x, g(x, rho, chi), color="#16324F", lw=2.5, label=r"$g(x)$")
    ax.axhline(0, color="0.3", lw=0.8)
    ax.axvline(xmax, color="#D1495B", ls="--", lw=1.2,
               label=fr"maximum: $x={xmax:.4f}$")
    ax.scatter([xmax], [g(xmax, rho, chi)], color="#D1495B", zorder=5)
    xzero = (-(1.0 - chi) - np.sqrt((1.0 - chi) ** 2 + 4.0 * chi * (rho + 1.0))) / (-2.0 * chi)
    ax.axvline(xzero, color="#777777", ls=":", lw=1.2,
               label=fr"$g=0$ at $x={xzero:.4f}$")
    # Mark equilibria for lambda=2: stable low root and unstable upper root.
    lam = 2.0
    coeff = [-chi, 1.0 - chi, rho + 1.0, -lam]
    eq = sorted(float(r.real) for r in np.roots(coeff)
                if abs(r.imag) < 1e-8 and r.real > 0)
    for z in eq:
        ax.scatter([z], [lam], s=55, color="#2A9D8F" if gprime(z, rho, chi) > 0 else "#E76F51",
                   zorder=6)
    ax.axhline(lam, color="#2A9D8F", ls="-.", lw=1.0,
               label=r"$λ=2$: stable/unstable equilibria")
    ax.axhline(5.2, color="#D1495B", ls="-.", lw=1.0,
               label=r"$λ=5.2>g_{max}$: overload")
    ax.set(xlabel=r"dimensionless load $x$", ylabel=r"$g(x)=x+\rho x/(1+x)-\chi x^2$",
           title=r"Scalar reduced model: $F(x)=\lambda-g(x)$")
    ax.set_ylim(-8, 6.5)
    ax.legend(frameon=False, fontsize=8, loc="lower left")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(OUT / "scalar_g_equilibria.png", dpi=220)
    fig.savefig(OUT / "scalar_g_equilibria.svg")
    plt.close(fig)

    # Panel 2: finite-pool versus ligand-excess approximation.
    M = np.linspace(0, 300, 1200)
    cf_exact, cb = finite_pool(M)
    cf_approx = approximate_free(M)
    v_exact = folding_rate(cf_exact)
    v_approx = folding_rate(cf_approx)
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.1))
    ax = axes[0]
    ax.plot(M, cf_exact, color="#16324F", lw=2.4, label=r"exact $C_f$")
    ax.plot(M, cf_approx, color="#D1495B", lw=2.0, ls="--", label=r"approx. $C_T/(1+M_T/K_d)$")
    ax.set(xlabel=r"total misfolded ligand $M_T$ ($\mu$M)", ylabel=r"free chaperone $C_f$ ($\mu$M)",
           title="Finite-pool mass balance")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.2)
    ax = axes[1]
    ax.plot(M, v_exact, color="#16324F", lw=2.4, label=r"exact $v_{fold}$")
    ax.plot(M, v_approx, color="#D1495B", lw=2.0, ls="--", label=r"approx. $v_{fold}$")
    ax.set(xlabel=r"total misfolded ligand $M_T$ ($\mu$M)", ylabel=r"relative folding rate ($k_{obs,max}=1$)",
           title="Rate consequence of binding closure")
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig.savefig(OUT / "finite_pool_vs_approx.png", dpi=220)
    fig.savefig(OUT / "finite_pool_vs_approx.svg")
    plt.close(fig)


def main():
    make_plots()
    xm = stationary_points()[0]
    print(f"stationary_x={xm:.12f}")
    print(f"g_max={g(xm):.12f}")
    xzero = (-(1.0 - 0.15) - np.sqrt((1.0 - 0.15) ** 2 + 4.0 * 0.15 * 5.0)) / (-2.0 * 0.15)
    print(f"positive_zero={xzero:.12f}")
    for M in (0, 1, 5, 10, 25, 50, 100, 150, 200, 300):
        ca = float(approximate_free(M))
        ce, cb = finite_pool(M)
        va, ve = float(folding_rate(ca)), float(folding_rate(ce))
        print(f"M={M:6.1f} Capprox={ca:.8f} Cexact={float(ce):.8f} "
              f"Cbound={float(cb):.8f} vratio={ve/va:.8f}")


if __name__ == "__main__":
    main()
