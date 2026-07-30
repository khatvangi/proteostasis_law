#!/usr/bin/env python3
"""
Figure 1 -- the proteostasis envelope and its regimes.

the reduced burden model, in nondimensional form:

    dx/dtau = lambda - x - rho*x/(1+x) + chi*x^2

    lambda : burden inflow          (translation-derived damage flux)
    x      : effective burden pool  (misfolded / damaged species)
    -x     : linear clearance
    -rho*x/(1+x) : saturable rescue (chaperone / protease capacity)
    +chi*x^2     : superlinear aggregation feedback

fixed points satisfy lambda = g(x) where g(x) = x + rho*x/(1+x) - chi*x^2.

g'(x) = 1 + rho/(1+x)^2 - 2*chi*x is positive at x=0 and negative as x grows
whenever chi > 0, so g has exactly one interior maximum. that maximum IS the
saddle-node: below it two fixed points exist (a stable low-burden node and an
unstable upper threshold); at it they merge; above it no bounded state remains.

NOTE ON A CORRECTION: the previous draft's legend described this transition as
leading to "bistability". that is wrong for this model -- there is no second
stable high-burden state, because +chi*x^2 dominates at large x and the burden
pool diverges. what the fold produces is loss of the low-burden state, i.e.
runaway, not a second attractor. the legend below says so.

this figure is a qualitative statement about regime STRUCTURE. parameters are
not fitted to any organism, and no claim is made that any organism sits near
the threshold.
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import figstyle as fs

ROOT = Path(__file__).resolve().parent.parent
FIGS = ROOT / "figures"

RHO = 4.0    # rescue capacity
CHI = 0.15   # aggregation feedback


def g(x, rho=RHO, chi=CHI):
    """inflow required to hold the burden pool at x (the fixed-point curve)."""
    return x + rho * x / (1.0 + x) - chi * x ** 2


def dxdt(x, lam, rho=RHO, chi=CHI):
    return lam - x - rho * x / (1.0 + x) + chi * x ** 2


def main():
    fs.setup()
    FIGS.mkdir(parents=True, exist_ok=True)

    # locate the fold: the interior maximum of g
    xs = np.linspace(1e-6, 12, 400_000)
    gs = g(xs)
    k = int(np.argmax(gs))
    x_fold, lam_crit = float(xs[k]), float(gs[k])

    fig, axes = plt.subplots(1, 2, figsize=(fs.COL_DOUBLE, 2.7))

    # ---------------- panel a: the fold ----------------
    ax = axes[0]
    stable = xs <= x_fold
    df_s = pd.DataFrame({"lam": gs[stable], "x": xs[stable]})
    df_u = pd.DataFrame({"lam": gs[~stable], "x": xs[~stable]})
    # only the part of the upper branch that is a real fixed point (lambda > 0)
    df_u = df_u[df_u.lam > 0]

    sns.lineplot(data=df_s, x="lam", y="x", ax=ax, color=fs.C["burden"],
                 linewidth=1.8, label="stable (low burden)")
    sns.lineplot(data=df_u, x="lam", y="x", ax=ax, color=fs.C["burden"],
                 linewidth=1.3, linestyle="--", label="unstable threshold")

    ax.axvline(lam_crit, color=fs.C["alert"], linewidth=1.0, zorder=1)
    ax.plot([lam_crit], [x_fold], marker="o", ms=4.5,
            color=fs.C["alert"], zorder=5)

    y_top = float(df_u.x.max())
    ax.axvspan(0, 0.62 * lam_crit, color=fs.C["capacity"], alpha=0.10, lw=0)
    ax.axvspan(0.62 * lam_crit, lam_crit, color="#e9c46a", alpha=0.22, lw=0)
    ax.axvspan(lam_crit, 1.35 * lam_crit, color=fs.C["alert"], alpha=0.10, lw=0)

    for xpos, lab in ((0.31 * lam_crit, "buffered"),
                      (0.81 * lam_crit, "vulnerable"),
                      (1.17 * lam_crit, "overload")):
        ax.text(xpos, y_top * 0.94, lab, ha="center", va="top",
                fontsize=6.8, style="italic", color=fs.C["muted"])

    ax.annotate(r"$\lambda_{\rm crit}$ (saddle-node)",
                xy=(lam_crit, x_fold), xytext=(0.42 * lam_crit, x_fold * 1.5),
                fontsize=6.8, color=fs.C["alert"],
                arrowprops=dict(arrowstyle="->", lw=0.7, color=fs.C["alert"]))

    ax.set_xlim(0, 1.35 * lam_crit)
    ax.set_ylim(0, y_top)
    ax.set_xlabel(r"burden inflow  $\lambda$")
    ax.set_ylabel(r"burden pool  $x^{*}$")
    ax.set_title("Saturable rescue + superlinear feedback\ngive a threshold")
    ax.legend(loc="lower right", fontsize=6.5)
    fs.panel_label(ax, "a")

    # ---------------- panel b: flow at three inflows ----------------
    ax = axes[1]
    xg = np.linspace(0, 1.15 * float(df_u.x.max()), 600)
    cases = [(0.45, "buffered", fs.C["capacity"]),
             (0.90, "vulnerable", "#d9a13a"),
             (1.10, "overload", fs.C["alert"])]

    for frac, lab, col in cases:
        lam = frac * lam_crit
        d = pd.DataFrame({"x": xg, "f": dxdt(xg, lam)})
        sns.lineplot(data=d, x="x", y="f", ax=ax, color=col, linewidth=1.5,
                     label=rf"$\lambda = {frac:.2f}\,\lambda_{{\rm crit}}$  ({lab})")
        # mark the stable fixed point where the curve crosses zero downward
        s = np.sign(d.f.to_numpy())
        cross = np.where(np.diff(s) < 0)[0]
        if len(cross):
            i = cross[0]
            ax.plot([xg[i]], [0], marker="o", ms=4, color=col, zorder=6)

    ax.axhline(0, color=fs.C["muted"], linewidth=0.7, zorder=1)
    ax.set_xlabel(r"burden pool  $x$")
    ax.set_ylabel(r"$dx/d\tau$")
    ax.set_title("Below threshold a stable state exists;\nabove it, none does")
    ax.legend(loc="upper left", fontsize=6.2)
    fs.panel_label(ax, "b")

    sns.despine(fig=fig)
    fig.tight_layout(w_pad=2.0)
    fs.save(fig, str(FIGS / "Fig1_envelope"))

    print(f"  rho = {RHO}, chi = {CHI}")
    print(f"  fold at x = {x_fold:.3f}, lambda_crit = {lam_crit:.3f}")
    print("  (illustrative parameters; not fitted to data)")


if __name__ == "__main__":
    main()
