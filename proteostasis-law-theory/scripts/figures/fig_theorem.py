"""Figure 1 -- the theorem, in two panels.

Panel (a) shows the Lagrange condition rather than asserting it: the aggregate
nullcline {G = 0}, contours of total removal R, and at the solved fold the two
gradients drawn as arrows that point the same way.

Panel (b) shows the equilibrium branch a*(j) with its VERTICAL tangent at the
saddle-node. Differentiating the equilibrium condition gives

    da*/dj = G_u / det J

so a vertical tangent is det J = 0 (the theorem) and a HORIZONTAL tangent is
G_u = 0, which is where infinitesimal homeostasis of a with respect to j would
sit. See `homeostasisPointExists` for why the second does not occur here.

This figure needs NO run root: it is one parameter point, computed from the
model, so a clean checkout reproduces it exactly.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

# manuscript figure number, by ORDER OF FIRST MENTION in bmb_v4.md.
# filenames are deliberately semantic so a reorder touches this line only.
FIGURE = "fig1"

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3",
           REPO_ROOT / "scripts" / "figures"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import _figstyle as F  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from proteostasis import model as M  # noqa: E402
import fold_theorem as FT  # noqa: E402
import boundary_structure as B  # noqa: E402
import dilution as D  # noqa: E402


def fields(p: M.Params, u_lo=2e-3, u_hi=2.0, a_lo=1e-4, a_hi=1.0, n=200):
    """G and R on one log-spaced grid, from a single flux evaluation per point.

    A LOG grid matters: the lower nullcline branch has `a` of order 1e-3 at small
    `u`, which a linear grid over [0, 1] cannot resolve, and the branch would be
    truncated at the low-influx end of panel (b).
    """
    ug = np.geomspace(u_lo, u_hi, n)
    ag = np.geomspace(a_lo, a_hi, n)
    UU, AA = np.meshgrid(ug, ag)
    GG = np.full(UU.shape, np.nan)
    RR = np.full(UU.shape, np.nan)
    for i in range(n):
        for k in range(n):
            try:
                f = M.fluxes(UU[i, k], AA[i, k], p)
            except (M.ModelError, OverflowError):
                continue
            GG[i, k] = (f["nucleate"] + f["grow"] - f["disaggregate"]
                        - f["degrade_a"])
            RR[i, k] = f["refold"] + f["degrade_u"] + f["degrade_a"]
    return UU, AA, GG, RR


def nullclineSegments(UU, AA, GG):
    """trace {G = 0} as a contour rather than by root-finding per column.

    Root-finding at fixed `u` loses the curve near its turning point, where the
    two roots merge: the branches come back as two disconnected pieces with the
    fold sitting in the gap between them. A contour follows the curve THROUGH the
    turn, which is the part of the figure that has to be right.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cs = ax.contour(UU, AA, GG, levels=[0.0])
    segs = [s for s in cs.allsegs[0] if len(s) > 2]
    plt.close(fig)
    return segs


def branchStats(p: M.Params, segs):
    """(u, a, j, G_u, det J) along the traced curve, ordered as the curve runs."""
    rows = []
    for seg in segs:
        for u, a in seg:
            if u <= 0.0 or a <= 0.0:
                continue
            try:
                Gu, _ = FT._centralGradient(FT.aggregateG, float(u), float(a), p)
                det = float(np.linalg.det(M.jacobian(float(u), float(a), p)))
                rows.append((float(u), float(a), FT.removalR(float(u), float(a), p),
                             Gu, det))
            except (M.ModelError, np.linalg.LinAlgError, OverflowError):
                continue
    return np.array(rows)


def turningPoint(p: M.Params, stats: np.ndarray):
    """the point where `G_a = 0`, located on the traced curve.

    `du*/dj = -G_a/det J` by Cramer on the second column, so the SOLUBLE pool has
    a horizontal tangent exactly where `G_a` vanishes. That is the nullcline's own
    turning point -- where `{G = 0}` runs vertical in the (u, a) plane -- and it is
    generic, since `det J = R_a . G_u` there rather than zero. It is also precisely
    where root-finding in `a` at fixed `u` loses the curve: the numerical
    difficulty and the horizontal tangent are one locus.
    """
    def Ga(u, a):
        return FT._centralGradient(FT.aggregateG, float(u), float(a), p)[1]

    ga = np.array([Ga(u, a) for u, a, _, _, _ in stats])
    idx = np.where(np.sign(ga[:-1]) != np.sign(ga[1:]))[0]
    if not len(idx):
        return None
    i = int(idx[0])
    u0, a0 = stats[i, 0], stats[i, 1]
    u1, a1 = stats[i + 1, 0], stats[i + 1, 1]
    try:
        s = brentq(lambda s: Ga(u0 + s * (u1 - u0), a0 + s * (a1 - a0)),
                   0.0, 1.0, xtol=1e-14)
    except ValueError:
        s = 0.5
    ut, at = u0 + s * (u1 - u0), a0 + s * (a1 - a0)
    det = float(np.linalg.det(M.jacobian(ut, at, p)))
    _, Ra = FT._centralGradient(FT.removalR, ut, at, p)
    Gu, _ = FT._centralGradient(FT.aggregateG, ut, at, p)
    return {"u": ut, "a": at, "j": FT.removalR(ut, at, p), "det_J": det,
            "Ra_times_Gu": Ra * Gu, "G_a": Ga(ut, at),
            "n_sign_changes": int(len(idx))}


def homeostasisPointExists(stats: np.ndarray) -> dict:
    """does da*/dj vanish anywhere, i.e. is there an infinitesimal-homeostasis point?

    `da*/dj = G_u/det J`, so this asks whether `G_u` changes sign. It does not,
    and the reason is structural rather than a property of these parameters:
    every term of `G = v_nuc + v_grow - v_dis - v_degA` pushes `G_u` the same way.
    Nucleation and aggregate growth rise with `u` directly; raising `u` also
    sequesters chaperone and protease, which LOWERS disaggregation and aggregate
    clearance, and those enter `G` with a minus sign. All four contributions are
    positive, so `G_u > 0` and no horizontal tangent exists.
    """
    Gu = stats[:, 3]
    return {"exists": bool((Gu < 0).any()),
            "min_G_u": float(Gu.min()), "max_G_u": float(Gu.max()),
            "n_points": int(len(Gu))}


def build():
    F.setStyle()
    p = M.Params().validate()
    jc, us, as_ = FT.foldSolve(p)
    ident = FT.determinantIdentity(us, as_, p)

    UU, AA, GG, RR = fields(p)
    segs = nullclineSegments(UU, AA, GG)
    stats = branchStats(p, segs)
    homeo = homeostasisPointExists(stats)
    turn = turningPoint(p, stats)

    fig, (ax_a, ax_b) = plt.subplots(
        1, 2, figsize=(F.W_SINGLE, 0.44 * F.W_SINGLE))

    # ---------------- panel (a): the Lagrange condition ----------------
    u_max, a_max = 1.5, 0.85
    cs = ax_a.contour(UU, AA, RR, levels=8, colors="0.78", linewidths=0.4)
    ax_a.clabel(cs, inline=True, fontsize=4.2, fmt="%.2f")
    for seg in segs:
        ax_a.plot(seg[:, 0], seg[:, 1], color="#1b3a6b", lw=1.5, zorder=3,
                  solid_capstyle="round")
    ax_a.plot([], [], color="#1b3a6b", lw=1.5, label=r"$\{G=0\}$")
    ax_a.plot([], [], color="0.78", lw=0.4, label=r"contours of $R$")

    gR = np.array(FT._centralGradient(FT.removalR, us, as_, p))
    gG = np.array(FT._centralGradient(FT.aggregateG, us, as_, p))
    # the two gradients are PARALLEL, so drawn from a common origin one hides the
    # other and the panel shows nothing. offset them along the shared normal so
    # the reader sees two arrows pointing the same way, which is the whole point.
    d = gR / np.linalg.norm(gR)
    nrm = np.array([-d[1], d[0]])
    L, off = 0.34, 0.035
    for vec, sgn, col, lab in ((gR, +1.0, "#b3341f", r"$\nabla R$"),
                               (gG, -1.0, "#0f7a5a", r"$\nabla G$")):
        v = vec / np.linalg.norm(vec) * L
        base = np.array([us, as_]) + sgn * off * nrm
        ax_a.annotate("", xy=base + v, xytext=base,
                      arrowprops=dict(arrowstyle="-|>", color=col, lw=1.25,
                                      shrinkA=0, shrinkB=0), zorder=5)
        ax_a.plot([], [], color=col, lw=1.25, label=lab)
    ax_a.plot([us], [as_], "o", ms=3.6, mfc="w", mec="k", mew=0.9, zorder=6)
    ax_a.annotate(rf"$j_{{\mathrm{{crit}}}}=R(u^*,a^*)={jc:.4f}$",
                  xy=(us, as_), xytext=(0.06, 0.62),
                  fontsize=6.0,
                  arrowprops=dict(arrowstyle="-", lw=0.4, color="0.35"))
    ax_a.set_xlim(0, u_max)
    ax_a.set_ylim(0, a_max)
    ax_a.set_xlabel(r"soluble misfolded burden $u$")
    ax_a.set_ylabel(r"aggregate burden $a$")
    ax_a.set_title("(a) the gradients align at the fold", loc="left")
    ax_a.legend(loc="upper right", handlelength=1.5, borderpad=0.2)

    # ---------------- panel (b): both coordinates, all computed ----------
    S = stats[np.argsort(stats[:, 2])]
    stable = S[S[:, 4] > 0]
    unstab = S[S[:, 4] < 0]
    C_A, C_U = "#1b3a6b", "#b3341f"
    for arr, ls, lw in ((stable, "-", 1.4), (unstab, (0, (3, 2)), 1.1)):
        if not len(arr):
            continue
        ax_b.plot(arr[:, 2], arr[:, 1], color=C_A, ls=ls, lw=lw)
        ax_b.plot(arr[:, 2], arr[:, 0], color=C_U, ls=ls, lw=lw)
    # label the curves inline: a legend would have to sit on top of one of them
    ax_b.plot([jc], [as_], "o", ms=3.4, mfc="w", mec=C_A, mew=0.9, zorder=6)
    ax_b.plot([jc], [us], "o", ms=3.4, mfc="w", mec=C_U, mew=0.9, zorder=6)
    ax_b.text(jc * 1.02, us, r"$u^*$", color=C_U, fontsize=7, va="center")
    ax_b.text(jc * 1.02, as_, r"$a^*$", color=C_A, fontsize=7, va="center")
    ax_b.text(0.004, 0.455, "dashed: unstable", fontsize=5.4,
              color="0.35", va="top")
    ax_b.set_xlim(0, jc * 1.10)
    ax_b.set_ylim(0, 0.62)
    ax_b.set_xlabel(r"damage influx $j$")
    ax_b.set_ylabel(r"equilibrium burden")
    ax_b.set_title("(b) two singularities, in order", loc="left")

    # the turn and the fold differ by ~0.1 % in j, so on the full axis they are
    # one point. this inset is a ZOOM on computed data, not a schematic.
    if turn is not None:
        ins = ax_b.inset_axes([0.30, 0.52, 0.42, 0.40])
        lo_j = jc * 0.9955
        for arr, ls, lw in ((stable[stable[:, 2] > lo_j], "-", 1.2),
                            (unstab[unstab[:, 2] > lo_j], (0, (2, 1.6)), 1.0)):
            if len(arr):
                ins.plot(arr[:, 2], arr[:, 1], color=C_A, ls=ls, lw=lw)
                ins.plot(arr[:, 2], arr[:, 0], color=C_U, ls=ls, lw=lw)
        ins.axvline(turn["j"], color=C_U, lw=0.5, ls=(0, (1, 2)))
        ins.axvline(jc, color="0.45", lw=0.5, ls=(0, (1, 2)))
        ins.plot([turn["j"]], [turn["u"]], "s", ms=3.0, mfc="w", mec=C_U, mew=0.8,
                 zorder=6)
        ins.plot([jc], [us], "o", ms=2.8, mfc="w", mec=C_U, mew=0.8, zorder=6)
        ins.plot([jc], [as_], "o", ms=2.8, mfc="w", mec=C_A, mew=0.8, zorder=6)
        ins.annotate(r"$G_a=0$", xy=(turn["j"], turn["u"]),
                     xytext=(turn["j"] - 3.2e-4, turn["u"] - 0.075),
                     fontsize=5.2, color=C_U, ha="center",
                     arrowprops=dict(arrowstyle="-|>", lw=0.4, color=C_U))
        ins.set_xlim(lo_j, jc * 1.0009)
        ins.set_ylim(0.18, 0.46)
        ins.set_xticks([turn["j"], jc])
        ins.set_xticklabels([r"$j_{\rm turn}$", r"$j_{\rm crit}$"], fontsize=5.0)
        ins.set_yticks([])
        ins.tick_params(pad=1.0, length=1.8)
        for s in ins.spines.values():
            s.set_linewidth(0.4)
            s.set_color("0.55")
        ins.set_title(r"zoom $\times$%d near the boundary" % round(1 / 0.0055),
                      fontsize=4.8, pad=1.5)
        # mark the zoomed strip on the main axes without connector lines
        ax_b.axvspan(lo_j, jc, color="0.85", alpha=0.55, lw=0, zorder=0)

    fig.tight_layout(pad=0.35)
    F.widthCheck(fig, F.W_SINGLE)
    hashes = F.save(fig, FIGURE)
    plt.close(fig)

    return {"j_crit": jc, "u_star": us, "a_star": as_,
            "sin_angle": ident["sin_angle"],
            "n_nullcline_points": int(len(stats)),
            "n_segments": len(segs),
            "homeostasis": homeo, "turn": turn, "hashes": hashes}


if __name__ == "__main__":
    out = build()
    print(f"Figure {FIGURE[3:]}")
    print("  j_crit = %.6f at (u*, a*) = (%.6f, %.6f)"
          % (out["j_crit"], out["u_star"], out["a_star"]))
    print("  sin(grad R, grad G) at the plotted fold = %.3e" % out["sin_angle"])
    tn = out["turn"]
    if tn:
        print("  turn (G_a=0): j=%.6f (j/j_crit=%.4f) u=%.6f a=%.6f"
              % (tn["j"], tn["j"] / out["j_crit"], tn["u"], tn["a"]))
        print("    det J = %.6g  vs  R_a.G_u = %.6g   (stable: %s)"
              % (tn["det_J"], tn["Ra_times_Gu"], tn["det_J"] > 0))
    h = out["homeostasis"]
    print("  nullcline points: %d | G_u in [%.4g, %.4g] | horizontal tangent: %s"
          % (h["n_points"], h["min_G_u"], h["max_G_u"], h["exists"]))
    for k, v in out["hashes"].items():
        print("  %-10s %s" % (k, v[:16]))
