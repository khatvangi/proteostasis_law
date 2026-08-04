"""Figure 2 -- where the boundary sits, and how widely.

The claim has two halves and the figure has to carry both: collapse happens FAR
BELOW saturation, and the spread is wide enough that one measurement discriminates
weakly. A bar chart of medians would show the first and hide the second, so the
distributions are drawn in full with the p5-p95 span marked explicitly.

Reads `data/figures/fig2_saturation.tsv`, never the run root.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "figures"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import _figstyle as F  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

DATA = REPO_ROOT / "data" / "figures" / "fig2_saturation.tsv"

# THE SCREEN THE STANDING LIMITATION ASKED FOR IS NOT APPLIED, AND THE DATA IS THE
# REASON. A screen is only defensible if the low-`s_a` draws form a distinct
# cluster. They do not: on a log axis the distribution runs smoothly from 1e-5 to
# 1 with no gap, and the median of the survivors slides continuously with the
# floor -- 0.090, 0.124, 0.162, 0.194, 0.227, 0.267, 0.309, 0.355 for floors of
# 1e-4 through 2e-2. Any floor is a free parameter that moves a load-bearing
# number by a factor of four, which is the degrees-of-freedom hazard this project
# keeps catching elsewhere. See D039.
#
# NO exclusion is applied either, and for the same discipline: the complete
# population is what §6 quotes, and dropping even the 47 numerically-zero draws
# shifts the medians to 0.1858 / 0.1602 / 0.0603 against the 0.175 / 0.155 / 0.056
# in the text. Introducing a second population to fix a cosmetic issue is how the
# 325-against-2884 confusion arose. The near-zero draws are noted in the caption
# and left in.
S_A_ZERO = 0.0


def build():
    F.setStyle()
    df = pd.read_csv(DATA, sep="\t")
    n_all = len(df)
    keep = df                       # complete population, no exclusion
    n_screened = n_all - len(keep)
    # the sensitivity that rules a screen out, recomputed here rather than quoted
    floors = (1e-4, 5e-4, 1e-3, 2e-3, 3e-3, 5e-3, 1e-2, 2e-2)
    sens = [(f, float(np.median(df.loc[df["s_a"] >= f, "s_a"]))) for f in floors]

    fig, ax = plt.subplots(figsize=(F.W_SINGLE, 0.36 * F.W_SINGLE))
    series = [("$s_{\\mathrm{ref}}$", keep["s_ref"].to_numpy(), "#1b3a6b"),
              ("$s_u$", keep["s_u"].to_numpy(), "#b3341f"),
              ("$s_a$", keep["s_a"].to_numpy(), "#0f7a5a")]

    stats = {}
    for i, (lab, v, col) in enumerate(series):
        y = len(series) - 1 - i
        # the full distribution, as a horizontal density strip
        h, edges = np.histogram(v, bins=60, range=(0.0, 1.0), density=True)
        h = h / h.max() * 0.34
        centres = 0.5 * (edges[:-1] + edges[1:])
        ax.fill_between(centres, y - h, y + h, color=col, alpha=0.30, lw=0)
        p5, p50, p95 = np.percentile(v, [5, 50, 95])
        ax.plot([p5, p95], [y, y], color=col, lw=1.1, solid_capstyle="butt")
        ax.plot([p50], [y], "o", ms=3.6, mfc="w", mec=col, mew=1.0, zorder=5)
        ax.text(-0.035, y, lab, ha="right", va="center", fontsize=7, color=col)
        ax.text(p95 + 0.02, y + 0.20,
                f"median {p50:.3f}   p5–p95 width {p95 - p5:.3f}",
                fontsize=5.2, color=col, va="center")
        stats[lab] = {"p5": p5, "median": p50, "p95": p95, "width": p95 - p5}

    ax.axvline(1.0, color="0.25", lw=0.9, ls=(0, (4, 2)))
    ax.text(1.0, len(series) - 0.42, "  capacity exhaustion\n  would predict $s = 1$",
            fontsize=5.4, va="top", ha="left", color="0.25")

    ax.set_xlim(-0.16, 1.28)
    ax.set_ylim(-0.62, len(series) - 0.30)
    ax.set_yticks([])
    ax.set_xticks([0.0, 0.25, 0.5, 0.75, 1.0])
    ax.set_xlabel("Michaelis saturation fraction at the collapse boundary")
    ax.spines["left"].set_visible(False)
    n_near_zero = int((df["s_a"] < 1e-9).sum())
    ax.set_title(f"all {n_all} folds, no screen and no exclusion "
                 f"({n_near_zero} sit at $s_a < 10^{{-9}}$ and are retained)",
                 loc="left", fontsize=5.8, color="0.35")

    # inset: why no screen. the median slides continuously with any floor chosen.
    ins = ax.inset_axes([0.60, 0.06, 0.28, 0.30])
    ins.plot([f for f, _ in sens], [m for _, m in sens], "o-", ms=2.0, lw=0.8,
             color="#0f7a5a")
    ins.set_xscale("log")
    ins.set_xlabel("screen floor on $s_a$", fontsize=4.6, labelpad=1.0)
    ins.set_ylabel("median $s_a$", fontsize=4.6, labelpad=1.0)
    ins.tick_params(labelsize=4.0, pad=1.0, length=1.8)
    for s in ins.spines.values():
        s.set_linewidth(0.4)
        s.set_color("0.55")
    ins.set_title("no natural break: any floor\nmoves the median", fontsize=4.4,
                  pad=1.5, color="0.35")

    fig.tight_layout(pad=0.35)
    F.widthCheck(fig, F.W_SINGLE)
    hashes = F.save(fig, "fig2")
    plt.close(fig)
    return {"n_all": n_all, "n_kept": len(keep), "n_screened": n_screened,
            "stats": stats, "sensitivity": sens, "hashes": hashes}


if __name__ == "__main__":
    o = build()
    print("Figure 2")
    print("  %d folds, %d excluded as numerically zero, %d plotted"
          % (o["n_all"], o["n_screened"], o["n_kept"]))
    print("  median s_a against screen floor (why no screen is applied):")
    for f, m in o["sensitivity"]:
        print("     floor %-8g -> median %.4f" % (f, m))
    for k, s in o["stats"].items():
        print("  %-16s median %.4f   p5-p95 %.4f-%.4f  width %.4f"
              % (k, s["median"], s["p5"], s["p95"], s["width"]))
    for k, v in o["hashes"].items():
        print("  %-10s %s" % (k, v[:16]))
