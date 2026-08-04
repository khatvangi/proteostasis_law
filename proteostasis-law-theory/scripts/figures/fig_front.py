"""Figure 3 -- the strategy front, and where along it the boundary actually binds.

The claim §7 makes is deflationary and easy to lose: the throughput optimum sits
EXACTLY on the feasibility boundary, and yet most of the front does not. A figure
showing only the optimum would restate the flattering half. So the panel carries
the front coloured by `j/j_crit`, and the inset plots `j/j_crit` along it, where
the 0.227-0.965 range is as legible as the optimum is.

The optimum is COMPUTED, not annotated from the text. The 469 feasible strategies
come from `data/figures/pareto.tsv`; the exact optimum comes from
`pareto.optimumOnBoundary()`, which parametrises the boundary rather than a grid.
Its `j/j_crit` of 1.000000 is not in the grid at all -- the best grid point reaches
0.8975, and that gap is discretisation, which is exactly why the grid cannot be
asked to supply this number.

Neither source touches the run root.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# manuscript figure number, by ORDER OF FIRST MENTION in bmb_v4.md.
# filenames are deliberately semantic so a reorder touches this line only.
FIGURE = "fig4"

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "figures",
           REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import _figstyle as F  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import pareto as P  # noqa: E402

DATA = REPO_ROOT / "data" / "figures" / "pareto.tsv"

CMAP = "viridis"


def build():
    F.setStyle()
    df = pd.read_csv(DATA, sep="\t")
    front = df[df["on_front"]].sort_values("throughput").reset_index(drop=True)

    # the exact optimum, solved on the boundary rather than read off the grid
    exact = P.optimumOnBoundary()
    if exact is None:
        raise RuntimeError("boundary optimisation returned nothing")
    exact_ratio = exact["j"] / exact["j_crit"]
    grid_best = df.loc[df["throughput"].idxmax()]

    fig, ax = plt.subplots(figsize=(F.W_DOUBLE, 0.80 * F.W_DOUBLE))

    # every feasible strategy, so the front reads as a boundary of something
    ax.scatter(df["throughput"], df["accuracy"], s=2.2, c="0.80", lw=0,
               zorder=1, rasterized=False)

    # the front, coloured by how close it runs to the feasibility boundary.
    # the connector is deliberately GREY and not colour-mapped: colour along a
    # segment would assert `j/j_crit` at points between the 13 solved strategies,
    # and nothing was solved there. the line is a reading guide; the markers are
    # the data.
    ax.plot(front["throughput"], front["accuracy"], "-", lw=0.7, color="0.55",
            zorder=3)
    sc = ax.scatter(front["throughput"], front["accuracy"], s=13,
                    c=front["j_over_jcrit"], cmap=CMAP,
                    norm=plt.Normalize(0.0, 1.0),
                    edgecolors="w", linewidths=0.35, zorder=4)

    # the optimum: on the boundary, and off the grid
    ax.plot([exact["throughput"]], [-exact["eps"]], marker="*", ms=8.5,
            mfc="#d62728", mec="w", mew=0.5, zorder=6, ls="none")
    ax.annotate(f"throughput optimum\n$j/j_{{\\mathrm{{crit}}}} = {exact_ratio:.6f}$",
                xy=(exact["throughput"], -exact["eps"]),
                xytext=(exact["throughput"] - 0.115, -exact["eps"] - 0.115),
                fontsize=5.2, color="#a01d24", ha="left", va="top",
                arrowprops=dict(arrowstyle="-", lw=0.5, color="#a01d24",
                                shrinkA=0.5, shrinkB=3.0))

    ax.set_xlabel("throughput")
    ax.set_ylabel("accuracy  ($-\\epsilon$)")
    ax.set_xlim(-0.02, 0.42)
    ax.set_ylim(-1.06, 0.06)

    cb = fig.colorbar(sc, ax=ax, pad=0.02, fraction=0.045)
    cb.set_label("$j/j_{\\mathrm{crit}}$ along the front", fontsize=5.6,
                 labelpad=2.0)
    cb.ax.tick_params(labelsize=5.0, length=1.8, pad=1.2)
    cb.outline.set_linewidth(0.4)

    # inset: the deflationary half, at the same visual weight as the optimum
    lo = float(front["j_over_jcrit"].min())
    hi = float(front["j_over_jcrit"].max())
    ins = ax.inset_axes([0.14, 0.14, 0.40, 0.33])
    ins.plot(front["throughput"], front["j_over_jcrit"], "-", lw=0.8,
             color="0.35", zorder=2)
    ins.scatter(front["throughput"], front["j_over_jcrit"], s=6,
                c=front["j_over_jcrit"], cmap=CMAP,
                norm=plt.Normalize(0.0, 1.0), lw=0, zorder=3)
    ins.axhline(1.0, color="#a01d24", lw=0.7, ls=(0, (3, 2)), zorder=1)
    ins.plot([exact["throughput"]], [exact_ratio], marker="*", ms=5.0,
             mfc="#d62728", mec="none", ls="none", zorder=4)
    ins.set_ylim(0.0, 1.18)
    ins.set_xlabel("throughput", fontsize=4.6, labelpad=1.0)
    ins.set_ylabel("$j/j_{\\mathrm{crit}}$", fontsize=4.6, labelpad=1.0)
    ins.tick_params(labelsize=4.0, pad=1.0, length=1.8)
    for s in ins.spines.values():
        s.set_linewidth(0.4)
        s.set_color("0.55")
    ins.set_title(f"the front spans {lo:.3f}–{hi:.3f}:\nonly its end presses on "
                  "the boundary", fontsize=4.4, pad=1.5, color="0.35")

    fig.tight_layout(pad=0.35)
    F.widthCheck(fig, F.W_DOUBLE)
    hashes = F.save(fig, FIGURE)
    plt.close(fig)
    return {"n_feasible": len(df), "n_front": len(front),
            "front_lo": lo, "front_hi": hi,
            "front_median": float(front["j_over_jcrit"].median()),
            "exact_ratio": float(exact_ratio),
            "exact_throughput": float(exact["throughput"]),
            "grid_best_ratio": float(grid_best["j_over_jcrit"]),
            "grid_best_throughput": float(grid_best["throughput"]),
            "hashes": hashes}


if __name__ == "__main__":
    o = build()
    print(f"Figure {FIGURE[3:]}")
    print("  %d feasible strategies, %d non-dominated"
          % (o["n_feasible"], o["n_front"]))
    print("  j/j_crit along the front : %.4f - %.4f (median %.4f)"
          % (o["front_lo"], o["front_hi"], o["front_median"]))
    print("  exact optimum            : T=%.6f, j/j_crit=%.6f"
          % (o["exact_throughput"], o["exact_ratio"]))
    print("  best grid strategy       : T=%.6f, j/j_crit=%.6f  <- discretisation"
          % (o["grid_best_throughput"], o["grid_best_ratio"]))
    for k, v in o["hashes"].items():
        print("  %-10s %s" % (k, v[:16]))
