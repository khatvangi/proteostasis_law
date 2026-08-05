"""Figure 4 -- the prediction, as a beta-indexed family.

The requirement on the old-pole aggregate load scales as `1/beta`, where `beta` is
the share of total aggregate held in the visible inclusion body. It is drawn as a
BAND rather than a curve because the match band itself is an interval, and on a
log axis the `1/beta` scaling is a straight line.

Two things this figure deliberately does NOT do:

  * it does not draw a vertical line at `beta = 0.145`. That figure combines
    Winkler's heat-stress molecule counts with a TWO-foci assumption, and Lindner's
    unstressed cells carry one focus in 46.5 % of cases and two in 1.2 %. The
    assumption does not transfer, so no lower bound on `beta` is defensible and
    the axis runs to the left edge without a stop (D033).
  * it does not present the wild-type aggregate fraction as a value. Tomoyasu et
    al. report it as undetected at 30 C, which is a bound.

Needs no run root.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# manuscript figure number, by ORDER OF FIRST MENTION in bmb_v4.md.
# filenames are deliberately semantic so a reorder touches this line only.
FIGURE = "fig5"

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3",
           REPO_ROOT / "scripts" / "figures"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import _figstyle as F  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import asymmetric_division as A  # noqa: E402

# damping is BETA-DEPENDENT and measured per beta (D033). An earlier version of
# this figure used a single 0.35, which put every marked value ~1.5% off the
# section 8.4 table -- a caption-against-text divergence with no cause but
# convenience. `A.dampingAtBeta` interpolates the measured points, so the three
# marked beta reproduce the table exactly.
damping = A.dampingAtBeta

# ONE range, owned here, used by the sweep AND by the manuscript's table rows.
#
# The table used to stop at beta = 0.25 while the figure plotted to 0.05. Prose
# written from the table then put the closest approach to the measured load at
# 16x, when over the figure's own range it is 3.19x -- the paper's only
# falsifiable prediction, understated fivefold. That is a structural fault, not
# an oversight: two ranges for one quantity will always let prose land in the
# wrong place. `tableRows()` below emits exactly the rows the manuscript prints,
# from the same endpoints the figure draws, so either source gives one answer.
BETA_LO, BETA_HI = 0.05, 1.0
TABLE_BETAS = (1.00, 0.75, 0.50, 0.25, 0.05)


def tableRows():
    """the Section 8.3 table, generated. one row per beta over the plotted range.

    `f_eff = (1 + beta)/2` is the partitioning fraction. The ratio columns pair
    the EXTREMES -- the smallest ratio is the measured lower bound against the
    largest requirement -- because that is the honest span, and pairing
    like-with-like would understate the closest approach.
    """
    rpo_lo, rpo_hi = A.TOMOYASU2001["rpoH_null_30C"]
    rows = []
    for b in TABLE_BETAS:
        lo, hi = A.requiredAggregateFractionBeta(b, damping(b))
        rows.append({
            "beta": b,
            "f_eff": (1.0 + b) / 2.0,
            "pct_lo": 100.0 * lo,
            "pct_hi": 100.0 * hi,
            "ratio_lo": rpo_lo / hi,
            "ratio_hi": rpo_hi / lo,
        })
    return rows


def build():
    F.setStyle()
    betas = np.geomspace(BETA_LO, BETA_HI, 200)
    lo = np.array([A.requiredAggregateFractionBeta(b, damping(b))[0] for b in betas])
    hi = np.array([A.requiredAggregateFractionBeta(b, damping(b))[1] for b in betas])

    fig, ax = plt.subplots(figsize=(F.W_DOUBLE, 0.80 * F.W_DOUBLE))

    rpo_lo, rpo_hi = A.TOMOYASU2001["rpoH_null_30C"]
    ax.axhspan(rpo_lo, 1.0, color="#c9d6e4", alpha=0.55, lw=0, zorder=0)
    ax.text(0.98, rpo_lo * 1.5, "aggregate load measured in $\\Delta rpoH$\n"
                                "(5–10 %); wild type reported only as\n"
                                "undetected, i.e. an upper bound",
            fontsize=4.8, ha="right", va="bottom", color="#33526f")

    ax.fill_between(betas, lo, hi, color="#b3341f", alpha=0.30, lw=0,
                    label="required, from the measured band")
    ax.plot(betas, lo, color="#b3341f", lw=0.9)
    ax.plot(betas, hi, color="#b3341f", lw=0.9)

    for b in (1.0, 0.5, 0.25):
        l, h = A.requiredAggregateFractionBeta(b, damping(b))
        ax.plot([b, b], [l, h], color="#7a1f0f", lw=1.6, solid_capstyle="butt",
                zorder=4)
        ax.annotate(f"{100*l:.2f}–{100*h:.2f} %", xy=(b, h),
                    xytext=(b, h * 1.9), fontsize=4.8, ha="center",
                    color="#7a1f0f")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(0.05, 1.06)
    ax.set_ylim(2e-4, 1.0)
    ax.set_xlabel(r"$\beta$, share of total aggregate held in the visible focus")
    ax.set_ylabel("old-pole aggregate required,\nas a fraction of total protein")
    ax.set_title(r"requirement $\propto 1/\beta$; no measurement bounds $\beta$",
                 loc="left", fontsize=6.4)
    ax.legend(loc="lower left", handlelength=1.6, borderpad=0.25)
    ax.grid(True, which="major", color="0.9", lw=0.35)
    ax.set_axisbelow(True)

    fig.tight_layout(pad=0.35)
    F.widthCheck(fig, F.W_DOUBLE)
    hashes = F.save(fig, FIGURE)
    plt.close(fig)
    # the closest and widest separations from the only measured load, over the
    # WHOLE plotted range rather than at the marked beta -- the prose quotes these
    return {"damping_by_beta": A.DAMPING_BY_BETA,
            "closest_ratio": float((rpo_lo / hi).min()),
            "widest_ratio": float((rpo_hi / lo).max()),
            "closest_at_beta": float(betas[(rpo_lo / hi).argmin()]),
            "at_beta_1": A.requiredAggregateFractionBeta(1.0, damping(1.0)),
            "at_beta_05": A.requiredAggregateFractionBeta(0.5, damping(0.5)),
            "at_beta_025": A.requiredAggregateFractionBeta(0.25, damping(0.25)),
            "rows": tableRows(),
            "rpoH": (rpo_lo, rpo_hi), "hashes": hashes}


if __name__ == "__main__":
    o = build()
    print(f"Figure {FIGURE[3:]}")
    print("  beta=1.00 -> %.4f%% - %.4f%% of proteome"
          % (100 * o["at_beta_1"][0], 100 * o["at_beta_1"][1]))
    print("  beta=0.25 -> %.4f%% - %.4f%% of proteome"
          % (100 * o["at_beta_025"][0], 100 * o["at_beta_025"][1]))
    for k, v in o["hashes"].items():
        print("  %-10s %s" % (k, v[:16]))
