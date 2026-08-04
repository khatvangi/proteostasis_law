"""Figure S1 -- the parallelism residual is bracket tolerance, not identity failure.

ONE PANEL, and it carries the POSITIVE claim: sin(angle) between the two gradients
tracks the recorded leading eigenvalue over four decades with a log-log correlation
of +0.996. The residual at the recorded states is not zero, and should not be,
because those states are bracketed approximations rather than exact folds. If the
identity were failing, the residual would not care how tight the bracket was.

The normalisation contrast belongs in the CAPTION, not in a second panel. It is a
methodological point about which statistic to quote, and giving it a panel would
put the paper's own retracted metric at the same visual weight as the result.

This script recomputes every number the caption uses and prints them, because the
caption's contrast previously had no generator -- see D041.

Reads `data/figures/identity.tsv`, never the run root.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# manuscript figure number, by ORDER OF FIRST MENTION in bmb_v4.md.
# filenames are deliberately semantic so a reorder touches this line only.
FIGURE = "figS1"

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "figures"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import _figstyle as F  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402


def captionNumbers(df: pd.DataFrame) -> dict:
    """every quantity the caption and section 5 quote, from one place.

    log-log is the reported convention because both quantities span decades and
    because the +0.996 statistic beside it is already defined that way. The raw
    Pearson is returned too, so a reader who correlates the columns directly finds
    their own number here rather than a discrepancy.
    """
    lg = np.log10
    out = {"n": len(df)}
    out["corr_parallelism"] = float(np.corrcoef(lg(df["sin_angle"]),
                                                lg(df["eig_abs"]))[0, 1])
    for key, col in (("max", "res_max_normalised"),
                     ("grad", "res_gradient_normalised")):
        v = df[col]
        out[f"{key}_median"] = float(v.median())
        out[f"{key}_p90"] = float(v.quantile(0.90))
        out[f"{key}_p99"] = float(v.quantile(0.99))
        out[f"{key}_max"] = float(v.max())
        out[f"{key}_corr_loglog"] = float(np.corrcoef(lg(v), lg(df["eig_abs"]))[0, 1])
        out[f"{key}_corr_raw"] = float(np.corrcoef(v, df["eig_abs"])[0, 1])
    # the same defect stated without a correlation: split at the median bracket
    m = df["eig_abs"].median()
    out["max_tighter_half"] = float(df.loc[df["eig_abs"] <= m, "res_max_normalised"].median())
    out["max_looser_half"] = float(df.loc[df["eig_abs"] > m, "res_max_normalised"].median())
    out["max_degradation"] = out["max_tighter_half"] / out["max_looser_half"]
    tightest = df.loc[df["eig_abs"].idxmin()]
    out["tightest_eig"] = float(tightest["eig_abs"])
    out["tightest_sin"] = float(tightest["sin_angle"])
    return out


def build():
    F.setStyle()
    df = pd.read_csv(REPO_ROOT / "data" / "figures" / "identity.tsv", sep="\t")
    c = captionNumbers(df)

    fig, ax = plt.subplots(figsize=(F.W_DOUBLE, 0.78 * F.W_DOUBLE))

    ax.scatter(df["eig_abs"], df["sin_angle"], s=5.0, c="#1b3a6b", alpha=0.45,
               lw=0, zorder=3)

    # the fitted line, in log space, drawn only across the data's own range
    lg = np.log10
    sl, ic = np.polyfit(lg(df["eig_abs"]), lg(df["sin_angle"]), 1)
    xs = np.array([df["eig_abs"].min(), df["eig_abs"].max()])
    ax.plot(xs, 10 ** (ic + sl * lg(xs)), "-", lw=0.8, color="#b3341f", zorder=4)

    # the tightest bracket in the population: the point that makes the argument
    ax.plot([c["tightest_eig"]], [c["tightest_sin"]], "o", ms=4.5, mfc="none",
            mec="#b3341f", mew=0.8, zorder=5)
    ax.annotate("tightest bracket in the full population:\n"
                f"$|\\lambda| = {c['tightest_eig']:.2e}$,  "
                f"$\\sin\\theta = {c['tightest_sin']:.2e}$",
                xy=(c["tightest_eig"], c["tightest_sin"]),
                xytext=(0.05, 0.90), textcoords="axes fraction",
                fontsize=5.2, color="#8d2a19", ha="left", va="top",
                arrowprops=dict(arrowstyle="-", lw=0.5, color="#8d2a19",
                                shrinkA=1.0, shrinkB=3.0))

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("$|\\lambda|$ at the recorded state  (bracket tightness)")
    ax.set_ylabel("$\\sin\\theta$ between $\\nabla R$ and $\\nabla G$")
    ax.set_title(f"$r = {c['corr_parallelism']:+.4f}$ over all "
                 f"{c['n']} folds, slope {sl:.2f}",
                 loc="left", fontsize=6.2, color="0.30")
    ax.grid(True, which="major", color="0.88", lw=0.4, zorder=0)

    fig.tight_layout(pad=0.35)
    F.widthCheck(fig, F.W_DOUBLE)
    c["slope"] = float(sl)
    c["hashes"] = F.save(fig, FIGURE)
    plt.close(fig)
    return c


if __name__ == "__main__":
    o = build()
    print(f"Figure {FIGURE[3:]}  (all %d folds of the load grid, no subsample)" % o["n"])
    print("  corr(log sin, log |eig|)      : %+.4f   slope %.3f"
          % (o["corr_parallelism"], o["slope"]))
    print("  tightest bracket              : |eig| %.2e -> sin %.2e"
          % (o["tightest_eig"], o["tightest_sin"]))
    print("  -- the caption's contrast, recomputed (D041) --")
    for key, lab in (("max", "max-normalised"), ("grad", "gradient-normalised")):
        print("  %-20s median %.3e  p90 %.3e  p99 %.3e  max %.3e"
              % (lab, o[f"{key}_median"], o[f"{key}_p90"], o[f"{key}_p99"],
                 o[f"{key}_max"]))
        print("  %-20s corr with |eig|: log-log %+.3f, raw Pearson %+.3f"
              % ("", o[f"{key}_corr_loglog"], o[f"{key}_corr_raw"]))
    print("  max-normalised degrades toward the fold: tighter half %.3e vs "
          "looser %.3e (%.2fx)"
          % (o["max_tighter_half"], o["max_looser_half"], o["max_degradation"]))
    for k, v in o["hashes"].items():
        print("  %-10s %s" % (k, v[:16]))
