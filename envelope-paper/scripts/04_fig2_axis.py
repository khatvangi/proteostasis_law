#!/usr/bin/env python3
"""
Figure 2 -- operational structure of the code in (mu, nu) space.

panel a: mu axis   -- synonyms cluster far more tightly than permutation nulls
panel b: nu axis   -- no detectable structure
panel c: how much of the variation in log mu is between amino acids

panel c replaces the old "combined (mu, nu)" panel. the combined statistic is a
weighted mixture of a strong signal and a null, so it carries no information the
first two panels do not already carry, and its value moves with the mixing
weights. the variance decomposition is the honest summary of the mu result: it
states the size of the amino-acid-level structure without attributing it to
selection, which the permutation test cannot establish.
"""
from pathlib import Path

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import figstyle as fs

ROOT = Path(__file__).resolve().parent.parent
COMP = ROOT / "data" / "computed"
FIGS = ROOT / "figures"


def main():
    fs.setup()
    FIGS.mkdir(parents=True, exist_ok=True)

    tests = pd.read_csv(COMP / "axis_tests.tsv", sep="\t")
    tests = tests[tests.null == "within_degeneracy"].set_index("axis")
    axes_df = pd.read_csv(COMP / "codon_axes.tsv", sep="\t")

    fig, axs = plt.subplots(1, 3, figsize=(fs.COL_DOUBLE, 2.9))

    # ---------------- panels a, b: nulls vs observed ----------------
    spec = [
        ("mu", "a", r"operational spread $\Delta$" "\n" r"$\mu$ (mistranslation)",
         "Synonyms cluster in error rate", "right"),
        ("nu", "b", r"operational spread $\Delta$" "\n" r"$\nu$ (supply, tAI)",
         "No structure on supply", "left"),
    ]
    for (axis, letter, xlab, title, box_side), ax in zip(spec, axs):
        null = pd.read_csv(COMP / f"null_{axis}.tsv", sep="\t")
        row = tests.loc[axis]

        sns.histplot(data=null, x="delta", bins=50, ax=ax,
                     color=fs.C["neutral"], edgecolor="white", linewidth=0.2)
        ax.axvline(row.null_mean, color=fs.C["muted"], linestyle="--",
                   linewidth=0.8)
        ax.axvline(row.observed, color=fs.C["alert"], linewidth=2.0, zorder=5)

        sig = ("***" if row.p_one_sided < 0.001 else
               "**" if row.p_one_sided < 0.01 else
               "*" if row.p_one_sided < 0.05 else "n.s.")
        ax.text(0.04 if box_side == "left" else 0.96, 0.96,
                f"observed\nz = {row.z:+.2f}\np = {row.p_one_sided:.4f} {sig}",
                transform=ax.transAxes, va="top", ha=box_side,
                fontsize=6.5, family="monospace",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor=fs.C["muted"], alpha=0.9))

        ax.set_xlabel(xlab)
        ax.set_ylabel("null codes (10,000)")
        ax.set_title(title)
        ax.set_yticks([])
        fs.panel_label(ax, letter)

    # ---------------- panel c: where log mu variation sits ----------------
    ax = axs[2]
    vd = pd.read_json(COMP / "mu_variance_decomposition.json", typ="series")
    order = axes_df.groupby("aa").log_mu.mean().sort_values().index.tolist()

    sns.stripplot(data=axes_df, x="log_mu", y="aa", order=order, ax=ax,
                  color=fs.C["burden"], size=3.0, alpha=0.8, jitter=0.10)
    means = axes_df.groupby("aa").log_mu.mean().reindex(order)
    ax.scatter(means.to_numpy(), range(len(order)), marker="|", s=150,
               color=fs.C["alert"], linewidths=1.3, zorder=6)
    # no legend: at this panel size the legend marker lands inside the data and
    # reads as a extra point. the red ticks are identified in the figure legend.

    # eta^2 goes in the title, not a floating box: at this panel size every
    # corner overlaps data
    ax.set_title("Most variation is between amino\n"
                 rf"acids, not within  ($\eta^2 = "
                 rf"{vd['eta_squared_log_mu_between_aa']:.2f}$)")
    ax.set_xlabel(r"$\log\,\mu$  (per codon)")
    ax.set_ylabel("amino acid")
    ax.tick_params(axis="y", labelsize=6.0, pad=1)
    fs.panel_label(ax, "c", dx=-0.20)

    sns.despine(fig=fig)
    fig.tight_layout(w_pad=1.6)
    fs.save(fig, str(FIGS / "Fig2_axis_structure"))

    for axis in ("mu", "nu", "2D"):
        r = tests.loc[axis]
        print(f"  {axis:>3}: observed={r.observed:.4f} z={r.z:+.2f} "
              f"p={r.p_one_sided:.4f} ({r.direction})")


if __name__ == "__main__":
    main()
