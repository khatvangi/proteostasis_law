#!/usr/bin/env python3
"""
Figure 4 -- the distinguishing prediction, tested in the model.

(a) margin ladder: observed joint damage against the additive expectation, as the
    starting margin is compressed. the two coincide at wild-type margin and
    separate only as the margin closes; past a point the combination collapses
    while each single perturbation is still survivable.
(b) interaction as a percentage of the additive expectation over the
    (error, capacity) perturbation grid at the observed rate. the greyed region is
    joint collapse, where the interaction is unbounded rather than large.
(c) the two capacity knobs. C_tot is nearly inert because the folding arm is
    97.9% saturated at the observed operating point; k_obs_max acts
    proportionally. the gap is a property of the parameterization, not of the
    framework.
"""
import json
from pathlib import Path

import numpy as np
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

    s1 = pd.read_csv(COMP / "supraadditivity_margin_sweep.tsv", sep="\t")
    s2 = pd.read_csv(COMP / "supraadditivity_effect_grid.tsv", sep="\t")
    s4 = pd.read_csv(COMP / "supraadditivity_knob_comparison.tsv", sep="\t")
    summ = json.loads((COMP / "supraadditivity_summary.json").read_text())

    fig, axs = plt.subplots(1, 3, figsize=(fs.COL_DOUBLE, 2.9))

    # ---------------- (a) margin ladder ----------------
    ax = axs[0]
    d = s1[np.isfinite(s1.D_error)].copy()
    d["additive"] = d.D_error + d.D_capacity
    d = d.sort_values("margin_baseline")

    sns.lineplot(data=d, x="margin_baseline", y="additive", ax=ax, marker="o",
                 ms=4.5, color=fs.C["muted"], linewidth=1.4)
    viable = d[~d.collapsed_both]
    sns.lineplot(data=viable, x="margin_baseline", y="D_both", ax=ax, marker="s",
                 ms=4.5, color=fs.C["burden"], linewidth=1.6)

    coll = d[d.collapsed_both]
    if len(coll):
        # collapse means unbounded margin loss, so these points get NO y-value --
        # drawing them on the additive line would imply a number they do not have.
        # they are marked at the top of the axis with an upward arrow instead.
        lo = min(viable.D_both.min(), d.additive.min())
        hi = max(viable.D_both.max(), d.additive.max())
        pad = 0.10 * (hi - lo)
        ax.set_ylim(lo - pad, hi + 3.2 * pad)
        y_top = ax.get_ylim()[1] - 0.5 * pad
        ax.scatter(coll.margin_baseline, [y_top] * len(coll), marker="X", s=52,
                   color=fs.C["alert"], zorder=8)
        for x in coll.margin_baseline:
            ax.annotate("", xy=(x, y_top + 0.42 * pad), xytext=(x, y_top - 0.7 * pad),
                        arrowprops=dict(arrowstyle="-|>", lw=0.9,
                                        color=fs.C["alert"]))
        ax.axvspan(coll.margin_baseline.max() + 0.155,
                   d.margin_baseline.min() - 0.05,
                   color=fs.C["alert"], alpha=0.09, lw=0, zorder=0)
        ax.text(coll.margin_baseline.min() - 0.02, y_top - 1.5 * pad,
                "each single perturbation\nsurvivable, combination\nlethal",
                ha="left", va="top", fontsize=5.8, color=fs.C["alert"])

    # mark the paper's own evaluation point -- Result 5 was previously anchored at
    # the window bottom, off the left of this ladder
    m_ec = summ["baseline_margin_log10"]
    ax.axvline(m_ec, color=fs.C["capacity"], linewidth=1.2, zorder=3)
    ax.annotate("E. coli\n(usage-wtd $\\bar\\mu$)", xy=(m_ec, ax.get_ylim()[0]),
                xytext=(3, 4), textcoords="offset points", fontsize=5.8,
                color=fs.C["capacity"], ha="left", va="bottom")

    ax.invert_xaxis()          # margin closes to the right
    ax.set_xlabel("starting margin (log$_{10}$)\n← more headroom      less →")
    ax.set_ylabel("margin lost (log$_{10}$ units)")
    ax.set_title("E. coli sits where the synergy\nis already appreciable")
    # series labelled in place; an in-plot legend puts its markers among the data
    a0 = d.loc[d.margin_baseline.idxmax()]
    ax.annotate("additive\nexpectation", xy=(a0.margin_baseline, a0.additive),
                xytext=(4, 10), textcoords="offset points", fontsize=5.8,
                color=fs.C["muted"], ha="left")
    v0 = viable.loc[viable.margin_baseline.idxmin()]
    ax.annotate("observed,\nboth together", xy=(v0.margin_baseline, v0.D_both),
                xytext=(6, -4), textcoords="offset points", fontsize=5.8,
                color=fs.C["burden"], ha="left")
    fs.panel_label(ax, "a")

    # ---------------- (b) effect-size grid ----------------
    ax = axs[1]
    # reindex to the full factor set: pivot_table drops rows whose interaction is
    # entirely NaN, and at this evaluation point the strongest perturbation rows
    # are fully collapsed -- dropping them would hide a third of the grid, which
    # is exactly what this panel is claiming
    efs = sorted(s2.error_factor.unique(), reverse=True)
    cfs = sorted(s2.capacity_factor.unique())
    piv = (s2.pivot_table(index="error_factor", columns="capacity_factor",
                          values="interaction_pct_of_additive")
             .reindex(index=efs, columns=cfs))
    collapsed = (s2.pivot_table(index="error_factor", columns="capacity_factor",
                                values="collapsed_both")
                   .reindex(index=efs, columns=cfs).fillna(False).astype(bool))
    # two different kinds of collapse, which must not look alike: cells where the
    # combination kills a system each perturbation alone leaves viable (the
    # prediction), versus cells where one perturbation alone already kills it
    synth = (s2.pivot_table(index="error_factor", columns="capacity_factor",
                            values="qualitative_supraadditive")
               .reindex(index=efs, columns=cfs).fillna(False).astype(bool))

    sns.heatmap(piv, ax=ax, cmap="rocket_r", vmin=0, annot=True, fmt=".1f",
                annot_kws={"fontsize": 5.6}, linewidths=0.4, linecolor="white",
                cbar_kws={"label": "% above additive", "pad": 0.02})
    # grey out the joint-collapse cells: there the interaction is unbounded
    for i, e in enumerate(piv.index):
        for j, c in enumerate(piv.columns):
            if not bool(collapsed.loc[e, c]):
                continue
            if bool(synth.loc[e, c]):
                ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor=fs.C["alert"],
                                           alpha=0.55, edgecolor="white",
                                           lw=0.4, zorder=3))
                ax.text(j + 0.5, i + 0.5, "SL", ha="center", va="center",
                        fontsize=5.4, color="white", fontweight="bold", zorder=4)
            else:
                ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="#d8d8d8",
                                           edgecolor="white", lw=0.4, zorder=3))
                ax.text(j + 0.5, i + 0.5, "1×", ha="center", va="center",
                        fontsize=5.0, color="#555", zorder=4)
    ax.set_xlabel("rescue throughput knocked down (÷)")  # at the consistent f
    ax.set_ylabel("error rate raised (×)")
    n_sl = int(synth.to_numpy().sum())
    ax.set_title(f"{n_sl} of {synth.size} combinations are\n"
                 "survivable alone, lethal together")
    ax.tick_params(labelsize=6)
    fs.panel_label(ax, "b", dx=-0.20)

    # ---------------- (c) knob comparison ----------------
    ax = axs[2]
    k = s4[np.isfinite(s4.D_capacity)].copy()
    k["knob_label"] = k.knob.map({"k_obs_max": "throughput\n$k_{obs,max}$",
                                  "C_tot": "pool size\n$C_{tot}$"})
    sns.barplot(data=k[k.f_base == k.f_base.min()], x="knob_label", y="D_capacity",
                ax=ax, hue="knob_label", legend=False,
                palette=[fs.C["capacity"], fs.C["neutral"]],
                edgecolor="white", linewidth=0.8)
    dref = k[k.f_base == k.f_base.min()].D_error.iloc[0]
    ax.axhline(dref, color=fs.C["alert"], linestyle="--", linewidth=1.0)
    ax.text(0.02, dref, "error ×3, for scale ", transform=ax.get_yaxis_transform(),
            ha="left", va="bottom", fontsize=6.0, color=fs.C["alert"])
    for i, r in enumerate(k[k.f_base == k.f_base.min()].itertuples()):
        ax.text(i, r.D_capacity + 0.012, f"{r.D_capacity:.3f}", ha="center",
                va="bottom", fontsize=6.5)
    ax.set_xlabel("")
    ax.set_ylabel("margin lost, capacity ÷3 (log$_{10}$)")
    ax.set_title("The pool-size knob is inert:\nthe folding arm is 97% saturated")
    fs.panel_label(ax, "c", dx=-0.22)

    sns.despine(fig=fig)
    fig.tight_layout(w_pad=2.0)
    fs.save(fig, str(FIGS / "Fig4_supraadditivity"))

    wt = d.loc[d.margin_baseline.idxmax()]
    print(f"  interaction at wild-type margin ({wt.margin_baseline:.2f} log10): "
          f"{wt.interaction_pct_of_additive:+.1f}% of additive")
    if len(coll):
        print(f"  joint collapse first appears at starting margin "
              f"{coll.margin_baseline.max():.2f} log10 "
              f"(x{10**coll.margin_baseline.max():.0f} headroom)")


if __name__ == "__main__":
    main()
