#!/usr/bin/env python3
"""
Figure 3 -- where the observed mistranslation rate sits inside the envelope.

three independent-ish estimates of the maximum tolerable per-codon error rate:
  observed    ~1e-4 to 1e-3 /codon   (literature consensus for E. coli)
  arithmetic  combinatorial bound: per-codon error x length x p_misfold,
              integrated over the E. coli proteome
  two-pool ODE  dynamical bound, set by the aggregated fraction reaching a
              lethal level (aggregation-death)

IMPORTANT -- this figure uses paired_mc_results.json, NOT the marginal MC runs.
the earlier draft compared the two bounds using independent Monte Carlo runs that
drew their SHARED parameters (N, p_misfold, S_syn) from different distributions.
proteostasis-P1/PAIRED_MC_TASK.md states plainly that the resulting statistic
"P(arith < ODE) = 0.654" is "not a valid paired statistic", and paired_mc.py was
written to replace it. the valid paired value is P(r < 1) = 0.768. the earlier
draft and figure carried the retracted number and the superseded backup run
(ODE MC median 1.7e-2 instead of the current 2.6e-2). both are corrected here.

panel a: paired marginal distributions of the two bounds vs the observed band
panel b: the paired ratio r = f_arith / f_ODE, with P(r < 1)
panel c: headroom at the observed rate, in both pools of the two-pool model
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import figstyle as fs

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
COMP = ROOT / "data" / "computed"
FIGS = ROOT / "figures"

OBS_LO, OBS_HI = 1e-4, 1e-3      # observed E. coli mistranslation rate window


def main():
    fs.setup()
    FIGS.mkdir(parents=True, exist_ok=True)

    paired = json.loads((RAW / "paired_mc_results.json").read_text())
    arith = json.loads((RAW / "arithmetic_results.json").read_text())
    twopool = json.loads((RAW / "two_pool_results.json").read_text())

    s = paired["samples"]
    f_arith = np.asarray(s["f_arith"], float)
    f_ode = np.asarray(s["f_ODE"], float)
    r = np.asarray(s["r"], float)

    det_arith = arith["A_reproduction"]["f_codon_crit_exact"]
    det_ode = twopool["B_compare"]["two_pool"]["f_codon_crit"]
    cross = twopool["E_ecoli_crosscheck"]

    fig, axs = plt.subplots(1, 3, figsize=(fs.COL_DOUBLE, 2.75))

    # ---------------- panel a: the two bounds vs observed ----------------
    ax = axs[0]
    long = pd.concat([
        pd.DataFrame({"bound": "arithmetic", "log10_f": np.log10(f_arith)}),
        pd.DataFrame({"bound": "two-pool ODE", "log10_f": np.log10(f_ode)}),
    ])
    sns.violinplot(data=long, x="bound", y="log10_f", ax=ax, hue="bound",
                   legend=False, cut=0, inner="quartile", linewidth=0.8,
                   palette={"arithmetic": fs.C["burden"],
                            "two-pool ODE": fs.C["ode"]}, alpha=0.75)

    ax.axhspan(np.log10(OBS_LO), np.log10(OBS_HI),
               color=fs.C["alert"], alpha=0.13, lw=0, zorder=0)
    ax.text(0.62, np.log10(OBS_LO) + 0.08, "observed E. coli rate",
            transform=ax.get_yaxis_transform(), ha="left", va="bottom",
            fontsize=6.2, color=fs.C["alert"])

    ax.scatter([0], [np.log10(det_arith)], marker="D", s=22, zorder=8,
               color=fs.C["burden"], edgecolor="white", linewidth=0.6)
    ax.scatter([1], [np.log10(det_ode)], marker="D", s=22, zorder=8,
               color=fs.C["ode"], edgecolor="white", linewidth=0.6)

    ax.set_ylabel(r"$\log_{10}$ tolerable rate  (/codon)")
    ax.set_xlabel("")
    ax.set_title("Bound medians lie 1-2 orders\nabove the observed rate")
    fs.panel_label(ax, "a")

    # ---------------- panel b: paired ratio ----------------
    ax = axs[1]
    p_arith_tighter = float((r < 1).mean())
    sns.histplot(x=np.log10(r), bins=55, ax=ax, color=fs.C["neutral"],
                 edgecolor="white", linewidth=0.2)
    ax.axvline(0, color=fs.C["alert"], linewidth=1.4, zorder=5)
    ax.axvline(np.log10(np.median(r)), color=fs.C["muted"], linestyle="--",
               linewidth=0.9, zorder=4)
    ax.text(0.04, 0.96,
            f"median r = {np.median(r):.2f}\n"
            f"P(r < 1) = {p_arith_tighter:.3f}",
            transform=ax.transAxes, va="top", ha="left",
            fontsize=6.5, family="monospace",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                      edgecolor=fs.C["muted"], alpha=0.9))
    ax.set_xlabel(r"$\log_{10}\,r$,   $r = f_{\rm arith}/f_{\rm ODE}$")
    ax.set_ylabel("paired MC draws (5,000)")
    ax.set_title("Paired: arithmetic is the\ntighter bound more often")
    ax.set_yticks([])
    fs.panel_label(ax, "b")

    # ---------------- panel c: headroom at the observed rate ----------------
    ax = axs[2]
    head = pd.DataFrame({
        "pool": ["misfolded\nmonomer $P$", "aggregated\npool $A$"],
        "headroom": [cross["P_headroom_ratio"], cross["A_headroom_ratio"]],
    })
    head["log10"] = np.log10(head.headroom)
    sns.barplot(data=head, x="pool", y="log10", ax=ax, hue="pool", legend=False,
                palette=[fs.C["burden"], fs.C["capacity"]], alpha=0.9,
                edgecolor="white", linewidth=0.8)
    for i, row in head.iterrows():
        ax.text(i, row["log10"] + 0.08, f"$\\times${row.headroom:,.0f}",
                ha="center", va="bottom", fontsize=7)
    ax.axhline(0, color=fs.C["alert"], linewidth=1.0)
    ax.text(0.02, 0.02, "collapse threshold", transform=ax.transAxes,
            fontsize=6.2, color=fs.C["alert"], va="bottom")
    ax.set_ylim(0, head["log10"].max() * 1.22)
    ax.set_ylabel(r"$\log_{10}$ headroom to threshold")
    ax.set_xlabel("")
    ax.set_title(r"At $f = 10^{-4}$ the system rests" "\n" "deep on the stable branch")
    fs.panel_label(ax, "c")

    sns.despine(fig=fig)
    fig.tight_layout(w_pad=1.9)
    fs.save(fig, str(FIGS / "Fig3_bounds"))

    # ---- emit the numbers the manuscript and tests consume ----
    summary = {
        "source": "paired_mc_results.json (valid paired MC), "
                  "arithmetic_results.json, two_pool_results.json",
        "observed_window_per_codon": [OBS_LO, OBS_HI],
        "arithmetic_deterministic": det_arith,
        "arithmetic_paired_median": float(np.median(f_arith)),
        "arithmetic_paired_ci95": [float(np.percentile(f_arith, 2.5)),
                                   float(np.percentile(f_arith, 97.5))],
        "ode_deterministic": det_ode,
        "ode_paired_median": float(np.median(f_ode)),
        "ode_paired_ci95": [float(np.percentile(f_ode, 2.5)),
                            float(np.percentile(f_ode, 97.5))],
        "paired_median_ratio_arith_over_ode": float(np.median(r)),
        "paired_P_arith_tighter": p_arith_tighter,
        "ode_over_arith_at_median": float(np.median(f_ode) / np.median(f_arith)),
        "headroom_P": cross["P_headroom_ratio"],
        "headroom_A": cross["A_headroom_ratio"],
        "mechanism_frac_aggregation_death":
            twopool["D_monte_carlo"]["mechanism"]["frac_aggregation_death"],
        "retracted_statistic_not_used": {
            "value": 0.654,
            "why": "PAIRED_MC_TASK.md: shared parameters were drawn from "
                   "different distributions in the two marginal runs, so it is "
                   "not a valid paired statistic",
            "replacement": p_arith_tighter,
        },
    }
    (COMP / "bounds_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"  wrote {COMP/'bounds_summary.json'}")
    for k in ("arithmetic_deterministic", "arithmetic_paired_median",
              "ode_deterministic", "ode_paired_median",
              "paired_median_ratio_arith_over_ode", "paired_P_arith_tighter",
              "ode_over_arith_at_median"):
        print(f"    {k:38s} {summary[k]:.4g}")


if __name__ == "__main__":
    main()
