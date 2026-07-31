"""
shared figure style. seaborn is the plotting interface throughout; matplotlib is
used only for figure handles and saving, per project guidelines.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

# jitter seed. seaborn's stripplot draws jitter from the numpy global RNG, so
# without this Fig 2c's points move on every rebuild -- the same data producing a
# visibly different panel, which a rebuild-to-verify workflow cannot check.
JITTER_SEED = 42

COL_SINGLE = 3.42   # inches, single-column
COL_DOUBLE = 7.09   # inches, double-column

C = {
    "burden":   "#1f4e79",   # deep blue  - burden / arithmetic
    "capacity": "#2a9d8f",   # teal       - capacity / headroom
    "alert":    "#c1121f",   # red        - observed / threshold
    "ode":      "#e07a5f",   # warm       - ODE
    "neutral":  "#b8b8b8",   # grey       - nulls / background
    "muted":    "#6b6b6b",
}


def setup():
    np.random.seed(JITTER_SEED)
    sns.set_theme(context="paper", style="ticks", font_scale=0.85)
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "font.family": "sans-serif",
        "axes.linewidth": 0.7,
        "xtick.major.width": 0.7,
        "ytick.major.width": 0.7,
        "axes.titlesize": 8,
        "axes.titleweight": "normal",
        "legend.frameon": False,
        # stable element ids in the svg. without a fixed salt matplotlib derives
        # them from object addresses, so every rebuild rewrites every <path id>
        "svg.hashsalt": "envelope-paper",
    })


def panel_label(ax, letter, dx=-0.16, dy=1.06):
    ax.text(dx, dy, letter, transform=ax.transAxes,
            fontsize=10, fontweight="bold", va="top", ha="left")


def save(fig, path_noext):
    # svg as well as png/pdf: 15_build_paper.py substitutes the svg when it builds
    # the html, so headless Chrome prints the figures as vector rather than
    # re-rastering a 300 dpi png. the markdown keeps the png, which every reader
    # and forge renders
    # drop the build timestamp: with it, every rebuild rewrites every figure and
    # `git status` after a verification rebuild is all noise, which hides a real
    # change among a dozen fake ones
    no_date = {"pdf": {"CreationDate": None}, "svg": {"Date": None}, "png": {}}
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{path_noext}.{ext}", metadata=no_date[ext])
    plt.close(fig)
    print(f"  wrote {path_noext}.png / .pdf / .svg")
