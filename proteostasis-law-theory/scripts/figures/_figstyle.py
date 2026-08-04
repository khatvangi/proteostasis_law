"""shared figure style and output for the BMB submission.

CONSTRAINTS THIS ENCODES (Bulletin of Mathematical Biology)
  * 84 mm wide for a double-column figure, 174 mm for single-column
  * no taller than 234 mm
  * vector output (PDF and SVG) plus a PNG for inspection

DELIBERATE CHOICES
  * NO seaborn, and no reliance on any local matplotlib configuration. This
    contradicts the project's usual plotting default and is a considered override
    for submission: a figure that renders differently on the typesetter's machine
    than on ours is a defect, so every rcParam that matters is set explicitly
    after `rcdefaults()` rather than inherited.
  * timestamps are stripped from SVG and PDF output, so regenerating a figure
    with unchanged inputs produces a byte-identical file and diffs stay clean.
  * one fixed seed, exported, for anything involving a draw. an unseeded jitter
    in an earlier project cost a rerun.
"""

from __future__ import annotations

import hashlib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
FIG_DIR = REPO_ROOT / "figures"
DATA_DIR = REPO_ROOT / "data" / "figures"

SEED = 20260804

MM = 1.0 / 25.4
W_DOUBLE = 84.0 * MM        # double-column width
W_SINGLE = 174.0 * MM       # single-column (full text width)
H_MAX = 234.0 * MM

# a fixed salt makes SVG element ids deterministic across runs
_HASHSALT = "proteostasis-collapse-threshold"


def setStyle() -> None:
    """explicit rcParams, so output does not depend on a local matplotlibrc."""
    plt.rcdefaults()
    matplotlib.rcParams.update({
        "svg.hashsalt": _HASHSALT,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans"],
        "font.size": 7.0,
        "axes.labelsize": 7.0,
        "axes.titlesize": 7.5,
        "xtick.labelsize": 6.5,
        "ytick.labelsize": 6.5,
        "legend.fontsize": 6.0,
        "legend.frameon": False,
        "axes.linewidth": 0.6,
        "grid.linewidth": 0.4,
        "lines.linewidth": 1.0,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,        # embed as TrueType, editable text
        "ps.fonttype": 42,
        "text.usetex": False,
    })


def save(fig, stem: str) -> dict:
    """write PDF + SVG + PNG with timestamps stripped; return sha256 per file."""
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    w_in, h_in = fig.get_size_inches()
    if h_in > H_MAX + 1e-9:
        raise ValueError(f"{stem}: height {h_in/MM:.1f} mm exceeds the 234 mm limit")
    out = {}
    fig.savefig(FIG_DIR / f"{stem}.pdf", format="pdf",
                metadata={"CreationDate": None, "Producer": None, "Creator": None})
    fig.savefig(FIG_DIR / f"{stem}.svg", format="svg", metadata={"Date": None})
    fig.savefig(FIG_DIR / f"{stem}.png", format="png", metadata={"Software": None})
    for ext in ("pdf", "svg", "png"):
        p = FIG_DIR / f"{stem}.{ext}"
        out[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
    return out


def widthCheck(fig, target: float) -> None:
    w_in, _ = fig.get_size_inches()
    if abs(w_in - target) > 1e-6:
        raise ValueError(f"figure width {w_in/MM:.2f} mm is not the required "
                         f"{target/MM:.0f} mm")
