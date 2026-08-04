"""Figure 5 -- under division, the enzymatic condition does not move; dilution does.

Corollary 3 splits `j_crit = C_enz . phi_enz / (1 - delta)`. The claim the figure
has to carry is a CONTRAST between two curves, so they share one linear y-axis:
`phi_enz` sits flat near 0.13 while `delta` climbs to 0.39 over the same sweep.
Twin axes would rescale the flat curve until its numerical noise looked like
structure, which is the opposite of the point.

Everything is a deterministic continuation at the base parameter point under
CONSTANT dilution (`k_mu = inf`). There is no sampling, so the reported band on
`phi_enz` is a complete enumeration over the sweep rather than an extremum over a
subsample -- the distinction that D036-D038 turned on. Computed here, not
extracted: no run root is involved.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

# manuscript figure number, by ORDER OF FIRST MENTION in bmb_v4.md.
# filenames are deliberately semantic so a reorder touches this line only.
FIGURE = "fig2"

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "figures",
           REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import _figstyle as F  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from proteostasis import model as M  # noqa: E402
import dilution as D  # noqa: E402
import boundary_structure as BS  # noqa: E402

MU_HI = 0.08          # the range section 4.2 quotes; the fold is gone by mu = 0.10
N_STEP = 33


def sweep():
    """continuation in mu at the base parameters, constant dilution."""
    p = M.Params()
    rows = []
    for mu in np.linspace(0.0, MU_HI, N_STEP):
        dec = BS.boundaryDecomposition(p, D.Growth(mu0=float(mu)))
        if dec is None:
            continue
        rows.append({"mu": float(mu), "phi_enz": dec["phi_enz"],
                     "delta": dec["delta"], "j_crit": dec["j_crit"],
                     "a_star": dec["a_star"],
                     "identity_rel_err": dec["identity_rel_err"]})
    return rows


def build():
    F.setStyle()
    rows = sweep()
    if len(rows) < N_STEP:
        raise RuntimeError(f"continuation lost the fold: {len(rows)}/{N_STEP}")
    mu = np.array([r["mu"] for r in rows])
    phi = np.array([r["phi_enz"] for r in rows])
    dlt = np.array([r["delta"] for r in rows])

    fig, ax = plt.subplots(figsize=(F.W_DOUBLE, 0.66 * F.W_DOUBLE))

    ax.plot(mu, dlt, "-", lw=1.3, color="#b3341f", zorder=3)
    ax.plot(mu, phi, "-", lw=1.3, color="#1b3a6b", zorder=3)

    # the band on phi_enz, drawn so "flat" is a measured statement not an eyeball one
    ax.fill_between(mu, phi.min(), phi.max(), color="#1b3a6b", alpha=0.13, lw=0,
                    zorder=1)

    # labels sit clear of both curves: delta above its own line on the left half,
    # phi_enz above its band on the right, where delta has already climbed past.
    ax.text(0.004, 0.425,
            r"$\delta$  (dilution share):  0 → " + f"{dlt.max():.2f}",
            fontsize=6.0, color="#b3341f", ha="left", va="top")
    ax.text(mu[-1] * 0.99, phi.max() + 0.030,
            r"$\varphi_{\mathrm{enz}}$  (enzymatic share):  "
            + f"{phi.min():.4f}–{phi.max():.4f}",
            fontsize=6.0, color="#1b3a6b", ha="right", va="bottom")

    ax.set_xlim(0.0, MU_HI)
    ax.set_ylim(0.0, 0.45)
    ax.set_xlabel(r"dilution rate $\mu$  (constant, $k_\mu = \infty$)")
    ax.set_ylabel("dimensionless share of the critical influx")
    ax.set_title(f"{len(rows)}-point continuation at the base parameters, "
                 "complete", loc="left", fontsize=5.8, color="0.35")
    ax.grid(True, axis="y", color="0.90", lw=0.4, zorder=0)

    fig.tight_layout(pad=0.35)
    F.widthCheck(fig, F.W_DOUBLE)
    hashes = F.save(fig, FIGURE)
    plt.close(fig)
    return {"n": len(rows), "mu_hi": float(mu.max()),
            "phi_lo": float(phi.min()), "phi_hi": float(phi.max()),
            "phi_spread_pct": float(100.0 * (phi.max() - phi.min())
                                    / ((phi.max() + phi.min()) / 2.0)),
            "delta_lo": float(dlt.min()), "delta_hi": float(dlt.max()),
            "identity_max": float(max(r["identity_rel_err"] for r in rows)),
            "identity_p99": float(np.percentile(
                [r["identity_rel_err"] for r in rows], 99)),
            "j_crit_lo": float(rows[0]["j_crit"]),
            "j_crit_hi": float(rows[-1]["j_crit"]),
            "a_star_lo": float(rows[0]["a_star"]),
            "a_star_hi": float(rows[-1]["a_star"]),
            "hashes": hashes}


if __name__ == "__main__":
    o = build()
    print(f"Figure {FIGURE[3:]}  (%d-point continuation, complete enumeration, no sampling)"
          % o["n"])
    print("  mu 0 -> %.3f" % o["mu_hi"])
    print("  phi_enz : %.4f - %.4f   (full width %.1f%% of the mean)"
          % (o["phi_lo"], o["phi_hi"], o["phi_spread_pct"]))
    print("  delta   : %.4f - %.4f" % (o["delta_lo"], o["delta_hi"]))
    print("  j_crit  : %.4f -> %.4f     a*: %.3f -> %.3f"
          % (o["j_crit_lo"], o["j_crit_hi"], o["a_star_lo"], o["a_star_hi"]))
    print("  Corollary 3 identity residual: p99 %.2e, max %.2e over the sweep"
          % (o["identity_p99"], o["identity_max"]))
    for k, v in o["hashes"].items():
        print("  %-10s %s" % (k, v[:16]))
