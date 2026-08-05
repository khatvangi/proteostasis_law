"""Figure 4 -- the fold is not always the first loss of stability.

Section 7's three claims, one panel each, chosen so that each panel carries a
DIFFERENT kind of evidence rather than three views of the same one:

  (a) the phenomenon.   `tr J` along the low-burden branch against `j/j_crit`,
      for a crossing network and a non-crossing one. This is the branch
      computation, and by itself it is exactly what a tracing artefact would
      also produce -- which is why it is not the whole figure.
  (b) the control.      Amplification of a perturbation at each network's own
      `tr J` maximum, crossers against non-crossers. The nonlinear integration
      shares no failure mode with (a). The panel's content is the SEPARATION,
      so both distributions are drawn even though one of them does nothing.
  (c) the region.       Where the crossers sit in `(kappa_a, rho_A)`. An artefact
      has no reason to concentrate on two kinetic parameters, and this is the
      transferable half of the result: the incidence rate belongs to the
      sampling box, the location belongs to the model.

Every number in the caption comes from `captionNumbers`, which reads the same
computed tables the panels do, so the caption cannot drift from the figure.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# manuscript figure number, by ORDER OF FIRST MENTION in the manuscript.
# filenames are deliberately semantic so a reorder touches this line only.
FIGURE = "fig4"

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "figures",
           REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import _figstyle as F  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import fold_theorem as FT  # noqa: E402
import genericity as GEN  # noqa: E402
import hopf_check as HC  # noqa: E402

COMPUTED = REPO_ROOT / "data" / "computed"

# the two networks drawn in panel (a). Chosen by RULE, not by eye.
#
# The rule carries one extra clause that the first version lacked, and the
# figure is why: picking purely on "crossing nearest the median position" chose
# a network whose traced branch begins at 0.83 j_crit, so the panel showed a
# 2 mm red stub with no visible approach and no visible crossing. A panel whose
# job is to show that the crossing PRECEDES the fold has to show the run-up.
# So the crossing exemplar is drawn from networks whose branch is resolved from
# below 0.3 j_crit, and among those it is the one nearest the median crossing.
_EXEMPLAR_RULE = ("median first-crossing position among branches resolved "
                  "below 0.3 j_crit; median max tr J")
_RESOLVED_BELOW = 0.30


def _tables():
    S = pd.read_csv(COMPUTED / "hopf_refined_kinetic_box.tsv", sep="\t")
    I = pd.read_csv(COMPUTED / "hopf_integration_kinetic_box.tsv", sep="\t")
    T = S[S["traced"] == True]  # noqa: E712
    clean = T[(T["tr_max"] >= 0.0) & (T["fold_is_j_max"] == 1)]
    return T, clean, I[I["tested"] == True]  # noqa: E712


def pickExemplars(T, clean):
    frac = clean["j_at_first_cross"] / clean["j_crit"]
    wide = clean[clean["j_lo"] / clean["j_crit"] < _RESOLVED_BELOW]
    pool, target = (wide, frac.median()) if len(wide) else (clean, frac.median())
    pool_frac = pool["j_at_first_cross"] / pool["j_crit"]
    cross_name = pool.loc[(pool_frac - target).abs().idxmin(), "name"]
    stable = T[T["tr_max"] < 0.0]
    stable_name = stable.loc[
        (stable["tr_max"] - stable["tr_max"].median()).abs().idxmin(), "name"]
    return str(cross_name), str(stable_name)


def branchOf(name: str):
    """recompute the branch profile for one network, from the run root."""
    run = FT.phase1RunDir()
    for nm, p, u, a in GEN.kineticBoxStates(run):
        if nm == name:
            out = HC.branchProfile(p, u, a)
            if out is None:
                return None
            B = out["branch"]
            j_crit = float(FT.removalR(u, a, p))
            return B.sort_values("j"), j_crit
    return None


def captionNumbers() -> dict:
    """every quantity the caption quotes, from the same tables the panels read."""
    T, clean, I = _tables()
    L = pd.read_csv(COMPUTED / "hopf_refined_load_grid.tsv", sep="\t")
    L = L[L["traced"] == True]  # noqa: E712
    cr = I[I["name"].str.startswith("+")]
    ct = I[I["name"].str.startswith("-")]
    g = pd.read_csv(COMPUTED / "hopf_parameter_corner.tsv", sep="\t", index_col=0)
    frac = clean["j_at_first_cross"] / clean["j_crit"]
    per = cr[cr["period_measured"].notna()]
    pr = (per["period_measured"] - per["period_predicted"]).abs() \
        / per["period_predicted"]
    fit = cr[cr["slope"].notna()]
    rel = (fit["slope"] - fit["lambda_max"]).abs() / fit["lambda_max"].abs()
    return {
        "n_traced": int(len(T)),
        "n_cross": int(len(clean)),
        "pct_cross": 100.0 * len(clean) / len(T),
        "n_load_grid": int(len(L)),
        "load_grid_tr_max": float(L["tr_max"].max()),
        "frac_median": float(frac.median()),
        "frac_min": float(frac.min()),
        "det_min_at_cross": float(clean["det_at_first_cross"].min()),
        "n_ctrl": int(len(ct)),
        "ctrl_escaped": int(ct["escaped"].sum()),
        "ctrl_ratio_median": float(ct["ratio_max"].median()),
        "cross_grew": int(cr["grew"].sum()),
        "cross_escaped": int(cr["escaped"].sum()),
        "cross_ratio_median": float(cr["ratio_max"].median()),
        "n_period": int(len(per)),
        "period_within5": int((pr < 0.05).sum()),
        "n_fit": int(len(fit)),
        "exp_within5": int((rel < 0.05).sum()),
        "rho_A_ratio": float(g.loc["p_rho_A", "ratio"]),
        "kappa_a_ratio": 1.0 / float(g.loc["p_kappa_a", "ratio"]),
    }


def build():
    F.setStyle()
    T, clean, I = _tables()
    cross_name, stable_name = pickExemplars(T, clean)

    fig, axes = plt.subplots(1, 3, figsize=(F.W_SINGLE, 52.0 * F.MM))
    ax_a, ax_b, ax_c = axes

    # ---- (a) tr J along the branch ---------------------------------------
    for name, colour, label in ((stable_name, "0.45", "no crossing"),
                                (cross_name, "#c1272d", "crossing")):
        got = branchOf(name)
        if got is None:
            continue
        B, j_crit = got
        ax_a.plot(B["j"] / j_crit, B["tr_J"], color=colour, lw=1.1, label=label)
    ax_a.axhline(0.0, color="0.2", lw=0.5, ls=(0, (3, 2)))
    ax_a.axvline(1.0, color="0.2", lw=0.5, ls=(0, (1, 2)))
    ax_a.set_xlabel(r"influx $j/j_{\mathrm{crit}}$")
    ax_a.set_ylabel(r"$\mathrm{tr}\,J$ on the low-burden branch")
    # symlog: the two exemplars differ by an order of magnitude in tr J, and on a
    # linear axis the crossing network is a flat line against the other's -6.5.
    # the linear window straddles zero so the crossing itself stays readable.
    ax_a.set_yscale("symlog", linthresh=0.1, linscale=0.5)
    # explicit ticks: symlog's own labelling put 10^-2 and -10^-2 on top of each
    # other across the zero line, which is where the panel's whole claim sits
    ax_a.set_yticks([-10.0, -1.0, 0.0, 1.0])
    ax_a.set_yticklabels(["-10", "-1", "0", "1"])
    ax_a.set_xlim(0.0, 1.04)
    ax_a.legend(loc="upper left")
    ax_a.text(0.02, 1.06, "a", transform=ax_a.transAxes, fontweight="bold")

    # ---- (b) the control ---------------------------------------------------
    cr = I[I["name"].str.startswith("+")]["ratio_max"]
    ct = I[I["name"].str.startswith("-")]["ratio_max"]
    bins = np.logspace(np.log10(min(ct.min(), cr.min(), 0.5)),
                       np.log10(max(cr.max(), ct.max()) * 1.4), 30)
    ax_b.hist(ct, bins=bins, color="0.45", alpha=0.85,
              label=f"no crossing (n={len(ct)})")
    ax_b.hist(cr, bins=bins, color="#c1272d", alpha=0.85,
              label=f"crossing (n={len(cr)})")
    ax_b.set_xscale("log")
    ax_b.set_xlabel(r"perturbation amplification $\max|\delta|/\delta_0$")
    ax_b.set_ylabel("networks")
    ax_b.legend(loc="upper center")
    ax_b.text(0.02, 1.06, "b", transform=ax_b.transAxes, fontweight="bold")

    # ---- (c) the parameter region -----------------------------------------
    run = FT.phase1RunDir()
    c = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    c = c[pd.to_numeric(c["fold_burden"], errors="coerce").notna()]
    c = c.assign(name=[f"draw{i}" for i in c.index])
    m = c.merge(T[["name", "tr_max", "fold_is_j_max"]], on="name")
    m["cross"] = (m["tr_max"] >= 0.0) & (m["fold_is_j_max"] == 1)
    ax_c.scatter(m.loc[~m["cross"], "p_kappa_a"], m.loc[~m["cross"], "p_rho_A"],
                 s=1.4, c="0.72", linewidths=0, rasterized=True)
    ax_c.scatter(m.loc[m["cross"], "p_kappa_a"], m.loc[m["cross"], "p_rho_A"],
                 s=4.0, c="#c1272d", linewidths=0)
    ax_c.set_xscale("log")
    ax_c.set_yscale("log")
    ax_c.set_xlabel(r"$\kappa_a$  (clearance saturation)")
    ax_c.set_ylabel(r"$\rho_A$  (clearance rate)")
    ax_c.text(0.02, 1.06, "c", transform=ax_c.transAxes, fontweight="bold")

    fig.tight_layout(pad=0.4, w_pad=1.4)
    F.widthCheck(fig, F.W_SINGLE)
    return fig, captionNumbers(), (cross_name, stable_name)


def main() -> int:
    fig, n, ex = build()
    digests = F.save(fig, FIGURE)
    plt.close(fig)
    print(f"  exemplars (by rule: {_EXEMPLAR_RULE}): "
          f"crossing {ex[0]}, stable {ex[1]}")
    for k, v in n.items():
        print(f"  {k:24s} {v}")
    for name, d in digests.items():
        print(f"  {name:14s} {d[:16]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
