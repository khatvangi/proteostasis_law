"""phase 3: does critical slowing survive cell division, and with what exponent?

`empirical/GATE4_PROPOSAL.md` §10.4 blocks the gate on this. Critical slowing was
established in Phase 1, which had no cell division; dilution changes both the
boundary location and the Jacobian (`J - mu.I` for constant `mu`), so the
eigenvalue governing recovery time is directly affected.

THE SHARPER FORM
----------------
"recovery time increases" is weak. At a saddle-node one eigenvalue passes through
zero, and the generic normal form gives

        |lambda|  ~  (j_crit - j)^(1/2)      so   tau = 1/|lambda| ~ (j_crit - j)^(-1/2)

The exponent 1/2 is **parameter-free**: it follows from the bifurcation being a
saddle-node, not from any rate constant. That makes it far more robust to the
parameter uncertainty that killed the saturation-fraction design (§9 of the
proposal), and it is measurable as a slope on a log-log plot rather than as an
absolute magnitude.

What this does and does not test: an exponent of 1/2 is evidence for a
saddle-node, which is what the fold theorem asserts the boundary IS. It is not
by itself evidence for the specific two-state model. The alternative it
discriminates against is a smooth decline with no bifurcation, where recovery
time stays bounded.

CLAIM LABELS
  Mathematical  : the 1/2 exponent, from the saddle-node normal form.
  Computational : every fitted slope below.
  Empirical     : nothing. no organism data is used.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
from scipy.optimize import root

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import dilution as D  # noqa: E402
import calibration as C  # noqa: E402
from fold_theorem import KINETIC_FIELDS, phase1RunDir  # noqa: E402


def stableEquilibrium(j: float, p: M.Params, g: D.Growth,
                      guess=(1e-3, 1e-6)) -> Optional[tuple]:
    """the low-burden equilibrium at influx j, found by damped Newton in logs."""
    pj = p.with_(j=float(j)).validate()

    def res(x):
        # bound the search: the nucleation term is u^m, which overflows long
        # before the solver would ever accept such a state
        lu = float(np.clip(x[0], -34.0, 6.0))
        la = float(np.clip(x[1], -34.0, 6.0))
        u, a = float(np.exp(lu)), float(np.exp(la))
        try:
            du, da = D.rhsDil(u, a, pj, g)
        except (M.ModelError, OverflowError, FloatingPointError):
            return [1e6, 1e6]
        if not (np.isfinite(du) and np.isfinite(da)):
            return [1e6, 1e6]
        return [du, da]

    s = root(res, [np.log(max(guess[0], 1e-12)), np.log(max(guess[1], 1e-14))],
             method="hybr", options={"xtol": 1e-13})
    if not s.success:
        return None
    u = float(np.exp(float(np.clip(s.x[0], -34.0, 6.0))))
    a = float(np.exp(float(np.clip(s.x[1], -34.0, 6.0))))
    if not (np.isfinite(u) and np.isfinite(a)):
        return None
    try:
        # reject anything that is not actually a root
        du, da = D.rhsDil(u, a, pj, g)
        if max(abs(du), abs(da)) > 1e-9 * max(1.0, u + a):
            return None
        ev = np.linalg.eigvals(D.jacobianDil(u, a, pj, g))
    except (M.ModelError, OverflowError, np.linalg.LinAlgError):
        return None
    lead = float(np.max(ev.real))
    if lead >= 0.0:                       # not the stable branch
        return None
    return u, a, lead


def slowingExponent(p: M.Params, g: D.Growth, n: int = 14,
                    lo: float = 1e-4, hi: float = 3e-2) -> Optional[Dict[str, float]]:
    """fit log|lambda| against log(j_crit - j) on the approach to the boundary.

    the saddle-node normal form predicts slope 1/2. distances are taken as
    RELATIVE offsets from j_crit so the ladder is comparable across networks.
    """
    try:
        out = D.foldSolveDil(p, g)
    except (M.ModelError, OverflowError, np.linalg.LinAlgError):
        return None
    if out is None:
        return None
    j_crit = out[0]

    # continue OUTWARD from the fold state. just below j_crit the stable
    # equilibrium sits near (u*, a*), so the fold is the only reliable seed;
    # starting from a fixed far guess fails to converge at every rung.
    _, u_f, a_f = out
    ds, lams = [], []
    guess = (u_f, max(a_f, 1e-12))
    for rel in np.geomspace(lo, hi, n):          # near -> far
        j = j_crit * (1.0 - rel)
        eq = stableEquilibrium(j, p, g, guess)
        if eq is None:
            continue
        u, a, lead = eq
        guess = (u, max(a, 1e-14))
        if lead >= 0.0 or not np.isfinite(lead):
            continue
        ds.append(j_crit - j)
        lams.append(abs(lead))
    if len(ds) < 6:
        return None
    x, y = np.log10(np.array(ds)), np.log10(np.array(lams))
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    ss = 1.0 - float(np.sum(resid ** 2) / np.sum((y - y.mean()) ** 2))
    return {"j_crit": j_crit, "n_points": len(ds), "slope": float(slope),
            "r2": ss, "tau_exponent": -float(slope)}


def sweep(k: int = 14, seed: int = 23) -> pd.DataFrame:
    """the exponent across networks, without dilution and with the measured law."""
    run = phase1RunDir()
    c = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    draws = c.sample(n=k, random_state=seed)

    rows: List[Dict[str, float]] = []
    for idx, (_, r) in enumerate(draws.iterrows()):
        kin = {f: float(r["p_" + f]) for f in KINETIC_FIELDS}
        chi = float(r["p_chi"])
        try:
            p = M.Params(**kin).with_(nu=float(r["p_nu"]), c_tot=chi,
                                      d_tot=1.0 - chi).validate()
        except M.ModelError:
            continue
        for label, g in (("no dilution", D.Growth(0.0)),
                         ("calibrated mu0=0.05", C.calibratedGrowth(0.05, 0.02)),
                         ("calibrated mu0=0.10", C.calibratedGrowth(0.10, 0.02))):
            try:
                s = slowingExponent(p, g)
            except (M.ModelError, ValueError, OverflowError,
                    np.linalg.LinAlgError):
                continue
            if s is None:
                continue
            s.update({"draw": idx, "regime": label})
            rows.append(s)
    return pd.DataFrame(rows)


def main() -> int:
    if not phase1RunDir().is_dir():
        print("SKIP: no phase 1 run root (results/ is gitignored).")
        return 0

    print("CRITICAL SLOWING: does tau ~ (j_crit - j)^(-1/2) survive dilution?")
    print("   the saddle-node normal form predicts eigenvalue slope 0.5 exactly.")
    print()
    df = sweep()
    if df.empty:
        print("   no networks yielded a usable ladder")
        return 0
    for regime, sub in df.groupby("regime"):
        print("%-22s networks=%-4d  slope median %.4f  (IQR %.4f-%.4f)  r2 median %.4f"
              % (regime, len(sub), sub["slope"].median(),
                 sub["slope"].quantile(.25), sub["slope"].quantile(.75),
                 sub["r2"].median()))
    print()
    near = df[(df["slope"] - 0.5).abs() < 0.05]
    print("   networks within 0.05 of the predicted 0.5 : %d of %d (%.1f%%)"
          % (len(near), len(df), 100 * len(near) / len(df)))
    print("   overall slope median %.4f, r2 median %.4f"
          % (df["slope"].median(), df["r2"].median()))
    print()
    print("   -> %s" % ("SURVIVES: the exponent is unchanged by dilution, so the "
                        "primary outcome stands"
                        if abs(df["slope"].median() - 0.5) < 0.05 else
                        "DOES NOT hold at 0.5; see the fitted values above"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
