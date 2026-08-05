"""Are the 26 folds with d2R/ds2 > 0 collapse points, or the other turning point?

Block 1a found `d2R/ds2 < 0` at 325 of 325 load-grid folds -- the fold is a
constrained local MAXIMUM of removal along the aggregate nullcline -- but at 26
of 2765 kinetic-box folds the second derivative is positive.

WHY THAT IS NOT A ROUNDING DETAIL
`j_crit = R(u*, a*)`, so the sign of `d2R/ds2` is the sign of `d2j/ds2` along the
curve. At a local MAXIMUM of `R`, solutions of `R = j` exist just below `j_crit`
and vanish just above: the equilibrium pair is DESTROYED as the influx rises,
which is collapse. At a local MINIMUM they are BORN as the influx rises. A fold
of the second kind is the lower turning point of a hysteresis loop, and calling
it a collapse threshold would be wrong.

WHAT THIS MODULE DOES
For each network it enumerates every `det J = 0` point on the traced nullcline
rather than the one `foldSolve` happened to return, polishes each, and classifies
it by the sign of `d2R/ds2`. If the 26 carry a second, oppositely oriented
singular point then folds come in two orientations and the model class exhibits
both -- a classification rather than an exception list.

The orientation is computed TWICE by different routes: the arclength second
difference of `R` from `genericity.py`, and the shape of `j` along the traced
profile near the point. They are independent up to the projection they share.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import fold_theorem as FT  # noqa: E402
import genericity as GEN  # noqa: E402
import hopf_check as HC  # noqa: E402

OUT_DIR = REPO_ROOT / "data" / "computed"


def singularCandidates(p: M.Params, u_star: float, a_star: float,
                       n: int = 150) -> Optional[List[Dict]]:
    """every det J = 0 point on the accessible nullcline, polished and classified.

    Seeds from each sign change of `det J` along the projected profile, so the
    candidate set comes from the curve rather than from whichever root a solver
    reached first -- which is the failure `polishFold` was written to avoid.
    """
    out = HC.branchProfile(p, u_star, a_star, n=n)
    if out is None:
        return None
    D = out["profile"]
    det = D["det_J"].to_numpy()
    js = D["j"].to_numpy()
    us = D["u"].to_numpy()
    as_ = D["a"].to_numpy()

    cands = []
    for k in np.nonzero(det[:-1] * det[1:] < 0.0)[0]:
        # seed at the sign change, polish onto {G = 0, det J = 0}
        s = GEN.polishFold(p, 0.5 * (us[k] + us[k + 1]),
                           0.5 * (as_[k] + as_[k + 1]))
        if s is None:
            cands.append({"polished": False, "seed_u": float(us[k]),
                          "seed_a": float(as_[k])})
            continue
        j_c, u_c, a_c = s
        g = GEN.genericityAt(u_c, a_c, p)
        # second route to the same orientation: the shape of j along the traced
        # profile in a window around the crossing. no arclength, no reprojection.
        w = slice(max(0, k - 6), min(len(js), k + 8))
        loc = js[w]
        by_profile = np.nan
        if loc.size >= 5:
            mid = np.argmin(np.abs(np.arange(w.start, w.stop) - k))
            by_profile = float(loc[mid] - 0.5 * (loc[0] + loc[-1]))
        cands.append({
            "polished": True, "u": u_c, "a": a_c, "j_crit": float(j_c),
            "d2R_signed": np.nan if g is None else g["d2R_signed"],
            "tr_J": np.nan if g is None else g["tr_J"],
            "det_J": np.nan if g is None else g["det_J"],
            # > 0 means j sits above its neighbours: a maximum, hence collapse
            "j_shape": by_profile,
        })

    # dedupe: two sign changes can polish to the same point
    keep = []
    for c in cands:
        if not c.get("polished"):
            keep.append(c)
            continue
        if any(k.get("polished")
               and abs(np.log(max(c["u"], 1e-300)) - np.log(max(k["u"], 1e-300))) < 1e-6
               and abs(np.log(max(c["a"], 1e-300)) - np.log(max(k["a"], 1e-300))) < 1e-6
               for k in keep):
            continue
        keep.append(c)
    return keep


def _worker(item):
    name, p, u, a = item
    try:
        cands = singularCandidates(p, u, a)
    except Exception:
        cands = None
    if cands is None:
        return {"name": name, "ok": False}
    pol = [c for c in cands if c.get("polished")]
    d2 = [c["d2R_signed"] for c in pol if np.isfinite(c.get("d2R_signed", np.nan))]
    shp = [c["j_shape"] for c in pol if np.isfinite(c.get("j_shape", np.nan))]
    return {
        "name": name, "ok": True,
        "n_candidates": len(cands),
        "n_polished": len(pol),
        "n_max": int(sum(1 for x in d2 if x < 0)),      # collapse orientation
        "n_min": int(sum(1 for x in d2 if x > 0)),      # birth orientation
        "has_both": bool(any(x < 0 for x in d2) and any(x > 0 for x in d2)),
        "j_at_max": float(max([c["j_crit"] for c in pol
                               if np.isfinite(c.get("d2R_signed", np.nan))
                               and c["d2R_signed"] < 0], default=np.nan)),
        "j_at_min": float(min([c["j_crit"] for c in pol
                               if np.isfinite(c.get("d2R_signed", np.nan))
                               and c["d2R_signed"] > 0], default=np.nan)),
        # do the two routes agree on orientation at every polished candidate?
        "orient_agree": int(sum(1 for c in pol
                                if np.isfinite(c.get("d2R_signed", np.nan))
                                and np.isfinite(c.get("j_shape", np.nan))
                                and (c["d2R_signed"] < 0) == (c["j_shape"] > 0))),
        "orient_n": int(min(len(d2), len(shp))),
    }


def main() -> int:
    run = FT.phase1RunDir()
    states = GEN.kineticBoxStates(run)
    byname = {nm: (p, u, a) for nm, p, u, a in states}

    G = pd.read_csv(OUT_DIR / "genericity_kinetic_box.tsv", sep="\t")
    G = G[G["ok"] == True]  # noqa: E712
    pos = G[G["d2R_signed"] > 0]["name"].tolist()
    neg = G[G["d2R_signed"] < 0]["name"].tolist()
    print(f"kinetic box: d2R/ds2 > 0 at {len(pos)} folds, < 0 at {len(neg)}")

    # the 26, plus a control block of the ordinary orientation
    ctrl = neg[::max(1, len(neg) // 150)]
    items = [(nm, *byname[nm]) for nm in pos + ctrl if nm in byname]
    labels = {nm: ("d2R>0" if nm in set(pos) else "d2R<0") for nm, *_ in items}

    import multiprocessing as mp
    workers = max(1, min(56, (os.cpu_count() or 2) - 8))
    ctx = mp.get_context("fork")
    with ctx.Pool(processes=workers) as pool:
        rows = pool.map(_worker, items, chunksize=2)
    R = pd.DataFrame(rows)
    R["group"] = R["name"].map(labels)
    R.to_csv(OUT_DIR / "fold_orientation.tsv", sep="\t", index=False)

    for grp, sub in R.groupby("group"):
        ok = sub[sub["ok"] == True]  # noqa: E712
        print(f"\n  {grp}: {len(sub)} networks, {len(ok)} traced")
        if ok.empty:
            continue
        print(f"    singular candidates per network: median "
              f"{ok['n_candidates'].median():.0f}  max {ok['n_candidates'].max():.0f}"
              f"  polished {ok['n_polished'].sum()}/{ok['n_candidates'].sum()}")
        print(f"    candidates of COLLAPSE orientation (d2R<0): median "
              f"{ok['n_max'].median():.0f}  none in {int((ok['n_max'] == 0).sum())}")
        print(f"    candidates of BIRTH orientation (d2R>0):    median "
              f"{ok['n_min'].median():.0f}")
        print(f"    carry BOTH orientations (hysteresis): "
              f"{int(ok['has_both'].sum())} of {len(ok)}")
        both = ok[ok["has_both"]]
        if not both.empty:
            print(f"      of those, j at the birth point < j at the collapse "
                  f"point in {int((both['j_at_min'] < both['j_at_max']).sum())} "
                  f"of {len(both)}")
        agr = ok["orient_agree"].sum(), ok["orient_n"].sum()
        print(f"    the two orientation routes agree at {agr[0]} of {agr[1]} "
              "polished candidates")
    return 0


if __name__ == "__main__":
    sys.exit(main())
