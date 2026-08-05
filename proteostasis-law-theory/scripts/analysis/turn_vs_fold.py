"""R1: do `j_turn` and `j_crit` lie in a fixed order, or only at one parameter set?

Section 3.1 says "The two lie in a fixed order and differ by one part in a
thousand." Both halves were read off ONE branch, the base parameter set.

They are different objects and nothing obviously relates them:

    j_turn   where the nullcline runs vertical in the burden plane, G_a = 0.
             Cramer gives du*/dj = -G_a/det J, so this is where the soluble
             coordinate has a horizontal tangent. det J != 0 there; it is a
             regular point of the branch.
    j_crit   where det J = 0, so R is stationary along the nullcline. Both
             coordinates turn together.

This walks every load-grid branch, locates each, and reports the ordering and the
separation over the whole population rather than over one member of it.
"""

from __future__ import annotations

import json
import multiprocessing as mp
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

from proteostasis import model as M  # noqa: E402
import fold_theorem as FT  # noqa: E402

OUT = REPO_ROOT / "data" / "computed" / "turn_vs_fold.json"
OUT_TSV = REPO_ROOT / "data" / "computed" / "turn_vs_fold.tsv"


def _worker(item):
    import hopf_check as HC
    name, p, u_star, a_star = item
    try:
        out = HC.branchProfile(p, u_star, a_star, n=150)
    except Exception:
        return {"name": name, "ok": False}
    if out is None:
        return {"name": name, "ok": False}

    B = out["branch"].sort_values("j").reset_index(drop=True)
    ga = []
    for _, r in B.iterrows():
        try:
            _, g_a = FT._centralGradient(FT.aggregateG, r["u"], r["a"], p)
        except (M.ModelError, OverflowError):
            g_a = np.nan
        ga.append(g_a)
    B = B.assign(G_a=ga).dropna(subset=["G_a"])
    if len(B) < 8:
        return {"name": name, "ok": False}

    g = B["G_a"].to_numpy(float)
    j = B["j"].to_numpy(float)
    idx = np.nonzero(np.sign(g[:-1]) * np.sign(g[1:]) < 0)[0]

    try:
        j_crit = float(FT.removalR(u_star, a_star, p))
    except Exception:
        return {"name": name, "ok": False}

    if idx.size == 0:
        return {"name": name, "ok": True, "n_turns": 0, "j_crit": j_crit,
                "j_turn": np.nan, "ratio": np.nan}

    # the turn nearest the fold, interpolated at the sign change
    js = []
    for k in idx:
        w = abs(g[k]) / max(abs(g[k]) + abs(g[k + 1]), 1e-300)
        js.append(j[k] + w * (j[k + 1] - j[k]))
    js = np.array(js, float)
    j_turn = float(js[np.argmin(np.abs(js - j_crit))])

    return {"name": name, "ok": True, "n_turns": int(idx.size),
            "j_crit": j_crit, "j_turn": j_turn,
            "ratio": float(j_turn / j_crit) if j_crit else np.nan}


def run(workers: int | None = None) -> dict:
    import genericity as GEN
    states = GEN.loadGridStates(FT.phase1RunDir())
    if workers is None:
        workers = max(1, min(56, (os.cpu_count() or 2) - 8))
    with mp.get_context("fork").Pool(processes=workers) as pool:
        rows = pool.map(_worker, states, chunksize=2)

    D = pd.DataFrame(rows)
    D.to_csv(OUT_TSV, sep="\t", index=False)
    ok = D[D["ok"] == True]  # noqa: E712
    have = ok[ok["j_turn"].notna()]

    out = {
        "n_networks": int(len(states)),
        "n_traced": int(len(ok)),
        "n_trace_failed": int(len(D) - len(ok)),
        "n_with_a_turn": int(len(have)),
        "n_without_a_turn": int((ok["n_turns"] == 0).sum()),
    }
    if not have.empty:
        r = have["ratio"]
        out.update({
            "n_turn_below_crit": int((r < 1.0).sum()),
            "n_turn_above_crit": int((r > 1.0).sum()),
            "ratio_min": float(r.min()),
            "ratio_median": float(r.median()),
            "ratio_max": float(r.max()),
            "separation_p99": float((1.0 - r).abs().quantile(0.99)),
            "separation_max": float((1.0 - r).abs().max()),
            "n_multiple_turns": int((have["n_turns"] > 1).sum()),
        })
    return out


def main() -> int:
    o = run()
    OUT.write_text(json.dumps(o, indent=2) + "\n")
    print(json.dumps(o, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
