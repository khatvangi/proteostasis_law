"""C0: reconcile the Hopf populations, then split every statistic by group.

TWO JOBS, and the first must land before the second means anything.

RECONCILIATION. Four numbers circulate and none of them was labelled:

    2767  fold states solved over the kinetic box
    2766  of those whose low-burden branch traces (1 trace failure)
     108  of those reaching tr J = 0 before the fold
     104  of THOSE whose branch's influx maximum is the fold

The last step is not a filter on quality but on meaning: where the branch's
j-maximum is not the fold, a second singular point terminates the run the walk
returned, so "the crossing precedes the fold" is about a different branch. Those
4 were excluded from the integration and from Section 7's headline, and they are
the whole of the 108-versus-104 gap. 104/2766 = 3.76%, 108/2766 = 3.90%.

THE SPLIT. Task B7 found that `tr J` at the fold is positive in only some of the
crossers. A branch going `tr J < 0 -> > 0 -> < 0` has at least TWO zeros, so
those networks do not describe a precursor to collapse at all: they describe an
INSTABILITY WINDOW interior to the stable branch, opened by one crossing and
closed by another, after which the branch is stable up to a fold that does
terminate a stable state.

    group "terminal"  tr J > 0 at the fold: one crossing, then the fold
    group "window"    tr J < 0 at the fold: two crossings, stability recovered

The integration was run at each network's tr J MAXIMUM. For a window network
that point lies BETWEEN the two crossings. So the escape statistics answer a
question they were not designed for and the answer decides what Section 7 may
say: if the window group largely fails to escape, its crossings are plausibly
supercritical and a stable limit cycle lives in the window -- for those networks
"sustained oscillation" is the correct word. If it escapes at the same rate as
the terminal group, both crossings are subcritical and the window is a hard
instability band with no cycle.

This is a regrouping of data already on disk. Nothing is recomputed.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
for _p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase3"):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

COMPUTED = REPO_ROOT / "data" / "computed"
OUT = COMPUTED / "hopf_groups.json"
OUT_TSV = COMPUTED / "hopf_groups.tsv"


def reconcile() -> tuple[dict, pd.DataFrame, pd.DataFrame]:
    S = pd.read_csv(COMPUTED / "hopf_refined_kinetic_box.tsv", sep="\t")
    T = S[S["traced"] == True]  # noqa: E712
    cross = T[T["tr_max"] >= 0.0]
    clean = cross[cross["fold_is_j_max"] == 1]

    counts = {
        "n_fold_states": int(len(S)),
        "n_traced": int(len(T)),
        "n_trace_failed": int(len(S) - len(T)),
        "n_cross": int(len(cross)),
        "n_cross_clean": int(len(clean)),
        "n_cross_ambiguous_by_multiplicity": int(len(cross) - len(clean)),
        "pct_cross_of_traced": float(100.0 * len(cross) / len(T)),
        "pct_clean_of_traced": float(100.0 * len(clean) / len(T)),
    }
    return counts, T, clean


def groups() -> dict:
    counts, T, clean = reconcile()

    # tr J at the SOLVED fold state, from the B3/B7 reduction
    G = pd.read_csv(COMPUTED / "genericity_structure.tsv", sep="\t")
    G = G[G["population"] == "kinetic_box"][["name", "tr_J", "det_J", "grad_G"]]
    G = G.rename(columns={"tr_J": "tr_J_at_fold"})

    C = clean.merge(G, on="name", how="left")
    n_missing = int(C["tr_J_at_fold"].isna().sum())
    C = C[C["tr_J_at_fold"].notna()].copy()
    C["group"] = np.where(C["tr_J_at_fold"] > 0.0, "terminal", "window")

    # the same split over ALL crossers, not only the clean ones
    A = T[T["tr_max"] >= 0.0].merge(G, on="name", how="left")
    A = A[A["tr_J_at_fold"].notna()]
    all_split = {
        "terminal": int((A["tr_J_at_fold"] > 0.0).sum()),
        "window": int((A["tr_J_at_fold"] <= 0.0).sum()),
        "n": int(len(A)),
    }

    # integration, run on the clean crossers ("+") and 205 controls ("-")
    I = pd.read_csv(COMPUTED / "hopf_integration_kinetic_box.tsv", sep="\t")
    I = I[I["tested"] == True].copy()  # noqa: E712
    I["is_cross"] = I["name"].str.startswith("+")
    I["name"] = I["name"].str.slice(1)
    X = I[I["is_cross"]].merge(C[["name", "group", "tr_J_at_fold", "j_crit",
                                  "j_at_first_cross", "tr_max"]],
                               on="name", how="left")
    n_unmatched = int(X["group"].isna().sum())

    def stats(sub: pd.DataFrame) -> dict:
        if sub.empty:
            return {"n": 0}
        fit = sub[sub["slope"].notna()]
        rel = ((fit["slope"] - fit["lambda_max"]).abs()
               / fit["lambda_max"].abs().clip(lower=1e-12)) if not fit.empty else None
        per = sub[sub["period_measured"].notna()]
        pr = ((per["period_measured"] - per["period_predicted"]).abs()
              / per["period_predicted"]) if not per.empty else None
        return {
            "n": int(len(sub)),
            "n_escaped": int(sub["escaped"].sum()),
            "frac_escaped": float(sub["escaped"].mean()),
            "n_grew_10x": int(sub["grew"].sum()),
            "amplification_median": float(sub["ratio_max"].median()),
            "amplification_p10": float(sub["ratio_max"].quantile(0.10)),
            "amplification_max": float(sub["ratio_max"].max()),
            "n_horizon_short": int(sub["horizon_short"].sum()),
            "n_complex_pair": int(sub["complex_pair"].sum()),
            "exponent_n_evaluable": int(len(fit)),
            "exponent_within_5pct": int((rel < 0.05).sum()) if rel is not None else 0,
            "period_n_measurable": int(len(per)),
            "period_within_5pct": int((pr < 0.05).sum()) if pr is not None else 0,
            "first_cross_frac_median": float(
                (sub["j_at_first_cross"] / sub["j_crit"]).median())
            if "j_crit" in sub else None,
            "tr_max_median": float(sub["tr_max"].median()) if "tr_max" in sub else None,
        }

    controls = I[~I["is_cross"]]
    out = {
        "reconciliation": counts,
        "n_clean_without_fold_trace": n_missing,
        "n_integrated_without_group": n_unmatched,
        "split_over_all_crossers": all_split,
        "split_over_clean_crossers": {
            "terminal": int((C["group"] == "terminal").sum()),
            "window": int((C["group"] == "window").sum()),
        },
        "integration": {
            "terminal": stats(X[X["group"] == "terminal"]),
            "window": stats(X[X["group"] == "window"]),
            "all_crossers": stats(X),
            "controls": {
                "n": int(len(controls)),
                "n_escaped": int(controls["escaped"].sum()),
                "amplification_median": float(controls["ratio_max"].median()),
            },
        },
    }
    X.to_csv(OUT_TSV, sep="\t", index=False)
    return out


# ---------------------------------------------------------------------------
# C0b tasks 1 and 3: how many times does tr J cross, and how wide is the window
# ---------------------------------------------------------------------------


def _zeroWorker(item):
    """re-trace one crossing network and count the zeros of tr J on its branch.

    The parity is the check. A branch ending with tr J > 0 must cross an ODD
    number of times; one ending with tr J < 0, having crossed at all, must cross
    an EVEN number >= 2. A network violating its own parity has been mis-traced,
    and is reported as such rather than counted.
    """
    import hopf_check as HC
    import fold_theorem as FT2
    from proteostasis import model as M2

    name, p, u, a = item
    try:
        out = HC.branchProfile(p, u, a, n=150)
    except Exception:
        return {"name": name, "ok": False}
    if out is None:
        return {"name": name, "ok": False}

    B = out["branch"]
    # order along the branch by influx: the accessible branch runs from low
    # burden up to the fold, and j is what the experiment varies.
    B = B.sort_values("j").reset_index(drop=True)
    tr = B["tr_J"].to_numpy(float)
    j = B["j"].to_numpy(float)
    sign = np.sign(tr)
    idx = np.nonzero(sign[:-1] * sign[1:] < 0)[0]

    j_cross = []
    for k in idx:
        # linear interpolation in j at the sign change
        t0, t1 = tr[k], tr[k + 1]
        w = abs(t0) / max(abs(t0) + abs(t1), 1e-300)
        j_cross.append(float(j[k] + w * (j[k + 1] - j[k])))

    try:
        j_crit = float(FT2.removalR(u, a, p))
        J = M2.jacobian(u, a, p)
        tr_fold = float(J[0, 0] + J[1, 1])
    except Exception:
        return {"name": name, "ok": False}

    return {
        "name": name, "ok": True,
        "n_zeros": int(len(j_cross)),
        "tr_at_fold": tr_fold,
        "tr_at_branch_end": float(tr[-1]),
        "j_crit": j_crit,
        "j_H1": j_cross[0] if j_cross else np.nan,
        "j_H2": j_cross[1] if len(j_cross) > 1 else np.nan,
        "j_last_cross": j_cross[-1] if j_cross else np.nan,
        "n_points": int(len(B)),
    }


def zeroCounts(workers: int | None = None) -> dict:
    """re-trace only the crossers -- 108 networks, not 2767."""
    import multiprocessing as mp
    import os
    import fold_theorem as FT2
    import genericity as GEN

    counts, T, _clean = reconcile()
    cross_names = set(T[T["tr_max"] >= 0.0]["name"])

    run = FT2.phase1RunDir()
    c = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    c = c[pd.to_numeric(c["fold_burden"], errors="coerce").notna()]

    items = []
    for i, r in c.iterrows():
        nm = f"draw{i}"
        if nm not in cross_names:
            continue
        try:
            p = FT2.paramsFromSampleRow(r)
            u0, a0 = FT2.foldStateFromSampleRow(r)
        except Exception:
            continue
        s = GEN.polishFold(p, u0, a0) or FT2.foldSolve(p)
        if s is None:
            continue
        items.append((nm, p, float(s[1]), float(s[2])))

    if workers is None:
        workers = max(1, min(56, (os.cpu_count() or 2) - 8))
    with mp.get_context("fork").Pool(processes=workers) as pool:
        rows = pool.map(_zeroWorker, items, chunksize=2)

    Z = pd.DataFrame(rows)
    ok = Z[Z["ok"] == True].copy()  # noqa: E712
    ok["group"] = np.where(ok["tr_at_fold"] > 0.0, "terminal", "window")
    ok["parity_ok"] = np.where(
        ok["group"] == "terminal", ok["n_zeros"] % 2 == 1,
        (ok["n_zeros"] % 2 == 0) & (ok["n_zeros"] >= 2))
    ok["window_width"] = (ok["j_H2"] - ok["j_H1"]) / ok["j_crit"]
    ok.to_csv(COMPUTED / "hopf_zero_counts.tsv", sep="\t", index=False)

    out = {"n_requested": len(items), "n_traced": int(len(ok)),
           "n_trace_failed": int(len(Z) - len(ok))}
    for g in ("terminal", "window"):
        sub = ok[ok["group"] == g]
        if sub.empty:
            out[g] = {"n": 0}
            continue
        d = {
            "n": int(len(sub)),
            "n_zeros_min": int(sub["n_zeros"].min()),
            "n_zeros_median": float(sub["n_zeros"].median()),
            "n_zeros_max": int(sub["n_zeros"].max()),
            "n_parity_violations": int((~sub["parity_ok"]).sum()),
            "parity_violation_names": sub[~sub["parity_ok"]]["name"].tolist()[:20],
            "first_cross_over_jcrit_median": float(
                (sub["j_H1"] / sub["j_crit"]).median()),
        }
        if g == "window":
            w = sub["window_width"].dropna()
            d.update({
                "n_with_two_crossings": int(len(w)),
                "window_width_min": float(w.min()) if len(w) else None,
                "window_width_median": float(w.median()) if len(w) else None,
                "window_width_max": float(w.max()) if len(w) else None,
                "j_H2_over_jcrit_median": float(
                    (sub["j_H2"] / sub["j_crit"]).median()),
            })
        out[g] = d
    return out


def main() -> int:
    o = groups()
    o["zero_counts"] = zeroCounts()
    OUT.write_text(json.dumps(o, indent=2) + "\n")
    print(json.dumps(o, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
