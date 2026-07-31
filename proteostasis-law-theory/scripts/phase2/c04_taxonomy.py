#!/usr/bin/env python
"""phase 2 step 4 -- sort every audited draw into the five outcome classes and
report the fractions against both denominators.

the classes, and what separates them:

  confirmed_multistable          >=2 stable roots found INDEPENDENTLY by all
                                 three finders, each surviving certification,
                                 separated in burden, and with demonstrated
                                 basin access (two nonempty grid basins, or a
                                 saddle whose unstable manifold reaches two
                                 different stable roots)
  locally_stable_no_basin        >=2 certified stable roots, but no trajectory
                                 evidence that the second one is reachable. it
                                 is a real root; it is not a demonstrated
                                 attractor
  root_finder_artifact           experiment C reported >=2 stable, the dense
                                 three-method search does not reproduce it, OR
                                 the "two" roots merge under a looser but still
                                 tight clustering tolerance
  threshold_only                 the second stable root sits within a decade of
                                 the search ceiling, where its existence cannot
                                 be separated from the box boundary
  numerical_failure              certification failed: non-unique rapid-
                                 equilibrium closure, broken pool balance,
                                 analytic/numerical jacobian disagreement,
                                 non-positive state, or an exception

precedence when a draw's two evaluation points disagree is generous to
multistability on purpose: experiment C itself took the maximum over evaluation
points, so the audit must not gain an advantage by taking the minimum.

usage:  PHASE2_AUDIT_DIR=... python scripts/phase2/c04_taxonomy.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _audit_common import auditRoot, loadSamples, runRoot, writeJson, writeTable  # noqa: E402

CEILING_BURDEN = 1.0e4

#: certification thresholds. a stable root must clear all of them.
CERT = dict(residual_cancel=1e-9, jac_rel_error=1e-4, pool_balance=1e-10)

#: the clustering tolerance at which two "distinct" stable roots merging is
#: treated as evidence they were one root all along. 1e-3 is a hundredfold
#: looser than experiment C's 1e-5 and still far tighter than any real
#: separation between physically distinct attractors.
MERGE_TOL = "0.001"

CLASSES = ("confirmed_multistable", "locally_stable_no_basin", "threshold_only",
           "root_finder_artifact", "numerical_failure", "single_attractor")

#: most-to-least multistable; used to reduce a draw's per-point classes to one
PRECEDENCE = {c: i for i, c in enumerate(CLASSES)}


def clopperPearson(k: int, n: int, alpha: float = 0.05):
    """exact binomial interval. used because the counts here are small enough
    that a normal approximation would give a lower bound below zero."""
    from scipy.stats import beta
    if n == 0:
        return (float("nan"), float("nan"))
    lo = 0.0 if k == 0 else float(beta.ppf(alpha / 2, k, n - k + 1))
    hi = 1.0 if k == n else float(beta.ppf(1 - alpha / 2, k + 1, n - k))
    return (lo, hi)


def certifiedStable(pt: dict) -> list:
    """the stable roots that all three finders found and that pass every check."""
    out = []
    for r in pt["roots"]:
        if r["n_methods"] < 3 or not r["stable"]:
            continue
        if not r.get("positive", False) or not r.get("closure_unique", False):
            continue
        if r.get("residual_cancel", np.inf) > CERT["residual_cancel"]:
            continue
        je = r.get("jac_rel_error")
        if je is None or je > CERT["jac_rel_error"]:
            continue
        if r.get("pool_balance_max", np.inf) > CERT["pool_balance"]:
            continue
        out.append(r)
    return sorted(out, key=lambda r: r["burden"])


def failsCertification(pt: dict) -> bool:
    """any unanimous stable root that exists but does not survive the checks."""
    unan_stable = [r for r in pt["roots"] if r["n_methods"] == 3 and r["stable"]]
    return len(unan_stable) > len(certifiedStable(pt))


def classifyPoint(rec: dict, pt: dict, basin: dict | None) -> tuple[str, dict]:
    ev: dict = {}
    if rec.get("error"):
        return "numerical_failure", dict(reason="task error")

    cs = certifiedStable(pt)
    ev["n_stable_unanimous"] = pt["n_stable_unanimous"]
    ev["n_stable_certified"] = len(cs)
    ev["n_stable_majority"] = pt["n_stable_majority"]
    ev["burdens"] = [r["burden"] for r in cs]

    if failsCertification(pt):
        return "numerical_failure", dict(ev, reason="unanimous stable root failed certification")

    if len(cs) < 2:
        if rec.get("phase1_max_stable", 0) and rec["phase1_max_stable"] >= 2:
            return "root_finder_artifact", dict(
                ev, reason="phase 1 reported >=2 stable; three-method dense search "
                           "certifies fewer")
        return "single_attractor", ev

    # two or more certified stable roots from here on
    ds = pt.get("dedupe_sensitivity", {})
    merged = ds.get(MERGE_TOL, {}).get("n_stable")
    ev["n_stable_at_merge_tol"] = merged
    if merged is not None and merged < 2:
        return "root_finder_artifact", dict(
            ev, reason=f"stable roots merge at clustering tolerance {MERGE_TOL}")

    if any(r["burden"] > CEILING_BURDEN for r in cs):
        return "threshold_only", dict(
            ev, reason=f"a stable root sits above burden {CEILING_BURDEN:g}, within "
                       "a decade of the search ceiling")

    if basin is None or basin.get("error"):
        return "locally_stable_no_basin", dict(ev, reason="no basin scan available")

    ev["n_nonempty_grid_basins"] = basin["n_stable_with_nonempty_grid_basin"]
    ev["saddle_separates"] = basin["saddle_separates"]
    ev["distinct_stable_reached_grid"] = basin["distinct_stable_reached_grid"]
    ev["n_saddle"] = basin["n_saddle"]
    if basin["any_positivity_violation"]:
        return "numerical_failure", dict(ev, reason="trajectory left the nonnegative orthant")
    if basin["n_stable_with_nonempty_grid_basin"] >= 2 or basin["saddle_separates"]:
        return "confirmed_multistable", dict(ev, reason="two nonempty basins and/or a "
                                                        "separating saddle")
    return "locally_stable_no_basin", dict(
        ev, reason="certified stable roots, but no trajectory reached the second one")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=None)
    ap.add_argument("--audit-dir", default=None)
    args = ap.parse_args()

    run = runRoot(args.run)
    out = auditRoot(args.audit_dir)
    df = loadSamples(run)

    with open(out / "c02_dense_multistability.json") as fh:
        c02 = json.load(fh)
    basins = {}
    bpath = out / "c03_basins.json"
    if bpath.exists():
        with open(bpath) as fh:
            for r in json.load(fh)["records"]:
                basins[(int(r["sample_index"]), float(r["frac"]))] = r

    rows = []
    for rec in c02["records"]:
        idx = int(rec["sample_index"])
        per_point = []
        for pt in rec.get("points", []):
            b = basins.get((idx, float(pt["frac"])))
            cls, ev = classifyPoint(rec, pt, b)
            per_point.append(dict(frac=pt["frac"], cls=cls, evidence=ev))
        if not per_point:
            per_point = [dict(frac=None, cls="numerical_failure",
                              evidence=dict(reason="no evaluation point completed"))]
        best = min(per_point, key=lambda q: PRECEDENCE[q["cls"]])
        rows.append(dict(
            sample_index=idx, group=rec["group"], cls=best["cls"],
            reason=best["evidence"].get("reason", ""),
            phase1_max_stable=rec.get("phase1_max_stable"),
            phase1_second_stable=rec.get("phase1_second_stable"),
            phase1_confirmed=rec.get("phase1_confirmed"),
            j_fold=rec.get("j_fold"),
            per_point=json.dumps(per_point),
        ))
    tab = pd.DataFrame(rows)
    writeTable(tab, out / "c04_taxonomy.tsv")

    # ---- counts ----------------------------------------------------------
    by_group = {g: {c: int((sub["cls"] == c).sum()) for c in CLASSES}
                for g, sub in tab.groupby("group")}

    n_total = int(len(df))
    n_fold = int(df["C1_fold_exists"].fillna(False).astype(bool).sum())
    n_single_pool = int((df["C1_fold_exists"].fillna(False).astype(bool)
                         & (df["max_stable_equilibria"].fillna(-1) == 1)).sum())

    cand = tab[tab.group == "candidate"]
    ctrl = tab[tab.group == "control"]
    zero = tab[tab.group == "zero_stable"]

    n_conf = int((cand["cls"] == "confirmed_multistable").sum())
    n_fn = int((ctrl["cls"].isin(("confirmed_multistable", "locally_stable_no_basin",
                                  "threshold_only"))).sum())
    fn_lo, fn_hi = clopperPearson(n_fn, len(ctrl))

    report = dict(
        classes=list(CLASSES),
        certification_thresholds=CERT,
        ceiling_burden=CEILING_BURDEN,
        merge_tolerance=MERGE_TOL,
        denominators=dict(all_draws=n_total, fold_evaluable=n_fold,
                          single_attractor_pool=n_single_pool),
        by_group=by_group,
        candidates=dict(
            n_audited=int(len(cand)),
            n_phase1_second_stable=int(cand["phase1_second_stable"].sum()),
            n_phase1_confirmed=int(cand["phase1_confirmed"].sum()),
            counts={c: int((cand["cls"] == c).sum()) for c in CLASSES},
        ),
        controls=dict(
            n_audited=int(len(ctrl)),
            counts={c: int((ctrl["cls"] == c).sum()) for c in CLASSES},
            false_negative_rate=float(n_fn / len(ctrl)) if len(ctrl) else None,
            false_negative_ci95=[fn_lo, fn_hi],
            implied_missed_in_single_pool=[fn_lo * n_single_pool, fn_hi * n_single_pool],
        ),
        zero_stable=dict(
            n_audited=int(len(zero)),
            counts={c: int((zero["cls"] == c).sum()) for c in CLASSES},
        ),
        headline=dict(
            confirmed_multistable_n=n_conf,
            frac_of_all_draws=n_conf / n_total,
            frac_of_fold_evaluable=n_conf / n_fold,
            ci95_of_fold_evaluable=list(clopperPearson(n_conf, n_fold)),
            phase1_reported_second_stable=int(
                df["C4_second_stable_attractor"].fillna(False).astype(bool).sum()),
            phase1_reported_confirmed=int(
                df["C4_second_attractor_confirmed"].fillna(False).astype(bool).sum()),
        ),
        note=("candidate counts are exhaustive: every draw experiment C reported "
              "with >=2 stable equilibria was re-audited. control counts are a "
              "seeded sample of the single-attractor pool, so the false-negative "
              "rate is an estimate with the interval given, not a census."),
    )
    writeJson(report, out / "c04_taxonomy.json")
    print(json.dumps({k: report[k] for k in
                      ("by_group", "candidates", "controls", "zero_stable", "headline")},
                     indent=2))
    print(f"-> {out/'c04_taxonomy.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
