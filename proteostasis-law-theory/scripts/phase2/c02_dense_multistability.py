#!/usr/bin/env python
"""phase 2 step 2 -- re-adjudicate every experiment C multistability candidate
with a denser, scale-adaptive search and three independent numerical methods.

which draws are re-examined:

  candidate   every draw experiment C recorded with >= 2 stable equilibria
              (the 68 cases whose status is under audit)
  zero_stable every draw where the attractor search found NO stable equilibrium
              at frac*j_fold even though the fold search had just tracked a
              stable low branch up to j_fold. that combination is internally
              inconsistent and is a solver-failure signature, so it is audited
              in full rather than assumed benign
  control     a seeded random sample of single-attractor draws, to measure the
              FALSE-NEGATIVE rate. re-auditing only the positives would measure
              how many claimed cases survive and say nothing about how many
              real cases were missed

per draw and per evaluation point (j = frac * j_fold, the exact points
experiment C used) the script runs:

  * MINPACK hybr on a dense scale-adaptive multi-start grid;
  * MINPACK lm (levenberg-marquardt) on the same grid;
  * a newton-free nullcline bracketing search (brentq + bisection);
  * pseudo-arclength continuation in j from the lowest stable root;

then certifies every consensus root (cancellation-aware residual, analytic vs
richardson jacobian, conserved-pool balance, knaster-tarski closure uniqueness,
positivity) and sweeps the duplicate-root clustering tolerance.

usage:  PHASE2_AUDIT_DIR=... python scripts/phase2/c02_dense_multistability.py --workers 20
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _audit_common import (auditRoot, loadProvenance, loadSamples, paramsFromRow,  # noqa: E402
                           parallelMap, runRoot, writeJson, writeTable)
from _roots import (certifyRoot, clusterRoots, continueBranch, dedupeSensitivity,  # noqa: E402
                    denseMultiStart, matchRoot, nullclineRoots)

#: burden above which a root sits close enough to the search ceiling that its
#: existence cannot be separated from the box. experiment C searched to 1e5 and
#: the log clip sits at exp(12)=1.6e5.
CEILING_BURDEN = 1.0e4

#: a stable root whose leading eigenvalue is this small relative to the
#: spectral radius is not safely distinguishable from marginal.
MARGINAL_EIG_REL = 1e-8

_CTX: dict = {}


def _auditPoint(p, frac: float, jf: float, anchors, cfg_audit: dict) -> dict:
    """audit one (draw, evaluation point) pair with all four methods."""
    pj = p.with_(j=frac * jf)
    out: dict = dict(frac=frac, j=frac * jf)

    n_grid = cfg_audit["n_grid"]
    raw_h = denseMultiStart(pj, n_grid=n_grid, lo=cfg_audit["lo"], hi=cfg_audit["hi"],
                            method="hybr", anchors=anchors)
    raw_l = denseMultiStart(pj, n_grid=n_grid, lo=cfg_audit["lo"], hi=cfg_audit["hi"],
                            method="lm", anchors=anchors)
    nl, nl_diag = nullclineRoots(pj, u_lo=cfg_audit["lo"], u_hi=cfg_audit["hi"],
                                 a_lo=cfg_audit["lo"], a_hi=cfg_audit["hi"],
                                 n_u=cfg_audit["n_u"], n_a=cfg_audit["n_a"])

    cl_h = clusterRoots(raw_h, 1e-5)
    cl_l = clusterRoots(raw_l, 1e-5)
    cl_n = clusterRoots(nl, 1e-5)

    for tag, rawlist, clist in (("hybr", raw_h, cl_h), ("lm", raw_l, cl_l),
                                ("nullcline", nl, cl_n)):
        out[f"{tag}_n_raw"] = len(rawlist)
        out[f"{tag}_n_eq"] = len(clist)
        out[f"{tag}_n_stable"] = int(sum(r.stable for r in clist))
    out["nullcline_diag"] = nl_diag

    # --- consensus: how many methods independently found each root ---------
    pooled = clusterRoots(list(raw_h) + list(raw_l) + list(nl), 1e-5)
    consensus = []
    for r in pooled:
        votes = [tag for tag, cl in (("hybr", cl_h), ("lm", cl_l), ("nullcline", cl_n))
                 if matchRoot(cl, r.u, r.a, rtol=1e-3) is not None]
        rec = r.asDict()
        rec["methods"] = votes
        rec["n_methods"] = len(votes)
        rec.update(certifyRoot(r.u, r.a, pj))
        consensus.append(rec)
    out["roots"] = consensus

    unan = [r for r in consensus if r["n_methods"] == 3]
    maj = [r for r in consensus if r["n_methods"] >= 2]
    out["n_eq_pooled"] = len(consensus)
    out["n_eq_unanimous"] = len(unan)
    out["n_stable_unanimous"] = int(sum(bool(r["stable"]) for r in unan))
    out["n_stable_majority"] = int(sum(bool(r["stable"]) for r in maj))
    out["n_saddle_unanimous"] = int(sum(r["kind"] == "saddle" for r in unan))

    # --- quality flags -----------------------------------------------------
    stable_unan = [r for r in unan if r["stable"]]
    out["max_residual_cancel_stable"] = (
        float(max(r["residual_cancel"] for r in stable_unan)) if stable_unan else None)
    out["max_jac_rel_error_stable"] = (
        float(max(r["jac_rel_error"] for r in stable_unan
                  if r["jac_rel_error"] is not None)) if stable_unan else None)
    out["max_pool_balance_stable"] = (
        float(max(r["pool_balance_max"] for r in stable_unan)) if stable_unan else None)
    out["all_closures_unique"] = bool(all(r["closure_unique"] for r in unan)) if unan else None
    out["all_positive"] = bool(all(r["positive"] for r in unan)) if unan else None
    out["any_stable_near_ceiling"] = bool(
        any(r["burden"] > CEILING_BURDEN for r in stable_unan))
    out["any_stable_marginal"] = bool(
        any(abs(r["eig_real_max"]) / max(abs(r["eig_real_min"]), 1e-300) < MARGINAL_EIG_REL
            for r in stable_unan))
    out["min_stable_separation_log10"] = None
    if len(stable_unan) >= 2:
        bs = sorted(r["burden"] for r in stable_unan)
        out["min_stable_separation_log10"] = float(
            min(np.log10(b2 / b1) for b1, b2 in zip(bs, bs[1:])))

    out["dedupe_sensitivity"] = dedupeSensitivity(list(raw_h) + list(raw_l) + list(nl))

    # --- pseudo-arclength continuation -------------------------------------
    cont = None
    if stable_unan:
        low = min(stable_unan, key=lambda r: r["burden"])
        try:
            c = continueBranch(pj, low["u"], low["a"], frac * jf, direction=+1,
                               max_steps=cfg_audit["cont_steps"])
            cont = {k: c[k] for k in ("status", "n_points", "n_folds", "j_min",
                                      "j_max", "burden_min", "burden_max")}
            cont["folds"] = c["folds"][:12]
            # does the traced branch actually visit the upper stable root?
            hi_r = max(stable_unan, key=lambda r: r["burden"])
            if len(stable_unan) >= 2:
                near = min((abs(np.log10(q["burden"] / hi_r["burden"]))
                            + abs(q["j"] - frac * jf) / max(frac * jf, 1e-300)
                            for q in c["points"]), default=np.inf)
                cont["reaches_upper_stable"] = bool(near < 0.05)
                cont["closest_approach_to_upper"] = float(near)
        except Exception as exc:                                  # noqa: BLE001
            cont = dict(status=f"error: {type(exc).__name__}: {exc}")
    out["continuation"] = cont
    return out


def _task(item):
    idx, group = item
    cfg = _CTX["cfg"]
    df = _CTX["df"]
    audit = _CTX["audit"]
    row = df.loc[df["sample_index"] == idx].iloc[0]
    rec = dict(sample_index=int(idx), group=group, error="")
    try:
        p = paramsFromRow(row, cfg["base_params"], cfg.get("rescue_total", 1.0))
        jf = float(row["j_fold"])
        rec["j_fold"] = jf
        rec["removal_ceiling"] = float(row["removal_ceiling"])
        rec["phase1_max_stable"] = (None if pd.isna(row["max_stable_equilibria"])
                                    else int(row["max_stable_equilibria"]))
        rec["phase1_max_eq"] = (None if pd.isna(row["max_equilibria"])
                                else int(row["max_equilibria"]))
        rec["phase1_second_stable"] = bool(row["C4_second_stable_attractor"])
        rec["phase1_confirmed"] = bool(row["C4_second_attractor_confirmed"])
        # anchor the multi-start on the scale of the fold equilibrium, which is
        # the only state experiment C recorded per draw
        anchors = []
        fb = row.get("fold_burden", np.nan)
        faf = row.get("fold_aggregate_fraction", np.nan)
        if np.isfinite(fb) and np.isfinite(faf) and fb > 0:
            anchors.append((max(fb * (1.0 - faf), 1e-12), max(fb * faf, 1e-12)))
        rec["points"] = [_auditPoint(p, float(f), jf, anchors, audit)
                         for f in cfg["attractor_fractions"]]
    except Exception as exc:                                      # noqa: BLE001
        import traceback
        rec["error"] = f"{type(exc).__name__}: {exc}\n{traceback.format_exc()[:600]}"
    return rec


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=None)
    ap.add_argument("--audit-dir", default=None)
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--n-controls", type=int, default=200)
    ap.add_argument("--control-seed", type=int, default=20260731)
    ap.add_argument("--n-grid", type=int, default=41)
    ap.add_argument("--n-u", type=int, default=400)
    ap.add_argument("--n-a", type=int, default=1000)
    ap.add_argument("--cont-steps", type=int, default=1500)
    ap.add_argument("--limit", type=int, default=0, help="debug: audit only N draws")
    args = ap.parse_args()

    run = runRoot(args.run)
    out = auditRoot(args.audit_dir)
    df = loadSamples(run)
    cfg = loadProvenance(run)["config"]

    fold_ok = df["C1_fold_exists"].fillna(False).astype(bool)
    cand = df.loc[df["max_stable_equilibria"].fillna(-1) >= 2, "sample_index"].tolist()
    zero = df.loc[fold_ok & (df["max_stable_equilibria"].fillna(-1) == 0),
                  "sample_index"].tolist()
    single_pool = df.loc[fold_ok & (df["max_stable_equilibria"].fillna(-1) == 1),
                         "sample_index"].tolist()
    rng = np.random.default_rng(args.control_seed)
    ctrl = sorted(int(i) for i in rng.choice(single_pool,
                                             size=min(args.n_controls, len(single_pool)),
                                             replace=False))

    items = ([(int(i), "candidate") for i in cand]
             + [(int(i), "zero_stable") for i in zero]
             + [(int(i), "control") for i in ctrl])
    if args.limit:
        items = items[:args.limit]

    audit_cfg = dict(n_grid=args.n_grid, lo=1e-12, hi=1.5e5, n_u=args.n_u,
                     n_a=args.n_a, cont_steps=args.cont_steps,
                     ceiling_burden=CEILING_BURDEN,
                     marginal_eig_rel=MARGINAL_EIG_REL,
                     control_seed=args.control_seed, n_controls=len(ctrl))
    _CTX.update(cfg=cfg, df=df, audit=audit_cfg)

    print(f"[c02] {len(items)} draws "
          f"({len(cand)} candidate / {len(zero)} zero_stable / {len(ctrl)} control), "
          f"{args.workers} workers", flush=True)
    t0 = time.time()
    recs = parallelMap(_task, items, args.workers, chunksize=1)
    print(f"[c02] done in {time.time()-t0:.1f}s", flush=True)

    writeJson(dict(audit_config=audit_cfg, groups=dict(
        candidate=cand, zero_stable=zero, control=ctrl), records=recs),
        out / "c02_dense_multistability.json")

    # flat per-(draw, point) table for the tests and the classifier
    rows = []
    for rec in recs:
        for pt in rec.get("points", []):
            r = {k: v for k, v in rec.items() if k not in ("points",)}
            r.update({k: v for k, v in pt.items()
                      if k not in ("roots", "dedupe_sensitivity", "continuation",
                                   "nullcline_diag")})
            c = pt.get("continuation") or {}
            r.update({f"cont_{k}": v for k, v in c.items() if k != "folds"})
            ds = pt.get("dedupe_sensitivity", {})
            for tol, v in ds.items():
                r[f"dedupe_{tol}_n_stable"] = v["n_stable"]
                r[f"dedupe_{tol}_n_eq"] = v["n_eq"]
            rows.append(r)
    tab = pd.DataFrame(rows)
    writeTable(tab, out / "c02_points.tsv")

    n_err = int((tab["error"].astype(str).str.strip() != "").sum()) if len(tab) else 0
    print(json.dumps(dict(
        n_draws=len(recs), n_points=len(tab), n_errors=n_err,
        candidates_with_2plus_stable_unanimous=int(
            tab.loc[tab.group == "candidate"]
            .groupby("sample_index")["n_stable_unanimous"].max().ge(2).sum()),
        controls_with_2plus_stable_unanimous=int(
            tab.loc[tab.group == "control"]
            .groupby("sample_index")["n_stable_unanimous"].max().ge(2).sum()),
        zero_stable_with_1plus_stable_unanimous=int(
            tab.loc[tab.group == "zero_stable"]
            .groupby("sample_index")["n_stable_unanimous"].max().ge(1).sum()),
    ), indent=2))
    print(f"-> {out/'c02_dense_multistability.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
