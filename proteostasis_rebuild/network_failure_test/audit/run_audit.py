"""AUDIT 2026-10-09 (new work). runs every audit part and writes audit/results/.
the original results/ directory is read, never written.

usage: python run_audit.py [part ...]     parts: repro struct cf pareto feedback failed
       NPROC (default 16) workers; every part is deterministic in its seeds.
"""
import json
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))

import cf_model as C  # noqa: E402
import feedback as fbm  # noqa: E402
import model as M  # noqa: E402
import run_all as RA  # noqa: E402
from pools import findPools  # noqa: E402

ORIG = HERE.parent / "results"
OUT = HERE / "results"
NPROC = int(os.environ.get("NPROC", 16))
P = findPools()


# ------------------------------------------------------------ repro
REPRO_SETTINGS = ["bio_measured_W2", "bio_blind_W0", "unstructured_measured_W2",
                  "bio_measured_W2_alpha0.3"]
REPRO_SEEDS = range(0, 30)


def runRepro(job):
    RA.initWorker(P)
    return RA.runInstance(job)


def partRepro(pool):
    """rerun original runInstance on a seed sample; compare every column."""
    st = {s["name"]: s for s in RA.SETTINGS}
    jobs = [(st[n], s) for n in REPRO_SETTINGS for s in REPRO_SEEDS]
    new = pd.DataFrame(pool.map(runRepro, jobs))
    old = pd.read_csv(ORIG / "instances.csv", low_memory=False)
    old = old.set_index(["setting", "seed"])
    rows = []
    for _, r in new.iterrows():
        o = old.loc[(r.setting, r.seed)]
        nbad, worst = 0, 0.0
        for col in new.columns:
            if col in ("setting", "seed"):
                continue
            a, b = r[col], o.get(col, np.nan)
            if isinstance(a, str) or isinstance(b, str):
                same = (str(a) if a == a else "") == (str(b) if b == b else "")
            else:
                same = (pd.isna(a) and pd.isna(b)) or np.isclose(a, b, rtol=1e-9, atol=1e-12)
                if not pd.isna(a) and not pd.isna(b) and b != 0:
                    worst = max(worst, abs(a / b - 1))
            nbad += not same
        rows.append({"setting": r.setting, "seed": r.seed, "n_columns": len(new.columns) - 2,
                     "n_mismatch": nbad, "max_rel_diff_numeric": worst,
                     "failed_new": r.get("failed", "")})
    pd.DataFrame(rows).to_csv(OUT / "repro_instances.csv", index=False)

    # feedback: first 8 seeds of each kind
    fold = pd.read_csv(ORIG / "feedback.csv").set_index(["kind", "seed"])
    fj = [(k, s) for k in ("none", "titration", "utilization") for s in range(8)]
    fnew = pool.map(runFeedbackRepro, fj)
    frows = []
    for r in fnew:
        o = fold.loc[(r["kind"], r["seed"])]
        frows.append({"kind": r["kind"], "seed": r["seed"],
                      "n_jumps_same": int(r["n_jumps"] == o.n_jumps),
                      "onset_same": int(np.isclose(r["onset_frac"], o.onset_frac)),
                      "hyst_rel_diff": abs(r["hysteresis"] - o.hysteresis) / max(o.hysteresis, 1e-12)})
    pd.DataFrame(frows).to_csv(OUT / "repro_feedback.csv", index=False)

    # pareto: exact rerun
    RA.initWorker(P)
    seed, _inst, pts = RA.findParetoInstance()
    po = pd.read_csv(ORIG / "pareto_front.csv")
    pn = pd.DataFrame(pts)
    same = bool(seed == 1 and np.allclose(po.values, pn.values, rtol=1e-9))
    return {"repro_instances": {"n": len(rows), "n_with_any_mismatch": int(sum(r["n_mismatch"] > 0 for r in rows)),
                                "max_rel_diff": float(max(r["max_rel_diff_numeric"] for r in rows))},
            "repro_feedback": {"n": len(frows), "n_jumps_all_same": bool(all(r["n_jumps_same"] for r in frows)),
                               "onset_all_same": bool(all(r["onset_same"] for r in frows)),
                               "max_hyst_rel_diff": float(max(r["hyst_rel_diff"] for r in frows))},
            "repro_pareto_identical": same}


def runFeedbackRepro(job):
    RA.initWorker(P)
    return RA.runFeedback(job)


# ------------------------------------------------------------ struct (no LP)
def partStruct():
    """structural single cuts depend on the route list only: no capacity,
    abundance, kinetics or load enters. large Monte Carlo of the generators."""
    N = 20000
    cfgs = [("orig", {}), ("orig_aggfree", {"agg_free": True}), ("nocut", {}),
            ("unstruct", {})] + [(f"density_q{q}", {"q": q}) for q in
                                 (0.0, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 1.0)]
    rows, mach = [], []
    for name, kw in cfgs:
        mode = name.split("_")[0]
        cnt = {m: 0 for m in M.MACHINES}
        per_class = {k: 0 for k in M.CLASSES}
        anycut = 0
        for s in range(N):
            rs, rx = np.random.default_rng(s), np.random.default_rng(s + 5_000_000)
            base, extra = C.findRoutes(rs, mode, kw.get("q"), kw.get("agg_free", False), rx)
            R = base + extra
            cut = C.findCutMachines(R)
            anycut += bool(cut)
            for m in cut:
                cnt[m] += 1
            for k in M.CLASSES:
                if C.findCutMachines([r for r in R if r[0] == k]):
                    per_class[k] += 1
        rows.append({"config": name, "n": N, "frac_any_single_cut": anycut / N,
                     **{f"frac_class_{k}_has_cut": v / N for k, v in per_class.items()}})
        mach += [{"config": name, "machine": m, "frac_is_single_cut": v / N} for m, v in cnt.items()]
    pd.DataFrame(rows).to_csv(OUT / "struct_single_cuts.csv", index=False)
    pd.DataFrame(mach).to_csv(OUT / "struct_cut_machines.csv", index=False)

    # analytic check of the original rule set: no single cut iff
    # GRO has the KJE alternative (0.2) AND HTP primary is HTPG alone with a KJE
    # alternative (0.5*0.5) AND DK has GRO or TF alternative (1-0.5*0.7)
    # AND DMG has >= 2 proteases AND AGG has a ClpB-free route (SEQ or LON)
    p_dmg2 = 1 - (0.3 ** 3 + 3 * 0.7 * 0.3 ** 2)
    p_agg = 1 - 0.3 * 0.7
    analytic = 1 - 0.2 * 0.25 * (1 - 0.5 * 0.7) * p_dmg2 * p_agg
    # independence across classes holds: each rule uses its own draw
    return {"orig_frac_any_single_cut_MC": rows[0]["frac_any_single_cut"],
            "orig_frac_any_single_cut_analytic": analytic,
            "orig_frac_single_cut_reported_by_original_run": 0.987}


# ------------------------------------------------------------ cf (LP ensembles)
NCF = 300
CF_SETTINGS = [
    ("orig_W2", dict()),
    ("orig_W2_kinresample", dict(kin_offset=100_000)),
    ("orig_W0", dict(W=0.0)),
    ("orig_W4", dict(W=4.0)),
    ("orig_blind_W2", dict(abundance="blind")),
    ("orig_shuffled_W0", dict(abundance="shuffled", W=0.0)),
    ("orig_aggfree_W2", dict(agg_free=True)),
    ("orig_aggshare0.01_W2", dict(agg_share=0.01)),
    ("nocut_W2", dict(routes="nocut")),
    ("nocut_aggshare0.01_W2", dict(routes="nocut", agg_share=0.01)),
    ("nocut_blind_W2", dict(routes="nocut", abundance="blind")),
    ("density_q0.2_W2", dict(routes="density", q=0.2)),
    ("density_q0.5_W2", dict(routes="density", q=0.5)),
    ("density_q0.8_W2", dict(routes="density", q=0.8)),
    ("density_q1.0_W2", dict(routes="density", q=1.0)),
    ("unstruct_W2", dict(routes="unstruct")),
    ("orig_W2_costflat", dict(cost="flat")),
    ("orig_W2_costreverse", dict(cost="reverse")),
]


def runCF(job):
    name, kw, seed = job
    kw = dict(kw)
    off = kw.pop("kin_offset", 0)
    try:
        inst = C.makeInstanceCF(seed, pools=P, kin_seed=seed + off, **kw)
        row = C.evalInstance(inst)
        row["failed"] = ""
    except Exception as e:
        row = {"failed": repr(e)[:200]}
    row.update({"setting": name, "seed": seed})
    return row


def partCF(pool):
    jobs = [(n, kw, s) for n, kw in CF_SETTINGS for s in range(NCF)]
    df = pd.DataFrame(pool.map(runCF, jobs, chunksize=2))
    df.to_csv(OUT / "cf_instances.csv", index=False)
    out = []
    for n, g in df.groupby("setting", sort=False):
        ok = g[g.failed == ""]
        e = {"setting": n, "n": len(g), "n_failed": int((g.failed != "").sum()),
             "frac_single_cut": float((ok.single_cuts > 0).mean()),
             "median_single_cuts": float(ok.single_cuts.median())}
        for rho in (0.25, 0.5, 0.9):
            e[f"frac_SL_any_rho{rho}"] = float(((ok[f"SL_struct_rho{rho}"] + ok[f"SL_load_rho{rho}"]) > 0).mean())
            e[f"frac_SL_struct_rho{rho}"] = float((ok[f"SL_struct_rho{rho}"] > 0).mean())
            e[f"frac_SL_load_rho{rho}"] = float((ok[f"SL_load_rho{rho}"] > 0).mean())
        e["frac_sat_events_ge2"] = float((ok.sat_events >= 2).mean())
        e["mean_sat_events"] = float(ok.sat_events.mean())
        bn = ok.bottlenecks.fillna("").str.split(";").explode()
        e["frac_bottleneck_ClpB"] = float(ok.bottlenecks.fillna("").str.contains("ClpB").mean())
        e["top_bottlenecks"] = ";".join(f"{k}:{v / len(ok):.2f}" for k, v in bn[bn != ""].value_counts().head(4).items())
        out.append(e)
    pd.DataFrame(out).to_csv(OUT / "cf_summary.csv", index=False)

    # paired comparisons on the same structures (orig_W2 vs kinetic resample;
    # orig vs flat / reverse costs): how much of each statistic is structure
    a = df[df.setting == "orig_W2"].set_index("seed")
    paired = {}
    for other in ("orig_W2_kinresample", "orig_W2_costflat", "orig_W2_costreverse", "orig_aggfree_W2"):
        b = df[df.setting == other].set_index("seed")
        paired[other] = {
            "frac_same_SL_load_rho0.5_presence": float((a["SL_load_rho0.5"] > 0).eq(b["SL_load_rho0.5"] > 0).mean()),
            "frac_same_sat_events": float(a.sat_events.eq(b.sat_events).mean()),
            "frac_same_first_saturating_machine": float(
                a.sat_order.fillna("").str.split("@").str[0].eq(b.sat_order.fillna("").str.split("@").str[0]).mean()),
            "frac_same_bottleneck_set": float(a.bottlenecks.eq(b.bottlenecks).mean()),
        }
    return {"paired": paired}


# ------------------------------------------------------------ pareto
def runPareto(seed):
    inst = M.makeInstance(seed, W=2.0, pools=P)
    pts = M.findParetoFront(inst, 0.5, n=40)
    atp = np.array([p["atp"] for p in pts])
    seq = np.array([p["sequestered_flux"] for p in pts])
    t = np.array([p["max_utilization"] for p in pts])
    # strictly non-dominated in (min atp, max reserve) among the computed points
    nd = sum(1 for i in range(len(pts)) if not any(
        (atp[j] <= atp[i] and t[j] <= t[i]) and (atp[j] < atp[i] - 1e-9 * atp[i] or t[j] < t[i])
        for j in range(len(pts))))
    knee = t[np.argmax(atp <= atp[-1] * (1 + 1e-9))]
    return {"seed": seed, "n_points": len(pts), "atp_max_reserve": atp[0], "atp_min": atp[-1],
            "rel_atp_span": atp[0] / atp[-1] - 1 if atp[-1] > 0 else np.nan,
            "n_nondominated": nd, "util_at_atp_min": float(knee),
            "seq_change": float(seq.max() - seq.min()),
            "has_seq_route": int(any(k == "AGG" and "SEQ" in m for k, m in inst["routes"]))}


def partPareto(pool):
    df = pd.DataFrame(pool.map(runPareto, range(NCF)))
    df.to_csv(OUT / "pareto_ensemble.csv", index=False)
    deg = df.rel_atp_span < 1e-6
    return {"n": len(df), "frac_degenerate_no_tradeoff": float(deg.mean()),
            "median_rel_atp_span_nondegenerate": float(df.rel_atp_span[~deg].median()),
            "q95_rel_atp_span": float(df.rel_atp_span.quantile(0.95)),
            "median_n_nondominated_of_40": float(df.n_nondominated.median()),
            "median_util_where_atp_floor_reached": float(df.util_at_atp_min[~deg].median()),
            "frac_seq_flux_changes_along_front": float((df.seq_change > 1e-9).mean()),
            "original_front_seed1": df[df.seed == 1].to_dict("records")}


# ------------------------------------------------------------ feedback (reanalysis)
def partFeedback():
    f = pd.read_csv(ORIG / "feedback.csv")
    out = {}
    for kind, g in f.groupby("kind"):
        conv = g[g.nonconverged_points == 0]
        e = {"n": len(g), "n_nonconverged_runs": int((g.nonconverged_points > 0).sum()),
             "frac_ge2_jumps_converged": float((conv.n_jumps >= 2).mean()),
             "frac_ge1_jump_converged": float((conv.n_jumps >= 1).mean()),
             "frac_ge1_jump_all": float((g.n_jumps >= 1).mean()),
             "frac_ge1_jump_among_nonconverged": float((g[g.nonconverged_points > 0].n_jumps >= 1).mean())
             if (g.nonconverged_points > 0).any() else None}
        if kind != "none":
            gb = np.log10(g.gamma * g.beta)
            e["median_log10_gamma_beta_jump"] = float(gb[g.n_jumps > 0].median()) if (g.n_jumps > 0).any() else None
            e["median_log10_gamma_beta_nojump"] = float(gb[g.n_jumps == 0].median())
        out[kind] = e
    return out


# ------------------------------------------------------------ failed instances
def partFailed():
    """the 6 original failures: re-solve lambda*_C with HiGHS presolve off."""
    from scipy.optimize import linprog
    s = json.loads((ORIG / "summary.json").read_text())
    rows = []
    for f in s["numerics"]["failed_detail"]:
        st = {x["name"]: x for x in RA.SETTINGS}[f["setting"]]
        inst = M.makeInstance(f["seed"], W=st["W"], wU=st["wU"], alpha=st["alpha"],
                              abundance=st["abundance"], structure=st["structure"], pools=P)
        U, cls, d, c = inst["U"], inst["cls"], inst["d"], inst["c"]
        R = len(U)
        Aeq = np.zeros((6, R + 1))
        Aeq[cls, np.arange(R)] = 1
        Aeq[:, R] = -d
        obj = np.zeros(R + 1)
        obj[R] = -1
        r = linprog(obj, A_ub=np.hstack([(U / c).T, np.zeros((15, 1))]), b_ub=np.ones(15),
                    A_eq=Aeq, b_eq=np.zeros(6), method="highs", options={"presolve": False})
        rows.append({"setting": f["setting"], "seed": f["seed"], "min_d": float(d.min()), "presolve_off_status": r.status,
                     "lamC_presolve_off": float(-r.fun) if r.status == 0 else None,
                     "n_struct_single_cuts": sum(len(x) == 1 for x in M.findStructuralCuts(inst, 1))})
    pd.DataFrame(rows).to_csv(OUT / "failed_instances_resolved.csv", index=False)
    return rows


def main(parts):
    OUT.mkdir(exist_ok=True)
    path = OUT / "audit_summary.json"
    summ = json.loads(path.read_text()) if path.exists() else {}
    summ["_label"] = "AUDIT 2026-10-09: new work, not the original run"
    with Pool(NPROC) as pool:
        if "repro" in parts:
            summ["repro"] = partRepro(pool)
        if "cf" in parts:
            summ["cf"] = partCF(pool)
        if "pareto" in parts:
            summ["pareto"] = partPareto(pool)
    if "struct" in parts:
        summ["struct"] = partStruct()
    if "feedback" in parts:
        summ["feedback"] = partFeedback()
    if "failed" in parts:
        summ["failed"] = partFailed()
    path.write_text(json.dumps(summ, indent=1, default=float))
    print(json.dumps(summ, indent=1, default=float)[:4000])


if __name__ == "__main__":
    main(sys.argv[1:] or ["struct", "feedback", "failed", "repro", "pareto", "cf"])
