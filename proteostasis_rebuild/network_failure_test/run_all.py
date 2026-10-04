"""run every ensemble, write results/. deterministic: every instance is a
pure function of (setting, seed).

survival criteria for the network hypothesis, fixed here BEFORE the full run
(the hypothesis: failure is distributed/cascading, not one pooled threshold)
  N1 pooling is a poor summary: median lambda*_C / lambda*_A < 0.5 in every
     bio setting with W >= 1 (abundance measured and blind)
  N2 redundancy-with-synergy exists: >= 25% of bio instances have >= 1
     synthetic-lethal machine pair (knockout, rho = 0.5) in C, robustly across
     W. (B has none by theorem; A has none unless two machines carry > half
     of all capacity)
  N3 multiple viability thresholds: >= 2 confirmed jumps in unserved load in
     >= 25% of converged feedback runs. 1 jump = a single cascade threshold.
     without feedback the LP admits exactly one threshold (theorem), so N3 can
     only come from feedback.
  any criterion that holds only in some settings is reported as conditional.
"""
import json
import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.dont_write_bytecode = True

import feedback as fbm
import model as M
from pools import findPools

OUT = Path(__file__).resolve().parent / "results"
RHOS = (0.25, 0.5, 0.9)
P = None

SETTINGS = (
    [dict(name=f"bio_measured_W{w}", W=w, abundance="measured", structure="bio", alpha=1.0, wU=1.0, n=1000)
     for w in (0, 1, 2, 3, 4)]
    + [dict(name=f"bio_blind_W{w}", W=w, abundance="blind", structure="bio", alpha=1.0, wU=1.0, n=1000)
       for w in (0, 2)]
    + [dict(name="bio_measured_W2_alpha0.3", W=2, abundance="measured", structure="bio", alpha=0.3, wU=1.0, n=500),
       dict(name="bio_measured_W2_alpha3", W=2, abundance="measured", structure="bio", alpha=3.0, wU=1.0, n=500),
       dict(name="bio_measured_W2_wU0", W=2, abundance="measured", structure="bio", alpha=1.0, wU=0.0, n=500),
       dict(name="unstructured_measured_W2", W=2, abundance="measured", structure="unstructured",
            alpha=1.0, wU=1.0, n=1000)]
)


def initWorker(pools):
    global P
    P = pools


def runInstance(args):
    st, seed = args
    inst = M.makeInstance(seed, W=st["W"], wU=st["wU"], alpha=st["alpha"],
                          abundance=st["abundance"], structure=st["structure"], pools=P)
    lamC, _x, y = M.solveLamC(inst, duals=True)
    row = {"setting": st["name"], "seed": seed, "lamA": M.solveLamA(inst),
           "lamB": M.solveLamB(inst), "lamC": lamC}
    row["C_over_A"] = row["lamC"] / row["lamA"]
    row["B_over_C"] = row["lamB"] / row["lamC"]
    yn = y / y.max() if y.max() > 0 else y
    row["bottlenecks"] = ";".join(M.MACHINES[i] for i in np.flatnonzero(yn > 1e-6))
    cuts = M.findStructuralCuts(inst, 3)
    for s in (1, 2, 3):
        row[f"C_struct_cuts_size{s}"] = sum(len(c) == s for c in cuts)
    row["B_struct_cuts_size1"] = len(M.findStructuralCutsB(inst))
    for kappa, tag in ((0.0, "ko"), (0.5, "kd50")):
        W1s = {}
        for arch in "ABC":
            _, W1, W2 = M.findKnockdowns(inst, arch, kappa)
            W1s[arch] = W1
            for rho in RHOS:
                row[f"{arch}_{tag}_single_lethal_rho{rho}"] = int(np.sum(W1 < rho))
                ps = M.findPairStats(W1, W2, rho)
                row[f"{arch}_{tag}_SL_pairs_rho{rho}"] = ps["synthetic_lethal_pairs"]
            ps = M.findPairStats(W1, W2, 0.5)
            row[f"{arch}_{tag}_pairs_beyond_weakest"] = ps["pairs_beyond_weakest_link"]
            row[f"{arch}_{tag}_synergistic_pairs"] = ps["synergistic_pairs"]
            row[f"{arch}_{tag}_min_epistasis"] = ps["min_epistasis"]
            row[f"{arch}_{tag}_singles"] = ";".join(f"{w:.6g}" for w in W1)
        if kappa == 0.5:
            rs = spearmanr(W1s["A"], W1s["C"]).statistic
            row["spearman_singles_A_vs_C_kd50"] = float(rs) if np.isfinite(rs) else np.nan
    ev = M.findSaturationEvents(inst, n=100)
    row["C_sat_events_before_failure"] = ev["n_events_before_failure"]
    row["C_machines_saturated_before_failure"] = ev["n_machines_saturated_before_failure"]
    row["C_first_event_frac"] = ev["first_event_frac"]
    row["C_sat_order"] = ";".join(f"{m}@{f:.2f}" for m, f in ev["order"])
    return row


def runFeedback(args):
    kind, seed = args
    rng = np.random.default_rng(10_000 + seed)
    inst = M.makeInstance(seed, W=2.0, pools=P)
    fb = fbm.makeFeedback(rng, kind) if kind != "none" else {"kind": "none", "gamma": 0.0, "beta": 0.0}
    r = fbm.findCascade(inst, fb)
    return {"kind": kind, "seed": seed, "gamma": fb["gamma"], "beta": fb["beta"],
            "theta": fb.get("theta", np.nan),
            "couplings": ";".join(f"{M.MACHINES[a]}->{M.MACHINES[b]}" for a, b in fb.get("couplings", [])),
            "onset_frac": r["onset_frac"], "n_jumps": r["n_jumps"],
            "jump_at": ";".join(f"{j['at_frac']:.4f}" for j in r["jumps"]),
            "jump_dZ": ";".join(f"{j['dZ']:.3f}" for j in r["jumps"]),
            "hysteresis": r["hysteresis"], "nonconverged_points": r["nonconverged_points"]}


def findParetoInstance():
    """first W=2 bio seed whose aggregate class has a sequestration route and a
    non-degenerate ATP front (illustrative only)."""
    for seed in range(1000):
        inst = M.makeInstance(seed, W=2.0, pools=P)
        if not any(k == "AGG" and "SEQ" in mods for k, mods in inst["routes"]):
            continue
        pts = M.findParetoFront(inst, 0.5, n=40)
        if pts[0]["atp"] > pts[-1]["atp"] * (1 + 1e-6):
            return seed, inst, pts
    return None, None, []


def main():
    OUT.mkdir(exist_ok=True)
    pools = findPools()
    nproc = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count()))
    jobs = [(st, s) for st in SETTINGS for s in range(st["n"])]
    with Pool(nproc, initializer=initWorker, initargs=(pools,)) as pool:
        rows = pool.map(runInstance, jobs, chunksize=4)
        fjobs = [("none", s) for s in range(100)] + [("titration", s) for s in range(600)] \
            + [("utilization", s) for s in range(600)]
        frows = pool.map(runFeedback, fjobs, chunksize=1)
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "instances.csv", index=False)
    fdf = pd.DataFrame(frows)
    fdf.to_csv(OUT / "feedback.csv", index=False)

    initWorker(pools)
    seed, inst, pts = findParetoInstance()
    pd.DataFrame(pts).to_csv(OUT / "pareto_front.csv", index=False)

    # per-machine knockout lethality (C, knockout, rho = 0.5) per setting
    mk = []
    for name, g in df.groupby("setting", sort=False):
        W = np.array([[float(v) for v in s.split(";")] for s in g["C_ko_singles"]])
        Wd = np.array([[float(v) for v in s.split(";")] for s in g["C_kd50_singles"]])
        for i, m in enumerate(M.MACHINES):
            mk.append({"setting": name, "machine": m,
                       "frac_ko_lethal_rho0.5": float(np.mean(W[:, i] < 0.5)),
                       "frac_ko_structural_zero": float(np.mean(W[:, i] == 0)),
                       "median_W_ko": float(np.median(W[:, i])),
                       "median_W_kd50": float(np.median(Wd[:, i])),
                       "frac_bottleneck": float(g["bottlenecks"].str.split(";").apply(lambda b: m in b).mean())})
    pd.DataFrame(mk).to_csv(OUT / "machine_knockouts.csv", index=False)

    summ = {"settings": {}, "feedback": {}, "pareto": {"seed": seed, "load_frac": 0.5,
                                                         "n_points": len(pts)}}
    q = lambda s: {"median": float(s.median()), "q05": float(s.quantile(0.05)),
                   "q95": float(s.quantile(0.95))}
    for name, g in df.groupby("setting", sort=False):
        e = {"n": len(g), "C_over_A": q(g.C_over_A), "B_over_C": q(g.B_over_C),
             "frac_C_struct_single_cut": float((g.C_struct_cuts_size1 > 0).mean()),
             "median_C_struct_cuts_size1": float(g.C_struct_cuts_size1.median()),
             "median_B_struct_cuts_size1": float(g.B_struct_cuts_size1.median()),
             "spearman_singles_A_vs_C_kd50": q(g.spearman_singles_A_vs_C_kd50.dropna()),
             "C_sat_events_before_failure": {str(k): int(v) for k, v in
                                             g.C_sat_events_before_failure.value_counts().sort_index().items()},
             "C_first_event_frac": q(g.C_first_event_frac.dropna())}
        for arch in "ABC":
            for tag in ("ko", "kd50"):
                for rho in RHOS:
                    e[f"{arch}_{tag}_frac_inst_with_SL_pair_rho{rho}"] = float((g[f"{arch}_{tag}_SL_pairs_rho{rho}"] > 0).mean())
                    e[f"{arch}_{tag}_mean_single_lethal_rho{rho}"] = float(g[f"{arch}_{tag}_single_lethal_rho{rho}"].mean())
                e[f"{arch}_{tag}_frac_inst_pairs_beyond_weakest"] = float((g[f"{arch}_{tag}_pairs_beyond_weakest"] > 0).mean())
                e[f"{arch}_{tag}_median_min_epistasis"] = float(g[f"{arch}_{tag}_min_epistasis"].median())
        bn = g["bottlenecks"].str.split(";").explode().value_counts(normalize=False) / len(g)
        e["bottleneck_frequency"] = {k: round(float(v), 4) for k, v in bn.items()}
        summ["settings"][name] = e
    for kind, g in fdf.groupby("kind", sort=False):
        conv = g[g.nonconverged_points == 0]
        summ["feedback"][kind] = {
            "n": len(g), "n_fully_converged": len(conv),
            "frac_ge1_jump_converged": float((conv.n_jumps >= 1).mean()) if len(conv) else None,
            "frac_ge2_jumps_converged": float((conv.n_jumps >= 2).mean()) if len(conv) else None,
            "n_jumps_counts_all": {str(k): int(v) for k, v in g.n_jumps.value_counts().sort_index().items()},
            "frac_onset_before_lamC": float((g.onset_frac < 1.0).mean()),
            "onset_frac": q(g.onset_frac.dropna()),
            "frac_hysteresis_gt_0.02_converged": float((conv.hysteresis > 0.02).mean()) if len(conv) else None}
    (OUT / "summary.json").write_text(json.dumps(summ, indent=1))
    print(json.dumps(summ, indent=1)[:6000])


if __name__ == "__main__":
    main()
