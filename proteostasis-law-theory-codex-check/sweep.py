#!/usr/bin/env python3
"""Deterministic Latin-hypercube falsification sweep."""
import argparse, csv, json, os, platform, sys, time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from model import Params, equilibria, classify_equilibrium, integrate


def sample_params(cfg):
    names = list(cfg["ranges"])
    engine = qmc.LatinHypercube(d=len(names), seed=int(cfg["seed"]), scramble=True)
    unit = engine.random(int(cfg["samples"]))
    rows = []
    base = Params().dict()
    for x in unit:
        vals = base.copy()
        for i, name in enumerate(names):
            spec = cfg["ranges"][name]
            lo, hi = spec["bounds"]
            vals[name] = float(np.exp(np.log(lo)+x[i]*(np.log(hi)-np.log(lo)))) if spec.get("scale") == "log" else float(lo+x[i]*(hi-lo))
        rows.append(Params(**vals))
    return rows


def classify(p, cfg):
    try:
        roots = equilibria(p, cfg["search_max"], cfg["root_grid"])
        info = [(r, *classify_equilibrium(r, p)) for r in roots]
        stable = [item for item in info if item[1] == "stable"]
        hu, ha = cfg["thresholds"]["U"], cfg["thresholds"]["A"]
        below = [item for item in stable if item[0][0] < hu and item[0][1] < ha]
        if below:
            chosen, label = below[0], "stable_subthreshold"
        elif stable:
            chosen, label = stable[0], "stable_overthreshold"
        else:
            trials = []
            for y0 in cfg["initial_conditions"]:
                sol, ok, escaped = integrate(p, y0, cfg["t_final"], cfg["escape"])
                trials.append((ok, escaped, sol.y[:, -1]))
            if all(ok and escaped for ok, escaped, _ in trials):
                return {"class": "no_bounded_attracting_state", "n_roots": len(roots), "Ueq": "", "Aeq": "", "max_real_eig": ""}
            return {"class": "numerical_failure", "n_roots": len(roots), "Ueq": "", "Aeq": "", "max_real_eig": ""}
        r, _, eig = chosen
        return {"class": label, "n_roots": len(roots), "Ueq": float(r[0]), "Aeq": float(r[1]), "max_real_eig": float(np.max(eig.real))}
    except Exception:
        return {"class": "numerical_failure", "n_roots": "", "Ueq": "", "Aeq": "", "max_real_eig": ""}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    cfg = json.loads(Path(args.config).read_text())
    out = Path(args.output); out.mkdir(parents=True, exist_ok=True)
    started = time.time(); records = []
    for i, p in enumerate(sample_params(cfg)):
        records.append({"sample": i, **p.dict(), **classify(p, cfg)})
    with (out/"classifications.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=records[0]); w.writeheader(); w.writerows(records)
    counts = {k: sum(r["class"] == k for r in records) for k in ["stable_subthreshold", "stable_overthreshold", "no_bounded_attracting_state", "numerical_failure"]}
    summary = {"counts": counts, "samples": len(records), "seed": cfg["seed"], "elapsed_seconds": time.time()-started,
               "python": sys.version, "platform": platform.platform(), "pid": os.getpid(), "config": cfg}
    (out/"summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

