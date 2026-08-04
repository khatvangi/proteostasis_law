"""Extract the reduced arrays the figures plot, from the Phase 1 run root.

WHY THIS EXISTS
The run root is gitignored, so a clean checkout cannot regenerate Figures 2, 3
and S1. Everything else in this repository is reproducible from the deposit and
the figures must not be the hole. This writes ONLY the reduced arrays each figure
needs -- never the run root -- into `data/figures/`, which is tracked.

PROVENANCE IS PART OF THE DATA
Every file carries which run root it came from, how many states it holds, and
whether it is a subsample. Section 5 previously headed a table "2884 fold states"
while the rows beneath were computed on a different population of 325, and three
of its numbers came from an undeclared 20-state subsample that flattered all
three (D036). An array without its provenance is how that happens. The manifest
records a sha256 per file and the test suite asserts them.

Run with the run root present:

    python scripts/figures/build_figure_data.py
"""

from __future__ import annotations

import hashlib
import json
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

OUT = REPO_ROOT / "data" / "figures"
MANIFEST = OUT / "MANIFEST.json"

# the two populations, kept distinct and named at every use
LOAD_GRID = "load_grid"      # experiment B: nu x chi at fixed kinetics
KINETIC_BOX = "kinetic_box"  # experiment C: latin hypercube over kinetics


def _write(name: str, df: pd.DataFrame, population: str, source: str,
           complete: bool, note: str) -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / name
    df.to_csv(path, sep="\t", index=False)
    return {
        "file": name,
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "population": population,
        "n_states": int(len(df)),
        "is_subsample": not complete,
        "source": source,
        "columns": list(df.columns),
        "note": note,
    }


def buildSaturation(run: Path) -> dict:
    """Fig. 2: the three saturation fractions at collapse, over the kinetic box."""
    c = pd.read_csv(run / "C" / "samples.tsv", sep="\t")
    c = c[c["C1_fold_exists"] == True]  # noqa: E712
    c = c[pd.to_numeric(c["fold_burden"], errors="coerce").notna()]
    rows = []
    for _, r in c.iterrows():
        try:
            p = FT.paramsFromSampleRow(r)
            u, a = FT.foldStateFromSampleRow(r)
            d = FT.phiDecomposition(u, a, p)
        except (M.ModelError, ValueError, KeyError):
            continue
        rows.append({"s_ref": d["s_ref"], "s_u": d["s_u"], "s_a": d["s_a"],
                     "phi": d["phi"], "cf_frac": d["cf_frac"],
                     "df_frac": d["df_frac"]})
    df = pd.DataFrame(rows)
    return _write("saturation.tsv", df, KINETIC_BOX, run.name,
                  complete=len(df) == len(c),
                  note="michaelis factors at the fold; every draw that admits a "
                       "fold and rebuilds, no subsampling")


def buildIdentity(run: Path) -> dict:
    """Fig. S1: parallelism and BOTH residual normalisations, over the load grid."""
    b = pd.read_csv(run / "B" / "fold_boundary.tsv", sep="\t")
    b = b[b["found"] == True]  # noqa: E712
    base = M.Params()
    rows = []
    for _, r in b.iterrows():
        try:
            p = M.allocationParams(base.with_(nu=float(r["nu"])),
                                   float(r["chi"])).validate()
            u, a = float(r["fold_u"]), float(r["fold_a"])
            detJ = float(np.linalg.det(M.jacobian(u, a, p)))
            Ru, Ra = FT._centralGradient(FT.removalR, u, a, p)
            Gu, Ga = FT._centralGradient(FT.aggregateG, u, a, p)
            cross = Ru * Ga - Ra * Gu
            err = abs(detJ + cross)
            nR, nG = float(np.hypot(Ru, Ra)), float(np.hypot(Gu, Ga))
        except (M.ModelError, np.linalg.LinAlgError, OverflowError):
            continue
        rows.append({
            "nu": float(r["nu"]), "chi": float(r["chi"]),
            "eig_abs": abs(float(r["fold_eig_real_max"])),
            "sin_angle": abs(cross) / max(nR * nG, 1e-300),
            "res_gradient_normalised": err / max(nR * nG, 1e-300),
            "res_max_normalised": err / max(abs(detJ), abs(cross), 1e-300),
        })
    df = pd.DataFrame(rows)
    return _write("identity.tsv", df, LOAD_GRID, run.name,
                  complete=len(df) == len(b),
                  note="all found folds; carries BOTH normalisations so the "
                       "degradation of the max-normalised one is shown, not asserted")


def buildPareto() -> dict:
    """Fig. 3: the strategy front. Computed, not extracted -- no run root needed."""
    import pareto as P
    grid = P.strategyGrid()
    front = P.paretoFront(grid)
    df = grid[grid["feasible"]].copy()
    df["accuracy"] = -df["eps"]
    df["j_over_jcrit"] = df["j"] / df["j_crit"]
    key = set(zip(front["alpha"].round(12), front["R"].round(12)))
    df["on_front"] = [(a, r) in key for a, r in
                      zip(df["alpha"].round(12), df["R"].round(12))]
    df = df[["alpha", "R", "throughput", "eps", "accuracy", "j", "j_crit",
             "j_over_jcrit", "on_front"]]
    return _write("pareto.tsv", df, "computed", "scripts/phase3/pareto.py",
                  complete=True,
                  note="all feasible strategies with an on_front flag, so the "
                       "figure shows the cloud and the front; deterministic, "
                       "independent of the run root")


def main() -> int:
    run = FT.phase1RunDir()
    if not (run / "B" / "fold_boundary.tsv").is_file():
        print(f"run root absent at {run}; nothing rebuilt. "
              "data/figures/ is tracked, so the figures still reproduce.")
        return 0
    entries = [buildSaturation(run), buildIdentity(run)]
    try:
        entries.append(buildPareto())
    except Exception as exc:                                  # noqa: BLE001
        print(f"  pareto skipped: {exc}")
    manifest = {
        "run_root": run.name,
        "populations": {
            LOAD_GRID: "experiment B: nascent occupancy x rescue allocation at "
                       "fixed kinetics",
            KINETIC_BOX: "experiment C: latin hypercube over kinetic parameters",
        },
        "files": entries,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n")
    for e in entries:
        print("  %-24s %6d states  %-12s %s"
              % (e["file"], e["n_states"], e["population"], e["sha256"][:16]))
    print("manifest ->", MANIFEST)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
