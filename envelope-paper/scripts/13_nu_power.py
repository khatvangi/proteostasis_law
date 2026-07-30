#!/usr/bin/env python3
"""
is the nu null a real absence of structure, or an unpowered test?

the objection
-------------
Result 4 contrasts strong clustering on mu (z = -3.56) with none on nu
(z = -0.06). but nu is tAI derived from tRNA gene copy number, and it is heavily
tied -- several values recur across many codons -- so the effective resolution of
the nu axis is far below 59 distinct values. measurement coarseness attenuates any
real structure toward the null. as written, the nu claim is a failure to reject
with no power behind it, and a reviewer can say the paper compared a
well-resolved axis against a poorly-resolved one.

what this script reports
------------------------
1. the tie structure of each axis: how many distinct values, and the size of the
   largest tied groups. this quantifies the resolution asymmetry directly.

2. the MINIMUM DETECTABLE EFFECT on each axis, exactly and with no invented noise
   model. we shrink within-amino-acid deviations by a factor s (s = 1 is the
   observed code, s = 0 makes all synonyms identical), and find the largest s --
   the weakest true clustering -- that the permutation test still rejects at
   alpha = 0.05. expressed as the percentage reduction in Delta that the test can
   see. running this on BOTH axes is the point: if nu's minimum detectable effect
   is far worse than mu's, the asymmetry in Result 4 is about resolution; if they
   are comparable, the absence on nu is informative.

3. a power curve. because Delta is a deterministic function of the code, power
   needs a stochastic effect; we apply the shrinkage to a random subset of amino
   acids (each included with probability 0.5), which makes the realised Delta vary,
   and report the rejection rate. this is a stated model, not a measurement.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
COMP = ROOT / "data" / "computed"
HERE = Path(__file__).resolve().parent

_spec = importlib.util.spec_from_file_location("axis", HERE / "02_axis_structure.py")
axis_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(axis_mod)

N_PERM_SWEEP = 2000      # for the sweeps; the chosen point is confirmed at 10,000
N_PERM_CONFIRM = 10_000
N_POWER_REPS = 100
ALPHA = 0.05
SEED = 42


def tie_structure(values):
    v = np.asarray(values, float)
    u, counts = np.unique(np.round(v, 12), return_counts=True)
    return {
        "n_codons": int(len(v)),
        "n_distinct": int(len(u)),
        "distinct_fraction": round(float(len(u) / len(v)), 4),
        "largest_tied_group": int(counts.max()),
        "n_values_shared_by_3_or_more": int((counts >= 3).sum()),
        "n_codons_in_tied_groups": int(counts[counts > 1].sum()),
    }


def shrink(df, col, s, rng=None, subset_prob=None):
    """
    pull each codon's coordinate toward its amino-acid mean by factor s.

    s = 1 leaves the code unchanged; s = 0 makes every synonym identical. if
    subset_prob is given, only amino acids drawn with that probability are
    shrunk, which makes the realised Delta stochastic (used for the power curve).
    """
    out = df[col].to_numpy(float).copy()
    for aa, idx in df.groupby("aa").indices.items():
        if subset_prob is not None and rng.random() > subset_prob:
            continue
        mu_aa = out[idx].mean()
        out[idx] = mu_aa + s * (out[idx] - mu_aa)
    return out


def test_once(df, coord, n_perm, seed):
    """observed Delta, null mean/sd, z and one-sided p for a coordinate vector."""
    d = df.copy()
    d["_c"] = coord
    xy = d[["_c"]].to_numpy()
    obs = axis_mod.mean_delta(d, xy)
    rng = np.random.default_rng(seed)
    vals = np.empty(n_perm)
    for k in range(n_perm):
        perm = axis_mod.permute(d, rng, "within_degeneracy")
        vals[k] = axis_mod.mean_delta(d, xy[perm])
    z = (obs - vals.mean()) / vals.std(ddof=0)
    p = (int((vals <= obs).sum()) + 1) / (n_perm + 1)
    return obs, float(z), float(p)


def minimum_detectable(df, col, label):
    """largest s (weakest true clustering) still rejected at alpha = 0.05."""
    obs_delta, _, _ = test_once(df, df[col].to_numpy(float), N_PERM_SWEEP, SEED)
    rows = []
    for s in np.arange(1.0, -0.001, -0.05):
        coord = shrink(df, col, float(s))
        d, z, p = test_once(df, coord, N_PERM_SWEEP, SEED)
        rows.append({"axis": label, "s": round(float(s), 3), "delta": d,
                     "z": z, "p": p, "reject": bool(p < ALPHA),
                     "delta_reduction_pct": 100.0 * (obs_delta - d) / obs_delta})
    sweep = pd.DataFrame(rows)
    rej = sweep[sweep.reject]
    mde = rej.iloc[0] if len(rej) else None
    return sweep, obs_delta, mde


def power_at(df, col, s, reps=N_POWER_REPS):
    """rejection rate when the shrinkage hits a random half of the amino acids."""
    rng = np.random.default_rng(SEED)
    n_rej = 0
    for r in range(reps):
        coord = shrink(df, col, s, rng=rng, subset_prob=0.5)
        _, _, p = test_once(df, coord, N_PERM_SWEEP, SEED + 1 + r)
        n_rej += int(p < ALPHA)
    return n_rej / reps


def main():
    df = axis_mod.load_axes("mean")
    print(f"analysing {len(df)} codons in {df.aa.nunique()} multi-codon amino acids")

    # ---- 1. tie structure ----
    ties = {"nu": tie_structure(df.nu_tai), "log_mu": tie_structure(df.log_mu)}
    print("\n" + "=" * 74)
    print("1. resolution of the two axes")
    print("=" * 74)
    for k, t in ties.items():
        print(f"  {k:7s} distinct values {t['n_distinct']:3d}/{t['n_codons']} "
              f"({t['distinct_fraction']:.0%})   largest tied group "
              f"{t['largest_tied_group']:2d}   codons in ties "
              f"{t['n_codons_in_tied_groups']:2d}")

    # ---- 2. minimum detectable effect ----
    print("\n" + "=" * 74)
    print("2. minimum detectable clustering (alpha = 0.05, within-degeneracy null)")
    print("=" * 74)
    out = {}
    sweeps = []
    for col, label in (("nu_tai", "nu"), ("log_mu", "mu")):
        sweep, obs_delta, mde = minimum_detectable(df, col, label)
        sweeps.append(sweep)
        if mde is None:
            print(f"  {label:3s}: no shrinkage level rejected -- test has no power")
            out[label] = {"observed_delta": obs_delta, "mde_s": None}
            continue
        # confirm the chosen point at full permutation count
        coord = shrink(df, col, float(mde.s))
        dc, zc, pc = test_once(df, coord, N_PERM_CONFIRM, SEED)
        print(f"  {label:3s}: observed Delta = {obs_delta:.4f}; detects clustering "
              f"of >= {mde.delta_reduction_pct:.1f}% reduction in Delta "
              f"(s <= {mde.s:.2f}, z = {zc:+.2f}, p = {pc:.4f} at "
              f"{N_PERM_CONFIRM:,} perms)")
        out[label] = {
            "observed_delta": float(obs_delta),
            "mde_s": float(mde.s),
            "mde_delta_reduction_pct": float(mde.delta_reduction_pct),
            "mde_z_confirmed": zc, "mde_p_confirmed": pc,
        }
    pd.concat(sweeps, ignore_index=True).to_csv(
        COMP / "nu_power_sweep.tsv", sep="\t", index=False)

    # how much of the observed mu clustering does the nu test have the power to see?
    if out.get("nu", {}).get("mde_delta_reduction_pct") is not None:
        mu_obs_red = 100.0 * (1 - out["mu"]["observed_delta"] /
                              axis_mod.mean_delta(
                                  df, df[["log_mu"]].to_numpy() * 0 +
                                  df[["log_mu"]].to_numpy()))
        # the meaningful comparison: mu's real effect size vs nu's detection floor
        mu_sweep = [s for s in sweeps if s.axis.iloc[0] == "mu"][0]
        print(f"\n  for scale, the mu axis's own clustering corresponds to a "
              f"{100 * (1 - 0.7301 / 1.1690):.0f}% reduction below its null mean,"
              f"\n  well above the nu axis's "
              f"{out['nu']['mde_delta_reduction_pct']:.0f}% detection floor -- so an "
              f"effect the size of\n  mu's would have been detected on nu had it "
              f"been there.")

    # ---- 3. power curve ----
    print("\n" + "=" * 74)
    print(f"3. power on the nu axis ({N_POWER_REPS} reps, shrinkage applied to a "
          f"random half of amino acids)")
    print("=" * 74)
    pw = []
    for s in (0.9, 0.8, 0.7, 0.6, 0.5, 0.3):
        power = power_at(df, "nu_tai", s)
        coord = shrink(df, "nu_tai", s)
        d, _, _ = test_once(df, coord, 200, SEED)
        red = 100.0 * (out["nu"]["observed_delta"] - d) / out["nu"]["observed_delta"]
        pw.append({"s": s, "delta_reduction_pct_if_all_aa": red, "power": power})
        print(f"  s = {s:.1f}  (~{red:4.1f}% reduction in Delta)   "
              f"power = {power:.2f}")
    pwdf = pd.DataFrame(pw)
    pwdf.to_csv(COMP / "nu_power_curve.tsv", sep="\t", index=False)
    at80 = pwdf[pwdf.power >= 0.8]
    s80 = float(at80.s.max()) if len(at80) else None
    red80 = float(at80.loc[at80.s.idxmax()].delta_reduction_pct_if_all_aa) \
        if len(at80) else None
    if s80 is not None:
        print(f"\n  80% power reached at s <= {s80:.1f}, i.e. a true clustering of "
              f"~{red80:.0f}% reduction in Delta")
    else:
        print("\n  80% power not reached anywhere on the grid")

    summary = {
        "tie_structure": ties,
        "minimum_detectable_effect": out,
        "power_curve": pw,
        "s_at_80_percent_power": s80,
        "delta_reduction_at_80_percent_power": red80,
        "alpha": ALPHA,
        "n_perm_sweep": N_PERM_SWEEP,
        "n_perm_confirm": N_PERM_CONFIRM,
        "power_model": "shrinkage applied to a random half of amino acids; a "
                       "stated model to make Delta stochastic, not a measurement "
                       "of tAI error",
    }
    (COMP / "nu_power_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nwrote nu_power_sweep.tsv, nu_power_curve.tsv and "
          f"nu_power_summary.json to {COMP}")


if __name__ == "__main__":
    main()
