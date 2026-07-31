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
import os
from concurrent.futures import ProcessPoolExecutor
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
# 40 replicates put a 95% interval of roughly +-0.14 on a power estimate, wide
# enough that "does not reach 0.80" covered 0.80. 400 replicates halve that twice
# over; the reps are independent, so they run on the 64 cores available
N_POWER_REPS = 400
N_PERM_POWER = 800
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
    return obs, float(z), float(p), float(vals.mean()), float(vals.std(ddof=0))


def minimum_detectable(df, col, label):
    """largest s (weakest true clustering) still rejected at alpha = 0.05."""
    obs_delta, _, _, obs_null, _ = test_once(df, df[col].to_numpy(float),
                                             N_PERM_SWEEP, SEED)
    rows = []
    for s in np.arange(1.0, -0.001, -0.05):
        coord = shrink(df, col, float(s))
        d, z, p, nm, nsd = test_once(df, coord, N_PERM_SWEEP, SEED)
        rows.append({"axis": label, "s": round(float(s), 3), "delta": d,
                     "null_mean": nm, "null_sd": nsd, "z": z, "p": p,
                     "reject": bool(p < ALPHA),
                     "delta_reduction_pct": 100.0 * (obs_delta - d) / obs_delta,
                     "pct_below_null": 100.0 * (nm - d) / nm})
    sweep = pd.DataFrame(rows)
    rej = sweep[sweep.reject]
    mde = rej.iloc[0] if len(rej) else None
    obs_pct_below_null = 100.0 * (obs_null - obs_delta) / obs_null
    return sweep, obs_delta, mde, obs_pct_below_null


def _power_rep(args):
    """
    one power replicate. seeded from (SEED, s, r) rather than from a shared
    generator advanced in a loop, so the result does not depend on execution order
    and the replicates can run in parallel.
    """
    df, col, s, r = args
    rng = np.random.default_rng([SEED, int(s * 1000), r])
    coord = shrink(df, col, s, rng=rng, subset_prob=0.5)
    _, _, p, _, _ = test_once(df, coord, N_PERM_POWER, SEED + 1 + r)
    return int(p < ALPHA)


def wilson(k, n, z=1.96):
    """95% CI on a proportion. at 40 reps the interval was +-0.14; report it."""
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def power_at(df, col, s, reps=N_POWER_REPS):
    """rejection rate when the shrinkage hits a random half of the amino acids."""
    with ProcessPoolExecutor(max_workers=min(32, (os.cpu_count() or 4) - 2)) as ex:
        hits = list(ex.map(_power_rep,
                           [(df, col, s, r) for r in range(reps)],
                           chunksize=4))
    k = int(sum(hits))
    lo, hi = wilson(k, reps)
    return k / reps, lo, hi, reps


def mu_shrinkage_pattern(df, col="mu_z"):
    """
    mu's observed per-amino-acid tightness, as a shrinkage factor per amino acid.

    for amino acid A, s_A is the observed root-mean-square within-A deviation
    divided by the deviation expected if values were assigned at random,
    sd * sqrt(1 - 1/n_A). so s_A = 1 means A's synonyms are exactly as spread as
    chance predicts, and s_A < 1 means they are tighter. this is the quantity a
    uniform shrinkage sweep averages away: mu's structure is strongly non-uniform
    across amino acids, and transferring the average understates what a
    concentrated effect would look like on nu.
    """
    v = df[col].to_numpy(float)
    sd = v.std(ddof=0)
    out = {}
    for aa, idx in df.groupby("aa").indices.items():
        x = v[idx]
        n = len(x)
        if n < 2:
            continue
        obs = np.sqrt(((x - x.mean()) ** 2).mean())
        expected = sd * np.sqrt(1 - 1 / n)
        out[aa] = float(obs / expected) if expected > 0 else 1.0
    return out


def apply_pattern(df, col, pattern):
    """shrink (or expand) each amino acid's deviations by its own factor."""
    out = df[col].to_numpy(float).copy()
    for aa, idx in df.groupby("aa").indices.items():
        s = pattern.get(aa, 1.0)
        m = out[idx].mean()
        out[idx] = m + s * (out[idx] - m)
    return out


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
        sweep, obs_delta, mde, obs_below_null = minimum_detectable(df, col, label)
        sweeps.append(sweep)
        if mde is None:
            print(f"  {label:3s}: no shrinkage level rejected -- test has no power")
            out[label] = {"observed_delta": obs_delta, "mde_s": None,
                          "observed_pct_below_null": obs_below_null}
            continue
        # confirm the chosen point at full permutation count
        coord = shrink(df, col, float(mde.s))
        dc, zc, pc, _, _ = test_once(df, coord, N_PERM_CONFIRM, SEED)
        print(f"  {label:3s}: observed Delta = {obs_delta:.4f}; detects clustering "
              f"of >= {mde.delta_reduction_pct:.1f}% reduction in Delta "
              f"(s <= {mde.s:.2f}, z = {zc:+.2f}, p = {pc:.4f} at "
              f"{N_PERM_CONFIRM:,} perms)")
        out[label] = {
            "observed_delta": float(obs_delta),
            "observed_pct_below_null": float(obs_below_null),
            "mde_s": float(mde.s),
            "mde_pct_below_null": float(mde.pct_below_null),
            "mde_delta_reduction_pct": float(mde.delta_reduction_pct),
            "mde_z_confirmed": zc, "mde_p_confirmed": pc,
        }
    pd.concat(sweeps, ignore_index=True).to_csv(
        COMP / "nu_power_sweep.tsv", sep="\t", index=False)

    # THE comparison that answers the objection: is mu's real effect larger than
    # nu's detection floor? if not, the asymmetry in Result 4 is about resolution.
    if out["nu"].get("mde_pct_below_null") is not None:
        mu_real = out["mu"]["observed_pct_below_null"]
        nu_floor = out["nu"]["mde_pct_below_null"]
        print(f"\n  mu's observed clustering sits {mu_real:.1f}% below its null "
              f"mean.")
        print(f"  the nu test's detection floor is {nu_floor:.1f}% below null "
              f"(alpha = {ALPHA}).")
        verdict = ("detectable" if mu_real > nu_floor else "NOT detectable")
        print(f"  so an effect of mu's magnitude on the nu axis would have been "
              f"{verdict}\n  -- by a margin of "
              f"{abs(mu_real - nu_floor):.1f} percentage points.")
        out["comparison"] = {
            "mu_observed_pct_below_null": mu_real,
            "nu_detection_floor_pct_below_null": nu_floor,
            "mu_effect_exceeds_nu_floor": bool(mu_real > nu_floor),
            "margin_percentage_points": float(mu_real - nu_floor),
        }

    # ---- 3. power curve ----
    print("\n" + "=" * 74)
    print(f"3. power on the nu axis ({N_POWER_REPS} reps, shrinkage applied to a "
          f"random half of amino acids)")
    print("=" * 74)
    pw = []
    for s in (0.8, 0.6, 0.4, 0.2):
        power, lo, hi, reps = power_at(df, "nu_tai", s)
        coord = shrink(df, "nu_tai", s)
        d, _, _, _, _ = test_once(df, coord, 200, SEED)
        red = 100.0 * (out["nu"]["observed_delta"] - d) / out["nu"]["observed_delta"]
        pw.append({"s": s, "delta_reduction_pct_if_all_aa": red, "power": power,
                   "ci95_low": lo, "ci95_high": hi, "n_reps": reps,
                   "ci_excludes_0.80": bool(hi < 0.80 or lo > 0.80)})
        print(f"  s = {s:.1f}  (~{red:4.1f}% reduction in Delta)   "
              f"power = {power:.3f}  95% CI [{lo:.3f}, {hi:.3f}]"
              f"{'' if (hi < 0.80 or lo > 0.80) else '   <-- interval covers 0.80'}")
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
    covered = [r for r in pw if not r["ci95_excludes_0.80"]] if False else [
        r for r in pw if not r["ci_excludes_0.80"]]
    if covered:
        print(f"  NOTE: {len(covered)} grid point(s) have an interval covering "
              "0.80; the wording must not claim more than that")
    else:
        print("  every grid point's 95% interval excludes 0.80, so the claim is "
              "exact at this grid")

    # ---- 4. mu's own non-uniform pattern, transferred to nu ----
    print("\n" + "=" * 74)
    print("4. mu's observed per-amino-acid pattern imposed on nu")
    print("=" * 74)
    pattern = mu_shrinkage_pattern(df, "mu_z")
    s_vals = np.array(sorted(pattern.values()))
    print(f"  mu's per-amino-acid tightness s_A: min {s_vals.min():.2f}, "
          f"median {np.median(s_vals):.2f}, max {s_vals.max():.2f}; "
          f"{int((s_vals > 1).sum())}/{len(s_vals)} amino acids above 1 "
          "(spread wider than chance)")

    nu_patterned = apply_pattern(df, "nu_z", pattern)
    d_np, z_np, p_np, nm_np, _ = test_once(df, nu_patterned, N_PERM_CONFIRM, SEED)
    print(f"  nu with mu's pattern: Delta = {d_np:.4f} vs null {nm_np:.4f}, "
          f"z = {z_np:+.2f}, p = {p_np:.4f}  ->  "
          f"{'DETECTED' if p_np < ALPHA else 'NOT detected'}")

    # calibration: the same transfer onto a structureless axis must reproduce
    # something close to mu's own z, or the transfer is not faithful
    rng_cal = np.random.default_rng(SEED)
    cal = []
    for k in range(20):
        shuffled = df.mu_z.to_numpy(float)[rng_cal.permutation(len(df))]
        d2 = df.copy(); d2["_sh"] = shuffled
        _, z_cal, _, _, _ = test_once(d2, apply_pattern(d2, "_sh", pattern),
                                      N_PERM_SWEEP, SEED + k)
        cal.append(z_cal)
    _, z_mu_own, _, _, _ = test_once(df, df.mu_z.to_numpy(float),
                                     N_PERM_SWEEP, SEED)
    print(f"  calibration -- same pattern on a structureless axis: "
          f"z = {np.mean(cal):+.2f} +- {np.std(cal):.2f} "
          f"(mu's own z is {z_mu_own:+.2f}); the transfer is faithful if these "
          "are comparable")

    nu_pattern_result = {
        "what_this_answers": "the minimum-detectable-effect sweep transfers a "
                             "UNIFORM shrinkage, but mu's effect is strongly "
                             "non-uniform across amino acids. this imposes mu's "
                             "observed per-amino-acid pattern on nu instead.",
        "mu_pattern_s_by_aa": pattern,
        "mu_pattern_s_min": float(s_vals.min()),
        "mu_pattern_s_median": float(np.median(s_vals)),
        "mu_pattern_s_max": float(s_vals.max()),
        "n_aa_wider_than_chance": int((s_vals > 1).sum()),
        "n_aa": int(len(s_vals)),
        "nu_with_mu_pattern": {"delta": float(d_np), "null_mean": float(nm_np),
                               "z": float(z_np), "p": float(p_np),
                               "n_perm": N_PERM_CONFIRM,
                               "detected": bool(p_np < ALPHA)},
        "calibration_on_structureless_axis": {"mean_z": float(np.mean(cal)),
                                              "sd_z": float(np.std(cal)),
                                              "n": len(cal),
                                              "mu_own_z_for_comparison":
                                                  float(z_mu_own)},
    }

    summary = {
        "tie_structure": ties,
        "minimum_detectable_effect": out,
        "power_curve": pw,
        "n_power_reps": N_POWER_REPS,
        "n_perm_per_power_rep": N_PERM_POWER,
        "power_ci_method": "Wilson 95%",
        "nu_under_mu_observed_pattern": nu_pattern_result,
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
