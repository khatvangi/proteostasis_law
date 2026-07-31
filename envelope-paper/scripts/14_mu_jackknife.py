#!/usr/bin/env python3
"""
is the mu clustering carried by a handful of codons, and is it biological at all?

two objections this addresses
----------------------------
1. LEVERAGE. mu_max is CCC at 2.03e-2, a rare codon. That single value sets the
   613-fold span on its own and gives Pro the largest Delta(mu) of any amino acid.
   The median-statistic sensitivity in Result 4 addresses skew WITHIN a codon's
   substitution set; it does nothing about a single codon dominating the whole
   axis. A leave-one-codon-out jackknife does.

2. DETECTABILITY. Landerer's per-codon mu is a mean over DETECTED substitutions,
   and detection probability varies with amino acid through ionization efficiency,
   whether the mass shift is resolvable, and peptide abundance. That manufactures
   amino-acid-level structure with no biological content, and it is a harder
   confound than the near-cognate-chemistry argument the paper already gives,
   because it predicts exactly the pattern we observe. Quantify the exposure:
   report the per-codon number of distinct detected substitutions, test whether mu
   is related to it, and check whether the clustering survives dropping the codons
   with the thinnest sampling.

   AND THEN ACTUALLY TEST IT. Dropping thin codons and reporting eta^2 for depth
   is exposure accounting, not a test: it cannot say whether the clustering
   survives removal of the depth-related component of mu. Regressing log mu on log
   depth and re-running the whole test on the residuals can. That is a reanalysis,
   not a new experiment, so conceding the point without running it concedes too
   early.

   the mediator objection to residualizing -- that a higher TRUE error rate
   generates more detected substitutions, so depth is a consequence of mu and
   residualizing on it would remove real signal -- is settled by the sign of the
   correlation. Spearman(depth, log mu) is NEGATIVE. Under error-driven detection
   it would be positive. A negative slope is the signature of thin-sampling
   inflation: a codon seen in few substitutions has its mean pulled up by whichever
   few were seen. Residualizing removes that, and takes no real signal with it
   unless the true relationship is negative, which no mechanism predicts.
"""
import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
COMP = ROOT / "data" / "computed"
HERE = Path(__file__).resolve().parent

_spec = importlib.util.spec_from_file_location("axis", HERE / "02_axis_structure.py")
axis_mod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(axis_mod)

N_PERM = 4000
SEED = 42


def test_mu(df, n_perm=N_PERM, seed=SEED):
    """within-degeneracy permutation test on the mu axis of a codon subset."""
    d = df.copy()
    v = d["log_mu"].to_numpy(float)
    d["_c"] = (v - v.mean()) / v.std(ddof=0)      # restandardize on the subset
    d["degeneracy"] = d.groupby("aa").codon.transform("size")
    d = d[d.degeneracy >= 2]
    if d.aa.nunique() < 5:
        return None
    xy = d[["_c"]].to_numpy()
    obs = axis_mod.mean_delta(d, xy)
    rng = np.random.default_rng(seed)
    vals = np.empty(n_perm)
    for k in range(n_perm):
        vals[k] = axis_mod.mean_delta(d, xy[axis_mod.permute(d, rng, "within_degeneracy")])
    z = (obs - vals.mean()) / vals.std(ddof=0)
    p = (int((vals <= obs).sum()) + 1) / (n_perm + 1)
    return {"n_codons": int(len(d)), "n_aa": int(d.aa.nunique()),
            "observed": float(obs), "null_mean": float(vals.mean()),
            "z": float(z), "p": float(p)}


def variance_partition(df, col):
    """eta^2 between amino acids, one-way, computed exactly as scripts/02 does."""
    groups = [g[col].to_numpy(float) for _, g in df.groupby("aa") if len(g) >= 2]
    allv = np.concatenate(groups)
    ss_b = sum(len(g) * (g.mean() - allv.mean()) ** 2 for g in groups)
    eta2 = float(ss_b / ((allv - allv.mean()) ** 2).sum())
    F, p = stats.f_oneway(*groups)
    return eta2, float(F), float(p)


def residualize_on_depth(df, n_perm=10_000, seed=SEED):
    """
    the test the exposure accounting could not do.

    regress log mu on log sampling depth, then re-run the whole analysis -- Delta,
    the within-degeneracy permutation test, and the variance partition -- on the
    residuals. if the amino-acid-level clustering survives, the detectability
    confound is not sufficient to explain it and the mu result stands on its own.
    if it collapses, the mu result is a property of the instrument and the paper
    has to say so.
    """
    d = df.dropna(subset=["n_subs"]).copy()
    x = np.log(d.n_subs.to_numpy(float))
    y = d.log_mu.to_numpy(float)
    fit = stats.linregress(x, y)
    d["resid"] = y - (fit.intercept + fit.slope * x)

    # the residual axis, standardized and tested exactly like log mu itself
    r = d.resid.to_numpy(float)
    d["_c"] = (r - r.mean()) / r.std(ddof=0)
    d["degeneracy"] = d.groupby("aa").codon.transform("size")
    d = d[d.degeneracy >= 2]
    xy = d[["_c"]].to_numpy()
    obs = axis_mod.mean_delta(d, xy)
    rng = np.random.default_rng(seed)
    vals = np.empty(n_perm)
    for k in range(n_perm):
        vals[k] = axis_mod.mean_delta(
            d, xy[axis_mod.permute(d, rng, "within_degeneracy")])
    z = float((obs - vals.mean()) / vals.std(ddof=0))
    p = (int((vals <= obs).sum()) + 1) / (n_perm + 1)

    eta2_resid, F_resid, p_resid = variance_partition(d, "resid")
    eta2_depth, _, _ = variance_partition(d, "n_subs")
    eta2_logmu, _, _ = variance_partition(d, "log_mu")
    return {
        "regression": {
            "model": "log mu ~ log(n_detected_substitutions)",
            "slope": float(fit.slope), "intercept": float(fit.intercept),
            "r": float(fit.rvalue), "r_squared": float(fit.rvalue ** 2),
            "p": float(fit.pvalue), "stderr": float(fit.stderr),
            "frac_log_mu_variance_explained_by_depth": float(fit.rvalue ** 2),
            "slope_sign_defuses_mediator_objection": bool(fit.slope < 0),
        },
        "residual_axis_test": {
            "n_codons": int(len(d)), "n_aa": int(d.aa.nunique()),
            "observed": float(obs), "null_mean": float(vals.mean()),
            "null_sd": float(vals.std(ddof=0)), "z": z, "p": float(p),
            "n_perm": int(n_perm),
            "significant_at_0.05": bool(p < 0.05),
            "pct_below_null": float(100.0 * (vals.mean() - obs) / vals.mean()),
        },
        "variance_partition_after_residualizing": {
            "eta2_between_aa_residual": eta2_resid,
            "F": F_resid, "p": p_resid,
            "eta2_between_aa_log_mu_same_subset": eta2_logmu,
            "eta2_between_aa_sampling_depth_same_subset": eta2_depth,
            "fraction_of_log_mu_eta2_retained": (eta2_resid / eta2_logmu
                                                 if eta2_logmu else None),
        },
    }


def main():
    df = axis_mod.load_axes("mean")

    # per-codon sampling depth, straight from the Landerer supplement
    ec = pd.read_excel(RAW / "Data_S2_error_detection_rate.xlsx", sheet_name="E. coli")
    depth = ec.rename(columns={"Codon": "codon",
                               "n_unique_aa_misincroporations": "n_subs"})
    df = df.merge(depth[["codon", "n_subs", "sd", "cv"]], on="codon", how="left")

    full = test_mu(df)
    print("=" * 76)
    print(f"full set: {full['n_codons']} codons, z = {full['z']:+.2f}, "
          f"p = {full['p']:.4f}")
    print("=" * 76)

    # ---- 1. leave-one-codon-out jackknife ----
    rows = []
    for c in df.codon:
        r = test_mu(df[df.codon != c])
        if r is None:
            continue
        rows.append({"dropped_codon": c, "aa": df.loc[df.codon == c, "aa"].iloc[0],
                     "mu": float(df.loc[df.codon == c, "mu"].iloc[0]),
                     **{k: r[k] for k in ("z", "p", "observed", "n_codons")}})
    jk = pd.DataFrame(rows).sort_values("z")
    jk.to_csv(COMP / "mu_jackknife.tsv", sep="\t", index=False)

    print("\n1. leave-one-codon-out jackknife")
    print(f"   z ranges {jk.z.min():+.2f} to {jk.z.max():+.2f}  "
          f"(full set {full['z']:+.2f})")
    print(f"   still significant at p < 0.05 after dropping any single codon: "
          f"{int((jk.p < 0.05).sum())}/{len(jk)}")
    print(f"   most influential (weakest z when dropped): "
          + ", ".join(f"{r.dropped_codon}/{r.aa} z={r.z:+.2f}"
                      for r in jk.tail(3).itertuples()))
    ccc = jk[jk.dropped_codon == "CCC"]
    if len(ccc):
        print(f"   dropping CCC (the 613-fold maximum): z = {ccc.iloc[0].z:+.2f}, "
              f"p = {ccc.iloc[0].p:.4f}")

    # span without the extreme codon
    no_ccc = df[df.codon != "CCC"]
    print(f"   mu span without CCC: {no_ccc.mu.max()/no_ccc.mu.min():.0f}-fold "
          f"(with CCC: {df.mu.max()/df.mu.min():.0f}-fold)")

    # ---- 2. detectability exposure ----
    print("\n2. detectability exposure")
    print(f"   per-codon detected substitutions: median {df.n_subs.median():.0f}, "
          f"range {df.n_subs.min():.0f}-{df.n_subs.max():.0f}")
    rho_mu, p_mu = stats.spearmanr(df.n_subs, df.log_mu)
    print(f"   Spearman(n_detected_substitutions, log mu) = {rho_mu:+.3f} "
          f"(p = {p_mu:.3g})")
    # does sampling depth itself carry amino-acid structure? if it does, the mu
    # clustering could be inherited from it
    groups = [g.n_subs.to_numpy() for _, g in df.groupby("aa") if len(g) >= 2]
    allv = np.concatenate(groups)
    ss_b = sum(len(g) * (g.mean() - allv.mean()) ** 2 for g in groups)
    eta2_depth = float(ss_b / ((allv - allv.mean()) ** 2).sum())
    F_d, p_d = stats.f_oneway(*groups)
    vd = json.loads((COMP / "mu_variance_decomposition.json").read_text())
    print(f"   eta^2 between amino acids: sampling depth {eta2_depth:.3f} "
          f"(F = {F_d:.2f}, p = {p_d:.3g})  vs  log mu "
          f"{vd['eta_squared_log_mu_between_aa']:.3f}")

    # clustering after dropping the most thinly sampled codons
    thin = []
    for q in (0, 10, 25):
        keep = df[df.n_subs >= np.percentile(df.n_subs, q)]
        r = test_mu(keep)
        thin.append({"drop_below_percentile": q, "n_codons": r["n_codons"],
                     "z": r["z"], "p": r["p"]})
        print(f"   dropping the thinnest {q:2d}% by detected substitutions "
              f"({r['n_codons']} codons): z = {r['z']:+.2f}, p = {r['p']:.4f}")
    pd.DataFrame(thin).to_csv(COMP / "mu_depth_sensitivity.tsv",
                              sep="\t", index=False)

    # ---- 3. the residualized test ----
    print("\n3. depth-residualized test (the one that decides it)")
    res = residualize_on_depth(df)
    rg, rt, vp = (res["regression"], res["residual_axis_test"],
                  res["variance_partition_after_residualizing"])
    print(f"   log mu ~ log depth: slope {rg['slope']:+.3f} "
          f"(p = {rg['p']:.3g}), R^2 = {rg['r_squared']:.3f}")
    print(f"   slope is {'negative -- thin-sampling inflation, not error-driven '
                        'detection' if rg['slope'] < 0 else 'POSITIVE -- mediator '
                        'objection applies, do not residualize'}")
    print(f"   residual axis: Delta = {rt['observed']:.4f} vs null "
          f"{rt['null_mean']:.4f}, z = {rt['z']:+.2f}, p = {rt['p']:.4f} "
          f"({rt['n_perm']:,} perms)")
    print(f"   eta^2 between amino acids: residual {vp['eta2_between_aa_residual']:.3f} "
          f"vs log mu {vp['eta2_between_aa_log_mu_same_subset']:.3f} "
          f"({100 * vp['fraction_of_log_mu_eta2_retained']:.0f}% retained)")
    print("   VERDICT: " + ("clustering SURVIVES residualization -- the confound "
                            "is not sufficient to explain it"
                            if rt["significant_at_0.05"] else
                            "clustering COLLAPSES -- the mu result is a property "
                            "of the instrument"))

    summary = {
        "full_set": full,
        "jackknife": {
            "z_min": float(jk.z.min()), "z_max": float(jk.z.max()),
            "n_significant": int((jk.p < 0.05).sum()), "n_total": int(len(jk)),
            "all_remain_significant": bool((jk.p < 0.05).all()),
            "most_influential_codons": jk.tail(3)[
                ["dropped_codon", "aa", "z", "p"]].to_dict("records"),
            "z_without_CCC": float(ccc.iloc[0].z) if len(ccc) else None,
            "p_without_CCC": float(ccc.iloc[0].p) if len(ccc) else None,
        },
        "mu_span_fold_with_CCC": float(df.mu.max() / df.mu.min()),
        "mu_span_fold_without_CCC": float(no_ccc.mu.max() / no_ccc.mu.min()),
        "detectability": {
            "n_detected_substitutions_median": float(df.n_subs.median()),
            "n_detected_substitutions_range": [float(df.n_subs.min()),
                                               float(df.n_subs.max())],
            "spearman_depth_vs_log_mu": float(rho_mu),
            "p_depth_vs_log_mu": float(p_mu),
            "eta2_between_aa_sampling_depth": eta2_depth,
            "F_depth": float(F_d), "p_depth": float(p_d),
            "eta2_between_aa_log_mu": vd["eta_squared_log_mu_between_aa"],
            "depth_sensitivity": thin,
            "residualized": res,
        },
        "caveat": (
            "the depth-residualized test is the substantive check: it removes the "
            "depth-related component of log mu and re-runs everything on the "
            "residuals. exposure accounting (eta^2 for depth, dropping thin "
            "codons) cannot settle the question and must not be quoted as if it "
            "had." if res["residual_axis_test"]["significant_at_0.05"] else
            "the clustering does not survive residualizing log mu on sampling "
            "depth. it must be reported as a property of the measured error "
            "landscape, not of decoding."),
    }
    (COMP / "mu_jackknife_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"\nwrote mu_jackknife.tsv, mu_depth_sensitivity.tsv and "
          f"mu_jackknife_summary.json to {COMP}")


if __name__ == "__main__":
    main()
