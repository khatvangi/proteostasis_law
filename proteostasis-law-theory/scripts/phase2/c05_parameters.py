#!/usr/bin/env python
"""phase 2 step 5 -- which parameters distinguish multistable draws?

this is DIAGNOSTIC, not causal. the latin hypercube varies 16 parameters
independently, so a parameter that separates the classes is a marker of where in
the box multistability lives; it is not a demonstration that changing it causes
multistability, and no claim of that kind is made here.

no p-value fishing. sixteen parameters times several tests would manufacture
significance by volume, so the design avoids the temptation structurally:

  * effect sizes are reported directly (rank AUC and standardised mean
    difference) with BOOTSTRAP intervals, not p-values;
  * the single inferential question -- 'does ANY parameter combination carry
    signal at all?' -- is asked once, of a cross-validated classifier, against a
    label-permutation null. one question, one null, no per-parameter testing;
  * the ranking comes from held-out permutation importance, which is computed on
    data the model did not fit.

necessary-looking conditions are reported as what they are: one-sided thresholds
that happen to contain every observed positive, with the number of positives
supporting them, and an explicit statement that containment in a finite sample
is not necessity.

usage:  PHASE2_AUDIT_DIR=... python scripts/phase2/c05_parameters.py
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from _audit_common import (auditRoot, loadProvenance, loadSamples, runRoot,  # noqa: E402
                           writeJson, writeTable)

#: the 16 sampled parameters and the scale they were drawn on. features are
#: modelled on the SAMPLING scale, so a uniform-in-log parameter enters as its
#: log -- otherwise the effect size would mostly measure the skew of the prior.
FEATURES = {
    "nu": "log", "chi": "linear", "rho_U": "log", "rho_A": "log",
    "alpha_n": "log", "alpha_g": "log", "alpha_d": "log", "m": "linear",
    "kappa_ref": "log", "kappa_u": "log", "kappa_a": "log", "kappa_dis": "log",
    "kappa_cu": "log", "kappa_ca": "log", "kappa_du": "log", "kappa_da": "log",
}


def buildFeatures(df: pd.DataFrame) -> pd.DataFrame:
    X = pd.DataFrame(index=df.index)
    for name, scale in FEATURES.items():
        v = df[f"p_{name}"].astype(float)
        X[name if scale == "linear" else f"log10_{name}"] = (
            v if scale == "linear" else np.log10(v))
    return X


def rankAuc(x: np.ndarray, y: np.ndarray) -> float:
    """P(x_positive > x_negative) + 0.5 P(tie); the Mann-Whitney effect size."""
    from scipy.stats import rankdata
    pos, neg = int(y.sum()), int((~y.astype(bool)).sum())
    if pos == 0 or neg == 0:
        return float("nan")
    r = rankdata(x)
    return float((r[y.astype(bool)].sum() - pos * (pos + 1) / 2.0) / (pos * neg))


def bootstrapAucCi(x: np.ndarray, y: np.ndarray, n_boot: int = 2000,
                   seed: int = 20260731):
    """stratified bootstrap interval for the rank AUC."""
    rng = np.random.default_rng(seed)
    ip = np.flatnonzero(y.astype(bool))
    inn = np.flatnonzero(~y.astype(bool))
    if len(ip) < 2 or len(inn) < 2:
        return (float("nan"), float("nan"))
    vals = np.empty(n_boot)
    for b in range(n_boot):
        s = np.concatenate([rng.choice(ip, len(ip), replace=True),
                            rng.choice(inn, len(inn), replace=True)])
        vals[b] = rankAuc(x[s], y[s])
    return (float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5)))


def cohensD(x: np.ndarray, y: np.ndarray) -> float:
    a, b = x[y.astype(bool)], x[~y.astype(bool)]
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1))
                 / (len(a) + len(b) - 2))
    return float((a.mean() - b.mean()) / sp) if sp > 0 else float("nan")


def crossValidated(X: pd.DataFrame, y: np.ndarray, seed: int = 20260731,
                   n_perm: int = 100) -> dict:
    """repeated stratified CV for two model families, plus a permutation null."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import average_precision_score, roc_auc_score
    from sklearn.model_selection import StratifiedKFold, cross_val_predict
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    def models():
        return {
            "logistic_l2": make_pipeline(
                StandardScaler(),
                LogisticRegression(max_iter=5000, C=1.0, class_weight="balanced",
                                   random_state=seed)),
            "hist_gbdt": HistGradientBoostingClassifier(
                max_depth=3, max_iter=300, learning_rate=0.06,
                l2_regularization=1.0, random_state=seed,
                class_weight="balanced"),
        }

    Xv = X.to_numpy(float)
    out: dict = {"n": int(len(y)), "n_positive": int(y.sum())}
    for name, mdl in models().items():
        # cross_val_predict has no repeated-splitter form, so the repeats are
        # explicit; reporting the spread across repeats is what keeps a single
        # lucky split from being mistaken for a result
        aucs, aps = [], []
        for rep in range(5):
            skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed + rep)
            pr = cross_val_predict(mdl, Xv, y, cv=skf, method="predict_proba")[:, 1]
            aucs.append(roc_auc_score(y, pr))
            aps.append(average_precision_score(y, pr))
        out[name] = dict(roc_auc_mean=float(np.mean(aucs)),
                         roc_auc_sd=float(np.std(aucs, ddof=1)),
                         avg_precision_mean=float(np.mean(aps)),
                         avg_precision_sd=float(np.std(aps, ddof=1)),
                         base_rate=float(y.mean()))

    # one permutation null, for the one inferential question
    rng = np.random.default_rng(seed)
    null = []
    for b in range(n_perm):
        yp = rng.permutation(y)
        skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed + b)
        pr = cross_val_predict(models()["hist_gbdt"], Xv, yp, cv=skf,
                               method="predict_proba")[:, 1]
        null.append(roc_auc_score(yp, pr))
    obs = out["hist_gbdt"]["roc_auc_mean"]
    out["permutation_null"] = dict(
        n_permutations=n_perm, model="hist_gbdt",
        null_auc_mean=float(np.mean(null)), null_auc_sd=float(np.std(null, ddof=1)),
        null_auc_p95=float(np.percentile(null, 95)),
        observed_auc=float(obs),
        observed_exceeds_all_permutations=bool(obs > max(null)),
        empirical_p_ge=float((np.sum(np.array(null) >= obs) + 1) / (n_perm + 1)),
    )

    # held-out permutation importance, averaged over folds
    from sklearn.inspection import permutation_importance
    imp = {c: [] for c in X.columns}
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    for tr, te in skf.split(Xv, y):
        mdl = models()["hist_gbdt"]
        mdl.fit(Xv[tr], y[tr])
        r = permutation_importance(mdl, Xv[te], y[te], n_repeats=20,
                                   random_state=seed, scoring="roc_auc")
        for k, c in enumerate(X.columns):
            imp[c].append(float(r.importances_mean[k]))
    out["permutation_importance"] = {
        c: dict(mean=float(np.mean(v)), sd=float(np.std(v, ddof=1))) for c, v in imp.items()}
    return out


def necessaryLooking(X: pd.DataFrame, y: np.ndarray, ranges: dict) -> list:
    """one-sided thresholds satisfied by every observed positive.

    reported with the number of positives behind them and how much of the
    negative mass they exclude. containing all K positives in a sample is not
    necessity: with K positives, a threshold this tight arises by chance with
    probability roughly (excluded negative fraction)^K under exchangeability, so
    the count K is the honest measure of how much the condition is worth.
    """
    yb = y.astype(bool)
    out = []
    for c in X.columns:
        v = X[c].to_numpy(float)
        pos, neg = v[yb], v[~yb]
        if len(pos) == 0:
            continue
        for side, thr, keep_pos, excl in (
            ("lower", float(pos.min()), True, float(np.mean(neg < pos.min()))),
            ("upper", float(pos.max()), True, float(np.mean(neg > pos.max()))),
        ):
            if excl < 0.05:            # excludes almost nothing -> not a condition
                continue
            out.append(dict(
                feature=c, side=side, threshold=thr, n_positives=int(yb.sum()),
                negatives_excluded_fraction=excl,
                approx_chance_probability=float((1.0 - excl) ** int(yb.sum())),
                statement=(f"every observed positive has {c} "
                           f"{'>=' if side == 'lower' else '<='} {thr:.4g}"),
            ))
    out.sort(key=lambda d: -d["negatives_excluded_fraction"])
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default=None)
    ap.add_argument("--audit-dir", default=None)
    ap.add_argument("--seed", type=int, default=20260731)
    args = ap.parse_args()

    run = runRoot(args.run)
    out = auditRoot(args.audit_dir)
    df = loadSamples(run)

    fold = df[df["C1_fold_exists"].fillna(False).astype(bool)].copy()
    X = buildFeatures(fold)

    labels = {"phase1_second_stable":
              fold["C4_second_stable_attractor"].fillna(False).astype(bool).to_numpy(),
              "phase1_confirmed":
              fold["C4_second_attractor_confirmed"].fillna(False).astype(bool).to_numpy()}

    # audit-corrected label: a phase 1 positive counts only if the phase 2 audit
    # confirmed it. controls that the audit reclassified upward are added.
    tax = out / "c04_taxonomy.tsv"
    if tax.exists():
        t = pd.read_csv(tax, sep="\t")
        conf = set(t.loc[t["cls"] == "confirmed_multistable", "sample_index"])
        rejected = set(t.loc[(t["group"] == "candidate")
                             & (t["cls"] != "confirmed_multistable"), "sample_index"])
        lab = fold["C4_second_stable_attractor"].fillna(False).astype(bool).copy()
        lab.loc[fold["sample_index"].isin(rejected).to_numpy()] = False
        lab.loc[fold["sample_index"].isin(conf).to_numpy()] = True
        labels["audit_confirmed"] = lab.to_numpy()

    report: dict = dict(
        design=dict(
            n_fold_evaluable=int(len(fold)),
            features=list(X.columns),
            note=("descriptive; the latin hypercube varies parameters "
                  "independently, so these are markers of where multistability "
                  "sits in the sampled box, not causes of it"),
        ),
        labels={},
    )
    tables = []
    for lname, y in labels.items():
        if y.sum() < 5:
            report["labels"][lname] = dict(skipped=f"only {int(y.sum())} positives")
            continue
        eff = []
        for c in X.columns:
            v = X[c].to_numpy(float)
            auc = rankAuc(v, y)
            lo, hi = bootstrapAucCi(v, y, seed=args.seed)
            eff.append(dict(label=lname, feature=c, rank_auc=auc,
                            rank_auc_lo95=lo, rank_auc_hi95=hi,
                            auc_distance_from_half=abs(auc - 0.5),
                            cohens_d=cohensD(v, y),
                            median_positive=float(np.median(v[y.astype(bool)])),
                            median_negative=float(np.median(v[~y.astype(bool)])),
                            ci_excludes_half=bool(lo > 0.5 or hi < 0.5)))
        eff.sort(key=lambda d: -d["auc_distance_from_half"])
        tables.extend(eff)
        cvr = crossValidated(X, y.astype(int), seed=args.seed)
        rank = sorted(cvr["permutation_importance"].items(),
                      key=lambda kv: -kv[1]["mean"])
        report["labels"][lname] = dict(
            n_positive=int(y.sum()), n=int(len(y)), base_rate=float(y.mean()),
            effect_sizes=eff,
            cross_validated=cvr,
            importance_ranking=[dict(rank=i + 1, feature=k, **v)
                                for i, (k, v) in enumerate(rank)],
            necessary_looking=necessaryLooking(X, y, {}),
        )

    writeTable(pd.DataFrame(tables), out / "c05_effect_sizes.tsv")
    writeJson(report, out / "c05_parameters.json")

    for lname, r in report["labels"].items():
        if "skipped" in r:
            print(f"{lname}: {r['skipped']}")
            continue
        cv = r["cross_validated"]
        print(f"\n== {lname}  (n={r['n']}, positives={r['n_positive']}, "
              f"base rate={r['base_rate']:.4f})")
        print(f"   gbdt CV ROC-AUC {cv['hist_gbdt']['roc_auc_mean']:.3f} "
              f"+-{cv['hist_gbdt']['roc_auc_sd']:.3f}  |  AP "
              f"{cv['hist_gbdt']['avg_precision_mean']:.3f} (base {r['base_rate']:.4f})")
        print(f"   logistic CV ROC-AUC {cv['logistic_l2']['roc_auc_mean']:.3f}")
        pn = cv["permutation_null"]
        print(f"   permutation null AUC {pn['null_auc_mean']:.3f}"
              f" (p95 {pn['null_auc_p95']:.3f}); empirical p>= {pn['empirical_p_ge']:.4f}")
        print("   top features by held-out permutation importance:")
        for d in r["importance_ranking"][:6]:
            print(f"     {d['rank']}. {d['feature']:<16} {d['mean']:+.4f} +-{d['sd']:.4f}")
        print("   largest rank-AUC effects:")
        for d in r["effect_sizes"][:6]:
            print(f"     {d['feature']:<16} AUC {d['rank_auc']:.3f} "
                  f"[{d['rank_auc_lo95']:.3f},{d['rank_auc_hi95']:.3f}] d={d['cohens_d']:+.2f}")
    print(f"\n-> {out/'c05_parameters.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
