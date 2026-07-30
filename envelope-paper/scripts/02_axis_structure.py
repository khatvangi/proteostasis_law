#!/usr/bin/env python3
"""
operational structure of the genetic code in (mu, nu) space.

mu = per-codon mistranslation rate (E. coli, Landerer et al. 2024, Data S2)
nu = translational supply (tAI), validated by scripts/01_validate_tai.py

for each multi-codon amino acid A we compute the operational spread

    Delta_A = mean pairwise Euclidean distance among its synonyms
              in standardized coordinates

and compare the mean over amino acids against permutation nulls.

what this script deliberately adds beyond the earlier version
-------------------------------------------------------------
1. it regenerates the observed deltas and nulls from RAW inputs. the previous
   version of this project shipped these as precomputed TSVs with no generating
   script, which is how a corrupt nu vector survived undetected.

2. `--mu-stat median` sensitivity. the Landerer per-codon distributions are
   heavily right-skewed (median mean/median ratio ~4.9), so the choice of
   summary statistic is a real analytic decision and its effect is reported
   rather than buried.

3. a variance decomposition of log mu into between- and within-amino-acid
   components. this matters for interpretation: the permutation null asks
   whether reassigning mu ACROSS amino acids would loosen within-amino-acid
   spread, and the answer is almost forced to be yes if mu carries any
   amino-acid-level structure at all. the decomposition reports how much
   structure there is, so the manuscript can claim "mu is organized at the
   amino-acid level" (supported) rather than "the code was optimized to
   minimize error cost" (not established by this test).
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "computed"

BASES = "TCAG"
_AA = ("FFLLSSSSYY**CC*W" "LLLLPPPPHHQQRRRR"
       "IIIMTTTTNNKKSSRR" "VVVVAAAADDEEGGGG")
CODON_AA = {
    b1 + b2 + b3: _AA[i]
    for i, (b1, b2, b3) in enumerate(
        (a, b, c) for a in BASES for b in BASES for c in BASES)
}

N_PERM = 10_000
SEED = 42


def load_axes(mu_stat="mean"):
    """
    return a dataframe of codons carrying both axes.

    mu is read from the Landerer supplement so that the choice of per-codon
    summary statistic (mean vs median) is explicit and switchable.
    """
    ec = pd.read_excel(RAW / "Data_S2_error_detection_rate.xlsx", sheet_name="E. coli")
    mu = ec[["Codon", mu_stat]].rename(columns={"Codon": "codon", mu_stat: "mu"})

    nu = pd.read_csv(OUT / "nu_tai_ecoli_validated.tsv", sep="\t")

    df = mu.merge(nu, on="codon")
    df["aa"] = df.codon.map(CODON_AA)
    df = df[(df.aa != "*") & (df.mu > 0)].copy()

    # log mu, because raw mu spans ~3 orders of magnitude and is right-skewed
    df["log_mu"] = np.log(df.mu)
    # standardize both axes across the analysed codon set
    for src, dst in (("log_mu", "mu_z"), ("nu_tai", "nu_z")):
        v = df[src].to_numpy(float)
        df[dst] = (v - v.mean()) / v.std(ddof=0)

    # keep only amino acids that actually have synonyms present
    df["degeneracy"] = df.groupby("aa").codon.transform("size")
    return df[df.degeneracy >= 2].reset_index(drop=True)


def coords(df, axis):
    if axis == "mu":
        return df[["mu_z"]].to_numpy()
    if axis == "nu":
        return df[["nu_z"]].to_numpy()
    return df[["mu_z", "nu_z"]].to_numpy()


def delta_per_aa(df, xy):
    """mean pairwise distance among synonyms, per amino acid."""
    out = {}
    for aa, idx in df.groupby("aa").indices.items():
        if len(idx) < 2:
            continue
        pts = xy[idx]
        d = [np.linalg.norm(pts[i] - pts[j])
             for i in range(len(pts)) for j in range(i + 1, len(pts))]
        out[aa] = float(np.mean(d))
    return out


def mean_delta(df, xy):
    return float(np.mean(list(delta_per_aa(df, xy).values())))


def permute(df, rng, scheme):
    """
    return a permutation of row indices, i.e. a reassignment of (mu, nu) values
    to codons.

    within_degeneracy : values may only move between codons whose amino acids
                        have the same degeneracy. preserves how many codons each
                        amino acid has and the value pool within each class.
    full              : values may move anywhere among the analysed codons.
                        preserves only the alphabet size.
    """
    n = len(df)
    perm = np.empty(n, dtype=int)
    if scheme == "full":
        perm[:] = rng.permutation(n)
        return perm
    for _, idx in df.groupby("degeneracy").indices.items():
        perm[idx] = rng.permutation(idx)
    return perm


def null_distribution(df, axis, scheme, n_perm=N_PERM, seed=SEED):
    xy = coords(df, axis)
    rng = np.random.default_rng(seed)
    vals = np.empty(n_perm)
    for k in range(n_perm):
        p = permute(df, rng, scheme)
        vals[k] = mean_delta(df, xy[p])
    return vals


def variance_decomposition(df):
    """
    how much of the variation in log mu sits between amino acids rather than
    within them. reported as eta^2 plus a one-way F test.

    this is the descriptive backbone of the mu claim: it quantifies
    amino-acid-level organization without attributing it to a cause.
    """
    from scipy import stats as st
    groups = [g.log_mu.to_numpy() for _, g in df.groupby("aa") if len(g) >= 2]
    allv = np.concatenate(groups)
    grand = allv.mean()
    ss_between = sum(len(g) * (g.mean() - grand) ** 2 for g in groups)
    ss_total = ((allv - grand) ** 2).sum()
    F, p = st.f_oneway(*groups)
    return {
        "eta_squared_log_mu_between_aa": round(float(ss_between / ss_total), 4),
        "F": round(float(F), 3),
        "p": float(f"{p:.3g}"),
        "n_amino_acids": len(groups),
        "n_codons": int(len(allv)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mu-stat", default="mean", choices=["mean", "median"])
    ap.add_argument("--n-perm", type=int, default=N_PERM)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    df = load_axes(args.mu_stat)
    print(f"analysing {len(df)} codons across "
          f"{df.aa.nunique()} multi-codon amino acids (mu = per-codon {args.mu_stat})")

    rows, nulls = [], {}
    for axis in ("mu", "nu", "2D"):
        xy = coords(df, axis)
        obs = mean_delta(df, xy)
        for scheme in ("within_degeneracy", "full"):
            v = null_distribution(df, axis, scheme, args.n_perm)
            z = (obs - v.mean()) / v.std(ddof=0)
            # one-sided empirical p in the observed direction, with the
            # standard +1 correction so p can never be reported as exactly 0
            if z < 0:
                p = (int((v <= obs).sum()) + 1) / (args.n_perm + 1)
            else:
                p = (int((v >= obs).sum()) + 1) / (args.n_perm + 1)
            rows.append({
                "axis": axis, "null": scheme, "mu_stat": args.mu_stat,
                "observed": obs, "null_mean": float(v.mean()),
                "null_sd": float(v.std(ddof=0)), "z": float(z),
                "p_one_sided": float(p),
                "direction": "clustered" if z < 0 else "spread",
            })
            if scheme == "within_degeneracy":
                nulls[axis] = v
            print(f"  {axis:>3} / {scheme:<17} observed={obs:.4f} "
                  f"null={v.mean():.4f}  z={z:+.2f}  p={p:.4f}")

    suffix = "" if args.mu_stat == "mean" else f"_{args.mu_stat}"
    pd.DataFrame(rows).to_csv(OUT / f"axis_tests{suffix}.tsv", sep="\t", index=False)

    for axis, v in nulls.items():
        pd.DataFrame({"delta": v}).to_csv(
            OUT / f"null_{axis}{suffix}.tsv", sep="\t", index=False)

    xy2 = coords(df, "2D")
    per_aa = pd.DataFrame(
        [{"aa": a, "degeneracy": int(df[df.aa == a].degeneracy.iloc[0]),
          "delta_mu": delta_per_aa(df, coords(df, "mu"))[a],
          "delta_nu": delta_per_aa(df, coords(df, "nu"))[a],
          "delta_2D": delta_per_aa(df, xy2)[a]}
         for a in sorted(df.aa.unique())])
    per_aa.to_csv(OUT / f"delta_per_aa{suffix}.tsv", sep="\t", index=False)

    df.to_csv(OUT / f"codon_axes{suffix}.tsv", sep="\t", index=False)

    if args.mu_stat == "mean":
        vd = variance_decomposition(df)
        (OUT / "mu_variance_decomposition.json").write_text(json.dumps(vd, indent=2))
        print(f"\nvariance decomposition of log mu: "
              f"eta^2 between amino acids = {vd['eta_squared_log_mu_between_aa']:.3f} "
              f"(F = {vd['F']}, p = {vd['p']})")

        span = df.mu.max() / df.mu.min()
        print(f"mu span across codons: {span:.0f}-fold "
              f"({df.mu.min():.2e} to {df.mu.max():.2e})")

    print(f"\nwrote analysis outputs to {OUT}")


if __name__ == "__main__":
    main()
