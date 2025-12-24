#!/usr/bin/env python

import argparse
import csv
import math
from collections import defaultdict

import numpy as np
import statsmodels.api as sm

TWO_CODON_AA = {
    "C": ["TGT","TGC"],
    "D": ["GAT","GAC"],
    "E": ["GAA","GAG"],
    "F": ["TTT","TTC"],
    "H": ["CAT","CAC"],
    "K": ["AAA","AAG"],
    "N": ["AAT","AAC"],
    "Q": ["CAA","CAG"],
    "Y": ["TAT","TAC"],
}

def main():
    ap = argparse.ArgumentParser(
        description="Fit logistic models: P(metal | μ, tAI) per 2-codon AA"
    )
    ap.add_argument("--residue_tsv", required=True,
                    help="errors/residue_kappa_table.tsv")
    args = ap.parse_args()

    # Collect data per AA
    data_by_aa = defaultdict(list)

    with open(args.residue_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            aa = row["aa"]
            if aa not in TWO_CODON_AA:
                continue
            codon = row["codon"]
            if codon not in TWO_CODON_AA[aa]:
                continue
            is_metal = int(row["is_metal"])
            try:
                mu = float(row["mu"])
                tai = float(row["w_tai"])
            except ValueError:
                continue
            if math.isnan(mu) or math.isnan(tai):
                continue
            data_by_aa[aa].append((mu, tai, is_metal))

    for aa, records in sorted(data_by_aa.items()):
        if len(records) < 50:
            print(f"[WARN] {aa}: too few residues ({len(records)})")
            continue
        print(f"\n=== AA {aa} (n={len(records)}) ===")

        X_mu = np.array([[r[0]] for r in records])       # μ
        X_tai = np.array([[r[1]] for r in records])      # tAI
        X_both = np.array([[r[0], r[1]] for r in records])
        y = np.array([r[2] for r in records])

        # Add intercept
        X_mu_ = sm.add_constant(X_mu)
        X_tai_ = sm.add_constant(X_tai)
        X_both_ = sm.add_constant(X_both)

        # μ-only
        model_mu = sm.Logit(y, X_mu_).fit(disp=0)
        # tAI-only
        model_tai = sm.Logit(y, X_tai_).fit(disp=0)
        # μ + tAI
        model_both = sm.Logit(y, X_both_).fit(disp=0)

        print("  μ-only:  coef(mu) = %.3g, p=%.3g, AIC=%.2f" %
              (model_mu.params[1], model_mu.pvalues[1], model_mu.aic))
        print("  tAI-only: coef(tAI)= %.3g, p=%.3g, AIC=%.2f" %
              (model_tai.params[1], model_tai.pvalues[1], model_tai.aic))
        print("  μ+tAI:   coef(mu)= %.3g, p=%.3g; coef(tAI)= %.3g, p=%.3g; AIC=%.2f" %
              (model_both.params[1], model_both.pvalues[1],
               model_both.params[2], model_both.pvalues[2],
               model_both.aic))

if __name__ == "__main__":
    main()

