#!/usr/bin/env python

import argparse
import csv
from collections import defaultdict
import numpy as np
from scipy.stats import fisher_exact


def load_mu_table(mu_tsv_path):
    """
    Load codon -> mu from a TSV with columns:
    codon    mu
    """
    codon_to_mu = {}
    with open(mu_tsv_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            codon = row["codon"].strip().upper()
            mu = float(row["mu"])
            codon_to_mu[codon] = mu
    return codon_to_mu


def load_codon_counts(summary_csv_path):
    """
    Load aa/codon ligand/background counts from
    metal_codon_bias_summary.csv
    """
    data = defaultdict(dict)
    with open(summary_csv_path, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            aa = row["aa"].strip().upper()
            codon = row["codon"].strip().upper()
            lig = int(row["ligand_count"])
            bg = int(row["background_count"])
            data[aa][codon] = (lig, bg)
    return data


def main():
    ap = argparse.ArgumentParser(
        description="Combine metal codon counts with codon-specific error rates (mu)."
    )
    ap.add_argument(
        "--summary_csv",
        required=True,
        help="metals/metal_codon_bias_summary.csv",
    )
    ap.add_argument(
        "--mu_tsv",
        required=True,
        help="errors/codon_error_rates.tsv (codon, mu)",
    )
    ap.add_argument(
        "--out_tsv",
        required=True,
        help="Output TSV with per-AA codon and weighted mu stats.",
    )

    args = ap.parse_args()

    print(f"[INFO] Loading codon error rates from {args.mu_tsv}")
    codon_to_mu = load_mu_table(args.mu_tsv)
    print(f"[INFO] Loaded mu for {len(codon_to_mu)} codons")

    print(f"[INFO] Loading codon counts from {args.summary_csv}")
    aa_data = load_codon_counts(args.summary_csv)

    rows_out = []

    for aa, codons in sorted(aa_data.items()):
        if len(codons) != 2:
            print(f"[WARN] Skipping {aa}: expected 2 synonymous codons, got {len(codons)}")
            continue

        # Sort codons alphabetically just for consistent output
        codon1, codon2 = sorted(codons.keys())

        lig1, bg1 = codons[codon1]
        lig2, bg2 = codons[codon2]

        # Retrieve mu; warn if missing
        if codon1 not in codon_to_mu or codon2 not in codon_to_mu:
            print(f"[WARN] Missing mu for {aa}: {codon1} or {codon2}")
            continue

        mu1 = codon_to_mu[codon1]
        mu2 = codon_to_mu[codon2]

        total_lig = lig1 + lig2
        total_bg = bg1 + bg2

        if total_lig == 0 or total_bg == 0:
            print(f"[WARN] No counts for {aa}, skipping.")
            continue

        # Weighted average mu at ligand and background sites
        mu_lig = (mu1 * lig1 + mu2 * lig2) / total_lig
        mu_bg = (mu1 * bg1 + mu2 * bg2) / total_bg

        # Enrichment OR (same as test_metal_codon_enrichment)
        table = np.array([[lig1, bg1],
                          [lig2, bg2]])
        OR, p_fisher = fisher_exact(table)

        # Figure out which codon is enriched at ligands
        # Compare ligand fraction vs background fraction
        frac1_lig = lig1 / total_lig
        frac1_bg = bg1 / total_bg

        if frac1_lig > frac1_bg:
            enriched_codon = codon1
            enriched_mu = mu1
            depleted_codon = codon2
            depleted_mu = mu2
        else:
            enriched_codon = codon2
            enriched_mu = mu2
            depleted_codon = codon1
            depleted_mu = mu1

        rows_out.append({
            "aa": aa,
            "codon1": codon1,
            "codon2": codon2,
            "ligand_codon1": lig1,
            "ligand_codon2": lig2,
            "background_codon1": bg1,
            "background_codon2": bg2,
            "mu_codon1": mu1,
            "mu_codon2": mu2,
            "mu_ligand_weighted": mu_lig,
            "mu_background_weighted": mu_bg,
            "OR_enrichment": OR,
            "p_fisher": p_fisher,
            "enriched_codon": enriched_codon,
            "enriched_mu": enriched_mu,
            "depleted_codon": depleted_codon,
            "depleted_mu": depleted_mu,
            "delta_mu_enriched_minus_depleted": enriched_mu - depleted_mu,
            "delta_mu_lig_minus_bg": mu_lig - mu_bg,
        })

    # Write TSV
    with open(args.out_tsv, "w", newline="") as f_out:
        fieldnames = [
            "aa",
            "codon1", "codon2",
            "ligand_codon1", "ligand_codon2",
            "background_codon1", "background_codon2",
            "mu_codon1", "mu_codon2",
            "mu_ligand_weighted", "mu_background_weighted",
            "OR_enrichment", "p_fisher",
            "enriched_codon", "enriched_mu",
            "depleted_codon", "depleted_mu",
            "delta_mu_enriched_minus_depleted",
            "delta_mu_lig_minus_bg",
        ]
        writer = csv.DictWriter(f_out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows_out:
            writer.writerow(row)

    print(f"[INFO] Wrote summary with mu to {args.out_tsv}")


if __name__ == "__main__":
    main()

