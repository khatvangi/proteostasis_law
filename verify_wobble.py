#!/usr/bin/env python3
"""
verify the wobble structuring claims:
  1. 2.5x Watson-Crick/wobble separation ratio (p=0.033)
  2. matched-null z = -1.74, p = 0.042

recomputed from scratch using the FRESH mode table.
"""
import csv
import random
from collections import defaultdict
import statistics

def main():
    random.seed(42)

    # load mode table
    records = []
    with open("errors/codon_modes_ecoli_FRESH.tsv") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            records.append({
                "aa": row["aa"],
                "codon": row["codon"],
                "mode": row["mode"],
                "mu_class": row["mu_class"],
                "tai_class": row["tai_class"],
                "mu": float(row["mu"]) if row["mu"] and row["mu"] != "None" else None,
                "tai": float(row["w_tai"]) if row["w_tai"] and row["w_tai"] != "None" else None,
            })

    # group by amino acid
    by_aa = defaultdict(list)
    for r in records:
        by_aa[r["aa"]].append(r)

    # === TEST 1: WC vs wobble mode separation for 2-codon AAs ===
    two_codon_aas = {aa: codons for aa, codons in by_aa.items() if len(codons) == 2}
    print(f"Two-codon amino acids: {len(two_codon_aas)}")
    print()

    n_separated = 0
    n_total = 0

    for aa, codons in sorted(two_codon_aas.items()):
        c1, c2 = codons
        if c1["mode"] == "unknown" or c2["mode"] == "unknown":
            print(f"  {aa}: SKIPPED (unknown mode)")
            continue

        n_total += 1
        separated = c1["mode"] != c2["mode"]
        if separated:
            n_separated += 1

        def wc_status(codon):
            last = codon[-1]
            return "WC" if last in "CG" else "Wob"

        print(f"  {aa}: {c1['codon']} ({wc_status(c1['codon'])}, {c1['mode']}) "
              f"vs {c2['codon']} ({wc_status(c2['codon'])}, {c2['mode']}) "
              f"-> {'SEPARATED' if separated else 'SAME MODE'}")

    if n_total > 0:
        obs_frac = n_separated / n_total
        print(f"\nSeparated: {n_separated}/{n_total} = {obs_frac:.3f}")
        print(f"NOTE: with n=2 within-family median split, separation is GUARANTEED")
        print(f"      (unless mu1==mu2 or tai1==tai2 exactly)")
        print(f"The '2.5x ratio' and p=0.033 need independent verification of what")
        print(f"exactly was being measured.")
    print()

    # === TEST 2: matched-null ensemble ===
    # observed total K (within-family, excluding Met)
    obs_K = 0
    for aa, codons in by_aa.items():
        modes = set(c["mode"] for c in codons if c["mode"] != "unknown")
        if modes:
            obs_K += len(modes)
        elif len(codons) == 1:
            obs_K += 1  # singleton

    print(f"=== MATCHED-NULL ANALYSIS ===")
    print(f"Observed K = {obs_K}")

    # null model: shuffle (mu, tai) VALUES across all 61 codons,
    # then recompute within-family medians and modes
    mu_tai_pairs = [(r["mu"], r["tai"]) for r in records if r["mu"] is not None and r["tai"] is not None]
    # also track which records have valid mu/tai
    valid_indices = [i for i, r in enumerate(records) if r["mu"] is not None and r["tai"] is not None]

    N_NULL = 10000
    null_Ks = []

    for trial in range(N_NULL):
        # shuffle mu-tai pairs
        shuffled_pairs = list(mu_tai_pairs)
        random.shuffle(shuffled_pairs)

        # reassign to records
        shuffled_records = []
        pair_idx = 0
        for i, r in enumerate(records):
            new_r = dict(r)
            if i in valid_indices:
                new_r["mu"], new_r["tai"] = shuffled_pairs[pair_idx]
                pair_idx += 1
            shuffled_records.append(new_r)

        # group by AA
        null_by_aa = defaultdict(list)
        for r in shuffled_records:
            null_by_aa[r["aa"]].append(r)

        # compute within-family medians and modes for each AA
        null_K = 0
        for aa, codons in null_by_aa.items():
            valid = [c for c in codons if c["mu"] is not None and c["tai"] is not None]
            if not valid:
                if len(codons) == 1:
                    null_K += 1
                continue

            if len(valid) == 1:
                null_K += 1
                continue

            # within-family median split
            mus = [c["mu"] for c in valid]
            tais = [c["tai"] for c in valid]
            mu_med = statistics.median(mus)
            tai_med = statistics.median(tais)

            modes = set()
            for c in valid:
                mu_cls = "low" if c["mu"] <= mu_med else "high"
                tai_cls = "high" if c["tai"] >= tai_med else "low"

                if mu_cls == "low" and tai_cls == "high":
                    mode = "safe_sprinter"
                elif mu_cls == "low" and tai_cls == "low":
                    mode = "safe_careful"
                elif mu_cls == "high" and tai_cls == "high":
                    mode = "risky_sprinter"
                else:
                    mode = "risky_careful"
                modes.add(mode)

            null_K += len(modes)

        null_Ks.append(null_K)

    null_mean = statistics.mean(null_Ks)
    null_std = statistics.stdev(null_Ks)
    z = (obs_K - null_mean) / null_std if null_std > 0 else 0

    # p-values
    p_lower = sum(1 for k in null_Ks if k <= obs_K) / N_NULL
    p_upper = sum(1 for k in null_Ks if k >= obs_K) / N_NULL

    print(f"Null mean K = {null_mean:.1f} +/- {null_std:.2f}")
    print(f"z-score = {z:.2f}")
    print(f"p (obs <= null) = {p_lower:.4f}")
    print(f"p (obs >= null) = {p_upper:.4f}")
    print()
    print(f"Manuscript claims: z = -1.74, p = 0.042")
    print(f"Recomputed:        z = {z:.2f}, p(lower) = {p_lower:.4f}")
    match = abs(z - (-1.74)) < 1.0
    print(f"Consistent? {'APPROXIMATELY' if match else '**MISMATCH**'}")
    print()

    # distribution summary
    print(f"Null K distribution:")
    print(f"  min = {min(null_Ks)}, max = {max(null_Ks)}")
    print(f"  5th percentile = {sorted(null_Ks)[500]}")
    print(f"  95th percentile = {sorted(null_Ks)[9500]}")

if __name__ == "__main__":
    main()
