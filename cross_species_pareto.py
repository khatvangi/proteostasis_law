#!/usr/bin/env python3
"""
Cross-species Pareto K test.

Since Landerer et al. only measured E. coli, this test uses E. coli mu values
for all species (same codon->mu mapping) but species-specific tAI values.

This is an explicit limitation: it tests whether the ARCHITECTURAL demand
(K_Pareto) is robust to changes in tRNA pools, not whether each species
independently generates the same mu landscape.
"""

import argparse
import csv
from collections import defaultdict

CODON_TABLE = {
    "TTT":"Phe","TTC":"Phe","TTA":"Leu","TTG":"Leu",
    "CTT":"Leu","CTC":"Leu","CTA":"Leu","CTG":"Leu",
    "ATT":"Ile","ATC":"Ile","ATA":"Ile","ATG":"Met",
    "GTT":"Val","GTC":"Val","GTA":"Val","GTG":"Val",
    "TCT":"Ser","TCC":"Ser","TCA":"Ser","TCG":"Ser",
    "CCT":"Pro","CCC":"Pro","CCA":"Pro","CCG":"Pro",
    "ACT":"Thr","ACC":"Thr","ACA":"Thr","ACG":"Thr",
    "GCT":"Ala","GCC":"Ala","GCA":"Ala","GCG":"Ala",
    "TAT":"Tyr","TAC":"Tyr",
    "CAT":"His","CAC":"His","CAA":"Gln","CAG":"Gln",
    "AAT":"Asn","AAC":"Asn","AAA":"Lys","AAG":"Lys",
    "GAT":"Asp","GAC":"Asp","GAA":"Glu","GAG":"Glu",
    "TGT":"Cys","TGC":"Cys","TGG":"Trp",
    "CGT":"Arg","CGC":"Arg","CGA":"Arg","CGG":"Arg",
    "AGA":"Arg","AGG":"Arg","AGT":"Ser","AGC":"Ser",
    "GGT":"Gly","GGC":"Gly","GGA":"Gly","GGG":"Gly",
}

def load_values(path, key_col, val_col):
    out = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = row[key_col].strip().upper().replace("U","T")
            try:
                out[key] = float(row[val_col].strip())
            except (ValueError, KeyError):
                pass
    return out

def compute_pareto_K(mu_map, tai_map):
    """Compute Pareto K for one species."""
    aa_codons = defaultdict(list)
    for codon, aa in CODON_TABLE.items():
        m = mu_map.get(codon)
        t = tai_map.get(codon)
        if m is not None and t is not None:
            aa_codons[aa].append((codon, m, t))

    total_K = 0
    aa_results = {}
    for aa in sorted(aa_codons.keys()):
        codons = aa_codons[aa]
        non_dom = []
        for i, (ci, mui, taii) in enumerate(codons):
            dominated = False
            for j, (cj, muj, taij) in enumerate(codons):
                if i == j: continue
                if muj <= mui and taij >= taii and (muj < mui or taij > taii):
                    dominated = True
                    break
            if not dominated:
                non_dom.append(ci)
        k_A = len(non_dom)
        total_K += k_A
        aa_results[aa] = (len(codons), k_A, non_dom)

    return total_K, aa_results

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mu_tsv", required=True, help="E. coli mu values (used for all species)")
    ap.add_argument("--tai_ecoli", required=True)
    ap.add_argument("--tai_bsub", required=True)
    ap.add_argument("--tai_scer", required=True)
    args = ap.parse_args()

    mu = load_values(args.mu_tsv, "codon", "mu")

    species = {
        "E. coli":       load_values(args.tai_ecoli, "codon", "w_tai"),
        "B. subtilis":   load_values(args.tai_bsub, "codon", "w_tai"),
        "S. cerevisiae": load_values(args.tai_scer, "codon", "w_tai"),
    }

    print("=" * 70)
    print("CROSS-SPECIES PARETO K TEST")
    print("=" * 70)
    print()
    print("mu values: E. coli (Landerer et al.) for all species")
    print("tAI values: species-specific from GtRNAdb tRNA gene counts")
    print()

    for sp_name, tai_map in species.items():
        K, results = compute_pareto_K(mu, tai_map)
        print(f"--- {sp_name} ---")
        print(f"{'AA':<5} {'Deg':<5} {'k_A':<5} {'Non-dominated'}")
        for aa in sorted(results.keys()):
            deg, k, nd = results[aa]
            print(f"{aa:<5} {deg:<5} {k:<5} {', '.join(nd)}")
        print(f"\nTotal Pareto K = {K}")
        print(f"Doublet impossible? {'YES' if K > 16 else 'NO'} ({K} vs 16)")
        print()

    # comparison table
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"{'Species':<20} {'K_Pareto':<12} {'vs 16':<10} {'Verdict'}")
    print("-" * 55)
    for sp_name, tai_map in species.items():
        K, _ = compute_pareto_K(mu, tai_map)
        verdict = "IMPOSSIBLE" if K > 16 else "FEASIBLE"
        cmp = ">" if K > 16 else "<="
        print(f"{sp_name:<20} {K:<12} {cmp} 16     {verdict}")

    print()
    print("NOTE: Same mu values used for all species.")
    print("Variation in K reflects different tRNA pools (tAI).")
    print("This tests architectural robustness, not independent mu measurement.")

if __name__ == "__main__":
    main()
