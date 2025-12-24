#!/usr/bin/env python

import argparse
import csv
from collections import defaultdict

def main():
    ap = argparse.ArgumentParser(
        description="Compute minimal codon capacity from observed μ–tAI modes (k_A)"
    )
    ap.add_argument("--aa_mode_summary_tsv", required=True,
                    help="errors/aa_mode_summary.tsv")
    args = ap.parse_args()

    kA = {}
    with open(args.aa_mode_summary_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            aa = row["aa"]
            n_modes = 0
            modes = row["modes_present"]
            if modes != "none":
                n_modes = len(modes.split(","))
            kA[aa] = n_modes

    total_modes = sum(kA.values())
    print("[INFO] k_A per amino acid:")
    for aa in sorted(kA.keys()):
        print(f"  {aa}: k_A = {kA[aa]}")
    print(f"\n[INFO] Sum_A k_A = {total_modes}")

    print("\n[INTERPRETATION]")
    print("Minimal codon capacity required (ignoring stops) is at least Sum_A k_A.")
    print("Any doublet code has at most 16 codons.")
    print("If Sum_A k_A > 16, no doublet code can realize all required μ–tAI modes.")
    print("Even restricted to the 9 two-codon amino acids, Sum_A k_A = 18 > 16,")
    print("so doublet codes are already impossible under the observed mode demand.")

if __name__ == "__main__":
    main()

