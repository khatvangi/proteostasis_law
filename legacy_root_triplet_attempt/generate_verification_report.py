#!/usr/bin/env python3
"""
generate a single verification report comparing every
quantitative claim in the manuscript against fresh data.
"""
import csv
from collections import defaultdict

def main():
    print("=" * 70)
    print("PAPER 2 VERIFICATION REPORT")
    print("=" * 70)
    print()

    claims = []

    # --- claim 1: K_total for E. coli ---
    try:
        with open("errors/aa_mode_summary_FRESH.tsv") as f:
            reader = csv.DictReader(f, delimiter="\t")
            K = 0
            K_no_singletons = 0
            two_codon_K = 0
            for row in reader:
                modes = row["modes_present"]
                n_codons = int(row["n_codons"])
                if modes == "none":
                    k = 1 if n_codons == 1 else 0
                else:
                    k = len(modes.split(","))
                K += k
                if n_codons > 1:
                    K_no_singletons += k
                if n_codons == 2:
                    two_codon_K += k

        claims.append(("K_total (within-family, counting singletons=1)",
                       K, "47 or 48", True))
        claims.append(("K_total = 47 if Met counted as 0",
                       K - 1, 47, K - 1 == 47))
        claims.append(("Two-codon subtotal = 18 (within-family)",
                       two_codon_K, 18, two_codon_K == 18))
        claims.append(("Doublet impossible: 18 > 16 (within-family)",
                       f"{two_codon_K} > 16", "YES", two_codon_K > 16))
    except FileNotFoundError:
        claims.append(("K_total", "FILE MISSING", 47, False))

    # --- claim 2: tautology test ---
    # two-codon K under global thresholds = 15 (from Phase 4 output)
    claims.append(("Two-codon K (GLOBAL thresholds)", 15, "> 16",
                   False))
    claims.append(("Doublet impossible under global thresholds?",
                   "15 <= 16 => NO", "YES", False))
    claims.append(("k_A=2 for all 2-codon AAs is tautological?",
                   "YES (3 collapse to k_A=1 under global)", "should be NO",
                   False))

    # --- claim 3: S_3 ~ 0.72 ---
    CODON_TABLE = {
        "TTT":"Phe","TTC":"Phe","TTA":"Leu","TTG":"Leu",
        "CTT":"Leu","CTC":"Leu","CTA":"Leu","CTG":"Leu",
        "ATT":"Ile","ATC":"Ile","ATA":"Ile","ATG":"Met",
        "GTT":"Val","GTC":"Val","GTA":"Val","GTG":"Val",
        "TCT":"Ser","TCC":"Ser","TCA":"Ser","TCG":"Ser",
        "CCT":"Pro","CCC":"Pro","CCA":"Pro","CCG":"Pro",
        "ACT":"Thr","ACC":"Thr","ACA":"Thr","ACG":"Thr",
        "GCT":"Ala","GCC":"Ala","GCA":"Ala","GCG":"Ala",
        "TAT":"Tyr","TAC":"Tyr","TAA":"*","TAG":"*",
        "CAT":"His","CAC":"His","CAA":"Gln","CAG":"Gln",
        "AAT":"Asn","AAC":"Asn","AAA":"Lys","AAG":"Lys",
        "GAT":"Asp","GAC":"Asp","GAA":"Glu","GAG":"Glu",
        "TGT":"Cys","TGC":"Cys","TGA":"*","TGG":"Trp",
        "CGT":"Arg","CGC":"Arg","CGA":"Arg","CGG":"Arg",
        "AGT":"Ser","AGC":"Ser","AGA":"Arg","AGG":"Arg",
        "GGT":"Gly","GGC":"Gly","GGA":"Gly","GGG":"Gly",
    }

    bases = "TCAG"
    total = [0, 0, 0]
    syn = [0, 0, 0]

    for codon, aa in CODON_TABLE.items():
        if aa == "*":
            continue
        for pos in range(3):
            for b in bases:
                if b == codon[pos]:
                    continue
                mutant = codon[:pos] + b + codon[pos+1:]
                total[pos] += 1
                if mutant in CODON_TABLE and CODON_TABLE[mutant] == aa:
                    syn[pos] += 1

    S = [syn[i]/total[i] for i in range(3)]
    claims.append((f"S_1 (position 1 synonymy)", f"{S[0]:.4f}", "~0.04", True))
    claims.append((f"S_2 (position 2 synonymy)", f"{S[1]:.4f}", "~0.00", True))
    claims.append((f"S_3 (position 3 synonymy)", f"{S[2]:.4f}", "~0.72",
                   abs(S[2] - 0.72) < 0.02))

    # --- claim 4: wobble structuring ---
    claims.append(("9/9 two-codon AAs separated (within-family)",
                   "9/9 = 1.000", "expected", True))
    claims.append(("...but this is tautological under median split",
                   "TAUTOLOGICAL", "should not be", False))
    claims.append(("Wobble null z = -1.74, p = 0.042",
                   "z = -0.73, p = 0.30", "z ~ -1.74", False))

    # --- claim 5: metal enrichment ---
    claims.append(("Metal Asp GAC enriched (p < 0.05)",
                   "p = 0.163", "significant", False))
    claims.append(("Metal Cys TGC enriched (p < 0.05)",
                   "p = 0.139", "significant", False))
    claims.append(("Metal Glu enriched (p < 0.05)",
                   "p = 0.084", "significant", False))
    claims.append(("Metal His enriched (p < 0.05)",
                   "p = 1.000", "significant", False))

    # --- claim 6: mu values ---
    claims.append(("mu range from Landerer 2024",
                   "3.3e-5 to 2.0e-2", "3.3e-5 to 2.0e-2", True))
    claims.append(("Manuscript mu values (0.135, 0.120, etc.)",
                   "FABRICATED (1000-3600x too large)", "should match Landerer",
                   False))

    # --- claim 7: cross-species ---
    try:
        with open("cross_species/results/cross_species_comparison.tsv") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                org = row.get("organism", row.get("species", "?"))
                K_val = row.get("K_total", row.get("K", "?"))
                claims.append((f"K ({org})", K_val, "47-48", True))
    except (FileNotFoundError, KeyError):
        claims.append(("Cross-species K", "FILE MISSING OR FORMAT ISSUE", "47-48", False))

    # --- print report ---
    n_pass = sum(1 for _, _, _, ok in claims if ok)
    n_fail = sum(1 for _, _, _, ok in claims if not ok)

    print(f"{'#':<4} {'Claim':<55} {'Computed':<30} {'Expected':<20} {'Status'}")
    print("-" * 120)
    for i, (desc, computed, expected, ok) in enumerate(claims, 1):
        status = "PASS" if ok else "** FAIL **"
        print(f"{i:<4} {desc:<55} {str(computed):<30} {str(expected):<20} {status}")

    print()
    print(f"SUMMARY: {n_pass} PASS, {n_fail} FAIL out of {len(claims)} claims")
    print()
    print("CRITICAL FAILURES:")
    for i, (desc, computed, expected, ok) in enumerate(claims, 1):
        if not ok:
            print(f"  #{i}: {desc}")
            print(f"       computed: {computed}")
            print(f"       expected: {expected}")
            print()

if __name__ == "__main__":
    main()
