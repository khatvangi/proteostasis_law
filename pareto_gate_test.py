#!/usr/bin/env python3
"""
GATEKEEPING TEST: Pareto non-dominance for all 9 two-codon AAs.
Run on boron with FRESH data files.

Usage:
  python3 pareto_gate_test.py \
    --mu_tsv errors/codon_error_rates_FRESH.tsv \
    --tai_tsv errors/ecoli_tai_ws_FRESH.tsv

If the two-codon Pareto K >= 17, the Doublet Impossibility Lemma
survives under a non-tautological method.
If K <= 16, mode-counting is dead.
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
    "AGA":"Arg","AGG":"Arg",
    "AGT":"Ser","AGC":"Ser",
    "GGT":"Gly","GGC":"Gly","GGA":"Gly","GGG":"Gly",
}

TWO_CODON_AAS = {
    "Phe": ("TTT", "TTC"),
    "Tyr": ("TAT", "TAC"),
    "Cys": ("TGT", "TGC"),
    "His": ("CAT", "CAC"),
    "Gln": ("CAA", "CAG"),
    "Asn": ("AAT", "AAC"),
    "Lys": ("AAA", "AAG"),
    "Asp": ("GAT", "GAC"),
    "Glu": ("GAA", "GAG"),
}

def load_values(path, key_col, val_col):
    """Load a TSV into {key: float_value}"""
    out = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = row[key_col].strip().upper()
            try:
                val = float(row[val_col].strip())
                out[key] = val
            except (ValueError, KeyError):
                pass
    return out

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mu_tsv", required=True)
    ap.add_argument("--tai_tsv", required=True)
    args = ap.parse_args()

    mu = load_values(args.mu_tsv, "codon", "mu")
    tai = load_values(args.tai_tsv, "codon", "w_tai")

    # normalize codon names (handle T/U)
    mu_norm = {}
    for k, v in mu.items():
        mu_norm[k.replace("U", "T")] = v
        mu_norm[k] = v
    tai_norm = {}
    for k, v in tai.items():
        tai_norm[k.replace("U", "T")] = v
        tai_norm[k] = v

    print("=" * 70)
    print("PARETO NON-DOMINANCE TEST: ALL 9 TWO-CODON AMINO ACIDS")
    print("=" * 70)
    print()
    print("Method: codon A dominates B if mu(A) <= mu(B) AND tAI(A) >= tAI(B)")
    print("with at least one strict inequality.")
    print("k_A = number of non-dominated codons (= 2 if tradeoff, 1 if dominated)")
    print()

    total_K = 0
    n_tradeoff = 0
    n_dominated = 0
    missing = []

    for aa in sorted(TWO_CODON_AAS.keys()):
        c1_dna, c2_dna = TWO_CODON_AAS[aa]

        mu1 = mu_norm.get(c1_dna)
        mu2 = mu_norm.get(c2_dna)
        tai1 = tai_norm.get(c1_dna)
        tai2 = tai_norm.get(c2_dna)

        if mu1 is None or mu2 is None or tai1 is None or tai2 is None:
            missing.append((aa, c1_dna, c2_dna, mu1, mu2, tai1, tai2))
            print(f"{aa}: MISSING DATA")
            print(f"  {c1_dna}: mu={mu1}, tAI={tai1}")
            print(f"  {c2_dna}: mu={mu2}, tAI={tai2}")
            print()
            continue

        # pareto dominance test
        c1_dom = (mu1 <= mu2 and tai1 >= tai2) and (mu1 < mu2 or tai1 > tai2)
        c2_dom = (mu2 <= mu1 and tai2 >= tai1) and (mu2 < mu1 or tai2 > tai1)

        if c1_dom:
            k_A = 1
            relation = f"{c1_dna} DOMINATES {c2_dna}"
            n_dominated += 1
        elif c2_dom:
            k_A = 1
            relation = f"{c2_dna} DOMINATES {c1_dna}"
            n_dominated += 1
        else:
            k_A = 2
            relation = "TRADEOFF (neither dominates)"
            n_tradeoff += 1

        total_K += k_A

        # which axis favors which codon?
        mu_winner = c1_dna if mu1 < mu2 else (c2_dna if mu2 < mu1 else "TIE")
        tai_winner = c1_dna if tai1 > tai2 else (c2_dna if tai2 > tai1 else "TIE")

        print(f"{aa}: k_A = {k_A}  [{relation}]")
        print(f"  {c1_dna}: mu = {mu1:.2e}, tAI = {tai1:.4f}")
        print(f"  {c2_dna}: mu = {mu2:.2e}, tAI = {tai2:.4f}")
        print(f"  Lower mu: {mu_winner} | Higher tAI: {tai_winner}")
        print()

    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"Tradeoff (k_A=2):  {n_tradeoff} / {n_tradeoff + n_dominated}")
    print(f"Dominated (k_A=1): {n_dominated} / {n_tradeoff + n_dominated}")
    if missing:
        print(f"Missing data:      {len(missing)}")
    print(f"Two-codon Pareto K = {total_K}")
    print(f"Doublet ceiling    = 16")
    print()

    if total_K >= 17:
        print("* RESULT: DOUBLET IMPOSSIBILITY LEMMA SURVIVES")
        print(f"  Two-codon K = {total_K} > 16 under non-tautological method.")
        print("  Proceed with Option 2: restructure manuscript around")
        print("  Pareto non-dominance.")
    elif total_K == 16:
        print("# RESULT: EXACTLY AT THRESHOLD")
        print(f"  Two-codon K = 16 = doublet ceiling.")
        print("  Lemma fails (need strict inequality).")
        print("  But adding 4/6-codon AA contributions pushes K above 16.")
        print("  Check total K across all AAs.")
    else:
        print("X RESULT: DOUBLET IMPOSSIBILITY LEMMA FAILS")
        print(f"  Two-codon K = {total_K} < 16 under non-tautological method.")
        print("  Mode-counting cannot establish doublet impossibility.")
        print("  Proceed with Option 1 (PEC formula) or Option 3 (shelve).")

    # also compute total K across ALL amino acids under Pareto
    print()
    print("=" * 70)
    print("FULL PARETO K (all 20 amino acids)")
    print("=" * 70)
    print()
    print("For amino acids with >2 codons, Pareto non-dominance")
    print("counts how many codons are NOT dominated by any synonym.")
    print()

    # group all codons by AA
    aa_codons = defaultdict(list)
    for codon_dna, aa in CODON_TABLE.items():
        m = mu_norm.get(codon_dna)
        t = tai_norm.get(codon_dna)
        if m is not None and t is not None:
            aa_codons[aa].append((codon_dna, m, t))

    full_K = 0
    print(f"{'AA':<5} {'Deg':<5} {'k_A':<5} {'Non-dominated codons'}")
    print("-" * 60)
    for aa in sorted(aa_codons.keys()):
        codons = aa_codons[aa]
        n = len(codons)

        # find Pareto front
        non_dominated = []
        for i, (ci, mui, taii) in enumerate(codons):
            dominated = False
            for j, (cj, muj, taij) in enumerate(codons):
                if i == j:
                    continue
                if muj <= mui and taij >= taii and (muj < mui or taij > taii):
                    dominated = True
                    break
            if not dominated:
                non_dominated.append(ci)

        k_A = len(non_dominated)
        full_K += k_A
        nd_str = ", ".join(non_dominated)
        print(f"{aa:<5} {n:<5} {k_A:<5} {nd_str}")

    print()
    print(f"TOTAL PARETO K = {full_K}")
    print(f"Doublet capacity = 16")
    print(f"Triplet capacity = 64")
    print(f"Doublet impossible (full Pareto)? {'YES' if full_K > 16 else 'NO'}")

if __name__ == "__main__":
    main()
