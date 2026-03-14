#!/usr/bin/env python3
"""
calculate_tai.py

Calculate tRNA Adaptation Index (tAI) for each codon.
Supports using E. coli tAI values as proxy for other organisms.

Rationale for E. coli proxy:
- Ribosome structure is highly conserved across bacteria
- Wobble pairing rules are universal
- Mode STRUCTURE depends more on relative tAI differences than absolute values

Usage:
    # Use E. coli proxy (recommended for cross-species analysis)
    python calculate_tai.py --use_ecoli_proxy --out data/bsub/tai_values.tsv

    # Use custom tRNA copy numbers
    python calculate_tai.py --trna data/bsub/trna_copy_numbers.tsv --out data/bsub/tai_values.tsv
"""

import argparse
import csv
import sys

# E. coli tAI proxy values (normalized, from stAIcalc)
# These capture the relative translation efficiency differences
ECOLI_TAI_PROXY = {
    'TTT': 0.295, 'TTC': 0.500, 'TTA': 0.098, 'TTG': 0.167,
    'TCT': 0.167, 'TCC': 0.295, 'TCA': 0.098, 'TCG': 0.167,
    'TAT': 0.295, 'TAC': 0.500, 'TGT': 0.098, 'TGC': 0.167,
    'TGG': 0.500,
    'CTT': 0.098, 'CTC': 0.098, 'CTA': 0.038, 'CTG': 1.000,
    'CCT': 0.098, 'CCC': 0.098, 'CCA': 0.167, 'CCG': 0.500,
    'CAT': 0.167, 'CAC': 0.295, 'CAA': 0.295, 'CAG': 0.500,
    'CGT': 0.795, 'CGC': 0.500, 'CGA': 0.038, 'CGG': 0.098,
    'ATT': 0.500, 'ATC': 0.795, 'ATA': 0.038, 'ATG': 0.500,
    'ACT': 0.295, 'ACC': 0.500, 'ACA': 0.098, 'ACG': 0.167,
    'AAT': 0.295, 'AAC': 0.500, 'AAA': 0.795, 'AAG': 0.167,
    'AGT': 0.098, 'AGC': 0.167, 'AGA': 0.038, 'AGG': 0.038,
    'GTT': 0.500, 'GTC': 0.295, 'GTA': 0.295, 'GTG': 0.295,
    'GCT': 0.295, 'GCC': 0.500, 'GCA': 0.500, 'GCG': 0.500,
    'GAT': 0.295, 'GAC': 0.500, 'GAA': 0.795, 'GAG': 0.295,
    'GGT': 0.500, 'GGC': 0.500, 'GGA': 0.098, 'GGG': 0.167,
}

# Standard codon table
CODON_TO_AA = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}


def write_ecoli_proxy(output_file):
    """Write E. coli tAI proxy values to output file."""
    with open(output_file, 'w') as f:
        f.write("codon\tw_tai\n")
        for codon in sorted(ECOLI_TAI_PROXY.keys()):
            f.write(f"{codon}\t{ECOLI_TAI_PROXY[codon]}\n")

    print(f"[INFO] Wrote E. coli tAI proxy values to {output_file}", file=sys.stderr)
    print(f"[INFO] Using E. coli as proxy (ribosome structure conserved)", file=sys.stderr)


def reverse_complement(seq):
    """Return reverse complement of DNA sequence."""
    comp = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return ''.join(comp.get(b, b) for b in reversed(seq))


def get_wobble_penalty(codon_pos3, anticodon_pos1):
    """
    Return wobble penalty s for codon position 3 vs anticodon position 1.

    Wobble penalties (dos Reis et al. 2004):
    - Watson-Crick pairing: s = 0
    - G:U wobble: s = 0.5
    - I:N wobble (from A->I editing): s = 0.5
    """
    # Watson-Crick: no penalty
    wc_pairs = {('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')}
    if (codon_pos3, anticodon_pos1) in wc_pairs:
        return 0.0
    # G:U wobble
    if (codon_pos3, anticodon_pos1) in {('T', 'G'), ('G', 'T')}:
        return 0.5
    # I:N wobble (inosine from A)
    if anticodon_pos1 == 'I':
        if codon_pos3 in ['T', 'C', 'A']:
            return 0.5
    return 1.0  # No recognition


def calculate_tai_from_trna(trna_file, output_file):
    """
    Calculate tAI for all codons given tRNA copy numbers.

    tAI formula (dos Reis et al. 2004):
      w_i = sum over anticodons j of: (1 - s_ij) * tRNA_j
    """
    # Load tRNA copy numbers (anticodon -> count)
    trna_counts = {}
    with open(trna_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            anticodon = row['anticodon'].upper()
            count = int(row['copy_number'])
            trna_counts[anticodon] = count

    # Calculate w_i for each codon
    codon_weights = {}
    for codon in CODON_TO_AA:
        if CODON_TO_AA[codon] == '*':
            continue

        w_i = 0.0
        codon_pos3 = codon[2]

        # Check all anticodons
        for anticodon, count in trna_counts.items():
            if len(anticodon) != 3:
                continue
            # Anticodon is 3'->5', codon is 5'->3'
            # Anticodon position 1 pairs with codon position 3
            anticodon_pos1 = anticodon[0]

            # Check if this anticodon can recognize this codon
            anticodon_rc = reverse_complement(anticodon)
            if anticodon_rc[:2] == codon[:2]:  # Simplified check
                s = get_wobble_penalty(codon_pos3, anticodon_pos1)
                w_i += (1 - s) * count

        codon_weights[codon] = w_i

    # Normalize to [0, 1]
    max_w = max(codon_weights.values()) if codon_weights else 1
    if max_w > 0:
        for codon in codon_weights:
            codon_weights[codon] /= max_w

    # Write output
    with open(output_file, 'w') as f:
        f.write("codon\tw_tai\n")
        for codon in sorted(codon_weights):
            f.write(f"{codon}\t{codon_weights[codon]:.6f}\n")

    print(f"[INFO] Wrote tAI values to {output_file}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Calculate tRNA Adaptation Index (tAI) for codons"
    )
    parser.add_argument("--trna",
                        help="tRNA copy number file (anticodon, copy_number)")
    parser.add_argument("--use_ecoli_proxy", action="store_true",
                        help="Use E. coli tAI values as proxy (recommended for cross-species)")
    parser.add_argument("--out", required=True,
                        help="Output tAI file")
    args = parser.parse_args()

    if args.use_ecoli_proxy:
        write_ecoli_proxy(args.out)
    elif args.trna:
        calculate_tai_from_trna(args.trna, args.out)
    else:
        print("Error: Must specify either --use_ecoli_proxy or --trna",
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
