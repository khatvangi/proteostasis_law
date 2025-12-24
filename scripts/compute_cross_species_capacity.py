#!/usr/bin/env python3
"""
compute_cross_species_capacity.py

Compute and compare mode capacity across organisms.
Shows that doublet impossibility is universal.

Usage:
    python compute_cross_species_capacity.py \
        --ecoli results/ecoli_aa_mode_summary.tsv \
        --bsub results/bsub_aa_mode_summary.tsv \
        --yeast results/yeast_aa_mode_summary.tsv \
        --out results/cross_species_comparison.tsv
"""

import argparse
import csv
import sys

def load_mode_summary(path):
    """Load AA mode summary and compute k_A for each amino acid."""
    aa_modes = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            aa = row['aa']
            modes_str = row.get('modes_present', '')
            if modes_str and modes_str != 'none' and modes_str.strip():
                modes = [m.strip() for m in modes_str.split(',') if m.strip()]
                k_A = len(modes)
            else:
                k_A = 0

            n_codons = int(row.get('n_codons', 0))
            aa_modes[aa] = {
                'k_A': k_A,
                'n_codons': n_codons,
                'modes': modes_str
            }
    return aa_modes

def analyze_organism(name, summary_path):
    """Analyze mode capacity for one organism."""
    try:
        aa_modes = load_mode_summary(summary_path)
    except FileNotFoundError:
        return None
    except Exception as e:
        print(f"Error loading {summary_path}: {e}", file=sys.stderr)
        return None

    total_k = sum(d['k_A'] for d in aa_modes.values())

    # Count 2-codon AAs
    two_codon_aas = [aa for aa, d in aa_modes.items() if d['n_codons'] == 2]
    two_codon_k = sum(aa_modes[aa]['k_A'] for aa in two_codon_aas)

    # By degeneracy class
    by_deg = {}
    for aa, d in aa_modes.items():
        deg = d['n_codons']
        if deg not in by_deg:
            by_deg[deg] = {'count': 0, 'total_k': 0, 'aas': []}
        by_deg[deg]['count'] += 1
        by_deg[deg]['total_k'] += d['k_A']
        by_deg[deg]['aas'].append(aa)

    return {
        'name': name,
        'total_k': total_k,
        'n_aa': len(aa_modes),
        'two_codon_k': two_codon_k,
        'two_codon_count': len(two_codon_aas),
        'by_degeneracy': by_deg,
        'aa_modes': aa_modes
    }

def print_results(results):
    """Print formatted comparison results."""

    print("=" * 70)
    print("CROSS-SPECIES MODE CAPACITY COMPARISON")
    print("=" * 70)
    print(f"\n{'Organism':<20} {'Total k_A':<12} {'2-codon k':<12} {'> 16?':<10}")
    print("-" * 54)

    for r in results:
        if r is None:
            continue
        exceeds = "YES" if r['two_codon_k'] > 16 else "NO"
        print(f"{r['name']:<20} {r['total_k']:<12} {r['two_codon_k']:<12} {exceeds:<10}")

    print("\n" + "=" * 70)
    print("DOUBLET IMPOSSIBILITY TEST")
    print("=" * 70)
    print("\nFor doublets to work, 2-codon amino acids must fit in <= 16 codons.")
    print("But 2-codon AAs alone require:\n")

    all_impossible = True
    for r in results:
        if r is None:
            continue
        status = "IMPOSSIBLE" if r['two_codon_k'] > 16 else "possible"
        if r['two_codon_k'] <= 16:
            all_impossible = False
        avg_modes = r['two_codon_k'] / r['two_codon_count'] if r['two_codon_count'] > 0 else 0
        print(f"  {r['name']}: {r['two_codon_k']} modes "
              f"({r['two_codon_count']} AAs x {avg_modes:.1f} modes/AA) "
              f"-> {status}")

    print("\n" + "=" * 70)
    print("MODE COUNT BY DEGENERACY CLASS")
    print("=" * 70)

    for r in results:
        if r is None:
            continue
        print(f"\n{r['name']}:")
        for deg in sorted(r['by_degeneracy'].keys()):
            d = r['by_degeneracy'][deg]
            avg_k = d['total_k'] / d['count'] if d['count'] > 0 else 0
            print(f"  {deg}-codon AAs: {d['count']} AAs, "
                  f"total k = {d['total_k']}, avg k = {avg_k:.2f}")

    print("\n" + "=" * 70)
    print("CONCLUSION")
    print("=" * 70)

    if all_impossible:
        print("\n** DOUBLET CODES ARE IMPOSSIBLE IN ALL ORGANISMS **")
        print("\nThe mode capacity constraint is UNIVERSAL:")
        print("  - All organisms require more modes than doublets can provide")
        print("  - This is not an E. coli-specific artifact")
        print("  - Triplet codes are the minimal viable architecture")
    else:
        print("\nWarning: Some organisms might support doublet codes (check data)")

    print()

def main():
    parser = argparse.ArgumentParser(
        description="Compare mode capacity across organisms"
    )
    parser.add_argument("--ecoli", help="E. coli AA mode summary TSV")
    parser.add_argument("--bsub", help="B. subtilis AA mode summary TSV")
    parser.add_argument("--yeast", help="S. cerevisiae AA mode summary TSV")
    parser.add_argument("--out", default="cross_species_comparison.tsv",
                        help="Output comparison TSV")
    args = parser.parse_args()

    results = []

    if args.ecoli:
        r = analyze_organism("E. coli", args.ecoli)
        if r:
            results.append(r)

    if args.bsub:
        r = analyze_organism("B. subtilis", args.bsub)
        if r:
            results.append(r)

    if args.yeast:
        r = analyze_organism("S. cerevisiae", args.yeast)
        if r:
            results.append(r)

    if not results:
        print("Error: No valid input files provided", file=sys.stderr)
        sys.exit(1)

    # Print results
    print_results(results)

    # Write output TSV
    with open(args.out, 'w') as f:
        f.write("organism\ttotal_k_A\ttwo_codon_k\ttwo_codon_count\texceeds_16\n")
        for r in results:
            if r is None:
                continue
            exceeds = "YES" if r['two_codon_k'] > 16 else "NO"
            f.write(f"{r['name']}\t{r['total_k']}\t{r['two_codon_k']}\t"
                    f"{r['two_codon_count']}\t{exceeds}\n")

    print(f"Wrote comparison to {args.out}")

if __name__ == "__main__":
    main()
