#!/usr/bin/env python3
"""
Decompose dominated codons into:
  (a) Wobble partners of non-dominated codons (structurally necessary)
  (b) Box fillers (same NN_ box as non-dominated codon)
  (c) Historical excess (extra codon blocks beyond Pareto demand)

Usage:
  python3 decompose_dominated.py \
    --mu_tsv errors/codon_error_rates_FRESH.tsv \
    --tai_tsv errors/ecoli_tai_ws_FRESH.tsv
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

def wobble_partner(codon):
    """Return the wobble partner at position 3.
    C <-> T (pyrimidine pair), A <-> G (purine pair)."""
    partners = {"C": "T", "T": "C", "A": "G", "G": "A"}
    return codon[:2] + partners[codon[2]]

def same_box(c1, c2):
    """Two codons are in the same box if they share positions 1-2."""
    return c1[:2] == c2[:2]

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

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mu_tsv", required=True)
    ap.add_argument("--tai_tsv", required=True)
    args = ap.parse_args()

    mu = load_values(args.mu_tsv, "codon", "mu")
    tai = load_values(args.tai_tsv, "codon", "w_tai")

    # group codons by AA, compute Pareto fronts
    aa_codons = defaultdict(list)
    for codon, aa in CODON_TABLE.items():
        m = mu.get(codon)
        t = tai.get(codon)
        if m is not None and t is not None:
            aa_codons[aa].append((codon, m, t))

    # find non-dominated codons per AA
    non_dom_all = set()
    dom_all = set()
    aa_fronts = {}

    for aa, codons in aa_codons.items():
        non_dom = []
        dom = []
        for i, (ci, mui, taii) in enumerate(codons):
            dominated = False
            for j, (cj, muj, taij) in enumerate(codons):
                if i == j: continue
                if muj <= mui and taij >= taii and (muj < mui or taij > taii):
                    dominated = True
                    break
            if dominated:
                dom.append(ci)
                dom_all.add(ci)
            else:
                non_dom.append(ci)
                non_dom_all.add(ci)
        aa_fronts[aa] = (non_dom, dom)

    print("=" * 75)
    print("DECOMPOSITION OF 61 SENSE CODONS")
    print("=" * 75)
    print()
    print(f"Non-dominated (Pareto front): {len(non_dom_all)}")
    print(f"Dominated:                    {len(dom_all)}")
    print(f"Total with data:              {len(non_dom_all) + len(dom_all)}")
    print()

    # classify each dominated codon
    wobble_partners = []
    box_fillers = []
    excess_block = []

    for codon in sorted(dom_all):
        aa = CODON_TABLE[codon]
        non_dom, _ = aa_fronts[aa]

        # is wobble partner non-dominated?
        wp = wobble_partner(codon)
        if wp in non_dom_all and CODON_TABLE.get(wp) == aa:
            wobble_partners.append((codon, aa, f"wobble partner of {wp}"))
            continue

        # is any codon in the same box non-dominated?
        box_has_nd = any(same_box(codon, nd) and nd != codon for nd in non_dom)
        if box_has_nd:
            box_fillers.append((codon, aa, f"box filler (same NN_ as non-dom)"))
            continue

        # no non-dominated codon in this box
        excess_block.append((codon, aa, "excess block (no non-dom in same box)"))

    print("=" * 75)
    print("CATEGORY 1: WOBBLE PARTNERS OF NON-DOMINATED CODONS")
    print(f"({len(wobble_partners)} codons)")
    print("=" * 75)
    print("These MUST exist because tRNA decoding links C/T and A/G at pos 3.")
    print()
    for codon, aa, reason in wobble_partners:
        print(f"  {codon} ({aa}): {reason}")

    print()
    print("=" * 75)
    print("CATEGORY 2: BOX FILLERS")
    print(f"({len(box_fillers)} codons)")
    print("=" * 75)
    print("Same NN_ box as a non-dominated codon. Exist because boxes are")
    print("decoded by tRNA families that read all 4 position-3 variants.")
    print()
    for codon, aa, reason in box_fillers:
        print(f"  {codon} ({aa}): {reason}")

    print()
    print("=" * 75)
    print("CATEGORY 3: EXCESS BLOCKS")
    print(f"({len(excess_block)} codons)")
    print("=" * 75)
    print("No non-dominated codon shares their NN_ box. These entire blocks")
    print("are operationally unnecessary -- likely historical assignments from")
    print("biosynthetic pathway coevolution.")
    print()
    for codon, aa, reason in excess_block:
        print(f"  {codon} ({aa}): {reason}")

    print()
    print("=" * 75)
    print("SUMMARY TABLE")
    print("=" * 75)
    print(f"  Non-dominated (Pareto core):  {len(non_dom_all):>3}")
    print(f"  Wobble partners (physics):    {len(wobble_partners):>3}")
    print(f"  Box fillers (wobble blocks):  {len(box_fillers):>3}")
    print(f"  Excess blocks (historical):   {len(excess_block):>3}")
    print(f"  ─────────────────────────────────")
    total = len(non_dom_all) + len(wobble_partners) + len(box_fillers) + len(excess_block)
    print(f"  Total:                        {total:>3}")
    print()
    print("This decomposition shows that the 61 sense codons partition into:")
    print("  - An operational core forced by proteostasis (Pareto)")
    print("  - Structural partners forced by decoding physics (wobble)")
    print("  - Historical residue from code expansion (coevolution)")

if __name__ == "__main__":
    main()
