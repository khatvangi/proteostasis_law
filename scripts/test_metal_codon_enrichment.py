#!/usr/bin/env python

import csv
import numpy as np
from collections import defaultdict
from scipy.stats import fisher_exact, chi2_contingency

in_csv = "metals/metal_codon_bias_summary.csv"

# aa -> codon -> (ligand_count, background_count)
data = defaultdict(dict)

with open(in_csv, "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        aa = row["aa"]
        codon = row["codon"]
        lig = int(row["ligand_count"])
        bg = int(row["background_count"])
        data[aa][codon] = (lig, bg)

for aa, codons in data.items():
    if len(codons) != 2:
        print(f"Skipping {aa}: not exactly 2 codons")
        continue

    codon1, codon2 = sorted(codons.keys())  # alphabetical just for consistency
    a, b = codons[codon1]  # ligand, background
    c, d = codons[codon2]

    table = np.array([[a, b],
                      [c, d]])

    OR, p_fisher = fisher_exact(table)
    chi2, p_chi, dof, exp = chi2_contingency(table)

    print(f"\nAA: {aa}")
    print(f"  Codon1: {codon1} (ligand={a}, bg={b})")
    print(f"  Codon2: {codon2} (ligand={c}, bg={d})")
    print(f"  Fisher OR = {OR:.3f}, p = {p_fisher:.3e}")
    print(f"  Chi2 = {chi2:.3f}, p = {p_chi:.3e}")

