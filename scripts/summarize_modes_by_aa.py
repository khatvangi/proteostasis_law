#!/usr/bin/env python

import argparse
import csv
from collections import defaultdict, Counter

def main():
    ap = argparse.ArgumentParser(
        description="Summarize μ–tAI modes per amino acid from codon_modes_ecoli.tsv"
    )
    ap.add_argument(
        "--codon_modes_tsv",
        required=True,
        help="errors/codon_modes_ecoli.tsv"
    )
    ap.add_argument(
        "--out_tsv",
        required=True,
        help="Output TSV: aa, n_codons, modes_present, n_safe_sprinter, n_safe_careful, n_risky_sprinter, n_risky_careful, n_unknown"
    )
    args = ap.parse_args()

    # aa -> list of mode labels for its codons
    aa_modes = defaultdict(list)

    with open(args.codon_modes_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            aa = row["aa"]
            mode = row["mode"]
            aa_modes[aa].append(mode)

    with open(args.out_tsv, "w") as out:
        cols = [
            "aa",
            "n_codons",
            "modes_present",
            "n_safe_sprinter",
            "n_safe_careful",
            "n_risky_sprinter",
            "n_risky_careful",
            "n_unknown",
        ]
        out.write("\t".join(cols) + "\n")

        for aa in sorted(aa_modes.keys()):
            modes = aa_modes[aa]
            n_codons = len(modes)
            counts = Counter(modes)

            n_safe_sprinter   = counts.get("safe_sprinter", 0)
            n_safe_careful    = counts.get("safe_careful", 0)
            n_risky_sprinter  = counts.get("risky_sprinter", 0)
            n_risky_careful   = counts.get("risky_careful", 0)
            n_unknown         = counts.get("unknown", 0)

            # Set of non-unknown modes present
            present = sorted(
                m for m in counts.keys()
                if m != "unknown" and counts[m] > 0
            )
            modes_present = ",".join(present) if present else "none"

            row_vals = [
                aa,
                str(n_codons),
                modes_present,
                str(n_safe_sprinter),
                str(n_safe_careful),
                str(n_risky_sprinter),
                str(n_risky_careful),
                str(n_unknown),
            ]
            out.write("\t".join(row_vals) + "\n")

    print(f"[INFO] Wrote AA mode summary to {args.out_tsv}")

if __name__ == "__main__":
    main()

