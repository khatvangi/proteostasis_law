#!/usr/bin/env python

import argparse
import pandas as pd

def main():
    ap = argparse.ArgumentParser(
        description="Extract E. coli codon-level error rates (μ) from Landerer Data_S2_error_detection_rate.xlsx"
    )
    ap.add_argument(
        "--xlsx",
        required=True,
        help="Path to Data_S2_error_detection_rate.xlsx"
    )
    ap.add_argument(
        "--out_tsv",
        required=True,
        help="Output TSV (codon\\tmu)"
    )
    args = ap.parse_args()

    # Read the Excel file, E. coli sheet
    df = pd.read_excel(args.xlsx, sheet_name="E. coli")

    # Expect columns: Codon, mean, median, sd, ...
    if "Codon" not in df.columns or "mean" not in df.columns:
        raise RuntimeError("Expected columns 'Codon' and 'mean' not found in E. coli sheet")

    # Keep only codon and mean error rate as μ
    out = df[["Codon", "mean"]].copy()
    out = out.rename(columns={"Codon": "codon", "mean": "mu"})

    # Sort codons for readability
    out = out.sort_values("codon")

    # Write TSV
    out.to_csv(args.out_tsv, sep="\t", index=False)
    print(f"[INFO] Wrote {len(out)} codon μ values to {args.out_tsv}")

if __name__ == "__main__":
    main()

