#!/usr/bin/env python

import os
import argparse
import csv
import pandas as pd


def load_sifts_segments(sifts_csv_path):
    """
    Load uniprot_segments_observed.csv.gz from SIFTS and
    keep only relevant columns.
    """
    df = pd.read_csv(
        sifts_csv_path,
        comment="#",
        compression="gzip"
    )

    # Normalize column names if needed
    # Expecting at least: PDB, CHAIN, SP_PRIMARY, RES_BEG, RES_END, SP_BEG, SP_END
    df = df.rename(columns={
        "PDB": "pdb_id",
        "CHAIN": "chain_id",
        "SP_PRIMARY": "uniprot_ac",
        "RES_BEG": "resseq_start",
        "RES_END": "resseq_end",
        "SP_BEG": "uniprot_start",
        "SP_END": "uniprot_end",
    })

    # Keep only those columns to avoid surprises
    cols = [
        "pdb_id", "chain_id", "uniprot_ac",
        "resseq_start", "resseq_end",
        "uniprot_start", "uniprot_end"
    ]
    df = df[cols]

    # Ensure ids formatted consistently
    df["pdb_id"] = df["pdb_id"].str.upper().str.strip()
    df["chain_id"] = df["chain_id"].astype(str).str.strip()

    return df


def map_row_to_uniprot(row, sifts_df):
    """
    Given a metal-contact row and SIFTS segments dataframe,
    find UniProt accession and residue index.

    Returns (uniprot_ac, uniprot_pos) or (None, None) if not mappable.
    """
    pdb_id = row["pdb_id"].upper()
    chain_id = row["chain_id"].strip()
    resseq = int(row["resseq"])
    icode = row["icode"]

    # For now, skip residues with insertion codes
    if isinstance(icode, str) and icode.strip() != "":
        return None, None

    # Subset to matching PDB + chain
    subset = sifts_df[
        (sifts_df["pdb_id"] == pdb_id) &
        (sifts_df["chain_id"] == chain_id) &
        (sifts_df["resseq_start"] <= resseq) &
        (sifts_df["resseq_end"] >= resseq)
    ]

    if subset.empty:
        return None, None

    # If multiple segments match (rare), take the first
    seg = subset.iloc[0]

    offset = resseq - int(seg["resseq_start"])
    uniprot_pos = int(seg["uniprot_start"]) + offset
    uniprot_ac = seg["uniprot_ac"]

    return uniprot_ac, uniprot_pos


def main():
    ap = argparse.ArgumentParser(
        description="Map metal-contacting residues to UniProt using SIFTS segments."
    )
    ap.add_argument("--metal_csv", required=True,
                    help="Input CSV from find_metal_sites.py (metal_sites_raw_ecoli.csv).")
    ap.add_argument("--sifts_segments", required=True,
                    help="Path to uniprot_segments_observed.csv.gz from SIFTS.")
    ap.add_argument("--out_csv", required=True,
                    help="Output CSV with UniProt mapping added.")

    args = ap.parse_args()

    print(f"[INFO] Loading SIFTS segments from {args.sifts_segments}")
    sifts_df = load_sifts_segments(args.sifts_segments)

    print(f"[INFO] Loading metal-contact residues from {args.metal_csv}")
    metals = []
    with open(args.metal_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            metals.append(row)

    print(f"[INFO] Mapping {len(metals)} metal-contact residues to UniProt")

    out_fields = list(metals[0].keys()) + ["uniprot_ac", "uniprot_pos"]

    mapped = 0
    unmapped = 0

    with open(args.out_csv, "w", newline="") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=out_fields)
        writer.writeheader()

        for row in metals:
            uniprot_ac, uniprot_pos = map_row_to_uniprot(row, sifts_df)

            if uniprot_ac is None:
                unmapped += 1
                row["uniprot_ac"] = ""
                row["uniprot_pos"] = ""
            else:
                mapped += 1
                row["uniprot_ac"] = uniprot_ac
                row["uniprot_pos"] = uniprot_pos

            writer.writerow(row)

    print(f"[INFO] Done. Mapped: {mapped}, Unmapped: {unmapped}")
    print(f"[INFO] Output written to {args.out_csv}")


if __name__ == "__main__":
    main()

