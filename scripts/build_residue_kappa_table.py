#!/usr/bin/env python

import argparse
import csv

def load_metal_sites(path):
    # metals/metal_sites_ecoli_with_uniprot.csv
    # columns: pdb_id,chain_id,resseq,icode,resname,atom_name,metal_element,metal_serial,distance,cif_path,uniprot_ac,uniprot_pos
    metal = set()
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            up = row["uniprot_ac"]
            pos = row["uniprot_pos"]
            if not up or not pos:
                continue
            metal.add((up, int(pos)))
    return metal

def load_mu(path):
    mu = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mu[row["codon"]] = float(row["mu"])
    return mu

def load_tai(path):
    tai = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            tai[row["codon"]] = float(row["w_tai"])
    return tai

def main():
    ap = argparse.ArgumentParser(
        description="Build residue-level table with κ (metal vs non-metal)"
    )
    ap.add_argument("--uniprot_pos_codon_tsv", required=True,
                    help="errors/uniprot_pos_codon_ecoli.tsv")
    ap.add_argument("--metal_uniprot_csv", required=True,
                    help="metals/metal_sites_ecoli_with_uniprot.csv")
    ap.add_argument("--mu_tsv", required=True,
                    help="errors/codon_error_rates.tsv")
    ap.add_argument("--tai_tsv", required=True,
                    help="errors/ecoli_tai_ws.tsv")
    ap.add_argument("--out_tsv", required=True,
                    help="Output TSV: residue-level κ table")
    args = ap.parse_args()

    metal_sites = load_metal_sites(args.metal_uniprot_csv)
    mu_map = load_mu(args.mu_tsv)
    tai_map = load_tai(args.tai_tsv)

    with open(args.uniprot_pos_codon_tsv) as f_in, \
         open(args.out_tsv, "w") as f_out:

        reader = csv.DictReader(f_in, delimiter="\t")
        cols = [
            "uniprot_ac","pos","aa","codon",
            "mu","w_tai","is_metal","kappa"
        ]
        writer = csv.DictWriter(f_out, fieldnames=cols, delimiter="\t")
        writer.writeheader()

        for row in reader:
            up = row["uniprot_ac"]
            pos = int(row["pos"])
            aa = row["aa"]
            codon = row["codon"]
            mu = mu_map.get(codon, float("nan"))
            tai = tai_map.get(codon, float("nan"))
            is_metal = 1 if (up, pos) in metal_sites else 0
            kappa = 1.0 if is_metal == 1 else 0.0  # first-pass κ

            writer.writerow({
                "uniprot_ac": up,
                "pos": pos,
                "aa": aa,
                "codon": codon,
                "mu": mu,
                "w_tai": tai,
                "is_metal": is_metal,
                "kappa": kappa,
            })

    print(f"[INFO] Wrote residue κ table to {args.out_tsv}")

if __name__ == "__main__":
    main()

