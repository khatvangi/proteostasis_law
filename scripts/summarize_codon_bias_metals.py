#!/usr/bin/env python

import argparse
import csv
import gzip
from collections import defaultdict

# Simple translation table for the codons we care about
CODON2AA = {
    "TGT": "CYS", "TGC": "CYS",
    "GAC": "ASP", "GAT": "ASP",
    "GAA": "GLU", "GAG": "GLU",
    "CAC": "HIS", "CAT": "HIS",
    # You can extend as needed
}

SIDECHAIN_ATOMS = {
    "CYS": {"SG"},
    "HIS": {"ND1", "NE2"},
    "ASP": {"OD1", "OD2"},
    "GLU": {"OE1", "OE2"},
}

TARGET_AAS = {"CYS", "HIS", "ASP", "GLU"}


def load_cds_by_gene(cds_fasta_path):
    """
    gene_name -> CDS nucleotide sequence
    using [gene=...] in header.
    """
    import re
    from collections import defaultdict

    gene_to_seq = defaultdict(str)
    open_func = gzip.open if cds_fasta_path.endswith(".gz") else open
    gene_pattern = re.compile(r"\[gene=([^\]]+)\]")

    with open_func(cds_fasta_path, "rt") as f:
        current_gene = None
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                current_gene = None
                m = gene_pattern.search(line[1:])
                if m:
                    current_gene = m.group(1).strip()
            else:
                if current_gene:
                    gene_to_seq[current_gene] += line.strip().upper()

    return gene_to_seq


def main():
    ap = argparse.ArgumentParser(
        description="Summarize codon bias at metal-ligand residues vs background for CYS/HIS/ASP/GLU."
    )
    ap.add_argument("--metal_codons_csv", required=True,
                    help="metal_sites_ecoli_with_codons.csv")
    ap.add_argument("--cds_fasta", required=True,
                    help="GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz")
    ap.add_argument("--out_csv", required=True,
                    help="Output summary CSV for AA/codon counts.")

    args = ap.parse_args()

    print(f"[INFO] Loading CDS by gene from {args.cds_fasta}")
    gene_to_seq = load_cds_by_gene(args.cds_fasta)
    print(f"[INFO] Loaded {len(gene_to_seq)} CDS sequences")

    # 1) Collect ligand codon counts
    ligand_counts = defaultdict(int)
    genes_with_ligands = set()

    print(f"[INFO] Reading metal sites from {args.metal_codons_csv}")
    with open(args.metal_codons_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            aa = row["resname"].strip().upper()
            if aa not in TARGET_AAS:
                continue

            atom = row["atom_name"].strip().upper()
            if atom not in SIDECHAIN_ATOMS[aa]:
                # restrict to sidechain donor atoms
                continue

            codon = row.get("codon", "").strip().upper()
            if len(codon) != 3:
                continue

            gene = row.get("gene_name", "").strip()
            if not gene:
                continue

            # sanity: codon translates to same AA
            if CODON2AA.get(codon) != aa:
                # skip weird mismatches
                continue

            key = (aa, codon)
            ligand_counts[key] += 1
            genes_with_ligands.add(gene)

    print(f"[INFO] Ligand genes: {len(genes_with_ligands)}")

    # 2) Background: all codons for those AAs in those genes
    background_counts = defaultdict(int)

    for gene in genes_with_ligands:
        seq = gene_to_seq.get(gene)
        if not seq:
            continue

        # step through codons in-frame
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3].upper()
            aa = CODON2AA.get(codon)
            if aa in TARGET_AAS:
                key = (aa, codon)
                background_counts[key] += 1

    # 3) Write summary
    with open(args.out_csv, "w", newline="") as f_out:
        fieldnames = ["aa", "codon", "ligand_count", "background_count"]
        writer = csv.DictWriter(f_out, fieldnames=fieldnames)
        writer.writeheader()

        all_keys = set(ligand_counts.keys()) | set(background_counts.keys())
        for aa, codon in sorted(all_keys):
            writer.writerow({
                "aa": aa,
                "codon": codon,
                "ligand_count": ligand_counts.get((aa, codon), 0),
                "background_count": background_counts.get((aa, codon), 0),
            })

    print(f"[INFO] Summary written to {args.out_csv}")


if __name__ == "__main__":
    main()

