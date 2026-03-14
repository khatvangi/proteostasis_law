#!/usr/bin/env python

import argparse
import csv
import gzip
from collections import defaultdict

def load_uniprot_gene_map(uniprot_tsv_path):
    """
    Build a mapping: UniProt accession -> primary gene name
    from UP000000625_ecoliK12_metadata.tsv.gz.

    Expected columns:
    Entry   Gene Names (primary)    Organism    Length
    """
    acc_to_gene = {}

    open_func = gzip.open if uniprot_tsv_path.endswith(".gz") else open
    with open_func(uniprot_tsv_path, "rt") as f:
        header = f.readline().rstrip("\n").split("\t")
        col_idx = {name: i for i, name in enumerate(header)}

        entry_col = col_idx.get("Entry")
        gene_col = col_idx.get("Gene Names (primary)")

        if entry_col is None or gene_col is None:
            raise ValueError(
                f"Could not find 'Entry' or 'Gene Names (primary)' in {uniprot_tsv_path} header: {header}"
            )

        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            acc = parts[entry_col].strip()
            gene_name = parts[gene_col].strip()

            if not acc or not gene_name:
                continue

            # UniProt sometimes has multiple gene names separated by ' '
            primary_gene = gene_name.split()[0]
            acc_to_gene[acc] = primary_gene

    return acc_to_gene


def load_cds_by_gene(cds_fasta_path):
    """
    Build mapping gene_name -> CDS nucleotide sequence from
    GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz using [gene=...] in header.

    Example header:
    >lcl|NC_000913.3_cds_NP_414542.1_1 [gene=thrL] [locus_tag=b0001] ...
    """
    import re
    gene_to_seq = defaultdict(str)

    open_func = gzip.open if cds_fasta_path.endswith(".gz") else open
    gene_pattern = re.compile(r"\[gene=([^\]]+)\]")

    with open_func(cds_fasta_path, "rt") as f:
        current_gene = None
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                current_gene = None
                header = line[1:]  # drop '>'
                m = gene_pattern.search(header)
                if m:
                    current_gene = m.group(1).strip()
            else:
                if current_gene:
                    gene_to_seq[current_gene] += line.strip().upper()

    return gene_to_seq


def main():
    ap = argparse.ArgumentParser(
        description="Map metal-contact UniProt residues to codons using MG1655 CDS."
    )
    ap.add_argument("--metal_uniprot_csv", required=True,
                    help="Input CSV with metal sites mapped to UniProt (metal_sites_ecoli_with_uniprot.csv).")
    ap.add_argument("--uniprot_tsv", required=True,
                    help="UP000000625_ecoliK12_metadata.tsv.gz (for UniProt->gene mapping).")
    ap.add_argument("--cds_fasta", required=True,
                    help="GCF_000005845.2_ASM584v2_cds_from_genomic.fna.gz (for gene->CDS).")
    ap.add_argument("--out_csv", required=True,
                    help="Output CSV with gene_name, codon, codon_start included.")

    args = ap.parse_args()

    print(f"[INFO] Loading UniProt->gene map from {args.uniprot_tsv}")
    acc_to_gene = load_uniprot_gene_map(args.uniprot_tsv)
    print(f"[INFO] Loaded {len(acc_to_gene)} UniProt->gene mappings")

    print(f"[INFO] Loading CDS by gene from {args.cds_fasta}")
    gene_to_seq = load_cds_by_gene(args.cds_fasta)
    print(f"[INFO] Loaded {len(gene_to_seq)} gene->CDS sequences")

    print(f"[INFO] Reading metal sites from {args.metal_uniprot_csv}")

    with open(args.metal_uniprot_csv, "r") as f_in:
        reader = csv.DictReader(f_in)
        rows = list(reader)

    out_fields = list(rows[0].keys()) + ["gene_name", "codon", "codon_start"]
    mapped = 0
    unmapped = 0
    skipped = 0

    with open(args.out_csv, "w", newline="") as f_out:
        writer = csv.DictWriter(f_out, fieldnames=out_fields)
        writer.writeheader()

        for row in rows:
            uniprot_ac = row.get("uniprot_ac", "").strip()
            uniprot_pos_str = row.get("uniprot_pos", "").strip()

            if not uniprot_ac or not uniprot_pos_str:
                row["gene_name"] = ""
                row["codon"] = ""
                row["codon_start"] = ""
                writer.writerow(row)
                unmapped += 1
                continue

            try:
                uniprot_pos = int(float(uniprot_pos_str))
            except ValueError:
                row["gene_name"] = ""
                row["codon"] = ""
                row["codon_start"] = ""
                writer.writerow(row)
                skipped += 1
                continue

            gene_name = acc_to_gene.get(uniprot_ac)
            if not gene_name:
                row["gene_name"] = ""
                row["codon"] = ""
                row["codon_start"] = ""
                writer.writerow(row)
                unmapped += 1
                continue

            cds_seq = gene_to_seq.get(gene_name)
            if not cds_seq:
                row["gene_name"] = gene_name
                row["codon"] = ""
                row["codon_start"] = ""
                writer.writerow(row)
                unmapped += 1
                continue

            # 1-based residue index -> 0-based nucleotide index
            codon_start = (uniprot_pos - 1) * 3
            codon_end = codon_start + 3

            if codon_end > len(cds_seq):
                row["gene_name"] = gene_name
                row["codon"] = ""
                row["codon_start"] = ""
                writer.writerow(row)
                skipped += 1
                continue

            codon = cds_seq[codon_start:codon_end]

            row["gene_name"] = gene_name
            row["codon"] = codon
            row["codon_start"] = codon_start
            writer.writerow(row)
            mapped += 1

    print(f"[INFO] Done. Codon mapped for {mapped} sites.")
    print(f"[INFO] Unmapped: {unmapped}, Skipped (weird indices): {skipped}")
    print(f"[INFO] Output written to {args.out_csv}")


if __name__ == "__main__":
    main()

