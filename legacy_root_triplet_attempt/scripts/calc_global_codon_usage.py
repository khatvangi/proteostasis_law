#!/usr/bin/env python

import argparse
import gzip
from collections import Counter, defaultdict

# Standard genetic code, DNA codons
CODON_TABLE = {
    # Phe
    "TTT": "F", "TTC": "F",
    # Leu
    "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    # Ile
    "ATT": "I", "ATC": "I", "ATA": "I",
    # Met
    "ATG": "M",
    # Val
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    # Ser
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "AGT": "S", "AGC": "S",
    # Pro
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    # Thr
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    # Ala
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    # Tyr
    "TAT": "Y", "TAC": "Y",
    # His
    "CAT": "H", "CAC": "H",
    # Gln
    "CAA": "Q", "CAG": "Q",
    # Asn
    "AAT": "N", "AAC": "N",
    # Lys
    "AAA": "K", "AAG": "K",
    # Asp
    "GAT": "D", "GAC": "D",
    # Glu
    "GAA": "E", "GAG": "E",
    # Cys
    "TGT": "C", "TGC": "C",
    # Trp
    "TGG": "W",
    # Arg
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGA": "R", "AGG": "R",
    # Gly
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    # Stops
    "TAA": "*", "TAG": "*", "TGA": "*",
}

def parse_fasta(handle):
    header = None
    seq_chunks = []
    for line in handle:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header, "".join(seq_chunks)
            header = line[1:]
            seq_chunks = []
        else:
            seq_chunks.append(line)
    if header is not None:
        yield header, "".join(seq_chunks)

def main():
    ap = argparse.ArgumentParser(
        description="Compute global codon usage from CDS FASTA"
    )
    ap.add_argument(
        "--cds_fasta",
        required=True,
        help="NCBI CDS FASTA (e.g. GCF_*_cds_from_genomic.fna.gz)",
    )
    ap.add_argument(
        "--out_tsv",
        required=True,
        help="Output TSV: aa, codon, count, freq_within_aa, freq_global",
    )
    args = ap.parse_args()

    # Open gz or plain
    if args.cds_fasta.endswith(".gz"):
        fh = gzip.open(args.cds_fasta, "rt")
    else:
        fh = open(args.cds_fasta, "r")

    codon_counts = Counter()
    aa_counts = Counter()
    total_codons = 0

    for hdr, seq in parse_fasta(fh):
        seq = seq.upper().replace("U", "T")
        # assume coding sequence, in frame
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if len(codon) != 3:
                continue
            if "N" in codon:
                continue
            aa = CODON_TABLE.get(codon)
            if aa is None:
                continue
            if aa == "*":
                # stop codon; skip
                continue
            codon_counts[codon] += 1
            aa_counts[aa] += 1
            total_codons += 1

    fh.close()

    # Write output
    with open(args.out_tsv, "w") as out:
        out.write("aa\tcodon\tcount\tfreq_within_aa\tfreq_global\n")
        for codon, c in sorted(codon_counts.items()):
            aa = CODON_TABLE[codon]
            within = c / aa_counts[aa] if aa_counts[aa] > 0 else 0.0
            global_freq = c / total_codons if total_codons > 0 else 0.0
            out.write(
                f"{aa}\t{codon}\t{c}\t{within:.8f}\t{global_freq:.8f}\n"
            )

    print(f"[INFO] Wrote global codon usage to {args.out_tsv}")
    print(f"[INFO] Total codons counted: {total_codons}")

if __name__ == "__main__":
    main()

