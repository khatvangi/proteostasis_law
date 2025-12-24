#!/usr/bin/env python

import argparse
import gzip
import re
from collections import defaultdict

CODON_TABLE = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","AGT":"S","AGC":"S","CCT":"P","CCC":"P",
    "CCA":"P","CCG":"P","ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A",
    "GCA":"A","GCG":"A","TAT":"Y","TAC":"Y","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R","AGA":"R",
    "AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
    "TAA":"*","TAG":"*","TGA":"*",
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
        description="Build UniProt position→codon table from CDS FASTA + GFF"
    )
    ap.add_argument("--cds_fasta", required=True,
                    help="GCF_*_cds_from_genomic.fna.gz")
    ap.add_argument("--gff", required=True,
                    help="GCF_*_genomic.gff.gz")
    ap.add_argument("--out_tsv", required=True,
                    help="Output TSV: uniprot_ac, pos, aa, codon")
    args = ap.parse_args()

    # 1) Parse GFF to map protein_id -> UniProt AC
    prot_to_uniprot = {}
    prot_to_strand = {}
    with gzip.open(args.gff, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            ftype = parts[2]
            if ftype != "CDS":
                continue
            attr = parts[8]
            # Dbxref=UniProtKB/Swiss-Prot:P0A7F3,...
            m_db = re.search(r"Dbxref=[^;]*UniProtKB/Swiss-Prot:([^,;]+)", attr)
            m_pid = re.search(r"protein_id=([^;]+)", attr)
            m_strand = parts[6]
            if m_db and m_pid:
                up = m_db.group(1)
                pid = m_pid.group(1)
                prot_to_uniprot[pid] = up
                prot_to_strand[pid] = m_strand

    print(f"[INFO] Mapped {len(prot_to_uniprot)} protein_ids to UniProt ACs")

    # 2) Parse CDS FASTA: header carries protein_id
    rows = []
    if args.cds_fasta.endswith(".gz"):
        fh = gzip.open(args.cds_fasta, "rt")
    else:
        fh = open(args.cds_fasta, "r")

    for hdr, seq in parse_fasta(fh):
        # header e.g. lcl|NC_000913.3_cds_NP_414542.1_1 [gene=thrL] ...
        m_pid = re.search(r"\[protein_id=([^]]+)\]", hdr)
        if not m_pid:
            continue
        pid = m_pid.group(1)
        if pid not in prot_to_uniprot:
            continue
        up = prot_to_uniprot[pid]
        strand = prot_to_strand.get(pid, "+")
        seq = seq.upper().replace("U", "T")

        # assume CDS in coding orientation; if strand is "-", reverse complement
        # BUT NCBI CDS FASTA is usually in coding (5'->3' mRNA-like) orientation,
        # so we do NOT reverse-complement here.
        pos = 0
        for i in range(0, len(seq) - 2, 3):
            codon = seq[i:i+3]
            if "N" in codon:
                aa = "X"
            else:
                aa = CODON_TABLE.get(codon, "X")
            if aa == "*":
                # stop
                break
            pos += 1  # 1-based AA index
            rows.append((up, pos, aa, codon))

    fh.close()
    print(f"[INFO] Built {len(rows)} UniProt position→codon rows")

    # 3) Write TSV
    with open(args.out_tsv, "w") as out:
        out.write("uniprot_ac\tpos\taa\tcodon\n")
        for up, pos, aa, codon in rows:
            out.write(f"{up}\t{pos}\t{aa}\t{codon}\n")

    print(f"[INFO] Wrote {len(rows)} rows to {args.out_tsv}")

if __name__ == "__main__":
    main()

