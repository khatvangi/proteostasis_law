#!/usr/bin/env python

import argparse
import csv
import random
from collections import defaultdict, Counter

def load_modes(path):
    # codon_modes_ecoli.tsv
    codon_to_mode = {}
    codon_to_aa = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            codon = row["codon"]
            aa = row["aa"]
            mode = row["mode"]
            codon_to_mode[codon] = mode
            codon_to_aa[codon] = aa
    return codon_to_mode, codon_to_aa

def load_aa_degeneracy(path):
    # From same file: count codons per AA
    degeneracy = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        by_aa = defaultdict(set)
        for row in reader:
            aa = row["aa"]
            codon = row["codon"]
            by_aa[aa].add(codon)
        for aa, codons in by_aa.items():
            degeneracy[aa] = len(codons)
    return degeneracy

def compute_n_modes_for_assignment(assign, codon_to_mode):
    # assign: aa -> list of codons
    # return: aa -> n_modes (non-unknown)
    result = {}
    for aa, codons in assign.items():
        modes = set()
        for c in codons:
            mode = codon_to_mode.get(c, "unknown")
            if mode != "unknown":
                modes.add(mode)
        result[aa] = len(modes)
    return result

def main():
    ap = argparse.ArgumentParser(
        description="Null model: random codon assignment vs n_modes-by-degeneracy"
    )
    ap.add_argument("--codon_modes_tsv", required=True,
                    help="errors/codon_modes_ecoli.tsv")
    ap.add_argument("--n_iter", type=int, default=1000,
                    help="Number of simulations")
    ap.add_argument("--out_tsv", required=True,
                    help="Output TSV summarizing null distributions")
    args = ap.parse_args()

    codon_to_mode, codon_to_aa_actual = load_modes(args.codon_modes_tsv)
    degeneracy = load_aa_degeneracy(args.codon_modes_tsv)

    all_codons = list(codon_to_mode.keys())
    aa_list = sorted(degeneracy.keys())

    # Observed n_modes
    observed_assign = defaultdict(list)
    for codon, aa in codon_to_aa_actual.items():
        observed_assign[aa].append(codon)
    observed_n_modes = compute_n_modes_for_assignment(observed_assign, codon_to_mode)

    # Null distributions: for each degeneracy value, collect mean n_modes
    null_by_deg = defaultdict(list)

    for it in range(args.n_iter):
        random.shuffle(all_codons)
        idx = 0
        assign = {}
        for aa in aa_list:
            n = degeneracy[aa]
            assign[aa] = all_codons[idx:idx+n]
            idx += n
        n_modes_sim = compute_n_modes_for_assignment(assign, codon_to_mode)
        for aa in aa_list:
            k = degeneracy[aa]
            null_by_deg[k].append(n_modes_sim[aa])

    # Summarize: per degeneracy, observed mean n_modes vs null distribution
    with open(args.out_tsv, "w") as out:
        out.write("n_codons\tobs_mean_n_modes\tnull_mean\tnull_sd\n")
        for k in sorted(set(degeneracy.values())):
            obs = [observed_n_modes[aa] for aa in aa_list if degeneracy[aa] == k]
            null_vals = null_by_deg[k]
            obs_mean = sum(obs)/len(obs)
            null_mean = sum(null_vals)/len(null_vals)
            null_sd = (sum((x-null_mean)**2 for x in null_vals)/(len(null_vals)-1))**0.5
            out.write(f"{k}\t{obs_mean:.3f}\t{null_mean:.3f}\t{null_sd:.3f}\n")

    print(f"[INFO] Wrote null summary to {args.out_tsv}")

if __name__ == "__main__":
    main()

