#!/usr/bin/env python3
"""
how much translation-derived damage the measured error rates actually imply.

the earlier draft asserted that "roughly one-fifth to one-quarter of synthesized
proteins are expected to carry at least one misincorporation". that figure was
not derived from the mu data the paper itself uses. computing it directly from
the E. coli codon usage and the Landerer per-codon rates gives ~16-18%, i.e.
about one in six. this script produces that number so the manuscript can cite a
computed artifact instead of a recollection.

model: errors at successive codons are treated as independent, so for a protein
of N codons the chance of carrying at least one misincorporation is
1 - exp(-mubar * N) with mubar the codon-usage-weighted mean per-codon rate.
independence is an approximation -- error rates are correlated along a transcript
via local tRNA availability -- and it is stated as such in the manuscript.
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
COMP = ROOT / "data" / "computed"


def main():
    mu = pd.read_csv(RAW / "codon_error_rates_ecoli.tsv", sep="\t")
    usage = pd.read_csv(RAW / "global_codon_usage_ecoli.tsv", sep="\t")
    prot = json.loads((RAW / "arithmetic_results.json").read_text())
    length = prot["length_distribution"]["summary"]

    m = usage.merge(mu, on="codon")
    w = m["count"] / m["count"].sum()
    mubar = float((w * m.mu).sum())

    out = {
        "n_codons_with_both_mu_and_usage": int(len(m)),
        "usage_weighted_mean_mu_per_codon": mubar,
        "unweighted_mean_mu_per_codon": float(mu.mu.mean()),
        "mu_min": float(mu.mu.min()),
        "mu_max": float(mu.mu.max()),
        "mu_fold_span": float(mu.mu.max() / mu.mu.min()),
        "proteome_source": prot["length_distribution"]["source"],
        "proteome_n": int(length["n"]),
        "median_length_aa": float(length["median"]),
        "mean_length_aa": float(length["mean"]),
    }
    for label, N in (("median_length", length["median"]),
                     ("mean_length", length["mean"])):
        lam = mubar * N
        out[f"expected_errors_per_protein_at_{label}"] = float(lam)
        out[f"frac_proteins_with_ge1_error_at_{label}"] = float(1 - np.exp(-lam))

    (COMP / "translation_burden.json").write_text(json.dumps(out, indent=2))

    print(f"usage-weighted mean mu = {mubar:.3e} /codon "
          f"({out['n_codons_with_both_mu_and_usage']} codons)")
    print(f"mu span = {out['mu_fold_span']:.0f}-fold "
          f"({out['mu_min']:.2e} to {out['mu_max']:.2e})")
    for label in ("median_length", "mean_length"):
        print(f"  at {label} ({out['median_length_aa' if label=='median_length' else 'mean_length_aa']:.0f} aa): "
              f"E[errors] = {out[f'expected_errors_per_protein_at_{label}']:.3f}, "
              f"P(>=1) = {out[f'frac_proteins_with_ge1_error_at_{label}']:.1%}")
    print(f"\nwrote {COMP/'translation_burden.json'}")


if __name__ == "__main__":
    main()
