#!/usr/bin/env python3
"""
validate the E. coli tAI (translational supply, nu) axis before anything uses it.

why this script exists
----------------------
the previous version of this project used `codon-deployment/data/computed/
ecoli_tai_ws.tsv` as nu. that file is corrupt: 45 of its 61 values are
bit-identical to the yeast file, and it assigns E. coli's dominant Leu codon
CTG a tAI of 0.059 (near the floor) while giving the minor Lys codon AAG the
maximum 1.00. the script that supposedly produced it does not reproduce it and
implements the dos Reis wobble rules with permuted indices.

so we do not trust any tAI vector on assertion. we validate it against an
independent, local, objective signal: translational selection.

the test
--------
tAI measures how well a codon is served by the tRNA pool. translational
selection predicts that highly expressed genes shift their synonymous codon
usage TOWARD well-served codons. therefore, for a valid tAI vector:

    delta_usage(codon) = f_within_aa(ribosomal genes) - f_within_aa(all genes)

must be positively correlated with tAI. ribosomal protein genes (rpl*/rps*/rpm*)
are the canonical high-expression set in E. coli and need no expression data.

this is a falsifiable pass/fail check, not a plausibility argument. it is run on
both the candidate reference vector and the known-corrupt vector, and the
corrupt vector is expected to fail. that comparison is what licenses the choice.
"""
import json
import re
from pathlib import Path

import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "computed"

BASES = "TCAG"
_AA = ("FFLLSSSSYY**CC*W" "LLLLPPPPHHQQRRRR"
       "IIIMTTTTNNKKSSRR" "VVVVAAAADDEEGGGG")
CODON_AA = {
    b1 + b2 + b3: _AA[i]
    for i, (b1, b2, b3) in enumerate(
        (a, b, c) for a in BASES for b in BASES for c in BASES)
}
SENSE = [c for c, a in CODON_AA.items() if a != "*"]

# ribosomal protein gene-name prefixes: the canonical E. coli high-expression set
RIBO_RE = re.compile(r"^(rpl|rps|rpm)", re.I)


def parse_cds(path):
    """yield (gene_name, sequence) for every CDS in an NCBI-style fasta."""
    gene, seq = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if gene is not None and seq:
                    yield gene, "".join(seq)
                m = re.search(r"\[gene=([^\]]+)\]", line)
                gene, seq = (m.group(1) if m else ""), []
            else:
                seq.append(line.strip().upper())
    if gene is not None and seq:
        yield gene, "".join(seq)


def codon_counts(path, only_ribosomal=False):
    """count sense codons over all CDS, or only over ribosomal protein genes."""
    counts = {c: 0 for c in SENSE}
    n_genes = 0
    for gene, s in parse_cds(path):
        if only_ribosomal and not RIBO_RE.match(gene):
            continue
        n_genes += 1
        for i in range(0, len(s) - 2, 3):
            c = s[i:i + 3]
            if c in counts:
                counts[c] += 1
    return counts, n_genes


def within_aa_freq(counts):
    """relative synonymous codon usage: fraction of its amino acid's codons."""
    df = pd.DataFrame({"codon": list(counts), "n": [counts[c] for c in counts]})
    df["aa"] = df.codon.map(CODON_AA)
    df["f"] = df.groupby("aa").n.transform(lambda x: x / x.sum())
    return df.set_index("codon").f


def load_tai(path):
    df = pd.read_csv(path, sep="\t")
    return dict(zip(df[df.columns[0]], df[df.columns[1]]))


def validate(name, tai, delta, f_all, f_ribo):
    """run the translational-selection test on one candidate tAI vector."""
    codons = [c for c in SENSE if c in tai]
    t = [tai[c] for c in codons]

    rho_delta, p_delta = stats.spearmanr(t, [delta[c] for c in codons])
    rho_all, _ = stats.spearmanr(t, [f_all[c] for c in codons])
    rho_ribo, _ = stats.spearmanr(t, [f_ribo[c] for c in codons])

    # per-amino-acid: does the higher-tAI synonym gain usage in ribosomal genes?
    gained = total = 0
    for aa in sorted(set(CODON_AA[c] for c in codons)):
        syn = [c for c in codons if CODON_AA[c] == aa]
        if len(syn) < 2:
            continue
        hi = max(syn, key=lambda c: tai[c])
        total += 1
        if delta[hi] > 0:
            gained += 1

    # sign test on that count
    #
    # NB on which statistic actually discriminates: the pooled correlation
    # rho(tAI, delta_usage) is ~+0.31 for BOTH the reference and the corrupt
    # vector, so on its own it is not diagnostic -- a vector can get the pooled
    # trend right while assigning the wrong synonym the higher weight inside
    # many amino acids. the per-amino-acid sign test is the discriminating
    # criterion, because Delta_A is computed strictly within amino acids and
    # that is exactly the resolution the downstream analysis depends on.
    p_sign = stats.binomtest(gained, total, 0.5, alternative="greater").pvalue
    verdict = "PASS" if (rho_delta > 0 and p_delta < 0.05 and p_sign < 0.05) else "FAIL"

    return {
        "vector": name,
        "rho_tai_vs_delta_usage": round(float(rho_delta), 4),
        "p_delta": float(f"{p_delta:.3g}"),
        "rho_tai_vs_usage_all_genes": round(float(rho_all), 4),
        "rho_tai_vs_usage_ribosomal": round(float(rho_ribo), 4),
        "n_aa_where_high_tai_synonym_gains": f"{gained}/{total}",
        "p_sign_test": float(f"{p_sign:.3g}"),
        "verdict": verdict,
    }


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cds = RAW / "ecoli_k12_cds.fna"

    all_counts, n_all = codon_counts(cds)
    ribo_counts, n_ribo = codon_counts(cds, only_ribosomal=True)
    print(f"CDS parsed: {n_all} genes total, {n_ribo} ribosomal protein genes")

    f_all = within_aa_freq(all_counts)
    f_ribo = within_aa_freq(ribo_counts)
    delta = (f_ribo - f_all).to_dict()

    candidates = {
        "reference (R tAI package, E. coli K-12)": RAW / "ecoli_tai_reference.tsv",
        "CORRUPT (codon-deployment/data/computed)": RAW / "ecoli_tai_CORRUPT_for_regression_test.tsv",
    }

    report = [validate(name, load_tai(p), delta, f_all, f_ribo)
              for name, p in candidates.items()]

    print("\n" + "=" * 78)
    print("translational-selection validation of candidate tAI vectors")
    print("=" * 78)
    for r in report:
        print(f"\n  {r['vector']}")
        print(f"    Spearman(tAI, usage shift in ribosomal genes) = "
              f"{r['rho_tai_vs_delta_usage']:+.3f}  (p = {r['p_delta']})")
        print(f"    Spearman(tAI, within-AA usage, all genes)     = "
              f"{r['rho_tai_vs_usage_all_genes']:+.3f}")
        print(f"    Spearman(tAI, within-AA usage, ribosomal)     = "
              f"{r['rho_tai_vs_usage_ribosomal']:+.3f}")
        print(f"    amino acids where the higher-tAI synonym gains usage: "
              f"{r['n_aa_where_high_tai_synonym_gains']}  (sign test p = {r['p_sign_test']})")
        print(f"    --> {r['verdict']}")

    chosen = report[0]
    if chosen["verdict"] != "PASS":
        raise SystemExit(
            "\nFATAL: the reference tAI vector failed the translational-selection "
            "test. Do not proceed -- no validated nu axis is available.")
    if report[1]["verdict"] == "PASS":
        raise SystemExit(
            "\nFATAL: the known-corrupt vector PASSED, so this test does not "
            "discriminate. Strengthen the test before relying on it.")

    # publish the validated vector under a neutral name for downstream scripts
    ref = pd.read_csv(candidates["reference (R tAI package, E. coli K-12)"], sep="\t")
    ref.columns = ["codon", "nu_tai"]
    ref = ref[ref.codon.isin(SENSE)].sort_values("codon")
    ref.to_csv(OUT / "nu_tai_ecoli_validated.tsv", sep="\t", index=False)

    (OUT / "tai_validation_report.json").write_text(json.dumps(report, indent=2))
    print(f"\nwrote {OUT/'nu_tai_ecoli_validated.tsv'} ({len(ref)} codons)")
    print(f"wrote {OUT/'tai_validation_report.json'}")

    top = ref.nlargest(4, "nu_tai")
    print("\n  sanity: highest-tAI codons = "
          + ", ".join(f"{r.codon}={r.nu_tai:.2f}" for r in top.itertuples()))


if __name__ == "__main__":
    main()
