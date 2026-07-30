#!/usr/bin/env python3
"""
reproduce the two results REMOVED from this paper, so that the Discussion's
account of why they were removed is itself checkable.

part A -- metal-binding-site codon enrichment
    recompute the four Fisher tests under two backgrounds:
      genome-wide  (as previously published)
      within-gene  (expression-matched: non-metal positions of the same genes)
    and report the site-annotation position concordance rate.

part B -- cross-species conservation of operational geometry
    show (i) that the published Delta_A used E. coli mu for ALL THREE species,
    by reproducing the published table exactly, and (ii) that the three
    "species-specific" tAI vectors are not independent.

caveat on part B: recomputing the cross-species correlations requires tRNA gene
counts for B. subtilis and S. cerevisiae. the only counts available locally are
the hardcoded vectors in the prior project's own script, which are themselves
unverified (its E. coli vector sums to 116 tRNA genes where K-12 has ~86). we
therefore recompute using those counts with the CANONICAL dos Reis formula and
label the result as indicative, not authoritative. what IS authoritative here is
the shared-mu artifact and the vector-identity count, both of which are exact.
"""
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
COMP = ROOT / "data" / "computed"

BASES = "TCAG"
_AA = ("FFLLSSSSYY**CC*W" "LLLLPPPPHHQQRRRR"
       "IIIMTTTTNNKKSSRR" "VVVVAAAADDEEGGGG")
CODONS_64 = [a + b + c for a in BASES for b in BASES for c in BASES]
CODON_AA = {CODONS_64[i]: _AA[i] for i in range(64)}
SENSE = [c for c in CODONS_64 if CODON_AA[c] != "*"]

AA3 = {"ASP": "D", "CYS": "C", "GLU": "E", "HIS": "H"}
PAIRS = {"D": ("GAC", "GAT"), "C": ("TGC", "TGT"),
         "E": ("GAA", "GAG"), "H": ("CAC", "CAT")}

# what the earlier draft reported, for side-by-side comparison
PUBLISHED = {"D": (1.349, 0.0070), "C": (1.282, 0.0298),
             "E": (1.356, 0.0289), "H": (1.399, 0.0035)}


# --------------------------------------------------------------------------
# part A
# --------------------------------------------------------------------------
def load_cds_by_uniprot():
    out, acc, seq = {}, None, []
    with open(RAW / "ecoli_k12_cds.fna") as fh:
        for line in fh:
            if line.startswith(">"):
                if acc and seq:
                    out[acc] = "".join(seq)
                m = re.search(r"UniProtKB/(?:Swiss-Prot|TrEMBL):(\w+)", line)
                acc, seq = (m.group(1) if m else None), []
            else:
                seq.append(line.strip().upper())
    if acc and seq:
        out[acc] = "".join(seq)
    return out


def part_a():
    cds = load_cds_by_uniprot()

    d = pd.read_csv(RAW / "metal_sites_ecoli_with_codons.csv")
    d = d[d.resname.isin(AA3)].copy()
    d["aa"] = d.resname.map(AA3)
    d = d[pd.to_numeric(d.uniprot_pos, errors="coerce").notna()]
    d["uniprot_pos"] = pd.to_numeric(d.uniprot_pos).astype(int)
    sites = d.drop_duplicates(["uniprot_ac", "uniprot_pos"])

    rows = []
    for r in sites.itertuples(index=False):
        s = cds.get(r.uniprot_ac)
        if s is None or 3 * r.uniprot_pos > len(s):
            continue
        c = s[3 * (r.uniprot_pos - 1):3 * r.uniprot_pos]
        if c not in CODON_AA:
            continue
        rows.append({"acc": r.uniprot_ac, "pos": r.uniprot_pos, "aa": r.aa,
                     "codon_annotation": r.codon, "codon_cds": c,
                     "aa_from_cds": CODON_AA[c]})
    sd = pd.DataFrame(rows)

    concordant = sd[sd.aa_from_cds == sd.aa].copy()
    diagnostics = {
        "unique_metal_residues_in_annotation": int(len(sites)),
        "resolvable_against_cds": int(len(sd)),
        "cds_codon_encodes_annotated_residue": int((sd.aa_from_cds == sd.aa).sum()),
        "position_concordance_frac": float((sd.aa_from_cds == sd.aa).mean()),
        "annotation_codon_matches_cds_codon_frac":
            float((sd.codon_annotation == sd.codon_cds).mean()),
        "cds_verified_sites_used": int(len(concordant)),
        "metalloprotein_genes": int(concordant.acc.nunique()),
    }

    site_keys = set(zip(concordant.acc, concordant.pos))
    metal_genes = set(concordant.acc)
    within, genomewide = {}, {}
    for acc, s in cds.items():
        in_metal_gene = acc in metal_genes
        for i in range(1, len(s) // 3 + 1):
            c = s[3 * (i - 1):3 * i]
            aa = CODON_AA.get(c)
            if aa not in PAIRS:
                continue
            genomewide[(aa, c)] = genomewide.get((aa, c), 0) + 1
            if in_metal_gene and (acc, i) not in site_keys:
                within[(aa, c)] = within.get((aa, c), 0) + 1

    results = []
    for aa, (enr, alt) in PAIRS.items():
        g = concordant[concordant.aa == aa]
        se, sa = int((g.codon_cds == enr).sum()), int((g.codon_cds == alt).sum())
        row = {"amino_acid": aa, "enriched_codon": enr, "alt_codon": alt,
               "n_sites_enriched": se, "n_sites_alt": sa,
               "published_or": PUBLISHED[aa][0], "published_p": PUBLISHED[aa][1]}
        for label, bg in (("genomewide", genomewide), ("within_gene", within)):
            be, ba = bg.get((aa, enr), 0), bg.get((aa, alt), 0)
            OR, p = fisher_exact([[se, sa], [be, ba]])
            row[f"{label}_or"] = float(OR)
            row[f"{label}_p"] = float(p)
            row[f"{label}_bg_n"] = int(be + ba)
        results.append(row)

    df = pd.DataFrame(results)
    df.to_csv(COMP / "removed_metal_site_test.tsv", sep="\t", index=False)

    print("=" * 92)
    print("PART A -- metal-binding-site enrichment under two backgrounds")
    print("=" * 92)
    print(f"  annotated D/C/E/H metal residues        : "
          f"{diagnostics['unique_metal_residues_in_annotation']}")
    print(f"  position concordance with the CDS       : "
          f"{diagnostics['position_concordance_frac']:.1%}  "
          f"(<-- ~40% of sites are mispositioned)")
    print(f"  annotation codon == CDS codon at that position: "
          f"{diagnostics['annotation_codon_matches_cds_codon_frac']:.1%}  "
          f"(so the bug is uniprot_pos, not the codon lookup)")
    print(f"  CDS-verified sites used                 : "
          f"{diagnostics['cds_verified_sites_used']} "
          f"in {diagnostics['metalloprotein_genes']} genes\n")
    print(f"  {'AA':<4}{'sites':>12} | {'genome-wide bg':^24} | "
          f"{'within-gene bg':^24} | published")
    for r in results:
        print(f"  {r['amino_acid']:<4}{r['n_sites_enriched']:>5}/"
              f"{r['n_sites_alt']:<6} | "
              f"OR={r['genomewide_or']:5.3f} p={r['genomewide_p']:7.4f}"
              f"{'*' if r['genomewide_p'] < .05 else ' '}      | "
              f"OR={r['within_gene_or']:5.3f} p={r['within_gene_p']:7.4f}"
              f"{'*' if r['within_gene_p'] < .05 else ' '}      | "
              f"{r['published_or']:.3f} / {r['published_p']:.4f}")
    n_sig = sum(r["within_gene_p"] < 0.05 for r in results)
    print(f"\n  significant under the expression-matched background: {n_sig} of 4")
    return diagnostics, df


# --------------------------------------------------------------------------
# part B
# --------------------------------------------------------------------------
S_WOBBLE = [0.0, 0.0, 0.0, 0.0, 0.41, 0.28, 0.9999, 0.68, 0.89]


def canonical_tai(trna, sking=1):
    """dos Reis et al. 2004 weights, matching R tAI::get.ws index conventions."""
    p = [1 - s for s in S_WOBBLE]
    W = np.zeros(64)
    for b in range(0, 64, 4):
        W[b + 0] = p[0] * trna[b + 0] + p[4] * trna[b + 1]   # NNT
        W[b + 1] = p[1] * trna[b + 1] + p[5] * trna[b + 0]   # NNC
        W[b + 2] = p[2] * trna[b + 2] + p[6] * trna[b + 1]   # NNA
        W[b + 3] = p[3] * trna[b + 3] + p[7] * trna[b + 2]   # NNG
    if sking == 1:
        W[34] = 1 - S_WOBBLE[8]                              # ATA, lysidine
    keep = [i for i in range(64) if CODON_AA[CODONS_64[i]] != "*" and i != 35]
    w = np.array([W[i] for i in keep])
    w = w / w.max()
    nz = w[w > 0]
    if len(nz) < len(w):
        w[w == 0] = np.exp(np.mean(np.log(nz)))
    return dict(zip([CODONS_64[i] for i in keep], w))


def tRNA_vectors_from_prior_script():
    """read the hardcoded tGCN vectors out of the prior project's script."""
    text = (RAW / "prior_tai_script_for_tRNA_counts.py.txt").read_text()
    out = {}
    for name, key in (("ecoli_trna", "ecoli"), ("bsub_trna", "bsubtilis"),
                      ("scer_trna", "scerevisiae")):
        m = re.search(rf"{name}\s*=\s*\[(.*?)\]", text, re.S)
        nums = [int(x) for x in re.findall(r"-?\d+", re.sub(r"#[^\n]*", "", m.group(1)))]
        assert len(nums) == 64, f"{name}: got {len(nums)} entries"
        out[key] = nums
    return out


def delta_per_aa(vals):
    """mean pairwise distance among synonyms given {codon: coordinate-vector}."""
    out = {}
    for aa in sorted({CODON_AA[c] for c in vals}):
        syn = [c for c in vals if CODON_AA[c] == aa]
        if len(syn) < 2:
            continue
        d = [np.linalg.norm(np.asarray(vals[x]) - np.asarray(vals[y]))
             for i, x in enumerate(syn) for y in syn[i + 1:]]
        out[aa] = float(np.mean(d))
    return out


def standardize(m):
    v = np.array(list(m.values()), float)
    v = (v - v.mean()) / v.std(ddof=0)
    return dict(zip(m.keys(), v))


def part_b():
    mu = pd.read_csv(RAW / "codon_error_rates_ecoli.tsv", sep="\t")
    mu = {r.codon: r.mu for r in mu.itertuples() if r.mu > 0}

    used = {
        "ecoli": RAW / "ecoli_tai_CORRUPT_for_regression_test.tsv",
        "bsubtilis": RAW / "bsubtilis_tai_AS_USED_previously.tsv",
        "scerevisiae": RAW / "scerevisiae_tai_AS_USED_previously.tsv",
    }
    tai_used = {}
    for k, p in used.items():
        t = pd.read_csv(p, sep="\t")
        tai_used[k] = dict(zip(t[t.columns[0]], t[t.columns[1]]))

    published = pd.read_csv(RAW / "delta_a_by_species_AS_PUBLISHED.tsv",
                            sep="\t").set_index("aa")

    # (i) reproduce the published Delta_A using E. coli mu for ALL species
    max_diff = 0.0
    for sp in ("ecoli", "bsubtilis", "scerevisiae"):
        codons = [c for c in SENSE if c in mu and c in tai_used[sp]]
        mz = standardize({c: np.log(mu[c]) for c in codons})
        nz = standardize({c: tai_used[sp][c] for c in codons})
        d = delta_per_aa({c: [mz[c], nz[c]] for c in codons})
        s = pd.Series(d)
        col = f"delta_{sp}"
        max_diff = max(max_diff,
                       float((s - published[col].reindex(s.index)).abs().max()))

    # (ii) how independent are the three "species-specific" vectors
    ec, sc, bs = tai_used["ecoli"], tai_used["scerevisiae"], tai_used["bsubtilis"]
    shared = [c for c in SENSE if c in ec and c in sc]
    identical_ec_sc = sum(1 for c in shared if ec[c] == sc[c])
    shared_b = [c for c in SENSE if c in ec and c in bs]
    identical_ec_bs = sum(1 for c in shared_b if ec[c] == bs[c])

    # (iii) indicative recomputation: nu-only Delta_A, independent canonical tAI
    counts = tRNA_vectors_from_prior_script()
    tai_canon = {"ecoli": canonical_tai(counts["ecoli"], 1),
                 "bsubtilis": canonical_tai(counts["bsubtilis"], 1),
                 "scerevisiae": canonical_tai(counts["scerevisiae"], 0)}

    def rho_table(table):
        deltas = {}
        for sp, t in table.items():
            codons = [c for c in SENSE if c in t and c in mu]
            nz = standardize({c: t[c] for c in codons})
            deltas[sp] = delta_per_aa({c: [nz[c]] for c in codons})
        aas = sorted(set.intersection(*[set(v) for v in deltas.values()]))
        out = {}
        for a, b in (("ecoli", "scerevisiae"), ("ecoli", "bsubtilis"),
                     ("bsubtilis", "scerevisiae")):
            r, p = stats.spearmanr([deltas[a][x] for x in aas],
                                   [deltas[b][x] for x in aas])
            out[f"{a}_vs_{b}"] = {"rho": round(float(r), 3), "p": float(f"{p:.3g}")}
        return out

    rho_used = rho_table(tai_used)
    rho_canon = rho_table(tai_canon)

    summary = {
        "published_delta_A_reproduced_using_ecoli_mu_for_all_species": {
            "max_abs_difference": max_diff,
            "exact": bool(max_diff < 1e-9),
            "meaning": "Delta_A is the 2D (mu, nu) distance and mu is E. coli's "
                       "for all three species, so the correlated vectors share "
                       "an identical coordinate by construction",
        },
        "tai_vector_independence": {
            "ecoli_vs_scerevisiae_bit_identical": identical_ec_sc,
            "ecoli_vs_bsubtilis_bit_identical": identical_ec_bs,
            "n_compared": len(shared),
        },
        "published_rho_2D": {
            k: round(float(stats.spearmanr(published[f"delta_{a}"],
                                           published[f"delta_{b}"])[0]), 3)
            for k, (a, b) in {
                "ecoli_vs_scerevisiae": ("ecoli", "scerevisiae"),
                "ecoli_vs_bsubtilis": ("ecoli", "bsubtilis"),
                "bsubtilis_vs_scerevisiae": ("bsubtilis", "scerevisiae")}.items()
        },
        "rho_nu_only_vectors_as_used": rho_used,
        "rho_nu_only_independent_canonical_tai": rho_canon,
        "caveat": "the canonical recomputation uses the prior script's own "
                  "hardcoded tGCN vectors, which are unverified (its E. coli "
                  "vector sums to 116 genes where K-12 has ~86). indicative "
                  "only; the shared-mu artifact and the identity counts are exact.",
    }
    (COMP / "removed_crossspecies_test.json").write_text(json.dumps(summary, indent=2))

    print("\n" + "=" * 92)
    print("PART B -- cross-species conservation of operational geometry")
    print("=" * 92)
    print(f"  published Delta_A reproduced with E. coli mu for all 3 species: "
          f"max|diff| = {max_diff:.2e}  -> {'EXACT' if max_diff < 1e-9 else 'no'}")
    print(f"  'species-specific' tAI vectors, bit-identical codons:")
    print(f"      E. coli vs S. cerevisiae : {identical_ec_sc} / {len(shared)}")
    print(f"      E. coli vs B. subtilis   : {identical_ec_bs} / {len(shared_b)}")
    print(f"\n  Spearman rho of Delta_A across species:")
    print(f"      published (2D, shared mu)      : "
          + ", ".join(f"{k.split('_vs_')[0][:2]}-{k.split('_vs_')[1][:2]}={v}"
                      for k, v in summary["published_rho_2D"].items()))
    print(f"      nu-only, vectors as used       : "
          + ", ".join(f"{k.split('_vs_')[0][:2]}-{k.split('_vs_')[1][:2]}="
                      f"{v['rho']}" for k, v in rho_used.items()))
    print(f"      nu-only, independent canonical : "
          + ", ".join(f"{k.split('_vs_')[0][:2]}-{k.split('_vs_')[1][:2]}="
                      f"{v['rho']}{'' if v['p'] < .05 else ' (ns)'}"
                      for k, v in rho_canon.items()))
    return summary


def main():
    COMP.mkdir(parents=True, exist_ok=True)
    diag, _ = part_a()
    part_b()
    (COMP / "removed_metal_site_diagnostics.json").write_text(json.dumps(diag, indent=2))
    print(f"\nwrote diagnostics to {COMP}")


if __name__ == "__main__":
    main()
