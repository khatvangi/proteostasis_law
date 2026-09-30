"""WO-06: every number in parameter_audit.tsv that is COMPUTED (not quoted) comes
from here and is written to audit_derived.json. test_wo06 re-reads this file.

checks:
  codon_usage   do the legacy genome codon counts reproduce from the local
                NC_000913.3 CDS file they were presumably built from?
  S_code        what the standard genetic code implies for the synonymous share
                of single-nucleotide misreads, if misreading were uniform. this
                is a property of the code, not a measurement of decoding errors.
  k_deg         half-lives implied by the legacy k_deg range vs the "1-10 hr"
                consensus the legacy cites for it.
  prot_tot      legacy 300 uM vs three independent records.
  lengths       median / mean of the legacy's own UniProt length file.
  k_clear       the legacy's own arithmetic -ln(0.05)/7200 s.
  pools         legacy C_tot vs measured DnaK and GroEL pools (schmidt_pools.json).
"""
import json
import math
from collections import Counter
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
EP = Path("/storage/kiran-stuff/proteostasis_law/envelope-paper/data/raw")
P1 = Path("/storage/kiran-stuff/proteostasis-P1")
N_A = 6.02214076e23

BASES = "TCAG"
AA = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
CODE = {a + b + c: AA[16 * i + 4 * j + k]
        for i, a in enumerate(BASES) for j, b in enumerate(BASES) for k, c in enumerate(BASES)}


def read_cds(path):
    seqs, cur = [], []
    for line in open(path):
        if line.startswith(">"):
            if cur:
                seqs.append("".join(cur))
            cur = []
        else:
            cur.append(line.strip().upper())
    if cur:
        seqs.append("".join(cur))
    return seqs


def codon_usage_check():
    usage = pd.read_csv(EP / "global_codon_usage_ecoli.tsv", sep="\t")
    cnt = Counter()
    seqs = read_cds(EP / "ecoli_k12_cds.fna")
    for s in seqs:
        if len(s) % 3:
            continue
        for i in range(0, len(s) - 3, 3):          # drop the stop codon
            cnt[s[i:i + 3]] += 1
    diff = {r.codon: int(r["count"]) - cnt.get(r.codon, 0) for _, r in usage.iterrows()}
    return {"n_cds_in_fna": len(seqs), "n_codons_in_tsv": int(len(usage)),
            "tsv_total": int(usage["count"].sum()),
            "fna_total_sense_excl_stop": int(sum(v for k, v in cnt.items() if CODE.get(k, "*") != "*")),
            "max_abs_count_diff": int(max(abs(v) for v in diff.values())),
            "n_codons_differing": int(sum(v != 0 for v in diff.values()))}


def synonymous_share(weights=None, positions=(0, 1, 2)):
    """share of single-nucleotide substitutions of sense codons that stay
    synonymous, excluding changes to stop. uniform over positions/bases."""
    num = den = 0.0
    for c, aa in CODE.items():
        if aa == "*":
            continue
        w = 1.0 if weights is None else weights.get(c, 0.0)
        for p in positions:
            for b in BASES:
                if b == c[p]:
                    continue
                m = c[:p] + b + c[p + 1:]
                if CODE[m] == "*":
                    continue
                den += w
                num += w * (CODE[m] == aa)
    return num / den


def run():
    usage = pd.read_csv(EP / "global_codon_usage_ecoli.tsv", sep="\t")
    wu = dict(zip(usage.codon, usage["count"].astype(float)))
    S = {"all_positions_unweighted": synonymous_share(),
         "all_positions_usage_weighted": synonymous_share(wu),
         "third_position_unweighted": synonymous_share(positions=(2,)),
         "third_position_usage_weighted": synonymous_share(wu, positions=(2,)),
         "legacy_value": 0.30, "legacy_ranges": {"two_pool": [0.25, 0.35], "arithmetic_paired": [0.20, 0.40]}}

    kd = {"legacy_baseline_per_s": 3e-4, "legacy_range_per_s": [1e-4, 1e-3],
          "half_life_min_at_baseline": math.log(2) / 3e-4 / 60,
          "half_life_min_at_range": [math.log(2) / 1e-3 / 60, math.log(2) / 1e-4 / 60],
          "cited_consensus_hours": [1, 10],
          "k_equiv_to_cited_consensus_per_s": [math.log(2) / (10 * 3600), math.log(2) / 3600],
          "christiano_median_half_life_h_Scer": 8.8,
          "k_from_christiano_median_per_s": math.log(2) / (8.8 * 3600)}
    kd["baseline_inside_cited_consensus"] = (kd["k_equiv_to_cited_consensus_per_s"][0]
                                            <= 3e-4 <= kd["k_equiv_to_cited_consensus_per_s"][1])

    pools = json.loads((HERE / "schmidt_pools.json").read_text())
    pc = {r["condition"]: r for r in pools["per_condition"]}
    milo_mM = [n / (N_A * 1e-15) * 1e3 for n in (2e6, 4e6)]
    pt = {"legacy_uM": 300.0, "legacy_range_uM": [250.0, 350.0],
          "milo2013_mM_from_2to4e6_per_um3": milo_mM,
          "bnid104726_mM": 4.0,
          "schmidt_glucose_mM_wholecell": pc["Glucose"]["total_protein_mM_wholecell"],
          "schmidt_exponential_range_mM": [pools["summary"]["total_protein_mM_wholecell"]["exp_min"],
                                           pools["summary"]["total_protein_mM_wholecell"]["exp_max"]]}
    pt["fold_low_vs_schmidt_glucose"] = pt["schmidt_glucose_mM_wholecell"] * 1e3 / 300.0
    pt["fold_low_vs_bnid"] = 4000.0 / 300.0

    lens = pd.read_csv(P1 / "ecoli_proteome_lengths.tsv", sep="\t")["Length"]
    ln = {"n": int(len(lens)), "median": float(lens.median()), "mean": float(lens.mean()),
          "legacy_N": 300.0, "legacy_label": "median protein length"}

    kc = {"legacy_k_clear": 4e-4, "recomputed_minus_ln_0p05_over_7200": -math.log(0.05) / 7200}

    def comb(cond):
        r = pc[cond]
        return {"dnaK_plus_groL_protomers_uM": r["dnaK_uM_protomer"] + r["groL_uM_protomer"],
                "dnaK_plus_groEL14_uM": r["dnaK_uM_protomer"] + r["groL_uM_oligomer"],
                "clpB6_uM": r["clpB_uM_oligomer"],
                "dnaK_over_dnaJ": r["dnaK_uM_protomer"] / r["dnaJ_uM_protomer"]}
    ch = {c: comb(c) for c in ("Glucose", "LB", "42°C glucose", "Stationary phase 1 day",
                               "Stationary phase 3 days")}
    ch["legacy_C_tot_uM"] = 50.0

    out = {"codon_usage": codon_usage_check(), "S_code": S, "k_deg": kd, "prot_tot": pt,
           "lengths": ln, "k_clear": kc, "pools": ch}
    (HERE / "audit_derived.json").write_text(json.dumps(out, indent=1))
    return out


if __name__ == "__main__":
    print(json.dumps(run(), indent=1))
