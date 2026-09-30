"""WO-07 derived inputs: synthesis weighting of chain length and codon usage.

why this exists: WO-06 flagged two quantity mismatches that are fixable from
files already on disk, without any new source.
  * N: the legacy used genome-encoded lengths (every gene counted once). the
    burden-relevant length is weighted by how many chains are MADE.
  * codon usage weights for aggregating per-codon error rates: the legacy used
    genomic counts (E05). the relevant weights count codons TRANSLATED.

in balanced exponential growth with negligible turnover, a protein's synthesis
rate is mu * copies, so copy number is a synthesis weight. this closure holds
ONLY in balanced growth: in stationary phase copies/cell are an abundance, not
a synthesis rate, and the stationary weights below are labelled as such.

inputs (read-only):
  Schmidt 2016 Table S6 copies/cell (WO-06/records, sha256 checked by WO-06)
  E. coli K-12 MG1655 CDS (NC_000913.3), mapped by UniProt accession
  the legacy UniProt length table (hashed in WO-00)
  Landerer 2024 Data_S4 per-dataset counts (via WO-05/landerer_s4.py)
"""
import hashlib
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "WO-06"))
sys.path.insert(0, str(ROOT / "WO-05"))
import schmidt_pools as sp  # noqa: E402
import landerer_s4 as L  # noqa: E402

CDS = Path("/storage/kiran-stuff/proteostasis_law/envelope-paper/data/raw/ecoli_k12_cds.fna")
CDS_SHA256 = "9bd477a9ecd7b84d00fe6c5a90f4f1a7cf6066c43f8e943dae9eb45fef78fbb5"
GENOME_LEN = Path("/storage/kiran-stuff/proteostasis-P1/ecoli_proteome_lengths.tsv")
GENOME_LEN_SHA256 = "8c472011e909f3f477b2b50d75f8451e02b0e9adf40b2275f0a81e206691bba7"
STOPS = {"TAA", "TAG", "TGA"}
EXP_COND, STAT_COND = "Glucose", "Stationary phase 1 day"


def sha256(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def read_cds():
    """accession -> in-frame sense-codon list (terminal stop removed)."""
    assert sha256(CDS) == CDS_SHA256, "CDS file changed"
    out, acc, seq = {}, None, []
    for line in CDS.read_text().splitlines() + [">"]:
        if line.startswith(">"):
            if acc and seq:
                s = "".join(seq).upper()
                cod = [s[i:i + 3] for i in range(0, len(s) - len(s) % 3, 3)]
                if cod and cod[-1] in STOPS:
                    cod = cod[:-1]
                out.setdefault(acc, cod)          # first CDS per accession
            m = re.search(r"UniProtKB/Swiss-Prot:(\w+)", line)
            acc, seq = (m.group(1) if m else None), []
        else:
            seq.append(line.strip())
    return out


def weighted_lengths(cond):
    """(lengths in codons, copy weights) for Schmidt proteins that map to a CDS."""
    copies, _, _, _ = sp.load()
    cds = read_cds()
    c = copies[cond].dropna()
    c = c[c > 0]
    keep = [a for a in c.index if a in cds]
    L_ = np.array([len(cds[a]) for a in keep], float)
    w = c.loc[keep].to_numpy(float)
    return L_, w, {"n_mapped": len(keep), "n_quantified": int(len(c)),
                   "copy_fraction_mapped": float(w.sum() / c.sum())}


def genome_lengths():
    assert sha256(GENOME_LEN) == GENOME_LEN_SHA256, "legacy length table changed"
    L_ = pd.read_csv(GENOME_LEN, sep="\t")["Length"].to_numpy(float)
    return L_, np.ones_like(L_)


def translated_codon_usage(cond):
    """codon weights = sum over proteins of copies x codon counts (61 sense codons)."""
    copies, _, _, _ = sp.load()
    cds = read_cds()
    c = copies[cond].dropna()
    tot = {}
    for a, n in c.items():
        if a in cds and n > 0:
            for k, v in pd.Series(cds[a]).value_counts().items():
                tot[k] = tot.get(k, 0.0) + n * v
    w = pd.Series(tot)
    w = w[[k for k in w.index if len(k) == 3 and set(k) <= set("ACGT") and k not in STOPS]]
    return w / w.sum()


def etel_aggregates(w):
    """the four WO-05 eTEL aggregates under codon weights w (MS, substitution level)."""
    cc, _ = L.load_s4("ecoli")
    r = L.per_codon_rates(cc)
    ds = L.per_dataset_usage_weighted(cc, w)
    return {"conditional_mean_DataS2": L.usage_weighted(r.conditional_mean, w),
            "zero_inclusive_mean": L.usage_weighted(r.zero_inclusive_mean, w),
            "pooled_ratio": L.usage_weighted(r.pooled_ratio, w),
            "dataset_median": float(ds.usage_weighted.median())}


def summarize(L_, w):
    return {"weighted_mean_codons": float((L_ * w).sum() / w.sum()),
            "unweighted_median_codons": float(np.median(L_))}


def run():
    Le, we, me = weighted_lengths(EXP_COND)
    Ls, ws, ms = weighted_lengths(STAT_COND)
    Lg, wg = genome_lengths()
    w_trans = translated_codon_usage(EXP_COND)
    w_gen = L.usage_weights()
    out = {
        "inputs_sha256": {"cds": CDS_SHA256, "genome_lengths": GENOME_LEN_SHA256,
                          "schmidt_xlsx": sp.XLSX_SHA256},
        "N_synthesis_weighted_exp_glucose": {**summarize(Le, we), **me,
            "schmidt_condition": EXP_COND,
            "closure": "copies = synthesis weight only in balanced growth with negligible turnover"},
        "N_abundance_weighted_stationary_1d": {**summarize(Ls, ws), **ms,
            "schmidt_condition": STAT_COND,
            "closure": "copy (standing-proteome) weight; NOT a synthesis weight in non-growing cells"},
        "N_genome_unweighted": {**summarize(Lg, wg), "n": int(len(Lg)),
                                "source": "UniProt UP000000625 length table (legacy, WO-06 L08b)"},
        "codon_usage_translated_exp_glucose": {"schmidt_condition": EXP_COND,
                                               "weights": w_trans.round(10).to_dict()},
        "etel_genomic_weights": etel_aggregates(w_gen),
        "etel_translated_weights_exp_glucose": etel_aggregates(w_trans),
    }
    return out


if __name__ == "__main__":
    r = run()
    (HERE / "derived.json").write_text(json.dumps(r, indent=1))
    print(json.dumps({k: v for k, v in r.items() if k != "codon_usage_translated_exp_glucose"},
                     indent=1))
