"""
Landerer et al. 2024 (MBE 41:msae048) per-dataset data, read-only.

source: the authors' supplementary archive msae048_supplementary_data.zip, which
contains Data_S2 (per-codon summary) and Data_S4_detected_misincorporations.zip
(per-dataset files). for every E. coli dataset (80) S4 has
  <PXD>_codon_counts.csv : codon, base_count (PSMs covering a position with this
                           codon), error_count (of those, PSMs carrying a
                           substitution there), detection_rate = error/base
  <PXD>_substitution_errors.csv : one row per substitution PSM
so the per-dataset numerators AND denominators, including covered-but-zero
cells, are available. nothing here is inferred from (sd/se)^2.

three cell types are kept apart:
  no coverage     codon row absent (or base_count == 0) -> no information
  covered, zero   base_count > 0, error_count == 0      -> an observed rate of 0
  covered, >0     base_count > 0, error_count  > 0

upstream mechanism (deTEL, git.mpi-cbg.de/tothpetroczylab/detel, main @ c3593d46,
eTEL/workflow/global_report.py:get_all_codon_count): detection_rate is mapped to
log10, -inf (rate 0) is replaced by NaN, and those rows are dropped. Data_S2
is reproduced exactly only under that drop (see reconstruct_s2).
"""
import hashlib
import io
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

SUPPL_ZIP = Path("/storage/kiran-stuff/triplet-proof/reviewer_response/position_rates/"
                 "raw/landerer2024/msae048_supplementary_data.zip")
USAGE_TSV = Path("/storage/kiran-stuff/proteostasis_law/envelope-paper/data/raw/"
                 "global_codon_usage_ecoli.tsv")
SHEET = {"ecoli": "E. coli", "yeast": "S. cerevisiae"}
STOPS = {"TAA", "TAG", "TGA"}


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def sha256_member(name):
    """sha256 of one file inside the authors' supplementary zip."""
    with zipfile.ZipFile(SUPPL_ZIP) as z:
        return hashlib.sha256(z.read(name)).hexdigest()


def _outer():
    return zipfile.ZipFile(SUPPL_ZIP)


def load_s2(org="ecoli"):
    with _outer() as z:
        return pd.read_excel(io.BytesIO(z.read("Data_S2_error_detection_rate.xlsx")),
                             sheet_name=SHEET[org])


def load_s4(org="ecoli"):
    """per-dataset tables: (codon_counts, substitution_errors), both with a 'ds' column."""
    with _outer() as z:
        inner = zipfile.ZipFile(io.BytesIO(z.read("Data_S4_detected_misincorporations.zip")))
    cc, se = [], []
    for n in sorted(inner.namelist()):
        if not n.startswith(f"{org}_filtered/") or not n.endswith(".csv"):
            continue
        ds = n.split("/")[-1].split("_")[0]
        d = pd.read_csv(io.BytesIO(inner.read(n)))
        d["ds"] = ds
        if n.endswith("_codon_counts.csv"):
            cc.append(d.drop(columns=[c for c in d.columns if c.startswith("Unnamed")]))
        elif n.endswith("_substitution_errors.csv"):
            se.append(d)
    return pd.concat(cc, ignore_index=True), pd.concat(se, ignore_index=True)


def cell_census(cc, n_datasets):
    """count no-coverage / covered-zero / covered-positive sense-codon cells."""
    sense = cc[~cc.codon.isin(STOPS)]
    covered = sense[sense.base_count > 0]
    return {
        "datasets": int(cc.ds.nunique()),
        "possible_sense_cells": 61 * n_datasets,
        "no_coverage_cells": int(61 * n_datasets - len(covered)),
        "covered_zero_cells": int((covered.error_count == 0).sum()),
        "covered_positive_cells": int((covered.error_count > 0).sum()),
        "stop_codon_rows": int(cc.codon.isin(STOPS).sum()),
        "max_abs_rate_minus_err_over_base": float(np.max(np.abs(
            covered.detection_rate - covered.error_count / covered.base_count))),
    }


def numerator_check(cc, se):
    """error_count in codon_counts equals the substitution-PSM rows per (ds, codon)."""
    k = se.groupby(["ds", "codon"]).size().rename("n_rows")
    j = cc.set_index(["ds", "codon"]).join(k).fillna({"n_rows": 0})
    return int((j.error_count != j.n_rows).sum())


def reconstruct_s2(cc, s2):
    """compare Data_S2 against two aggregation rules. returns max abs deviations."""
    sense = cc[~cc.codon.isin(STOPS) & (cc.base_count > 0)]
    s2 = s2.set_index("Codon")
    n_s2 = ((s2.sd / s2.se) ** 2).round().astype(int)
    out = {}
    for label, sub in (("zero_inclusive_covering", sense),
                       ("detected_only_error_count_gt_0", sense[sense.error_count > 0])):
        g = sub.groupby("codon").detection_rate.agg(["mean", "std", "median", "count"])
        g = g.reindex(s2.index)
        out[label] = {
            "max_abs_dev_mean": float(np.nanmax(np.abs(g["mean"] - s2["mean"]))),
            "max_abs_dev_sd": float(np.nanmax(np.abs(g["std"] - s2["sd"]))),
            "max_abs_dev_median": float(np.nanmax(np.abs(g["median"] - s2["median"]))),
            "n_equals_sd_over_se_sq": bool((g["count"] == n_s2).all()),
        }
    return out


def usage_weights():
    """the legacy weights: E. coli genome codon counts (scripts/06), 61 sense codons."""
    u = pd.read_csv(USAGE_TSV, sep="\t")
    u = u[~u.codon.isin(STOPS)]
    return (u.set_index("codon")["count"] / u["count"].sum())


def per_codon_rates(cc):
    """three E. coli per-codon eTEL statistics (all MS-detected, substitution level)."""
    sense = cc[~cc.codon.isin(STOPS) & (cc.base_count > 0)]
    g = sense.groupby("codon")
    return pd.DataFrame({
        # Data_S2 'mean': mean of per-dataset rates over datasets with >= 1 detection
        "conditional_mean": sense[sense.error_count > 0].groupby("codon").detection_rate.mean(),
        # mean of per-dataset rates over every dataset covering the codon (zeros kept)
        "zero_inclusive_mean": g.detection_rate.mean(),
        # PSM-pooled ratio: sum of errors / sum of covering PSMs, across datasets
        "pooled_ratio": g.error_count.sum() / g.base_count.sum(),
        "n_covering": g.size(),
        "n_detected": sense[sense.error_count > 0].groupby("codon").size(),
        "sum_base": g.base_count.sum(),
        "sum_err": g.error_count.sum(),
    }).fillna({"n_detected": 0})


def usage_weighted(series, w):
    s = series.reindex(w.index)
    assert s.notna().all(), "a sense codon lacks a rate"
    return float((w * s).sum())


def per_dataset_usage_weighted(cc, w):
    """each dataset's own usage-weighted rate over the codons it covers (weights
    renormalised over covered codons), plus its PSM-pooled global rate."""
    sense = cc[~cc.codon.isin(STOPS) & (cc.base_count > 0)]
    rows = []
    for ds, d in sense.groupby("ds"):
        d = d.set_index("codon")
        wd = w.reindex(d.index)
        rows.append({"ds": ds, "usage_weighted": float((wd * d.detection_rate).sum() / wd.sum()),
                     "pooled": float(d.error_count.sum() / d.base_count.sum()),
                     "psm": int(d.base_count.sum()), "n_codons": len(d)})
    return pd.DataFrame(rows)
