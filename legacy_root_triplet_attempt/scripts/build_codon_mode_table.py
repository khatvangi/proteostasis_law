#!/usr/bin/env python

import argparse
import csv
from collections import defaultdict
import math

def load_tsv_map(path, key_col, val_cols):
    """
    Load a TSV with a key column and a set of value columns into a dict.
    key_col: name of column to use as key
    val_cols: list of column names to keep
    Returns: dict[key] = {col: value, ...}
    """
    out = {}
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            key = row[key_col].strip()
            if not key:
                continue
            rec = {}
            for col in val_cols:
                val = row[col].strip()
                rec[col] = val
            out[key] = rec
    return out

def median(values):
    vals = sorted(v for v in values if v is not None and not math.isnan(v))
    if not vals:
        return float("nan")
    n = len(vals)
    mid = n // 2
    if n % 2 == 1:
        return vals[mid]
    else:
        return 0.5 * (vals[mid-1] + vals[mid])

def main():
    ap = argparse.ArgumentParser(
        description="Join global codon usage with mu and tAI; define modes per AA"
    )
    ap.add_argument(
        "--usage_tsv",
        required=True,
        help="errors/global_codon_usage.tsv"
    )
    ap.add_argument(
        "--mu_tsv",
        required=True,
        help="errors/codon_error_rates.tsv (codon, mu)"
    )
    ap.add_argument(
        "--tai_tsv",
        required=True,
        help="errors/ecoli_tai_ws.tsv (codon, w_tai)"
    )
    ap.add_argument(
        "--out_tsv",
        required=True,
        help="Output TSV with modes"
    )
    args = ap.parse_args()

    # Load mu and tAI maps
    mu_map = load_tsv_map(args.mu_tsv, "codon", ["mu"])
    tai_map = load_tsv_map(args.tai_tsv, "codon", ["w_tai"])

    # Load usage and attach mu/tai
    records = []
    by_aa = defaultdict(list)

    with open(args.usage_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            aa = row["aa"]
            codon = row["codon"]
            count = int(row["count"])
            freq_within = float(row["freq_within_aa"])
            freq_global = float(row["freq_global"])

            mu = None
            if codon in mu_map:
                try:
                    mu = float(mu_map[codon]["mu"])
                except ValueError:
                    mu = None

            tai = None
            if codon in tai_map:
                try:
                    tai = float(tai_map[codon]["w_tai"])
                except ValueError:
                    tai = None

            rec = {
                "aa": aa,
                "codon": codon,
                "count": count,
                "freq_within_aa": freq_within,
                "freq_global": freq_global,
                "mu": mu,
                "w_tai": tai,
            }
            records.append(rec)
            by_aa[aa].append(rec)

    # Compute per-AA medians for mu and tAI
    aa_medians = {}
    for aa, recs in by_aa.items():
        mu_vals = [r["mu"] for r in recs if r["mu"] is not None]
        tai_vals = [r["w_tai"] for r in recs if r["w_tai"] is not None]
        aa_medians[aa] = {
            "mu_med": median(mu_vals),
            "tai_med": median(tai_vals),
        }

    # Assign classes
    for rec in records:
        aa = rec["aa"]
        mu = rec["mu"]
        tai = rec["w_tai"]
        mu_med = aa_medians[aa]["mu_med"]
        tai_med = aa_medians[aa]["tai_med"]

        # mu class
        if mu is None or math.isnan(mu_med):
            rec["mu_class"] = "unknown"
        else:
            rec["mu_class"] = "low_mu" if mu <= mu_med else "high_mu"

        # tAI class
        if tai is None or math.isnan(tai_med):
            rec["tai_class"] = "unknown"
        else:
            rec["tai_class"] = "high_tai" if tai >= tai_med else "low_tai"

        # combined mode
        if rec["mu_class"] == "unknown" or rec["tai_class"] == "unknown":
            rec["mode"] = "unknown"
        else:
            if rec["mu_class"] == "low_mu" and rec["tai_class"] == "high_tai":
                rec["mode"] = "safe_sprinter"
            elif rec["mu_class"] == "low_mu" and rec["tai_class"] == "low_tai":
                rec["mode"] = "safe_careful"
            elif rec["mu_class"] == "high_mu" and rec["tai_class"] == "high_tai":
                rec["mode"] = "risky_sprinter"
            elif rec["mu_class"] == "high_mu" and rec["tai_class"] == "low_tai":
                rec["mode"] = "risky_careful"
            else:
                rec["mode"] = "unknown"

    # Write output
    with open(args.out_tsv, "w") as out:
        cols = [
            "aa", "codon", "count",
            "freq_within_aa", "freq_global",
            "mu", "w_tai",
            "mu_class", "tai_class", "mode",
        ]
        out.write("\t".join(cols) + "\n")
        for rec in sorted(records, key=lambda r: (r["aa"], r["codon"])):
            vals = []
            for c in cols:
                v = rec[c]
                if isinstance(v, float):
                    vals.append(f"{v:.8g}")
                else:
                    vals.append(str(v))
            out.write("\t".join(vals) + "\n")

    print(f"[INFO] Wrote codon mode table to {args.out_tsv}")

if __name__ == "__main__":
    main()

