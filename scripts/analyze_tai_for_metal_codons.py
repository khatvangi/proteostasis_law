#!/usr/bin/env python
import argparse
import csv
from collections import defaultdict
from math import isfinite

from scipy.stats import fisher_exact


def load_tai(tai_tsv):
    codon2tai = {}
    with open(tai_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            codon = row["codon"].strip().upper()
            w = row["w_tai"]
            try:
                w = float(w)
            except ValueError:
                continue
            codon2tai[codon] = w
    return codon2tai


def load_summary(summary_csv):
    aa_map = defaultdict(list)
    with open(summary_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            aa = row["aa"].strip().upper()
            codon = row["codon"].strip().upper()
            lig = int(row["ligand_count"])
            bg = int(row["background_count"])
            aa_map[aa].append(
                {
                    "codon": codon,
                    "ligand": lig,
                    "background": bg,
                }
            )
    return aa_map


def main():
    ap = argparse.ArgumentParser(
        description="Analyze tAI for metal-binding codons (ligand vs background)."
    )
    ap.add_argument(
        "--summary_csv",
        required=True,
        help="Input CSV with aa,codon,ligand_count,background_count.",
    )
    ap.add_argument(
        "--tai_tsv",
        required=True,
        help="TSV with codon and w_tai (from tAI).",
    )
    ap.add_argument(
        "--out_tsv",
        required=True,
        help="Output TSV with per-AA tAI summary.",
    )
    args = ap.parse_args()

    codon2tai = load_tai(args.tai_tsv)
    print(f"[INFO] Loaded tAI weights for {len(codon2tai)} codons")

    aa_map = load_summary(args.summary_csv)
    print(f"[INFO] Loaded codon counts for {len(aa_map)} amino acids")

    out_fields = [
        "aa",
        "codon1",
        "codon2",
        "ligand_codon1",
        "ligand_codon2",
        "background_codon1",
        "background_codon2",
        "tai_codon1",
        "tai_codon2",
        "tai_ligand_weighted",
        "tai_background_weighted",
        "OR_enrichment",
        "p_fisher",
        "enriched_codon",
        "enriched_tai",
        "depleted_codon",
        "depleted_tai",
        "delta_tai_enriched_minus_depleted",
        "delta_tai_lig_minus_bg",
    ]

    with open(args.out_tsv, "w") as out_f:
        writer = csv.DictWriter(out_f, fieldnames=out_fields, delimiter="\t")
        writer.writeheader()

        for aa in sorted(aa_map.keys()):
            rows = aa_map[aa]
            if len(rows) != 2:
                print(f"[WARN] Skipping aa={aa}, expected 2 codons, found {len(rows)}")
                continue

            r1, r2 = rows[0], rows[1]
            c1, c2 = r1["codon"], r2["codon"]
            lig1, lig2 = r1["ligand"], r2["ligand"]
            bg1, bg2 = r1["background"], r2["background"]

            t1 = codon2tai.get(c1)
            t2 = codon2tai.get(c2)
            if t1 is None or t2 is None:
                print(f"[WARN] Missing tAI for aa={aa}, codons={c1},{c2}")
                continue

            total_lig = lig1 + lig2
            total_bg = bg1 + bg2

            if total_lig > 0:
                tai_lig = (lig1 * t1 + lig2 * t2) / total_lig
            else:
                tai_lig = float("nan")

            if total_bg > 0:
                tai_bg = (bg1 * t1 + bg2 * t2) / total_bg
            else:
                tai_bg = float("nan")

            # 2×2 contingency table for enrichment (same as before)
            table = [[lig1, lig2], [bg1, bg2]]
            try:
                OR, p = fisher_exact(table)
            except Exception:
                OR, p = float("nan"), float("nan")

            # Decide enriched/depleted by OR
            if isfinite(OR) and OR > 1.0:
                enriched_codon, enriched_tai = c1, t1
                depleted_codon, depleted_tai = c2, t2
            elif isfinite(OR) and OR < 1.0:
                enriched_codon, enriched_tai = c2, t2
                depleted_codon, depleted_tai = c1, t1
            else:
                enriched_codon, enriched_tai = "", float("nan")
                depleted_codon, depleted_tai = "", float("nan")

            delta_enriched_minus_depleted = (
                enriched_tai - depleted_tai
                if (isfinite(enriched_tai) and isfinite(depleted_tai))
                else float("nan")
            )
            delta_lig_minus_bg = (
                tai_lig - tai_bg
                if (isfinite(tai_lig) and isfinite(tai_bg))
                else float("nan")
            )

            writer.writerow(
                {
                    "aa": aa,
                    "codon1": c1,
                    "codon2": c2,
                    "ligand_codon1": lig1,
                    "ligand_codon2": lig2,
                    "background_codon1": bg1,
                    "background_codon2": bg2,
                    "tai_codon1": t1,
                    "tai_codon2": t2,
                    "tai_ligand_weighted": tai_lig,
                    "tai_background_weighted": tai_bg,
                    "OR_enrichment": OR,
                    "p_fisher": p,
                    "enriched_codon": enriched_codon,
                    "enriched_tai": enriched_tai,
                    "depleted_codon": depleted_codon,
                    "depleted_tai": depleted_tai,
                    "delta_tai_enriched_minus_depleted": delta_enriched_minus_depleted,
                    "delta_tai_lig_minus_bg": delta_lig_minus_bg,
                }
            )

    print(f"[INFO] Wrote tAI summary to {args.out_tsv}")


if __name__ == "__main__":
    main()

