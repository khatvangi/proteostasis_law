#!/usr/bin/env python3
"""
recompute mode assignments using GLOBAL thresholds instead of
within-family medians. this tests whether k_A = 2 for 2-codon
AAs is an empirical finding or a methodological tautology.
"""
import csv
from collections import defaultdict
import statistics

def main():
    # load the mode table (which has μ and tAI per codon)
    records = []
    with open("errors/codon_modes_ecoli_FRESH.tsv") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            mu = float(row["mu"]) if row["mu"] and row["mu"] != "None" else None
            tai = float(row["w_tai"]) if row["w_tai"] and row["w_tai"] != "None" else None
            records.append({
                "aa": row["aa"],
                "codon": row["codon"],
                "mu": mu,
                "tai": tai,
                "within_mode": row["mode"],
            })

    # compute GLOBAL medians across all sense codons
    all_mu = [r["mu"] for r in records if r["mu"] is not None]
    all_tai = [r["tai"] for r in records if r["tai"] is not None]

    global_mu_med = statistics.median(all_mu)
    global_tai_med = statistics.median(all_tai)

    print(f"Global μ median:   {global_mu_med:.2e}")
    print(f"Global tAI median: {global_tai_med:.4f}")
    print()

    # assign modes using global thresholds
    by_aa = defaultdict(list)
    for r in records:
        mu = r["mu"]
        tai = r["tai"]

        if mu is not None:
            mu_class = "low" if mu <= global_mu_med else "high"
        else:
            mu_class = "?"

        if tai is not None:
            tai_class = "high" if tai >= global_tai_med else "low"
        else:
            tai_class = "?"

        if mu_class != "?" and tai_class != "?":
            if mu_class == "low" and tai_class == "high":
                mode = "safe_sprinter"
            elif mu_class == "low" and tai_class == "low":
                mode = "safe_careful"
            elif mu_class == "high" and tai_class == "high":
                mode = "risky_sprinter"
            else:
                mode = "risky_careful"
        else:
            mode = "unknown"

        r["global_mode"] = mode
        by_aa[r["aa"]].append(r)

    # compute k_A under global thresholds
    print(f"{'AA':<5} {'Deg':<5} {'k_A(within)':<14} {'k_A(global)':<14} "
          f"{'Global modes':<45} {'Changed?'}")
    print("-" * 100)

    total_K_within = 0
    total_K_global = 0
    two_codon_K_within = 0
    two_codon_K_global = 0
    changes = []

    for aa in sorted(by_aa.keys()):
        codons = by_aa[aa]
        deg = len(codons)

        # global k_A
        global_modes = set(c["global_mode"] for c in codons if c["global_mode"] != "unknown")
        k_global = len(global_modes) if global_modes else (1 if deg == 1 else 0)

        # within-family k_A
        within_modes = set(c["within_mode"] for c in codons if c["within_mode"] != "unknown")
        k_within = len(within_modes) if within_modes else (1 if deg == 1 else 0)

        changed = "YES" if k_global != k_within else ""
        if changed:
            changes.append((aa, deg, k_within, k_global,
                            sorted(within_modes), sorted(global_modes)))

        total_K_within += k_within
        total_K_global += k_global
        if deg == 2:
            two_codon_K_within += k_within
            two_codon_K_global += k_global

        modes_str = ", ".join(sorted(global_modes)) if global_modes else "unknown"
        print(f"{aa:<5} {deg:<5} {k_within:<14} {k_global:<14} {modes_str:<45} {changed}")

    print()
    print(f"Total K (within-family): {total_K_within}")
    print(f"Total K (global):        {total_K_global}")
    print()
    print(f"Two-codon K (within):    {two_codon_K_within}")
    print(f"Two-codon K (global):    {two_codon_K_global}")
    print(f"Doublet ceiling:         16")
    print(f"Doublet impossible (within)? {'YES' if two_codon_K_within > 16 else 'NO'}")
    print(f"Doublet impossible (global)? {'YES' if two_codon_K_global > 16 else 'NO'}")
    print()

    if changes:
        print(f"CHANGED amino acids ({len(changes)}):")
        for aa, deg, kw, kg, wm, gm in changes:
            print(f"  {aa}: {deg} codons, k_A {kw} → {kg}")
            print(f"    within: {wm}")
            print(f"    global: {gm}")
    else:
        print("No changes — global and within-family agree for all AAs")

    print()
    # show the 2-codon AAs in detail
    print("=== TWO-CODON AAs: DETAILED GLOBAL MODE ASSIGNMENTS ===")
    for aa in sorted(by_aa.keys()):
        codons = by_aa[aa]
        if len(codons) != 2:
            continue
        c1, c2 = codons
        print(f"  {aa}: {c1['codon']} (μ={c1['mu']:.2e}, tAI={c1['tai']:.3f}) → {c1['global_mode']}")
        print(f"       {c2['codon']} (μ={c2['mu']:.2e}, tAI={c2['tai']:.3f}) → {c2['global_mode']}")
        same = c1['global_mode'] == c2['global_mode']
        print(f"       → {'SAME MODE (k_A=1)' if same else 'DIFFERENT MODES (k_A=2)'}")

if __name__ == "__main__":
    main()
