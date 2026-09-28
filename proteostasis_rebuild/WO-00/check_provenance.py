#!/usr/bin/env python3
"""
WO-00 provenance checks.

  python check_provenance.py              # G0.1 + G0.2 on CLAIM_REGISTER.tsv
  python check_provenance.py --hash       # G0.3 write legacy_hashes.tsv
  python check_provenance.py --verify     # re-hash and compare against manifest

paths in the register are relative to /storage/kiran-stuff. the anchor text of
each claim must occur on exactly the cited line, so a register row cannot
drift away from the file it describes.
"""
import csv
import hashlib
import sys
from pathlib import Path

BASE = Path("/storage/kiran-stuff")
HERE = Path(__file__).resolve().parent
REGISTER = HERE.parent / "CLAIM_REGISTER.tsv"
MANIFEST = HERE / "legacy_hashes.tsv"

STATUSES = {"REPRODUCED", "CORRECTED", "REJECTED", "CONDITIONAL",
            "UNVERIFIED", "OPEN"}
WOS = {f"WO-{i:02d}" for i in range(11)}

# every legacy file this run reads. hashed once at WO-00, re-verified at the end.
LEGACY = [
    "proteostasis-P1/two_pool_ode.py",
    "proteostasis-P1/LITERATURE_ANCHORS.md",
    "proteostasis-P1/arithmetic_stress_test.py",
    "proteostasis-P1/arithmetic_summary.md",
    "proteostasis-P1/paired_mc.py",
    "proteostasis-P1/paired_mc_summary.md",
    "proteostasis-P1/two_pool_summary.md",
    "proteostasis-P1/two_pool_results.json",
    "proteostasis-P1/ecoli_proteome_lengths.tsv",
    "proteostasis-P1/ecoli_proteome_with_genes.tsv",
    "proteostasis_law/envelope-paper/manuscript/MANUSCRIPT.md",
    "proteostasis_law/envelope-paper/PLAN_P1_REPAIR.md",
    "proteostasis_law/envelope-paper/scripts/vendor/two_pool_ode.py",
    "proteostasis_law/envelope-paper/scripts/06_translation_burden.py",
    "proteostasis_law/envelope-paper/scripts/09_supraadditivity.py",
    "proteostasis_law/envelope-paper/scripts/11_headroom_sensitivity.py",
    "proteostasis_law/envelope-paper/scripts/12_chaperone_availability.py",
    "proteostasis_law/envelope-paper/data/raw/codon_error_rates_ecoli.tsv",
    "proteostasis_law/envelope-paper/data/raw/global_codon_usage_ecoli.tsv",
    "proteostasis_law/envelope-paper/data/raw/Data_S2_error_detection_rate.xlsx",
    "proteostasis_law/envelope-paper/data/raw/arithmetic_results.json",
    "proteostasis_law/envelope-paper/data/computed/translation_burden.json",
    "proteostasis_law/envelope-paper/data/computed/headroom_sensitivity_summary.json",
    "proteostasis_law/envelope-paper/data/computed/supraadditivity_summary.json",
    "proteostasis_law/envelope-paper/data/computed/chaperone_availability_summary.json",
    "proteostasis_law/investigation_2026-09-26/AUDIT.md",
    "proteostasis_law/investigation_2026-09-26/METHODS.md",
    "proteostasis_law/investigation_2026-09-26/EVIDENCE_AND_FALSIFICATION.md",
    "proteostasis_law/investigation_2026-09-26/redraw.py",
    "proteostasis_law/investigation_2026-09-26/bound_identity_checks.py",
]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_register():
    with open(REGISTER, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def check_register():
    """return a list of failure strings; empty means G0.1 and G0.2 pass."""
    fails = []
    rows = read_register()
    ids = [r["claim_id"] for r in rows]
    if len(ids) != len(set(ids)):
        fails.append("duplicate claim ids")
    for r in rows:
        cid = r["claim_id"]
        path = BASE / r["legacy_file"]
        if not path.is_file():
            fails.append(f"{cid}: missing file {path}")
            continue
        lines = path.read_text(encoding="utf-8").splitlines()
        n = int(r["line"])
        if not (1 <= n <= len(lines)):
            fails.append(f"{cid}: line {n} out of range")
        elif r["anchor_text"] not in lines[n - 1]:
            fails.append(f"{cid}: anchor not on line {n} of {r['legacy_file']}")
        if r["status"] not in STATUSES:
            fails.append(f"{cid}: bad status {r['status']!r}")
        if r["owner_wo"] not in WOS:
            fails.append(f"{cid}: bad owner {r['owner_wo']!r}")
    return fails, len(rows)


def write_manifest():
    with open(MANIFEST, "w") as f:
        f.write("path\tsha256\tbytes\n")
        for rel in LEGACY:
            p = BASE / rel
            f.write(f"{rel}\t{sha256(p)}\t{p.stat().st_size}\n")
    print(f"wrote {MANIFEST} ({len(LEGACY)} files)")


def verify_manifest():
    """return list of files whose hash changed since WO-00."""
    changed = []
    with open(MANIFEST) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if sha256(BASE / row["path"]) != row["sha256"]:
                changed.append(row["path"])
    return changed


if __name__ == "__main__":
    if "--hash" in sys.argv:
        write_manifest()
    elif "--verify" in sys.argv:
        ch = verify_manifest()
        print("legacy files unchanged" if not ch else f"CHANGED: {ch}")
        sys.exit(1 if ch else 0)
    else:
        fails, n = check_register()
        for x in fails:
            print("FAIL", x)
        print(f"{n} claims checked, {len(fails)} failures")
        sys.exit(1 if fails else 0)
