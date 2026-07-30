#!/usr/bin/env python3
"""
regenerate everything: analyses, figures, tables. then run the test suite.

    python scripts/run_all.py              # full rebuild (~5 min, dominated by permutations)
    python scripts/run_all.py --fast       # skip the two 10,000-permutation runs
    python scripts/run_all.py --no-tests   # rebuild only

the pipeline is ordered and fails closed: if 01_validate_tai.py rejects the
supply axis it exits non-zero and nothing downstream runs.
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parent
ROOT = SCRIPTS.parent

# (script, args, label, slow?)
STEPS = [
    ("01_validate_tai.py", [], "validate the nu (tAI) axis", False),
    ("02_axis_structure.py", [], "axis structure, mu = per-codon mean", True),
    ("02_axis_structure.py", ["--mu-stat", "median"],
     "axis structure, mu = per-codon median (sensitivity)", True),
    ("06_translation_burden.py", [], "translation burden magnitude", False),
    ("07_removed_results.py", [], "diagnostics for the two removed results", False),
    ("03_fig1_envelope.py", [], "Fig 1  envelope regimes", False),
    ("04_fig2_axis.py", [], "Fig 2  axis structure", False),

    ("09_supraadditivity.py", [], "supraadditivity test (2x2 on the two-pool ODE)", False),
    ("11_headroom_sensitivity.py", [], "headroom sensitivity (evaluation point x anchoring)", False),
    ("12_chaperone_availability.py", [], "chaperone availability sweep (theta)", False),
    ("13_nu_power.py", [], "axis power: minimum detectable effect + power curve", True),
    ("14_mu_jackknife.py", [], "mu jackknife and detectability exposure", True),
    ("05_fig3_bounds.py", [], "Fig 3  bounds and headroom sensitivity", False),
    ("10_fig4_supraadditivity.py", [], "Fig 4  supraadditivity", False),
    ("08_make_tables.py", [], "Tables 1-5, 7 and S1-S9", False),
]


def run(script, args, cwd):
    t0 = time.time()
    r = subprocess.run([sys.executable, str(SCRIPTS / script), *args],
                       cwd=cwd, capture_output=True, text=True)
    return r, time.time() - t0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true",
                    help="skip the 10,000-permutation runs (reuses existing outputs)")
    ap.add_argument("--no-tests", action="store_true")
    args = ap.parse_args()

    print(f"rebuilding {ROOT.name}/\n")
    failed = []
    for script, extra, label, slow in STEPS:
        if slow and args.fast:
            print(f"  skip  {label}")
            continue
        # the figure scripts import figstyle as a sibling module
        cwd = SCRIPTS if script.startswith(("03_", "04_", "05_", "10_")) else ROOT
        r, dt = run(script, extra, cwd)
        if r.returncode == 0:
            print(f"  ok    {label}  ({dt:.0f}s)")
        else:
            print(f"  FAIL  {label}  ({dt:.0f}s)")
            print("        " + (r.stderr.strip().splitlines() or ["(no stderr)"])[-1])
            failed.append(label)
            # 01 is a gate: if the supply axis is invalid, stop
            if script.startswith("01_"):
                print("\n  the nu axis failed validation -- halting, "
                      "downstream results would be meaningless")
                return 1

    if failed:
        print(f"\n{len(failed)} step(s) failed: {', '.join(failed)}")
        return 1

    if args.no_tests:
        print("\nrebuild complete (tests skipped)")
        return 0

    print("\nrunning the test suite")
    r = subprocess.run([sys.executable, "-m", "unittest", "discover",
                        "-s", "tests"], cwd=ROOT, capture_output=True, text=True)
    tail = r.stderr.strip().splitlines()[-3:]
    for line in tail:
        print("  " + line)
    if r.returncode != 0:
        print("\nTESTS FAILED -- the manuscript and the computed data disagree")
        print(r.stderr)
        return 1
    print("\nrebuild complete and verified")
    return 0


if __name__ == "__main__":
    sys.exit(main())
