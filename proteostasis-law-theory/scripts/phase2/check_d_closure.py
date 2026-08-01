#!/usr/bin/env python
"""compact tracked check of the experiment D closure against its gitignored output.

WHY THIS EXISTS.  the numbers in `EXPERIMENT_D_FINAL.md` and
`PHASE2_CLOSURE_FINAL.md` are computed in
`results/phase2/closure_*/D_final/run_d_final.py`, which is gitignored along with
everything it writes.  a clean checkout therefore carries the prose and the
statistics (`scripts/phase2/d_final.py`) but not the evidence, and nothing would
notice if the prose drifted from the run.  this module is the tracked assertion
that closes that gap.

WHAT IT DOES NOT DO.  it does not re-run the 10,000-replicate grouped bootstrap
and it contains no copy of the estimator -- that lives once, in `d_final.py`.
It reads the shipped estimates and checks them against the values the tracked
documents state, plus the hashes that make those estimates attributable to a
specific run of a specific tree.

SKIP SEMANTICS.  when the D run root or the D_final output directory is absent,
`runChecks` returns `None` and `main` prints an explicit SKIP line and exits 0.
That is the same convention `tests/phase2` already uses: a clean checkout must
not fail for lacking gitignored results, but it must say so out loud rather than
silently reporting success.

usage: python scripts/phase2/check_d_closure.py
"""

from __future__ import annotations

import csv
import json
import math
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
if str(REPO / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO / "scripts"))

from proteostasis.provenance import hashFile, hashObject, loadConfig   # noqa: E402

PINNED = json.loads((REPO / "scripts" / "phase2" / "D_RUN_HASHES.json").read_text())
RUN_ROOT = REPO / PINNED["run_root"]

#: the two tracked documents whose numbers this module holds to the output
DOCS = ("EXPERIMENT_D_FINAL.md", "PHASE2_CLOSURE_FINAL.md")

#: (pair, null) -> (n backgrounds with a within-background majority worse than
#: null, verdict on the PRIMARY subset).  46 usable backgrounds throughout.
#: chaperone-only additive and bliss are False here on purpose: they hold only on
#: the prespecified damaging/informative subsets and are conditional, not
#: headline.
HEADLINE_MAJORITY = {
    ("influx_x_total_capacity", "additive"): (46, True),
    ("influx_x_total_capacity", "multiplicative"): (44, True),
    ("influx_x_total_capacity", "bliss"): (41, True),
    ("influx_x_chaperone_only", "additive"): (28, False),
    ("influx_x_chaperone_only", "multiplicative"): (35, True),
    ("influx_x_chaperone_only", "bliss"): (28, False),
    ("nascent_x_total_capacity", "additive"): (37, True),
    ("nascent_x_total_capacity", "multiplicative"): (38, True),
    ("nascent_x_total_capacity", "bliss"): (37, True),
}

#: the one chaperone-only null that IS supported on the primary set, with the
#: majority-fraction interval the documents quote for it
CHAPERONE_MULTIPLICATIVE_MAJORITY_CI = (0.6304347826086957, 0.8695652173913043)

#: pair -> (n backgrounds with >=1 collapse, per-background median fraction and
#: its CI, cell-level fraction and its CI)
HEADLINE_COLLAPSE = {
    "overall": (43, 0.226666666667, 0.206666666667, 0.266666666667,
                0.220579710145, 0.189855072464, 0.250434782609),
    "influx_x_total_capacity": (43, 0.4, 0.36, 0.48,
                                0.39652173913, 0.347826086957, 0.442608695652),
    "influx_x_chaperone_only": (36, 0.2, 0.2, 0.2405,
                                0.22347826087, 0.179130434783, 0.267826086957),
    "nascent_x_total_capacity": (13, 0.0, 0.0, 0.0,
                                 0.0417391304348, 0.0208695652174, 0.0660869565217),
}

#: numeric strings each document must state, checked after markdown emphasis is
#: stripped.  these are the counts and the two headline rates; if a document
#: stops saying them, it has stopped reporting this run.
DOC_TOKENS = {
    "EXPERIMENT_D_FINAL.md": ("60", "58", "46", "12", "19 and 37", "4968", "3450",
                              "36 / 36", "43 / 46 = 0.9348", "0.4000", "20260801",
                              "0.2206"),
    "PHASE2_CLOSURE_FINAL.md": ("58 completed", "46 usable", "12", "19 and 37",
                                "4968", "3450", "36/36", "43 / 46 = 0.9348",
                                "0.2206", "20260801", "71.7"),
}


def resultsPresent() -> bool:
    """whether there is anything on this host to assert against."""
    return RUN_ROOT.is_dir() and dFinalDir() is not None


def dFinalDir() -> Path | None:
    """the newest `results/phase2/closure_*/D_final` that has a results json.

    resolved rather than hardcoded so a re-run under a new timestamp is picked up
    without editing tracked source.
    """
    root = REPO / "results" / "phase2"
    if not root.is_dir():
        return None
    for c in sorted(root.glob("closure_*/D_final"), reverse=True):
        if (c / "d_final_results.json").exists():
            return c
    return None


def readTsv(path: Path) -> list:
    with path.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def near(got, want, name: str) -> tuple:
    ok = math.isclose(float(got), float(want), rel_tol=1e-9, abs_tol=1e-12)
    return (name, ok, f"{float(got):.12g} vs {float(want):.12g}")


def normalise(text: str) -> str:
    """strip markdown emphasis and collapse whitespace, so a token match is about
    the number and not about how the sentence was wrapped."""
    return re.sub(r"\s+", " ", text.replace("*", "").replace("`", ""))


# --------------------------------------------------------------------------

def checkTrackedSources() -> list:
    """the D result is attributable to this tree only while these match."""
    out = [(f"source_hash {rel}", hashFile(REPO / rel) == want, want[:12])
           for rel, want in PINNED["source_files"].items()]
    out.append(("combined_source_hash",
                hashObject(PINNED["source_files"]) == PINNED["source_hash"],
                PINNED["source_hash"][:12]))
    cfg = loadConfig(REPO / PINNED["config_path"])
    out.append(("config_hash", hashObject(cfg) == PINNED["config_hash"],
                PINNED["config_hash"][:12]))
    return out


def checkRunArtifacts() -> list:
    return [(f"artifact_hash {name}", hashFile(RUN_ROOT / name) == want, want[:12])
            for name, want in PINNED["artifact_sha256"].items()]


def checkValidationAndCounts(d_final: Path) -> list:
    """validate_d_run.py's own verdict, and its counts against the pinned manifest."""
    v = json.loads((d_final / "validation.json").read_text())
    c, p = v["counts"], PINNED["counts"]
    out = [("validation_all_passed", v["all_passed"] is True, str(v["all_passed"])),
           ("validation_36_of_36", (v["n_checks"], v["n_passed"]) == (36, 36),
            f'{v["n_passed"]}/{v["n_checks"]}'),
           ("no_check_failed", all(x["passed"] for x in v["checks"]),
            f'{len(v["checks"])} checks')]
    for key, pinned_key in (("n_requested", "n_backgrounds_requested"),
                            ("n_completed", "n_backgrounds_completed"),
                            ("n_usable", "n_usable_backgrounds"),
                            ("n_unusable", "n_unusable_backgrounds"),
                            ("n_unresolved_timeout", "n_unresolved_timeout"),
                            ("n_failed_process", "n_failed_process"),
                            ("n_cells", "n_cells"),
                            ("n_double_cells", "n_double_cells")):
        out.append((f"count {key}", c[key] == p[pinned_key],
                    f"{c[key]} vs {p[pinned_key]}"))
    # the count identities the closure prose rests on
    out.append(("46 + 12 == 58", c["n_usable"] + c["n_unusable"] == c["n_completed"],
                f'{c["n_usable"]}+{c["n_unusable"]}={c["n_completed"]}'))
    out.append(("58 + 2 + 0 == 60",
                c["n_completed"] + c["n_unresolved_timeout"] + c["n_failed_process"]
                == c["n_requested"], str(c["n_requested"])))
    out.append(("timeouts_are_not_failures", c["n_failed_process"] == 0
                and c["n_unresolved_timeout"] == 2, "0 failed, 2 unresolved"))
    return out


def checkHeadlineEstimates(d_final: Path) -> list:
    """majority counts, verdicts and the collapse rates the documents report."""
    out = []
    est = {(r["pair"], r["null"]): r
           for r in readTsv(d_final / "pair_null_estimates.tsv") if r["subset"] == "all"}
    for key, (n_major, supported) in HEADLINE_MAJORITY.items():
        r = est[key]
        tag = f"{key[0]}/{key[1]}"
        out.append((f"majority {tag}", int(r["n_backgrounds_majority_worse"]) == n_major,
                    f'{r["n_backgrounds_majority_worse"]} of {r["n_backgrounds"]}'))
        out.append((f"n_backgrounds {tag}", int(r["n_backgrounds"]) == 46,
                    r["n_backgrounds"]))
        out.append((f"primary_supported {tag}",
                    (r["supported"] == "True") is supported, r["supported"]))

    # the two chaperone-only nulls fail for the stated reason, not some other one:
    # the median CI touches or crosses zero AND the majority CI includes one half
    for null in ("additive", "bliss"):
        r = est[("influx_x_chaperone_only", null)]
        lo, hi = float(r["bg_median_excess_lo"]), float(r["bg_median_excess_hi"])
        out.append((f"chaperone {null} median CI includes zero", lo <= 0.0 <= hi,
                    f"[{lo:.6g}, {hi:.6g}]"))
        out.append((f"chaperone {null} majority CI includes half",
                    float(r["frac_backgrounds_majority_worse_lo"]) <= 0.5
                    <= float(r["frac_backgrounds_majority_worse_hi"]),
                    f'[{r["frac_backgrounds_majority_worse_lo"]}, '
                    f'{r["frac_backgrounds_majority_worse_hi"]}]'))
        out.append((f"chaperone {null} failing criteria recorded",
                    r["failing_criteria"] != "", r["failing_criteria"]))

    r = est[("influx_x_chaperone_only", "multiplicative")]
    lo, hi = CHAPERONE_MULTIPLICATIVE_MAJORITY_CI
    out.append(near(r["frac_backgrounds_majority_worse_lo"], lo,
                    "chaperone multiplicative majority CI lo"))
    out.append(near(r["frac_backgrounds_majority_worse_hi"], hi,
                    "chaperone multiplicative majority CI hi"))

    coll = {r["pair"]: r for r in readTsv(d_final / "collapse_estimates.tsv")}
    for pair, (n_bg, med, med_lo, med_hi, cell, cell_lo, cell_hi) in HEADLINE_COLLAPSE.items():
        r = coll[pair]
        out.append((f"collapse n_backgrounds {pair}",
                    int(r["n_backgrounds_with_collapse"]) == n_bg,
                    f'{r["n_backgrounds_with_collapse"]} of {r["n_backgrounds"]}'))
        for col, want in (("bg_median_frac_collapse", med),
                          ("bg_median_frac_collapse_lo", med_lo),
                          ("bg_median_frac_collapse_hi", med_hi),
                          ("cell_frac_collapse", cell),
                          ("cell_frac_collapse_lo", cell_lo),
                          ("cell_frac_collapse_hi", cell_hi)):
            out.append(near(r[col], want, f"collapse {pair} {col}"))
    return out


def checkSensitivityArithmetic(d_final: Path) -> list:
    """every bound is k over a stated denominator, and nothing is imputed."""
    out = []
    rows = readTsv(d_final / "sensitivity_bounds.tsv")
    out.append(("sensitivity rows", len(rows) == 13, str(len(rows))))
    for r in rows:
        k, tag = int(r["k"]), f'{r["quantity"][:24]} {r["pair"]}/{r["null"]}'
        out.append((f"denominators {tag}",
                    (int(r["n_usable"]), int(r["n_unresolved"]), int(r["n_unusable"]),
                     int(r["n_requested"])) == (46, 2, 12, 60), r["n_requested"]))
        out.append(near(r["conditional"], k / 46, f"conditional {tag}"))
        out.append(near(r["requested_lo"], k / 60, f"requested_lo {tag}"))
        out.append(near(r["requested_hi"], (k + 2) / 60, f"requested_hi {tag}"))
        out.append((f"usable bounds bracket conditional {tag}",
                    float(r["usable_lo"]) <= float(r["conditional"]) + 1e-12
                    and float(r["usable_hi"]) >= float(r["conditional"]) - 1e-12,
                    f'[{r["usable_lo"]}, {r["usable_hi"]}]'))
    return out


def checkDocuments() -> list:
    out = []
    for name, tokens in DOC_TOKENS.items():
        text = normalise((REPO / name).read_text())
        for t in tokens:
            out.append((f"{name} states {t!r}", t in text, ""))
    return out


# --------------------------------------------------------------------------

def runChecks():
    """every check as (name, ok, detail), or None when there is nothing to check."""
    d_final = dFinalDir()
    if not RUN_ROOT.is_dir() or d_final is None:
        return None
    return (checkTrackedSources() + checkRunArtifacts()
            + checkValidationAndCounts(d_final) + checkHeadlineEstimates(d_final)
            + checkSensitivityArithmetic(d_final) + checkDocuments())


def main() -> int:
    checks = runChecks()
    if checks is None:
        missing = RUN_ROOT if not RUN_ROOT.is_dir() else "results/phase2/closure_*/D_final"
        print(f"SKIP: experiment D closure output absent ({missing}); "
              f"nothing to assert on this host. tracked prose and statistics "
              f"were not checked against a run.")
        return 0
    bad = [c for c in checks if not c[1]]
    for name, ok, detail in checks:
        print(f"{'ok  ' if ok else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ""))
    print(f"\n{len(checks) - len(bad)}/{len(checks)} checks pass "
          f"against {dFinalDir()}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
