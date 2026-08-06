"""D2: archive the baseline source of each phase-1 experiment, and be exact
about what that is and is not.

`data/PROVENANCE.md` records that the working tree was dirty at launch in all
four experiments, and `provenance.json` stores no patch. The launch-time diff is
therefore GONE and cannot be reconstructed; a diff between the recorded commit
and the next one would be a plausible-looking guess, not the code that ran, and
is not produced here.

What is recoverable is exact: the tracked source at the commit each experiment
recorded as its baseline. This writes that source out per experiment, with a
SHA-256 for every file, so the deposit carries a fixed artifact instead of a
twelve-hex commit prefix the reader has to resolve.

The scope of the gap is worth stating precisely rather than generally. Phase 1
produced the raw ensembles -- the load grid and the kinetic box. Every derived
number in the manuscript is recomputed from those raw outputs by the tracked
code in `scripts/`, under the test suite. So the provenance gap sits between the
source and the raw ensembles, not between the ensembles and the paper.
"""

from __future__ import annotations

import hashlib
import json
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
OUT = REPO_ROOT / "data" / "provenance_snapshots"

# experiment -> (baseline commit, the driver it ran), from data/PROVENANCE.md
EXPERIMENTS = {
    "A": ("73e1c0ab341530f7b2d369c459008674897f6287",
          "scripts/run_experiment_a.py"),
    "B": ("a17dfafdd966537088dc87aaf95faca893e935a3",
          "scripts/run_experiment_b.py"),
    "C": ("a17dfafdd966537088dc87aaf95faca893e935a3",
          "scripts/run_experiment_c.py"),
    "D": ("850726c76afec231948fdbd1fd1241b04e89f2b6",
          "scripts/run_experiment_d.py"),
}

# the model package the drivers import; carried with every experiment because a
# driver alone does not determine the field it integrated
PACKAGE_PREFIX = "scripts/proteostasis/"


def _git(*args: str) -> str:
    return subprocess.run(["git", *args], cwd=REPO_ROOT, check=True,
                          capture_output=True, text=True).stdout


def _treeFiles(commit: str, prefix: str) -> list[str]:
    out = _git("ls-tree", "-r", "--name-only", commit, prefix)
    return [ln for ln in out.splitlines() if ln.endswith(".py")]


def _addedAt(path: str) -> str | None:
    """the commit that first introduced a path, or None."""
    out = _git("log", "--diff-filter=A", "--format=%H", "--", path).split()
    return out[-1] if out else None


def snapshotOne(exp: str, commit: str, driver: str) -> dict:
    dest = OUT / exp
    dest.mkdir(parents=True, exist_ok=True)

    # Experiment A recorded commit 73e1c0ab as its baseline. That commit holds
    # 17 files, all markdown: NO Python existed in the repository at all, so A
    # ran entirely from an uncommitted tree and its recorded commit is not a
    # baseline in any useful sense. The source is taken from the earliest commit
    # that contains the driver, and both facts are recorded rather than one.
    recorded, fell_back = commit, False
    if not _treeFiles(commit, PACKAGE_PREFIX):
        alt = _addedAt(driver)
        if alt:
            commit, fell_back = alt, True

    wanted = [driver] + _treeFiles(commit, PACKAGE_PREFIX)

    files = {}
    for rel in wanted:
        try:
            blob = subprocess.run(["git", "show", f"{commit}:{rel}"],
                                  cwd=REPO_ROOT, check=True,
                                  capture_output=True).stdout
        except subprocess.CalledProcessError:
            # a file that did not exist at that commit is recorded as absent
            files[rel] = None
            continue
        target = dest / rel.replace("scripts/", "", 1)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(blob)
        files[rel] = hashlib.sha256(blob).hexdigest()

    return {"experiment": exp, "recorded_commit": recorded,
            "source_taken_from": commit, "driver": driver,
            "recorded_commit_held_no_source": fell_back,
            "tree_dirty_at_launch": True,
            "launch_diff_captured": False,
            "files": files}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    index = {e: snapshotOne(e, c, d) for e, (c, d) in EXPERIMENTS.items()}
    (OUT / "index.json").write_text(json.dumps(index, indent=2) + "\n")

    n = sum(len([v for v in r["files"].values() if v]) for r in index.values())
    print(f"{len(index)} experiments, {n} source files snapshotted")
    for e, r in index.items():
        miss = [k for k, v in r["files"].items() if v is None]
        note = ""
        if r["recorded_commit_held_no_source"]:
            note = (f"  <- recorded {r['recorded_commit'][:12]} held NO source; "
                    f"taken from {r['source_taken_from'][:12]}")
        print(f"  {e}  {r['source_taken_from'][:12]}  "
              f"{len(r['files']) - len(miss)} files"
              + (f"  ABSENT: {miss}" if miss else "") + note)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
