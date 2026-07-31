#!/usr/bin/env python3
"""
build the Zenodo deposition tarball, then verify it by using it.

a deposition that has not been extracted and run somewhere else is a guess. this
script builds the archive, extracts it to a scratch directory, and runs the
pipeline and the test suite THERE, so what is uploaded is known to be
self-contained rather than assumed to be.

    python scripts/16_make_deposit.py                # build + verify
    python scripts/16_make_deposit.py --no-verify    # build only
    python scripts/16_make_deposit.py --full         # verify with the full rebuild

outputs, into dist/:
    <name>.tar.gz          the archive
    <name>.tar.gz.sha256   its checksum, for the Zenodo record
    <name>.contents.txt    file listing with sizes
and inside the archive: SHA256SUMS, COMMIT, BUILT.
"""
import argparse
import hashlib
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DIST = ROOT / "dist"
VERSION = "v1.0.0"
NAME = f"proteostasis-envelope-{VERSION}"

# what goes in. everything else -- dist/, __pycache__, .gitignore -- stays out.
INCLUDE = ["README.md", "DEPOSIT.md", "manuscript", "scripts", "data", "tables",
           "figures", "tests"]
EXCLUDE_NAMES = {"__pycache__", ".gitignore", ".DS_Store", "dist"}
EXCLUDE_SUFFIX = {".pyc", ".pyo"}


def keep(path: Path) -> bool:
    if any(part in EXCLUDE_NAMES for part in path.parts):
        return False
    return path.suffix not in EXCLUDE_SUFFIX


def collect():
    files = []
    for entry in INCLUDE:
        p = ROOT / entry
        if not p.exists():
            sys.exit(f"ERROR: {entry} is missing; refusing to build a partial "
                     "deposition")
        if p.is_file():
            files.append(p)
        else:
            files += [q for q in sorted(p.rglob("*")) if q.is_file() and keep(q)]
    return files


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def git(*args, default=""):
    try:
        return subprocess.run(["git", *args], cwd=ROOT, capture_output=True,
                              text=True, check=True).stdout.strip()
    except Exception:
        return default


def build(built_at):
    DIST.mkdir(exist_ok=True)
    files = collect()

    # the archive states which commit it came from and whether that commit was
    # clean. an archive built from a dirty tree is not reproducible from the
    # repository, and saying so is more useful than hiding it.
    head = git("rev-parse", "HEAD", default="unknown")
    dirty = [ln for ln in git("status", "--porcelain").splitlines()
             if ln and "proteostasis-paper" not in ln]
    commit_txt = (f"{head}\n"
                  f"repository: https://github.com/khatvangi/proteostasis_law\n"
                  f"path: envelope-paper/\n"
                  f"tree state at build: "
                  f"{'DIRTY -- ' + str(len(dirty)) + ' uncommitted change(s)'
                     if dirty else 'clean'}\n")
    if dirty:
        commit_txt += "".join(f"  {ln}\n" for ln in dirty)

    sums = "\n".join(f"{sha256(f)}  {f.relative_to(ROOT).as_posix()}"
                     for f in files) + "\n"

    tar_path = DIST / f"{NAME}.tar.gz"
    with tarfile.open(tar_path, "w:gz") as tar:
        for f in files:
            tar.add(f, arcname=f"{NAME}/{f.relative_to(ROOT).as_posix()}")
        for fname, text in (("SHA256SUMS", sums), ("COMMIT", commit_txt),
                            ("BUILT", built_at + "\n")):
            data = text.encode()
            info = tarfile.TarInfo(f"{NAME}/{fname}")
            info.size = len(data)
            info.mtime = 0                      # keep the archive reproducible
            import io
            tar.addfile(info, io.BytesIO(data))

    digest = sha256(tar_path)
    (DIST / f"{NAME}.tar.gz.sha256").write_text(
        f"{digest}  {tar_path.name}\n")
    listing = "\n".join(
        f"{f.stat().st_size:>10}  {f.relative_to(ROOT).as_posix()}" for f in files)
    (DIST / f"{NAME}.contents.txt").write_text(listing + "\n")

    print(f"built {tar_path.relative_to(ROOT)}")
    print(f"  {len(files)} files, {tar_path.stat().st_size / 1e6:.1f} MB")
    print(f"  sha256 {digest}")
    print(f"  commit {head[:12]}  ({'DIRTY' if dirty else 'clean'} tree)")
    return tar_path, bool(dirty)


def verify(tar_path, full):
    """extract elsewhere and run it, so 'self-contained' is a finding not a hope."""
    with tempfile.TemporaryDirectory(prefix="deposit-verify-") as tmp:
        with tarfile.open(tar_path) as tar:
            tar.extractall(tmp, filter="data")
        root = Path(tmp) / NAME
        print(f"\nverifying the extracted archive in {root}")

        # 1. checksums inside the archive must match what was extracted
        bad = []
        for line in (root / "SHA256SUMS").read_text().splitlines():
            want, rel = line.split("  ", 1)
            got = sha256(root / rel)
            if got != want:
                bad.append(rel)
        print(f"  checksums: {'OK' if not bad else 'MISMATCH ' + str(bad[:3])}")
        if bad:
            return False

        # 2. the pipeline must run from the archive alone
        cmd = [sys.executable, "scripts/run_all.py"] + ([] if full else ["--fast"])
        r = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
        tail = [ln for ln in r.stdout.strip().splitlines() if ln.strip()][-6:]
        for ln in tail:
            print("  " + ln)
        if r.returncode != 0:
            print("  FAILED: the archive does not reproduce on its own")
            print(r.stdout[-2000:])
            return False

        # 3. and the suite must pass there, not just here
        t = subprocess.run([sys.executable, "-m", "unittest", "discover",
                            "-s", "tests"], cwd=root, capture_output=True,
                           text=True)
        print("  test suite in the extracted copy: "
              + (t.stderr.strip().splitlines() or ["?"])[-1])
        return t.returncode == 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-verify", action="store_true")
    ap.add_argument("--full", action="store_true",
                    help="verify with the full rebuild rather than --fast")
    ap.add_argument("--built", default=None,
                    help="timestamp to record (defaults to now)")
    a = ap.parse_args()

    # passed in rather than read from the clock inside the build, so the archive
    # is reproducible from a known commit + timestamp pair
    built = a.built or subprocess.run(["date", "-u", "+%Y-%m-%dT%H:%M:%SZ"],
                                      capture_output=True, text=True).stdout.strip()
    tar_path, dirty = build(built)
    if a.no_verify:
        print("\nskipped verification (--no-verify)")
        return 0
    ok = verify(tar_path, a.full)
    print("\n" + ("DEPOSITION VERIFIED: extracts, rebuilds and passes its own "
                  "tests in a clean directory" if ok else
                  "DEPOSITION FAILED VERIFICATION -- do not upload"))
    if ok and dirty:
        print("NOTE: built from a dirty tree; commit first if the archive is "
              "meant to correspond to a public commit")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
