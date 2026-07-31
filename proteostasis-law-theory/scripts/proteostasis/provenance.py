"""configuration loading, provenance stamping and content-addressed outputs.

every experiment writes a `provenance.json` recording the git commit, the
config hash, package versions, the seed, and a sha256 of each output file. the
reproducibility test reruns a small deterministic slice and compares hashes, so
an accidental dependence on wall-clock time, dict ordering or unseeded rng
becomes a test failure rather than a silent drift.
"""

from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]


def loadConfig(path: str | os.PathLike) -> Dict[str, Any]:
    with open(path, "r") as fh:
        cfg = json.load(fh)
    if "experiment" not in cfg or "seed" not in cfg:
        raise ValueError(f"config {path} must declare 'experiment' and 'seed'")
    return cfg


def canonicalJson(obj: Any) -> str:
    """stable serialization: sorted keys, fixed separators, no wall-clock data."""
    def default(o):
        if is_dataclass(o):
            return asdict(o)
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, Path):
            return str(o)
        raise TypeError(f"not serializable: {type(o)}")
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=default)


def hashObject(obj: Any) -> str:
    return hashlib.sha256(canonicalJson(obj).encode()).hexdigest()


def hashFile(path: str | os.PathLike) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def gitCommit(repo: Path = REPO_ROOT) -> Dict[str, Any]:
    def run(*args):
        try:
            return subprocess.run(args, cwd=repo, capture_output=True, text=True,
                                  check=True).stdout.strip()
        except Exception:
            return None
    return {
        "commit": run("git", "rev-parse", "HEAD"),
        "branch": run("git", "rev-parse", "--abbrev-ref", "HEAD"),
        "dirty": bool(run("git", "status", "--porcelain")),
        "remote": run("git", "remote", "-v") or "",
    }


def environmentInfo() -> Dict[str, Any]:
    import scipy
    return {
        "python": sys.version.split()[0],
        "executable": sys.executable,
        "platform": platform.platform(),
        "node": platform.node(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pandas": pd.__version__,
    }


def writeTable(df: pd.DataFrame, path: str | os.PathLike) -> str:
    """write a tidy TSV deterministically (fixed column order, fixed float format)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    out = df.reindex(sorted(df.columns), axis=1)
    out.to_csv(path, sep="\t", index=False, float_format="%.12g", lineterminator="\n")
    return hashFile(path)


def writeJson(obj: Any, path: str | os.PathLike) -> str:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(json.dumps(obj, indent=2, sort_keys=True,
                            default=lambda o: canonicalJson(o) and json.loads(canonicalJson(o))))
        fh.write("\n")
    return hashFile(path)


def writeProvenance(outdir: str | os.PathLike, config: Dict[str, Any],
                    outputs: Dict[str, str], extra: Optional[Dict] = None) -> Path:
    """record everything needed to decide whether a rerun reproduced a result."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    prov = {
        "experiment": config.get("experiment"),
        "config": config,
        "config_hash": hashObject(config),
        "git": gitCommit(),
        "environment": environmentInfo(),
        "slurm": {k: v for k, v in os.environ.items() if k.startswith("SLURM_")},
        "output_sha256": outputs,
        "extra": extra or {},
    }
    path = outdir / "provenance.json"
    with open(path, "w") as fh:
        fh.write(json.dumps(prov, indent=2, sort_keys=True))
        fh.write("\n")
    return path
