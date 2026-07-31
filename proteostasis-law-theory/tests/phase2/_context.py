"""put `scripts/` and `scripts/phase2/` on the import path for the phase 2 tests."""

import os
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
for p in (REPO_ROOT / "scripts", REPO_ROOT / "scripts" / "phase2"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))


def auditDir() -> Path:
    """the phase 2 audit directory to assert against.

    resolution order: the PHASE2_AUDIT_DIR environment variable, then the newest
    `results/phase2/audit_*` on disk. results are gitignored, so on a clean
    checkout there is nothing to assert against and the artefact-dependent tests
    skip rather than fail -- the model-level tests below still run and are the
    ones that actually pin the science.
    """
    env = os.environ.get("PHASE2_AUDIT_DIR")
    if env and Path(env).is_dir():
        return Path(env)
    root = REPO_ROOT / "results" / "phase2"
    cands = sorted(root.glob("audit_*")) if root.is_dir() else []
    return cands[-1] if cands else root / "__missing__"
