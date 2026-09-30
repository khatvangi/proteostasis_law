"""
run every permanent check for stages 0-4 and write results.json.

    cd checks && PYTHONDONTWRITEBYTECODE=1 python run_checks.py

regenerates ../EQUATION_TERM_LEDGER.tsv first (write_ledger.py), then runs
check_conservation, check_baseline and check_stage4_fold. exit 1 if any fails.
every number quoted in files 00-04 and STATUS.md is read from results.json.
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import write_ledger  # noqa: E402
import check_conservation  # noqa: E402
import check_baseline  # noqa: E402
import check_stage4_fold  # noqa: E402


def clean(o):
    if isinstance(o, dict):
        return {str(k): clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.floating, np.integer, np.bool_)):
        return o.item()
    return o


def main():
    write_ledger.main()
    out = {}
    for name, mod in (("conservation", check_conservation), ("baseline", check_baseline),
                      ("stage4", check_stage4_fold)):
        t = time.time()
        out[name] = mod.run()
        print(f"{name}: {time.time() - t:.0f} s")
    passes = {f"{g}.{k}": bool(v["pass"]) for g, r in out.items() for k, v in r.items()}
    out["summary"] = {"all_pass": all(passes.values()), "n_checks": len(passes), "checks": passes}
    (HERE / "results.json").write_text(json.dumps(clean(out), indent=1, default=str) + "\n")
    for k, v in passes.items():
        print(("PASS " if v else "FAIL ") + k)
    return 0 if out["summary"]["all_pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
