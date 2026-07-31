#!/bin/bash
set -euo pipefail
cd /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory-codex-check
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
exec /home/kiran/miniforge3/bin/python3 sweep.py --config config/phase1.json --output results/phase1
