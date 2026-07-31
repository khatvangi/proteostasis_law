#!/bin/bash
#SBATCH --job-name=prot-law-p1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:45:00
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
set -euo pipefail
cd /storage/kiran-stuff/proteostasis_law/proteostasis-law-theory-codex-check
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
python3 sweep.py --config config/phase1.json --output results/phase1

