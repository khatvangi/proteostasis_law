#!/bin/bash
#SBATCH -J nft_audit
#SBATCH -p compute
#SBATCH -c 16
#SBATCH --mem=12G
#SBATCH -t 03:00:00
#SBATCH -o results/slurm_%j.log
# AUDIT 2026-10-09: heavy parts of the network_failure_test audit
cd /storage/kiran-stuff/proteostasis_law/proteostasis_rebuild/network_failure_test/audit
export NPROC=16 PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1
python run_audit.py pareto cf
python run_audit.py repro
