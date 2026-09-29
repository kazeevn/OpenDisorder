#!/usr/bin/env bash
#SBATCH --job-name=rt-select
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
set -euo pipefail
python "$BENCH/select_entries.py" --csv "$COD_CSV" --cod-db "$COD_DB" --procs "$SLURM_CPUS_PER_TASK" \
  --limit "${LIMIT:-0}" --out "$OUT/selected.csv" --rejected "$OUT/rejected.csv"
