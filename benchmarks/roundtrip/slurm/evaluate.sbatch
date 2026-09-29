#!/usr/bin/env bash
#SBATCH --job-name=rt-evaluate
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%A_%a.out
# Shards are disjoint sets of parents; each shard builds its own small working index for
# NoveltyChecker (entries/keys ~ a few hundred MB + candidate CIFs), so use node-local scratch
# via WORK_DIR if the shared file system is slow.
set -euo pipefail
export OMP_NUM_THREADS=1
python "$BENCH/evaluate.py" --selected "$OUT/selected.csv" \
  --structures "$OUT"/orderings/orderings.shard*.jsonl "$OUT"/relaxed/relaxed_*.jsonl \
  --cod-db "$COD_DB" --out-dir "$OUT/results" --procs "$SLURM_CPUS_PER_TASK" --novelty \
  --work-index "${WORK_DIR:-$OUT/results}/work_index.shard$SLURM_ARRAY_TASK_ID.sqlite" \
  --shard "slurm/$NSHARDS"
