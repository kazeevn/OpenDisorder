#!/usr/bin/env bash
#SBATCH --job-name=rt-relax
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=24:00:00
#SBATCH --output=%x-%A_%a.out
# One GPU per shard; the job can be resubmitted after a time-out: finished structures are
# skipped (resumable).  For CPU-only partitions drop --gres and set DEVICE=cpu.
set -euo pipefail
BACKEND="${BACKEND:-mace}"
python "$BENCH/relax.py" --in "$OUT"/orderings/orderings.shard*.jsonl --backend "$BACKEND" \
  ${MODEL:+--model "$MODEL"} --device "${DEVICE:-cuda}" --threads "$SLURM_CPUS_PER_TASK" \
  --shard "slurm/$NSHARDS" --out "$OUT/relaxed/relaxed_$BACKEND.shard$SLURM_ARRAY_TASK_ID.jsonl"
