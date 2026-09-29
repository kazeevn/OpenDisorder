#!/usr/bin/env bash
#SBATCH --job-name=rt-generate
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=%x-%A_%a.out
# enum/enum_ewald need enumlib: e.g. `conda install -c conda-forge enumlib` in the job's env
set -euo pipefail
export OMP_NUM_THREADS=1
python "$BENCH/generate_orderings.py" --selected "$OUT/selected.csv" --cod-db "$COD_DB" \
  --out "$OUT/orderings/orderings.shard$SLURM_ARRAY_TASK_ID.jsonl" --shard "slurm/$NSHARDS" \
  --methods "${METHODS:-random,ewald,sqs,enum}" --sizes auto --n-sizes 3 --max-atoms 120 \
  --shapes 2 --n-per-size 5 --controls 3 --enum-mem-gb 6
