#!/usr/bin/env bash
# SLURM version of run_all.sh: submits the four steps as dependent (array) jobs.
#
#   COD_CSV=/data/cod/all.csv COD_DB=/data/cod/cod.sqlite OUT=/scratch/$USER/roundtrip \
#   BACKENDS="mace orb chgnet" bash slurm/pipeline.sh
#
# Partition names, accounts, modules/conda environments are site specific: edit the
# #SBATCH lines in the *.sbatch files (or pass extra sbatch options in SBATCH_EXTRA).
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
export COD_DB OUT COD_CSV BENCH="$(dirname "$HERE")"
: "${COD_DB:?}" "${OUT:?}" "${COD_CSV:?}"
BACKENDS="${BACKENDS:-mace}"
NGEN="${NGEN:-16}"      # generation shards (CPU)
NRELAX="${NRELAX:-32}"  # relaxation shards per backend (GPU)
NEVAL="${NEVAL:-32}"    # evaluation shards (CPU)
mkdir -p "$OUT/logs"
X=(${SBATCH_EXTRA:-})

sel=$(sbatch --parsable "${X[@]}" "$HERE/select.sbatch")
gen=$(sbatch --parsable "${X[@]}" --dependency=afterok:$sel --array=0-$((NGEN-1)) --export=ALL,NSHARDS=$NGEN "$HERE/generate.sbatch")
deps=""
for b in $BACKENDS; do
  r=$(sbatch --parsable "${X[@]}" --dependency=afterok:$gen --array=0-$((NRELAX-1)) \
      --export=ALL,BACKEND=$b,NSHARDS=$NRELAX "$HERE/relax.sbatch")
  deps="$deps:$r"
done
ev=$(sbatch --parsable "${X[@]}" --dependency=afterok$deps --array=0-$((NEVAL-1)) --export=ALL,NSHARDS=$NEVAL "$HERE/evaluate.sbatch")
sbatch "${X[@]}" --dependency=afterok:$ev --job-name=rt-summary --time=00:30:00 --cpus-per-task=1 \
  --output="$OUT/logs/summary.out" --wrap="python $BENCH/evaluate.py --out-dir $OUT/results --summarize-only"
echo "submitted: select $sel, generate $gen, relax$deps, evaluate $ev"
