#!/usr/bin/env bash
# Round-trip benchmark, whole pipeline on one machine.  See README.md; for a cluster use the
# SLURM scripts in slurm/ (same steps, sharded).
#
#   COD_CSV=/data/cod/all.csv COD_DB=~/.cache/opendisorder/cod.sqlite OUT=/data/roundtrip \
#   BACKENDS="mace chgnet" DEVICE=cuda PROCS=16 ./run_all.sh
#
# Environment variables (defaults in brackets):
#   COD_CSV    COD CSV metadata dump, e.g. from
#              https://www.crystallography.net/cod/result?format=csv  (or the one opendisorder
#              build-index downloaded)                                        [required]
#   COD_DB     opendisorder index with CIFs (build-index + import-archive); opened read-only
#   OUT        output directory                                               [./roundtrip_out]
#   METHODS    ordering methods                                               [random,ewald,sqs,enum]
#   BACKENDS   MLIP backends for relax.py (space separated)                   [mace]
#   DEVICE     cpu | cuda                                                      [cpu]
#   PROCS      worker processes for selection / evaluation                    [4]
#   LIMIT      stratified sample of parents (0 = all that pass the filter)    [0]
#   SMOKE=1    tiny run: 3 parents, few orderings, noop relaxation
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
COD_DB="${COD_DB:-$HOME/.cache/opendisorder/cod.sqlite}"
OUT="${OUT:-./roundtrip_out}"
METHODS="${METHODS:-random,ewald,sqs,enum}"
BACKENDS="${BACKENDS:-mace}"
DEVICE="${DEVICE:-cpu}"
PROCS="${PROCS:-4}"
LIMIT="${LIMIT:-0}"
mkdir -p "$OUT"

if [[ "${SMOKE:-0}" == 1 ]]; then
  SEL_ARGS=(--ids "${SMOKE_IDS:-1005047,1008793,1008967}")
  GEN_ARGS=(--methods random,ewald,sqs --n-sizes 2 --max-atoms 60 --n-per-size 3 --controls 2 --sqs-iterations 20000)
  BACKENDS="noop"
else
  SEL_ARGS=(--csv "${COD_CSV:?set COD_CSV}" --limit "$LIMIT" --procs "$PROCS")
  GEN_ARGS=(--methods "$METHODS" --sizes auto --n-sizes 3 --max-atoms 120 --shapes 2 --n-per-size 5 --controls 3)
fi

# 1. quality-filtered, compositionally disordered parents
python "$HERE/select_entries.py" "${SEL_ARGS[@]}" --cod-db "$COD_DB" \
  --out "$OUT/selected.csv" --rejected "$OUT/rejected.csv"

# 2. ordered realizations (+ controls); enum needs enum.x/makestr.x (conda-forge enumlib) on PATH
python "$HERE/generate_orderings.py" --selected "$OUT/selected.csv" --cod-db "$COD_DB" \
  --out "$OUT/orderings.jsonl" "${GEN_ARGS[@]}"

# 3. MLIP relaxations, one file per backend (resumable)
for b in $BACKENDS; do
  python "$HERE/relax.py" --in "$OUT/orderings.jsonl" --backend "$b" --device "$DEVICE" \
    --out "$OUT/relaxed_$b.jsonl"
done

# 4. match back to COD: parent recall, decoy false positives, NoveltyChecker end to end
python "$HERE/evaluate.py" --selected "$OUT/selected.csv" \
  --structures "$OUT/orderings.jsonl" "$OUT"/relaxed_*.jsonl \
  --cod-db "$COD_DB" --out-dir "$OUT/results" --procs "$PROCS" --novelty

echo "summary: $OUT/results/summary.md"
