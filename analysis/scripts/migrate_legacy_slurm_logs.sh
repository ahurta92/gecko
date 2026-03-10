#!/usr/bin/env bash
set -euo pipefail

DEST_ROOT=${1:-/gpfs/projects/rjh/adrian/development/gecko/analysis/slurm_logs/legacy}
mkdir -p "$DEST_ROOT"

sources=(
  /gpfs/scratch/ahurtado/project_data/freq_dev
  /gpfs/scratch/ahurtado/project_data/response_scaling_dev
  /gpfs/scratch/ahurtado/project_data/data/raman_paper
  /gpfs/scratch/ahurtado/project_data/data/NLO
  /gpfs/scratch/ahurtado/project_data/excited_state_dev
)

echo "Destination: $DEST_ROOT"

total=0
copied=0
skipped=0

for src in "${sources[@]}"; do
  if [[ ! -d "$src" ]]; then
    continue
  fi

  src_tag=$(basename "$src")
  while IFS= read -r -d '' f; do
    total=$((total + 1))
    rel=${f#"$src"/}
    dst="$DEST_ROOT/$src_tag/$rel"
    dstdir=$(dirname "$dst")
    mkdir -p "$dstdir"

    if [[ -f "$dst" ]]; then
      skipped=$((skipped + 1))
      continue
    fi

    cp -p "$f" "$dst"
    copied=$((copied + 1))
  done < <(find "$src" -type f \( -name '*.out' -o -name '*.err' -o -name 'slurm-*.out' \) -print0)
done

echo "Total found: $total"
echo "Copied:      $copied"
echo "Skipped:     $skipped"
