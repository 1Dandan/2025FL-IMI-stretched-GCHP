#!/bin/bash
set -euo pipefail

Dir="/n/holylfs06/LABS/jacob_lab2/Lab/dzhang8/imi-gchp-test"
imiDir="imi-gchp-precomputedK"
outdirsubname="configs_faces_for_permian"
outdir="${Dir}/${imiDir}/${outdirsubname}" 
outputDir="output-gchp-stretch-soil"
runname_base="Test_Stretch_1month"
target_coords_path="${Dir}/supportData/target_coords.csv"
out="reference_rundir_path_1month.txt"

cd "${Dir}/${imiDir}"

# Count data rows (exclude header)
n_face=$(awk -F, 'NR>1 && NF>=2 {c++} END {print c+0}' "$target_coords_path")
if [[ "$n_face" -eq 0 ]]; then
  echo "No data rows found in $target_coords_path" >&2
  exit 1
fi

# Zero-padding width based on number of rows (e.g., 220 -> 3)
pad=${#n_face}

# Optional tag file: accept absolute/relative or ${Dir}/${imiDir}/<name>
tag_file_arg="${1:-}"
tag_file=""
if [[ -n "$tag_file_arg" ]]; then
  if   [[ -f "$tag_file_arg" ]]; then tag_file="$tag_file_arg"
  elif [[ -f "${Dir}/${imiDir}/$tag_file_arg" ]]; then tag_file="${Dir}/${imiDir}/$tag_file_arg"
  else
    echo "Warning: tag file '$tag_file_arg' not found; proceeding with all indices." >&2
  fi
fi

# Build indices
indices=()
if [[ -n "$tag_file" ]]; then
  while read -r raw; do
    tag=$(echo "$raw" | tr -d '\r' | xargs)   # trim
    [[ -z "$tag" ]] && continue
    idx="${tag#T}"; idx=$((10#$idx))          # strip leading T, force base-10
    indices+=( "$idx" )
  done < "$tag_file"
else
  mapfile -t indices < <(seq "$n_face")
fi


tmp=$(mktemp)
for idx in "${indices[@]}"; do
  # zero-padded tag, e.g., 1 -> 001
  printf -v tag "%0${pad}d" "$idx"

  rundir_path="${Dir}/${outputDir}/${runname_base}_T${tag}"
  echo ${rundir_path} >> "$tmp"
done

mv -f "$tmp" "$out"