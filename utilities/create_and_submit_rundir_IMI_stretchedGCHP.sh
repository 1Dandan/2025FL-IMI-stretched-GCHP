#!/bin/bash
set -euo pipefail

Dir="/n/holylfs06/LABS/jacob_lab2/Lab/dzhang8/imi-gchp-test"
imiDir="imi-stretch-gchp-cp2"
templatefname="config-soil-1month.yml"
template="${Dir}/${imiDir}/${templatefname}"
target_coords_path="${Dir}/supportData/target_coords.csv"
outdirsubname="configs_faces_for_permian"
outdir="${Dir}/${imiDir}/${outdirsubname}"
gchp_exec_dir="${Dir}/output-gchp-stretch-soil/Test_Stretch_1month/GEOSChem_build"
runname_base="Test_Stretch_1month"

cd "${Dir}/${imiDir}"
mkdir -p "$outdir"

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

# Preload lon/lat arrays once (0-based arrays; row 1 -> index 0)
mapfile -t lon_arr < <(awk -F, 'NR>1 {gsub(/\r/,""); gsub(/^[ \t]+|[ \t]+$/, "", $1); print $1}' "$target_coords_path")
mapfile -t lat_arr < <(awk -F, 'NR>1 {gsub(/\r/,""); gsub(/^[ \t]+|[ \t]+$/, "", $2); print $2}' "$target_coords_path")

make_config() {
  local idx="$1" lon="$2" lat="$3"

  # zero-padded tag, e.g., 1 -> 001
  printf -v tag "%0${pad}d" "$idx"

  local outyml="${outdir}/config-soil_T${tag}.yml"

  cp "$template" "$outyml"
  # Update fields
  sed -i "s/^RunName:.*/RunName: \"${runname_base}_T${tag}\"/" "$outyml"
  sed -i "s/^TARGET_LAT:.*/TARGET_LAT: ${lat}/" "$outyml"
  sed -i "s/^TARGET_LON:.*/TARGET_LON: ${lon}/" "$outyml"

  # Prepare run dir + symlink
  local run_root="${Dir}/output-gchp-stretch-soil/${runname_base}_T${tag}"
  mkdir -p "$run_root"
  ln -nsf "${gchp_exec_dir}" "${run_root}/GEOSChem_build"

  # submit if needed:
  # sbatch -p serial_requeue,unrestricted,shared -t 0-12:00 --mem 16G -c 1 \
  #   -o "imi_output_T${tag}.log" run_imi.sh "$outyml"
}

# Main loop
for idx in "${indices[@]}"; do
  # Range guard
  if (( idx < 1 || idx > ${#lon_arr[@]} )); then
    echo "Skipping idx=$idx (out of range 1..${#lon_arr[@]})" >&2
    continue
  fi
  lon="${lon_arr[idx-1]}"
  lat="${lat_arr[idx-1]}"
  if [[ -z "${lon}" || -z "${lat}" ]]; then
    echo "Missing lon/lat for idx=$idx; skipping" >&2
    continue
  fi
  make_config "$idx" "$lon" "$lat"
done