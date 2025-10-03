#!/bin/bash
set -euo pipefail

Dir="/n/holylfs06/LABS/jacob_lab2/Lab/dzhang8/imi-gchp-test"
imiDir="imi-gchp-precomputedK"
outputDir="output-gchp-stretch-soil"
target_coords_path="${Dir}/supportData/target_coords.csv"
out="reference_rundir_path.txt"

cd "${Dir}/${imiDir}"

# Count data rows (exclude header). Robust if last line lacks trailing newline.
n_face=$(awk -F, 'NR>1 && NF>=2 {c++} END {print c+0}' "$target_coords_path")
if [[ "$n_face" -eq 0 ]]; then
  echo "No data rows found in $target_coords_path" >&2
  exit 1
fi
pad=${#n_face}

tmp=$(mktemp)
awk -F, -v base="${Dir}/${outputDir}/Test_Stretch_1day_T" -v pad="$pad" '
  NR==1 {next}                 # skip header
  NF<2 {next}                  # skip incomplete/blank lines
  {
    # strip Windows CRs and any stray spaces on first two fields
    for(i=1;i<=2;i++){ gsub(/\r/,"",$i); sub(/^[ \t]+/,"",$i); sub(/[ \t]+$/,"",$i) }
    tag = sprintf("%0*d", pad, NR-1)
    print base tag
  }
' "$target_coords_path" > "$tmp"

mv -f "$tmp" "$out"