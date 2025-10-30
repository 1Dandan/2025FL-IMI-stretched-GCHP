#!/bin/bash
set -euo pipefail

read -rp "An example successful IMI run tag (only 3-digit padded number): " SuccessIMIRunTag

# -------------------------
# User-configurable section
# -------------------------
Dir="/n/holylfs06/LABS/jacob_lab2/Lab/dzhang8/imi-gchp-test"
imiDir="imi-stretch-gchp-cp2"
outdirsubname="configs_faces_for_permian"
outdir="${Dir}/${imiDir}/${outdirsubname}" 
runname_base="Test_Stretch_1month"

# Where to save results
failed_list="${Dir}/${imiDir}/failed_runs.txt"
: > "$failed_list"

SuccessInvRun="${Dir}/output-gchp-stretch-soil/${runname_base}_T${SuccessIMIRunTag}/inversion"
THRESHOLD=$(find "${SuccessInvRun}/data_converted" -type f | wc -l)


# What counts as "enough work done" inside inversion/ ?
# Choose the counting mode and threshold:
COUNT_MODE="${COUNT_MODE:-dirs}"   # dirs | files
THRESHOLD="${THRESHOLD:-10}"       # integer

# Optional filter: only test tags that match this regex (e.g., '^T0(01|02)$')
ONLY_TAGS_REGEX="${ONLY_TAGS_REGEX:-}"

# -------------------------
# Start
# -------------------------
tag_file="${1:-}"

# -------------------------
# Build config list
# -------------------------
configs=()
if [[ -n "${Dir}/${imiDir}/$tag_file" && -f "${Dir}/${imiDir}/$tag_file" ]]; then
  while read -r raw; do
    tag="$(echo "$raw" | xargs)"               # trim
    [[ -z "$tag" ]] && continue
    cfg="${outdir}/config-soil_${tag}.yml"
    [[ -f "$cfg" ]] && configs+=("$cfg")
  done < "${Dir}/${imiDir}/$tag_file"
else
  shopt -s nullglob
  configs=( "${outdir}"/config-soil_T*.yml )
  shopt -u nullglob
fi

if (( ${#configs[@]} == 0 )); then
  echo "No configs found in ${outdir} (pattern config-soil_T*.yml). Nothing to do."
  exit 0
fi

for yml in "${configs[@]}"; do
  base="$(basename "$yml")"
  if [[ ! "$base" =~ _T([0-9]+)\.yml$ ]]; then
    echo "Skipping non-matching file: $base"
    continue
  fi
  tag="T${BASH_REMATCH[1]}"

  # Regex filter
  if [[ -n "$ONLY_TAGS_REGEX" ]] && ! [[ "$tag" =~ $ONLY_TAGS_REGEX ]]; then
    continue
  fi

  run_root="${Dir}/output-gchp-stretch-soil/${runname_base}_${tag}"
  invdir="${run_root}/inversion"

  # Condition: EITHER inversion does not exist OR item count < THRESHOLD
  should_flag=false
  if [[ ! -d "$invdir" ]]; then
    should_flag=true
  else
    item_count=$(find "$invdir/data_converted" -type f | wc -l)
    if (( item_count < THRESHOLD )); then
      should_flag=true
    fi
  fi

  if $should_flag; then
    # Append once per unique tag
    if ! grep -qx "$tag" "$failed_list"; then
      echo "$tag" >> "$failed_list"
    fi
  fi

done