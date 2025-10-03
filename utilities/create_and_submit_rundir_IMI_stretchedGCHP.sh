#!/bin/bash
set -euo pipefail

Dir="/n/holylfs06/LABS/jacob_lab2/Lab/dzhang8/imi-gchp-test"
imiDir="imi-stretch-gchp-cp2"
template="${Dir}/${imiDir}/config-soil.yml"              # template to copy from
target_coords_path="${Dir}/supportData/target_coords.csv"
outdir="${Dir}/${imiDir}/configs_per_target"        # output folder for generated configs
gchp_exec_dir="${Dir}/output-gchp-stretch-soil/Test_Stretch_1month/GEOSChem_build"

cd "${Dir}/${imiDir}"

mkdir -p "$outdir"

# Count data rows (exclude header), robust even if last line has no trailing \n
n_face=$(awk -F, 'NR>1 && NF>=2 {c++} END {print c+0}' "$target_coords_path")
if [[ "$n_face" -eq 0 ]]; then
  echo "No data rows found in $target_coords_path"; exit 1
fi

# Zero-padding width based on number of rows (e.g., 220 -> width=3 => T001..T220)
pad=${#n_face}

# Iterate rows; awk prints: index lon lat
awk -F, '
  NR==1 {next}          # skip header
  NF<2 {next}           # skip incomplete/blank lines
  {
    # strip trailing CR (Windows) from each field
    gsub(/\r/,"",$1); gsub(/\r/,"",$2)

    # trim leading/trailing spaces
    lon=$1; lat=$2
    sub(/^[ \t]+/,"",lon); sub(/[ \t]+$/,"",lon)
    sub(/^[ \t]+/,"",lat); sub(/[ \t]+$/,"",lat)

    print NR-1, lon, lat
  }' "$target_coords_path" \
| while read -r idx lon lat; do
    # Build zero-padded tag, e.g., 1 -> 001
    printf -v tag "%0${pad}d" "$idx"

    # Output filename (you can change naming if you prefer lat/lon in the name)
    outyml="configs_per_target/config-soil_T${tag}.yml"

    # Copy template
    cp "$template" "$outyml"

    # RunName
    sed -i "s/^RunName:.*/RunName: \"Test_Stretch_1day_T${tag}\"/" "$outyml"
    # TARGET_LAT / TARGET_LON
    sed -i "s/^TARGET_LAT:.*/TARGET_LAT: ${lat}/" "$outyml"
    sed -i "s/^TARGET_LON:.*/TARGET_LON: ${lon}/" "$outyml"

    # echo "Created ${outyml}  (lon=${lon}, lat=${lat}, tag=T${tag})"

    mkdir -p "${Dir}/output-gchp-stretch-soil/Test_Stretch_1day_T${tag}"
    ln -nsf "${gchp_exec_dir}" "${Dir}/output-gchp-stretch-soil/Test_Stretch_1day_T${tag}/GEOSChem_build"
    # submit IMI
    sbatch -p serial_requeue,unrestricted,shared \
    -t 0-12:00 --mem 16G -c 1 -o imi_output_T${tag}.log run_imi.sh "$outyml"
done
