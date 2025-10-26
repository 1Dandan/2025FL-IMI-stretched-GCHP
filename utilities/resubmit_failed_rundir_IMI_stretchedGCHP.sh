#!/bin/bash
set -euo pipefail

# -------------------------
# User-configurable section
# -------------------------
Dir="/n/holylfs06/LABS/jacob_lab2/Lab/dzhang8/imi-gchp-test"
imiDir="imi-stretch-gchp-cp2"
outdir="${Dir}/${imiDir}/configs_per_target"
gchp_exec_dir="${Dir}/output-gchp-stretch-soil/Test_Stretch_1month/GEOSChem_build"
failed_list="${Dir}/${imiDir}/failed_runs.txt"

# Slurm knobs
SBATCH_PARTITIONS="sapphire,huce_cascade,seas_compute,shared,unrestricted"
SBATCH_TIME="0-12:00"
SBATCH_MEM="16G"
SBATCH_CPUS="1"
IMI_LAUNCHER="run_imi.sh"

# -------------------------
# Helper
# -------------------------

tweak_yaml_for_redo() {
  local yml="$1"
  [[ -f "$yml" ]] || { log "ERROR: YAML not found: $yml"; return 1; }

  sed -i -E \
    -e 's/^RunSetup:.*/RunSetup: true/' \
    -e 's/^SetupTemplateRundir:.*/SetupTemplateRundir: false/' \
    -e 's/^SetupSpinupRun:.*/SetupSpinupRun: true/' \
    -e 's/^SetupJacobianRuns:.*/SetupJacobianRuns: true/' \
    -e 's/^SetupInversion:.*/SetupInversion: true/' \
    -e 's/^SetupPosteriorRun:.*/SetupPosteriorRun: false/' \
    -e 's/^DoHemcoPriorEmis:.*/DoHemcoPriorEmis: true/' \
    -e 's/^DoSpinup:.*/DoSpinup: true/' \
    -e 's/^DoJacobian:.*/DoJacobian: true/' \
    -e 's/^ReDoJacobian:.*/ReDoJacobian: true/' \
    -e 's/^DoInversion:.*/DoInversion: true/' \
    -e 's/^DoPosterior:.*/DoPosterior: false/' \
    -e 's/^DoPreview:.*/DoPreview: true/' \
    "$yml"
}

submit_job() {
  local tag="$1" yml="$2"
  local logf="imi_output_${tag}.log"

  sbatch -p "${SBATCH_PARTITIONS}" \
         -t "${SBATCH_TIME}" \
         --mem "${SBATCH_MEM}" \
         -c "${SBATCH_CPUS}" \
         -o "${logf}" \
         "${IMI_LAUNCHER}" "${yml}"
}

# -------------------------
# Main
# -------------------------


[[ -f "$failed_list" ]] || { echo "No failed list found: $failed_list"; exit 1; }

cd ${Dir}/${imiDir}

while IFS= read -r tag; do
    # Skip blanks or comments
    [[ -z "$tag" || "$tag" =~ ^# ]] && continue

    yml="configs_per_target/config-soil_${tag}.yml"
    logf="imi_output_${tag}.log"

    if [[ ! -f "$yml" ]]; then
        echo "Missing YAML for ${tag}: $yml"
        continue
    fi

    tweak_yaml_for_redo "$yml"

    echo "Resubmitting ${tag}"
    submit_job "$tag" "$yml"
    
done < "$failed_list"