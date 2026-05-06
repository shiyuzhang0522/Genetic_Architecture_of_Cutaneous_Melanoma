#!/usr/bin/env bash

set -euo pipefail

############################################################
# gsMap pipeline
# Step 4: Run spatial LDSC
#
# Usage:
#   bash Step4.spatial.LDSC.sh run <sample_name>
#   bash Step4.spatial.LDSC.sh run-all
#
# Examples:
#   bash Step4.spatial.LDSC.sh run MEL131
#   bash Step4.spatial.LDSC.sh run-all
############################################################

############################################################
# Common paths
############################################################

BASE_WD="/path/to/gsMap_pipeline"
SUMSTATS_FILE="${BASE_WD}/melanoma_gwas_for_gsmap.tsv"

TRAIT_NAME="Melanoma"
N_PROCESSES=16

############################################################
# Sample list
############################################################

SAMPLES=(
  MEL109
  MEL126
  MEL131
  MEL162
  MEL170
  MEL176
  MEL201
  MEL222
  MEL239
  MEL256
  MEL271
  MEL74
  MEL75
  MEL99
)

############################################################
# Usage
############################################################

usage() {
  cat <<EOF
Usage:
  bash $0 run <sample_name>
  bash $0 run-all

Examples:
  bash $0 run MEL131
  bash $0 run-all
EOF
  exit 1
}

############################################################
# Activate environment
############################################################

activate_env() {
  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate gsMap_new
  else
    echo "[WARN] conda command not found; assuming gsmap is already on PATH" >&2
  fi
}

############################################################
# Run spatial LDSC for one sample
############################################################

run_one_sample() {
  local sample_name="$1"
  local working_dir="${BASE_WD}/wd_${sample_name}"
  local w_prefix="${working_dir}/${sample_name}/generate_ldscore/w_ld/weights."

  echo "[INFO] Running spatial LDSC"
  echo "[INFO] Sample:        ${sample_name}"
  echo "[INFO] Working dir:   ${working_dir}"
  echo "[INFO] Sumstats file: ${SUMSTATS_FILE}"
  echo "[INFO] Weight prefix: ${w_prefix}"
  echo "[INFO] Processes:     ${N_PROCESSES}"

  gsmap run_spatial_ldsc \
    --workdir "${working_dir}" \
    --sample_name "${sample_name}" \
    --trait_name "${TRAIT_NAME}" \
    --sumstats_file "${SUMSTATS_FILE}" \
    --w_file "${w_prefix}" \
    --num_processes "${N_PROCESSES}"

  echo "[INFO] Finished spatial LDSC for ${sample_name}"
}

############################################################
# Run spatial LDSC for all samples
############################################################

run_all_samples() {
  for sample in "${SAMPLES[@]}"; do
    run_one_sample "${sample}"
  done

  echo "[INFO] Finished spatial LDSC for all samples."
}

############################################################
# Parse mode
############################################################

[[ $# -ge 1 ]] || usage

mode="$1"

activate_env

case "${mode}" in
  run)
    [[ $# -eq 2 ]] || usage
    run_one_sample "$2"
    ;;

  run-all)
    [[ $# -eq 1 ]] || usage
    run_all_samples
    ;;

  *)
    usage
    ;;
esac