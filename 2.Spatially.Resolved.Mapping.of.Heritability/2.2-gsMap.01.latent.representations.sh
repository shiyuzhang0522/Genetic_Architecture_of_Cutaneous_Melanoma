#!/usr/bin/env bash
############################################################
## gsMap pipeline
## Step 1) Find latent representations
##
## This script contains two usage modes:
##   Part A: run one sample directly
##   Part B: run all samples in batch
##
## Usage:
##   bash Step1.latent.representations.sh single <sample_name>
##   bash Step1.latent.representations.sh batch
##
## Examples:
##   bash Step1.latent.representations.sh single MEL131
##   bash Step1.latent.representations.sh batch
##
## Author: Shelley
############################################################

set -euo pipefail

############################################################
## Common paths
############################################################
PIPE_BASE="/public/home/hpc8301200407/WGS/GWAS_ST/gsMap_pipeline"
LOG_DIR="${PIPE_BASE}/logs_step1"
mkdir -p "${LOG_DIR}"

############################################################
## Part A) Run latent representation step for one sample
############################################################
run_one_sample() {
  local sample_name="$1"
  local working_dir="${PIPE_BASE}/wd_${sample_name}"
  local h5ad_file="${PIPE_BASE}/ST_h5ad/${sample_name}.h5ad"

  mkdir -p "${working_dir}"

  if [ ! -f "${h5ad_file}" ]; then
    echo "[ERROR] h5ad file not found: ${h5ad_file}" >&2
    return 1
  fi

  echo "[INFO] Sample:       ${sample_name}"
  echo "[INFO] Working dir:  ${working_dir}"
  echo "[INFO] Input h5ad:   ${h5ad_file}"

  gsmap run_find_latent_representations \
    --workdir "${working_dir}" \
    --sample_name "${sample_name}" \
    --input_hdf5_path "${h5ad_file}" \
    --data_layer "count" \
    --epochs 1000

  echo "[INFO] Finished run_find_latent_representations for ${sample_name}"
}

############################################################
## Part B) Run Step 1 for all samples
############################################################
run_batch() {
  local samples=(
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

  for s in "${samples[@]}"; do
    echo "[$(date)] >>> Step 1 for ${s}"
    out_log="${LOG_DIR}/${s}.step1.out"
    err_log="${LOG_DIR}/${s}.step1.err"

    bash "$0" single "${s}" > "${out_log}" 2> "${err_log}"

    echo "[$(date)] >>> Finished ${s}"
    echo "[INFO] Logs: ${out_log}, ${err_log}"
    echo
  done
}

############################################################
## Parse usage
############################################################
usage() {
  cat <<EOF
Usage:
  bash $0 single <sample_name>
  bash $0 batch

Examples:
  bash $0 single MEL131
  bash $0 batch
EOF
  exit 1
}

[ $# -ge 1 ] || usage

mode="$1"

if [ "${mode}" = "single" ]; then
  [ $# -eq 2 ] || usage
  run_one_sample "$2"
elif [ "${mode}" = "batch" ]; then
  run_batch
else
  usage
fi