#!/usr/bin/env bash
############################################################
## gsMap pipeline
## Step 2) Generate gene specificity scores
##
## This script contains two usage modes:
##   Part A: run one sample directly
##   Part B: run all samples in batch
##
## Usage:
##   bash Step2.Generate.GSS.sh single <sample_name>
##   bash Step2.Generate.GSS.sh batch
##
## Examples:
##   bash Step2.Generate.GSS.sh single MEL131
##   bash Step2.Generate.GSS.sh batch
##
## Author: Shelley
############################################################

set -euo pipefail

############################################################
## Common paths
############################################################
PIPE_BASE="/public/home/hpc8301200407/WGS/GWAS_ST/gsMap_pipeline"
LOG_DIR="${PIPE_BASE}/logs_step2"
mkdir -p "${LOG_DIR}"

############################################################
## Part A) Run step 2 for one sample
############################################################
run_one_sample() {
  local sample_name="$1"
  local working_dir="${PIPE_BASE}/wd_${sample_name}"

  echo "[INFO] Sample:       ${sample_name}"
  echo "[INFO] Working dir:  ${working_dir}"

  gsmap run_latent_to_gene \
    --workdir "${working_dir}" \
    --sample_name "${sample_name}" \
    --latent_representation "latent_GVAE"

  echo "[INFO] Finished step 2: run_latent_to_gene for ${sample_name}"
}

############################################################
## Part B) Run step 2 for all samples
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
    echo "[$(date)] >>> Step 2 for ${s}"
    out_log="${LOG_DIR}/${s}.step2.out"
    err_log="${LOG_DIR}/${s}.step2.err"

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