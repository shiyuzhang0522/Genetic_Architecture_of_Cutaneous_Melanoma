#!/usr/bin/env bash
#SBATCH --job-name=spatial_LDSC
#SBATCH --output=/public/home/hpc8301200407/WGS/GWAS_ST/gsMap_pipeline/logs_step4/%x.%j.out
#SBATCH --error=/public/home/hpc8301200407/WGS/GWAS_ST/gsMap_pipeline/logs_step4/%x.%j.err
#SBATCH --partition=cpuQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --qos=cpuq
#SBATCH --account=pi_dengguangtong

############################################################
## gsMap pipeline
## Step 4) Run spatial LDSC
##
## This script contains two usage modes:
##   run    : run spatial LDSC for one sample
##   submit : submit jobs for all samples
##
## Usage:
##   sbatch Step4.spatial.LDSC.sh run <sample_name>
##   bash   Step4.spatial.LDSC.sh submit
##
## Examples:
##   sbatch Step4.spatial.LDSC.sh run MEL131
##   bash   Step4.spatial.LDSC.sh submit
##
## Author: Shelley
############################################################

set -euo pipefail

############################################################
## Common paths
############################################################
BASE_WD="/public/home/hpc8301200407/WGS/GWAS_ST/gsMap_pipeline"
LOG_DIR="${BASE_WD}/logs_step4"
SUMSTATS_FILE="${BASE_WD}/melanoma_gwas_for_gsmap.tsv"

mkdir -p "${LOG_DIR}"

############################################################
## Sample list
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
## Usage
############################################################
usage() {
  cat <<EOF
Usage:
  sbatch $0 run <sample_name>
  bash   $0 submit

Examples:
  sbatch $0 run MEL131
  bash   $0 submit
EOF
  exit 1
}

############################################################
## Part A) Run spatial LDSC for one sample
############################################################
run_one_sample() {
  local sample_name="$1"
  local working_dir="${BASE_WD}/wd_${sample_name}"
  local w_prefix="${working_dir}/${sample_name}/generate_ldscore/w_ld/weights."
  local num_processes="${SLURM_CPUS_PER_TASK:-16}"

  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate gsMap_new || true
  else
    echo "[WARN] conda command not found; assuming gsmap is on PATH" >&2
  fi

  echo "[INFO] Sample:        ${sample_name}"
  echo "[INFO] Working dir:   ${working_dir}"
  echo "[INFO] Sumstats file: ${SUMSTATS_FILE}"
  echo "[INFO] Weight prefix: ${w_prefix}"
  echo "[INFO] Processes:     ${num_processes}"

  gsmap run_spatial_ldsc \
    --workdir "${working_dir}" \
    --sample_name "${sample_name}" \
    --trait_name "Melanoma" \
    --sumstats_file "${SUMSTATS_FILE}" \
    --w_file "${w_prefix}" \
    --num_processes "${num_processes}"

  echo "[INFO] Finished Step 4: spatial LDSC for ${sample_name}"
}

############################################################
## Part B) Submit jobs for all samples
############################################################
submit_all_jobs() {
  for s in "${SAMPLES[@]}"; do
    echo "[SUBMIT] spatial LDSC for sample: ${s}"
    sbatch "$0" run "${s}"
  done

  echo "[INFO] All Step 4 jobs submitted."
}

############################################################
## Parse mode
############################################################
[[ $# -ge 1 ]] || usage

mode="$1"

case "${mode}" in
  run)
    [[ $# -eq 2 ]] || usage
    run_one_sample "$2"
    ;;
  submit)
    [[ $# -eq 1 ]] || usage
    submit_all_jobs
    ;;
  *)
    usage
    ;;
esac