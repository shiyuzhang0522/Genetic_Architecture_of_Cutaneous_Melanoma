#!/usr/bin/env bash

set -euo pipefail

############################################################
# gsMap pipeline
# Step 3: Generate LD scores
#
# Usage:
#   bash Step3.ldscore.sh run <sample_name> <chr>
#   bash Step3.ldscore.sh run-sample <sample_name>
#   bash Step3.ldscore.sh run-all
#   bash Step3.ldscore.sh check
#
# Examples:
#   bash Step3.ldscore.sh run MEL131 1
#   bash Step3.ldscore.sh run-sample MEL131
#   bash Step3.ldscore.sh run-all
#   bash Step3.ldscore.sh check
############################################################

############################################################
# Common paths
############################################################

BASE_WD="/path/to/gsMap_pipeline"
GTF_FILE="/path/to/gencode.annotation.gtf"
LDREF_ROOT="${BASE_WD}/LD_ref.panel"

N_THREADS=1

############################################################
# Sample list
############################################################

SAMPLES=(
  MEL126
  MEL162
  MEL176
  MEL222
  MEL256
  MEL74
  MEL99
  MEL109
  MEL131
  MEL170
  MEL201
  MEL239
  MEL271
  MEL75
)

############################################################
# Usage
############################################################

usage() {
  cat <<EOF
Usage:
  bash $0 run <sample_name> <chr>
  bash $0 run-sample <sample_name>
  bash $0 run-all
  bash $0 check

Examples:
  bash $0 run MEL131 1
  bash $0 run-sample MEL131
  bash $0 run-all
  bash $0 check
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
# Run one sample-chromosome job
############################################################

run_one_job() {
  local sample_name="$1"
  local chr="$2"
  local working_dir="${BASE_WD}/wd_${sample_name}"
  local bfile_root="${LDREF_ROOT}/chr${chr}.UKBB_LD10k_MAF0.05"

  echo "[INFO] Generating LD score"
  echo "[INFO] Sample:      ${sample_name}"
  echo "[INFO] Chromosome:  ${chr}"
  echo "[INFO] Working dir: ${working_dir}"
  echo "[INFO] LD panel:    ${bfile_root}"
  echo "[INFO] GTF file:    ${GTF_FILE}"

  export OMP_NUM_THREADS="${N_THREADS}"

  gsmap run_generate_ldscore \
    --workdir "${working_dir}" \
    --sample_name "${sample_name}" \
    --chrom "${chr}" \
    --bfile_root "${bfile_root}" \
    --gtf_annotation_file "${GTF_FILE}" \
    --gene_window_size 50000

  echo "[INFO] Finished LD-score generation for ${sample_name}, chr${chr}"
}

############################################################
# Run all chromosomes for one sample
############################################################

run_one_sample() {
  local sample_name="$1"

  for chr in $(seq 1 22); do
    run_one_job "${sample_name}" "${chr}"
  done

  echo "[INFO] Finished all chromosomes for ${sample_name}"
}

############################################################
# Run all samples and chromosomes
############################################################

run_all_jobs() {
  for sample in "${SAMPLES[@]}"; do
    for chr in $(seq 1 22); do
      run_one_job "${sample}" "${chr}"
    done
  done

  echo "[INFO] Finished all sample-chromosome jobs."
}

############################################################
# Check completion
############################################################

check_done_files() {
  echo "[CHECK] Looking for *_generate_ldscore.done files"

  for sample in "${SAMPLES[@]}"; do
    done_file="${BASE_WD}/wd_${sample}/${sample}/generate_ldscore/${sample}_generate_ldscore.done"

    if [[ -f "${done_file}" ]]; then
      echo "[OK]   ${sample}: found $(basename "${done_file}")"
    else
      echo "[MISS] ${sample}: DONE file not found"
    fi
  done
}

############################################################
# Parse mode
############################################################

[[ $# -ge 1 ]] || usage

mode="$1"

activate_env

case "${mode}" in
  run)
    [[ $# -eq 3 ]] || usage
    run_one_job "$2" "$3"
    ;;

  run-sample)
    [[ $# -eq 2 ]] || usage
    run_one_sample "$2"
    ;;

  run-all)
    [[ $# -eq 1 ]] || usage
    run_all_jobs
    ;;

  check)
    [[ $# -eq 1 ]] || usage
    check_done_files
    ;;

  *)
    usage
    ;;
esac