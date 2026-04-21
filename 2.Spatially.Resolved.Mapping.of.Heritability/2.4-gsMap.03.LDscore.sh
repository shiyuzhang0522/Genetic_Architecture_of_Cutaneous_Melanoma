#!/usr/bin/env bash
#SBATCH --job-name=ldscore
#SBATCH --output=/public/home/hpc8301200407/WGS/GWAS_ST/gsMap_pipeline/logs_step3/%x.%j.out
#SBATCH --error=/public/home/hpc8301200407/WGS/GWAS_ST/gsMap_pipeline/logs_step3/%x.%j.err
#SBATCH --partition=cpuQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --qos=cpuq
#SBATCH --account=pi_dengguangtong

############################################################
## gsMap pipeline
## Step 3) Generate LD scores
##
## This script contains three usage modes:
##   run     : run one sample-chromosome job
##   submit  : submit all sample × chromosome jobs
##   check   : check whether per-sample done files exist
##
## Usage:
##   sbatch Step3.ldscore.sh run <sample_name> <chr>
##   bash   Step3.ldscore.sh submit
##   bash   Step3.ldscore.sh check
##
## Examples:
##   sbatch Step3.ldscore.sh run MEL131
##   bash   Step3.ldscore.sh submit
##   bash   Step3.ldscore.sh check
##
## Author: Shelley
############################################################

set -euo pipefail

############################################################
## Common paths
############################################################
BASE_WD="/public/home/hpc8301200407/WGS/GWAS_ST/gsMap_pipeline"
LOG_DIR="${BASE_WD}/logs_step3"
GTF_FILE="/public/home/hpc8301200407/software/gsMap_resource/gencode.v49.basic.annotation.gtf"
LDREF_ROOT="${BASE_WD}/LD_ref.panel"

mkdir -p "${LOG_DIR}"

############################################################
## Sample list
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
## Usage
############################################################
usage() {
  cat <<EOF
Usage:
  sbatch $0 run <sample_name> <chr>
  bash   $0 submit
  bash   $0 check

Examples:
  sbatch $0 run MEL131
  bash   $0 submit
  bash   $0 check
EOF
  exit 1
}

############################################################
## Part A) Run one sample-chromosome job
############################################################
run_one_job() {
  local sample_name="$1"
  local chr="$2"
  local working_dir="${BASE_WD}/wd_${sample_name}"

  echo "[INFO] Generating LD score for sample ${sample_name}, chr${chr}..."
  echo "[INFO] Sample:       ${sample_name}"
  echo "[INFO] Chromosome:   ${chr}"
  echo "[INFO] Working dir:  ${working_dir}"

  if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate gsMap_new || true
  else
    echo "[WARN] conda command not found; assuming gsmap is on PATH" >&2
  fi

  export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-1}"

  gsmap run_generate_ldscore \
    --workdir "${working_dir}" \
    --sample_name "${sample_name}" \
    --chrom "${chr}" \
    --bfile_root "${LDREF_ROOT}/chr${chr}.UKBB_LD10k_MAF0.05" \
    --gtf_annotation_file "${GTF_FILE}" \
    --gene_window_size 50000

  echo "[INFO] Finished step 3: run_generate_ldscore for ${sample_name} chr${chr}"
}

############################################################
## Part B) Submit all jobs
############################################################
submit_all_jobs() {
  mapfile -t IDLE_NODES < <(sinfo -h -N -p cpuQ -t idle -o '%N')

  if [[ ${#IDLE_NODES[@]} -eq 0 ]]; then
    echo "[ERROR] No idle nodes found in partition cpuQ" >&2
    exit 1
  fi

  echo "[INFO] Idle nodes:"
  printf '  %s\n' "${IDLE_NODES[@]}"

  local job_idx=0

  for sample in "${SAMPLES[@]}"; do
    for chr in $(seq 1 22); do
      local node="${IDLE_NODES[$((job_idx % ${#IDLE_NODES[@]}))]}"

      echo "[SUBMIT] sample=${sample}, chr=${chr} -> node=${node}"

      sbatch \
        --nodelist="${node}" \
        "$0" run "${sample}" "${chr}"

      ((job_idx+=1))
    done
  done

  echo "[INFO] All submissions done."
}

############################################################
## Part C) Check completion
############################################################
check_done_files() {
  echo "[CHECK] Looking for *_generate_ldscore.done files"

  for s in "${SAMPLES[@]}"; do
    done_file="${BASE_WD}/wd_${s}/${s}/generate_ldscore/${s}_generate_ldscore.done"
    if [[ -f "${done_file}" ]]; then
      echo "[OK]   ${s}: found $(basename "${done_file}")"
    else
      echo "[MISS] ${s}: DONE file not found"
    fi
  done
}

############################################################
## Parse mode
############################################################
[[ $# -ge 1 ]] || usage

mode="$1"

case "${mode}" in
  run)
    [[ $# -eq 3 ]] || usage
    run_one_job "$2" "$3"
    ;;
  submit)
    [[ $# -eq 1 ]] || usage
    submit_all_jobs
    ;;
  check)
    [[ $# -eq 1 ]] || usage
    check_done_files
    ;;
  *)
    usage
    ;;
esac