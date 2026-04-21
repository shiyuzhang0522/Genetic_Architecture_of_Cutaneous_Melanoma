#!/usr/bin/env bash
############################################################
## BHR (burden heritability regression), step-1
## Calculate ALT allele counts in UKBB WGS data on UKB-RAP
##
## Modes:
##   1) Run one job:
##        bash alt_counts_ukbrap.sh run <chr> <case|control>
##
##   2) Submit all chromosomes for both groups:
##        bash alt_counts_ukbrap.sh submit_all
##
## Example:
##   bash alt_counts_ukbrap.sh run 20 case
##   bash alt_counts_ukbrap.sh submit_all
##
## Author: Shelley
############################################################

set -euo pipefail

############################################################
## Step 1) Usage
############################################################
usage() {
  cat <<EOF
Usage:
  bash $0 run <chr> <case|control>
  bash $0 submit_all

Examples:
  bash $0 run 20 case
  bash $0 run 20 control
  bash $0 submit_all
EOF
  exit 1
}

[[ $# -ge 1 ]] || usage
mode="$1"

############################################################
## Step 2) Shared configuration
############################################################

## UKB-RAP PLINK2 pfile prefix
## Keep the escaped quotes because the RAP path contains spaces.
pfile_prefix="/mnt/project/Bulk/\"DRAGEN WGS\"/\"DRAGEN population level WGS variants, PLINK format [500k release]\"/ukb24308_c"

## RAP output folder
output_folder="project-GxVyGzQJ5V1qVf14k938ZFf8:/BHR/counts/allele_counts"

############################################################
## Step 3) Run one job
############################################################
run_one() {
  [[ $# -eq 2 ]] || usage

  local chr="$1"
  local group="$2"

  if ! [[ "${chr}" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
    echo "[ERROR] chr must be 1-22. Got: ${chr}" >&2
    exit 1
  fi

  if [[ "${group}" != "case" && "${group}" != "control" ]]; then
    echo "[ERROR] group must be 'case' or 'control'. Got: ${group}" >&2
    exit 1
  fi

  echo "[INFO] Running for chromosome ${chr}; sample group = ${group}"

  local variant_list="/mnt/project/BHR/counts/variant.list/chr${chr}.variant_list.keep.txt"
  local keep_file

  if [[ "${group}" == "case" ]]; then
    keep_file="/mnt/project/BHR/counts/FID_IID/cases.keep"
  else
    keep_file="/mnt/project/BHR/counts/FID_IID/controls.keep"
  fi

  local output_prefix="chr${chr}.${group}.ALT_counts"

  local cmd
  cmd="plink2 \
    --no-pheno \
    --pfile ${pfile_prefix}${chr}_b0_v1 \
    --keep ${keep_file} \
    --extract ${variant_list} \
    --freq counts \
    --out ${output_prefix}"

  echo "[INFO] Command to run inside swiss-army-knife:"
  echo "${cmd}"

  dx run app-swiss-army-knife \
    -icmd="${cmd}" \
    --folder="${output_folder}" \
    --instance-type=mem2_ssd1_v2_x8 \
    --name="ALTcount_chr${chr}_${group}" \
    --priority=high \
    -y
}

############################################################
## Step 4) Submit all chromosomes for case and control
############################################################
submit_all() {
  for chr in $(seq 1 22); do
    for group in case control; do
      echo "[SUBMIT] chr${chr} ${group}"
      run_one "${chr}" "${group}"
    done
  done

  echo "[INFO] Submitted chr1-22 for both case and control"
}

############################################################
## Step 5) Dispatch
############################################################
case "${mode}" in
  run)
    [[ $# -eq 3 ]] || usage
    run_one "$2" "$3"
    ;;
  submit_all)
    [[ $# -eq 1 ]] || usage
    submit_all
    ;;
  *)
    usage
    ;;
esac