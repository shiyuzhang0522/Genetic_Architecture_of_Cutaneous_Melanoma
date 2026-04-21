#!/usr/bin/env bash
#SBATCH --job-name=PRS-CS
#SBATCH --output=/public/home/hpc8301200407/WGS/Melanoma_WGS/PRS/My.PRS/PRScs/logs/%x_%A_%a.out
#SBATCH --error=/public/home/hpc8301200407/WGS/Melanoma_WGS/PRS/My.PRS/PRScs/logs/%x_%A_%a.err
#SBATCH --partition=fatQ
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --qos=fatq
#SBATCH --account=pi_dengguangtong
#SBATCH --array=1-22

############################################################
## Run PRS-CS by chromosome
##
## Purpose:
##   Estimate posterior SNP effect sizes using PRS-CS with
##   the UK Biobank European LD reference panel.
##
## Important note:
##   The GWAS summary-statistics input file must be lifted
##   over to hg19/GRCh37 before running PRS-CS, to match the
##   genome build of the PRS-CS LD reference panel.
##
## Input:
##   - GWAS.PRSCS.input.txt   (hg19/GRCh37)
##   - snpinfo_ukbb_hm3.bim
##   - UKBB LD blocks for EUR
##
## Output:
##   - Chromosome-specific PRS-CS posterior effect estimates
##
## Author: Shelley
## Reference:https://github.com/getian107/PRScs
############################################################

set -euo pipefail

############################################################
## Step 1) Get chromosome index from SLURM array
############################################################
chr="${SLURM_ARRAY_TASK_ID}"

############################################################
## Step 2) Define input and output paths
############################################################
ref_dir="/public/home/hpc8301200407/software/PRScs/LD.Ref.Panel/UKBB/ldblk_ukbb_eur"
bim_prefix="/public/home/hpc8301200407/WGS/Melanoma_WGS/PRS/My.PRS/PRScs/snpinfo_ukbb_hm3"
sst_file="/public/home/hpc8301200407/WGS/Melanoma_WGS/PRS/My.PRS/PRScs/GWAS.PRSCS.input.txt"
out_dir="/public/home/hpc8301200407/WGS/Melanoma_WGS/PRS/My.PRS/PRScs/outputs/chr${chr}"

mkdir -p "${out_dir}"

############################################################
## Step 3) Activate environment
############################################################
eval "$(~/miniconda-fresh/condabin/conda shell.bash hook)"
conda activate PRS.CS

############################################################
## Step 4) Print run information
############################################################
echo "[INFO] Running PRS-CS for chr ${chr} ..."
echo "[INFO] ref_dir    = ${ref_dir}"
echo "[INFO] bim_prefix = ${bim_prefix}"
echo "[INFO] sst_file   = ${sst_file}"
echo "[INFO] out_dir    = ${out_dir}"
echo "[INFO] NOTE: sst_file must already be lifted over to hg19/GRCh37."

############################################################
## Step 5) Set thread environment variables
############################################################
threads="${SLURM_CPUS_PER_TASK:-8}"
export MKL_NUM_THREADS="${threads}"
export NUMEXPR_NUM_THREADS="${threads}"
export OMP_NUM_THREADS="${threads}"

############################################################
## Step 6) Sanity checks
############################################################
[[ -d "${ref_dir}" ]] || { echo "[ERROR] ref_dir not found: ${ref_dir}" >&2; exit 1; }
[[ -f "${sst_file}" ]] || { echo "[ERROR] sst_file not found: ${sst_file}" >&2; exit 1; }
[[ -f "${bim_prefix}.bim" ]] || { echo "[ERROR] BIM not found: ${bim_prefix}.bim" >&2; exit 1; }

############################################################
## Step 7) Run PRS-CS
############################################################
python /public/home/hpc8301200407/software/PRScs/PRScs.py \
  --ref_dir="${ref_dir}" \
  --bim_prefix="${bim_prefix}" \
  --sst_file="${sst_file}" \
  --n_gwas=200679 \
  --out_dir="${out_dir}" \
  --n_iter=10000 \
  --n_burnin=5000 \
  --chrom="${chr}" \
  --seed=1220

echo "[INFO] PRS-CS successfully finished for chr ${chr}!"