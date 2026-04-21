#!/bin/bash
#SBATCH --job-name=clump_meta
#SBATCH --qos=high
#SBATCH --partition=GPU2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --array=1-22
#SBATCH --output=/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/post-GWAS/Clump/log/clump_%A_%a.out
#SBATCH --error=/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/post-GWAS/Clump/log/clump_%A_%a.err

set -euo pipefail

chr="${SLURM_ARRAY_TASK_ID}"

############################################################
## Input files and software
############################################################
ref_dir="/hpc/home/zhangsy/0223_DATA_TRANSFER/WGS_data"
bfile="${ref_dir}/chr${chr}/plink_QCed_chr${chr}"

sumstats="/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/Meta.GWAS.melanoma1.QCed.txt"

outdir="/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/post-GWAS/Clump"
outpref="${outdir}/clumped.meta.GWAS.chr${chr}"

plink2_exec="/hpc/home/zhangsy/software/plink2"

############################################################
## Sanity checks
############################################################
[[ -f "${sumstats}" ]] || { echo "ERROR: summary statistics not found: ${sumstats}" >&2; exit 1; }
[[ -f "${bfile}.bed" ]] || { echo "ERROR: reference BED not found: ${bfile}.bed" >&2; exit 1; }
[[ -f "${bfile}.bim" ]] || { echo "ERROR: reference BIM not found: ${bfile}.bim" >&2; exit 1; }
[[ -f "${bfile}.fam" ]] || { echo "ERROR: reference FAM not found: ${bfile}.fam" >&2; exit 1; }
[[ -x "${plink2_exec}" ]] || { echo "ERROR: plink2 executable not found: ${plink2_exec}" >&2; exit 1; }

mkdir -p "${outdir}"

echo "[$(date)] Starting clumping for chr${chr}"
echo "  Reference panel : ${bfile}"
echo "  Summary stats   : ${sumstats}"
echo "  Threads         : ${SLURM_CPUS_PER_TASK}"

"${plink2_exec}" \
  --bfile "${bfile}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --clump "${sumstats}" \
  --clump-p1 5e-8 \
  --clump-p2 1e-5 \
  --clump-r2 0.1 \
  --clump-kb 1000 \
  --clump-snp-field MarkerName \
  --clump-field P.value \
  --out "${outpref}"

echo "[$(date)] Finished clumping for chr${chr}"
echo "  Output: ${outpref}.clumped"