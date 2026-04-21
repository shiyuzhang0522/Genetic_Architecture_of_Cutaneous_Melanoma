#!/usr/bin/env bash
#SBATCH --job-name=cojo_meta
#SBATCH --qos=high
#SBATCH --partition=GPU1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --array=1-22
#SBATCH --output=/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/post-GWAS/COJO/log/cojo_%A_%a.out
#SBATCH --error=/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/post-GWAS/COJO/log/cojo_%A_%a.err

set -euo pipefail

chr="${SLURM_ARRAY_TASK_ID}"

############################################################
## Paths
############################################################
base_dir="/hpc/home/zhangsy/Meta-GWAS/GWAS_summary/Meta_analysis/post-GWAS/COJO"

bfile="/hpc/home/zhangsy/COJO_LD_Panel/chr${chr}_LD_panel_for_COJO"
sumstats="${base_dir}/GWAS_input/GWAS_meta_chr${chr}.txt"
outdir="${base_dir}/results"

gcta_exec="/hpc/home/zhangsy/software/gcta-1.94.1-linux-kernel-3-x86_64/gcta64"

mkdir -p "${outdir}" "${outdir}/log"

outpref="${outdir}/COJO.meta.GWAS.chr${chr}"

############################################################
## Sanity checks
############################################################
[[ -f "${sumstats}" ]] || { echo "[ERROR] Missing sumstats: ${sumstats}" >&2; exit 1; }

[[ -f "${bfile}.bed" && -f "${bfile}.bim" && -f "${bfile}.fam" ]] || {
  echo "[ERROR] Missing PLINK files for bfile prefix: ${bfile}" >&2
  exit 1
}

[[ -x "${gcta_exec}" ]] || {
  echo "[ERROR] GCTA executable not found: ${gcta_exec}" >&2
  exit 1
}

############################################################
## Run COJO (stepwise selection)
############################################################
echo "[$(date)] Starting COJO on chr${chr}"
echo "  LD reference : ${bfile}"
echo "  Sumstats     : ${sumstats}"
echo "  Threads      : ${SLURM_CPUS_PER_TASK}"
echo "  Output       : ${outpref}"

"${gcta_exec}" \
  --bfile "${bfile}" \
  --chr "${chr}" \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --cojo-file "${sumstats}" \
  --cojo-slct \
  --diff-freq 0.3 \
  --out "${outpref}"

echo "[$(date)] Finished COJO on chr${chr}"