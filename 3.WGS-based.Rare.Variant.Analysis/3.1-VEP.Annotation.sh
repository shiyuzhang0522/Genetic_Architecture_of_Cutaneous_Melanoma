#!/usr/bin/env bash
#SBATCH --job-name=vep_anno
#SBATCH --partition=cpuQ
#SBATCH --qos=cpuq
#SBATCH --account=pi_dengguangtong
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=190464
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

############################################################
## VEP annotation of WGS variant lists (site-only VCF)
## GRCh38 with LoFTEE and dbNSFP under Singularity
##
## Description:
##   Annotates site-level variants (no sample genotypes)
##   using Ensembl VEP with LoFTEE (LoF prediction) and
##   dbNSFP (functional scores such as REVEL and CADD).
##
## Input:
##   - chr*.vcf (site-only)
##
## Usage:
##   sbatch vep_annotate.sh <chr>
##
## Author: Shelley
############################################################

set -euo pipefail

############################################################
## Step 1) Environment
############################################################
module purge
module load singularity

############################################################
## Step 2) Parse command-line argument
############################################################
if [[ $# -ne 1 ]]; then
  echo "ERROR: Expected 1 argument: <chromosome_number>" >&2
  echo "Usage: $0 <chromosome_number>" >&2
  exit 1
fi

CHR="$1"
if ! [[ "${CHR}" =~ ^([1-9]|1[0-9]|2[0-2])$ ]]; then
  echo "ERROR: chromosome_number must be an integer in [1..22]. Got: ${CHR}" >&2
  exit 1
fi

############################################################
## Step 3) Define paths
############################################################
SCRATCH_DIR="/public/home/hpc8301200407/Shelley"

INPUT_DIR="${SCRATCH_DIR}/plink2"
OUTPUT_DIR="${SCRATCH_DIR}/variant.annotation/VEP_new"
LOG_DIR="${OUTPUT_DIR}/log"

CACHE_DIR="/public/home/hpc8301200407/tool/ensembl-vep/cache.annotations"
SIF_IMAGE="/public/home/hpc8301200407/tool/vep.sif"

REF_FASTA_DIR="${OUTPUT_DIR}/FASTA"
REF_FASTA="${REF_FASTA_DIR}/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
SAMTOOLS_BIN="${HOME}/samtools/bin/samtools"

############################################################
## Step 4) Define input and output files
############################################################
INPUT_VCF="${INPUT_DIR}/chr${CHR}.copied.pvar.vcf"
OUTPUT_VCF="${OUTPUT_DIR}/chr${CHR}.VEP.anno.vcf"

mkdir -p "${OUTPUT_DIR}" "${LOG_DIR}" "logs"

############################################################
## Step 5) Sanity checks
############################################################
[[ -f "${INPUT_VCF}" ]] || { echo "ERROR: Input VCF not found: ${INPUT_VCF}" >&2; exit 1; }
[[ -f "${SIF_IMAGE}" ]] || { echo "ERROR: Singularity image not found: ${SIF_IMAGE}" >&2; exit 1; }
[[ -d "${CACHE_DIR}" ]] || { echo "ERROR: VEP cache directory not found: ${CACHE_DIR}" >&2; exit 1; }
[[ -f "${REF_FASTA}" ]] || { echo "ERROR: Reference FASTA not found: ${REF_FASTA}" >&2; exit 1; }
[[ -f "${REF_FASTA}.fai" ]] || { echo "ERROR: Missing FASTA index: ${REF_FASTA}.fai" >&2; exit 1; }

GZI_PATH="${REF_FASTA%.gz}.gz.gzi"
[[ -f "${GZI_PATH}" ]] || { echo "ERROR: Missing FASTA gzi index: ${GZI_PATH}" >&2; exit 1; }

############################################################
## Step 6) Build VEP command
############################################################
LOFTEE_DIR="${CACHE_DIR}/loftee"

VEP_CMD=(
  vep
  -i "${INPUT_VCF}"
  --assembly GRCh38
  --vcf
  --format vcf
  --cache
  --dir_cache /cache
  --fasta /ref.fa.gz
  -o "${OUTPUT_VCF}"
  --plugin "LoF,loftee_path:/plugins,\
human_ancestor_fa:/plugins/vep_data/human_ancestor.fa.gz,\
conservation_file:/plugins/vep_data/loftee.sql,\
gerp_bigwig:/plugins/vep_data/gerp_conservation_scores.homo_sapiens.GRCh38.bw"
  --plugin "dbNSFP,/plugins/dbNSFP5.1a.txt.gz,REVEL_score,CADD_phred"
  --everything
  --offline
  --force_overwrite
)

############################################################
## Step 7) Build Singularity command
############################################################
SING_CMD=(
  singularity exec --cleanenv
  -W "${HOME}"
  --bind "${SCRATCH_DIR}"
  --bind "${CACHE_DIR}":/cache
  --bind "${LOFTEE_DIR}":/plugins
  --bind "${REF_FASTA}":/ref.fa.gz
  --bind "${REF_FASTA}.fai":/ref.fa.gz.fai
  --bind "${GZI_PATH}":/ref.fa.gz.gzi
)

if [[ -x "${SAMTOOLS_BIN}" ]]; then
  SING_CMD+=( --bind "${SAMTOOLS_BIN}":/usr/local/bin/samtools )
fi

SING_CMD+=( "${SIF_IMAGE}" )

############################################################
## Step 8) Run annotation
############################################################
echo "[INFO] CHR=${CHR}"
echo "[INFO] INPUT_VCF=${INPUT_VCF}"
echo "[INFO] OUTPUT_VCF=${OUTPUT_VCF}"
echo "[INFO] CACHE_DIR=${CACHE_DIR}"
echo "[INFO] REF_FASTA=${REF_FASTA}"
echo

echo "[INFO] Running command:"
printf "  %q " "${SING_CMD[@]}" "${VEP_CMD[@]}"
echo
echo

"${SING_CMD[@]}" "${VEP_CMD[@]}"

echo "[INFO] Done: ${OUTPUT_VCF}"