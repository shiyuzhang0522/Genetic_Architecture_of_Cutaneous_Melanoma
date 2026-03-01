#!/usr/bin/env bash
###############################################################################
# Pipeline: LDSC (Liftover hg38->hg19 + h2 estimation)
#
# Purpose:
#   1) Subset SNPs (exclude INDELs) from GWAS summary table
#   2) Convert to BED (hg38), remove unstable hg38 CUPs regions
#   3) LiftOver hg38 -> hg19 (GRCh37)
#   4) Run LDSC h2 on liability and observed scale
#
# Author: Shiyu Zhang (Shelley)
# Date:   2025-12-17
###############################################################################
set -euo pipefail

############################################
# Config
############################################

# Input GWAS summary table (hg38 coordinates assumed here)
IN_GWAS="Meta.GWAS.melanoma1.QCed.txt"

# Working directory for intermediate files
OUTDIR="/public/home/hpc8301200407/WGS/Melanoma_WGS/BHR/LDSC/stable.liftover"

# CUPs BED (hg38) to exclude unstable regions
CUPS_BED="${OUTDIR}/FASTA_BED.ALL_GRCh38.novel_CUPs.bed"

# LiftOver chain (hg38 -> hg19)
CHAIN="/public/home/hpc8301200407/software/liftover/hg38ToHg19.over.chain.gz"

# Tools (edit if not in PATH)
LIFTOVER_BIN="liftOver"
LDSC_PY="/public/home/hpc8301200407/software/ldsc/ldsc.py"

# LD score reference
REF_LD="/public/home/hpc8301200407/WGS/Melanoma_WGS/BHR/LDSC/LDscore/UKBB.EUR"
W_LD="/public/home/hpc8301200407/WGS/Melanoma_WGS/BHR/LDSC/LDscore/UKBB.EUR"

# LDSC input (must be hg19 and formatted for ldsc.py --h2)
# (This file is produced by a separate R/Python formatting step, not included here.)
LDSC_SUMSTATS_HG19="${OUTDIR}/melanoma.GWAS.hg19.LDSC.input.txt"

# Prevalence (edit to your setting)
SAMP_PREV="0.03934525"
POP_PREV="0.0147"

mkdir -p "$OUTDIR"

############################################
# Step 0) Build hg38 SNP-only BED (exclude INDELs)
############################################
# Output: melanoma1.snps.noindel.bed
BED_SNP="${OUTDIR}/melanoma1.snps.noindel.bed"

awk 'BEGIN{FS=OFS="\t"}
NR==1 {next}
{
  chr=$1; pos=$2; ref=$3; alt=$4;

  # keep SNPs only
  if (length(ref)==1 && length(alt)==1) {
    chr="chr"chr;
    start=pos-1;
    end=pos;
    name=chr":"pos":"ref":"alt;  # stable ID for mapping
    print chr, start, end, name
  }
}' "$IN_GWAS" > "$BED_SNP"

echo "[INFO] SNP BED written: $BED_SNP"
echo "[INFO] SNP BED lines  : $(wc -l < "$BED_SNP")"
head -n 3 "$BED_SNP"

############################################
# Step 1) Remove SNPs overlapping hg38 CUPs
############################################
# Output: melanoma1.snps.noindel.noCUPs.bed
test -s "$CUPS_BED"

CUPS_SORTED="${OUTDIR}/FASTA_BED.ALL_GRCh38.novel_CUPs.sorted.bed"
SNP_SORTED="${OUTDIR}/melanoma1.snps.noindel.sorted.bed"
BED_NO_CUPS="${OUTDIR}/melanoma1.snps.noindel.noCUPs.bed"

sort -k1,1 -k2,2n "$CUPS_BED" > "$CUPS_SORTED"
sort -k1,1 -k2,2n "$BED_SNP"   > "$SNP_SORTED"

bedtools intersect -v \
  -a "$SNP_SORTED" \
  -b "$CUPS_SORTED" \
  > "$BED_NO_CUPS"

n_in=$(wc -l < "$SNP_SORTED")
n_out=$(wc -l < "$BED_NO_CUPS")
echo "[INFO] CUPs filter input SNPs : $n_in"
echo "[INFO] CUPs filter kept SNPs  : $n_out"
echo "[INFO] CUPs filter removed    : $(( n_in - n_out ))"
head -n 3 "$BED_NO_CUPS"

############################################
# Step 2) LiftOver hg38 -> hg19
############################################
# Output: melanoma1.snps.noindel.noCUPs.hg19.bed (+ unmapped)
test -s "$CHAIN"
command -v "$LIFTOVER_BIN" >/dev/null 2>&1

BED_NO_CUPS_SORTED="${OUTDIR}/melanoma1.snps.noindel.noCUPs.sorted.bed"
BED_HG19="${OUTDIR}/melanoma1.snps.noindel.noCUPs.hg19.bed"
UNMAPPED="${OUTDIR}/melanoma1.snps.noindel.noCUPs.hg19.unmapped"

sort -k1,1 -k2,2n "$BED_NO_CUPS" > "$BED_NO_CUPS_SORTED"

"$LIFTOVER_BIN" \
  "$BED_NO_CUPS_SORTED" \
  "$CHAIN" \
  "$BED_HG19" \
  "$UNMAPPED"

echo "[INFO] LiftOver output BED : $BED_HG19"
echo "[INFO] LiftOver mapped     : $(wc -l < "$BED_HG19")"
echo "[INFO] LiftOver unmapped   : $(wc -l < "$UNMAPPED")"
head -n 3 "$BED_HG19"

############################################
# Step 3) Run LDSC h2
############################################
# Assumes LDSC_SUMSTATS_HG19 already exists (hg19 + munged/compatible).
test -s "$LDSC_SUMSTATS_HG19"
test -e "$LDSC_PY"

# Liability scale
python "$LDSC_PY" \
  --h2 "$LDSC_SUMSTATS_HG19" \
  --ref-ld "$REF_LD" \
  --w-ld "$W_LD" \
  --samp-prev "$SAMP_PREV" \
  --pop-prev "$POP_PREV" \
  --out "${OUTDIR}/LDSC.h2.melanoma.liability"

# Observed scale
python "$LDSC_PY" \
  --h2 "$LDSC_SUMSTATS_HG19" \
  --ref-ld "$REF_LD" \
  --w-ld "$W_LD" \
  --out "${OUTDIR}/LDSC.h2.melanoma.observed"

echo "[INFO] Done."