#!/usr/bin/env bash
############################################################
## Prepare LDSC input and estimate SNP heritability
##
## Modules:
##   Step 1) Extract autosomal SNPs and convert to BED
##   Step 2) Remove hg38 conversion-unstable positions (CUPs)
##   Step 3) Lift SNP coordinates from hg38 to hg19
##   Step 4) Run LDSC heritability analysis
##
## Notes:
##   - This script prepares the SNP coordinate set only.
##   - Formatting the final LDSC summary-statistics input file
##     is handled separately in an R script.
##
## Author: Shelley
## Date:  2025-12-17
############################################################

set -euo pipefail

############################################################
## Step 0) Define inputs and outputs
############################################################
in_sumstats="Meta.GWAS.melanoma1.QCed.txt"
outdir="/public/home/hpc8301200407/WGS/Melanoma_WGS/BHR/LDSC/stable.liftover"

cups_bed="${outdir}/FASTA_BED.ALL_GRCh38.novel_CUPs.bed"
chain="/public/home/hpc8301200407/software/liftover/hg38ToHg19.over.chain.gz"

ldsc_py="/public/home/hpc8301200407/software/ldsc/ldsc.py"
ref_ld="/public/home/hpc8301200407/WGS/Melanoma_WGS/BHR/LDSC/LDscore/UKBB.EUR"
w_ld="/public/home/hpc8301200407/WGS/Melanoma_WGS/BHR/LDSC/LDscore/UKBB.EUR"

ldsc_input="${outdir}/melanoma.GWAS.hg19.LDSC.input.txt"

mkdir -p "${outdir}"

############################################################
## Step 1) Extract SNPs only and convert to BED
############################################################
awk 'BEGIN{FS=OFS="\t"}
NR==1 {next}
{
  chr=$1; pos=$2; ref=$3; alt=$4;
  if (length(ref)==1 && length(alt)==1) {
    chr="chr"chr;
    start=pos-1;
    end=pos;
    name=chr":"pos":"ref":"alt;
    print chr, start, end, name
  }
}' "${in_sumstats}" > "${outdir}/melanoma1.snps.noindel.bed"

echo "[INFO] BED written: ${outdir}/melanoma1.snps.noindel.bed"
echo "[INFO] Number of SNPs extracted: $(wc -l < "${outdir}/melanoma1.snps.noindel.bed")"
head -n 3 "${outdir}/melanoma1.snps.noindel.bed"

############################################################
## Step 2) Remove hg38 CUPs
############################################################
## CUPs = conversion-unstable positions
############################################################
sort -k1,1 -k2,2n "${cups_bed}" \
  > "${outdir}/FASTA_BED.ALL_GRCh38.novel_CUPs.sorted.bed"

sort -k1,1 -k2,2n "${outdir}/melanoma1.snps.noindel.bed" \
  > "${outdir}/melanoma1.snps.noindel.sorted.bed"

bedtools intersect -v \
  -a "${outdir}/melanoma1.snps.noindel.sorted.bed" \
  -b "${outdir}/FASTA_BED.ALL_GRCh38.novel_CUPs.sorted.bed" \
  > "${outdir}/melanoma1.snps.noindel.noCUPs.bed"

echo "[INFO] Input SNPs : $(wc -l < "${outdir}/melanoma1.snps.noindel.sorted.bed")"
echo "[INFO] Kept SNPs  : $(wc -l < "${outdir}/melanoma1.snps.noindel.noCUPs.bed")"
echo "[INFO] Removed    : $(( $(wc -l < "${outdir}/melanoma1.snps.noindel.sorted.bed") - $(wc -l < "${outdir}/melanoma1.snps.noindel.noCUPs.bed") ))"
head -n 3 "${outdir}/melanoma1.snps.noindel.noCUPs.bed"

############################################################
## Step 3) Lift hg38 SNP coordinates to hg19
############################################################
in_bed="${outdir}/melanoma1.snps.noindel.noCUPs.bed"
sorted_bed="${outdir}/melanoma1.snps.noindel.noCUPs.sorted.bed"
out_bed="${outdir}/melanoma1.snps.noindel.noCUPs.hg19.bed"
unmapped="${outdir}/melanoma1.snps.noindel.noCUPs.hg19.unmapped"

test -s "${in_bed}"
test -s "${chain}"

sort -k1,1 -k2,2n "${in_bed}" > "${sorted_bed}"

liftOver \
  "${sorted_bed}" \
  "${chain}" \
  "${out_bed}" \
  "${unmapped}"

echo "[INFO] LiftOver complete"
echo "[INFO] Lifted BED   : ${out_bed}"
echo "[INFO] Unmapped BED : ${unmapped}"
echo "[INFO] Successfully lifted: $(wc -l < "${out_bed}")"
echo "[INFO] Unmapped variants  : $(grep -vc '^#' "${unmapped}" || true)"

############################################################
## Step 4) Run LDSC heritability analysis
############################################################
## The LDSC-formatted summary-statistics file is assumed to
## have already been prepared separately:
##   ${ldsc_input}
############################################################

test -s "${ldsc_input}"
test -f "${ldsc_py}"

## liability-scale heritability
python "${ldsc_py}" \
  --h2 "${ldsc_input}" \
  --ref-ld "${ref_ld}" \
  --w-ld "${w_ld}" \
  --samp-prev 0.03934525 \
  --pop-prev 0.0147 \
  --out "${outdir}/LDSC.h2.melanoma.liability"

## observed-scale heritability
python "${ldsc_py}" \
  --h2 "${ldsc_input}" \
  --ref-ld "${ref_ld}" \
  --w-ld "${w_ld}" \
  --out "${outdir}/LDSC.h2.melanoma.observed"

echo "[INFO] LDSC finished"
echo "[INFO] Liability-scale output: ${outdir}/LDSC.h2.melanoma.liability"
echo "[INFO] Observed-scale output : ${outdir}/LDSC.h2.melanoma.observed"