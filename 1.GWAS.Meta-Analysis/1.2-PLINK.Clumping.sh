#!/bin/bash

set -euo pipefail

############################################################
# GWAS LD clumping using PLINK2
#
# Usage:
#   bash clump_gwas.sh <chromosome>
#
# Example:
#   bash clump_gwas.sh 1
############################################################

chr="$1"

# Input files
bfile="/path/to/reference_panel/chr${chr}"
sumstats="/path/to/gwas_summary_statistics.txt"

# Output
outpref="/path/to/output/clumped.gwas.chr${chr}"

# Run LD clumping
plink2 \
  --bfile "${bfile}" \
  --threads 16 \
  --clump "${sumstats}" \
  --clump-p1 5e-8 \
  --clump-p2 1e-5 \
  --clump-r2 0.1 \
  --clump-kb 1000 \
  --clump-snp-field MarkerName \
  --clump-field P.value \
  --out "${outpref}"