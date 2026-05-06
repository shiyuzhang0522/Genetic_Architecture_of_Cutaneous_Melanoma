#!/usr/bin/env bash

set -euo pipefail

############################################################
# GCTA-COJO stepwise selection by chromosome
#
# Usage:
#   bash run_cojo.sh <chromosome>
#
# Example:
#   bash run_cojo.sh 1
############################################################

chr="$1"

# Input files
bfile="/path/to/LD_reference_panel/chr${chr}"
sumstats="/path/to/COJO_input/GWAS_chr${chr}.txt"

# Output
outpref="/path/to/output/COJO.gwas.chr${chr}"

# Run COJO stepwise selection
gcta64 \
  --bfile "${bfile}" \
  --chr "${chr}" \
  --threads 32 \
  --cojo-file "${sumstats}" \
  --cojo-slct \
  --diff-freq 0.3 \
  --out "${outpref}"