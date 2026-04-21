############################################################
## Harmonize and summarize prioritized protein-coding genes
## from GWAS and rare-variant collapsing analysis
## for the melanoma biologic network
##
## Purpose:
##   - Read gene-prioritization results from:
##       (1) GWAS gene prioritization (FUMA SNP2GENE)
##       (2) Rare-variant protein-coding associations
##   - Harmonize the two sources into a unified gene table
##   - Restrict to protein-coding genes
##   - Check for missing or inconsistent gene symbols
##   - Summarize source-level evidence for each unique gene
##
## Inputs:
##   - FUMA genes.txt
##   - protein_coding.FDR0.05.associations.txt
##
## Outputs:
##   - all_genes_long_pc_dt.txt
##   - gene_evidence_dt.txt
##   - missing_symbol_dt.txt
##   - multi_symbol_dt.txt
##
## Author: Shelley Zhang
## Date:   2026-03-06
############################################################

############################################################
## Step 1) Load package
############################################################
suppressPackageStartupMessages({
  library(data.table)
})

cat("[INFO] Package loaded: data.table\n")

############################################################
## Step 2) Define input and output paths
############################################################
fuma_file <- "/Users/shelleyz/Projects/Melanoma-WGS/FUMA_input/FUMA_job674132-2_SNP2GENE/genes.txt"
rare_file <- "/Users/shelleyz/Projects/Melanoma-WGS/Updated_RareVariantAnalysis/FDR005_Associations/protein_coding.FDR0.05.associations.txt"

out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/Updated.Omnigenic.Architecture/Vertices"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(fuma_file))
stopifnot(file.exists(rare_file))

############################################################
## Step 3) Read input data
############################################################

############################
## Step 3.1) FUMA genes
############################
fuma_dt <- fread(fuma_file)

cat("[INFO] Loaded FUMA genes: ", fuma_file, "\n", sep = "")
cat("[INFO] fuma_dt: ", nrow(fuma_dt), " rows x ", ncol(fuma_dt), " cols\n", sep = "")
print(head(fuma_dt, 3))

############################
## Step 3.2) Rare-variant genes
############################
rare_dt <- fread(rare_file)

cat("[INFO] Loaded rare-variant genes: ", rare_file, "\n", sep = "")
cat("[INFO] rare_dt: ", nrow(rare_dt), " rows x ", ncol(rare_dt), " cols\n", sep = "")
print(head(rare_dt, 3))

############################################################
## Step 4) Harmonize each source to a common schema
## Common columns:
##   - ensg
##   - symbol
##   - type
##   - source
############################################################

############################
## Step 4.1) GWAS genes (FUMA)
############################
stopifnot(all(c("ensg", "symbol", "type") %in% names(fuma_dt)))

fuma_prior_dt <- fuma_dt[, .(
  ensg   = as.character(ensg),
  symbol = as.character(symbol),
  type   = as.character(type),
  source = "GWAS"
)]

cat("[INFO] fuma_prior_dt: ", nrow(fuma_prior_dt), " rows x ", ncol(fuma_prior_dt), " cols\n", sep = "")
print(head(fuma_prior_dt, 5))

############################
## Step 4.2) Rare-variant protein-coding genes
############################
stopifnot(all(c("Region", "gene_symbol") %in% names(rare_dt)))

rare_prior_dt <- rare_dt[, .(
  ensg   = as.character(Region),
  symbol = as.character(gene_symbol)
)]

rare_prior_dt[, type := "protein_coding"]
rare_prior_dt[, source := "RV_PC"]

## drop blank symbols
rare_prior_dt <- rare_prior_dt[!is.na(symbol) & trimws(symbol) != ""]

## deduplicate
rare_prior_dt <- unique(rare_prior_dt)

cat("[INFO] rare_prior_dt: ", nrow(rare_prior_dt), " rows x ", ncol(rare_prior_dt), " cols\n", sep = "")
print(head(rare_prior_dt, 5))

############################################################
## Step 5) Combine sources and restrict to protein-coding
############################################################
all_genes_long_dt <- rbindlist(
  list(
    fuma_prior_dt[, .(ensg, symbol, type, source)],
    rare_prior_dt[, .(ensg, symbol, type, source)]
  ),
  use.names = TRUE,
  fill = TRUE
)

all_genes_long_pc_dt <- all_genes_long_dt[type == "protein_coding"]

cat("[INFO] all_genes_long_pc_dt: ", nrow(all_genes_long_pc_dt), " rows x ",
    ncol(all_genes_long_pc_dt), " cols\n", sep = "")
cat("[INFO] Unique protein-coding ENSG: ", uniqueN(all_genes_long_pc_dt$ensg), "\n", sep = "")
print(head(all_genes_long_pc_dt, 5))

############################################################
## Step 6) Sanity checks
##   - missing symbols
##   - ENSG IDs mapping to multiple symbols
############################################################

############################
## Step 6.1) Missing symbols
############################
missing_symbol_dt <- all_genes_long_pc_dt[
  is.na(symbol) | trimws(symbol) == "",
  .(ensg, symbol, type, source)
][order(ensg, source)]

cat("[INFO] Rows with missing symbol: ", nrow(missing_symbol_dt), "\n", sep = "")
cat("[INFO] Unique ENSG with missing symbol: ", uniqueN(missing_symbol_dt$ensg), "\n", sep = "")
print(missing_symbol_dt)

############################
## Step 6.2) ENSG with multiple distinct symbols
############################
tmp_dt <- copy(all_genes_long_pc_dt)
tmp_dt[, symbol_norm := toupper(trimws(symbol))]

multi_symbol_ensg <- tmp_dt[
  !is.na(symbol_norm) & symbol_norm != "",
  .(n_symbol = uniqueN(symbol_norm)),
  by = ensg
][n_symbol > 1, ensg]

multi_symbol_dt <- tmp_dt[
  ensg %in% multi_symbol_ensg,
  .(ensg, symbol, type, source)
][order(ensg, source, symbol)]

cat("[INFO] ENSG with >1 distinct symbols: ", length(multi_symbol_ensg), "\n", sep = "")
print(multi_symbol_dt)

############################################################
## Step 7) Save long-format protein-coding gene table and QC
############################################################
out_pc_file      <- file.path(out_dir, "all_genes_long_pc_dt.txt")
out_missing_file <- file.path(out_dir, "missing_symbol_dt.txt")
out_multi_file   <- file.path(out_dir, "multi_symbol_dt.txt")

fwrite(all_genes_long_pc_dt, out_pc_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(missing_symbol_dt,    out_missing_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(multi_symbol_dt,      out_multi_file, sep = "\t", quote = FALSE, na = "NA")

cat("[INFO] Saved long-format and QC tables to: ", out_dir, "\n", sep = "")

############################################################
## Step 8) Collapse to one row per gene and summarize
## evidence across sources
############################################################

############################
## Step 8.1) Helper to choose a representative symbol
############################
pick_first_nonmissing <- function(x) {
  x <- x[!is.na(x) & trimws(x) != ""]
  if (length(x) == 0) return(NA_character_)
  x[1]
}

dt0 <- copy(all_genes_long_pc_dt)

gene_evidence_dt <- dt0[
  ,
  .(
    symbol = pick_first_nonmissing(symbol),
    type = pick_first_nonmissing(type),
    evidence_types = paste(sort(unique(source)), collapse = "+"),
    evidence_n = uniqueN(source)
  ),
  by = ensg
]

setorder(gene_evidence_dt, -evidence_n, evidence_types, symbol)

gene_evidence_dt[, evidence_class := fifelse(
  evidence_n == 1,
  evidence_types,
  paste0("multi_", evidence_n)
)]

cat("[INFO] gene_evidence_dt: ", nrow(gene_evidence_dt), " rows x ",
    ncol(gene_evidence_dt), " cols\n", sep = "")
print(head(gene_evidence_dt, 10))

cat("[INFO] Evidence pattern counts:\n")
print(
  gene_evidence_dt[
    ,
    .N,
    by = .(evidence_n, evidence_types)
  ][order(-evidence_n, -N)]
)

############################################################
## Step 9) Save collapsed evidence table
############################################################
out_evidence_file <- file.path(out_dir, "gene_evidence_dt.txt")

fwrite(
  gene_evidence_dt,
  out_evidence_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("[INFO] Saved gene_evidence_dt -> ", out_evidence_file, "\n", sep = "")
cat("[INFO] Finished.\n")