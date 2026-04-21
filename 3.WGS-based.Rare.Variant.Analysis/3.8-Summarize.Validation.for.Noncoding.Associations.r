############################################################
## Summarize and visualize GEL replication results for
## noncoding rare-variant associations
## (ncRNA and melanocyte-specific cCREs)
##
## Discovery cohort:
##   - UK Biobank (UKBB), WGS
##
## Replication cohort:
##   - Genomics England (GEL), WGS
##
## Purpose:
##   - Aggregate GEL replication results for ncRNA and cCRE units
##     aligned to UKBB discovery signals
##   - Generate clean summary tables and figures
##
## Author: Shelley
## Date:   2026-02-24
############################################################

############################################################
## Step 1) Load packages and define inputs
############################################################
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

UKBB_cCRE_file  <- "/Users/shelleyz/Projects/Melanoma-WGS/Updated_RareVariantAnalysis/FDR005_Associations/cCRE.FDR0.05.associations.txt"
UKBB_ncRNA_file <- "/Users/shelleyz/Projects/Melanoma-WGS/Updated_RareVariantAnalysis/FDR005_Associations/ncRNA.FDR0.05.associations.txt"

GEL_cCRE_file   <- "/Users/shelleyz/Projects/Melanoma-WGS/Validation/100kGP-GEL/Results/Firth_cCRE_all_masks_results.tsv"
GEL_ncRNA_file  <- "/Users/shelleyz/Projects/Melanoma-WGS/Validation/100kGP-GEL/Results/Firth_ncRNA_all_masks_results.tsv"

stopifnot(file.exists(UKBB_cCRE_file))
stopifnot(file.exists(UKBB_ncRNA_file))
stopifnot(file.exists(GEL_cCRE_file))
stopifnot(file.exists(GEL_ncRNA_file))

cat("[INFO] Step 1 complete: all input files exist.\n")

############################################################
## Step 2) Read data
############################################################
UKBB_cCRE_dt  <- fread(UKBB_cCRE_file)
UKBB_ncRNA_dt <- fread(UKBB_ncRNA_file)
GEL_cCRE_dt   <- fread(GEL_cCRE_file)
GEL_ncRNA_dt  <- fread(GEL_ncRNA_file)

cat("[INFO] UKBB_cCRE_dt : ", nrow(UKBB_cCRE_dt),  " rows x ", ncol(UKBB_cCRE_dt),  " cols\n", sep = "")
cat("[INFO] UKBB_ncRNA_dt: ", nrow(UKBB_ncRNA_dt), " rows x ", ncol(UKBB_ncRNA_dt), " cols\n", sep = "")
cat("[INFO] GEL_cCRE_dt  : ", nrow(GEL_cCRE_dt),   " rows x ", ncol(GEL_cCRE_dt),   " cols\n", sep = "")
cat("[INFO] GEL_ncRNA_dt : ", nrow(GEL_ncRNA_dt),  " rows x ", ncol(GEL_ncRNA_dt),  " cols\n", sep = "")

str(UKBB_cCRE_dt)
str(UKBB_ncRNA_dt)
str(GEL_cCRE_dt)
str(GEL_ncRNA_dt)

############################################################
## Step 3) Harmonize UKBB discovery tables
############################################################

############################
## Step 3.1) UKBB cCRE
## Region format: chr12_PLS_EH38E1615598
############################
stopifnot(all(c(
  "Region", "Group", "max_MAF",
  "Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
  "BETA_Burden", "SE_Burden", "FDR"
) %in% names(UKBB_cCRE_dt)))

UKBB_cCRE_h <- copy(UKBB_cCRE_dt)

UKBB_cCRE_h[
  , c("chr", "cCRE_type", "accession_id") := tstrsplit(Region, "_", fixed = TRUE)
]

UKBB_cCRE_h[, chr := as.integer(sub("^chr", "", chr))]

UKBB_cCRE_h <- UKBB_cCRE_h[
  , .(
    Region,
    chr,
    cCRE_type,
    accession_id,
    Group,
    max_MAF,
    Pvalue,
    Pvalue_Burden,
    Pvalue_SKAT,
    BETA_Burden,
    SE_Burden,
    FDR
  )
]

cat("[INFO] UKBB_cCRE_h created: ", nrow(UKBB_cCRE_h), " rows x ", ncol(UKBB_cCRE_h), " cols\n", sep = "")

############################
## Step 3.2) UKBB ncRNA
## Region format: RP11-460I13.2_ENSG00000227050_antisense
############################
stopifnot(all(c(
  "Region", "Group", "max_MAF",
  "Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
  "BETA_Burden", "SE_Burden", "FDR"
) %in% names(UKBB_ncRNA_dt)))

UKBB_ncRNA_h <- copy(UKBB_ncRNA_dt)

UKBB_ncRNA_h[
  , c("gene_symbol", "ensembl_id", "ncRNA_type") := tstrsplit(Region, "_", fixed = TRUE)
]

UKBB_ncRNA_h <- UKBB_ncRNA_h[
  , .(
    Region,
    gene_symbol,
    ensembl_id,
    ncRNA_type,
    Group,
    max_MAF,
    Pvalue,
    Pvalue_Burden,
    Pvalue_SKAT,
    BETA_Burden,
    SE_Burden,
    FDR
  )
]

cat("[INFO] UKBB_ncRNA_h created: ", nrow(UKBB_ncRNA_h), " rows x ", ncol(UKBB_ncRNA_h), " cols\n", sep = "")

############################################################
## Step 4) Harmonize GEL replication tables
############################################################

############################
## Step 4.1) GEL cCRE
############################
stopifnot(all(c("predictor", "coef", "se", "lower95", "upper95", "Chisq", "p") %in% names(GEL_cCRE_dt)))

GEL_cCRE_h <- copy(GEL_cCRE_dt)

GEL_cCRE_h[, c("region_part", "mask_part") := tstrsplit(predictor, "__", fixed = TRUE)]

GEL_cCRE_h[, chr := sub("^chr([0-9]+|X|Y)_.*$", "\\1", region_part)]
GEL_cCRE_h[, chr := fifelse(chr %in% c("X", "Y"), chr, as.character(as.integer(chr)))]

GEL_cCRE_h[, accession_id := sub("^chr([0-9]+|X|Y)_([^_]+)_.*$", "\\2", region_part)]
GEL_cCRE_h[, cCRE_type    := sub("^chr([0-9]+|X|Y)_[^_]+_(.*)$", "\\2", region_part)]

GEL_cCRE_h <- GEL_cCRE_h[
  , .(
    predictor,
    chr,
    accession_id,
    cCRE_type,
    mask_part,
    coef,
    se,
    lower95,
    upper95,
    Chisq,
    p
  )
]

cat("[INFO] GEL_cCRE_h created: ", nrow(GEL_cCRE_h), " rows x ", ncol(GEL_cCRE_h), " cols\n", sep = "")

############################
## Step 4.2) GEL ncRNA
############################
stopifnot(all(c("predictor", "coef", "se", "lower95", "upper95", "Chisq", "p") %in% names(GEL_ncRNA_dt)))

GEL_ncRNA_h <- copy(GEL_ncRNA_dt)
GEL_ncRNA_h[, c("region_part", "mask_part") := tstrsplit(predictor, "__", fixed = TRUE)]

GEL_ncRNA_h[, chr := sub("^chr([0-9]+|X|Y)_.*$", "\\1", region_part)]
GEL_ncRNA_h[, chr := fifelse(chr %in% c("X", "Y"), chr, as.character(as.integer(chr)))]

GEL_ncRNA_h[, gene_symbol := sub("^chr([0-9]+|X|Y)_", "", region_part)]

GEL_ncRNA_h <- GEL_ncRNA_h[
  , .(
    predictor,
    chr,
    gene_symbol,
    mask_part,
    coef,
    se,
    lower95,
    upper95,
    Chisq,
    p
  )
]

cat("[INFO] GEL_ncRNA_h created: ", nrow(GEL_ncRNA_h), " rows x ", ncol(GEL_ncRNA_h), " cols\n", sep = "")

############################################################
## Step 5) Select the best GEL mask per unit
############################################################

############################
## Step 5.1) cCRE: best mask per accession_id
############################
stopifnot(all(c("accession_id", "p") %in% names(GEL_cCRE_h)))

GEL_cCRE_best <- GEL_cCRE_h[
  !is.na(p) & is.finite(p)
][
  order(accession_id, p)
][
  , .SD[1], by = accession_id
]

cat("[INFO] GEL_cCRE_best created: ", nrow(GEL_cCRE_best),
    " rows (unique accession_id: ", uniqueN(GEL_cCRE_h$accession_id), ")\n", sep = "")

############################
## Step 5.2) ncRNA: best mask per gene_symbol
############################
stopifnot(all(c("gene_symbol", "p") %in% names(GEL_ncRNA_h)))

GEL_ncRNA_best <- GEL_ncRNA_h[
  !is.na(p) & is.finite(p)
][
  order(gene_symbol, p)
][
  , .SD[1], by = gene_symbol
]

cat("[INFO] GEL_ncRNA_best created: ", nrow(GEL_ncRNA_best),
    " rows (unique gene_symbol: ", uniqueN(GEL_ncRNA_h$gene_symbol), ")\n", sep = "")

############################################################
## Step 6) Restrict GEL best hits to units present in UKBB discovery
############################################################

############################
## Step 6.1) cCRE
############################
ukbb_cCRE_ids <- unique(UKBB_cCRE_h$accession_id)
GEL_cCRE_best_inUKBB <- GEL_cCRE_best[accession_id %in% ukbb_cCRE_ids]

cat("[INFO] GEL_cCRE_best total: ", nrow(GEL_cCRE_best), "\n", sep = "")
cat("[INFO] UKBB cCRE unique accession_id: ", length(ukbb_cCRE_ids), "\n", sep = "")
cat("[INFO] GEL_cCRE_best_inUKBB: ", nrow(GEL_cCRE_best_inUKBB), "\n", sep = "")

############################
## Step 6.2) ncRNA
############################
ukbb_ncRNA_syms <- unique(UKBB_ncRNA_h$gene_symbol)
GEL_ncRNA_best_inUKBB <- GEL_ncRNA_best[gene_symbol %in% ukbb_ncRNA_syms]

cat("[INFO] GEL_ncRNA_best total: ", nrow(GEL_ncRNA_best), "\n", sep = "")
cat("[INFO] UKBB ncRNA unique gene_symbol: ", length(ukbb_ncRNA_syms), "\n", sep = "")
cat("[INFO] GEL_ncRNA_best_inUKBB: ", nrow(GEL_ncRNA_best_inUKBB), "\n", sep = "")

############################################################
## Step 7) Annotate GEL best hits with UKBB discovery statistics
############################################################

############################
## Step 7.1) cCRE: join by accession_id
############################
ukbb_cCRE_anno <- unique(
  UKBB_cCRE_h[
    , .(
      accession_id,
      UKBB_Group         = Group,
      UKBB_max_MAF       = max_MAF,
      UKBB_Pvalue        = Pvalue,
      UKBB_Pvalue_Burden = Pvalue_Burden,
      UKBB_Pvalue_SKAT   = Pvalue_SKAT,
      UKBB_BETA_Burden   = BETA_Burden,
      UKBB_SE_Burden     = SE_Burden,
      UKBB_FDR           = FDR
    )
  ],
  by = c("accession_id", "UKBB_Group", "UKBB_max_MAF")
)

GEL_cCRE_best_annot <- merge(
  GEL_cCRE_best_inUKBB,
  ukbb_cCRE_anno,
  by = "accession_id",
  all.x = TRUE,
  allow.cartesian = TRUE
)

cat("[INFO] GEL_cCRE_best_annot: ", nrow(GEL_cCRE_best_annot), " rows x ",
    ncol(GEL_cCRE_best_annot), " cols\n", sep = "")

############################
## Step 7.2) ncRNA: join by gene_symbol
############################
ukbb_ncRNA_anno <- unique(
  UKBB_ncRNA_h[
    , .(
      gene_symbol,
      UKBB_ensembl_id    = ensembl_id,
      UKBB_ncRNA_type    = ncRNA_type,
      UKBB_Group         = Group,
      UKBB_max_MAF       = max_MAF,
      UKBB_Pvalue        = Pvalue,
      UKBB_Pvalue_Burden = Pvalue_Burden,
      UKBB_Pvalue_SKAT   = Pvalue_SKAT,
      UKBB_BETA_Burden   = BETA_Burden,
      UKBB_SE_Burden     = SE_Burden,
      UKBB_FDR           = FDR
    )
  ],
  by = c("gene_symbol", "UKBB_Group", "UKBB_max_MAF")
)

GEL_ncRNA_best_annot <- merge(
  GEL_ncRNA_best_inUKBB,
  ukbb_ncRNA_anno,
  by = "gene_symbol",
  all.x = TRUE,
  allow.cartesian = TRUE
)

cat("[INFO] GEL_ncRNA_best_annot: ", nrow(GEL_ncRNA_best_annot), " rows x ",
    ncol(GEL_ncRNA_best_annot), " cols\n", sep = "")

############################################################
## Step 8) Save annotated GEL replication tables
############################################################
out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/Validation/Updated_Validation_20260224"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_cCRE <- file.path(out_dir, "Validation.GEL_cCRE_best_annot.tsv")
fwrite(GEL_cCRE_best_annot, out_cCRE, sep = "\t")
cat("[INFO] Saved GEL_cCRE_best_annot:\n  ", out_cCRE, "\n", sep = "")

out_ncRNA <- file.path(out_dir, "Validation.GEL_ncRNA_best_annot.tsv")
fwrite(GEL_ncRNA_best_annot, out_ncRNA, sep = "\t")
cat("[INFO] Saved GEL_ncRNA_best_annot:\n  ", out_ncRNA, "\n", sep = "")