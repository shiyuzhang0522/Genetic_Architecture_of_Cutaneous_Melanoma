############## R code to run BHR (aggregate only)
############## Purpose: Aggregate pLoF heritability across freq bins
##############          with fixed genes = set 2 (45 validated genes)
############## Author: Shelley
############## Date: 2025-12-18 (updated 2026-01-06) 
############## Aggregate pLoF heritability across freq bins:
##############   ultra_rare_MAC_le_10 + MAC_gt_10_ALT_AF_le_1e-4

suppressPackageStartupMessages({
  library(data.table)
  library(bhr)
})

############################
# Step 1) fixed inputs
############################
ref_file <- "/public/home/hpc8301200407/software/BHR/bhr/reference_files/ms_baseline_oe5.txt"
sum_dir  <- "/public/home/hpc8301200407/WGS/Melanoma_WGS/BHR/pipeline/input/runs_by_freqbin"

## gene set 2 annotation (45 validated genes)
gene_set_dir <- "/public/home/hpc8301200407/WGS/Melanoma_WGS/BHR/pipeline/Gene.Set.Enrichment.Analysis/Gene.Set.Annotations"
gs2_file <- file.path(gene_set_dir, "gene_set_2_validated45.txt")

stopifnot(file.exists(ref_file))
stopifnot(dir.exists(sum_dir))
stopifnot(file.exists(gs2_file))

baseline_model <- read.table(ref_file)
stopifnot("gene" %in% names(baseline_model))

## read gene set 2 and extract fixed genes (membership == 1)
gs2_dt <- fread(gs2_file)
stopifnot("gene" %in% names(gs2_dt))

mem_cols <- setdiff(names(gs2_dt), "gene")
if (length(mem_cols) != 1) {
  stop("[ERROR] Expected exactly 1 membership column in ", gs2_file,
       " but found: ", paste(mem_cols, collapse = ", "))
}
mem_col <- mem_cols[1]

fixed_genes_set2 <- gs2_dt[get(mem_col) == 1, gene]
fixed_genes_set2 <- fixed_genes_set2[!is.na(fixed_genes_set2) & fixed_genes_set2 != ""]
fixed_genes_set2 <- sub("\\.\\d+$", "", fixed_genes_set2)  # drop ENSG version suffix if present
fixed_genes_set2 <- unique(fixed_genes_set2)

cat("[INFO] Loaded fixed genes (set2_validated45):", length(fixed_genes_set2), "\n")

############################
# Step 2) define aggregate inputs
############################
freq_bins <- c("ultra_rare_MAC_le_10", "MAC_gt_10_ALT_AF_le_1e-4")
mask <- "pLoF"

sum_files <- file.path(sum_dir, freq_bins, paste0("bhr_input_", freq_bins, "_", mask, ".txt"))
names(sum_files) <- freq_bins

missing <- sum_files[!file.exists(sum_files)]
if (length(missing) > 0) {
  stop("[ERROR] Missing input files:\n", paste(missing, collapse = "\n"))
}

cat("[INFO] Aggregate freq bins:\n")
print(sum_files)

## read each freq-bin sumstats into a list of data.frames
ss_list <- lapply(sum_files, function(f) as.data.frame(fread(f)))

############################
# Step 3) run BHR in aggregate mode (with fixed genes)
############################
agg_fit_plof <- BHR(
  mode        = "aggregate",
  ss_list     = ss_list,
  trait_list  = list("melanoma"),
  annotations = list(baseline_model),
  fixed_genes = fixed_genes_set2
)

####### print the results #######
str(agg_fit_plof)
# List of 2
# $ aggregated_mixed_model_h2  : num 0.00341
# $ aggregated_mixed_model_h2se: num 0.000943

############################
# liability conversion + save (aggregate output)
############################

## observed-scale h2 from aggregate fit
h2_obs    <- agg_fit_plof$aggregated_mixed_model_h2
h2_obs_se <- agg_fit_plof$aggregated_mixed_model_h2se

## ---- fill these ----
K <- 0.0147
n_cases <- 3816
n_controls <- 339903
## -------------------

P <- n_cases / (n_cases + n_controls)

## Lee 2011 Eq.23-style scaling factor
t <- qnorm(1 - K)
Z <- dnorm(t)
scale_factor <- (K * (1 - K) / (Z^2)) * (K * (1 - K) / (P * (1 - P)))

## liability-scale
h2_liab    <- h2_obs * scale_factor
h2_liab_se <- h2_obs_se * scale_factor

## quick print
cat("[RESULT] Aggregate pLoF (fixed set2) observed h2 =", h2_obs, " (", h2_obs_se, ")\n", sep = "")
cat("[RESULT] Aggregate pLoF (fixed set2) liability h2 =", h2_liab, " (", h2_liab_se, ")\n", sep = "")

## small results table (nice for downstream reporting)
res_dt <- data.table(
  analysis      = "BHR_aggregate_fixed_set2",
  fixed_set     = "gene_set_2_validated45",
  n_fixed_genes = length(fixed_genes_set2),
  mask          = mask,
  freq_bins     = paste(freq_bins, collapse = " + "),
  K             = K,
  P             = P,
  scale_factor  = scale_factor,
  h2_obs        = h2_obs,
  h2_obs_se     = h2_obs_se,
  h2_liab       = h2_liab,
  h2_liab_se    = h2_liab_se
)

############################
# Save
############################
output_dir <- "/public/home/hpc8301200407/WGS/Melanoma_WGS/BHR/updated.pipeline.Fixed.Genes/results/aggregated"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

txt_file <- file.path(output_dir, "bhr_aggregate_pLoF_fixed_set2_validated45.h2.txt")
fwrite(res_dt, file = txt_file, sep = "\t", quote = FALSE)
cat("[INFO] Saved -> ", txt_file, "\n", sep = "")