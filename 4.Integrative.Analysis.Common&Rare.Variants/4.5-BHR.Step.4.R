############################################################
## BHR - Step3. Aggregate rare-variant heritability estimate
## Run aggregate BHR for pLoF variants with fixed genes
##
## Purpose:
##   Estimate aggregate rare-variant heritability across
##   selected frequency bins, with predefined fixed-gene
##   sets included in the BHR model.
##
## Fixed-effect gene sets:
##   - gene_set_2_validated45
##
## Aggregate frequency bins:
##   - ultra_rare_MAC_le_10
##   - MAC_gt_10_ALT_AF_le_1e-4
##
## Author: Shelley
## Date:   2025-12-18
## Updated: 2026-01-06
############################################################

############################################################
## Step 1) Load packages
############################################################
suppressPackageStartupMessages({
  library(data.table)
  library(bhr)
})

############################################################
## Step 2) Define inputs
############################################################
ref_file     <- "/path/BHR/bhr/reference_files/ms_baseline_oe5.txt"
sum_dir      <- "/path/WGS/Melanoma_WGS/BHR/pipeline/input/runs_by_freqbin"
gene_set_dir <- "/path/WGS/Melanoma_WGS/BHR/pipeline/Gene.Set.Enrichment.Analysis/Gene.Set.Annotations"
output_dir   <- "/path/WGS/Melanoma_WGS/BHR/updated.pipeline.Fixed.Genes/results/aggregated"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(ref_file))
stopifnot(dir.exists(sum_dir))
stopifnot(dir.exists(gene_set_dir))

baseline_model <- read.table(ref_file)
stopifnot("gene" %in% names(baseline_model))

############################################################
## Step 3) Define analysis settings
############################################################
freq_bins <- c(
  "ultra_rare_MAC_le_10",
  "MAC_gt_10_ALT_AF_le_1e-4"
)

mask <- "pLoF"

gene_set_files <- c(
  gene_set_2_validated45 = file.path(gene_set_dir, "gene_set_2_validated45.txt")
)

stopifnot(all(file.exists(gene_set_files)))

sum_files <- file.path(sum_dir, freq_bins, paste0("bhr_input_", freq_bins, "_", mask, ".txt"))
names(sum_files) <- freq_bins

missing <- sum_files[!file.exists(sum_files)]
if (length(missing) > 0) {
  stop("[ERROR] Missing input files:\n", paste(missing, collapse = "\n"))
}

cat("[INFO] Aggregate freq bins:\n")
print(sum_files)

############################################################
## Step 4) Helper functions
############################################################

############################
## Step 4.1) Read one fixed-gene set
############################
read_fixed_genes <- function(gs_file) {
  dt <- fread(gs_file)
  stopifnot("gene" %in% names(dt))

  mem_cols <- setdiff(names(dt), "gene")
  if (length(mem_cols) != 1) {
    stop("[ERROR] Expected exactly 1 membership column in ", gs_file,
         " but found: ", paste(mem_cols, collapse = ", "))
  }

  mem_col <- mem_cols[1]

  genes <- dt[get(mem_col) == 1, gene]
  genes <- genes[!is.na(genes) & genes != ""]
  genes <- sub("\\.\\d+$", "", genes)

  unique(genes)
}

############################
## Step 4.2) Liability conversion factor
############################
obs2lia_factor <- function(K, P) {
  t <- qnorm(1 - K)
  Z <- dnorm(t)
  (K * (1 - K) / (Z^2)) * (K * (1 - K) / (P * (1 - P)))
}

############################
## Step 4.3) Run one aggregate BHR analysis
############################
run_aggregate_bhr <- function(fixed_set_name, gs_file) {
  cat("\n=================================================\n")
  cat("[INFO] Running aggregate BHR for: ", fixed_set_name, "\n", sep = "")

  fixed_genes <- read_fixed_genes(gs_file)
  cat("[INFO] Loaded fixed genes: ", length(fixed_genes), "\n", sep = "")

  ss_list <- lapply(sum_files, function(f) as.data.frame(fread(f)))

  fit <- BHR(
    mode        = "aggregate",
    ss_list     = ss_list,
    trait_list  = list("melanoma"),
    annotations = list(baseline_model),
    fixed_genes = fixed_genes
  )

  str(fit)

  h2_obs    <- fit$aggregated_mixed_model_h2
  h2_obs_se <- fit$aggregated_mixed_model_h2se

  K <- 0.0147
  n_cases <- 3816
  n_controls <- 339903
  P <- n_cases / (n_cases + n_controls)

  scale_factor <- obs2lia_factor(K, P)

  h2_liab    <- h2_obs * scale_factor
  h2_liab_se <- h2_obs_se * scale_factor

  cat("[RESULT] ", fixed_set_name, " observed h2 = ", h2_obs,
      " (", h2_obs_se, ")\n", sep = "")
  cat("[RESULT] ", fixed_set_name, " liability h2 = ", h2_liab,
      " (", h2_liab_se, ")\n", sep = "")

  res_dt <- data.table(
    analysis      = paste0("BHR_aggregate_", fixed_set_name),
    fixed_set     = fixed_set_name,
    n_fixed_genes = length(fixed_genes),
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

  txt_file <- file.path(output_dir, paste0("bhr_aggregate_", fixed_set_name, ".h2.txt"))
  fwrite(res_dt, file = txt_file, sep = "\t", quote = FALSE)

  cat("[INFO] Saved -> ", txt_file, "\n", sep = "")

  list(
    fit = fit,
    summary = res_dt
  )
}

############################################################
## Step 5) Run all requested fixed-gene sets
############################################################
agg_results <- lapply(names(gene_set_files), function(nm) {
  run_aggregate_bhr(
    fixed_set_name = nm,
    gs_file = gene_set_files[[nm]]
  )
})
names(agg_results) <- names(gene_set_files)

############################################################
## Step 6) Save combined summary table
############################################################
summary_dt <- rbindlist(
  lapply(agg_results, function(x) x$summary),
  use.names = TRUE,
  fill = TRUE
)

summary_file <- file.path(output_dir, "bhr_aggregate_fixed_gene_sets.summary.txt")
fwrite(summary_dt, file = summary_file, sep = "\t", quote = FALSE)

cat("\n[INFO] Combined summary saved -> ", summary_file, "\n", sep = "")
cat("[DONE] Aggregate BHR analyses completed.\n")