############################################################
## BHR - Step3
## Run BHR with significant genes treated as fixed effects
##
## Purpose:
##   Estimate rare-variant heritability with predefined
##   significant-gene sets included as fixed effects in BHR.
##
## Fixed-effect gene sets:
##   Set 1: 128 significant genes from UK Biobank (FDR < 0.05)
##   Set 2: 45 genes validated in >=1 validation cohort
##   Set 3: 10 high-confidence genes validated in >=2 cohorts
##   Set 4: Mendelian melanoma genes (OCA2, CDKN2A, MITF, MC1R)
##
## Author: Shelley
## Date: 2026-01-06
############################################################

############################################################
## Step 1) Load packages
############################################################
suppressPackageStartupMessages({
  library(data.table)
  library(bhr)
  library(stringr)
})

cat("[INFO] Packages loaded.\n")

############################################################
## Step 2) Define inputs
############################################################
ref_file     <- "/path/BHR/bhr/reference_files/ms_baseline_oe5.txt"
sum_dir      <- "/path/WGS/Melanoma_WGS/BHR/pipeline/input/runs_by_freqbin"
gene_set_dir <- "/path/WGS/Melanoma_WGS/BHR/pipeline/Gene.Set.Enrichment.Analysis/Gene.Set.Annotations"
out_dir      <- "/path/WGS/Melanoma_WGS/BHR/updated.pipeline.Fixed.Genes/results"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(ref_file))
stopifnot(dir.exists(sum_dir))
stopifnot(dir.exists(gene_set_dir))

############################################################
## Step 3) Read baseline reference and fixed-effect gene sets
############################################################

############################
## Step 3.1) Baseline reference
############################
ref_dt <- read.table(ref_file)
stopifnot("gene" %in% names(ref_dt))

cat("[INFO] Loaded baseline reference: ", nrow(ref_dt), " rows x ", ncol(ref_dt), " cols\n", sep = "")

############################
## Step 3.2) Fixed-effect gene sets
############################
gs_files <- c(
  gene_set_1_UKBB128     = file.path(gene_set_dir, "gene_set_1_UKBB128.txt"),
  gene_set_2_validated45 = file.path(gene_set_dir, "gene_set_2_validated45.txt"),
  gene_set_3_highconf10  = file.path(gene_set_dir, "gene_set_3_highconf10.txt"),
  gene_set_4_core4       = file.path(gene_set_dir, "gene_set_4_core4.txt")
)
stopifnot(all(file.exists(gs_files)))

fixed_gene_sets <- lapply(names(gs_files), function(nm) {
  dt <- fread(gs_files[[nm]])
  stopifnot("gene" %in% names(dt))

  mem_cols <- setdiff(names(dt), "gene")
  if (length(mem_cols) != 1) {
    stop("[ERROR] Expected exactly 1 membership column in ", nm,
         " but found: ", paste(mem_cols, collapse = ", "))
  }

  mem_col <- mem_cols[1]
  genes <- dt[get(mem_col) == 1, gene]
  genes <- genes[!is.na(genes) & genes != ""]
  genes <- sub("\\.\\d+$", "", genes)

  unique(genes)
})
names(fixed_gene_sets) <- names(gs_files)

cat("[INFO] Fixed gene set sizes:\n")
print(sapply(fixed_gene_sets, length))

############################################################
## Step 4) Define BHR run grid
############################################################
freq_bins <- c(
  "ultra_rare_MAC_le_10",
  "MAC_gt_10_ALT_AF_le_1e-4",
  "ALT_AF_gt_1e-4_to_le_1e-3"
)

masks <- c(
  "pLoF",
  "damaging_missense_or_protein_altering",
  "synonymous"
)

fixed_set_names <- names(fixed_gene_sets)

############################################################
## Step 5) Run BHR across fixed-gene sets, frequency bins, and masks
############################################################
bhr_results <- list()

for (fs in fixed_set_names) {
  fixed_genes_now <- fixed_gene_sets[[fs]]

  cat("\n=================================================\n")
  cat("[INFO] Fixed gene set: ", fs, " (n = ", length(fixed_genes_now), ")\n", sep = "")

  for (fb in freq_bins) {
    for (mk in masks) {
      cat("\n-------------------------------------------------\n")
      cat("[INFO] Running BHR for:\n")
      cat("       fixed_set = ", fs, "\n", sep = "")
      cat("       freq_bin  = ", fb, "\n", sep = "")
      cat("       mask      = ", mk, "\n", sep = "")

      sum_file <- file.path(sum_dir, fb, paste0("bhr_input_", fb, "_", mk, ".txt"))

      if (!file.exists(sum_file)) {
        cat("[WARN] Missing file, skipping: ", sum_file, "\n", sep = "")
        next
      }

      sum_df <- as.data.frame(fread(sum_file))

      fit <- BHR(
        trait1_sumstats = sum_df,
        annotations     = list(ref_dt),
        num_blocks      = 100,
        mode            = "univariate",
        fixed_genes     = fixed_genes_now
      )

      key <- paste(fs, fb, mk, sep = "__")
      bhr_results[[key]] <- fit

      cat("[DONE] Finished: ", key, "\n", sep = "")

      if (!is.null(fit$mixed_model) && !is.null(fit$mixed_model$heritabilities)) {
        print(fit$mixed_model$heritabilities)
      }
    }
  }
}

cat("\n[INFO] Total successful fits stored: ", length(bhr_results), "\n", sep = "")

############################################################
## Step 6) Save raw BHR fit objects
############################################################
ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_rdata <- file.path(out_dir, paste0("bhr_results_fixed_genes_", ts, ".RData"))

save(bhr_results, file = out_rdata)

cat("[INFO] Saved BHR results to:\n", out_rdata, "\n", sep = "")

############################################################
## Step 7) Summarize heritability estimates
############################################################
stopifnot(length(bhr_results) > 0)

############################
## Step 7.1) Helper to extract h2 table
############################
extract_h2_table <- function(fit) {
  h2mat <- fit$mixed_model$heritabilities
  if (is.null(h2mat)) return(NULL)

  dt_long <- as.data.table(as.table(h2mat))
  if (ncol(dt_long) != 3) {
    stop("Unexpected heritability table shape; ncol = ", ncol(dt_long))
  }

  setnames(dt_long, names(dt_long), c("metric", "component", "value"))
  dt_wide <- dcast(dt_long, component ~ metric, value.var = "value")

  if ("h2" %in% names(dt_wide)) {
    dt_wide[, h2 := as.numeric(h2)]
  }
  if ("h2_se" %in% names(dt_wide)) {
    dt_wide[, h2_se := as.numeric(h2_se)]
  }

  dt_wide[]
}

############################
## Step 7.2) Extract h2 across all runs
############################
h2_list <- lapply(names(bhr_results), function(key) {
  fit <- bhr_results[[key]]
  if (is.null(fit) || is.null(fit$mixed_model)) return(NULL)

  out <- extract_h2_table(fit)
  if (is.null(out)) return(NULL)

  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  fixed_set <- ifelse(length(parts) >= 1, parts[1], NA_character_)
  freq_bin  <- ifelse(length(parts) >= 2, parts[2], NA_character_)
  mask      <- ifelse(length(parts) >= 3, parts[3], NA_character_)

  out[, `:=`(
    run       = key,
    fixed_set = fixed_set,
    freq_bin  = freq_bin,
    mask      = mask
  )]

  setcolorder(
    out,
    c("run", "fixed_set", "freq_bin", "mask", "component",
      setdiff(names(out), c("run", "fixed_set", "freq_bin", "mask", "component")))
  )

  out
})

h2_dt <- rbindlist(h2_list, use.names = TRUE, fill = TRUE)

cat("[INFO] h2_dt rows: ", nrow(h2_dt), "\n", sep = "")
print(h2_dt[1:min(10, .N)])

############################
## Step 7.3) Total h2 summary
############################
h2_total <- h2_dt[component == "total"]

if (nrow(h2_total) > 0 && all(c("h2", "h2_se") %in% names(h2_total))) {
  h2_total[, z := h2 / h2_se]
  h2_total[, p_two_sided := 2 * pnorm(-abs(z))]
}

cat("[INFO] h2_total rows: ", nrow(h2_total), "\n", sep = "")
print(h2_total[order(fixed_set, freq_bin, mask)])

############################################################
## Step 8) Save observed-scale heritability summaries
############################################################
out_h2_all   <- file.path(out_dir, paste0("bhr_h2_all_components_", ts, ".txt"))
out_h2_total <- file.path(out_dir, paste0("bhr_h2_total_", ts, ".txt"))

fwrite(
  h2_dt,
  file = out_h2_all,
  sep = "\t",
  quote = FALSE,
  col.names = TRUE
)

fwrite(
  h2_total,
  file = out_h2_total,
  sep = "\t",
  quote = FALSE,
  col.names = TRUE
)

cat("[INFO] h2 tables written:\n",
    " - ", out_h2_all, "\n",
    " - ", out_h2_total, "\n", sep = "")

############################################################
## Step 9) Convert total h2 from observed scale to liability scale
############################################################
K <- 0.0147
n_cases    <- 3816
n_controls <- 339903
P <- n_cases / (n_cases + n_controls)

obs2lia_factor <- function(K, P) {
  t <- qnorm(1 - K)
  Z <- dnorm(t)
  (K * (1 - K) / (Z^2)) * (K * (1 - K) / (P * (1 - P)))
}

scale_factor <- obs2lia_factor(K, P)

cat("[INFO] Liability conversion:\n")
cat("       K = ", K, "\n", sep = "")
cat("       P = ", P, "\n", sep = "")
cat("       scale_factor = ", scale_factor, "\n", sep = "")

stopifnot(all(c("h2", "h2_se") %in% names(h2_total)))

h2_total[, `:=`(
  h2_liab    = h2 * scale_factor,
  h2_se_liab = h2_se * scale_factor
)]

h2_total[, z_liab := h2_liab / h2_se_liab]
h2_total[, p_two_sided_liab := 2 * pnorm(-abs(z_liab))]

out_h2_total_liab <- file.path(out_dir, paste0("bhr_h2_total_liability_", ts, ".txt"))

fwrite(
  h2_total,
  file = out_h2_total_liab,
  sep = "\t",
  quote = FALSE,
  col.names = TRUE
)

cat("[INFO] Saved liability-scale h2 table to:\n", out_h2_total_liab, "\n", sep = "")

############################################################
## Step 10) Split BHR fit results by fixed-gene set
############################################################
make_set_results <- function(bhr_results, fixed_set_name) {
  out <- bhr_results[
    grepl(paste0("^", fixed_set_name, "__"), names(bhr_results), perl = TRUE)
  ]
  cat("[INFO] ", fixed_set_name, ": ", length(out), " fits\n", sep = "")
  out
}

bhr_results_set1 <- make_set_results(bhr_results, "gene_set_1_UKBB128")
bhr_results_set2 <- make_set_results(bhr_results, "gene_set_2_validated45")
bhr_results_set3 <- make_set_results(bhr_results, "gene_set_3_highconf10")
bhr_results_set4 <- make_set_results(bhr_results, "gene_set_4_core4")

print(sapply(
  list(
    set1 = bhr_results_set1,
    set2 = bhr_results_set2,
    set3 = bhr_results_set3,
    set4 = bhr_results_set4
  ),
  length
))

############################################################
## Step 11) Summarize significant-gene outputs
############################################################

############################
## Step 11.1) Parse run key
############################
parse_key <- function(key) {
  parts <- strsplit(key, "__", fixed = TRUE)[[1]]
  data.table(
    run       = key,
    fixed_set = ifelse(length(parts) >= 1, parts[1], NA_character_),
    freq_bin  = ifelse(length(parts) >= 2, parts[2], NA_character_),
    mask      = ifelse(length(parts) >= 3, parts[3], NA_character_)
  )
}

############################
## Step 11.2) Extract run-level significant-gene summary
############################
extract_sig_summary <- function(fit) {
  sg <- fit$significant_genes
  if (is.null(sg)) return(NULL)

  out <- data.table(
    number_significant_genes = if (!is.null(sg$number_significant_genes)) {
      sg$number_significant_genes
    } else {
      NA_integer_
    },
    fraction_burdenh2_significant = if (!is.null(sg$fraction_burdenh2_significant)) {
      as.numeric(sg$fraction_burdenh2_significant[1, 1])
    } else {
      NA_real_
    },
    fraction_burdenh2_significant_se = if (!is.null(sg$fraction_burdenh2_significant_se)) {
      as.numeric(sg$fraction_burdenh2_significant_se[1, 1])
    } else {
      NA_real_
    }
  )

  if (!is.null(sg$significant_genes)) {
    out[, significant_genes := paste(sg$significant_genes, collapse = ",")]
  } else {
    out[, significant_genes := NA_character_]
  }

  out
}

############################
## Step 11.3) Extract gene-level significant table
############################
extract_sig_table <- function(fit) {
  sg <- fit$significant_genes
  if (is.null(sg) || is.null(sg$sig_table)) return(NULL)

  dt <- as.data.table(sg$sig_table)
  stopifnot("gene" %in% names(dt))
  dt[]
}

############################
## Step 11.4) Summarize one fixed-gene result list
############################
summarize_bhr_sig_results <- function(res_list) {
  stopifnot(is.list(res_list), length(res_list) > 0)

  summary_list <- lapply(names(res_list), function(key) {
    fit  <- res_list[[key]]
    meta <- parse_key(key)
    sig  <- extract_sig_summary(fit)
    if (is.null(sig)) return(NULL)
    cbind(meta, sig)
  })
  sig_summary_dt <- rbindlist(summary_list, use.names = TRUE, fill = TRUE)

  sigtable_list <- lapply(names(res_list), function(key) {
    fit  <- res_list[[key]]
    meta <- parse_key(key)
    dt   <- extract_sig_table(fit)
    if (is.null(dt)) return(NULL)

    dt[, `:=`(
      run       = meta$run,
      fixed_set = meta$fixed_set,
      freq_bin  = meta$freq_bin,
      mask      = meta$mask
    )]

    setcolorder(
      dt,
      c("run", "fixed_set", "freq_bin", "mask",
        setdiff(names(dt), c("run", "fixed_set", "freq_bin", "mask")))
    )

    dt
  })
  sig_table_dt <- rbindlist(sigtable_list, use.names = TRUE, fill = TRUE)

  list(
    sig_summary_dt = sig_summary_dt,
    sig_table_dt   = sig_table_dt
  )
}

############################################################
## Step 12) Run significant-gene summaries for all four sets
############################################################
sets_list <- list(
  set1_UKBB128     = bhr_results_set1,
  set2_validated45 = bhr_results_set2,
  set3_highconf10  = bhr_results_set3,
  set4_core4       = bhr_results_set4
)

for (nm in names(sets_list)) {
  res_list <- sets_list[[nm]]

  if (is.null(res_list) || length(res_list) == 0) {
    cat("[WARN] Skipping empty results list: ", nm, "\n", sep = "")
    next
  }

  cat("\n=================================================\n")
  cat("[INFO] Summarizing: ", nm, " (", length(res_list), " fits)\n", sep = "")

  out <- summarize_bhr_sig_results(res_list)
  sig_summary_dt <- out$sig_summary_dt
  sig_table_dt   <- out$sig_table_dt

  cat("[INFO] sig_summary rows: ", nrow(sig_summary_dt), "\n", sep = "")
  cat("[INFO] sig_table rows  : ", nrow(sig_table_dt), "\n", sep = "")

  f_sum <- file.path(out_dir, paste0("bhr_sig_summary_", nm, "_", ts, ".txt"))
  f_tab <- file.path(out_dir, paste0("bhr_sig_table_",   nm, "_", ts, ".txt"))

  fwrite(sig_summary_dt, f_sum, sep = "\t", quote = FALSE, col.names = TRUE)
  fwrite(sig_table_dt,   f_tab, sep = "\t", quote = FALSE, col.names = TRUE)

  cat("[INFO] Saved:\n")
  cat("  - ", f_sum, "\n", sep = "")
  cat("  - ", f_tab, "\n", sep = "")
}

cat("\n[DONE] All fixed-gene sets summarized and saved. Timestamp: ", ts, "\n", sep = "")