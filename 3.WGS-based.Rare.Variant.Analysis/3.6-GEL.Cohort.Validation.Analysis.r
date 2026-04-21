############################################################
## Firth logistic regression for rare-variant burden masks
## Supports: protein_coding, cCRE, ncRNA
## Validation for results with FDR<0.05 in the UKB discovery set
## Replication cohort: GEL
##
## Usage:
##   Rscript firth_burden_regression.R <burden_type> [chunk_id] [n_chunks]
##
## Examples:
##   Rscript firth_burden_regression.R protein_coding 1 100
##   Rscript firth_burden_regression.R cCRE 1 10
##   Rscript firth_burden_regression.R ncRNA 1 10
##
## Author: Shelley
## Date:   2025-12-10
############################################################

############################################################
## Step 1) Load packages
############################################################
suppressPackageStartupMessages({
  library(data.table)
  library(logistf)
})

############################################################
## Step 2) Parse command-line arguments
############################################################
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1L) {
  stop(
    "Usage: Rscript firth_burden_regression.R <burden_type> [chunk_id] [n_chunks]\n",
    "  burden_type: protein_coding | cCRE | ncRNA"
  )
}

burden_type <- args[1]
valid_types <- c("protein_coding", "cCRE", "ncRNA")

if (!(burden_type %in% valid_types)) {
  stop(
    "[ERROR] Invalid burden_type: ", burden_type,
    ". Must be one of: ", paste(valid_types, collapse = ", ")
  )
}

default_n_chunks <- switch(
  burden_type,
  protein_coding = 100L,
  cCRE = 10L,
  ncRNA = 10L
)

chunk_id <- if (length(args) >= 2L) as.integer(args[2]) else 1L
n_chunks <- if (length(args) >= 3L) as.integer(args[3]) else default_n_chunks

if (is.na(chunk_id) || is.na(n_chunks) ||
    chunk_id < 1L || chunk_id > n_chunks) {
  stop(
    "[ERROR] Invalid chunk_id / n_chunks: chunk_id = ", chunk_id,
    ", n_chunks = ", n_chunks
  )
}

cat("[INFO] Burden type: ", burden_type, "\n", sep = "")
cat("[INFO] Running chunk ", chunk_id, " of ", n_chunks, "\n", sep = "")

############################################################
## Step 3) Define input and output paths
############################################################
project_dir <- "/re_gecip/cancer_melanoma/Shelley_project/Project.noncoding.WGS"

burden_file <- switch(
  burden_type,
  protein_coding = file.path(project_dir, "Mask.Burden/protein.coding/protein_coding_burden_dt.txt"),
  cCRE           = file.path(project_dir, "Mask.Burden/cCRE/cCRE_burden_dt.txt"),
  ncRNA          = file.path(project_dir, "Mask.Burden/ncRNA/ncRNA_burden_dt.txt")
)

pheno_file <- file.path(
  project_dir,
  "Validation_cohort/cleaned.phenotypes/Validation_pheno_cleaned_EUR_unrelated.txt"
)

base_out_dir <- switch(
  burden_type,
  protein_coding = file.path(project_dir, "Firth.Regression.Results/protein.coding"),
  cCRE           = file.path(project_dir, "Firth.Regression.Results/cCRE"),
  ncRNA          = file.path(project_dir, "Firth.Regression.Results/ncRNA")
)

out_dir <- file.path(base_out_dir, sprintf("chunk%03d", chunk_id))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

cat("[INFO] Burden file: ", burden_file, "\n", sep = "")
cat("[INFO] Phenotype file: ", pheno_file, "\n", sep = "")
cat("[INFO] Output directory: ", out_dir, "\n", sep = "")

############################################################
## Step 4) Read burden and phenotype data
############################################################
cat("[INFO] Reading burden data...\n")
burden_dt <- fread(burden_file)

cat("[INFO] Burden data loaded: ",
    nrow(burden_dt), " rows x ", ncol(burden_dt), " columns\n", sep = "")
str(burden_dt, max.level = 1)

cat("[INFO] Reading phenotype data...\n")
pheno_dt <- fread(pheno_file)

cat("[INFO] Phenotype data loaded: ",
    nrow(pheno_dt), " rows x ", ncol(pheno_dt), " columns\n", sep = "")
str(pheno_dt, max.level = 1)

############################################################
## Step 5) Harmonize IDs and merge burden with phenotype
############################################################
burden_dt[, IID := as.character(IID)]
pheno_dt[, Platekey := as.character(Platekey)]

n_burden_ids <- uniqueN(burden_dt$IID)
n_pheno_ids  <- uniqueN(pheno_dt$Platekey)
n_overlap    <- length(intersect(burden_dt$IID, pheno_dt$Platekey))

cat("[INFO] Unique IDs in burden_dt: ", n_burden_ids, "\n", sep = "")
cat("[INFO] Unique IDs in pheno_dt : ", n_pheno_ids, "\n", sep = "")
cat("[INFO] Overlapping IDs        : ", n_overlap, "\n", sep = "")

merged_dt <- merge(
  x = burden_dt,
  y = pheno_dt,
  by.x = "IID",
  by.y = "Platekey",
  all = FALSE
)

cat("[INFO] After inner join, merged_dt has ",
    nrow(merged_dt), " rows x ", ncol(merged_dt), " columns\n", sep = "")
cat("[INFO] Missing case_control after merge: ",
    sum(is.na(merged_dt$case_control)), "\n", sep = "")

str(merged_dt, max.level = 1)

############################################################
## Step 6) Remove masks with MAC < 5
############################################################
burden_cols <- setdiff(names(burden_dt), c("FID", "IID"))
cat("[INFO] Number of burden mask columns: ", length(burden_cols), "\n", sep = "")

stopifnot(all(burden_cols %in% names(merged_dt)))

mac_vec <- colSums(merged_dt[, ..burden_cols], na.rm = TRUE)

cat("[INFO] MAC summary across burden columns:\n")
print(summary(mac_vec))

keep_burden_cols <- burden_cols[mac_vec >= 5]
drop_burden_cols <- setdiff(burden_cols, keep_burden_cols)

cat("[INFO] Masks with MAC >= 5 kept   : ", length(keep_burden_cols), "\n", sep = "")
cat("[INFO] Masks with MAC  < 5 dropped: ", length(drop_burden_cols), "\n", sep = "")

if (length(keep_burden_cols) == 0L) {
  cat("[WARN] No masks with MAC >= 5; exiting.\n")
  quit(save = "no", status = 0)
}

keep_all_cols <- c(
  setdiff(names(merged_dt), burden_cols),
  keep_burden_cols
)

merged_filt_dt <- merged_dt[, ..keep_all_cols]

cat("[INFO] After MAC filtering, merged_filt_dt has ",
    nrow(merged_filt_dt), " rows x ", ncol(merged_filt_dt), " columns\n", sep = "")

str(merged_filt_dt, max.level = 1)

############################################################
## Step 7) Build phenotype and covariates
############################################################
analysis_dt <- merged_filt_dt[case_control %in% c("case", "control")]

cat("[INFO] Samples with valid case_control: ", nrow(analysis_dt), "\n", sep = "")

analysis_dt[, pheno_case := fifelse(case_control == "case", 1L, 0L)]
analysis_dt[, sex_bin := fifelse(`Participant Phenotypic Sex` == "Male", 1L, 0L)]

analysis_dt[, age := age_at_recruitment]
analysis_dt[, age_sq := age^2]

pc_cols <- paste0("Pc", 1:20)
stopifnot(all(pc_cols %in% names(analysis_dt)))

cat("[INFO] PC columns used: ", paste(pc_cols, collapse = ", "), "\n", sep = "")

str(
  analysis_dt[, c("IID", "pheno_case", "sex_bin", "age", "age_sq", pc_cols), with = FALSE],
  max.level = 1
)

############################################################
## Step 8) Assign masks to the current chunk
############################################################
n_masks <- length(keep_burden_cols)
cat("[INFO] Total masks after MAC filter: ", n_masks, "\n", sep = "")

mask_indices_list <- split(
  seq_len(n_masks),
  cut(
    seq_len(n_masks),
    breaks = n_chunks,
    labels = FALSE
  )
)

current_idx <- mask_indices_list[[chunk_id]]

if (length(current_idx) == 0L) {
  cat("[WARN] No masks assigned to chunk ", chunk_id, "; exiting.\n", sep = "")
  quit(save = "no", status = 0)
}

all_masks <- keep_burden_cols[current_idx]
cat("[INFO] Masks in this chunk: ", length(all_masks), "\n", sep = "")

############################################################
## Step 9) Run Firth logistic regression for this chunk
############################################################
firth_summary_list <- vector("list", length(all_masks))
names(firth_summary_list) <- all_masks

failed_masks  <- character(0)
failed_reason <- character(0)

fit_ctrl <- logistf.control(maxit = 2000)
pl_ctrl  <- logistpl.control(maxit = 2000)

for (i in seq_along(all_masks)) {
  m <- all_masks[i]

  cat("[INFO] (", i, "/", length(all_masks), ") Running Firth for mask: ",
      m, "\n", sep = "")

  ## Wrap mask names in backticks to safely handle non-syntactic names
  mask_term <- sprintf("`%s`", m)

  form_str <- paste0(
    "pheno_case ~ ", mask_term,
    " + age + age_sq + sex_bin",
    " + age:sex_bin + age_sq:sex_bin",
    " + ", paste(pc_cols, collapse = " + ")
  )
  form <- as.formula(form_str)

  res <- tryCatch(
    {
      fit_m <- logistf(
        formula = form,
        data = analysis_dt,
        control = fit_ctrl,
        plcontrol = pl_ctrl
      )

      sum_m <- summary(fit_m)

      out_txt <- file.path(out_dir, paste0("Firth_", m, "_summary.txt"))
      capture.output(summary(fit_m), file = out_txt)

      sum_m
    },
    error = function(e) {
      cat("[WARN] Firth failed for mask ", m, ":\n     ",
          conditionMessage(e), "\n", sep = "")
      failed_masks  <<- c(failed_masks, m)
      failed_reason <<- c(failed_reason, conditionMessage(e))
      NULL
    }
  )

  firth_summary_list[[m]] <- res
}

cat("[INFO] Firth loop finished for ", burden_type, " chunk ", chunk_id, "\n", sep = "")
cat("[INFO] Successful masks in this chunk: ",
    sum(!vapply(firth_summary_list, is.null, logical(1))), "\n", sep = "")
cat("[INFO] Failed masks in this chunk    : ",
    length(failed_masks), "\n", sep = "")

############################################################
## Step 10) Save outputs for this chunk
############################################################
rds_prefix <- switch(
  burden_type,
  protein_coding = "protein_coding_firth",
  cCRE           = "cCRE_firth",
  ncRNA          = "ncRNA_firth"
)

saveRDS(
  firth_summary_list,
  file = file.path(
    out_dir,
    sprintf("%s_chunk%03d_summaries.rds", rds_prefix, chunk_id)
  )
)

if (length(failed_masks) > 0) {
  failed_dt <- data.table(mask = failed_masks, reason = failed_reason)
  fwrite(
    failed_dt,
    file = file.path(
      out_dir,
      sprintf("%s_chunk%03d_failed_masks.txt", rds_prefix, chunk_id)
    ),
    sep = "\t"
  )
}

cat("[INFO] Done. ", burden_type, " chunk ", chunk_id,
    " finished at ", format(Sys.time()), "\n", sep = "")