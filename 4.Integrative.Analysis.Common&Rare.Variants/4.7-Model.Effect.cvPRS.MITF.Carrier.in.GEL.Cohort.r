############################################################
## Estimate odds ratios for cvPRS and MITF missense
## carrier status in the GEL cohort
##
## Purpose:
##   - Fit logistic regression models to estimate odds ratios
##     for:
##       (1) common-variant PRS (cvPRS)
##       (2) MITF missense-variant carrier status
##     on melanoma case-control risk in the GEL cohort.
##   - Quantify model improvement using Nagelkerke R2.
##
## Author: Shelley (Shiyu Zhang)
## Date:   2026-03-05
############################################################

############################################################
## Step 1) Load packages
############################################################
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(broom)
  library(rcompanion)
})

############################################################
## Step 2) Define input and output paths
############################################################
carrier_raw <- "/re_gecip/cancer_melanoma/Shelley_project/Project.noncoding.WGS/GT.dosage/MITF/chr3_69964940_G_A.raw"
pheno_tsv   <- "/re_gecip/cancer_melanoma/Shelley_project/Project.noncoding.WGS/PRS/performance/general.cohort/prs_eur_unrelated_qc_dt.tsv"
dx_age_file <- "/re_gecip/cancer_melanoma/Shelley_project/Project.noncoding.WGS/PRS/rvBurden/Trial.1/results/melanoma.cases.age.at.diagnosis.info.txt"

out_dir <- "/re_gecip/cancer_melanoma/Shelley_project/Project.noncoding.WGS/Updated.rvBurden/MITF.only"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(carrier_raw))
stopifnot(file.exists(pheno_tsv))
stopifnot(file.exists(dx_age_file))

############################################################
## Step 3) Read MITF carrier-status data
############################################################
carrier_dt <- fread(carrier_raw)

## PLINK2 .raw dosage column is expected to be the last column
dosage_col <- names(carrier_dt)[ncol(carrier_dt)]
message("[INFO] Dosage column detected: ", dosage_col)

carrier_dt <- carrier_dt[, .(
  FID,
  IID,
  mitf_var_dosage = as.numeric(get(dosage_col))
)]

## Carrier definition: at least one ALT allele
carrier_dt[, mitf_missense_carrier := fifelse(
  is.na(mitf_var_dosage), NA_integer_,
  fifelse(mitf_var_dosage >= 1, 1L, 0L)
)]

############################################################
## Step 4) Read phenotype / covariate data
############################################################
pheno_dt <- fread(pheno_tsv)

message("[INFO] Loaded carrier_dt: n = ", nrow(carrier_dt), ", cols = ", ncol(carrier_dt))
message("[INFO] Loaded pheno_dt:   n = ", nrow(pheno_dt),   ", cols = ", ncol(pheno_dt))

############################################################
## Step 5) Merge phenotype table with MITF carrier status
##   Join key: pheno_dt$platekey <-> carrier_dt$IID
############################################################
stopifnot("platekey" %in% names(pheno_dt))
stopifnot("IID" %in% names(carrier_dt))

message("[INFO] Unique platekey in pheno_dt: ", uniqueN(pheno_dt$platekey), " / n = ", nrow(pheno_dt))
message("[INFO] Unique IID in carrier_dt:    ", uniqueN(carrier_dt$IID),     " / n = ", nrow(carrier_dt))

merged_dt <- merge(
  pheno_dt,
  carrier_dt[, .(IID, mitf_var_dosage, mitf_missense_carrier)],
  by.x = "platekey",
  by.y = "IID",
  all.x = TRUE,
  sort = FALSE
)

message("[INFO] merged_dt: n = ", nrow(merged_dt), ", cols = ", ncol(merged_dt))
message("[INFO] Non-missing carrier calls: ",
        sum(!is.na(merged_dt$mitf_missense_carrier)), " / ", nrow(merged_dt))
message("[INFO] Carrier count: ",
        sum(merged_dt$mitf_missense_carrier == 1L, na.rm = TRUE))

############################################################
## Step 6) QC checks before model fitting
############################################################

############################
## Step 6.1) PRS distribution
############################
message("[INFO] Summary of PRSCS_PRS_z:")
print(summary(merged_dt$PRSCS_PRS_z))

message("[INFO] Quantiles of PRSCS_PRS_z:")
print(quantile(
  merged_dt$PRSCS_PRS_z,
  probs = c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1),
  na.rm = TRUE
))

message("[INFO] SD of PRSCS_PRS_z:")
print(sd(merged_dt$PRSCS_PRS_z, na.rm = TRUE))

############################
## Step 6.2) MITF carrier status by case/control
############################
tab_mitf <- table(
  merged_dt$MM_label,
  merged_dt$mitf_missense_carrier,
  useNA = "ifany"
)

message("[INFO] MITF missense carrier status by MM_label:")
print(tab_mitf)

message("[INFO] Carrier proportion within MM groups:")
print(prop.table(tab_mitf, margin = 1))

############################################################
## Step 7) Build analysis dataset
##   - remove samples with missing MITF carrier information
##   - recalculate standardized cvPRS from PRSCS_PRS_SCORE1_SUM
############################################################
analysis_dt <- merged_dt[!is.na(mitf_missense_carrier)]

message("[INFO] N before filtering: ", nrow(merged_dt))
message("[INFO] N after removing NA carrier status: ", nrow(analysis_dt))

analysis_dt[, cvPRS_z := as.numeric(scale(PRSCS_PRS_SCORE1_SUM))]
analysis_dt[, mitf_carrier_f := factor(mitf_missense_carrier, levels = c(0, 1))]

message("[INFO] Summary of recalculated cvPRS_z:")
print(summary(analysis_dt$cvPRS_z))
message("[INFO] Mean cvPRS_z:")
print(mean(analysis_dt$cvPRS_z, na.rm = TRUE))
message("[INFO] SD cvPRS_z:")
print(sd(analysis_dt$cvPRS_z, na.rm = TRUE))

message("[INFO] Case/control by MITF carrier:")
print(table(
  analysis_dt$MM_label,
  analysis_dt$mitf_missense_carrier,
  useNA = "ifany"
))

############################################################
## Step 8) Define covariates for logistic models
############################################################
pc_terms <- paste0("PC", 1:20)
covar_terms <- c(
  pc_terms,
  "age_at_consent", "age2", "sex_bin", "age_sex", "age2_sex"
)

############################################################
## Step 9) Fit logistic regression models for case-control risk
############################################################

############################
## Step 9.1) Helper
############################
fit_logistic_model <- function(formula_str, data) {
  glm(
    formula = as.formula(formula_str),
    data = data,
    family = binomial()
  )
}

############################
## Step 9.2) Model formulas
############################
formula_m1 <- paste(
  "MM_bin ~ cvPRS_z +",
  paste(covar_terms, collapse = " + ")
)

formula_m2 <- paste(
  "MM_bin ~ mitf_carrier_f +",
  paste(covar_terms, collapse = " + ")
)

formula_m3 <- paste(
  "MM_bin ~ cvPRS_z + mitf_carrier_f +",
  paste(covar_terms, collapse = " + ")
)

############################
## Step 9.3) Fit models
############################
fit_logit_m1 <- fit_logistic_model(formula_m1, analysis_dt)
fit_logit_m2 <- fit_logistic_model(formula_m2, analysis_dt)
fit_logit_m3 <- fit_logistic_model(formula_m3, analysis_dt)

############################################################
## Step 10) Save logistic regression results on OR scale
############################################################
add_or_cols <- function(dt, model_label) {
  dt <- as.data.table(dt)
  dt[, model := model_label]
  dt[, OR := exp(estimate)]
  dt[, lower95 := exp(conf.low)]
  dt[, upper95 := exp(conf.high)]
  dt
}

tidy_m1 <- add_or_cols(
  broom::tidy(fit_logit_m1, conf.int = TRUE),
  "Model1_cvPRS"
)
tidy_m2 <- add_or_cols(
  broom::tidy(fit_logit_m2, conf.int = TRUE),
  "Model2_MITF_carrier"
)
tidy_m3 <- add_or_cols(
  broom::tidy(fit_logit_m3, conf.int = TRUE),
  "Model3_joint"
)

tidy_all <- rbindlist(list(tidy_m1, tidy_m2, tidy_m3), use.names = TRUE, fill = TRUE)

fwrite(tidy_m1, file.path(out_dir, "model1_cvPRS_OR_results.tsv"), sep = "\t")
fwrite(tidy_m2, file.path(out_dir, "model2_MITF_carrier_OR_results.tsv"), sep = "\t")
fwrite(tidy_m3, file.path(out_dir, "model3_joint_OR_results.tsv"), sep = "\t")
fwrite(tidy_all, file.path(out_dir, "model1-3_combined_OR_results.tsv"), sep = "\t")

message("[INFO] Logistic regression results saved to: ", out_dir)

############################################################
## Step 11) Save analysis dataset
############################################################
fwrite(
  analysis_dt,
  file = file.path(out_dir, "analysis_dt.tsv"),
  sep = "\t"
)

message("[INFO] Saved analysis_dt.tsv")

############################################################
## Step 12) Evaluate model improvement using Nagelkerke R2
############################################################

############################
## Step 12.1) Fit nested models
############################
fit_r2_m0 <- glm(
  formula = as.formula(paste("MM_bin ~", paste(covar_terms, collapse = " + "))),
  data = analysis_dt,
  family = binomial()
)

fit_r2_m1 <- fit_logit_m1
fit_r2_m2 <- fit_logit_m2
fit_r2_m3 <- fit_logit_m3

############################
## Step 12.2) Compute Nagelkerke R2 for nested comparisons
############################
nk_m1_vs_m0 <- nagelkerke(fit_r2_m1, null = fit_r2_m0, restrictNobs = TRUE)
nk_m2_vs_m0 <- nagelkerke(fit_r2_m2, null = fit_r2_m0, restrictNobs = TRUE)
nk_m3_vs_m1 <- nagelkerke(fit_r2_m3, null = fit_r2_m1, restrictNobs = TRUE)

extract_nk <- function(nk_obj, comparison_name) {
  r2 <- nk_obj$Pseudo.R.squared.for.model.vs.null[
    "Nagelkerke (Cragg and Uhler)", "Pseudo.R.squared"
  ]

  chisq <- nk_obj$Likelihood.ratio.test[1, "Chisq"]
  pval  <- nk_obj$Likelihood.ratio.test[1, "p.value"]

  data.table(
    comparison = comparison_name,
    Nagelkerke_R2 = r2,
    LR_chisq = chisq,
    LR_pvalue = pval
  )
}

r2_summary_dt <- rbindlist(list(
  extract_nk(nk_m1_vs_m0, "Model1_vs_Model0 (cvPRS effect)"),
  extract_nk(nk_m2_vs_m0, "Model2_vs_Model0 (MITF effect)"),
  extract_nk(nk_m3_vs_m1, "Model3_vs_Model1 (MITF beyond PRS)")
))

print(r2_summary_dt)

############################
## Step 12.3) Save nested Nagelkerke outputs
############################
fwrite(
  r2_summary_dt,
  file.path(out_dir, "Nagelkerke_R2_summary.tsv"),
  sep = "\t"
)

nk_outfile <- file.path(out_dir, "Nagelkerke_outputs_full.txt")
con <- file(nk_outfile, open = "wt")
sink(con)
sink(con, type = "message")

cat("############################################################\n")
cat("# Nagelkerke outputs (rcompanion)\n")
cat("# Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("############################################################\n\n")

cat("===== nk_m1_vs_m0: Model1 vs Model0 (cvPRS effect) =====\n")
print(nk_m1_vs_m0)
cat("\n\n")

cat("===== nk_m2_vs_m0: Model2 vs Model0 (MITF effect) =====\n")
print(nk_m2_vs_m0)
cat("\n\n")

cat("===== nk_m3_vs_m1: Model3 vs Model1 (MITF beyond PRS) =====\n")
print(nk_m3_vs_m1)
cat("\n\n")

sink(type = "message")
sink()
close(con)

message("[INFO] Saved Nagelkerke summary and full outputs")

############################################################
## Step 13) Absolute Nagelkerke R2 for each model
############################################################
fit_null0 <- glm(MM_bin ~ 1, data = analysis_dt, family = binomial())

get_abs_nk_r2 <- function(fit, null_fit, model_name) {
  nk <- nagelkerke(fit, null = null_fit, restrictNobs = TRUE)
  r2 <- nk$Pseudo.R.squared.for.model.vs.null[
    "Nagelkerke (Cragg and Uhler)", "Pseudo.R.squared"
  ]

  data.table(
    model = model_name,
    Nagelkerke_R2 = as.numeric(r2)
  )
}

r2_models_dt <- rbindlist(list(
  get_abs_nk_r2(fit_r2_m0, fit_null0, "Model0_covariates"),
  get_abs_nk_r2(fit_r2_m1, fit_null0, "Model1_covariates+PRS"),
  get_abs_nk_r2(fit_r2_m2, fit_null0, "Model2_covariates+MITF"),
  get_abs_nk_r2(fit_r2_m3, fit_null0, "Model3_joint")
))

print(r2_models_dt)

fwrite(
  r2_models_dt,
  file.path(out_dir, "Nagelkerke_R2_models_absolute.tsv"),
  sep = "\t"
)

message("[INFO] Saved absolute Nagelkerke R2 table")