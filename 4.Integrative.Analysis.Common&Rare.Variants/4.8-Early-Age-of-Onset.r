############################################################
## Association of common polygenic risk and MITF rare-
## variant carrier status with age at melanoma diagnosis
##
## Purpose:
##   - Restrict to melanoma cases with available diagnosis age.
##   - Evaluate whether:
##       (1) standardized cvPRS is associated with age at diagnosis
##       (2) MITF rare-variant carrier status is associated with age
##           at diagnosis
##       (3) a prevalence-matched high-cvPRS group shows earlier
##           diagnosis
##   - Save reproducible summary tables, model outputs, and a
##     publication-ready density plot.
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
  library(ggplot2)
})

############################################################
## Step 2) Define input and output paths
############################################################
in_dt <- "/re_gecip/cancer_melanoma/Shelley_project/Project.noncoding.WGS/Updated.rvBurden/MITF.only/analysis_dx_dt.tsv"
out_dir <- "/re_gecip/cancer_melanoma/Shelley_project/Project.noncoding.WGS/Updated.rvBurden/MITF.only"

stopifnot(file.exists(in_dt))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## Step 3) Read input data
############################################################
analysis_dx_dt <- fread(in_dt)

message("[INFO] Loaded analysis_dx_dt: n = ", nrow(analysis_dx_dt), ", cols = ", ncol(analysis_dx_dt))
message(
  "[INFO] Key columns present: ",
  paste(
    intersect(
      c("platekey", "MM_bin", "cvPRS_z", "mitf_missense_carrier", "dx_age_from_event"),
      names(analysis_dx_dt)
    ),
    collapse = ", "
  )
)

str(analysis_dx_dt)

############################################################
## Step 4) Estimate MITF carrier frequency in the full cohort
############################################################
stopifnot("mitf_missense_carrier" %in% names(analysis_dx_dt))

carrier_tab <- table(analysis_dx_dt$mitf_missense_carrier, useNA = "ifany")
print(carrier_tab)

n_total <- nrow(analysis_dx_dt)
n_carrier <- sum(analysis_dx_dt$mitf_missense_carrier == 1, na.rm = TRUE)
carrier_freq <- n_carrier / n_total

message("[INFO] Total samples: ", n_total)
message("[INFO] MITF carriers: ", n_carrier)
message("[INFO] Carrier frequency: ", signif(carrier_freq, 4))

############################################################
## Step 5) Define a prevalence-matched high-cvPRS group
############################################################
stopifnot("cvPRS_z" %in% names(analysis_dx_dt))

prs_threshold <- quantile(
  analysis_dx_dt$cvPRS_z,
  probs = 1 - carrier_freq,
  na.rm = TRUE
)

message("[INFO] PRS threshold equivalent to carrier rate: ", prs_threshold)

analysis_dx_dt[, cvPRS_high_risk := as.integer(cvPRS_z >= prs_threshold)]

message("[INFO] Individuals in PRS high-risk group: ",
        sum(analysis_dx_dt$cvPRS_high_risk == 1, na.rm = TRUE))
print(table(analysis_dx_dt$cvPRS_high_risk, useNA = "ifany"))

############################################################
## Step 6) Restrict to melanoma cases with diagnosis age
############################################################
stopifnot(all(c("MM_bin", "dx_age_from_event", "mitf_missense_carrier", "cvPRS_high_risk") %in% names(analysis_dx_dt)))

dx_case_dt <- analysis_dx_dt[
  MM_bin == 1 & !is.na(dx_age_from_event)
]

message("[INFO] Cases with diagnosis age: ", nrow(dx_case_dt))

message("[INFO] MITF carriers among cases: ",
        sum(dx_case_dt$mitf_missense_carrier == 1, na.rm = TRUE))
print(table(dx_case_dt$mitf_missense_carrier, useNA = "ifany"))

message("[INFO] PRS high-risk individuals among cases: ",
        sum(dx_case_dt$cvPRS_high_risk == 1, na.rm = TRUE))
print(table(dx_case_dt$cvPRS_high_risk, useNA = "ifany"))

############################################################
## Step 7) Wilcoxon rank-sum tests for diagnosis age
############################################################

############################
## Step 7.1) MITF carriers vs non-carriers
############################
carrier_summary_dt <- dx_case_dt[
  ,
  .(
    n = .N,
    mean_dx_age = mean(dx_age_from_event),
    median_dx_age = median(dx_age_from_event),
    sd_dx_age = sd(dx_age_from_event)
  ),
  by = mitf_missense_carrier
]

print(carrier_summary_dt)

wilcox_mitf <- wilcox.test(
  dx_age_from_event ~ mitf_missense_carrier,
  data = dx_case_dt
)

print(wilcox_mitf)
message("[INFO] MITF Wilcoxon p-value: ", wilcox_mitf$p.value)

############################
## Step 7.2) High-risk PRS vs others
############################
prs_summary_dt <- dx_case_dt[
  ,
  .(
    n = .N,
    mean_dx_age = mean(dx_age_from_event),
    median_dx_age = median(dx_age_from_event),
    sd_dx_age = sd(dx_age_from_event)
  ),
  by = cvPRS_high_risk
]

print(prs_summary_dt)

wilcox_prs <- wilcox.test(
  dx_age_from_event ~ cvPRS_high_risk,
  data = dx_case_dt
)

print(wilcox_prs)
message("[INFO] PRS high-risk Wilcoxon p-value: ", wilcox_prs$p.value)

############################################################
## Step 8) Save descriptive summary tables
############################################################
fwrite(
  carrier_summary_dt,
  file = file.path(out_dir, "carrier_dx_age_summary.txt"),
  sep = "\t"
)
message("[INFO] Saved carrier summary")

fwrite(
  prs_summary_dt,
  file = file.path(out_dir, "prs_highrisk_dx_age_summary.txt"),
  sep = "\t"
)
message("[INFO] Saved PRS summary")

############################################################
## Step 9) Linear models for diagnosis age adjusting for
##         sex and ancestry PCs
############################################################

############################
## Step 9.1) Define model formulas
############################
pc_terms <- paste0("PC", 1:20)
covar_terms <- c("sex_bin", pc_terms)

formula_m1 <- paste(
  "dx_age_from_event ~ cvPRS_z +",
  paste(covar_terms, collapse = " + ")
)

formula_m2 <- paste(
  "dx_age_from_event ~ cvPRS_high_risk +",
  paste(covar_terms, collapse = " + ")
)

formula_m3 <- paste(
  "dx_age_from_event ~ mitf_missense_carrier +",
  paste(covar_terms, collapse = " + ")
)

formula_m4 <- paste(
  "dx_age_from_event ~ cvPRS_z + mitf_missense_carrier +",
  paste(covar_terms, collapse = " + ")
)

############################
## Step 9.2) Fit models
############################
fit_lm_m1 <- lm(as.formula(formula_m1), data = dx_case_dt)
fit_lm_m2 <- lm(as.formula(formula_m2), data = dx_case_dt)
fit_lm_m3 <- lm(as.formula(formula_m3), data = dx_case_dt)
fit_lm_m4 <- lm(as.formula(formula_m4), data = dx_case_dt)

############################
## Step 9.3) Extract and save tidy results
############################
tidy_m1 <- as.data.table(broom::tidy(fit_lm_m1, conf.int = TRUE))
tidy_m2 <- as.data.table(broom::tidy(fit_lm_m2, conf.int = TRUE))
tidy_m3 <- as.data.table(broom::tidy(fit_lm_m3, conf.int = TRUE))
tidy_m4 <- as.data.table(broom::tidy(fit_lm_m4, conf.int = TRUE))

fwrite(
  tidy_m1,
  file = file.path(out_dir, "lm_model1_cvPRS_dx_age_results.txt"),
  sep = "\t"
)

fwrite(
  tidy_m2,
  file = file.path(out_dir, "lm_model2_cvPRS_highrisk_dx_age_results.txt"),
  sep = "\t"
)

fwrite(
  tidy_m3,
  file = file.path(out_dir, "lm_model3_MITF_carrier_dx_age_results.txt"),
  sep = "\t"
)

fwrite(
  tidy_m4,
  file = file.path(out_dir, "lm_model4_joint_dx_age_results.txt"),
  sep = "\t"
)

message("[INFO] Linear model results saved to: ", out_dir)

############################################################
## Step 10) Density plot of diagnosis age by MITF carrier
##          status
############################################################
stopifnot(all(c("dx_age_from_event", "mitf_missense_carrier") %in% names(dx_case_dt)))

plot_dt <- copy(dx_case_dt)
plot_dt <- plot_dt[!is.na(dx_age_from_event) & !is.na(mitf_missense_carrier)]
plot_dt[, group := fifelse(mitf_missense_carrier == 1L, "Carrier", "Non-carrier")]

############################
## Step 10.1) Group counts and medians
############################
n_dt <- plot_dt[, .N, by = group]
n_carrier_case <- n_dt[group == "Carrier", N]
n_noncar_case  <- n_dt[group == "Non-carrier", N]

med_dt <- plot_dt[, .(median_dx_age = median(dx_age_from_event)), by = group]

############################
## Step 10.2) Wilcoxon p-value annotation
############################
p_wilcox <- wilcox.test(dx_age_from_event ~ group, data = plot_dt)$p.value
p_lab <- format.pval(p_wilcox, digits = 2, eps = 1e-4)

############################
## Step 10.3) Plot settings
############################
col_noncar <- "#4C78A8"
col_car    <- "#D94343"

lab_noncar <- paste0("Non-carrier (n=", n_noncar_case, ")")
lab_car    <- paste0("Carrier (n=", n_carrier_case, ")")

plot_dt[, group_lab := fifelse(group == "Carrier", lab_car, lab_noncar)]
plot_dt[, group_lab := factor(group_lab, levels = c(lab_noncar, lab_car))]

med_dt[, group_lab := fifelse(group == "Carrier", lab_car, lab_noncar)]
med_dt[, group_lab := factor(group_lab, levels = c(lab_noncar, lab_car))]

fill_vals  <- setNames(c(col_noncar, col_car), c(lab_noncar, lab_car))
color_vals <- setNames(c(col_noncar, col_car), c(lab_noncar, lab_car))

############################
## Step 10.4) Build plot
############################
p <- ggplot(plot_dt, aes(x = dx_age_from_event, fill = group_lab, color = group_lab)) +
  geom_density(alpha = 0.25, linewidth = 0.45, adjust = 1.15) +
  geom_vline(
    data = med_dt,
    aes(xintercept = median_dx_age, color = group_lab),
    linetype = "dashed",
    linewidth = 0.45,
    show.legend = FALSE
  ) +
  labs(
    title = "Age at diagnosis",
    subtitle = paste0("Carrier vs non-carrier (Wilcoxon p = ", p_lab, ")"),
    x = "Age (years)",
    y = "Density",
    fill = NULL,
    color = NULL
  ) +
  scale_fill_manual(values = fill_vals) +
  scale_color_manual(values = color_vals) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    plot.title = element_text(face = "bold", size = 7.2, margin = margin(b = 1.5)),
    plot.subtitle = element_text(size = 6.6, margin = margin(b = 2.5)),
    axis.title = element_text(size = 6.8),
    axis.text = element_text(size = 6.6),
    axis.line = element_line(linewidth = 0.35),
    axis.ticks = element_line(linewidth = 0.30),
    legend.position = "bottom",
    legend.text = element_text(size = 6.6),
    legend.key.width = unit(0.9, "lines"),
    legend.spacing.x = unit(0.5, "lines"),
    plot.margin = margin(3, 3, 3, 3)
  )

print(p)

############################################################
## Step 11) Save density plot
############################################################
ggsave(
  filename = file.path(out_dir, "Fig_dx_age_density_MITF_carrier_square.pdf"),
  plot = p,
  width = 3.1,
  height = 3.1,
  units = "in",
  device = cairo_pdf
)

message("[INFO] Saved density plot")