############################################################
## Summarize and visualize validation results for
## protein-coding rare-variant associations
## across three independent cohorts
##
## Discovery cohort:
##   - UK Biobank (UKBB), WGS
##
## Validation cohorts:
##   - All of Us (AoU), WGS
##   - Mass General Brigham Biobank (MGB), WES
##   - Genomics England (GEL), WGS
##
## Author: Shelley
## Date:   2026-02-24
############################################################

############################################################
## Step 1) Load packages
############################################################
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
  library(ComplexUpset)
  library(scales)
  library(grid)
  library(yorkregression)
  library(ggrepel)
})

############################################################
## Step 2) Read data
############################################################

## discovery
discovery_file <- "/Users/shelleyz/Projects/Melanoma-WGS/Updated_RareVariantAnalysis/FDR005_Associations/protein_coding.FDR0.05.associations.txt"
stopifnot(file.exists(discovery_file))

disc_dt <- fread(discovery_file)
cat("[INFO] Discovery loaded: ", discovery_file, "\n", sep = "")
cat("[INFO] disc_dt: ", nrow(disc_dt), " rows x ", ncol(disc_dt), " cols\n", sep = "")
print(head(disc_dt, 3))
cat("[INFO] Discovery columns:\n")
print(names(disc_dt))

## validation
aou_file <- "/Users/shelleyz/Projects/Melanoma-WGS/Validation/20251211_proteincoding/AllofUS.tsv"
mgb_file <- "/Users/shelleyz/Projects/Melanoma-WGS/Validation/20251211_proteincoding/MGB_53k.tsv"
gel_file <- "/Users/shelleyz/Projects/Melanoma-WGS/Validation/100kGP-GEL/Results/Firth_protein_coding_all_masks_results.tsv"

stopifnot(file.exists(aou_file))
stopifnot(file.exists(mgb_file))
stopifnot(file.exists(gel_file))

aou_dt <- fread(aou_file)
mgb_dt <- fread(mgb_file)
gel_dt <- fread(gel_file)

cat("[INFO] AoU loaded: ", aou_file, "\n", sep = "")
cat("[INFO] aou_dt: ", nrow(aou_dt), " rows x ", ncol(aou_dt), " cols\n", sep = "")
print(head(aou_dt, 3))
cat("[INFO] AoU columns:\n")
print(names(aou_dt))

cat("[INFO] MGB loaded: ", mgb_file, "\n", sep = "")
cat("[INFO] mgb_dt: ", nrow(mgb_dt), " rows x ", ncol(mgb_dt), " cols\n", sep = "")
print(head(mgb_dt, 3))
cat("[INFO] MGB columns:\n")
print(names(mgb_dt))

cat("[INFO] GEL loaded: ", gel_file, "\n", sep = "")
cat("[INFO] gel_dt: ", nrow(gel_dt), " rows x ", ncol(gel_dt), " cols\n", sep = "")
print(head(gel_dt, 3))
cat("[INFO] GEL columns:\n")
print(names(gel_dt))

############################################################
## Step 3) Harmonize discovery and validation tables
############################################################

############################
## Step 3.1) Discovery (UKBB)
############################
req_disc <- c(
  "ensembl_gene_id", "gene_symbol", "Group", "max_MAF",
  "Pvalue", "Pvalue_Burden", "Pvalue_SKAT",
  "BETA_Burden", "SE_Burden", "FDR"
)
missing_disc <- setdiff(req_disc, names(disc_dt))
stopifnot(length(missing_disc) == 0)

disc_h <- copy(disc_dt)[
  ,
  .(
    ENSG         = as.character(ensembl_gene_id),
    gene_symbol  = as.character(gene_symbol),
    Group        = as.character(Group),
    max_MAF      = as.numeric(max_MAF),
    Pvalue       = as.numeric(Pvalue),
    Pvalue_Burden = as.numeric(Pvalue_Burden),
    Pvalue_SKAT   = as.numeric(Pvalue_SKAT),
    BETA_Burden  = as.numeric(BETA_Burden),
    SE_Burden    = as.numeric(SE_Burden),
    FDR          = as.numeric(FDR)
  )
]

############################
## Step 3.2) AoU
############################
disc_genes <- unique(disc_h$gene_symbol)

aou_sub <- aou_dt[
  Gene %in% disc_genes | `Search Gene` %in% disc_genes
]

cat("[INFO] AoU subset: ", nrow(aou_sub), " rows x ", ncol(aou_sub), " cols\n", sep = "")

req_aou <- c("Gene", "Best Mask", "P-Value: Best Mask", "Cases", "Controls", "Beta")
missing_aou <- setdiff(req_aou, names(aou_sub))
stopifnot(length(missing_aou) == 0)

aou_h <- copy(aou_sub)[
  ,
  .(
    gene_symbol          = as.character(Gene),
    AoU_Best_Mask        = as.character(`Best Mask`),
    AoU_Pvalue_Best_Mask = as.numeric(`P-Value: Best Mask`),
    Cases_AoU            = as.integer(Cases),
    Controls_AoU         = as.integer(Controls),
    BETA_AoU             = as.numeric(Beta)
  )
]

############################
## Step 3.3) MGB
############################
mgb_sub <- mgb_dt[
  Gene %in% disc_genes | `Search Gene` %in% disc_genes
]

cat("[INFO] MGB subset: ", nrow(mgb_sub), " rows x ", ncol(mgb_sub), " cols\n", sep = "")

req_mgb <- c("Gene", "Best Mask", "P-Value: Best Mask", "Cases", "Controls", "Beta")
missing_mgb <- setdiff(req_mgb, names(mgb_sub))
stopifnot(length(missing_mgb) == 0)

mgb_h <- copy(mgb_sub)[
  ,
  .(
    gene_symbol          = as.character(Gene),
    MGB_Best_Mask        = as.character(`Best Mask`),
    MGB_Pvalue_Best_Mask = as.numeric(`P-Value: Best Mask`),
    Cases_MGB            = as.integer(Cases),
    Controls_MGB         = as.integer(Controls),
    BETA_MGB             = as.numeric(Beta)
  )
]

############################
## Step 3.4) GEL
############################
req_gel <- c("mask", "p", "coef", "se")
missing_gel <- setdiff(req_gel, names(gel_dt))
stopifnot(length(missing_gel) == 0)

gel_h <- copy(gel_dt)
gel_h[, c("ENSG", "GEL_Best_Mask") := tstrsplit(mask, "__", fixed = TRUE)]
gel_h[, GEL_Pvalue := as.numeric(p)]
gel_h[, BETA_GEL   := as.numeric(coef)]
gel_h[, SE_GEL     := as.numeric(se)]

gel_best <- gel_h[
  is.finite(GEL_Pvalue) & !is.na(GEL_Pvalue) & GEL_Pvalue > 0 & GEL_Pvalue <= 1
][
  order(ENSG, GEL_Pvalue)
][
  , .SD[1], by = ENSG
]

gel_h2 <- copy(gel_best)[
  ,
  .(
    ENSG                 = as.character(ENSG),
    GEL_Best_Mask        = as.character(GEL_Best_Mask),
    GEL_Pvalue_Best_Mask = as.numeric(GEL_Pvalue),
    GEL_BETA             = as.numeric(BETA_GEL),
    GEL_SE               = as.numeric(SE_GEL)
  )
]

ensg2sym <- unique(disc_h[, .(ENSG, gene_symbol)])

gel_h3 <- merge(
  gel_h2,
  ensg2sym,
  by = "ENSG",
  all.x = TRUE
)

cat("[INFO] Missing gene_symbol after GEL merge: ", sum(is.na(gel_h3$gene_symbol)), "\n", sep = "")

############################################################
## Step 4) Define validation criteria
############################################################
n_targets <- 128L
p_nominal <- 0.05
p_strong  <- 0.05 / n_targets

disc_lead <- copy(disc_h)[
  is.finite(Pvalue) & !is.na(Pvalue)
][
  order(gene_symbol, Pvalue)
][
  , .SD[1], by = gene_symbol
][
  , .(
    gene_symbol,
    ENSG,
    DISC_Pvalue = as.numeric(Pvalue),
    DISC_BETA   = as.numeric(BETA_Burden)
  )
]

add_validation_flags <- function(dt, cohort_prefix, p_col, beta_col) {
  stopifnot("gene_symbol" %in% names(dt))
  stopifnot(all(c(p_col, beta_col) %in% names(dt)))

  x <- merge(dt, disc_lead, by = "gene_symbol", all.x = TRUE)

  x[, beta_align := {
    b1 <- get(beta_col)
    b0 <- DISC_BETA
    is.finite(b1) & !is.na(b1) &
      is.finite(b0) & !is.na(b0) &
      (sign(b1) == sign(b0))
  }]

  x[, nominal_validation := beta_align & is.finite(get(p_col)) & (get(p_col) < p_nominal)]
  x[, strong_validation  := beta_align & is.finite(get(p_col)) & (get(p_col) < p_strong)]

  setnames(
    x,
    old = c("beta_align", "nominal_validation", "strong_validation"),
    new = c(
      paste0(cohort_prefix, "_beta_align"),
      paste0(cohort_prefix, "_nominal_validation"),
      paste0(cohort_prefix, "_strong_validation")
    )
  )

  x
}

aou_v <- add_validation_flags(
  dt = aou_h,
  cohort_prefix = "AoU",
  p_col = "AoU_Pvalue_Best_Mask",
  beta_col = "BETA_AoU"
)

mgb_v <- add_validation_flags(
  dt = mgb_h,
  cohort_prefix = "MGB",
  p_col = "MGB_Pvalue_Best_Mask",
  beta_col = "BETA_MGB"
)

gel_v <- add_validation_flags(
  dt = gel_h3,
  cohort_prefix = "GEL",
  p_col = "GEL_Pvalue_Best_Mask",
  beta_col = "GEL_BETA"
)

cat("[INFO] Strong threshold: ", signif(p_strong, 3), " (0.05/", n_targets, ")\n", sep = "")

############################################################
## Step 5) Save validation tables
############################################################
out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/Validation/Updated_Validation_20260224"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(aou_v, file.path(out_dir, "AoU.validation_table.tsv"), sep = "\t")
fwrite(mgb_v, file.path(out_dir, "MGB.validation_table.tsv"), sep = "\t")
fwrite(gel_v, file.path(out_dir, "GEL.validation_table.tsv"), sep = "\t")

cat("[INFO] Saved validation tables to:\n  ", out_dir, "\n", sep = "")

############################################################
## Step 6) UpSet plot for nominal validation
############################################################
up_nominal <- merge(
  aou_v[, .(gene_symbol, AoU = as.logical(AoU_nominal_validation))],
  mgb_v[, .(gene_symbol, MGB = as.logical(MGB_nominal_validation))],
  by = "gene_symbol", all = TRUE
)
up_nominal <- merge(
  up_nominal,
  gel_v[, .(gene_symbol, GEL = as.logical(GEL_nominal_validation))],
  by = "gene_symbol", all = TRUE
)

up_nominal[is.na(AoU), AoU := FALSE]
up_nominal[is.na(MGB), MGB := FALSE]
up_nominal[is.na(GEL), GEL := FALSE]

up_nominal <- up_nominal[AoU | MGB | GEL]

sets <- c("AoU", "MGB", "GEL")
p_up_nominal <- upset(
  as.data.frame(up_nominal),
  intersect = sets,
  name = NULL,
  sort_intersections_by = "cardinality",
  base_annotations = list(
    "Intersection size" = intersection_size(
      counts = TRUE,
      text = list(size = 2.0),
      width = 0.55
    )
  ),
  matrix = intersection_matrix(
    geom = geom_point(size = 1.3)
  ),
  set_sizes = upset_set_size(),
  width_ratio = 0.18,
  min_size = 1
) +
  labs(title = "Nominal validation (p<0.05) with effect alignment") +
  theme_classic(base_size = 6) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    axis.text.x  = element_text(size = 6, color = "black"),
    axis.text.y  = element_text(size = 6, color = "black"),
    plot.title = element_text(size = 5, hjust = 0),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

out_pdf <- file.path(out_dir, "Fig_upset_validation.nominal.smallpanel.pdf")
ggsave(
  filename = out_pdf,
  plot = p_up_nominal,
  width = 80, height = 70, units = "mm",
  useDingbats = FALSE
)

cat("[INFO] Saved:\n  ", out_pdf, "\n", sep = "")

############################################################
## Step 7) Print genes with nominal validation in >=2 cohorts
############################################################
nom2_dt <- merge(
  aou_v[, .(gene_symbol, AoU_nominal = as.logical(AoU_nominal_validation))],
  mgb_v[, .(gene_symbol, MGB_nominal = as.logical(MGB_nominal_validation))],
  by = "gene_symbol", all = TRUE
)
nom2_dt <- merge(
  nom2_dt,
  gel_v[, .(gene_symbol, GEL_nominal = as.logical(GEL_nominal_validation))],
  by = "gene_symbol", all = TRUE
)

nom2_dt[is.na(AoU_nominal), AoU_nominal := FALSE]
nom2_dt[is.na(MGB_nominal), MGB_nominal := FALSE]
nom2_dt[is.na(GEL_nominal), GEL_nominal := FALSE]

nom2_dt[, n_valid_cohorts := as.integer(AoU_nominal) + as.integer(MGB_nominal) + as.integer(GEL_nominal)]
nom2_dt2 <- nom2_dt[n_valid_cohorts >= 2][order(-n_valid_cohorts, gene_symbol)]

cat("[INFO] Genes with nominal validation in >=2 cohorts: ", nrow(nom2_dt2), "\n", sep = "")
print(nom2_dt2)

############################################################
## Step 8) Directional concordance by cohort
############################################################
dir_concordance_test <- function(dt, cohort_name, beta_col, disc_beta_col = "DISC_BETA") {
  stopifnot(all(c(beta_col, disc_beta_col) %in% names(dt)))

  x <- as.data.table(copy(dt))
  x <- x[
    is.finite(get(beta_col)) & !is.na(get(beta_col)) &
      is.finite(get(disc_beta_col)) & !is.na(get(disc_beta_col))
  ]
  x <- x[get(disc_beta_col) != 0 & get(beta_col) != 0]
  x[, dir_match := (sign(get(beta_col)) == sign(get(disc_beta_col)))]

  N <- nrow(x)
  K <- x[dir_match == TRUE, .N]

  bt <- binom.test(K, N, p = 0.5, alternative = "greater")

  out <- list(
    cohort = cohort_name,
    N = N,
    K = K,
    concordance_rate = K / N,
    p_value = unname(bt$p.value),
    conf_low = unname(bt$conf.int[1]),
    conf_high = unname(bt$conf.int[2])
  )

  invisible(list(summary = out, detail = x))
}

aou_dir <- dir_concordance_test(aou_v, cohort_name = "AoU", beta_col = "BETA_AoU")
mgb_dir <- dir_concordance_test(mgb_v, cohort_name = "MGB", beta_col = "BETA_MGB")
gel_dir <- dir_concordance_test(gel_v, cohort_name = "GEL", beta_col = "GEL_BETA")

dir_sum_dt <- rbindlist(list(aou_dir$summary, mgb_dir$summary, gel_dir$summary), fill = TRUE)
print(dir_sum_dt)

############################################################
## Step 9) Visualize directional concordance
############################################################
format_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  formatC(p, format = "e", digits = 2)
}

summarize_dir_two_sided <- function(dt, cohort_label, align_col, p_col, p_thr = 0.05) {
  stopifnot(all(c(align_col, p_col) %in% names(dt)))
  x <- as.data.table(copy(dt))

  all_dt <- x[!is.na(get(align_col))]
  N_all <- nrow(all_dt)
  K_all <- all_dt[get(align_col) == TRUE, .N]
  prop_all <- if (N_all == 0) NA_real_ else (K_all / N_all)

  bt_two <- if (N_all == 0) NULL else binom.test(K_all, N_all, p = 0.5, alternative = "two.sided")
  p_two  <- if (is.null(bt_two)) NA_real_ else unname(bt_two$p.value)
  ci_two <- if (is.null(bt_two)) c(NA_real_, NA_real_) else unname(bt_two$conf.int)

  sub_dt <- x[!is.na(get(p_col)) & get(p_col) < p_thr & !is.na(get(align_col))]
  N_sub <- nrow(sub_dt)
  K_sub <- if (N_sub == 0) 0L else sub_dt[get(align_col) == TRUE, .N]
  prop_sub <- if (N_sub == 0) NA_real_ else (K_sub / N_sub)

  rbindlist(list(
    data.table(
      cohort = cohort_label,
      group = "All genes",
      N = N_all, K = K_all,
      prop = prop_all,
      ci_low = ci_two[1], ci_high = ci_two[2],
      p_value = p_two
    ),
    data.table(
      cohort = cohort_label,
      group = "Validation p<0.05",
      N = N_sub, K = K_sub,
      prop = prop_sub,
      ci_low = NA_real_, ci_high = NA_real_,
      p_value = NA_real_
    )
  ), use.names = TRUE, fill = TRUE)
}

plot_dt <- rbindlist(list(
  summarize_dir_two_sided(aou_v, "AoU", "AoU_beta_align", "AoU_Pvalue_Best_Mask"),
  summarize_dir_two_sided(mgb_v, "MGB", "MGB_beta_align", "MGB_Pvalue_Best_Mask"),
  summarize_dir_two_sided(gel_v, "GEL", "GEL_beta_align", "GEL_Pvalue_Best_Mask")
), use.names = TRUE, fill = TRUE)

plot_dt[, cohort := factor(cohort, levels = c("GEL", "MGB", "AoU"))]
plot_dt[, group := factor(group, levels = c("All genes", "Validation p<0.05"))]
plot_dt[, lab_kn := sprintf("%d/%d", K, N)]
plot_dt[group == "All genes", p_lab := paste0("P=", vapply(p_value, format_p, character(1)))]

dir_sum_dt2 <- plot_dt[group == "All genes",
                      .(cohort, N, K,
                        concordance_rate = prop,
                        p_value,
                        conf_low = ci_low,
                        conf_high = ci_high)]

warm2 <- c(
  "All genes" = "#e6bd5cff",
  "Validation p<0.05" = "#e79653ff"
)

dodge_w <- 0.68
bar_w   <- 0.56

err_dt <- plot_dt[group == "All genes"]
err_dt[, y_err := as.numeric(cohort) - dodge_w / 4]

p_dir_panel <- ggplot(plot_dt, aes(x = prop, y = cohort, fill = group)) +
  geom_col(position = position_dodge(width = dodge_w), width = bar_w) +
  geom_errorbarh(
    data = err_dt,
    aes(y = y_err, xmin = ci_low, xmax = ci_high),
    inherit.aes = FALSE,
    height = 0.18,
    linewidth = 0.35
  ) +
  geom_text(
    aes(label = lab_kn),
    position = position_dodge(width = dodge_w),
    size = 2.6,
    hjust = -0.05
  ) +
  geom_text(
    data = err_dt,
    aes(y = y_err, x = pmin(ci_high + 0.04, 0.98), label = p_lab),
    inherit.aes = FALSE,
    size = 2.6,
    hjust = 0
  ) +
  scale_fill_manual(values = warm2) +
  scale_x_continuous(
    limits = c(0, 1.02),
    breaks = c(0, 0.25, 0.5, 0.75, 1.0),
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  labs(
    x = "Proportion aligned (sign concordance)",
    y = NULL,
    fill = NULL
  ) +
  theme_classic(base_size = 8) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 6.2),
    legend.key.size = unit(3.0, "mm"),
    legend.spacing.x = unit(1.5, "mm"),
    axis.text.y = element_text(size = 8),
    axis.text.x = element_text(size = 7),
    axis.title.x = element_text(size = 8),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

out_pdf2 <- file.path(out_dir, "Fig_directional_concordance_all_vs_validationP005.twosided.pdf")
ggsave(out_pdf2, p_dir_panel, width = 65, height = 70, units = "mm", useDingbats = FALSE)

fwrite(dir_sum_dt2, file.path(out_dir, "DirectionalConcordance_summary_allgenes.twosided.tsv"), sep = "\t")
fwrite(plot_dt, file.path(out_dir, "DirectionalConcordance_plotdata.long.tsv"), sep = "\t")

############################################################
## Step 10) UKBB vs GEL York regression
############################################################
gel_tmp <- copy(gel_best)

if (!("BETA_GEL" %in% names(gel_tmp)) && ("coef" %in% names(gel_tmp))) {
  gel_tmp[, BETA_GEL := as.numeric(coef)]
}
if (!("SE_GEL" %in% names(gel_tmp)) && ("se" %in% names(gel_tmp))) {
  gel_tmp[, SE_GEL := as.numeric(se)]
}

ukb_gel_dt2 <- merge(
  disc_h[, .(
    ENSG,
    gene_symbol,
    x    = as.numeric(BETA_Burden),
    sd_x = as.numeric(SE_Burden)
  )],
  gel_tmp[, .(
    ENSG,
    y    = as.numeric(BETA_GEL),
    sd_y = as.numeric(SE_GEL)
  )],
  by = "ENSG",
  all = FALSE
)

ukb_gel_dt2 <- ukb_gel_dt2[
  is.finite(x) & is.finite(y) &
    is.finite(sd_x) & is.finite(sd_y) &
    sd_x > 0 & sd_y > 0
]

fit_york <- york(
  x = ukb_gel_dt2$x,
  y = ukb_gel_dt2$y,
  sd_x = ukb_gel_dt2$sd_x,
  sd_y = ukb_gel_dt2$sd_y,
  r_xy_errors = 0
)

fit_sum <- summary(fit_york)
coef_mat <- fit_sum$Coefficients

reg_line <- fit_sum$Regression_Test[1]
chisq_val <- suppressWarnings(as.numeric(gsub(".*Chisq-statistic:\\s*([0-9.]+).*", "\\1", reg_line)))
p_val <- suppressWarnings(as.numeric(gsub(".*p-value:\\s*([0-9.eE+-]+).*", "\\1", reg_line)))
if (is.na(p_val) && grepl("p-value:\\s*0\\b", reg_line)) p_val <- 0

york_sum_dt <- data.table(
  cohort_pair  = "UKBB_vs_GEL",
  N            = nrow(ukb_gel_dt2),
  intercept    = unname(coef_mat["intercept", "Estimate"]),
  intercept_se = unname(coef_mat["intercept", "Std_Error"]),
  slope        = unname(coef_mat["slope", "Estimate"]),
  slope_se     = unname(coef_mat["slope", "Std_Error"]),
  chisq        = chisq_val,
  df           = nrow(ukb_gel_dt2) - 2,
  p_value      = p_val
)

print(york_sum_dt)

############################################################
## Step 11) Plot York regression
############################################################
dtp <- as.data.table(ukb_gel_dt2)
b0    <- unname(coef_mat["intercept", "Estimate"])
b0_se <- unname(coef_mat["intercept", "Std_Error"])
b1    <- unname(coef_mat["slope", "Estimate"])
b1_se <- unname(coef_mat["slope", "Std_Error"])

z_slope <- b1 / b1_se
p_slope <- 2 * pnorm(-abs(z_slope))

x_rng <- range(dtp$x, na.rm = TRUE)
line_dt <- data.table(
  x = seq(x_rng[1], x_rng[2], length.out = 200)
)
line_dt[, y := b0 + b1 * x]

annot_txt <- sprintf(
  "York slope = %.3f (SE %.3f)\nP(slope) = %.2e\nIntercept = %.3f (SE %.3f)\nN = %d",
  b1, b1_se, p_slope, b0, b0_se, nrow(dtp)
)

x_annot <- x_rng[1] + 0.02 * diff(x_rng)
y_rng <- range(dtp$y, na.rm = TRUE)
y_annot <- y_rng[2] - 0.02 * diff(y_rng)

theme_smallpanel <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      legend.position = "none",
      axis.title = element_text(size = base_size + 1),
      axis.text  = element_text(size = base_size),
      plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "mm")
    )
}

p_scatter_york <- ggplot(dtp, aes(x = x, y = y)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey35") +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.35, color = "grey35") +
  geom_errorbarh(
    aes(xmin = x - sd_x, xmax = x + sd_x),
    height = 0, linewidth = 0.20, alpha = 0.55, color = "grey40"
  ) +
  geom_errorbar(
    aes(ymin = y - sd_y, ymax = y + sd_y),
    width = 0, linewidth = 0.20, alpha = 0.55, color = "grey40"
  ) +
  geom_point(aes(color = x), size = 1.25, alpha = 0.9) +
  scale_color_gradientn(
    colors = c("#fff5eb", "#fdd0a2", "#fdae6b", "#fd8d3c", "#e6550d", "#a63603"),
    guide = "none"
  ) +
  geom_line(
    data = line_dt, aes(x = x, y = y),
    inherit.aes = FALSE,
    linewidth = 0.7, color = "black"
  ) +
  labs(
    x = "UKBB effect (beta)",
    y = "GEL effect (beta)"
  ) +
  annotate(
    "text",
    x = x_annot, y = y_annot,
    label = annot_txt,
    hjust = 0, vjust = 1,
    size = 2.7
  ) +
  theme_smallpanel(base_size = 8)

out_pdf3 <- file.path(out_dir, "Fig_York_UKBB_vs_GEL.smallpanel.pdf")
ggsave(
  filename = out_pdf3,
  plot = p_scatter_york,
  width = 75,
  height = 70,
  units = "mm",
  useDingbats = FALSE
)

############################################################
## Step 12) Merge AoU, MGB, and GEL validation tables
############################################################
aou_all <- as.data.table(copy(aou_v))
mgb_all <- as.data.table(copy(mgb_v))
gel_all <- as.data.table(copy(gel_v))

aou_all[, gene_symbol := as.character(gene_symbol)]
mgb_all[, gene_symbol := as.character(gene_symbol)]
gel_all[, gene_symbol := as.character(gene_symbol)]

merged_all <- merge(
  aou_all, mgb_all,
  by = "gene_symbol", all = TRUE,
  suffixes = c(".AoUtbl", ".MGBtbl")
)
merged_all <- merge(
  merged_all, gel_all,
  by = "gene_symbol", all = TRUE,
  suffixes = c("", ".GELtbl")
)

nom_cols <- c("AoU_nominal_validation", "MGB_nominal_validation", "GEL_nominal_validation")
for (cc in nom_cols) {
  if (cc %in% names(merged_all)) {
    merged_all[, (cc) := as.logical(get(cc))]
  }
}

if (all(nom_cols %in% names(merged_all))) {
  merged_all[
    , n_nominal_validated :=
        as.integer(AoU_nominal_validation %in% TRUE) +
        as.integer(MGB_nominal_validation %in% TRUE) +
        as.integer(GEL_nominal_validation %in% TRUE)
  ]
  setorder(merged_all, -n_nominal_validated, gene_symbol)
} else {
  setorder(merged_all, gene_symbol)
}

out_tsv <- file.path(out_dir, "ValidationTables_AoU_MGB_GEL.merged_allcols.tsv")
fwrite(merged_all, out_tsv, sep = "\t", quote = FALSE, na = "NA")

cat("[INFO] Saved merged table:\n  ", out_tsv, "\n", sep = "")