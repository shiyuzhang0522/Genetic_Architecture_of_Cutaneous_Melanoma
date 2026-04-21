############################################################
## gsMap signal heterogeneity across tumor spots
## Adjusting for spot-level UMI depth
## Author: Shelley
## Date:   2026-02-12
############################################################

############################################################
## Step 1) Load packages
############################################################
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
})

############################################################
## Step 2) Read per-sample gsMap and QC tables
############################################################

## sample IDs
sample_ids <- c(
  "MEL109",
  "MEL126",
  "MEL131",
  "MEL162",
  "MEL170",
  "MEL176",
  "MEL201",
  "MEL222",
  "MEL239",
  "MEL256",
  "MEL271",
  "MEL74",
  "MEL75",
  "MEL99"
)

cat("[INFO] Number of samples: ", length(sample_ids), "\n")

## base directories
gsmap_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/gsMAP"
qc_base   <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/SpaCET.outputs"

## build per-sample file paths
gsmap_paths <- file.path(gsmap_dir, paste0(sample_ids, "_gsMap_SpaCETgrid_table.txt"))
qc_paths    <- file.path(qc_base, sample_ids, "tables", paste0(sample_ids, ".spot_QC_metrics.tsv"))

## sanity checks
missing_gsmap <- gsmap_paths[!file.exists(gsmap_paths)]
missing_qc    <- qc_paths[!file.exists(qc_paths)]

if (length(missing_gsmap) > 0) {
  cat("[WARN] Missing gsMap files:\n")
  cat("  - ", paste(missing_gsmap, collapse = "\n  - "), "\n")
}
if (length(missing_qc) > 0) {
  cat("[WARN] Missing QC metric files:\n")
  cat("  - ", paste(missing_qc, collapse = "\n  - "), "\n")
}

stopifnot(length(missing_gsmap) == 0, length(missing_qc) == 0)

## read into lists
gsmap_list <- setNames(vector("list", length(sample_ids)), sample_ids)
qc_list    <- setNames(vector("list", length(sample_ids)), sample_ids)

for (sid in sample_ids) {
  gsmap_path <- file.path(gsmap_dir, paste0(sid, "_gsMap_SpaCETgrid_table.txt"))
  qc_path    <- file.path(qc_base, sid, "tables", paste0(sid, ".spot_QC_metrics.tsv"))

  dt_gs <- fread(gsmap_path)
  dt_gs[, sample_id := sid]

  dt_qc <- fread(qc_path)
  dt_qc[, sample_id := sid]

  gsmap_list[[sid]] <- dt_gs
  qc_list[[sid]]    <- dt_qc

  cat("[INFO] Loaded ", sid, " | gsMap: ", nrow(dt_gs), " rows; QC: ", nrow(dt_qc), " rows\n")
}

## quick peek
sid0 <- sample_ids[1]
cat("\n[INFO] Example gsMap columns (first 25):\n")
print(head(names(gsmap_list[[sid0]]), 25))

cat("\n[INFO] Example QC columns (first 25):\n")
print(head(names(qc_list[[sid0]]), 25))

cat("\n[INFO] Example gsMap rows:\n")
print(gsmap_list[[sid0]][1:3])

cat("\n[INFO] Example QC rows:\n")
print(qc_list[[sid0]][1:3])

############################################################
## Step 3) Merge gsMap signal with spot-level QC metrics
############################################################

merged_list <- vector("list", length(sample_ids))
names(merged_list) <- sample_ids

for (sid in sample_ids) {
  dt_gs <- copy(gsmap_list[[sid]])
  dt_qc <- copy(qc_list[[sid]])

  ## harmonize spot identifier
  setnames(dt_gs, "spot", "barcode")

  ## define gsMap signal
  stopifnot("logp" %in% names(dt_gs))
  dt_gs <- dt_gs[, .(
    barcode,
    gsMap_log10p = as.numeric(logp),
    p_value      = as.numeric(p),
    z            = as.numeric(z),
    beta         = as.numeric(beta),
    se           = as.numeric(se),
    SpaCET_position,
    sample_id
  )]

  ## QC table: keep UMI and coordinates
  stopifnot(all(c("barcode", "UMI") %in% names(dt_qc)))
  keep_qc <- intersect(
    c("barcode", "UMI", "log1p_UMI", "nGene", "pixel_row", "pixel_col",
      "array_row", "array_col", "coordinate_x_um", "coordinate_y_um"),
    names(dt_qc)
  )
  dt_qc <- dt_qc[, ..keep_qc]

  ## merge
  dt_m <- merge(dt_gs, dt_qc, by = "barcode", all.x = TRUE)

  miss_umi <- mean(is.na(dt_m$UMI))
  if (miss_umi > 0) {
    cat("[WARN] ", sid, ": missing UMI for ", sprintf("%.2f", 100 * miss_umi), "% spots\n", sep = "")
  }

  merged_list[[sid]] <- dt_m
  cat("[INFO] Merged ", sid, ": ", nrow(dt_m), " spots\n", sep = "")
}

############################################################
## Step 4) Per-sample UMI adjustment
############################################################

resid_list   <- vector("list", length(sample_ids))
fit_sum_list <- vector("list", length(sample_ids))
names(resid_list)   <- sample_ids
names(fit_sum_list) <- sample_ids

for (sid in sample_ids) {
  dt <- copy(merged_list[[sid]])

  stopifnot(all(c("gsMap_log10p", "UMI") %in% names(dt)))

  if ("log1p_UMI" %in% names(dt) && !all(is.na(dt$log1p_UMI))) {
    dt[, logUMI := as.numeric(log1p_UMI)]
  } else {
    dt[, logUMI := log1p(as.numeric(UMI))]
  }

  dt_fit <- dt[is.finite(gsMap_log10p) & is.finite(logUMI)]
  if (nrow(dt_fit) < 50) {
    warning("[WARN] Too few spots after filtering for sample: ", sid)
    next
  }

  fit <- lm(gsMap_log10p ~ logUMI, data = dt_fit)

  dt[, gsMap_resid := NA_real_]
  dt[dt_fit, gsMap_resid := resid(fit), on = "barcode"]

  dt[, gsMap_fit := NA_real_]
  dt[dt_fit, gsMap_fit := fitted(fit), on = "barcode"]

  s <- summary(fit)
  fit_sum_list[[sid]] <- data.table(
    sample_id    = sid,
    n_spots_fit  = nrow(dt_fit),
    slope_logUMI = unname(coef(fit)["logUMI"]),
    intercept    = unname(coef(fit)["(Intercept)"]),
    r2           = unname(s$r.squared),
    adj_r2       = unname(s$adj.r.squared),
    p_logUMI     = unname(coef(s)["logUMI", "Pr(>|t|)"]),
    cor_spearman = suppressWarnings(cor(dt_fit$gsMap_log10p, dt_fit$logUMI, method = "spearman")),
    cor_pearson  = suppressWarnings(cor(dt_fit$gsMap_log10p, dt_fit$logUMI, method = "pearson"))
  )

  resid_list[[sid]] <- dt

  cat("[INFO] Step 4 done: ", sid,
      " | n=", nrow(dt_fit),
      " | slope=", signif(fit_sum_list[[sid]]$slope_logUMI, 3),
      " | R2=", signif(fit_sum_list[[sid]]$r2, 3),
      "\n", sep = "")
}

fit_sum_dt <- rbindlist(fit_sum_list, use.names = TRUE, fill = TRUE)

cat("\n[INFO] Fit summary table (first 10 rows):\n")
print(fit_sum_dt[1:min(10, .N)])

############################################################
## Save per-sample fit summaries
############################################################

out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/5.gsMap.signal.heterogeneity"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

fit_sum_path_tsv <- file.path(out_dir, "gsMap_UMI_adjustment_fit_summary.tsv")
fwrite(fit_sum_dt, fit_sum_path_tsv, sep = "\t", quote = FALSE)

############################################################
## Step 5) Read SpaCET labels and cell-type proportions
############################################################

spacet_base <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/SpaCET.outputs"

label_list <- setNames(vector("list", length(sample_ids)), sample_ids)
prop_list  <- setNames(vector("list", length(sample_ids)), sample_ids)

for (sid in sample_ids) {
  label_path <- file.path(spacet_base, sid, "tables", paste0(sid, ".SpaCET.interface_labels.tsv"))
  prop_path  <- file.path(spacet_base, sid, "tables", paste0(sid, ".SpaCET.celltype_prop.tsv"))

  stopifnot(file.exists(label_path), file.exists(prop_path))

  label_list[[sid]] <- fread(label_path)
  prop_list[[sid]]  <- fread(prop_path)

  cat("[INFO] Loaded ", sid,
      " | labels: ", nrow(label_list[[sid]]), " x ", ncol(label_list[[sid]]),
      " | props: ", nrow(prop_list[[sid]]), " x ", ncol(prop_list[[sid]]),
      "\n", sep = "")
}

sid0 <- sample_ids[1]
cat("\n[INFO] Example label columns (first 25):\n")
print(head(names(label_list[[sid0]]), 25))
cat("\n[INFO] Example prop columns (first 25):\n")
print(head(names(prop_list[[sid0]]), 25))

############################################################
## Step 6) Merge labels and cell-type proportions into resid_list
############################################################

for (sid in sample_ids) {
  dt_main  <- copy(resid_list[[sid]])
  dt_label <- copy(label_list[[sid]])
  dt_prop  <- copy(prop_list[[sid]])

  stopifnot("SpaCET_position" %in% names(dt_main))
  stopifnot("spot_id" %in% names(dt_label))
  stopifnot("spot_id" %in% names(dt_prop))

  if (anyDuplicated(dt_label$spot_id)) {
    warning("[WARN] duplicated spot_id in label file for: ", sid)
  }
  if (anyDuplicated(dt_prop$spot_id)) {
    warning("[WARN] duplicated spot_id in prop file for: ", sid)
  }

  dt_main <- merge(
    dt_main,
    dt_label,
    by.x = "SpaCET_position",
    by.y = "spot_id",
    all.x = TRUE
  )

  dt_main <- merge(
    dt_main,
    dt_prop,
    by.x = "SpaCET_position",
    by.y = "spot_id",
    all.x = TRUE
  )

  miss_lab <- mean(is.na(dt_main$interface_label))
  miss_mal <- if ("Malignant" %in% names(dt_main)) mean(is.na(dt_main$Malignant)) else NA_real_

  cat("[INFO] ", sid,
      " | merged rows: ", nrow(dt_main),
      " | missing interface_label: ", sprintf("%.2f%%", 100 * miss_lab),
      " | missing Malignant: ", ifelse(is.na(miss_mal), "NA", sprintf("%.2f%%", 100 * miss_mal)),
      "\n", sep = "")

  resid_list[[sid]] <- dt_main
}

cat("\n[INFO] Step 6 complete: merged interface_label and cell-type proportions by SpaCET_position.\n")

sid0 <- sample_ids[1]
print(resid_list[[sid0]][1:3, .(SpaCET_position, barcode, interface_label, Malignant, CAF, UMI, gsMap_log10p, gsMap_resid)])

############################################################
## Step 7) Define tumor purity flags
############################################################

cutoff_loose  <- 0.5
cutoff_strict <- 0.99

for (sid in sample_ids) {
  dt <- copy(resid_list[[sid]])

  stopifnot(all(c("interface_label", "Malignant") %in% names(dt)))

  dt[, tumor_pure_loose  := (interface_label == "Tumor" & Malignant >= cutoff_loose)]
  dt[, tumor_pure_strict := (interface_label == "Tumor" & Malignant >= cutoff_strict)]

  chk <- dt[, .N, by = .(interface_label, tumor_pure_loose, tumor_pure_strict)]
  setorder(chk, interface_label, tumor_pure_loose, tumor_pure_strict)

  cat("\n[INFO] ", sid, " | counts by interface_label and purity flags:\n", sep = "")
  print(chk)

  resid_list[[sid]] <- dt
}

cat("\n[INFO] Step 7 complete: added tumor_pure_loose and tumor_pure_strict.\n")

############################################################
## Step 8) Summary table of gsMap_log10p heterogeneity
############################################################

het_list <- vector("list", length(sample_ids))
names(het_list) <- sample_ids

for (sid in sample_ids) {
  dt <- resid_list[[sid]]

  dt_A <- dt[interface_label == "Tumor" & is.finite(gsMap_log10p)]
  dt_B <- dt[interface_label == "Tumor" & Malignant >= 0.99 & is.finite(gsMap_log10p)]

  if (nrow(dt_A) < 30) next

  het_list[[sid]] <- data.table(
    sample_id = sid,
    n_tumor_all    = nrow(dt_A),
    n_tumor_strict = nrow(dt_B),

    mean_all   = mean(dt_A$gsMap_log10p),
    median_all = median(dt_A$gsMap_log10p),
    sd_all     = sd(dt_A$gsMap_log10p),
    IQR_all    = IQR(dt_A$gsMap_log10p),
    MAD_all    = mad(dt_A$gsMap_log10p),
    var_all    = var(dt_A$gsMap_log10p),

    mean_strict   = ifelse(nrow(dt_B) > 10, mean(dt_B$gsMap_log10p), NA_real_),
    median_strict = ifelse(nrow(dt_B) > 10, median(dt_B$gsMap_log10p), NA_real_),
    sd_strict     = ifelse(nrow(dt_B) > 10, sd(dt_B$gsMap_log10p), NA_real_),
    IQR_strict    = ifelse(nrow(dt_B) > 10, IQR(dt_B$gsMap_log10p), NA_real_),
    MAD_strict    = ifelse(nrow(dt_B) > 10, mad(dt_B$gsMap_log10p), NA_real_),
    var_strict    = ifelse(nrow(dt_B) > 10, var(dt_B$gsMap_log10p), NA_real_)
  )
}

het_dt <- rbindlist(het_list, use.names = TRUE, fill = TRUE)
print(het_dt)

############################################################
## Step 9) Summary table of gsMap_resid heterogeneity
############################################################

het_resid_list <- vector("list", length(sample_ids))
names(het_resid_list) <- sample_ids

for (sid in sample_ids) {
  dt <- resid_list[[sid]]

  dt_A <- dt[
    interface_label == "Tumor" &
      is.finite(gsMap_resid)
  ]

  dt_B <- dt[
    interface_label == "Tumor" &
      Malignant >= 0.99 &
      is.finite(gsMap_resid)
  ]

  if (nrow(dt_A) < 30) next

  het_resid_list[[sid]] <- data.table(
    sample_id = sid,
    n_tumor_all    = nrow(dt_A),
    n_tumor_strict = nrow(dt_B),

    mean_all_resid   = mean(dt_A$gsMap_resid),
    median_all_resid = median(dt_A$gsMap_resid),
    sd_all_resid     = sd(dt_A$gsMap_resid),
    IQR_all_resid    = IQR(dt_A$gsMap_resid),
    MAD_all_resid    = mad(dt_A$gsMap_resid),
    var_all_resid    = var(dt_A$gsMap_resid),

    mean_strict_resid   = ifelse(nrow(dt_B) > 10, mean(dt_B$gsMap_resid), NA_real_),
    median_strict_resid = ifelse(nrow(dt_B) > 10, median(dt_B$gsMap_resid), NA_real_),
    sd_strict_resid     = ifelse(nrow(dt_B) > 10, sd(dt_B$gsMap_resid), NA_real_),
    IQR_strict_resid    = ifelse(nrow(dt_B) > 10, IQR(dt_B$gsMap_resid), NA_real_),
    MAD_strict_resid    = ifelse(nrow(dt_B) > 10, mad(dt_B$gsMap_resid), NA_real_),
    var_strict_resid    = ifelse(nrow(dt_B) > 10, var(dt_B$gsMap_resid), NA_real_)
  )
}

het_resid_dt <- rbindlist(het_resid_list, use.names = TRUE, fill = TRUE)
print(het_resid_dt)

############################################################
## Save summary tables
############################################################

out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/5.gsMap.signal.heterogeneity"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

raw_tsv <- file.path(out_dir, "gsMap_logp_heterogeneity_summary.tsv")
fwrite(het_dt, raw_tsv, sep = "\t", quote = FALSE)

resid_tsv <- file.path(out_dir, "gsMap_residual_heterogeneity_summary.tsv")
fwrite(het_resid_dt, resid_tsv, sep = "\t", quote = FALSE)

############################################################
## Step 10) A4 plotting
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

make_two_panel_boxplot <- function(dt_all,
                                   y_col,
                                   y_label,
                                   main_title_prefix,
                                   order_by = c("strict", "all"),
                                   strict_cutoff = 0.99,
                                   min_strict_n = 30,
                                   palette_fill = c("#BFD7EA", "#D8C7F0"),
                                   palette_line = c("#1F4E79", "#5A2A82"),
                                   dashed_zero = FALSE,
                                   y_symmetric = FALSE) {

  order_by <- match.arg(order_by)

  stopifnot(all(c("sample_id", "interface_label", "Malignant", y_col) %in% names(dt_all)))

  dtA <- dt_all[interface_label == "Tumor" & is.finite(get(y_col))]
  dtB <- dt_all[interface_label == "Tumor" & Malignant >= strict_cutoff & is.finite(get(y_col))]

  ordB <- dtB[, .(var_strict = var(get(y_col), na.rm = TRUE), n_strict = .N), by = sample_id]
  ordA <- dtA[, .(var_all = var(get(y_col), na.rm = TRUE), n_all = .N), by = sample_id]
  ord  <- merge(ordA, ordB, by = "sample_id", all = TRUE)

  ord[, var_for_order := fifelse(!is.na(var_strict) & n_strict >= min_strict_n, var_strict, var_all)]
  setorder(ord, -var_for_order)
  lvls <- ord$sample_id

  dtA[, sample_id := factor(sample_id, levels = lvls)]
  dtB[, sample_id := factor(sample_id, levels = lvls)]

  yA <- dtA[[y_col]]
  if (y_symmetric) {
    q <- quantile(yA, probs = c(0.01, 0.99), na.rm = TRUE)
    y_abs <- max(abs(q))
    y_abs <- max(1, ceiling(y_abs * 2) / 2)
    y_lim <- c(-y_abs, y_abs)
  } else {
    y_max <- quantile(yA, probs = 0.995, na.rm = TRUE)
    y_max <- max(6, ceiling(y_max * 2) / 2)
    y_lim <- c(0, y_max)
  }

  theme_ng <- theme_classic(base_size = 10.5) +
    theme(
      plot.title   = element_text(face = "bold", size = 11.5, hjust = 0),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 10.5),
      axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
      panel.grid   = element_blank(),
      plot.margin  = margin(3, 5, 2, 5, "mm")
    )

  pA <- ggplot(dtA, aes(x = sample_id, y = .data[[y_col]])) +
    (if (dashed_zero) geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) else NULL) +
    geom_boxplot(
      width = 0.55,
      outlier.shape = 16, outlier.size = 0.6, outlier.alpha = 0.20,
      fill = palette_fill[1], color = palette_line[1], linewidth = 0.32
    ) +
    geom_point(
      position = position_jitter(width = 0.18, height = 0),
      size = 0.45, alpha = 0.10, color = palette_line[1]
    ) +
    coord_cartesian(ylim = y_lim) +
    labs(
      title = paste0("A  |  ", main_title_prefix, " in all Tumor spots"),
      y = y_label
    ) +
    theme_ng +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  pB <- ggplot(dtB, aes(x = sample_id, y = .data[[y_col]])) +
    (if (dashed_zero) geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) else NULL) +
    geom_boxplot(
      width = 0.55,
      outlier.shape = 16, outlier.size = 0.6, outlier.alpha = 0.20,
      fill = palette_fill[2], color = palette_line[2], linewidth = 0.32
    ) +
    geom_point(
      position = position_jitter(width = 0.18, height = 0),
      size = 0.45, alpha = 0.10, color = palette_line[2]
    ) +
    coord_cartesian(ylim = y_lim) +
    labs(
      title = paste0("B  |  ", main_title_prefix, " in Tumor spots with Malignant ≥ ", strict_cutoff),
      y = y_label
    ) +
    theme_ng

  pA / pB + plot_layout(heights = c(1, 1))
}

dt_all <- rbindlist(resid_list, use.names = TRUE, fill = TRUE)
stopifnot(all(c("sample_id", "interface_label", "Malignant", "gsMap_log10p", "gsMap_resid") %in% names(dt_all)))

p_raw <- make_two_panel_boxplot(
  dt_all = dt_all,
  y_col = "gsMap_log10p",
  y_label = "gsMap logp",
  main_title_prefix = "gsMap logp",
  strict_cutoff = 0.99,
  min_strict_n = 30,
  dashed_zero = FALSE,
  y_symmetric = FALSE
)

p_resid <- make_two_panel_boxplot(
  dt_all = dt_all,
  y_col = "gsMap_resid",
  y_label = "gsMap residual",
  main_title_prefix = "UMI-adjusted gsMap (residual)",
  strict_cutoff = 0.99,
  min_strict_n = 30,
  dashed_zero = TRUE,
  y_symmetric = TRUE
)

library(patchwork)

p_A4 <- wrap_plots(p_raw, p_resid, ncol = 1, guides = "collect") +
  plot_layout(heights = c(1, 1))

print(p_A4)

out_pdf <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/5.gsMap.signal.heterogeneity/Fig_gsMap_raw_and_residual_distribution.pdf"
ggsave(out_pdf, p_A4, width = 190, height = 240, units = "mm")

############################################################
## Step 11) Moran's I test for spatial heterogeneity
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(spdep)
})

moran_test_one_sample <- function(
  dt,
  sid,
  tumor_filter_expr,
  coord_x = "array_col",
  coord_y = "array_row",
  k = 6,
  n_perm = 999
) {
  dt_sub <- dt[
    eval(tumor_filter_expr) &
      is.finite(gsMap_resid) &
      is.finite(get(coord_x)) &
      is.finite(get(coord_y))
  ]

  if (nrow(dt_sub) < 50) return(NULL)

  coords <- as.matrix(dt_sub[, .(get(coord_x), get(coord_y))])
  resids <- dt_sub$gsMap_resid

  knn <- spdep::knearneigh(coords, k = k)
  nb  <- spdep::knn2nb(knn)
  lw  <- spdep::nb2listw(nb, style = "W", zero.policy = TRUE)

  mor <- spdep::moran.mc(
    x = resids,
    listw = lw,
    nsim = n_perm,
    zero.policy = TRUE,
    alternative = "greater"
  )

  data.table(
    sample_id   = sid,
    tumor_def   = NA_character_,
    n_spots     = nrow(dt_sub),
    k_neighbors = k,
    moran_I     = unname(mor$statistic),
    p_value     = unname(mor$p.value),
    expectation = unname(mor$estimate[["Expectation"]]),
    variance    = unname(mor$estimate[["Variance"]])
  )
}

k_use <- 6
nperm_use <- 999

moran_results <- list()

for (sid in sample_ids) {
  dt <- resid_list[[sid]]

  stopifnot(all(c("gsMap_resid", "interface_label", "Malignant", "array_row", "array_col") %in% names(dt)))

  resA <- moran_test_one_sample(
    dt = dt,
    sid = sid,
    tumor_filter_expr = quote(interface_label == "Tumor"),
    coord_x = "array_col",
    coord_y = "array_row",
    k = k_use,
    n_perm = nperm_use
  )
  if (!is.null(resA)) {
    resA[, tumor_def := "Tumor_all"]
    moran_results[[paste0(sid, "_A")]] <- resA
  }

  strict_expr <- if ("tumor_pure_strict" %in% names(dt)) {
    quote(tumor_pure_strict == TRUE)
  } else {
    quote(interface_label == "Tumor" & Malignant >= 0.99)
  }

  resB <- moran_test_one_sample(
    dt = dt,
    sid = sid,
    tumor_filter_expr = strict_expr,
    coord_x = "array_col",
    coord_y = "array_row",
    k = k_use,
    n_perm = nperm_use
  )
  if (!is.null(resB)) {
    resB[, tumor_def := "Tumor_strict"]
    moran_results[[paste0(sid, "_B")]] <- resB
  }

  cat("[INFO] Moran's I done: ", sid,
      " | all: ", ifelse(is.null(resA), "NA", signif(resA$moran_I, 3)),
      " p=", ifelse(is.null(resA), "NA", signif(resA$p_value, 3)),
      " | strict: ", ifelse(is.null(resB), "NA", signif(resB$moran_I, 3)),
      " p=", ifelse(is.null(resB), "NA", signif(resB$p_value, 3)),
      "\n", sep = "")
}

moran_dt <- rbindlist(moran_results, use.names = TRUE, fill = TRUE)
moran_dt[, p_FDR := p.adjust(p_value, method = "BH")]

setorder(moran_dt, tumor_def, p_FDR)
print(moran_dt)

############################################################
## Save Moran's I summary table
############################################################

out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/5.gsMap.signal.heterogeneity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_path <- file.path(out_dir, "gsMap_MoransI_spatial_heterogeneity_summary.tsv")

fwrite(
  moran_dt,
  file = out_path,
  sep = "\t",
  quote = FALSE
)

cat("[INFO] Moran's I summary table saved to:\n", out_path, "\n")