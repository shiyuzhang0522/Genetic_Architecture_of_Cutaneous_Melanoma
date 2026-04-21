############################################################
## Spatial-context heterogeneity in tumor regions
## and its association with gsMap signals
##
## Author: Shelley
## Date:   2026-02-12
############################################################

############################################################
## Step 1) Load packages
############################################################
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(fs)
})

############################################################
## Step 2) Define sample IDs and read per-sample data
##   Inputs:
##     - gsMap grid table
##     - SpaCET interface labels
##     - SpaCET cell-type proportions
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
gsmap_dir  <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/gsMAP"
spacet_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/SpaCET.outputs"

## helper: safe fread with clear message
safe_fread <- function(path, ...) {
  if (!file.exists(path)) {
    warning("[WARN] Missing file: ", path, call. = FALSE)
    return(NULL)
  }
  fread(path, data.table = TRUE, ...)
}

## initialize lists
gsmap_list <- setNames(vector("list", length(sample_ids)), sample_ids)
iface_list <- setNames(vector("list", length(sample_ids)), sample_ids)
prop_list  <- setNames(vector("list", length(sample_ids)), sample_ids)

## read files into lists
for (sid in sample_ids) {
  gsmap_path <- file.path(gsmap_dir, paste0(sid, "_gsMap_SpaCETgrid_table.txt"))
  iface_path <- file.path(spacet_dir, sid, "tables", paste0(sid, ".SpaCET.interface_labels.tsv"))
  prop_path  <- file.path(spacet_dir, sid, "tables", paste0(sid, ".SpaCET.celltype_prop.tsv"))

  cat("[INFO] Reading: ", sid, "\n")
  cat("       gsMap:   ", gsmap_path, "\n")
  cat("       iface:   ", iface_path, "\n")
  cat("       prop:    ", prop_path, "\n")

  gsmap_list[[sid]] <- safe_fread(gsmap_path)
  iface_list[[sid]] <- safe_fread(iface_path)
  prop_list[[sid]]  <- safe_fread(prop_path)
}

## quick QC summary
qc_dt <- data.table(
  sample_id = sample_ids,
  gsmap_ok  = vapply(gsmap_list, \(x) !is.null(x), logical(1)),
  iface_ok  = vapply(iface_list, \(x) !is.null(x), logical(1)),
  prop_ok   = vapply(prop_list,  \(x) !is.null(x), logical(1))
)

cat("\n[INFO] Read summary (TRUE = loaded):\n")
print(qc_dt)

############################################################
## Step 3) Merge gsMap, interface labels, and cell-type
##         proportions for each sample
##
## Keys:
##   - gsMap: SpaCET_position
##   - iface: spot_id
##   - props: spot_id
############################################################

merge_one_sample <- function(sid, gsmap_list, iface_list, prop_list) {
  g <- gsmap_list[[sid]]
  i <- iface_list[[sid]]
  p <- prop_list[[sid]]

  if (is.null(g) || is.null(i) || is.null(p)) {
    warning("[WARN] Skip (NULL input): ", sid, call. = FALSE)
    return(NULL)
  }

  ## sanity checks
  stopifnot("SpaCET_position" %in% names(g))
  stopifnot("spot_id" %in% names(i))
  stopifnot("spot_id" %in% names(p))

  g <- as.data.table(g)
  i <- as.data.table(i)
  p <- as.data.table(p)

  ## de-duplicate keys if needed
  if (anyDuplicated(g$SpaCET_position)) {
    warning("[WARN] Duplicated gsMap SpaCET_position in ", sid, " (keeping first).", call. = FALSE)
    g <- g[!duplicated(SpaCET_position)]
  }
  if (anyDuplicated(i$spot_id)) {
    warning("[WARN] Duplicated iface spot_id in ", sid, " (keeping first).", call. = FALSE)
    i <- i[!duplicated(spot_id)]
  }
  if (anyDuplicated(p$spot_id)) {
    warning("[WARN] Duplicated prop spot_id in ", sid, " (keeping first).", call. = FALSE)
    p <- p[!duplicated(spot_id)]
  }

  ## rename keys to match gsMap
  setnames(i, "spot_id", "SpaCET_position")
  setnames(p, "spot_id", "SpaCET_position")

  ## set keys
  setkey(g, SpaCET_position)
  setkey(i, SpaCET_position)
  setkey(p, SpaCET_position)

  ## merge
  m <- merge(g, i, by = "SpaCET_position", all.x = TRUE)
  m <- merge(m, p, by = "SpaCET_position", all.x = TRUE)

  ## add sample_id
  m[, sample_id := sid]

  ## quick merge QC
  n_all   <- nrow(m)
  n_iface <- sum(!is.na(m$interface_label))
  cat(sprintf("[INFO] Merged %-35s  n=%d | iface=%d (%.1f%%)\n",
              sid, n_all, n_iface, 100 * n_iface / n_all))

  m[]
}

merged_list <- setNames(vector("list", length(sample_ids)), sample_ids)

for (sid in sample_ids) {
  merged_list[[sid]] <- merge_one_sample(sid, gsmap_list, iface_list, prop_list)
}

############################################################
## Step 4) Compute minimum distance to the nearest interface
##         for tumor spots
############################################################

add_min_dist_to_interface_one_sample <- function(dt,
                                                 pos_col = "SpaCET_position",
                                                 verbose = FALSE) {
  dt <- copy(as.data.table(dt))
  stopifnot("interface_label" %in% names(dt), pos_col %in% names(dt))

  ## standardized label column
  dt[, interface := trimws(as.character(interface_label))]

  ## parse "76x22" -> array_row and array_col
  tmp <- tstrsplit(dt[[pos_col]], "x", fixed = TRUE)
  dt[, array_row := as.integer(tmp[[1]])]
  dt[, array_col := as.integer(tmp[[2]])]

  tumor_dt <- dt[interface == "Tumor"]
  int_dt   <- dt[interface == "Interface"]

  dt[, min_dist_to_interface := NA_real_]

  if (verbose) {
    cat("[QC] labels: ", paste(sort(unique(dt$interface)), collapse = " | "), "\n", sep = "")
    cat("[QC] n_tumor: ", nrow(tumor_dt), "  n_interface: ", nrow(int_dt), "\n", sep = "")
  }

  if (nrow(tumor_dt) > 0L && nrow(int_dt) > 0L) {
    tumor_coords <- as.matrix(tumor_dt[, .(array_row, array_col)])
    int_coords   <- as.matrix(int_dt[, .(array_row, array_col)])

    dmin <- vapply(seq_len(nrow(tumor_coords)), function(i) {
      dr <- tumor_coords[i, 1] - int_coords[, 1]
      dc <- tumor_coords[i, 2] - int_coords[, 2]
      min(sqrt(dr^2 + dc^2))
    }, numeric(1))

    tumor_dt[, min_dist_to_interface := dmin]
    dt[tumor_dt, min_dist_to_interface := i.min_dist_to_interface, on = pos_col]

  } else {
    if (verbose) {
      if (nrow(tumor_dt) == 0L) cat("[WARN] No Tumor spots. Distances remain NA.\n")
      if (nrow(int_dt) == 0L)   cat("[WARN] No Interface spots. Distances remain NA.\n")
    }
  }

  dt[]
}

for (sid in sample_ids) {
  if (is.null(merged_list[[sid]])) {
    warning("[WARN] merged_list[[", sid, "]] is NULL. Skipping.", call. = FALSE)
    next
  }

  merged_list[[sid]] <- add_min_dist_to_interface_one_sample(
    merged_list[[sid]],
    verbose = FALSE
  )

  dt <- merged_list[[sid]]
  n_tumor <- dt[interface == "Tumor", .N]
  n_intf  <- dt[interface == "Interface", .N]
  med_d   <- dt[interface == "Tumor", median(min_dist_to_interface, na.rm = TRUE)]
  max_d   <- dt[interface == "Tumor", max(min_dist_to_interface, na.rm = TRUE)]

  cat(sprintf("[INFO] %-35s  Tumor=%4d  Interface=%4d  medianDist=%.3f  maxDist=%.3f\n",
              sid, n_tumor, n_intf, med_d, max_d))
}

############################################################
## Step 5) Per-sample Spearman summary
##   Tumor spots: distance to interface vs gsMap signal
############################################################

suppressPackageStartupMessages({
  library(data.table)
})

spearman_summary <- function(x, y, conf_level = 0.95) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  n <- length(x)

  if (n < 10L) {
    return(list(n = n, rho = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, p = NA_real_))
  }

  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  rho <- unname(ct$estimate)
  p   <- ct$p.value

  alpha <- 1 - conf_level
  zcrit <- qnorm(1 - alpha / 2)

  rho_clip <- max(min(rho, 0.999999), -0.999999)
  z <- atanh(rho_clip)
  se <- 1 / sqrt(n - 3)

  z_lo <- z - zcrit * se
  z_hi <- z + zcrit * se

  ci_lo <- tanh(z_lo)
  ci_hi <- tanh(z_hi)

  list(n = n, rho = rho, ci_lo = ci_lo, ci_hi = ci_hi, p = p)
}

## choose gsMap signal column
signal_col <- "logp"

corr_sum_dt <- rbindlist(lapply(sample_ids, function(sid) {
  dt <- merged_list[[sid]]
  if (is.null(dt)) return(NULL)

  if (!("interface" %in% names(dt))) {
    dt[, interface := trimws(as.character(interface_label))]
  }

  tumor_dt <- dt[
    interface == "Tumor" &
      !is.na(min_dist_to_interface) &
      !is.na(get(signal_col))
  ]

  x <- tumor_dt$min_dist_to_interface
  y <- tumor_dt[[signal_col]]

  s <- spearman_summary(x, y, conf_level = 0.95)

  data.table(
    sample_id = sid,
    signal_col = signal_col,
    n_tumor_used = s$n,
    spearman_rho = s$rho,
    spearman_ci95_lo = s$ci_lo,
    spearman_ci95_hi = s$ci_hi,
    spearman_p = s$p
  )
}), use.names = TRUE, fill = TRUE)

corr_sum_dt[, spearman_p_BH := p.adjust(spearman_p, method = "BH")]

corr_sum_dt[
  ,
  direction := fifelse(
    is.na(spearman_rho),
    NA_character_,
    fifelse(spearman_rho > 0, "positive", "negative")
  )
]

setorder(corr_sum_dt, spearman_p_BH)
corr_sum_dt

############################################################
## Step 6) Per-sample linear model
##   logp ~ min_dist_to_interface + Malignant
##   Tumor spots only
############################################################

suppressPackageStartupMessages({
  library(data.table)
})

signal_col <- "logp"

lm_sum_dt <- rbindlist(lapply(sample_ids, function(sid) {
  dt <- merged_list[[sid]]
  if (is.null(dt)) return(NULL)

  if (!("interface" %in% names(dt))) {
    dt[, interface := trimws(as.character(interface_label))]
  }

  dtt <- dt[
    interface == "Tumor" &
      !is.na(min_dist_to_interface) &
      !is.na(get(signal_col)) &
      !is.na(Malignant)
  ]

  n <- nrow(dtt)
  if (n < 20L) {
    return(data.table(
      sample_id = sid,
      n_tumor_used = n,
      beta_dist = NA_real_,
      ci95_lo = NA_real_,
      ci95_hi = NA_real_,
      p_dist = NA_real_
    ))
  }

  fit <- lm(
    dtt[[signal_col]] ~ min_dist_to_interface + Malignant,
    data = dtt
  )

  sm <- summary(fit)
  ci <- confint(fit)

  data.table(
    sample_id = sid,
    n_tumor_used = n,
    beta_dist = sm$coefficients["min_dist_to_interface", "Estimate"],
    ci95_lo   = ci["min_dist_to_interface", 1],
    ci95_hi   = ci["min_dist_to_interface", 2],
    p_dist    = sm$coefficients["min_dist_to_interface", "Pr(>|t|)"]
  )
}), use.names = TRUE, fill = TRUE)

lm_sum_dt[, p_dist_BH := p.adjust(p_dist, method = "BH")]

lm_sum_dt[
  ,
  direction := fifelse(
    is.na(beta_dist),
    NA_character_,
    fifelse(beta_dist > 0, "positive", "negative")
  )
]

setorder(lm_sum_dt, p_dist_BH)
lm_sum_dt

############################################################
## Step 7) Save result tables
############################################################

suppressPackageStartupMessages({
  library(data.table)
})

out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/6.spatial.heterogeneity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

corr_out <- file.path(out_dir, "Tumor_dist_vs_gsMap.corr_summary.tsv")
lm_out   <- file.path(out_dir, "Tumor_dist_vs_gsMap.lm_summary.adjustMalignantProp.tsv")

stopifnot(
  exists("corr_sum_dt"),
  exists("lm_sum_dt")
)

fwrite(corr_sum_dt, file = corr_out, sep = "\t", na = "NA")
fwrite(lm_sum_dt, file = lm_out, sep = "\t", na = "NA")

cat("[INFO] Saved result tables:\n")
cat("  - ", corr_out, "\n", sep = "")
cat("  - ", lm_out,   "\n", sep = "")

############################################################
## Step 8) Visualization: per-sample scatter and fitted line
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

signal_col <- "logp"

malignant_cols <- c("#FDE0DD", "#FA9FB5", "#C51B8A", "#7A0177")

plot_dist_vs_gsmap_one_sample <- function(dt, sid) {
  dt <- as.data.table(dt)

  if (!("interface" %in% names(dt))) {
    dt[, interface := trimws(as.character(interface_label))]
  }

  dtt <- dt[
    interface == "Tumor" &
      !is.na(min_dist_to_interface) &
      !is.na(get(signal_col)) &
      !is.na(Malignant)
  ]

  if (nrow(dtt) < 20L) {
    warning("[WARN] Too few tumor spots for plotting: ", sid, call. = FALSE)
    return(NULL)
  }

  ct <- suppressWarnings(
    cor.test(
      dtt$min_dist_to_interface,
      dtt[[signal_col]],
      method = "spearman",
      exact = FALSE
    )
  )
  rho  <- unname(ct$estimate)
  pval <- ct$p.value

  ann_lab <- sprintf("Spearman \u03C1 = %.3f\np = %.2g", rho, pval)

  ggplot(
    dtt,
    aes(
      x = min_dist_to_interface,
      y = .data[[signal_col]],
      color = Malignant
    )
  ) +
    geom_point(size = 0.6, alpha = 0.75) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.5, color = "black") +
    scale_color_gradientn(colors = malignant_cols, name = "Malignant\nproportion") +
    labs(
      title = sid,
      x = "Distance to nearest interface (grid units)",
      y = expression("gsMap signal  (" * -log[10] * "(p))")
    ) +
    annotate(
      "text",
      x = Inf, y = Inf,
      label = ann_lab,
      hjust = 1.05, vjust = 1.15,
      size = 2,
      color = "black"
    ) +
    theme_classic(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 7, hjust = 0),
      axis.title = element_text(size = 6),
      axis.text  = element_text(size = 6),
      legend.title = element_text(size = 6),
      legend.text  = element_text(size = 6),
      legend.key.height = unit(0.3, "cm")
    )
}

plot_list <- setNames(vector("list", length(sample_ids)), sample_ids)

for (sid in sample_ids) {
  dt <- merged_list[[sid]]
  if (is.null(dt)) next
  plot_list[[sid]] <- plot_dist_vs_gsmap_one_sample(dt, sid)
}

############################################################
## Step 9) Combine plots with custom ordering
##   - negatives: ascending
##   - positives: descending
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(patchwork)
})

stopifnot(exists("corr_sum_dt"), exists("plot_list"))

ord_dt <- corr_sum_dt[!is.na(spearman_rho), .(sample_id, spearman_rho)]

neg_ids <- ord_dt[spearman_rho < 0, sample_id][order(ord_dt[spearman_rho < 0, spearman_rho])]
pos_ids <- ord_dt[spearman_rho >= 0, sample_id][order(-ord_dt[spearman_rho >= 0, spearman_rho])]

ord_samples <- c(neg_ids, pos_ids)

cat("[INFO] Custom sample order:\n")
print(corr_sum_dt[match(ord_samples, sample_id), .(sample_id, spearman_rho)])

ordered_plots <- plot_list[ord_samples]
ordered_plots <- ordered_plots[!vapply(ordered_plots, is.null, logical(1))]

n_col <- 4

combined_plot <-
  wrap_plots(ordered_plots, ncol = n_col) +
  plot_annotation(
    title = "Spatial association between gsMap signal and tumor depth",
    subtitle = "Negatives (most negative→0) followed by positives (largest→smallest)",
    theme = theme(
      plot.title = element_text(face = "bold", size = 8),
      plot.subtitle = element_text(size = 6)
    )
  )

combined_plot

############################################################
## Step 10) Save combined figure
############################################################

suppressPackageStartupMessages({
  library(patchwork)
  library(ggplot2)
})

out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/6.spatial.heterogeneity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_pdf <- file.path(out_dir, "Fig_TumorDist_vs_gsMap_logp.bySample.order_custom.A4panel.pdf")

a4_w_mm <- 200
panel_h_mm <- 180

ggsave(
  filename = out_pdf,
  plot = combined_plot,
  width = a4_w_mm,
  height = panel_h_mm,
  units = "mm",
  useDingbats = FALSE
)