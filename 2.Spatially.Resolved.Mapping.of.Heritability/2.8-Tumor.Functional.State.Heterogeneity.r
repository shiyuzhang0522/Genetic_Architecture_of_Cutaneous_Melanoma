############################################################
## Cell functional state heterogeneity
## Association between gsMap signals & tumor-spot cell status
## 
## Author: Shelley
## Last update:   2026-03-18
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
##   Inputs per sample:
##     1) gsMap grid table
##     2) SpaCET interface labels
##     3) SpaCET gene-set score table
############################################################

## samples
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

## helper: construct per-sample file paths
make_paths_one_sample <- function(sid, gsmap_dir, spacet_dir) {
  list(
    sid = sid,
    gsmap_txt = file.path(gsmap_dir, paste0(sid, "_gsMap_SpaCETgrid_table.txt")),
    interface_tsv = file.path(spacet_dir, sid, "tables", paste0(sid, ".SpaCET.interface_labels.tsv")),
    status_tsv = file.path(spacet_dir, sid, "tables", paste0(sid, ".SpaCET.GeneSetScore.TLS.tsv"))
  )
}

paths_list <- lapply(
  sample_ids,
  make_paths_one_sample,
  gsmap_dir = gsmap_dir,
  spacet_dir = spacet_dir
)
names(paths_list) <- sample_ids

## helper: safe fread with existence check
safe_fread <- function(fp, ...) {
  if (!file.exists(fp)) {
    stop("[ERROR] File not found: ", fp)
  }
  data.table::fread(fp, ...)
}

## read per-sample tables
gsmap_list     <- vector("list", length(sample_ids)); names(gsmap_list) <- sample_ids
interface_list <- vector("list", length(sample_ids)); names(interface_list) <- sample_ids
status_list    <- vector("list", length(sample_ids)); names(status_list) <- sample_ids

for (sid in sample_ids) {
  p <- paths_list[[sid]]

  gsmap_list[[sid]]     <- safe_fread(p$gsmap_txt)
  interface_list[[sid]] <- safe_fread(p$interface_tsv)
  status_list[[sid]]    <- safe_fread(p$status_tsv)

  cat(
    "[INFO] Loaded: ", sid, "\n",
    "  - gsMap:      ", basename(p$gsmap_txt), "\n",
    "  - interface:  ", basename(p$interface_tsv), "\n",
    "  - status:     ", basename(p$status_tsv), "\n",
    sep = ""
  )
}

cat(
  "[INFO] Done. Objects in memory:\n",
  "  - gsmap_list\n",
  "  - interface_list\n",
  "  - status_list\n",
  sep = ""
)

str(gsmap_list[[1]])
str(interface_list[[1]])
str(status_list[[1]])

############################################################
## Step 3) Merge per-sample tables
##   Key mapping:
##     - gsMap:     SpaCET_position
##     - interface: spot_id
##     - status:    spot_id
############################################################

merge_one_sample <- function(sid, gs_dt, interface_dt, status_dt) {
  gs <- copy(gs_dt)
  it <- copy(interface_dt)
  st <- copy(status_dt)

  ## sanity checks
  stopifnot(all(c("spot", "SpaCET_position") %in% names(gs)))
  stopifnot(all(c("spot_id", "interface_label") %in% names(it)))
  stopifnot("spot_id" %in% names(st))

  ## standardize join keys
  gs[, SpaCET_position := as.character(SpaCET_position)]
  it[, spot_id := as.character(spot_id)]
  st[, spot_id := as.character(spot_id)]

  ## enforce uniqueness in annotation tables
  if (anyDuplicated(it$spot_id)) {
    warning("[WARN] Duplicate spot_id in interface table for ", sid,
            " (keeping first occurrence).")
    it <- it[!duplicated(spot_id)]
  }
  if (anyDuplicated(st$spot_id)) {
    warning("[WARN] Duplicate spot_id in status table for ", sid,
            " (keeping first occurrence).")
    st <- st[!duplicated(spot_id)]
  }

  ## merge gsMap + interface
  m1 <- merge(
    x = gs,
    y = it,
    by.x = "SpaCET_position",
    by.y = "spot_id",
    all.x = TRUE,
    sort = FALSE
  )

  ## merge + status
  merged <- merge(
    x = m1,
    y = st,
    by.x = "SpaCET_position",
    by.y = "spot_id",
    all.x = TRUE,
    sort = FALSE
  )

  ## add sample_id
  merged[, sample_id := sid]
  setcolorder(merged, c("sample_id", "spot", "SpaCET_position"))

  ## quick QC
  n0 <- nrow(gs)
  n1 <- nrow(merged)
  if (n0 != n1) {
    warning("[WARN] Row count changed after merge for ", sid,
            ": gsMap n=", n0, " -> merged n=", n1,
            ". This usually indicates duplicated keys in gsMap.")
  }

  miss_it <- merged[is.na(interface_label), .N]
  if (miss_it > 0) {
    warning("[WARN] Missing interface_label for ", sid, ": ", miss_it, " spots.")
  }

  merged[]
}

merged_list <- vector("list", length(sample_ids))
names(merged_list) <- sample_ids

for (sid in sample_ids) {
  merged_list[[sid]] <- merge_one_sample(
    sid = sid,
    gs_dt = gsmap_list[[sid]],
    interface_dt = interface_list[[sid]],
    status_dt = status_list[[sid]]
  )
  cat("[INFO] Merged: ", sid, " | n=", nrow(merged_list[[sid]]), "\n", sep = "")
}

cat("[INFO] Step 3 done. Output object: merged_list\n")

############################################################
## Step 4) Explore the association between gsMap logp and
##         tumor-spot status scores
############################################################

############################################################
## Step 4.1) Subset tumor spots
############################################################

tumor_list <- vector("list", length(sample_ids))
names(tumor_list) <- sample_ids

for (sid in sample_ids) {
  dt <- copy(merged_list[[sid]])
  tumor_dt <- dt[interface_label == "Tumor"]

  tumor_list[[sid]] <- tumor_dt

  cat(
    "[INFO] ", sid,
    " | total spots: ", nrow(dt),
    " | tumor spots: ", nrow(tumor_dt), "\n",
    sep = ""
  )
}

n_tumor_total <- sum(sapply(tumor_list, nrow))
cat("[INFO] Total tumor spots across all samples: ", n_tumor_total, "\n")

############################################################
## Step 4.2) Prepare hallmark feature list and correlation helper
############################################################

min_tumor_spots <- 30

## identify hallmark columns shared across samples
hallmark_cols_list <- lapply(
  tumor_list,
  function(dt) grep("^HALLMARK_", names(dt), value = TRUE)
)
hallmark_cols <- Reduce(intersect, hallmark_cols_list)
hallmark_cols <- sort(hallmark_cols)

cat("[INFO] Hallmark features found (intersection across samples): ",
    length(hallmark_cols), "\n")
if (length(hallmark_cols) > 0) {
  cat("[INFO] Example hallmark columns: ",
      paste(head(hallmark_cols, 5), collapse = ", "), "\n")
}

gsmap_y <- "logp"

has_logp <- sapply(tumor_list, function(dt) gsmap_y %in% names(dt))
if (!all(has_logp)) {
  bad <- names(has_logp)[!has_logp]
  stop("[ERROR] Missing '", gsmap_y, "' in these samples: ", paste(bad, collapse = ", "))
}

safe_spearman <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 5) {
    return(list(rho = NA_real_, p = NA_real_, n = length(x)))
  }
  if (stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(list(rho = NA_real_, p = NA_real_, n = length(x)))
  }

  ct <- suppressWarnings(stats::cor.test(x, y, method = "spearman", exact = FALSE))
  list(rho = unname(ct$estimate), p = ct$p.value, n = length(x))
}

############################################################
## Step 4.3) Per-sample Spearman correlation
##   gsMap logp vs Hallmark scores in tumor spots
############################################################

hallmark_cor_list <- vector("list", length(sample_ids))
names(hallmark_cor_list) <- sample_ids

for (sid in sample_ids) {
  dt <- tumor_list[[sid]]

  if (nrow(dt) < min_tumor_spots) {
    warning("[WARN] ", sid, ": tumor spots < min_tumor_spots (n=",
            nrow(dt), "). Returning NA results.")
    tmp <- data.table(
      sample_id = sid,
      feature = hallmark_cols,
      rho = NA_real_,
      p = NA_real_,
      n = nrow(dt)
    )
    hallmark_cor_list[[sid]] <- tmp
    next
  }

  res_list <- lapply(hallmark_cols, function(feat) {
    out <- safe_spearman(dt[[gsmap_y]], dt[[feat]])
    data.table(
      sample_id = sid,
      feature = feat,
      rho = out$rho,
      p = out$p,
      n = out$n
    )
  })

  res_dt <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
  hallmark_cor_list[[sid]] <- res_dt

  cat("[INFO] ", sid, ": computed ", nrow(res_dt),
      " hallmark correlations (n tumor spots=", nrow(dt), ")\n", sep = "")
}

hallmark_cor_dt <- rbindlist(hallmark_cor_list, use.names = TRUE, fill = TRUE)

cat("[INFO] Combined hallmark correlation table (raw p only): ",
    nrow(hallmark_cor_dt), " rows\n", sep = "")

cat("[INFO] Strongest raw associations (top 10 by |rho|, excluding NA):\n")
print(hallmark_cor_dt[!is.na(rho)][order(-abs(rho))][1:10])

############################################################
## Step 4.4) Global BH correction across all samples and pathways
############################################################

hallmark_cor_dt[, p_bh := stats::p.adjust(p, method = "BH")]

############################################################
## Step 5) Extract significant positive associations
############################################################

hallmark_pos_sig_dt <- hallmark_cor_dt[
  !is.na(rho) & !is.na(p_bh) &
    rho > 0 & p_bh < 0.05
]

cat("[INFO] Number of significant positive associations: ",
    nrow(hallmark_pos_sig_dt), "\n", sep = "")

for (sid in sample_ids) {
  tmp <- hallmark_pos_sig_dt[sample_id == sid]

  cat("\n============================================================\n")
  cat("[INFO] Sample: ", sid, "\n", sep = "")
  cat("[INFO] Number of significant positive pathways: ", nrow(tmp), "\n", sep = "")

  if (nrow(tmp) == 0) {
    cat("[INFO] No significant positive pathways detected.\n")
  } else {
    print(tmp[, .(feature, rho, p, p_bh, n)][order(-rho, p_bh)])
  }
}

############################################################
## Step 5.1) Summarize significant positive associations
##           across samples for each Hallmark pathway
############################################################

hallmark_summary_dt <- hallmark_pos_sig_dt[
  ,
  .(
    n_samples_sig = uniqueN(sample_id),
    n_total = .N,
    mean_rho = mean(rho, na.rm = TRUE),
    median_rho = median(rho, na.rm = TRUE),
    sd_rho = sd(rho, na.rm = TRUE),
    min_rho = min(rho, na.rm = TRUE),
    max_rho = max(rho, na.rm = TRUE)
  ),
  by = feature
]

setorder(hallmark_summary_dt, -n_samples_sig, -median_rho)

cat("[INFO] Per-pathway summary of significant positive associations:\n")
print(hallmark_summary_dt)

############################################################
## Step 5.2) Visualize significant positive Hallmark correlations
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

############################################################
## Step 5.2.1) Settings
############################################################

min_sig_samples <- 6
max_samples <- length(sample_ids)

hall_pal <- c("#C4D9EC", "#9ECAE1", "#6BAED6", "#3182BD", "#08306B")

############################################################
## Step 5.2.2) Helpers
############################################################

get_rho_limits <- function(x, pad = 0.12) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(-0.1, 0.1))
  lo <- min(x)
  hi <- max(x)
  rng <- hi - lo
  if (!is.finite(rng) || rng == 0) rng <- max(abs(c(lo, hi)), 0.1)
  lo <- lo - pad * rng
  hi <- hi + pad * rng
  lo <- min(lo, 0)
  hi <- max(hi, 0)
  c(lo, hi)
}

clean_hallmark_label <- function(x) {
  x <- gsub("^HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  tools::toTitleCase(tolower(x))
}

############################################################
## Step 5.2.3) Build per-pathway summary table for plotting
############################################################

stopifnot(all(c("sample_id", "feature", "rho", "p", "p_bh", "n") %in% names(hallmark_pos_sig_dt)))

hallmark_plot_long_dt <- copy(hallmark_pos_sig_dt)

hallmark_plot_summary_dt <- hallmark_plot_long_dt[
  ,
  .(
    n_sig_samples = uniqueN(sample_id),
    mean_rho_sig = mean(rho, na.rm = TRUE),
    median_rho_sig = median(rho, na.rm = TRUE),
    sd_rho_sig = sd(rho, na.rm = TRUE),
    q1_rho_sig = as.numeric(quantile(rho, 0.25, na.rm = TRUE)),
    q3_rho_sig = as.numeric(quantile(rho, 0.75, na.rm = TRUE)),
    min_rho_sig = min(rho, na.rm = TRUE),
    max_rho_sig = max(rho, na.rm = TRUE)
  ),
  by = feature
]

hallmark_plot_summary_dt <- hallmark_plot_summary_dt[n_sig_samples >= min_sig_samples]

cat("[INFO] Number of Hallmark pathways passing min_sig_samples = ",
    min_sig_samples, ": ", nrow(hallmark_plot_summary_dt), "\n", sep = "")

if (nrow(hallmark_plot_summary_dt) == 0) {
  stop("[ERROR] No Hallmark pathways passed the minimum sample threshold.")
}

hallmark_plot_long_dt <- hallmark_plot_long_dt[feature %in% hallmark_plot_summary_dt$feature]

############################################################
## Step 5.2.4) Clean labels and define ordering
############################################################

hallmark_plot_summary_dt[, feature_label := clean_hallmark_label(feature)]
hallmark_plot_long_dt[, feature_label := clean_hallmark_label(feature)]

setorder(hallmark_plot_summary_dt, -n_sig_samples, -median_rho_sig)
feat_levels <- rev(hallmark_plot_summary_dt$feature_label)

hallmark_plot_summary_dt[, feature_label := factor(feature_label, levels = feat_levels)]
hallmark_plot_long_dt[, feature_label := factor(feature_label, levels = feat_levels)]

hallmark_plot_long_dt <- merge(
  hallmark_plot_long_dt,
  hallmark_plot_summary_dt[, .(feature, median_rho_sig)],
  by = "feature",
  all.x = TRUE,
  sort = FALSE
)

############################################################
## Step 5.2.5) Define scales
############################################################

col_lim <- range(hallmark_plot_summary_dt$median_rho_sig, na.rm = TRUE)
if (!all(is.finite(col_lim)) || diff(col_lim) == 0) {
  col_lim <- c(min(0, col_lim[1]), max(0.3, col_lim[2]))
}

rho_lim <- get_rho_limits(hallmark_plot_long_dt$rho, pad = 0.14)

max_median_rho <- max(hallmark_plot_summary_dt$median_rho_sig, na.rm = TRUE)
cat("[INFO] Maximum median rho across Hallmark pathways: ",
    sprintf("%.2f", max_median_rho), "\n", sep = "")

############################################################
## Step 5.2.6) Left panel: number of significant samples
############################################################

p_bar <- ggplot(
  hallmark_plot_summary_dt,
  aes(x = -n_sig_samples, y = feature_label, fill = median_rho_sig)
) +
  geom_col(width = 0.70) +
  geom_vline(
    xintercept = -min_sig_samples,
    linetype = "dashed",
    linewidth = 0.40,
    colour = "grey50"
  ) +
  scale_x_continuous(
    limits = c(-max_samples, 0),
    breaks = -seq(0, max_samples, by = 2),
    labels = function(x) abs(x),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_fill_gradientn(
    colours = hall_pal,
    limits = col_lim,
    oob = scales::squish,
    name = "Median rho\nacross sig. samples"
  ) +
  labs(
    title = "Hallmark",
    x = "N samples",
    y = NULL
  ) +
  theme_minimal(base_size = 6) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 6, face = "bold", hjust = 0, colour = "black"),
    axis.title.x = element_text(size = 6, margin = margin(t = 5), colour = "black"),
    axis.text.y  = element_text(size = 6, colour = "black"),
    axis.text.x  = element_text(size = 6, colour = "black"),
    axis.ticks   = element_line(linewidth = 0.26, colour = "black"),
    axis.ticks.length = unit(2, "pt"),
    legend.position = "none",
    plot.margin = margin(2, 1, 2, 2)
  )

############################################################
## Step 5.2.7) Right panel: rho distribution across samples
############################################################

p_rho <- ggplot() +
  geom_linerange(
    data = hallmark_plot_summary_dt,
    aes(
      y = feature_label,
      xmin = q1_rho_sig,
      xmax = q3_rho_sig,
      colour = median_rho_sig
    ),
    linewidth = 0.55,
    alpha = 0.70
  ) +
  geom_segment(
    data = hallmark_plot_summary_dt,
    aes(
      y = feature_label,
      yend = feature_label,
      x = median_rho_sig,
      xend = median_rho_sig,
      colour = median_rho_sig
    ),
    linewidth = 1.00,
    lineend = "round"
  ) +
  geom_point(
    data = hallmark_plot_long_dt,
    aes(
      y = feature_label,
      x = rho,
      colour = median_rho_sig
    ),
    position = position_jitter(width = 0, height = 0.11, seed = 1),
    size = 0.9,
    alpha = 0.80
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    colour = "grey70"
  ) +
  scale_colour_gradientn(
    colours = hall_pal,
    limits = col_lim,
    oob = scales::squish,
    name = "Median rho\nacross sig. samples",
    guide = guide_colorbar(
      title.position = "top",
      barheight = unit(28, "mm"),
      barwidth = unit(2.8, "mm"),
      frame.colour = "black",
      ticks.colour = "black"
    )
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.08))
  ) +
  coord_cartesian(xlim = rho_lim, clip = "on") +
  labs(
    x = "Spearman rho",
    y = NULL
  ) +
  theme_minimal(base_size = 6) +
  theme(
    panel.grid = element_blank(),
    axis.title.x = element_text(size = 6, margin = margin(t = 5), colour = "black"),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_text(size = 5.5, colour = "black"),
    axis.ticks   = element_line(linewidth = 0.26, colour = "black"),
    axis.ticks.length = unit(2, "pt"),
    axis.line    = element_line(linewidth = 0.26, colour = "black"),
    legend.position = "right",
    legend.title = element_text(size = 5.5, colour = "black"),
    legend.text  = element_text(size = 5.0, colour = "black"),
    plot.margin = margin(2, 2, 2, 0)
  )

############################################################
## Step 5.2.8) Combine and print
############################################################

p_hallmark_only <- p_bar + p_rho + patchwork::plot_layout(widths = c(0.60, 1.2))
print(p_hallmark_only)

############################################################
## Step 5.3) Save Hallmark-only figure
############################################################

out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsMap_0318/cell.status.heterogeneity"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_pdf <- file.path(
  out_dir,
  "Fig_gsmap_sigpos_hallmark_summary.only.withLegend.pdf"
)

ggsave(
  filename = out_pdf,
  plot = p_hallmark_only,
  width = 100,
  height = 100,
  units = "mm"
)

cat("[INFO] Saved:\n  ", out_pdf, "\n", sep = "")