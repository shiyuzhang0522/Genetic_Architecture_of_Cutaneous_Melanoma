############################################################
## gsMap + ST analysis
## Global pattern of spatially resolved mapping of heritability
## Regional Cauchy p-values and cell-type correlation analysis
## Author:  Shelley
## Last updated:  2026-02-11
############################################################

############################################################
## Step 1) Load packages
############################################################
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(forcats)
})

cat("[INFO] Packages loaded successfully.\n")

############################################################
## Step 2) Read gsMap and SpaCET outputs
############################################################
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

gsmap_dir  <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/gsMAP"
spacet_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/SpaCET.outputs"

stopifnot(dir.exists(gsmap_dir), dir.exists(spacet_dir))

read_spacet_interface <- function(sample_id, spacet_dir) {
  f <- file.path(
    spacet_dir, sample_id, "tables",
    paste0(sample_id, ".SpaCET.interface_labels.tsv")
  )
  if (!file.exists(f)) {
    warning("[WARN] Missing SpaCET interface file: ", f)
    return(NULL)
  }
  dt <- fread(f)
  dt[, sample_id := sample_id]
  dt
}

read_gsmap_pvals <- function(sample_id, gsmap_dir) {
  f <- file.path(gsmap_dir, paste0(sample_id, "_gsMap_SpaCETgrid_table.txt"))
  if (!file.exists(f)) {
    warning("[WARN] Missing gsMap grid table: ", f)
    return(NULL)
  }
  dt <- fread(f)
  dt[, sample_id := sample_id]
  dt
}

spacet_interface_list <- setNames(vector("list", length(sample_ids)), sample_ids)
gsmap_list <- setNames(vector("list", length(sample_ids)), sample_ids)

for (sid in sample_ids) {
  cat("[INFO] Reading sample: ", sid, "\n")
  spacet_interface_list[[sid]] <- read_spacet_interface(sid, spacet_dir)
  gsmap_list[[sid]] <- read_gsmap_pvals(sid, gsmap_dir)
}

n_spacet <- sum(!vapply(spacet_interface_list, is.null, logical(1)))
n_gsmap  <- sum(!vapply(gsmap_list, is.null, logical(1)))

cat("[INFO] SpaCET interface files loaded: ", n_spacet, "/", length(sample_ids), "\n")
cat("[INFO] gsMap grid tables loaded: ", n_gsmap, "/", length(sample_ids), "\n")

############################################################
## Step 3) Merge interface labels with gsMap spot-level p-values
############################################################
merge_one_sample <- function(sid, spacet_interface_list, gsmap_list) {
  sp_dt <- spacet_interface_list[[sid]]
  gs_dt <- gsmap_list[[sid]]

  if (is.null(sp_dt) || is.null(gs_dt)) {
    warning("[WARN] Skip merge (missing input): ", sid)
    return(NULL)
  }

  stopifnot(all(c("spot_id", "interface_label", "sample_id") %in% names(sp_dt)))
  stopifnot(all(c("SpaCET_position", "p", "sample_id") %in% names(gs_dt)))

  sp_dt2 <- copy(sp_dt)
  gs_dt2 <- copy(gs_dt)

  sp_dt2[, SpaCET_position := spot_id]
  sp_dt2[, SpaCET_position := str_trim(SpaCET_position)]
  gs_dt2[, SpaCET_position := str_trim(SpaCET_position)]

  if (any(duplicated(sp_dt2$SpaCET_position))) {
    warning("[WARN] Duplicated SpaCET_position in SpaCET table: ", sid, " (keeping first occurrence)")
    sp_dt2 <- sp_dt2[!duplicated(SpaCET_position)]
  }

  if (any(duplicated(gs_dt2$SpaCET_position))) {
    warning("[WARN] Duplicated SpaCET_position in gsMap table: ", sid, " (keeping first occurrence)")
    gs_dt2 <- gs_dt2[!duplicated(SpaCET_position)]
  }

  m_dt <- merge(
    x = gs_dt2,
    y = sp_dt2[, .(SpaCET_position, interface_label)],
    by = "SpaCET_position",
    all.x = TRUE,
    sort = FALSE
  )

  n_miss <- sum(is.na(m_dt$interface_label))
  if (n_miss > 0) {
    warning("[WARN] Missing interface_label after merge: ", sid, " (", n_miss, " spots)")
  }

  m_dt[, interface_label := factor(interface_label, levels = c("Tumor", "Interface", "Stroma"))]
  m_dt[]
}

merged_list <- setNames(vector("list", length(sample_ids)), sample_ids)

for (sid in sample_ids) {
  cat("[INFO] Merging sample: ", sid, "\n")
  merged_list[[sid]] <- merge_one_sample(sid, spacet_interface_list, gsmap_list)
}

############################################################
## Step 4) Region-level Cauchy combined p-values
############################################################
cauchy_combine_p <- function(p, w = NULL, eps = 1e-15) {
  p <- as.numeric(p)
  p <- p[is.finite(p) & !is.na(p)]
  if (length(p) == 0) return(NA_real_)

  p <- pmin(pmax(p, eps), 1 - eps)

  if (is.null(w)) {
    w <- rep(1 / length(p), length(p))
  } else {
    w <- as.numeric(w)
    stopifnot(length(w) == length(p))
    w <- w / sum(w)
  }

  T <- sum(w * tan((0.5 - p) * pi))
  p_region <- 0.5 - atan(T) / pi
  p_region <- pmin(pmax(p_region, 0), 1)

  p_region
}

region_cauchy_list <- setNames(vector("list", length(sample_ids)), sample_ids)

for (sid in sample_ids) {
  dt <- merged_list[[sid]]

  if (is.null(dt)) {
    warning("[WARN] Skip (missing merged dt): ", sid)
    region_cauchy_list[[sid]] <- NULL
    next
  }

  stopifnot(all(c("interface_label", "p") %in% names(dt)))

  dt2 <- dt[!is.na(interface_label) & !is.na(p)]

  out_dt <- dt2[, .(
    n_spots = .N,
    p_cauchy = cauchy_combine_p(p)
  ), by = .(interface_label)]

  out_dt[, interface_label := factor(interface_label, levels = c("Tumor", "Interface", "Stroma"))]
  setorder(out_dt, interface_label)
  out_dt[, neglog10_p := -log10(pmax(p_cauchy, .Machine$double.xmin))]
  out_dt[, sample_id := sid]
  setcolorder(out_dt, c("sample_id", "interface_label", "n_spots", "p_cauchy", "neglog10_p"))

  region_cauchy_list[[sid]] <- out_dt
  cat("[INFO] Done: ", sid, " | regions: ", nrow(out_dt), "\n")
}

region_cauchy_summary_dt <- rbindlist(region_cauchy_list, use.names = TRUE, fill = TRUE)
cat("[INFO] Region-level summary rows: ", nrow(region_cauchy_summary_dt), "\n")

############################################################
## Step 5) Plot regional Cauchy signals
############################################################
stopifnot(all(c("sample_id", "interface_label", "neglog10_p") %in% names(region_cauchy_summary_dt)))

region_levels <- c("Tumor", "Interface", "Stroma")
region_cauchy_summary_dt[, interface_label := factor(interface_label, levels = region_levels)]

region_cols <- c(
  "Tumor"     = "#b20000ff",
  "Interface" = "#E69F00",
  "Stroma"    = "#4b85b4ff"
)

thr_dt <- region_cauchy_summary_dt[!is.na(neglog10_p), .(
  med_neglog10_p = median(neglog10_p, na.rm = TRUE)
), by = .(interface_label)]

thr_dt[, interface_label := factor(interface_label, levels = region_levels)]
setorder(thr_dt, interface_label)

cat("[INFO] Region median thresholds (-log10 p):\n")
print(thr_dt)

out_pdf <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/1.Cauchy/Fig_region_cauchy_smallpanel.pdf"

ymax <- max(region_cauchy_summary_dt$neglog10_p, na.rm = TRUE)
thr_dt[, y_lab := sprintf("%.2f", med_neglog10_p)]
thr_dt[, y_text := pmin(med_neglog10_p + 0.06, ymax - 0.02)]

p_region <- ggplot(
  region_cauchy_summary_dt,
  aes(x = interface_label, y = neglog10_p, group = sample_id)
) +
  geom_line(linewidth = 0.35, alpha = 0.35, colour = "grey60") +
  geom_point(aes(colour = interface_label), size = 2.0, alpha = 0.95) +
  geom_hline(
    data = thr_dt,
    aes(yintercept = med_neglog10_p, colour = interface_label),
    linetype = "dashed",
    linewidth = 0.65,
    alpha = 0.95,
    inherit.aes = FALSE
  ) +
  geom_text(
    data = thr_dt,
    aes(x = 1.0, y = y_text, label = y_lab, colour = interface_label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 0,
    size = 2.6
  ) +
  scale_colour_manual(values = region_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(x = "Region", y = expression(-log[10](Cauchy~p))) +
  theme_minimal(base_size = 8.5) +
  theme(
    legend.position = "none",
    panel.grid.major = element_line(linewidth = 0.32, colour = "grey88"),
    panel.grid.minor = element_blank(),
    axis.title.x = element_text(size = 8.5, margin = margin(t = 5), colour = "black"),
    axis.title.y = element_text(size = 8.5, margin = margin(r = 5), colour = "black"),
    axis.text.x  = element_text(size = 8.5, colour = "black"),
    axis.text.y  = element_text(size = 7.5, colour = "black"),
    axis.ticks   = element_line(linewidth = 0.30, colour = "black"),
    axis.ticks.length = unit(2, "pt"),
    axis.line    = element_line(linewidth = 0.28, colour = "black"),
    plot.margin  = margin(3, 5, 3, 3)
  ) +
  coord_cartesian(clip = "off")

ggsave(out_pdf, p_region, width = 70, height = 60, units = "mm", useDingbats = FALSE)
cat("[INFO] Saved:\n  ", out_pdf, "\n")

############################################################
## Step 6) Save region-level summary table
############################################################
out_dir  <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/1.Cauchy"
out_file <- file.path(out_dir, "region_cauchy_summary_dt.txt")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
  cat("[INFO] Created directory: ", out_dir, "\n")
}

fwrite(
  region_cauchy_summary_dt,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("[INFO] Saved region_cauchy_summary_dt →\n  ", out_file, "\n")

############################################################
## Step 7) Read SpaCET cell-type proportions
############################################################
spacet_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/SpaCET.outputs"
stopifnot(dir.exists(spacet_dir))

read_spacet_celltype_prop <- function(sample_id, spacet_dir) {
  f <- file.path(
    spacet_dir, sample_id, "tables",
    paste0(sample_id, ".SpaCET.celltype_prop.tsv")
  )

  if (!file.exists(f)) {
    warning("[WARN] Missing SpaCET celltype_prop file: ", f)
    return(NULL)
  }

  dt <- fread(f)
  dt[, sample_id := sample_id]

  has_spot <- any(c("spot_id", "SpaCET_position") %in% names(dt))
  if (!has_spot) {
    warning("[WARN] No spot id column found in: ", f,
            "\n[WARN] Columns are: ", paste(names(dt), collapse = ", "))
  }

  dt[]
}

spacet_celltype_list <- setNames(vector("list", length(sample_ids)), sample_ids)

for (sid in sample_ids) {
  cat("[INFO] Reading SpaCET celltype_prop: ", sid, "\n")
  spacet_celltype_list[[sid]] <- read_spacet_celltype_prop(sid, spacet_dir)
}

n_loaded <- sum(!vapply(spacet_celltype_list, is.null, logical(1)))
cat("[INFO] SpaCET celltype_prop files loaded: ", n_loaded, "/", length(sample_ids), "\n")

############################################################
## Step 8) Merge gsMap with cell-type proportions
############################################################
merge_gsmap_with_celltype <- function(sid, merged_list, spacet_celltype_list) {
  gdt <- merged_list[[sid]]
  cdt <- spacet_celltype_list[[sid]]

  if (is.null(gdt) || is.null(cdt)) {
    warning("[WARN] Skip merge (missing input): ", sid)
    return(NULL)
  }

  stopifnot("SpaCET_position" %in% names(gdt))
  stopifnot("spot_id" %in% names(cdt))

  gdt2 <- copy(gdt)
  cdt2 <- copy(cdt)

  gdt2[, spot_id := str_trim(SpaCET_position)]
  cdt2[, spot_id := str_trim(spot_id)]

  if (any(duplicated(gdt2$spot_id))) {
    warning("[WARN] Duplicated spot_id in gsMap table: ", sid, " (keeping first)")
    gdt2 <- gdt2[!duplicated(spot_id)]
  }
  if (any(duplicated(cdt2$spot_id))) {
    warning("[WARN] Duplicated spot_id in celltype table: ", sid, " (keeping first)")
    cdt2 <- cdt2[!duplicated(spot_id)]
  }

  if ("sample_id" %in% names(cdt2)) cdt2[, sample_id := NULL]

  m_dt <- merge(
    x = gdt2,
    y = cdt2,
    by = "spot_id",
    all.x = TRUE,
    sort = FALSE
  )

  celltype_cols <- setdiff(names(cdt2), "spot_id")
  n_miss <- sum(!complete.cases(m_dt[, ..celltype_cols]))
  if (n_miss > 0) {
    warning("[WARN] Some spots missing cell-type proportions after merge: ",
            sid, " (spots with any NA celltype values = ", n_miss, ")")
  }

  m_dt[]
}

gsmap_celltype_merged_list <- setNames(vector("list", length(sample_ids)), sample_ids)

for (sid in sample_ids) {
  cat("[INFO] Merging gsMap + celltype: ", sid, "\n")
  gsmap_celltype_merged_list[[sid]] <- merge_gsmap_with_celltype(
    sid = sid,
    merged_list = merged_list,
    spacet_celltype_list = spacet_celltype_list
  )
}

n_ok <- sum(!vapply(gsmap_celltype_merged_list, is.null, logical(1)))
cat("[INFO] gsMap + celltype merged: ", n_ok, "/", length(sample_ids), "\n")

############################################################
## Step 9) Cell-type proportion and gsMap signal correlation
############################################################
celltype_panels <- list(
  Major = c(
    "Malignant",
    "CAF", "Endothelial",
    "Plasma", "B cell", "T CD4", "T CD8", "NK", "cDC", "pDC",
    "Macrophage", "Mast", "Neutrophil"
  ),
  Bcell = c(
    "B cell naive", "B cell non-switched memory",
    "B cell switched memory", "B cell exhausted"
  ),
  TCD4 = c(
    "T CD4 naive", "Th1", "Th2", "Th17", "Tfh", "Treg"
  ),
  TCD8 = c(
    "T CD8 naive", "T CD8 central memory",
    "T CD8 effector memory", "T CD8 effector", "T CD8 exhausted"
  ),
  cDC = c(
    "cDC1 CLEC9A", "cDC2 CD1C", "cDC3 LAMP3"
  ),
  Macrophage_sub = c(
    "Macrophage M1", "Macrophage M2"
  )
)

celltypes_to_test <- unique(unlist(celltype_panels))

celltype_cor_list <- list()

for (sid in sample_ids) {
  cat("[INFO] Processing: ", sid, "\n")

  dt <- gsmap_celltype_merged_list[[sid]]
  if (is.null(dt)) next

  dt[, signal := -log10(pmax(p, .Machine$double.xmin))]
  available_ct <- intersect(celltypes_to_test, names(dt))

  res_list <- lapply(available_ct, function(ct) {
    x <- dt[[ct]]
    y <- dt$signal

    if (sd(x, na.rm = TRUE) == 0 || all(is.na(x))) return(NULL)

    test <- suppressWarnings(cor.test(x, y, method = "spearman"))

    data.table(
      sample_id = sid,
      cell_type = ct,
      rho = unname(test$estimate),
      pval = test$p.value
    )
  })

  celltype_cor_list[[sid]] <- rbindlist(res_list)
}

celltype_cor_dt <- rbindlist(celltype_cor_list)
celltype_cor_dt[, FDR := p.adjust(pval, method = "BH")]

############################################################
## Step 10) Prepare cell-type visualization data
############################################################
suppressPackageStartupMessages({
  library(grid)
  library(patchwork)
})

stopifnot(all(c("sample_id", "cell_type", "rho", "pval") %in% names(celltype_cor_dt)))

dt_viz <- copy(celltype_cor_dt)
dt_viz[pval <= 0 | is.na(pval), pval := .Machine$double.xmin]
dt_viz[, FDR := p.adjust(pval, method = "BH")]

dt_viz[, major_lineage := fcase(
  cell_type == "Malignant", "Malignant",
  cell_type %in% c("CAF", "Endothelial"), cell_type,
  cell_type %in% c("Plasma", "B cell", "T CD4", "T CD8", "NK", "cDC", "pDC", "Macrophage", "Mast", "Neutrophil"), cell_type,
  cell_type %in% c("B cell naive", "B cell non-switched memory", "B cell switched memory", "B cell exhausted"), "B cell",
  cell_type %in% c("T CD4 naive", "Th1", "Th2", "Th17", "Tfh", "Treg"), "T CD4",
  cell_type %in% c("T CD8 naive", "T CD8 central memory", "T CD8 effector memory", "T CD8 effector", "T CD8 exhausted"), "T CD8",
  cell_type %in% c("cDC1 CLEC9A", "cDC2 CD1C", "cDC3 LAMP3"), "cDC",
  cell_type %in% c("Macrophage M1", "Macrophage M2"), "Macrophage",
  default = "Other"
)]

print(dt_viz[, .N, by = major_lineage][order(-N)])

dt_plot <- dt_viz[rho > 0 & FDR < 0.05]

cat("[INFO] Rows kept (rho > 0 and FDR < 0.05): ", nrow(dt_plot), "\n")
cat("[INFO] Unique cell types kept: ", uniqueN(dt_plot$cell_type), "\n")
cat("[INFO] Unique samples represented: ", uniqueN(dt_plot$sample_id), "\n\n")

keep_ct_dt <- dt_plot[, .(
  n_sig_samples = uniqueN(sample_id),
  median_rho = median(rho, na.rm = TRUE)
), by = .(cell_type, major_lineage)][order(-n_sig_samples, -median_rho)]

print(keep_ct_dt)
cat("\n[INFO] Surviving cell types by major_lineage:\n")
print(keep_ct_dt[, .N, by = major_lineage][order(-N)])

############################################################
## Step 11) Reorder cell-type blocks
############################################################
block_list <- list(
  Malignant = c("Malignant"),
  Plasma    = c("Plasma"),
  Bcell     = c("B cell", "B cell naive", "B cell non-switched memory", "B cell switched memory", "B cell exhausted"),
  TCD4      = c("T CD4", "T CD4 naive", "Th1", "Th2", "Th17", "Tfh", "Treg"),
  TCD8      = c("T CD8", "T CD8 naive", "T CD8 central memory", "T CD8 effector memory", "T CD8 effector", "T CD8 exhausted"),
  NK        = c("NK"),
  cDC       = c("cDC", "cDC1 CLEC9A", "cDC2 CD1C", "cDC3 LAMP3"),
  pDC       = c("pDC"),
  Macro     = c("Macrophage", "Macrophage M1", "Macrophage M2"),
  Mast      = c("Mast")
)

survivors <- keep_ct_dt$cell_type

order_one_block <- function(v) {
  v <- v[v %in% survivors]
  if (length(v) <= 1) return(v)
  major <- v[1]
  subs  <- setdiff(v, major)
  if (length(subs) == 0) return(major)

  subs_ord <- keep_ct_dt[cell_type %in% subs][order(-n_sig_samples, -median_rho), cell_type]
  c(major, subs_ord)
}

block_score_dt <- rbindlist(lapply(names(block_list), function(bn) {
  ct_vec <- block_list[[bn]]
  ct_vec <- ct_vec[ct_vec %in% survivors]
  if (length(ct_vec) == 0) return(NULL)

  tmp <- keep_ct_dt[cell_type %in% ct_vec]
  data.table(
    block = bn,
    max_n = max(tmp$n_sig_samples),
    max_median_rho = max(tmp$median_rho)
  )
}), fill = TRUE)

setorder(block_score_dt, -max_n, -max_median_rho)
print(block_score_dt)

block_order <- block_score_dt$block

x_levels_final <- unlist(
  lapply(block_order, function(bn) order_one_block(block_list[[bn]])),
  use.names = FALSE
)
x_levels_final <- unique(x_levels_final)

cat("[INFO] # cell types to plot: ", length(x_levels_final), "\n")
print(x_levels_final)

dt_plot[, cell_type := factor(cell_type, levels = x_levels_final)]

count_dt <- copy(keep_ct_dt[cell_type %in% x_levels_final])
count_dt[, cell_type := factor(cell_type, levels = x_levels_final)]
count_dt[, n_total_samples := uniqueN(dt_viz$sample_id)]

############################################################
## Step 12) Plot cell-type correlation results
############################################################
major_cols <- c(
  "Malignant"   = "#7A1F2B",
  "CAF"         = "#7A6A58",
  "Endothelial" = "#355C7D",
  "Plasma"      = "#8E6C8A",
  "B cell"      = "#3D6EA6",
  "T CD4"       = "#3B8C6E",
  "T CD8"       = "#B77A1A",
  "NK"          = "#6A9FB5",
  "cDC"         = "#9C4F3D",
  "pDC"         = "#7F7F7F",
  "Macrophage"  = "#5C7A3A",
  "Mast"        = "#6E4E3A",
  "Neutrophil"  = "#3F3F3F"
)

ct_to_block <- rbindlist(lapply(names(block_list), function(bn) {
  data.table(block = bn, cell_type = block_list[[bn]])
}))
ct_to_block <- ct_to_block[cell_type %in% x_levels_final]
ct_to_block[, cell_type := factor(cell_type, levels = x_levels_final)]
setorder(ct_to_block, cell_type)

rle_block <- rle(as.character(ct_to_block$block))
end_idx <- cumsum(rle_block$lengths)
vline_pos <- end_idx[-length(end_idx)] + 0.5

p_bar <- ggplot(count_dt, aes(x = cell_type, y = n_sig_samples, fill = major_lineage)) +
  geom_col(width = 0.78, colour = NA, alpha = 0.95) +
  geom_hline(
    yintercept = count_dt$n_total_samples[1],
    linewidth = 0.35, linetype = "dashed", colour = "grey35"
  ) +
  geom_vline(xintercept = vline_pos, colour = "grey92", linewidth = 0.35) +
  scale_fill_manual(values = major_cols, guide = "none") +
  scale_y_continuous(
    limits = c(0, count_dt$n_total_samples[1]),
    breaks = sort(unique(c(0, 5, 10, count_dt$n_total_samples[1]))),
    expand = expansion(mult = c(0, 0.03))
  ) +
  labs(x = NULL, y = "# samples") +
  theme_minimal(base_size = 8.0) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.28, colour = "grey88"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 8.0, margin = margin(r = 4), colour = "black"),
    axis.text.y  = element_text(size = 7.2, colour = "black"),
    plot.margin  = margin(3, 5, 0, 3)
  )

rho_sum_dt <- dt_plot[, .(
  n_sig_samples = uniqueN(sample_id),
  median_rho = median(rho, na.rm = TRUE),
  q1 = quantile(rho, 0.25, na.rm = TRUE),
  q3 = quantile(rho, 0.75, na.rm = TRUE)
), by = .(cell_type, major_lineage)]

rho_sum_dt[, cell_type := factor(cell_type, levels = levels(dt_plot$cell_type))]
setorder(rho_sum_dt, cell_type)

p_rho <- ggplot(dt_plot, aes(x = cell_type, y = rho)) +
  geom_vline(xintercept = vline_pos, colour = "grey92", linewidth = 0.35) +
  geom_linerange(
    data = rho_sum_dt,
    aes(x = cell_type, ymin = q1, ymax = q3, colour = major_lineage),
    linewidth = 1.05,
    alpha = 0.45,
    inherit.aes = FALSE
  ) +
  geom_segment(
    data = rho_sum_dt,
    aes(x = cell_type, xend = cell_type, y = median_rho, yend = median_rho, colour = major_lineage),
    linewidth = 1.1,
    lineend = "round",
    inherit.aes = FALSE
  ) +
  geom_point(
    aes(colour = major_lineage),
    position = position_jitter(width = 0.12, height = 0, seed = 1),
    size = 1,
    alpha = 0.90
  ) +
  scale_colour_manual(values = major_cols, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  labs(x = "Cell type", y = "Spearman rho") +
  theme_minimal(base_size = 8.0) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(linewidth = 0.28, colour = "grey88"),
    axis.title.x = element_text(size = 8.0, margin = margin(t = 6), colour = "black"),
    axis.title.y = element_text(size = 8.0, margin = margin(r = 6), colour = "black"),
    axis.text.x  = element_text(size = 6.9, angle = 45, hjust = 1, vjust = 1, colour = "black"),
    axis.text.y  = element_text(size = 7.2, colour = "black"),
    axis.ticks   = element_line(linewidth = 0.26, colour = "black"),
    axis.ticks.length = unit(2, "pt"),
    axis.line    = element_line(linewidth = 0.26, colour = "black"),
    plot.margin  = margin(0, 5, 3, 3)
  )

p_combined <- p_bar / p_rho + plot_layout(heights = c(1, 4))

out_pdf <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/2.cell.type/Fig_celltype_rho_with_counts.smallpanel.pdf"

ggsave(
  filename = out_pdf,
  plot = p_combined,
  width = 140,
  height = 75,
  units = "mm",
  useDingbats = FALSE
)

cat("[INFO] Saved:\n  ", out_pdf, "\n")

############################################################
## Step 13) Save relevant tables
############################################################
out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/2.cell.type"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

fwrite(
  celltype_cor_dt,
  file = file.path(out_dir, "Table_celltype_cor_dt.raw.txt"),
  sep = "\t"
)

fwrite(
  rho_sum_dt,
  file = file.path(out_dir, "Table_rho_sum_dt.sigPos_FDRlt0.05.median_IQR.txt"),
  sep = "\t"
)

cat("[INFO] Saved tables to:\n  ", out_dir, "\n")