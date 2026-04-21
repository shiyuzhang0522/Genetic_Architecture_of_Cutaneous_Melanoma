############################################################
## Immune cell colocalization and enrichment at the boundary
##
## Author: Shelley
## Last Updated:   2026-03-18
############################################################

############################################################
## Step 1) Load packages
############################################################
suppressPackageStartupMessages({
  library(data.table)
})

############################################################
## Step 2) Read full data tables
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
gsmap_dir     <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/gsMAP"
spacet_outdir <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/SpaCET.outputs"

## read full files into lists
gsmap_list <- setNames(lapply(sample_ids, function(sid) {
  fread(file.path(gsmap_dir, paste0(sid, "_gsMap_SpaCETgrid_table.txt")))
}), sample_ids)

iface_list <- setNames(lapply(sample_ids, function(sid) {
  fread(file.path(spacet_outdir, sid, "tables", paste0(sid, ".SpaCET.interface_labels.tsv")))
}), sample_ids)

prop_list <- setNames(lapply(sample_ids, function(sid) {
  fread(file.path(spacet_outdir, sid, "tables", paste0(sid, ".SpaCET.celltype_prop.tsv")))
}), sample_ids)

############################################################
## Step 3) Merge gsMap, interface labels, and cell-type
##         proportions into one table per sample
##
## Merge key:
##   gsMap$SpaCET_position  <->  interface$spot_id
##                         <->  prop$spot_id
############################################################

merge_one_sample <- function(gsmap_dt, iface_dt, prop_dt, sid = "sample") {
  dt_g <- copy(gsmap_dt)
  dt_i <- copy(iface_dt)
  dt_p <- copy(prop_dt)

  stopifnot("SpaCET_position" %in% names(dt_g))
  stopifnot(all(c("spot_id", "interface_label") %in% names(dt_i)))
  stopifnot("spot_id" %in% names(dt_p))

  setnames(dt_i, old = "spot_id", new = "SpaCET_position")
  setnames(dt_p, old = "spot_id", new = "SpaCET_position")

  merged_dt <- merge(dt_g, dt_i, by = "SpaCET_position", all.x = TRUE, sort = FALSE)
  merged_dt <- merge(merged_dt, dt_p, by = "SpaCET_position", all.x = TRUE, sort = FALSE)

  merged_dt[, sample_id := sid]
  merged_dt[]
}

merged_list <- setNames(vector("list", length(sample_ids)), sample_ids)

for (sid in sample_ids) {
  merged_list[[sid]] <- merge_one_sample(
    gsmap_dt = gsmap_list[[sid]],
    iface_dt = iface_list[[sid]],
    prop_dt  = prop_list[[sid]],
    sid      = sid
  )
}

## quick peek
sid0 <- sample_ids[1]
cat("[INFO] Merged columns (first 30) for ", sid0, ":\n", sep = "")
print(head(names(merged_list[[sid0]]), 30))
cat("[INFO] Merged rows (first 3):\n")
print(merged_list[[sid0]][1:3])

############################################################
## Step 4) Define boundary set (Interface plus 6-neighbors)
############################################################

get_hex_neighbors <- function(spot_vec) {
  if (length(spot_vec) == 0) return(character(0))

  xy <- t(vapply(strsplit(spot_vec, "x"), function(v) as.numeric(v), numeric(2)))
  sx <- xy[, 1]
  sy <- xy[, 2]

  c(
    paste0(sx,   "x", sy - 2),
    paste0(sx,   "x", sy + 2),
    paste0(sx-1, "x", sy - 1),
    paste0(sx+1, "x", sy + 1),
    paste0(sx-1, "x", sy + 1),
    paste0(sx+1, "x", sy - 1)
  )
}

add_boundary_flags_one <- function(dt,
                                   interface_label_name = "Interface",
                                   spot_col = "SpaCET_position",
                                   label_col = "interface_label") {
  stopifnot(all(c(spot_col, label_col) %in% names(dt)))

  dt <- copy(dt)
  all_spots <- dt[[spot_col]]

  iface_spots <- dt[
    !is.na(get(label_col)) & get(label_col) == interface_label_name,
    get(spot_col)
  ]

  neigh <- unique(get_hex_neighbors(iface_spots))
  neigh <- intersect(neigh, all_spots)

  boundary_set <- unique(c(iface_spots, neigh))
  boundary_set <- intersect(boundary_set, all_spots)

  dt[, is_interface := (!is.na(get(label_col))) & (get(label_col) == interface_label_name)]
  dt[, is_neighbor_of_interface := get(spot_col) %in% neigh]
  dt[, is_boundary := get(spot_col) %in% boundary_set]

  dt[]
}

merged_list2 <- merged_list
for (sid in names(merged_list2)) {
  merged_list2[[sid]] <- add_boundary_flags_one(
    dt = merged_list2[[sid]],
    interface_label_name = "Interface"
  )
}

## summarize boundary composition
boundary_comp_dt <- rbindlist(lapply(names(merged_list2), function(sid) {
  dt_all <- merged_list2[[sid]]
  dt_b <- dt_all[is_boundary == TRUE]

  if (nrow(dt_b) == 0) {
    return(data.table(
      sample_id = sid,
      n_boundary = 0L,
      n_interface = 0L,
      n_neighbor_only = 0L,
      interface_label = NA_character_,
      N = 0L,
      frac_spots = NA_real_
    ))
  }

  dt_comp <- dt_b[, .N, by = interface_label][order(-N)]
  dt_comp[, frac_spots := N / sum(N)]
  dt_comp[, sample_id := sid]

  dt_comp[, n_boundary := nrow(dt_b)]
  dt_comp[, n_interface := sum(dt_all$is_interface, na.rm = TRUE)]
  dt_comp[, n_neighbor_only := sum(dt_all$is_neighbor_of_interface & !dt_all$is_interface, na.rm = TRUE)]

  setcolorder(
    dt_comp,
    c("sample_id", "n_boundary", "n_interface", "n_neighbor_only",
      "interface_label", "N", "frac_spots")
  )
  dt_comp
}), use.names = TRUE, fill = TRUE)

print(boundary_comp_dt)

############################################################
## Step 4.1) Sanity check interface flag
############################################################

interface_check_dt <- rbindlist(lapply(names(merged_list2), function(sid) {
  dt <- merged_list2[[sid]]

  n_is_interface <- dt[, sum(is_interface, na.rm = TRUE)]
  n_label_interface <- dt[, sum(!is.na(interface_label) & interface_label == "Interface")]

  data.table(
    sample_id = sid,
    n_is_interface = n_is_interface,
    n_label_interface = n_label_interface,
    is_match = (n_is_interface == n_label_interface)
  )
}), use.names = TRUE, fill = TRUE)

print(interface_check_dt)

if (!all(interface_check_dt$is_match)) {
  warning("[WARN] Mismatch detected between is_interface and interface_label == 'Interface' in some samples.")
} else {
  cat("[INFO] Sanity check passed: for all samples, n(is_interface == TRUE) matches n(interface_label == 'Interface').\n")
}

############################################################
## Step 5) Define lineage column sets
############################################################

col_malignant <- "Malignant"

immune_major_cols <- c(
  "Plasma", "B cell", "T CD4", "T CD8", "NK", "cDC", "pDC",
  "Macrophage", "Mast", "Neutrophil"
)

stromal_cols <- c("CAF", "Endothelial")

immune_sub_cols <- c(
  "B cell naive", "B cell non-switched memory", "B cell switched memory", "B cell exhausted",
  "T CD4 naive", "Th1", "Th2", "Th17", "Tfh", "Treg",
  "T CD8 naive", "T CD8 central memory", "T CD8 effector memory", "T CD8 effector", "T CD8 exhausted",
  "cDC1 CLEC9A", "cDC2 CD1C", "cDC3 LAMP3",
  "Macrophage M1", "Macrophage M2"
)

cat("[INFO] Major immune lineages: ", length(immune_major_cols), "\n", sep = "")
cat("[INFO] Immune sublineages: ", length(immune_sub_cols), "\n", sep = "")

############################################################
## Step 6) Colocalization analysis
##   Spearman rho between malignant proportion and each
##   major immune lineage at:
##     1) exact interface
##     2) boundary region
############################################################

############################################################
## Step 6.1) Helper for robust Spearman correlation
############################################################

safe_spearman <- function(x, y, min_n = 10) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  n_ok <- length(x)

  if (n_ok < min_n) {
    return(list(rho = NA_real_, p = NA_real_, n = n_ok))
  }

  if (stats::sd(x) == 0 || stats::sd(y) == 0) {
    return(list(rho = NA_real_, p = NA_real_, n = n_ok))
  }

  ct <- suppressWarnings(
    stats::cor.test(x, y, method = "spearman", exact = FALSE)
  )

  list(
    rho = unname(ct$estimate),
    p   = ct$p.value,
    n   = n_ok
  )
}

############################################################
## Step 6.2) Sanity check required columns
############################################################

required_cols <- c("is_interface", "is_boundary", col_malignant, immune_major_cols)

col_check_dt <- rbindlist(lapply(names(merged_list2), function(sid) {
  dt <- merged_list2[[sid]]

  data.table(
    sample_id = sid,
    missing_cols = paste(setdiff(required_cols, names(dt)), collapse = "; ")
  )
}), use.names = TRUE, fill = TRUE)

col_check_dt[, has_all_required := missing_cols == ""]

print(col_check_dt)

if (!all(col_check_dt$has_all_required)) {
  stop("[ERROR] Some samples are missing required columns. Please check `col_check_dt`.")
} else {
  cat("[INFO] All required columns are present in every sample.\n")
}

############################################################
## Step 6.3) Define region settings
############################################################

region_def_dt <- data.table(
  region = c("interface", "boundary"),
  flag_col = c("is_interface", "is_boundary")
)

print(region_def_dt)

############################################################
## Step 6.4) Compute per-sample colocalization results
############################################################

coloc_major_dt <- rbindlist(lapply(names(merged_list2), function(sid) {
  dt <- merged_list2[[sid]]

  rbindlist(lapply(seq_len(nrow(region_def_dt)), function(i) {
    region_i <- region_def_dt$region[i]
    flag_col_i <- region_def_dt$flag_col[i]

    dt_sub <- dt[get(flag_col_i) == TRUE]

    rbindlist(lapply(immune_major_cols, function(ic) {
      out <- safe_spearman(
        x = dt_sub[[col_malignant]],
        y = dt_sub[[ic]],
        min_n = 10
      )

      data.table(
        sample_id = sid,
        region = region_i,
        immune_type = ic,
        rho = out$rho,
        p = out$p,
        n = out$n
      )
    }), use.names = TRUE, fill = TRUE)

  }), use.names = TRUE, fill = TRUE)

}), use.names = TRUE, fill = TRUE)

############################################################
## Step 6.5) Quick summaries
############################################################

cat("[INFO] Colocalization result table dimensions: ",
    nrow(coloc_major_dt), " rows x ", ncol(coloc_major_dt), " columns\n", sep = "")

cat("[INFO] First few rows of colocalization results:\n")
print(coloc_major_dt[1:10])

coloc_n_dt <- unique(coloc_major_dt[, .(sample_id, region, n)])
cat("[INFO] Number of spots used per sample and region:\n")
print(coloc_n_dt[order(region, sample_id)])

cat("[INFO] Top raw positive correlations by rho:\n")
print(coloc_major_dt[!is.na(rho)][order(-rho, p)][1:20])

############################################################
## Step 7) Multiple-testing adjustment
##   1) BH within each spatial region
##   2) global BH across all tests
############################################################

############################################################
## Step 7.1) Sanity check input table
############################################################

stopifnot(all(c("sample_id", "region", "immune_type", "rho", "p", "n") %in% names(coloc_major_dt)))

cat("[INFO] Number of total colocalization tests: ", nrow(coloc_major_dt), "\n", sep = "")
cat("[INFO] Regions present: ", paste(unique(coloc_major_dt$region), collapse = ", "), "\n", sep = "")

############################################################
## Step 7.2) BH within each region
############################################################

coloc_major_dt[, p_bh_region := stats::p.adjust(p, method = "BH"), by = region]

############################################################
## Step 7.3) Global BH across all tests
############################################################

coloc_major_dt[, p_bh_global := stats::p.adjust(p, method = "BH")]

############################################################
## Step 7.4) Inspect top signals
############################################################

cat("[INFO] Top associations ranked by within-region BH-adjusted p-value:\n")
print(
  coloc_major_dt[order(p_bh_region)][
    1:20,
    .(sample_id, region, immune_type, rho, p, p_bh_region, p_bh_global, n)
  ]
)

cat("[INFO] Top associations ranked by global BH-adjusted p-value:\n")
print(
  coloc_major_dt[order(p_bh_global)][
    1:20,
    .(sample_id, region, immune_type, rho, p, p_bh_region, p_bh_global, n)
  ]
)

############################################################
## Step 8) Extract significant positive colocalization
##         results based on region-level FDR
############################################################

############################################################
## Step 8.1) Define output directory
############################################################

out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsMap_0318/spatial.context.heterogeneity"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

############################################################
## Step 8.2) Subset significant positive results
############################################################

coloc_sigpos_region_dt <- coloc_major_dt[
  !is.na(rho) &
    !is.na(p_bh_region) &
    rho > 0 &
    p_bh_region < 0.05
]

cat("[INFO] Total significant positive colocalization results (region-level FDR < 0.05): ",
    nrow(coloc_sigpos_region_dt), "\n", sep = "")

coloc_sigpos_interface_dt <- coloc_sigpos_region_dt[region == "interface"]
coloc_sigpos_boundary_dt  <- coloc_sigpos_region_dt[region == "boundary"]

############################################################
## Step 8.3) Print region-specific summaries
############################################################

summarize_region_sigpos <- function(dt, region_name) {
  cat("\n============================================================\n")
  cat("[INFO] Region: ", region_name, "\n", sep = "")
  cat("[INFO] Number of significant positive results: ", nrow(dt), "\n", sep = "")

  if (nrow(dt) == 0) {
    cat("[INFO] No significant positive colocalization results detected.\n")
    return(invisible(NULL))
  }

  cat("[INFO] Number of unique samples: ", uniqueN(dt$sample_id), "\n", sep = "")
  cat("[INFO] Number of unique immune cell types: ", uniqueN(dt$immune_type), "\n", sep = "")

  immune_summary_dt <- dt[
    ,
    .(
      n_sig_samples = uniqueN(sample_id),
      mean_rho = mean(rho, na.rm = TRUE),
      median_rho = median(rho, na.rm = TRUE),
      min_p_bh_region = min(p_bh_region, na.rm = TRUE)
    ),
    by = immune_type
  ][order(-n_sig_samples, -median_rho)]

  cat("[INFO] Summary by immune cell type:\n")
  print(immune_summary_dt)

  cat("[INFO] Top individual signals:\n")
  print(
    dt[order(p_bh_region, -rho)][
      ,
      .(sample_id, immune_type, rho, p, p_bh_region, p_bh_global, n)
    ][1:min(20, .N)]
  )

  invisible(immune_summary_dt)
}

interface_summary_dt <- summarize_region_sigpos(
  dt = coloc_sigpos_interface_dt,
  region_name = "interface"
)

boundary_summary_dt <- summarize_region_sigpos(
  dt = coloc_sigpos_boundary_dt,
  region_name = "boundary"
)

############################################################
## Step 8.4) Save TSV files
############################################################

out_interface <- file.path(
  out_dir,
  "Table_coloc_majorImmune_sigpos.regionFDR0.05.interface.tsv"
)

out_boundary <- file.path(
  out_dir,
  "Table_coloc_majorImmune_sigpos.regionFDR0.05.boundary.tsv"
)

fwrite(
  coloc_sigpos_interface_dt,
  file = out_interface,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

fwrite(
  coloc_sigpos_boundary_dt,
  file = out_boundary,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("\n[INFO] Saved TSV files:\n")
cat("  ", out_interface, "\n", sep = "")
cat("  ", out_boundary, "\n", sep = "")

############################################################
## Step 9) Estimate overall T-cell proportion per spot
##   T_cell_overall = T CD4 + T CD8
############################################################

############################################################
## Step 9.1) Sanity check required columns
############################################################

required_cols_step9 <- c("T CD4", "T CD8")

col_check_step9_dt <- rbindlist(lapply(names(merged_list2), function(sid) {
  dt <- merged_list2[[sid]]

  data.table(
    sample_id = sid,
    missing_cols = paste(setdiff(required_cols_step9, names(dt)), collapse = "; ")
  )
}), use.names = TRUE, fill = TRUE)

col_check_step9_dt[, has_required := missing_cols == ""]

print(col_check_step9_dt)

if (!all(col_check_step9_dt$has_required)) {
  stop("[ERROR] Missing T-cell columns in some samples.")
} else {
  cat("[INFO] All samples contain T CD4 and T CD8 columns.\n")
}

############################################################
## Step 9.2) Add overall T-cell proportion
############################################################

merged_list3 <- merged_list2

for (sid in names(merged_list3)) {
  dt <- merged_list3[[sid]]

  dt[, T_cell_overall := rowSums(
    .SD,
    na.rm = TRUE
  ), .SDcols = c("T CD4", "T CD8")]

  merged_list3[[sid]] <- dt
}

############################################################
## Step 9.3) Quick sanity check
############################################################

sid0 <- names(merged_list3)[1]

cat("[INFO] Preview of T_cell_overall for sample: ", sid0, "\n", sep = "")
print(
  merged_list3[[sid0]][
    ,
    .(
      min_T = min(T_cell_overall, na.rm = TRUE),
      max_T = max(T_cell_overall, na.rm = TRUE),
      mean_T = mean(T_cell_overall, na.rm = TRUE)
    )
  ]
)

check_overflow_dt <- rbindlist(lapply(names(merged_list3), function(sid) {
  dt <- merged_list3[[sid]]

  data.table(
    sample_id = sid,
    n_over_1 = sum(dt$T_cell_overall > 1, na.rm = TRUE)
  )
}))

cat("[INFO] Spots with T_cell_overall > 1 (sanity check):\n")
print(check_overflow_dt)

############################################################
## Step 9.4) Optional summary by region
############################################################

tcell_summary_dt <- rbindlist(lapply(names(merged_list3), function(sid) {
  dt <- merged_list3[[sid]]

  rbindlist(list(
    dt[is_interface == TRUE,
       .(sample_id = sid,
         region = "interface",
         mean_T = mean(T_cell_overall, na.rm = TRUE),
         median_T = median(T_cell_overall, na.rm = TRUE),
         n = .N)],

    dt[is_boundary == TRUE,
       .(sample_id = sid,
         region = "boundary",
         mean_T = mean(T_cell_overall, na.rm = TRUE),
         median_T = median(T_cell_overall, na.rm = TRUE),
         n = .N)]
  ))

}), use.names = TRUE, fill = TRUE)

cat("[INFO] T-cell abundance summary by region:\n")
print(tcell_summary_dt)

############################################################
## Step 9.5) Summarize T-cell distribution by spatial region
############################################################

############################################################
## Step 9.5.1) Helper function
############################################################

summarize_T <- function(x) {
  x <- x[is.finite(x)]

  if (length(x) == 0) {
    return(list(
      n_spots = 0L,
      mean_T = NA_real_,
      median_T = NA_real_,
      q1_T = NA_real_,
      q3_T = NA_real_,
      IQR_T = NA_real_,
      p90_T = NA_real_
    ))
  }

  list(
    n_spots = length(x),
    mean_T = mean(x),
    median_T = median(x),
    q1_T = as.numeric(quantile(x, 0.25)),
    q3_T = as.numeric(quantile(x, 0.75)),
    IQR_T = IQR(x),
    p90_T = as.numeric(quantile(x, 0.90))
  )
}

############################################################
## Step 9.5.2) Compute summaries
############################################################

tcell_region_summary_dt <- rbindlist(lapply(names(merged_list3), function(sid) {
  dt <- merged_list3[[sid]]

  dt_interface <- dt[is_interface == TRUE]
  dt_boundary  <- dt[is_boundary == TRUE]
  dt_tumor     <- dt[interface_label == "Tumor"]
  dt_stroma    <- dt[interface_label == "Stroma"]

  rbindlist(list(
    as.data.table(c(list(sample_id = sid, region = "interface"),
                    summarize_T(dt_interface$T_cell_overall))),
    as.data.table(c(list(sample_id = sid, region = "boundary"),
                    summarize_T(dt_boundary$T_cell_overall))),
    as.data.table(c(list(sample_id = sid, region = "tumor"),
                    summarize_T(dt_tumor$T_cell_overall))),
    as.data.table(c(list(sample_id = sid, region = "stroma"),
                    summarize_T(dt_stroma$T_cell_overall)))
  ), use.names = TRUE, fill = TRUE)

}), use.names = TRUE, fill = TRUE)

############################################################
## Step 9.5.3) Output
############################################################

cat("[INFO] T-cell distribution summary across regions:\n")
print(tcell_region_summary_dt)

############################################################
## Step 9.5.4) Save T-cell region summary
############################################################

out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsMap_0318/spatial.context.heterogeneity"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_file <- file.path(
  out_dir,
  "Table_Tcell_distribution_by_region_per_sample.tsv"
)

fwrite(
  tcell_region_summary_dt,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("[INFO] Saved:\n  ", out_file, "\n", sep = "")

############################################################
## Step 10) Estimate per-sample T-cell enrichment at the
##          interface and boundary using bootstrap
############################################################

############################################################
## Step 10.1) Settings
############################################################

eps <- 1e-4
B   <- 2000

############################################################
## Step 10.2) Helper to compute observed enrichment ratio
############################################################

calc_region_enrichment <- function(dt, target_flag, region_name, eps = 1e-4) {
  stopifnot(all(c(target_flag, "T_cell_overall") %in% names(dt)))

  x_target <- dt[get(target_flag) == TRUE,  T_cell_overall]
  x_bg     <- dt[get(target_flag) == FALSE, T_cell_overall]

  x_target <- x_target[is.finite(x_target)]
  x_bg     <- x_bg[is.finite(x_bg)]

  n_target <- length(x_target)
  n_bg     <- length(x_bg)

  med_target <- if (n_target > 0) median(x_target) else NA_real_
  med_bg     <- if (n_bg > 0) median(x_bg) else NA_real_

  mean_target <- if (n_target > 0) mean(x_target) else NA_real_
  mean_bg     <- if (n_bg > 0) mean(x_bg) else NA_real_

  q1_target <- if (n_target > 0) as.numeric(quantile(x_target, 0.25)) else NA_real_
  q3_target <- if (n_target > 0) as.numeric(quantile(x_target, 0.75)) else NA_real_
  q1_bg     <- if (n_bg > 0) as.numeric(quantile(x_bg, 0.25)) else NA_real_
  q3_bg     <- if (n_bg > 0) as.numeric(quantile(x_bg, 0.75)) else NA_real_

  ratio_median <- if (n_target > 0 && n_bg > 0) {
    (med_target + eps) / (med_bg + eps)
  } else {
    NA_real_
  }

  log_ratio_median <- if (is.finite(ratio_median)) log(ratio_median) else NA_real_

  data.table(
    region = region_name,
    n_target = n_target,
    n_background = n_bg,
    median_target = med_target,
    median_background = med_bg,
    mean_target = mean_target,
    mean_background = mean_bg,
    q1_target = q1_target,
    q3_target = q3_target,
    q1_background = q1_bg,
    q3_background = q3_bg,
    ratio_median = ratio_median,
    log_ratio_median = log_ratio_median
  )
}

############################################################
## Step 10.3) Helper to bootstrap enrichment ratio
############################################################

bootstrap_ratio <- function(x_target, x_bg, B = 2000, eps = 1e-4, min_n = 5) {
  x_target <- x_target[is.finite(x_target)]
  x_bg     <- x_bg[is.finite(x_bg)]

  n_t <- length(x_target)
  n_b <- length(x_bg)

  if (n_t < min_n || n_b < min_n) {
    return(list(
      boot_median = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      prop_gt1 = NA_real_,
      p_boot_two_sided = NA_real_
    ))
  }

  ratios <- numeric(B)

  for (i in seq_len(B)) {
    bt_t <- sample(x_target, size = n_t, replace = TRUE)
    bt_b <- sample(x_bg, size = n_b, replace = TRUE)

    med_t <- median(bt_t)
    med_b <- median(bt_b)

    ratios[i] <- (med_t + eps) / (med_b + eps)
  }

  log_ratios <- log(ratios)

  p_left  <- (sum(log_ratios <= 0) + 1) / (B + 1)
  p_right <- (sum(log_ratios >= 0) + 1) / (B + 1)
  p_two_sided <- min(1, 2 * min(p_left, p_right))

  list(
    boot_median = median(ratios),
    ci_low = as.numeric(quantile(ratios, 0.025)),
    ci_high = as.numeric(quantile(ratios, 0.975)),
    prop_gt1 = mean(ratios > 1),
    p_boot_two_sided = p_two_sided
  )
}

############################################################
## Step 10.4) Compute observed enrichment and bootstrap CI
############################################################

tcell_enrich_boot_dt <- rbindlist(lapply(names(merged_list3), function(sid) {
  dt <- merged_list3[[sid]]

  ## interface vs non-interface
  obs_if <- calc_region_enrichment(
    dt = dt,
    target_flag = "is_interface",
    region_name = "interface",
    eps = eps
  )

  x_t_if <- dt[is_interface == TRUE,  T_cell_overall]
  x_b_if <- dt[is_interface == FALSE, T_cell_overall]

  boot_if <- bootstrap_ratio(
    x_target = x_t_if,
    x_bg = x_b_if,
    B = B,
    eps = eps
  )

  obs_if[, `:=`(
    ratio_boot_median = boot_if$boot_median,
    ci_low = boot_if$ci_low,
    ci_high = boot_if$ci_high,
    prop_gt1 = boot_if$prop_gt1,
    p_boot_two_sided = boot_if$p_boot_two_sided
  )]

  ## boundary vs non-boundary
  obs_bd <- calc_region_enrichment(
    dt = dt,
    target_flag = "is_boundary",
    region_name = "boundary",
    eps = eps
  )

  x_t_bd <- dt[is_boundary == TRUE,  T_cell_overall]
  x_b_bd <- dt[is_boundary == FALSE, T_cell_overall]

  boot_bd <- bootstrap_ratio(
    x_target = x_t_bd,
    x_bg = x_b_bd,
    B = B,
    eps = eps
  )

  obs_bd[, `:=`(
    ratio_boot_median = boot_bd$boot_median,
    ci_low = boot_bd$ci_low,
    ci_high = boot_bd$ci_high,
    prop_gt1 = boot_bd$prop_gt1,
    p_boot_two_sided = boot_bd$p_boot_two_sided
  )]

  out <- rbindlist(list(obs_if, obs_bd), use.names = TRUE, fill = TRUE)
  out[, sample_id := sid]
  setcolorder(out, c("sample_id", "region"))
  out

}), use.names = TRUE, fill = TRUE)

############################################################
## Step 10.5) Quick summaries
############################################################

cat("[INFO] Per-sample T-cell enrichment with bootstrap CI and two-sided P value:\n")
print(tcell_enrich_boot_dt)

tcell_enrich_summary_dt <- tcell_enrich_boot_dt[
  ,
  .(
    n_samples = sum(is.finite(ratio_median)),
    n_ratio_gt_1 = sum(ratio_median > 1, na.rm = TRUE),
    n_ci_excludes_1 = sum(
      is.finite(ci_low) & is.finite(ci_high) &
        (ci_low > 1 | ci_high < 1),
      na.rm = TRUE
    ),
    n_p_boot_lt_0.05 = sum(
      is.finite(p_boot_two_sided) & p_boot_two_sided < 0.05,
      na.rm = TRUE
    ),
    mean_ratio = mean(ratio_median, na.rm = TRUE),
    median_ratio = median(ratio_median, na.rm = TRUE),
    q1_ratio = as.numeric(quantile(ratio_median, 0.25, na.rm = TRUE)),
    q3_ratio = as.numeric(quantile(ratio_median, 0.75, na.rm = TRUE)),
    mean_prop_gt1 = mean(prop_gt1, na.rm = TRUE),
    median_prop_gt1 = median(prop_gt1, na.rm = TRUE)
  ),
  by = region
]

cat("[INFO] Cross-sample summary of T-cell enrichment:\n")
print(tcell_enrich_summary_dt)

############################################################
## Step 10.6) Wilcoxon signed-rank test on observed
##            log-ratio across samples
############################################################

wilcox_region <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 3) {
    return(list(n = length(x), p = NA_real_))
  }

  wt <- suppressWarnings(
    wilcox.test(x = x, mu = 0, alternative = "greater", exact = FALSE)
  )

  list(n = length(x), p = wt$p.value)
}

tcell_enrich_test_dt <- rbindlist(lapply(unique(tcell_enrich_boot_dt$region), function(reg) {
  x <- tcell_enrich_boot_dt[region == reg, log_ratio_median]
  out <- wilcox_region(x)

  data.table(
    region = reg,
    n_samples_tested = out$n,
    p_wilcox_greater_0 = out$p
  )
}), use.names = TRUE, fill = TRUE)

cat("[INFO] Wilcoxon signed-rank test on observed log-ratio:\n")
print(tcell_enrich_test_dt)

############################################################
## Step 10.7) Save outputs
############################################################

out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsMap_0318/spatial.context.heterogeneity"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_file_1 <- file.path(
  out_dir,
  "Table_Tcell_enrichment_ratio.bootstrapCI.by_sample.tsv"
)

out_file_2 <- file.path(
  out_dir,
  "Table_Tcell_enrichment_ratio.bootstrapCI.summary.tsv"
)

out_file_3 <- file.path(
  out_dir,
  "Table_Tcell_enrichment_ratio.wilcox.tsv"
)

fwrite(tcell_enrich_boot_dt, out_file_1, sep = "\t", quote = FALSE, na = "NA")
fwrite(tcell_enrich_summary_dt, out_file_2, sep = "\t", quote = FALSE, na = "NA")
fwrite(tcell_enrich_test_dt, out_file_3, sep = "\t", quote = FALSE, na = "NA")

cat("[INFO] Saved:\n")
cat("  ", out_file_1, "\n", sep = "")
cat("  ", out_file_2, "\n", sep = "")
cat("  ", out_file_3, "\n", sep = "")

############################################################
## Step 11.1) Read gsMap-distance correlation summary
############################################################

corr_file <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsmap_ST_analysis/6.spatial.heterogeneity/Tumor_dist_vs_gsMap.corr_summary.tsv"

if (!file.exists(corr_file)) {
  stop("[ERROR] File not found: ", corr_file)
}

corr_dt <- fread(corr_file, sep = "\t")

cat("[INFO] Loaded gsMap-distance correlation summary:\n")
print(corr_dt)

cat("[INFO] Columns in corr_dt:\n")
print(names(corr_dt))

if (!"sample_id" %in% names(corr_dt)) {
  stop("[ERROR] `sample_id` column is missing in corr_dt.")
}

############################################################
## Step 11.2) Subset boundary T-cell enrichment results
############################################################

tcell_boundary_dt <- tcell_enrich_boot_dt[
  region == "boundary",
  .(
    sample_id,
    region,
    n_target,
    n_background,
    median_target,
    median_background,
    mean_target,
    mean_background,
    ratio_median,
    log_ratio_median,
    ci_low,
    ci_high,
    p_boot_two_sided
  )
]

cat("[INFO] Boundary-only T-cell enrichment table:\n")
print(tcell_boundary_dt)

############################################################
## Step 11.3) Merge gsMap-distance pattern with boundary
##            T-cell enrichment by sample_id
############################################################

step11_merge_dt <- merge(
  x = corr_dt,
  y = tcell_boundary_dt,
  by = "sample_id",
  all = FALSE,
  sort = FALSE
)

cat("[INFO] Merged Step 11 table:\n")
print(step11_merge_dt)

cat("[INFO] Dimensions of merged table: ",
    nrow(step11_merge_dt), " rows x ", ncol(step11_merge_dt), " columns\n", sep = "")

############################################################
## Step 11.4) Correlate gsMap-distance pattern with
##            boundary T-cell abundance
############################################################

dt_test <- step11_merge_dt[
  is.finite(spearman_rho) & is.finite(median_target)
]

cat("[INFO] Number of samples used: ", nrow(dt_test), "\n", sep = "")

cor_test <- cor.test(
  x = dt_test$spearman_rho,
  y = dt_test$median_target,
  method = "spearman",
  exact = FALSE
)

cat("[INFO] Spearman correlation results:\n")
print(cor_test)

############################################################
## Step 11.4B) Correlate gsMap-distance pattern with
##             boundary T-cell log-enrichment
############################################################

dt_test_ratio <- step11_merge_dt[
  is.finite(spearman_rho) & is.finite(log_ratio_median)
]

cat("[INFO] Number of samples used: ", nrow(dt_test_ratio), "\n", sep = "")

cor_test_logratio <- cor.test(
  x = dt_test_ratio$spearman_rho,
  y = dt_test_ratio$log_ratio_median,
  method = "spearman",
  exact = FALSE
)

cat("[INFO] Spearman correlation results (gsMap-distance rho vs boundary T-cell log-enrichment):\n")
print(cor_test_logratio)

############################################################
## Step 11.5) Visualize the association between gsMap-distance
##            pattern and boundary T-cell enrichment
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(grid)
})

############################################################
## Step 11.5.1) Prepare plotting table
############################################################

plot_dt <- copy(step11_merge_dt)

plot_dt <- plot_dt[
  is.finite(spearman_rho) &
    is.finite(log_ratio_median) &
    is.finite(median_target)
]

plot_dt[, median_target_pct := median_target * 100]
plot_dt[, color_value := log10(median_target_pct + 0.05)]

cat("[INFO] Plotting samples: ", nrow(plot_dt), "\n", sep = "")

############################################################
## Step 11.5.2) Correlation annotation text
############################################################

cor_test_logratio <- suppressWarnings(
  cor.test(
    x = plot_dt$spearman_rho,
    y = plot_dt$log_ratio_median,
    method = "spearman",
    exact = FALSE
  )
)

rho_lab <- sprintf("%.2f", unname(cor_test_logratio$estimate))
p_lab <- if (cor_test_logratio$p.value < 0.001) {
  "< 0.001"
} else {
  sprintf("= %.3f", cor_test_logratio$p.value)
}

annot_text <- paste0("Spearman \u03C1 = ", rho_lab, "\n", "P ", p_lab)

############################################################
## Step 11.5.3) Optional sample labels
############################################################

label_samples <- c(
  "MEL75",
  "MEL201"
)

plot_dt[, label := ifelse(sample_id %in% label_samples, sample_id, NA_character_)]

plot_dt[, label_short := fifelse(
  sample_id == "MEL201_CA", "MEL201",
  fifelse(sample_id == "MEL75_CA", "MEL75",
  fifelse(sample_id == "MEL239_CA", "MEL239",
  fifelse(sample_id == "MEL256_CA", "MEL256", NA_character_)))
)]

############################################################
## Step 11.5.4) Plot
############################################################

col_pal <- c("#f7dddd", "#e5a9a9", "#C97373", "#8F2D2D")

p_step11 <- ggplot(
  plot_dt,
  aes(x = spearman_rho, y = log_ratio_median, colour = color_value)
) +
  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    linewidth = 0.5,
    colour = "black",
    fill = "grey85",
    alpha = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    colour = "grey65"
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    colour = "grey65"
  ) +
  geom_point(
    size = 2.4,
    alpha = 0.95
  ) +
  geom_text_repel(
    data = plot_dt[!is.na(label)],
    aes(label = label_short),
    size = 2.1,
    min.segment.length = 0,
    segment.size = 0.25,
    box.padding = 0.18,
    point.padding = 0.10,
    seed = 1,
    colour = "black"
  ) +
  annotate(
    "text",
    x = Inf, y = Inf,
    label = annot_text,
    hjust = 1.02, vjust = 1.15,
    size = 2.5,
    colour = "black"
  ) +
  scale_colour_gradientn(
    colours = col_pal,
    name = "Boundary T-cell\nproportion (%)",
    breaks = log10(c(0.1, 0.5, 1, 2, 4) + 0.05),
    labels = c("0.1", "0.5", "1", "2", "4"),
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barheight = unit(18, "mm"),
      barwidth = unit(2.0, "mm"),
      frame.colour = NA,
      ticks.colour = "white"
    )
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0.03, 0.03))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.04, 0.06))
  ) +
  labs(
    x = "gsMap-distance Spearman \u03C1",
    y = "Boundary T-cell enrichment (log ratio)"
  ) +
  theme_classic(base_size = 7) +
  theme(
    axis.title = element_text(size = 7, colour = "black"),
    axis.text = element_text(size = 6, colour = "black"),
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks.length = unit(2, "pt"),
    legend.position = "right",
    legend.title = element_text(size = 5.5, colour = "black", lineheight = 0.9),
    legend.text = element_text(size = 5, colour = "black"),
    legend.key.height = unit(14, "mm"),
    legend.key.width  = unit(2, "mm"),
    plot.margin = margin(3, 4, 2, 2)
  )

print(p_step11)

############################################################
## Step 11.5.5) Save as PDF
############################################################

out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/ST/Updated_gsMap_0318/spatial.context.heterogeneity"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_pdf <- file.path(
  out_dir,
  "Fig_step11_gsMapDist_vs_boundaryTcellEnrichment.smallpanel.pdf"
)

ggsave(
  filename = out_pdf,
  plot = p_step11,
  width = 80,
  height = 54,
  units = "mm"
)

cat("[INFO] Saved:\n  ", out_pdf, "\n", sep = "")