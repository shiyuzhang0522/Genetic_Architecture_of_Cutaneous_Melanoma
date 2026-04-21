############################################################
## Add metadata and normalize evidence weights
##
## Purpose:
##   - Read linked PPI and TR edge tables for the curated melanoma network
##   - Inspect and resolve duplicated edges within each layer
##   - Collapse reciprocal duplicates to undirected layer-specific edges
##   - Combine PPI and TR layers into one network-ready edge table
##   - Derive robust evidence-based edge weights from n_references
##
## Notes:
##   - De-duplication is performed within each interaction layer first
##   - PPI and TR are then combined without collapsing cross-layer overlaps
##   - Therefore, the same gene pair may appear more than once if supported
##     by different interaction layers
##
## Inputs:
##   - edges_ppi_linked_dt.txt
##   - edges_tr_linked_dt.txt
##   - all_genes_long_pc_dt.txt
##
## Outputs:
##   - layer-specific deduplicated edge tables
##   - combined edge table with processed weights
##
## Author: Shelley
## Date:   2026-01-21
## Updated: 2026-03-06
############################################################

############################################################
## Step 1) Load packages
############################################################
suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

cat("[INFO] Packages loaded: data.table, igraph, ggplot2, patchwork, scales\n")

############################################################
## Step 2) Read input tables
############################################################

network_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/Updated.Omnigenic.Architecture/Network"
vertex_file <- "/Users/shelleyz/Projects/Melanoma-WGS/Updated.Omnigenic.Architecture/Vertices/all_genes_long_pc_dt.txt"

ppi_file <- file.path(network_dir, "edges_ppi_linked_dt.txt")
tr_file  <- file.path(network_dir, "edges_tr_linked_dt.txt")

stopifnot(file.exists(ppi_file))
stopifnot(file.exists(tr_file))
stopifnot(file.exists(vertex_file))

edges_ppi_dt <- fread(ppi_file)
edges_tr_dt  <- fread(tr_file)
vertices_dt  <- fread(vertex_file)

cat("[INFO] Loaded PPI edges:", ppi_file, "\n")
cat("[INFO] edges_ppi_dt:", nrow(edges_ppi_dt), "rows x", ncol(edges_ppi_dt), "cols\n")
print(head(edges_ppi_dt, 3))

cat("[INFO] Loaded TR edges:", tr_file, "\n")
cat("[INFO] edges_tr_dt:", nrow(edges_tr_dt), "rows x", ncol(edges_tr_dt), "cols\n")
print(head(edges_tr_dt, 3))

cat("[INFO] Loaded vertices:", vertex_file, "\n")
cat("[INFO] vertices_dt:", nrow(vertices_dt), "rows x", ncol(vertices_dt), "cols\n")
print(head(vertices_dt, 3))

stopifnot(nrow(edges_ppi_dt) > 0)
stopifnot(nrow(edges_tr_dt) > 0)
stopifnot(nrow(vertices_dt) > 0)

############################################################
## Step 3) Sanity check duplicated directed edges
############################################################

required_edge_cols <- c(
  "from_ensg_filled", "to_ensg_filled",
  "from_ensg", "to_ensg",
  "from_symbol", "to_symbol",
  "n_references", "curation_effort", "n_resources",
  "consensus_direction", "consensus_stimulation", "consensus_inhibition",
  "type", "interaction_type",
  "from_ensg_missing", "to_ensg_missing",
  "from_match_method", "to_match_method"
)

stopifnot(all(required_edge_cols %in% names(edges_ppi_dt)))
stopifnot(all(required_edge_cols %in% names(edges_tr_dt)))

cat("\n[PPI] unique directed edges:\n")
ppi_unique_edges <- unique(edges_ppi_dt[, .(from_ensg_filled, to_ensg_filled)])
cat("[INFO] n_unique_edges:", nrow(ppi_unique_edges), "\n")

cat("\n[TR] unique directed edges:\n")
tr_unique_edges <- unique(edges_tr_dt[, .(from_ensg_filled, to_ensg_filled)])
cat("[INFO] n_unique_edges:", nrow(tr_unique_edges), "\n")

cat("\n[PPI] duplicated directed rows:",
    sum(duplicated(edges_ppi_dt[, .(from_ensg_filled, to_ensg_filled)])), "\n")
cat("[TR] duplicated directed rows:",
    sum(duplicated(edges_tr_dt[, .(from_ensg_filled, to_ensg_filled)])), "\n")

ppi_dup_dt <- edges_ppi_dt[
  duplicated(edges_ppi_dt[, .(from_ensg_filled, to_ensg_filled)]) |
    duplicated(edges_ppi_dt[, .(from_ensg_filled, to_ensg_filled)], fromLast = TRUE)
][order(from_ensg_filled, to_ensg_filled)]

tr_dup_dt <- edges_tr_dt[
  duplicated(edges_tr_dt[, .(from_ensg_filled, to_ensg_filled)]) |
    duplicated(edges_tr_dt[, .(from_ensg_filled, to_ensg_filled)], fromLast = TRUE)
][order(from_ensg_filled, to_ensg_filled)]

cat("\n[PPI] duplicated edge rows subset:", nrow(ppi_dup_dt), "\n")
if (nrow(ppi_dup_dt) > 0) {
  print(ppi_dup_dt[, .(
    from_ensg_filled, to_ensg_filled,
    from_symbol, to_symbol,
    interaction_type, type,
    n_references, n_resources
  )][1:min(50, .N)])
}

cat("\n[TR] duplicated edge rows subset:", nrow(tr_dup_dt), "\n")
if (nrow(tr_dup_dt) > 0) {
  print(tr_dup_dt[, .(
    from_ensg_filled, to_ensg_filled,
    from_symbol, to_symbol,
    interaction_type, type,
    n_references, n_resources
  )][1:min(50, .N)])
}

cat("[INFO] Duplicate inspection complete.\n")

############################################################
## Step 4) Build ENSG -> symbol dictionary
############################################################

stopifnot(all(c("ensg", "symbol") %in% names(vertices_dt)))

ensg2symbol_dt <- unique(
  vertices_dt[
    !is.na(ensg) & ensg != "" &
      !is.na(symbol) & symbol != "",
    .(ensg, symbol)
  ]
)

ensg_symbol_multi_dt <- ensg2symbol_dt[, .N, by = ensg][N > 1]

cat("[INFO] ENSG with multiple symbols in dictionary:", nrow(ensg_symbol_multi_dt), "\n")
if (nrow(ensg_symbol_multi_dt) > 0) {
  print(ensg2symbol_dt[ensg %in% ensg_symbol_multi_dt$ensg][order(ensg, symbol)])
  stop("[ERROR] ENSG-to-symbol dictionary is not one-to-one.")
}

setkey(ensg2symbol_dt, ensg)

############################################################
## Step 5) Helper functions for edge de-duplication
############################################################

dedup_directed_edges <- function(dt) {
  dt2 <- copy(dt)

  dt2[, n_references := as.numeric(n_references)]
  dt2[, curation_effort := as.numeric(curation_effort)]
  dt2[, n_resources := as.numeric(n_resources)]

  setorder(
    dt2,
    from_ensg_filled, to_ensg_filled,
    -n_references, -curation_effort, -n_resources
  )

  dt2 <- dt2[, .SD[1], by = .(from_ensg_filled, to_ensg_filled)]
  dt2[]
}

dedup_undirected_edges <- function(dt, ensg2symbol_dt) {
  dt2 <- copy(dt)

  dt2[, n_references := as.numeric(n_references)]
  dt2[, curation_effort := as.numeric(curation_effort)]
  dt2[, n_resources := as.numeric(n_resources)]
  dt2[, consensus_direction := as.numeric(consensus_direction)]
  dt2[, consensus_stimulation := as.numeric(consensus_stimulation)]
  dt2[, consensus_inhibition := as.numeric(consensus_inhibition)]

  dt2[, gene1 := pmin(from_ensg_filled, to_ensg_filled)]
  dt2[, gene2 := pmax(from_ensg_filled, to_ensg_filled)]

  setorder(dt2, gene1, gene2, -n_references, -curation_effort, -n_resources)

  out_dt <- dt2[
    ,
    .(
      from_ensg_filled = first(from_ensg_filled),
      to_ensg_filled   = first(to_ensg_filled),
      from_ensg        = first(from_ensg),
      to_ensg          = first(to_ensg),
      from_symbol      = first(from_symbol),
      to_symbol        = first(to_symbol),

      interaction_type = first(interaction_type),
      type             = first(type),

      n_references    = max(n_references, na.rm = TRUE),
      n_resources     = max(n_resources, na.rm = TRUE),
      curation_effort = max(curation_effort, na.rm = TRUE),

      consensus_direction   = max(consensus_direction, na.rm = TRUE),
      consensus_stimulation = max(consensus_stimulation, na.rm = TRUE),
      consensus_inhibition  = max(consensus_inhibition, na.rm = TRUE),

      from_ensg_missing = any(from_ensg_missing, na.rm = TRUE),
      to_ensg_missing   = any(to_ensg_missing,   na.rm = TRUE),

      from_match_method = paste(sort(unique(from_match_method)), collapse = ";"),
      to_match_method   = paste(sort(unique(to_match_method)),   collapse = ";"),

      n_dir_rows    = .N,
      bidirectional = uniqueN(paste(from_ensg_filled, to_ensg_filled, sep = "->")) > 1
    ),
    by = .(gene1, gene2)
  ]

  out_dt[ensg2symbol_dt, symbol1 := i.symbol, on = c(gene1 = "ensg")]
  out_dt[ensg2symbol_dt, symbol2 := i.symbol, on = c(gene2 = "ensg")]

  for (cc in c(
    "n_references", "n_resources", "curation_effort",
    "consensus_direction", "consensus_stimulation", "consensus_inhibition"
  )) {
    set(out_dt, i = which(!is.finite(out_dt[[cc]])), j = cc, value = NA_real_)
  }

  setcolorder(
    out_dt,
    c(
      "from_ensg_filled", "to_ensg_filled",
      "from_symbol", "to_symbol",
      "gene1", "gene2", "symbol1", "symbol2",
      "interaction_type", "type",
      "n_references", "n_resources", "curation_effort",
      "consensus_direction", "consensus_stimulation", "consensus_inhibition",
      "bidirectional", "n_dir_rows",
      "from_ensg", "to_ensg",
      "from_ensg_missing", "to_ensg_missing",
      "from_match_method", "to_match_method"
    )
  )

  out_dt[]
}

############################################################
## Step 6) De-duplicate edges within each layer
############################################################

cat("\n[Step 6a] De-duplicating exact directed PPI edges...\n")
edges_ppi_dir_dedup_dt <- dedup_directed_edges(edges_ppi_dt)
cat("[INFO] PPI rows:", nrow(edges_ppi_dt), "->", nrow(edges_ppi_dir_dedup_dt), "\n")

cat("\n[Step 6a] De-duplicating exact directed TR edges...\n")
edges_tr_dir_dedup_dt <- dedup_directed_edges(edges_tr_dt)
cat("[INFO] TR rows:", nrow(edges_tr_dt), "->", nrow(edges_tr_dir_dedup_dt), "\n")

cat("\n[Step 6b] Collapsing undirected duplicate PPI edges...\n")
edges_ppi_dedup_dt <- dedup_undirected_edges(edges_ppi_dir_dedup_dt, ensg2symbol_dt)
cat("[INFO] PPI rows after undirected collapse:",
    nrow(edges_ppi_dir_dedup_dt), "->", nrow(edges_ppi_dedup_dt), "\n")

cat("\n[Step 6b] Collapsing undirected duplicate TR edges...\n")
edges_tr_dedup_dt <- dedup_undirected_edges(edges_tr_dir_dedup_dt, ensg2symbol_dt)
cat("[INFO] TR rows after undirected collapse:",
    nrow(edges_tr_dir_dedup_dt), "->", nrow(edges_tr_dedup_dt), "\n")

cat("\n[PPI] duplicated undirected pairs after collapse:",
    edges_ppi_dedup_dt[, sum(duplicated(.SD)), .SDcols = c("gene1", "gene2")], "\n")
cat("[TR] duplicated undirected pairs after collapse:",
    edges_tr_dedup_dt[, sum(duplicated(.SD)), .SDcols = c("gene1", "gene2")], "\n")

cat("[INFO] PPI bidirectional links collapsed:",
    sum(edges_ppi_dedup_dt$bidirectional, na.rm = TRUE), "\n")
cat("[INFO] TR bidirectional links collapsed:",
    sum(edges_tr_dedup_dt$bidirectional, na.rm = TRUE), "\n")

############################################################
## Step 7) Combine PPI and TR layers
##
## Cross-layer duplicate gene pairs are retained intentionally
############################################################

all_cols <- union(names(edges_ppi_dedup_dt), names(edges_tr_dedup_dt))

ppi2 <- copy(edges_ppi_dedup_dt)[, ..all_cols]
tr2  <- copy(edges_tr_dedup_dt)[,  ..all_cols]

edges_all_dt <- rbindlist(list(ppi2, tr2), use.names = TRUE, fill = TRUE)

cat("[INFO] Combined edges_all_dt:", nrow(edges_all_dt), "rows x", ncol(edges_all_dt), "cols\n")

cat("[INFO] duplicated directed pairs after combine:",
    sum(duplicated(edges_all_dt[, .(from_ensg_filled, to_ensg_filled)])), "\n")

dup_links_dt <- edges_all_dt[
  duplicated(edges_all_dt[, .(from_ensg_filled, to_ensg_filled)]) |
    duplicated(edges_all_dt[, .(from_ensg_filled, to_ensg_filled)], fromLast = TRUE)
][order(from_ensg_filled, to_ensg_filled, interaction_type, type)]

cat("[INFO] duplicated links shared across PPI/TR:",
    uniqueN(paste(dup_links_dt$from_ensg_filled, dup_links_dt$to_ensg_filled, sep = "->")), "\n")

if (nrow(dup_links_dt) > 0) {
  print(
    dup_links_dt[, .(
      from_ensg_filled,
      to_ensg_filled,
      from_symbol,
      to_symbol,
      interaction_type,
      type,
      n_references,
      n_resources,
      curation_effort
    )]
  )
}

############################################################
## Step 8) Visualize raw evidence support (n_references)
############################################################

stopifnot(all(c("type", "n_references") %in% names(edges_all_dt)))

edges_all_dt[, type := factor(type, levels = sort(unique(as.character(type))))]

theme_ng <- theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    axis.title = element_text(size = 11),
    axis.text  = element_text(size = 10, color = "black"),
    axis.line  = element_line(linewidth = 0.7, color = "black"),
    axis.ticks = element_line(linewidth = 0.6, color = "black"),
    legend.title = element_text(size = 10),
    legend.text  = element_text(size = 9),
    legend.position = "right"
  )

type_levels <- levels(edges_all_dt$type)
pal <- setNames(hue_pal(l = 55, c = 85)(length(type_levels)), type_levels)

p_den_type_raw <- ggplot(edges_all_dt, aes(x = n_references, color = type)) +
  geom_density(linewidth = 0.9, adjust = 1.0) +
  scale_color_manual(values = pal) +
  labs(
    title = "Raw n_references by interaction type",
    x = "Number of references",
    y = "Density",
    color = "Type"
  ) +
  theme_ng

p_den_all_raw <- ggplot(edges_all_dt, aes(x = n_references)) +
  geom_density(linewidth = 0.9, adjust = 1.0, color = "black") +
  labs(
    title = "Raw n_references across all edges",
    x = "Number of references",
    y = "Density"
  ) +
  theme_ng +
  theme(legend.position = "none")

p_raw_panel <- p_den_type_raw + p_den_all_raw + plot_layout(ncol = 2)
print(p_raw_panel)

############################################################
## Step 9) Compute robust evidence-based edge weights
############################################################

stopifnot(all(c("n_references", "type", "gene1", "gene2") %in% names(edges_all_dt)))

cat("[INFO] duplicated undirected pairs before weighting:",
    edges_all_dt[, sum(duplicated(.SD)), .SDcols = c("type", "gene1", "gene2")], "\n")

cap_q <- 0.99
eps   <- 0.01

edges_all_dt[, weight_raw := log1p(n_references)]

cap_val <- as.numeric(quantile(edges_all_dt$weight_raw, probs = cap_q, na.rm = TRUE))
edges_all_dt[, weight_raw_cap := pmin(weight_raw, cap_val)]

w_min <- min(edges_all_dt$weight_raw_cap, na.rm = TRUE)
w_max <- max(edges_all_dt$weight_raw_cap, na.rm = TRUE)

if (isTRUE(all.equal(w_min, w_max))) {
  edges_all_dt[, weight_01 := 0]
} else {
  edges_all_dt[, weight_01 := (weight_raw_cap - w_min) / (w_max - w_min)]
}

edges_all_dt[, weight_final := eps + (1 - eps) * weight_01]

cat("[INFO] Final weight range:",
    paste(range(edges_all_dt$weight_final, na.rm = TRUE), collapse = " - "), "\n")

############################################################
## Step 10) Visualize processed edge weights
############################################################

p_den_type_w <- ggplot(edges_all_dt, aes(x = weight_final, color = type)) +
  geom_density(linewidth = 0.9, adjust = 1.0) +
  geom_vline(xintercept = eps, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  annotate("text", x = eps, y = Inf, label = paste0("ε = ", eps),
           vjust = 1.3, hjust = -0.1, size = 3.2, color = "grey40") +
  scale_color_manual(values = pal) +
  labs(
    title = "Processed edge weight by interaction type",
    x = "Edge weight",
    y = "Density",
    color = "Type"
  ) +
  theme_ng

p_den_all_w <- ggplot(edges_all_dt, aes(x = weight_final)) +
  geom_density(linewidth = 0.9, adjust = 1.0, color = "black") +
  geom_vline(xintercept = eps, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  annotate("text", x = eps, y = Inf, label = paste0("ε = ", eps),
           vjust = 1.3, hjust = -0.1, size = 3.2, color = "grey40") +
  labs(
    title = "Processed edge weight across all edges",
    x = "Edge weight",
    y = "Density"
  ) +
  theme_ng +
  theme(legend.position = "none")

p_w_panel <- p_den_type_w + p_den_all_w + plot_layout(ncol = 2)
print(p_w_panel)

############################################################
## Step 11) Save outputs
############################################################

ppi_out <- file.path(network_dir, "edges_ppi_dedup_dt.txt")
tr_out  <- file.path(network_dir, "edges_tr_dedup_dt.txt")
all_out <- file.path(network_dir, "edges_all_dt_with_weights.txt")

fwrite(edges_ppi_dedup_dt, file = ppi_out, sep = "\t", quote = FALSE, na = "NA")
fwrite(edges_tr_dedup_dt,  file = tr_out,  sep = "\t", quote = FALSE, na = "NA")
fwrite(edges_all_dt,       file = all_out, sep = "\t", quote = FALSE, na = "NA")

cat("[INFO] Saved:\n")
cat("  - ", ppi_out, "\n", sep = "")
cat("  - ", tr_out,  "\n", sep = "")
cat("  - ", all_out, "\n", sep = "")

cat("[INFO] File size (bytes): ", file.info(all_out)$size, "\n", sep = "")
cat("[INFO] Finished.\n")