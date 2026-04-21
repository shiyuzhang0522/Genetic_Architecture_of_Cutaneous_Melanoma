############################################################
## Construct the curated melanoma interaction network
##
## Purpose:
##   - Read curated melanoma vertices and OmniPath-derived edge tables
##   - Harmonize endpoint metadata across PPI and TR layers
##   - Link curated melanoma genes to molecular interactions
##   - Identify genes not connected in the curated network
##   - Benchmark observed network connectivity against random
##     protein-coding gene sets from the Ensembl universe
##
## Inputs:
##   - Vertices:
##       all_genes_long_pc_dt.txt
##   - Edges:
##       OmniPath_PPI_expanded_with_ENSG.txt
##       OmniPath_TR_expanded_with_ENSG.txt
##
## Outputs:
##   - Linked PPI edge table
##   - Linked TR edge table
##   - Reports of vertices missing from PPI / TR / any linked edge
##   - Ensembl protein-coding universe
##   - Permutation benchmark summary
##
## Author: Shelley
## Date:   2026-01-20
############################################################

############################################################
## Step 1) Load packages
############################################################
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(biomaRt)
  library(igraph)
})

cat("[INFO] Libraries loaded: data.table, dplyr, biomaRt, igraph\n")

############################################################
## Step 2) Read input tables
############################################################

## 2.1 Vertices
vertices_file <- "/Users/shelleyz/Projects/Melanoma-WGS/Updated.Omnigenic.Architecture/Vertices/all_genes_long_pc_dt.txt"
stopifnot(file.exists(vertices_file))

vertices_dt <- fread(vertices_file)

cat("[INFO] Loaded vertices:", vertices_file, "\n")
cat("[INFO] vertices_dt:", nrow(vertices_dt), "rows x", ncol(vertices_dt), "cols\n")
print(head(vertices_dt, 3))

## 2.2 OmniPath PPI edges
ppi_file <- "/Users/shelleyz/Projects/Melanoma-WGS/Omnigenic.Architecture/Edges/OmniPath_PPI_expanded_with_ENSG.txt"
stopifnot(file.exists(ppi_file))

edges_ppi_dt <- fread(ppi_file)

cat("[INFO] Loaded PPI edges:", ppi_file, "\n")
cat("[INFO] edges_ppi_dt:", nrow(edges_ppi_dt), "rows x", ncol(edges_ppi_dt), "cols\n")
print(head(edges_ppi_dt, 3))

## 2.3 OmniPath TR edges
tr_file <- "/Users/shelleyz/Projects/Melanoma-WGS/Omnigenic.Architecture/Edges/OmniPath_TR_expanded_with_ENSG.txt"
stopifnot(file.exists(tr_file))

edges_tr_dt <- fread(tr_file)

cat("[INFO] Loaded TR edges:", tr_file, "\n")
cat("[INFO] edges_tr_dt:", nrow(edges_tr_dt), "rows x", ncol(edges_tr_dt), "cols\n")
print(head(edges_tr_dt, 3))

############################################################
## Step 3) Standardize edge tables to a unified schema
############################################################

required_edge_cols <- c(
  "source_ensg", "target_ensg",
  "source_genesymbol", "target_genesymbol",
  "type",
  "consensus_direction", "consensus_stimulation", "consensus_inhibition",
  "curation_effort", "n_references", "n_resources"
)

missing_ppi <- setdiff(required_edge_cols, names(edges_ppi_dt))
missing_tr  <- setdiff(required_edge_cols, names(edges_tr_dt))

if (length(missing_ppi) > 0) {
  stop("[ERROR] edges_ppi_dt is missing: ", paste(missing_ppi, collapse = ", "))
}
if (length(missing_tr) > 0) {
  stop("[ERROR] edges_tr_dt is missing: ", paste(missing_tr, collapse = ", "))
}

edges_ppi_clean_dt <- edges_ppi_dt[
  ,
  .(
    from_ensg = source_ensg,
    from_symbol = source_genesymbol,
    to_ensg = target_ensg,
    to_symbol = target_genesymbol,
    consensus_direction = consensus_direction,
    consensus_stimulation = consensus_stimulation,
    consensus_inhibition = consensus_inhibition,
    type = type,
    curation_effort = curation_effort,
    n_references = n_references,
    n_resources = n_resources
  )
][, interaction_type := "PPI"]

edges_tr_clean_dt <- edges_tr_dt[
  ,
  .(
    from_ensg = source_ensg,
    from_symbol = source_genesymbol,
    to_ensg = target_ensg,
    to_symbol = target_genesymbol,
    consensus_direction = consensus_direction,
    consensus_stimulation = consensus_stimulation,
    consensus_inhibition = consensus_inhibition,
    type = type,
    curation_effort = curation_effort,
    n_references = n_references,
    n_resources = n_resources
  )
][, interaction_type := "TR"]

cat("[INFO] Standardized edge tables created.\n")
cat("[INFO] edges_ppi_clean_dt:", nrow(edges_ppi_clean_dt), "rows\n")
cat("[INFO] edges_tr_clean_dt :", nrow(edges_tr_clean_dt),  "rows\n")

############################################################
## Step 4) Prepare curated vertex set
############################################################

stopifnot(all(c("ensg", "symbol") %in% names(vertices_dt)))

vertices_dt <- as.data.table(copy(vertices_dt))
vertices_dt <- vertices_dt[!is.na(ensg) & ensg != ""]

v_ensg <- unique(vertices_dt$ensg)
v_sym  <- unique(vertices_dt$symbol[!is.na(vertices_dt$symbol) & vertices_dt$symbol != ""])

cat("[INFO] Curated vertex ENSG count:", length(v_ensg), "\n")

############################################################
## Step 5) Inspect duplicated gene symbols and apply manual fixes
############################################################

dup_symbol_dt <- vertices_dt[
  !is.na(symbol) & symbol != "",
  .N,
  by = symbol
][N > 1][order(-N, symbol)]

cat("[INFO] Duplicated symbols before manual fix:", nrow(dup_symbol_dt), "\n")
print(dup_symbol_dt)

## manual symbol disambiguation
vertices_dt[ensg == "ENSG00000250264", symbol := "TAP2-HLA-DOB"]
vertices_dt[ensg == "ENSG00000198211", symbol := "MC1R-TUBB3"]

cat("[INFO] Manual symbol fixes applied.\n")

dup_symbol_dt <- vertices_dt[
  !is.na(symbol) & symbol != "",
  .N,
  by = symbol
][N > 1][order(-N, symbol)]

cat("[INFO] Duplicated symbols after manual fix:", nrow(dup_symbol_dt), "\n")
print(dup_symbol_dt)

if (nrow(dup_symbol_dt) > 0) {
  dup_symbol_rows_dt <- vertices_dt[symbol %in% dup_symbol_dt$symbol][order(symbol, ensg)]
  cat("[INFO] Full rows for remaining duplicated symbols:\n")
  print(dup_symbol_rows_dt)
}

############################################################
## Step 6) Build deterministic symbol -> ENSG map
##
## Used only as a fallback when an edge endpoint ENSG is missing
############################################################

sym2ensg_dt <- unique(
  vertices_dt[
    !is.na(symbol) & symbol != "" &
      !is.na(ensg) & ensg != "",
    .(symbol, ensg)
  ]
)

setkey(sym2ensg_dt, symbol)

cat("[INFO] sym2ensg_dt rows:", nrow(sym2ensg_dt), "\n")

############################################################
## Step 7) Fill missing edge endpoints using symbol fallback
############################################################

fill_endpoint_ensg <- function(edge_dt,
                               sym2ensg_dt,
                               from_ensg_col = "from_ensg",
                               to_ensg_col   = "to_ensg",
                               from_sym_col  = "from_symbol",
                               to_sym_col    = "to_symbol") {
  dt <- copy(edge_dt)
  stopifnot(all(c(from_ensg_col, to_ensg_col, from_sym_col, to_sym_col) %in% names(dt)))

  dt[, from_ensg_missing := is.na(get(from_ensg_col)) | get(from_ensg_col) == ""]
  dt[, to_ensg_missing   := is.na(get(to_ensg_col))   | get(to_ensg_col) == ""]

  dt[from_ensg_missing == TRUE,
     from_ensg_filled := sym2ensg_dt[match(get(from_sym_col), symbol), ensg]]
  dt[from_ensg_missing == FALSE,
     from_ensg_filled := get(from_ensg_col)]

  dt[to_ensg_missing == TRUE,
     to_ensg_filled := sym2ensg_dt[match(get(to_sym_col), symbol), ensg]]
  dt[to_ensg_missing == FALSE,
     to_ensg_filled := get(to_ensg_col)]

  dt[, from_match_method := fifelse(from_ensg_missing, "symbol_fallback", "ensg_direct")]
  dt[, to_match_method   := fifelse(to_ensg_missing,   "symbol_fallback", "ensg_direct")]

  dt[]
}

edges_ppi_filled_dt <- fill_endpoint_ensg(
  edge_dt = edges_ppi_clean_dt,
  sym2ensg_dt = sym2ensg_dt
)

edges_tr_filled_dt <- fill_endpoint_ensg(
  edge_dt = edges_tr_clean_dt,
  sym2ensg_dt = sym2ensg_dt
)

############################################################
## Step 8) Subset to interactions internal to curated vertices
############################################################

subset_edges_to_vertices <- function(edge_dt_filled, vertex_ensg) {
  dt <- copy(edge_dt_filled)

  dt_valid <- dt[
    !is.na(from_ensg_filled) & from_ensg_filled != "" &
      !is.na(to_ensg_filled) & to_ensg_filled != ""
  ]

  dt_linked <- dt_valid[
    from_ensg_filled %in% vertex_ensg &
      to_ensg_filled %in% vertex_ensg
  ]

  dt_linked[]
}

edges_ppi_linked_dt <- subset_edges_to_vertices(edges_ppi_filled_dt, v_ensg)
edges_tr_linked_dt  <- subset_edges_to_vertices(edges_tr_filled_dt,  v_ensg)

cat("[INFO] PPI edges linked to curated vertices:", nrow(edges_ppi_linked_dt), "\n")
cat("[INFO] TR edges linked to curated vertices :", nrow(edges_tr_linked_dt),  "\n")

cat("[INFO] PPI linked unique vertices:",
    uniqueN(c(edges_ppi_linked_dt$from_ensg_filled, edges_ppi_linked_dt$to_ensg_filled)), "\n")
cat("[INFO] TR linked unique vertices:",
    uniqueN(c(edges_tr_linked_dt$from_ensg_filled, edges_tr_linked_dt$to_ensg_filled)), "\n")

cat("[INFO] PPI endpoint match methods:\n")
print(list(
  from = table(edges_ppi_linked_dt$from_match_method, useNA = "ifany"),
  to   = table(edges_ppi_linked_dt$to_match_method,   useNA = "ifany")
))

cat("[INFO] TR endpoint match methods:\n")
print(list(
  from = table(edges_tr_linked_dt$from_match_method, useNA = "ifany"),
  to   = table(edges_tr_linked_dt$to_match_method,   useNA = "ifany")
))

############################################################
## Step 9) Identify curated vertices missing from linked network
############################################################

ppi_linked_ensg <- unique(c(edges_ppi_linked_dt$from_ensg_filled, edges_ppi_linked_dt$to_ensg_filled))
tr_linked_ensg  <- unique(c(edges_tr_linked_dt$from_ensg_filled,  edges_tr_linked_dt$to_ensg_filled))

ppi_linked_ensg <- ppi_linked_ensg[!is.na(ppi_linked_ensg) & ppi_linked_ensg != ""]
tr_linked_ensg  <- tr_linked_ensg[!is.na(tr_linked_ensg) & tr_linked_ensg != ""]

missing_in_ppi <- setdiff(v_ensg, ppi_linked_ensg)
missing_in_tr  <- setdiff(v_ensg, tr_linked_ensg)
missing_in_any <- setdiff(v_ensg, union(ppi_linked_ensg, tr_linked_ensg))

missing_in_ppi_dt <- vertices_dt[ensg %in% missing_in_ppi][order(symbol, ensg)]
missing_in_tr_dt  <- vertices_dt[ensg %in% missing_in_tr ][order(symbol, ensg)]
missing_in_any_dt <- vertices_dt[ensg %in% missing_in_any][order(symbol, ensg)]

cat("[INFO] Vertices missing in linked PPI:", nrow(missing_in_ppi_dt), "\n")
cat("[INFO] Vertices missing in linked TR :", nrow(missing_in_tr_dt),  "\n")
cat("[INFO] Vertices missing in linked ANY:", nrow(missing_in_any_dt), "\n")

############################################################
## Step 10) Save linked edge tables and missing-vertex reports
############################################################

out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/Updated.Omnigenic.Architecture/Network"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ppi_out <- file.path(out_dir, "edges_ppi_linked_dt.txt")
tr_out  <- file.path(out_dir, "edges_tr_linked_dt.txt")

miss_ppi_out <- file.path(out_dir, "vertices_missing_in_PPI.txt")
miss_tr_out  <- file.path(out_dir, "vertices_missing_in_TR.txt")
miss_any_out <- file.path(out_dir, "vertices_missing_in_any_OP_edge.txt")

fwrite(edges_ppi_linked_dt, ppi_out, sep = "\t", quote = FALSE, na = "NA")
fwrite(edges_tr_linked_dt,  tr_out,  sep = "\t", quote = FALSE, na = "NA")

fwrite(missing_in_ppi_dt, miss_ppi_out, sep = "\t", quote = FALSE, na = "NA")
fwrite(missing_in_tr_dt,  miss_tr_out,  sep = "\t", quote = FALSE, na = "NA")
fwrite(missing_in_any_dt, miss_any_out, sep = "\t", quote = FALSE, na = "NA")

cat("[INFO] Saved linked edge tables and missing-vertex reports.\n")

############################################################
## Step 11) Build Ensembl protein-coding universe for benchmarking
############################################################

mart <- useEnsembl(
  biomart = "genes",
  dataset = "hsapiens_gene_ensembl"
)

attr_dt <- listAttributes(mart)
attr_names <- attr_dt$name

if ("gene_biotype" %in% attr_names) {
  biotype_attr <- "gene_biotype"
} else if ("transcript_biotype" %in% attr_names) {
  biotype_attr <- "transcript_biotype"
} else {
  stop("[ERROR] Neither gene_biotype nor transcript_biotype is available in this BioMart dataset.")
}

cat("[INFO] Using BioMart biotype attribute:", biotype_attr, "\n")

bm_raw <- getBM(
  attributes = c("ensembl_gene_id", "hgnc_symbol", biotype_attr),
  mart = mart
)

bm_dt <- as.data.table(bm_raw)
setnames(bm_dt, "ensembl_gene_id", "ensg")
setnames(bm_dt, "hgnc_symbol", "symbol")
setnames(bm_dt, biotype_attr, "biotype")

bm_dt <- bm_dt[!is.na(ensg) & ensg != ""]

pc_universe_dt <- unique(bm_dt[biotype == "protein_coding", .(ensg, symbol, biotype)])
pc_universe <- unique(pc_universe_dt$ensg)

cat("[INFO] Ensembl protein-coding universe size:", length(pc_universe), "\n")

pc_universe_out <- file.path(out_dir, "ensembl_pc_universe.txt")
fwrite(pc_universe_dt, pc_universe_out, sep = "\t", quote = FALSE, na = "NA")

cat("[INFO] Saved protein-coding universe ->", pc_universe_out, "\n")

############################################################
## Step 12) Sanity check: observed genes must be in protein-coding universe
############################################################

missing_in_pc_universe <- setdiff(v_ensg, pc_universe)
cat("[INFO] Observed ENSG not found in pc_universe:", length(missing_in_pc_universe), "\n")

if (length(missing_in_pc_universe) > 0) {
  print(missing_in_pc_universe)
  stop("[ERROR] Some observed ENSG are not in the Ensembl protein-coding universe.")
}

############################################################
## Step 13) Define permutation edge pools
##
## Use the filled edge tables, not the already-linked subset
############################################################

edges_ppi_perm_dt <- edges_ppi_filled_dt[
  !is.na(from_ensg_filled) & from_ensg_filled != "" &
    !is.na(to_ensg_filled) & to_ensg_filled != ""
]

edges_tr_perm_dt <- edges_tr_filled_dt[
  !is.na(from_ensg_filled) & from_ensg_filled != "" &
    !is.na(to_ensg_filled) & to_ensg_filled != ""
]

setkey(edges_ppi_perm_dt, from_ensg_filled, to_ensg_filled)
setkey(edges_tr_perm_dt,  from_ensg_filled, to_ensg_filled)

cat("[INFO] PPI permutation edge pool:", nrow(edges_ppi_perm_dt), "\n")
cat("[INFO] TR permutation edge pool :", nrow(edges_tr_perm_dt),  "\n")

############################################################
## Step 14) Helper to compute undirected network metrics
############################################################

calc_network_metrics <- function(gset, ppi_dt, tr_dt) {
  gset <- unique(gset)
  gset <- gset[!is.na(gset) & gset != ""]

  ## PPI as undirected
  ppi_sub <- ppi_dt[
    from_ensg_filled %in% gset & to_ensg_filled %in% gset,
    .(
      gene1 = pmin(from_ensg_filled, to_ensg_filled),
      gene2 = pmax(from_ensg_filled, to_ensg_filled)
    )
  ]
  ppi_sub <- unique(ppi_sub[gene1 != gene2])

  ## TR converted to undirected for summary metrics
  tr_sub <- tr_dt[
    from_ensg_filled %in% gset & to_ensg_filled %in% gset,
    .(
      gene1 = pmin(from_ensg_filled, to_ensg_filled),
      gene2 = pmax(from_ensg_filled, to_ensg_filled)
    )
  ]
  tr_sub <- unique(tr_sub[gene1 != gene2])

  edge_any_dt <- unique(rbindlist(list(ppi_sub, tr_sub), use.names = TRUE, fill = TRUE))
  vertex_df <- data.frame(name = gset, stringsAsFactors = FALSE)

  if (nrow(edge_any_dt) == 0) {
    g <- make_empty_graph(n = length(gset), directed = FALSE)
    g <- set_vertex_attr(g, "name", value = gset)
  } else {
    g <- graph_from_data_frame(
      d = as.data.frame(edge_any_dt),
      directed = FALSE,
      vertices = vertex_df
    )
  }

  n_genes <- gorder(g)
  n_edges <- gsize(g)
  network_density <- if (n_genes > 1) edge_density(g, loops = FALSE) else NA_real_

  deg <- degree(g, mode = "all", loops = FALSE)
  mean_degree <- mean(deg)

  comp <- components(g, mode = "weak")
  largest_component <- if (length(comp$csize) > 0) max(comp$csize) else 0L

  linked_genes <- sum(deg > 0)
  unlinked_genes <- sum(deg == 0)

  data.table(
    n_genes = n_genes,
    n_edges = n_edges,
    linked_genes = linked_genes,
    unlinked_genes = unlinked_genes,
    network_density = network_density,
    mean_degree = mean_degree,
    largest_component = largest_component
  )
}

############################################################
## Step 15) Compute observed network metrics
############################################################

obs_res <- calc_network_metrics(
  gset = v_ensg,
  ppi_dt = edges_ppi_perm_dt,
  tr_dt  = edges_tr_perm_dt
)

cat("[INFO] Observed network metrics:\n")
print(obs_res)

############################################################
## Step 16) Permutation benchmark
############################################################

set.seed(20260306)

n_perm <- 5000L
perm_res_list <- vector("list", n_perm)

for (b in seq_len(n_perm)) {
  gset_b <- sample(pc_universe, size = length(v_ensg), replace = FALSE)

  perm_res_list[[b]] <- cbind(
    calc_network_metrics(
      gset = gset_b,
      ppi_dt = edges_ppi_perm_dt,
      tr_dt  = edges_tr_perm_dt
    ),
    data.table(perm = b)
  )

  if (b %% 100L == 0L) {
    cat("[INFO] Permutation ", b, "/", n_perm, "\n", sep = "")
  }
}

perm_res_dt <- rbindlist(perm_res_list, use.names = TRUE, fill = TRUE)

############################################################
## Step 17) Empirical p-values
############################################################

emp_p_ge <- function(null_vec, obs_val) {
  (sum(null_vec >= obs_val, na.rm = TRUE) + 1) / (sum(!is.na(null_vec)) + 1)
}

p_density <- emp_p_ge(perm_res_dt$network_density,   obs_res$network_density)
p_degree  <- emp_p_ge(perm_res_dt$mean_degree,       obs_res$mean_degree)
p_lcc     <- emp_p_ge(perm_res_dt$largest_component, obs_res$largest_component)

cat("[INFO] Empirical p for network_density  :", signif(p_density, 3), "\n")
cat("[INFO] Empirical p for mean_degree      :", signif(p_degree,  3), "\n")
cat("[INFO] Empirical p for largest_component:", signif(p_lcc,     3), "\n")

############################################################
## Step 18) Summarize permutation results
############################################################

perm_summary_dt <- data.table(
  metric = c("network_density", "mean_degree", "largest_component"),
  observed = c(
    obs_res$network_density,
    obs_res$mean_degree,
    obs_res$largest_component
  ),
  null_mean = c(
    mean(perm_res_dt$network_density, na.rm = TRUE),
    mean(perm_res_dt$mean_degree, na.rm = TRUE),
    mean(perm_res_dt$largest_component, na.rm = TRUE)
  ),
  null_sd = c(
    sd(perm_res_dt$network_density, na.rm = TRUE),
    sd(perm_res_dt$mean_degree, na.rm = TRUE),
    sd(perm_res_dt$largest_component, na.rm = TRUE)
  ),
  null_q05 = c(
    quantile(perm_res_dt$network_density, 0.05, na.rm = TRUE),
    quantile(perm_res_dt$mean_degree, 0.05, na.rm = TRUE),
    quantile(perm_res_dt$largest_component, 0.05, na.rm = TRUE)
  ),
  null_q50 = c(
    quantile(perm_res_dt$network_density, 0.50, na.rm = TRUE),
    quantile(perm_res_dt$mean_degree, 0.50, na.rm = TRUE),
    quantile(perm_res_dt$largest_component, 0.50, na.rm = TRUE)
  ),
  null_q95 = c(
    quantile(perm_res_dt$network_density, 0.95, na.rm = TRUE),
    quantile(perm_res_dt$mean_degree, 0.95, na.rm = TRUE),
    quantile(perm_res_dt$largest_component, 0.95, na.rm = TRUE)
  ),
  emp_p_ge = c(p_density, p_degree, p_lcc)
)

cat("[INFO] Permutation summary:\n")
print(perm_summary_dt)

############################################################
## Step 19) Save permutation summary
############################################################

perm_out <- file.path(out_dir, "permutation_summary.txt")

fwrite(
  perm_summary_dt,
  perm_out,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("[INFO] Saved permutation summary ->", perm_out, "\n")
cat("[INFO] Finished.\n")