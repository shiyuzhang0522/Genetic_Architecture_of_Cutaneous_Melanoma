############################################################
## Refine the melanoma interaction network and construct the
## integrated undirected graph
##
## Purpose:
##   - Collapse layer-specific edges to unique gene-pair links
##   - Integrate evidence from PPI and TR for the same gene pair
##   - Build a graph object with curated vertex metadata
##   - Summarize global and node-level network topology
##
## Notes:
##   - PPI/TR evidence is first retained at the edge-row level
##   - Gene-pair integration is then performed across duplicated
##     gene1-gene2 pairs
##   - For gene pairs supported by both PPI and TR, the final edge
##     weight is computed using a probabilistic union:
##         1 - prod(1 - weight_final)
##
## Author: Shelley
## Date:   2026-03-06
############################################################

############################################################
## Step 1) Load packages
############################################################
suppressPackageStartupMessages({
  library(data.table)
  library(fs)
  library(igraph)
})

cat("[INFO] Packages loaded: data.table, fs, igraph\n")

############################################################
## Step 2) Read input tables
############################################################
network_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/Updated.Omnigenic.Architecture/Network"
vertex_file <- "/Users/shelleyz/Projects/Melanoma-WGS/Updated.Omnigenic.Architecture/Vertices/all_genes_long_pc_dt.txt"
edge_file   <- file.path(network_dir, "edges_all_dt_with_weights.txt")

stopifnot(file_exists(edge_file))
stopifnot(file_exists(vertex_file))

edges_all_dt <- fread(edge_file)
all_genes_long_pc_dt <- fread(vertex_file)

cat("[INFO] Loaded edges_all_dt:", nrow(edges_all_dt), "rows x", ncol(edges_all_dt), "cols\n")
cat("[INFO] Loaded all_genes_long_pc_dt:", nrow(all_genes_long_pc_dt), "rows x", ncol(all_genes_long_pc_dt), "cols\n")

############################################################
## Step 3) Subset columns of interest
############################################################
keep_cols <- c(
  "gene1", "gene2",
  "symbol1", "symbol2",
  "interaction_type", "type",
  "n_references", "n_resources", "curation_effort",
  "weight_raw", "weight_raw_cap", "weight_01", "weight_final"
)

stopifnot(all(keep_cols %in% names(edges_all_dt)))

edges_sub_dt <- copy(edges_all_dt)[, ..keep_cols]

cat("[INFO] edges_sub_dt:", nrow(edges_sub_dt), "rows x", ncol(edges_sub_dt), "cols\n")

############################################################
## Step 4) Sanity checks on edge table
############################################################

## unique genes / symbols
unique_ensg <- unique(c(edges_sub_dt$gene1, edges_sub_dt$gene2))
unique_symbols <- unique(c(edges_sub_dt$symbol1, edges_sub_dt$symbol2))

cat("[INFO] Number of unique ENSG IDs:", length(unique_ensg), "\n")
cat("[INFO] Number of unique gene symbols:", length(unique_symbols), "\n")

## interaction_type vs type consistency
tab <- table(edges_sub_dt$interaction_type, edges_sub_dt$type)
print(tab)

mismatch_dt <- edges_sub_dt[
  (interaction_type == "PPI" & type != "post_translational") |
  (interaction_type == "TR"  & type != "transcriptional")
]

cat("[INFO] Number of inconsistent rows:", nrow(mismatch_dt), "\n")
if (nrow(mismatch_dt) > 0) {
  print(head(mismatch_dt))
  stop("[ERROR] Found inconsistent interaction_type vs type mappings.")
}

## duplicated gene-pair rows before collapse
dup_pairs_dt <- edges_sub_dt[
  , .N, by = .(gene1, gene2)
][N > 1, .(gene1, gene2, dup_n = N)][order(-dup_n, gene1, gene2)]

dup_rows_dt <- edges_sub_dt[
  dup_pairs_dt,
  on = .(gene1, gene2),
  nomatch = 0
][order(gene1, gene2)]

cat("[INFO] Number of duplicated gene1-gene2 pairs:", nrow(dup_pairs_dt), "\n")
cat("[INFO] Number of rows in duplicated subset:", nrow(dup_rows_dt), "\n")

if (nrow(dup_rows_dt) > 0) {
  print(dup_rows_dt)
}

############################################################
## Step 5) Collapse to gene-pair level and integrate evidence
############################################################
pair_final_dt <- edges_sub_dt[
  ,
  {
    has_ppi <- any(interaction_type == "PPI")
    has_tr  <- any(interaction_type == "TR")

    type_final <- fifelse(
      has_ppi & has_tr, "both",
      fifelse(has_ppi, "post_translational", "transcriptional")
    )

    .(
      symbol1 = first(na.omit(symbol1)),
      symbol2 = first(na.omit(symbol2)),
      interaction_type_combined = paste(sort(unique(interaction_type)), collapse = ";"),
      type_combined = paste(sort(unique(type)), collapse = ";"),
      n_evidence_rows = .N,
      n_references_sum = sum(n_references, na.rm = TRUE),
      n_resources_sum = sum(n_resources, na.rm = TRUE),
      curation_effort_sum = sum(curation_effort, na.rm = TRUE),
      type_final = type_final,
      weight_final_integrated = 1 - prod(1 - weight_final)
    )
  },
  by = .(gene1, gene2)
][order(gene1, gene2)]

cat("[INFO] pair_final_dt:", nrow(pair_final_dt), "rows x", ncol(pair_final_dt), "cols\n")

## uniqueness check after collapse
pair_dup_check_dt <- pair_final_dt[
  , .N, by = .(gene1, gene2)
][N > 1]

cat("[INFO] Number of duplicated gene1-gene2 pairs after integration:", nrow(pair_dup_check_dt), "\n")
if (nrow(pair_dup_check_dt) > 0) {
  print(pair_dup_check_dt)
  stop("[ERROR] Duplicated gene-pair rows remain after integration.")
} else {
  cat("[INFO] All gene1-gene2 combinations are unique in pair_final_dt.\n")
}

############################################################
## Step 6) Save integrated gene-pair edge table
############################################################
pair_out_file <- file.path(network_dir, "Integrated_gene_pair_network_edges.txt")
fwrite(pair_final_dt, pair_out_file, sep = "\t", quote = FALSE, na = "NA")

cat("[INFO] Saved integrated gene-pair table:\n  ", pair_out_file, "\n", sep = "")
cat("[INFO] Rows:", nrow(pair_final_dt), " Columns:", ncol(pair_final_dt), "\n")

############################################################
## Step 7) Prepare vertex metadata
############################################################
stopifnot(all(c("ensg", "symbol", "type", "source") %in% names(all_genes_long_pc_dt)))

vertices_dt <- all_genes_long_pc_dt[
  ,
  {
    src_u <- sort(unique(source))

    source_final <- fifelse(
      all(c("GWAS", "RV_PC") %in% src_u), "GWAS+RV",
      fifelse("GWAS" %in% src_u, "GWAS", "RV_PC")
    )

    .(
      symbol = first(na.omit(symbol)),
      type = first(na.omit(type)),
      source_final = source_final
    )
  },
  by = .(name = ensg)
][order(name)]

cat("[INFO] Number of unique vertices in metadata:", nrow(vertices_dt), "\n")
print(table(vertices_dt$source_final))

############################################################
## Step 8) Keep only vertices present in the edge list
############################################################
edge_nodes <- unique(c(pair_final_dt$gene1, pair_final_dt$gene2))
vertices_sub_dt <- vertices_dt[name %in% edge_nodes][order(name)]

cat("[INFO] Number of unique nodes in edge list:", length(edge_nodes), "\n")
cat("[INFO] Number of vertices retained for graph:", nrow(vertices_sub_dt), "\n")

############################################################
## Step 9) Prepare edge table for graph construction
############################################################
edges_for_graph_dt <- unique(
  pair_final_dt[, .(
    from = gene1,
    to = gene2,
    symbol1,
    symbol2,
    weight = weight_final_integrated,
    type_final,
    interaction_type_combined,
    type_combined,
    n_evidence_rows,
    n_references_sum,
    n_resources_sum,
    curation_effort_sum
  )]
)

cat("[INFO] Number of edges for graph:", nrow(edges_for_graph_dt), "\n")
print(head(edges_for_graph_dt))

## vertex coverage check
missing_nodes <- setdiff(
  unique(c(edges_for_graph_dt$from, edges_for_graph_dt$to)),
  vertices_sub_dt$name
)

cat("[INFO] Number of edge nodes missing in vertex table:", length(missing_nodes), "\n")
if (length(missing_nodes) > 0) {
  print(missing_nodes)
  stop("[ERROR] Some edge nodes are not represented in vertices_sub_dt.")
}

############################################################
## Step 10) Create graph object
############################################################
g_net <- graph_from_data_frame(
  d = edges_for_graph_dt,
  vertices = vertices_sub_dt,
  directed = FALSE
)

print(g_net)
cat("[INFO] Number of vertices in graph:", vcount(g_net), "\n")
cat("[INFO] Number of edges in graph:", ecount(g_net), "\n")
cat("[INFO] Is weighted graph:", "weight" %in% edge_attr_names(g_net), "\n")
cat("[INFO] Vertex attributes:", paste(vertex_attr_names(g_net), collapse = ", "), "\n")
cat("[INFO] Edge attributes:", paste(edge_attr_names(g_net), collapse = ", "), "\n")

############################################################
## Step 11) Save graph object and final tables
############################################################
graph_out_file   <- file.path(network_dir, "g_net_integrated_undirected.RData")
vertex_out_file  <- file.path(network_dir, "vertices_sub_dt_for_graph.txt")
edgegraph_out_file <- file.path(network_dir, "edges_for_graph_dt.txt")

dir_create(network_dir)

save(g_net, file = graph_out_file)
fwrite(vertices_sub_dt, vertex_out_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(edges_for_graph_dt, edgegraph_out_file, sep = "\t", quote = FALSE, na = "NA")

cat("[INFO] Saved graph object:\n  ", graph_out_file, "\n", sep = "")
cat("[INFO] Saved graph vertices:\n  ", vertex_out_file, "\n", sep = "")
cat("[INFO] Saved graph edges:\n  ", edgegraph_out_file, "\n", sep = "")

############################################################
## Step 12) Global network topology
############################################################

n_nodes <- vcount(g_net)
n_edges <- ecount(g_net)

edge_type_dt <- data.table(
  type_final = names(table(E(g_net)$type_final)),
  n_edges = as.integer(table(E(g_net)$type_final))
)[order(type_final)]

edge_type_dt[, prop_edges := n_edges / sum(n_edges)]

net_density <- edge_density(g_net, loops = FALSE)
avg_degree  <- mean(degree(g_net, mode = "all"))

comp_obj <- components(g_net, mode = "weak")

comp_size_dt <- data.table(
  component_id = seq_along(comp_obj$csize),
  component_size = as.integer(comp_obj$csize)
)[order(-component_size)]

gcc_size <- max(comp_obj$csize)
gcc_id   <- which.max(comp_obj$csize)
gcc_prop <- gcc_size / n_nodes

network_basic_dt <- data.table(
  n_nodes = n_nodes,
  n_edges = n_edges,
  density = net_density,
  average_degree = avg_degree,
  n_connected_components = comp_obj$no,
  gcc_id = gcc_id,
  gcc_size = gcc_size,
  gcc_prop = gcc_prop
)

cat("Basic network summary:\n")
print(network_basic_dt)

cat("\nEdge counts by type:\n")
print(edge_type_dt)

cat("\nConnected component sizes:\n")
print(comp_size_dt)

cat("\nCheck: any isolated nodes? ", sum(degree(g_net) == 0), "\n")

############################################################
## Step 13) Node-level topology metrics
############################################################
node_topology_dt <- data.table(
  name = V(g_net)$name,
  symbol = V(g_net)$symbol,
  type = V(g_net)$type,
  source_final = V(g_net)$source_final,
  degree = degree(g_net, mode = "all"),
  strength = strength(g_net, mode = "all", weights = E(g_net)$weight),
  coreness = coreness(g_net, mode = "all")
)

setorder(node_topology_dt, -coreness, -degree, -strength, symbol)

cat("Coreness summary:\n")
print(summary(node_topology_dt$coreness))

cat("\nTop nodes by coreness:\n")
print(head(node_topology_dt, 20))

############################################################
## Step 14) Save topology outputs
############################################################
analysis_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/Updated.Omnigenic.Architecture/Analysis"
dir_create(analysis_dir)

node_out_file <- file.path(analysis_dir, "Node_Topology_Metrics.txt")
net_out_file  <- file.path(analysis_dir, "Network_Basic_Summary.txt")
etype_out_file <- file.path(analysis_dir, "Network_Edge_Type_Summary.txt")
comp_out_file <- file.path(analysis_dir, "Network_Component_Sizes.txt")

fwrite(node_topology_dt, node_out_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(network_basic_dt, net_out_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(edge_type_dt, etype_out_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(comp_size_dt, comp_out_file, sep = "\t", quote = FALSE, na = "NA")

cat("[INFO] Saved node topology table:\n  ", node_out_file, "\n", sep = "")
cat("[INFO] Saved network summary table:\n  ", net_out_file, "\n", sep = "")
cat("[INFO] Saved edge type summary table:\n  ", etype_out_file, "\n", sep = "")
cat("[INFO] Saved component size table:\n  ", comp_out_file, "\n", sep = "")