############################################################
## R script: Identify and summarize network modules in the
## melanoma genetic architecture network
##
## Purpose:
##   Part I
##   - Detect community structure using consensus Louvain
##     clustering on the unweighted topology of the integrated
##     undirected network
##   - Save module assignments and per-module gene lists
##
##   Part II
##   - Generate a compact module-size panel directly from the
##     module summary table
##
## Notes:
##   - Community detection is performed on the unweighted graph
##     topology by setting all edge weights to 1
##   - This prioritizes modular structure in connectivity rather
##     than weighted evidence strength
##
## Author: Shelley (Shiyu Zhang)
## Date:   2026-03-09
############################################################

############################################################
## Part I. Consensus Louvain module detection
############################################################

############################################################
## Step 1) Load required packages
############################################################
suppressPackageStartupMessages({
  library(data.table)
  library(igraph)
  library(EGAnet)
  library(ggplot2)
  library(fs)
})

cat("[INFO] Packages loaded successfully.\n")

############################################################
## Step 2) Read network object into R
############################################################
net_file <- "/Users/shelleyz/Projects/Melanoma-WGS/Updated.Omnigenic.Architecture/Network/g_net_integrated_undirected.RData"
stopifnot(file.exists(net_file))

load(net_file)

g <- g_net

stopifnot(exists("g"))
stopifnot(is.igraph(g))

cat("[INFO] Network object loaded.\n")
cat("[INFO] |V| =", gorder(g), "nodes\n")
cat("[INFO] |E| =", gsize(g), "edges\n")

############################################################
## Step 3) Sanity checks for graph structure and attributes
############################################################
n_isolates <- sum(degree(g) == 0)
comp <- components(g)

cat("[INFO] # isolated nodes:", n_isolates, "\n")
cat("[INFO] # connected components:", comp$no, "\n")
cat("[INFO] Degree summary:\n")
print(summary(degree(g)))

cat("[INFO] Node attributes:\n")
print(vertex_attr_names(g))

cat("[INFO] Edge attributes:\n")
print(edge_attr_names(g))

## node metadata preview
node_attr_dt <- data.table(name = V(g)$name)
for (att in setdiff(vertex_attr_names(g), "name")) {
  node_attr_dt[[att]] <- vertex_attr(g, att)
}

cat("[INFO] Node attribute preview:\n")
print(head(node_attr_dt, 10))
cat("[INFO] Node attribute dimensions:", nrow(node_attr_dt), "x", ncol(node_attr_dt), "\n")

## edge metadata preview
edge_ends_mat <- ends(g, es = E(g), names = TRUE)
edge_attr_dt <- data.table(
  from = edge_ends_mat[, 1],
  to   = edge_ends_mat[, 2]
)
for (att in edge_attr_names(g)) {
  edge_attr_dt[[att]] <- edge_attr(g, att)
}

cat("[INFO] Edge attribute preview:\n")
print(head(edge_attr_dt, 10))
cat("[INFO] Edge attribute dimensions:", nrow(edge_attr_dt), "x", ncol(edge_attr_dt), "\n")

cat("[INFO] Missing values in node attributes:\n")
print(colSums(is.na(node_attr_dt)))

cat("[INFO] Missing values in edge attributes:\n")
print(colSums(is.na(edge_attr_dt)))

if ("weight" %in% edge_attr_names(g)) {
  cat("[INFO] Edge weight summary:\n")
  print(summary(E(g)$weight))
}

############################################################
## Step 4) Create effectively unweighted graph
############################################################
g_unw <- g
E(g_unw)$weight <- 1

stopifnot("weight" %in% edge_attr_names(g_unw))
stopifnot(all(E(g_unw)$weight == 1))

cat("[INFO] Effectively unweighted network created.\n")
cat("[INFO] |V| =", gorder(g_unw), "\n")
cat("[INFO] |E| =", gsize(g_unw), "\n")

############################################################
## Step 5) Run consensus Louvain clustering
############################################################
set.seed(1)

consensus_iter <- 1000
order_use <- "higher"
consensus_method_use <- "iterative"
resolution_use <- 1

mem_cons <- community.consensus(
  network          = g_unw,
  order            = order_use,
  resolution       = resolution_use,
  consensus.method = consensus_method_use,
  consensus.iter   = consensus_iter,
  allow.singleton  = FALSE,
  membership.only  = TRUE
)

stopifnot(length(mem_cons) == gorder(g_unw))

V(g_unw)$module_louvain_cons <- as.integer(mem_cons)

cat("[INFO] Consensus Louvain completed successfully.\n")
cat("[INFO] Number of modules detected:",
    length(unique(V(g_unw)$module_louvain_cons)), "\n")

############################################################
## Step 6) Summarize module assignments
############################################################
node_module_dt <- data.table(
  name                = V(g_unw)$name,
  symbol              = V(g_unw)$symbol,
  type                = V(g_unw)$type,
  source_final        = V(g_unw)$source_final,
  module_louvain_cons = V(g_unw)$module_louvain_cons
)

setorder(node_module_dt, module_louvain_cons, symbol)

cat("[INFO] node_module_dt created successfully.\n")
cat("[INFO] Dimensions:", nrow(node_module_dt), "x", ncol(node_module_dt), "\n")
print(head(node_module_dt, 10))

module_summary_dt <- node_module_dt[
  ,
  .(n_genes = .N),
  by = module_louvain_cons
][order(-n_genes)]

cat("[INFO] Module summary table:\n")
print(module_summary_dt)

mod_sizes <- sort(table(V(g_unw)$module_louvain_cons), decreasing = TRUE)
cat("[INFO] Module size summary:\n")
print(mod_sizes)

############################################################
## Step 7) Save module outputs
############################################################
out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/Updated.Omnigenic.Architecture/Analysis/Modules"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## node-module table
node_module_file <- file.path(out_dir, "node_module_consensus_louvain.txt")
fwrite(node_module_dt, file = node_module_file, sep = "\t", quote = FALSE, na = "NA")
cat("[INFO] Saved node_module_dt ->", node_module_file, "\n")

## module size summary
module_size_file <- file.path(out_dir, "module_size_summary.txt")
fwrite(module_summary_dt, file = module_size_file, sep = "\t", quote = FALSE, na = "NA")
cat("[INFO] Saved module summary ->", module_size_file, "\n")

## one full gene list per module
module_ids <- sort(unique(node_module_dt$module_louvain_cons))
for (mod in module_ids) {
  mod_dt <- unique(
    node_module_dt[module_louvain_cons == mod, .(
      name, symbol, type, source_final
    )]
  )

  mod_file <- file.path(out_dir, paste0("module_", mod, "_gene_list.txt"))
  fwrite(mod_dt, file = mod_file, sep = "\t", quote = FALSE, na = "NA")

  cat("[INFO] Saved module", mod, "gene list:",
      nrow(mod_dt), "genes ->", mod_file, "\n")
}

## one symbol-only gene list per module
for (mod in module_ids) {
  mod_symbols <- node_module_dt[
    module_louvain_cons == mod & !is.na(symbol) & symbol != "",
    unique(symbol)
  ]

  mod_file <- file.path(out_dir, paste0("module_", mod, "_symbols.txt"))
  writeLines(text = mod_symbols, con = mod_file)

  cat("[INFO] Saved module", mod, "symbol-only list:",
      length(mod_symbols), "genes ->", mod_file, "\n")
}

## save updated graph with module labels
graph_mod_file <- file.path(out_dir, "g_net_integrated_undirected_with_modules.RData")
save(g_unw, file = graph_mod_file)
cat("[INFO] Saved graph with module labels ->", graph_mod_file, "\n")

############################################################
## Part II. Visualize module-size summary
############################################################

############################################################
## Step 8) Read module summary for plotting
############################################################
plot_dt <- copy(module_summary_dt)
setnames(plot_dt, c("module_louvain_cons", "n_genes"), c("Module_ID", "N_genes"))

## If you already have curated human-readable module labels,
## you can merge them here. Otherwise, keep generic labels.
## Example generic labels:
plot_dt[, Module_label := paste0("Module ", Module_ID)]

## order by module size
setorder(plot_dt, -N_genes, Module_ID)

plot_dt[, Module_label_plot := factor(
  Module_label,
  levels = rev(Module_label)
)]

cat("[INFO] Plotting table:\n")
print(plot_dt)

############################################################
## Step 9) Define colors
############################################################
## Generic palette for module-size display.
## Replace with your final biology-based labels/colors if desired.
module_cols <- setNames(
  c(
    rgb( 79,  27, 192, maxColorValue = 255),
    rgb(105,  60, 204, maxColorValue = 255),
    rgb(130,  91, 215, maxColorValue = 255),
    rgb(  6,  69, 157, maxColorValue = 255),
    rgb( 29,  90, 173, maxColorValue = 255),
    rgb( 67, 125, 200, maxColorValue = 255),
    rgb(202,  80, 134, maxColorValue = 255),
    rgb(227, 136, 175, maxColorValue = 255)
  )[seq_len(nrow(plot_dt))],
  plot_dt$Module_label
)

############################################################
## Step 10) Plot module-size panel
############################################################
p_module_size <- ggplot(
  plot_dt,
  aes(x = N_genes, y = Module_label_plot, fill = Module_label)
) +
  geom_col(width = 0.68, color = NA) +
  geom_text(
    aes(label = N_genes),
    hjust = -0.20,
    size = 2.8,
    family = "sans"
  ) +
  scale_fill_manual(
    values = module_cols,
    breaks = plot_dt$Module_label,
    guide = "none"
  ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.12)),
    breaks = NULL
  ) +
  labs(
    title = "Module size",
    x = NULL,
    y = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(
      hjust = 0,
      size = 10,
      face = "bold",
      margin = margin(b = 6)
    ),
    axis.text.y = element_text(
      size = 8.5,
      color = "black",
      face = "plain"
    ),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.5
    ),
    plot.margin = margin(6, 18, 6, 6)
  )

print(p_module_size)

############################################################
## Step 11) Save module-size panel
############################################################
out_pdf <- file.path(out_dir, "Module.Size.Panel.v1.pdf")

ggsave(
  filename = out_pdf,
  plot     = p_module_size,
  width    = 80,
  height   = 40,
  units    = "mm"
)

cat("[INFO] Module-size panel saved to:\n")
cat("       ", out_pdf, "\n")