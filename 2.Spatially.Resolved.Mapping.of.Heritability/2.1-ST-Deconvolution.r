############################################################
## Spatial transcriptomics deconvolution (SpaCET pipeline)
## Reference: https://github.com/data2intelligence/SpaCET
## Author: Shiyu Zhang
## Last Updated:   2026-02-10
############################################################

############################################################
## Step 1) Load packages
############################################################
suppressPackageStartupMessages({
  library(SpaCET)
  library(data.table)
  library(Matrix)
})

############################################################
## Step 2) Read ST data and create SpaCET object
############################################################
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript run_spacet.R <sample_ST_ID>\n")
  cat("Example: Rscript run_spacet.R MEL162\n")
  quit(save = "no", status = 1)
}

sample_ST_ID <- args[1]

base_dir <- "/public/home/hpc8301200407/MEL.ST.10x.Visium"
visiumPath <- file.path(base_dir, sample_ST_ID)

## output directory
BASE_OUTDIR <- "/public/home/hpc8301200407/WGS/GWAS_ST/SpaCET_process_updated_20260210"
outdir <- file.path(BASE_OUTDIR, sample_ST_ID)

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(outdir, "tables"), showWarnings = FALSE)
dir.create(file.path(outdir, "rds"), showWarnings = FALSE)
dir.create(file.path(outdir, "log"), showWarnings = FALSE)
dir.create(file.path(outdir, "fig"), showWarnings = FALSE)

## sanity checks
if (!dir.exists(visiumPath)) {
  stop("[ERROR] visiumPath not found: ", visiumPath)
}

cat("[INFO] sample_ST_ID: ", sample_ST_ID, "\n", sep = "")
cat("[INFO] visiumPath:   ", visiumPath, "\n", sep = "")
cat("[INFO] outdir:       ", outdir, "\n", sep = "")

## create SpaCET object from 10X Visium folder
SpaCET_obj <- create.SpaCET.object.10X(visiumPath = visiumPath)

## robust accessors
stopifnot("input" %in% slotNames(SpaCET_obj))
stopifnot(all(c("counts", "spotCoordinates", "metaData") %in% names(SpaCET_obj@input)))

counts <- SpaCET_obj@input$counts
coords <- SpaCET_obj@input$spotCoordinates
meta   <- SpaCET_obj@input$metaData

stopifnot(inherits(counts, "dgCMatrix"))
stopifnot(ncol(counts) == nrow(coords))

cat("[INFO] SpaCET object created for sample: ", sample_ST_ID, "\n", sep = "")
cat("[INFO] Genes x spots: ", nrow(counts), " x ", ncol(counts), "\n", sep = "")

############################################################
## Step 3) Extract per-spot QC metrics
############################################################
spot_ids <- colnames(counts)

## UMI per spot
umi_vec <- Matrix::colSums(counts)

## nGene per spot
ngene_vec <- Matrix::colSums(counts > 0)

qc_dt <- data.table(
  spot_id   = spot_ids,
  barcode   = if ("barcode" %in% colnames(meta)) meta$barcode else NA_character_,
  UMI       = as.numeric(umi_vec),
  nGene     = as.numeric(ngene_vec),
  log1p_UMI = log1p(as.numeric(umi_vec))
)

## add coordinates
coord_dt <- data.table::as.data.table(coords, keep.rownames = FALSE)
if (nrow(coord_dt) == nrow(qc_dt)) {
  qc_dt <- cbind(qc_dt, coord_dt)
}

qc_outfile <- file.path(outdir, "tables", paste0(sample_ST_ID, ".spot_QC_metrics.tsv"))
data.table::fwrite(qc_dt, qc_outfile, sep = "\t")

cat("[INFO] Saved spot QC metrics: ", qc_outfile, "\n", sep = "")

############################################################
## Step 4) QC filtering and visualization
############################################################
SpaCET_obj <- SpaCET.quality.control(SpaCET_obj, min.genes = 100)
cat("[INFO] QC filtering done (min.genes = 100).\n")

fig.QC <- SpaCET.visualize.spatialFeature(
  SpaCET_obj,
  spatialType = "QualityControl",
  spatialFeatures = c("UMI", "Gene"),
  imageBg = TRUE,
  imageSize = "CaptureArea",
  pointSize = 0.4
)

qc_pdf <- file.path(outdir, "fig", paste0(sample_ST_ID, ".QC_UMI_Gene.minGenes100.pdf"))

pdf(qc_pdf, width = 7.2, height = 2.4, useDingbats = FALSE)
print(fig.QC)
dev.off()

cat("[INFO] Saved QC PDF: ", qc_pdf, "\n", sep = "")

############################################################
## Step 5) Deconvolution and output export
############################################################
SpaCET_obj <- SpaCET.deconvolution(SpaCET_obj, cancerType = "SKCM", coreNo = 8)
cat("[INFO] SpaCET deconvolution finished.\n")

stopifnot("results" %in% slotNames(SpaCET_obj))
stopifnot("deconvolution" %in% names(SpaCET_obj@results))

deconv <- SpaCET_obj@results$deconvolution
stopifnot("propMat" %in% names(deconv))

propMat <- deconv$propMat

stopifnot(is.matrix(propMat))
stopifnot(!is.null(rownames(propMat)))
stopifnot(!is.null(colnames(propMat)))

cat("[INFO] propMat dim (cell types x spots): ",
    nrow(propMat), " x ", ncol(propMat), "\n", sep = "")

## save cell-type proportion matrix
prop_dt <- data.table::as.data.table(t(propMat), keep.rownames = "spot_id")
prop_outfile <- file.path(outdir, "tables", paste0(sample_ST_ID, ".SpaCET.celltype_prop.tsv"))
data.table::fwrite(prop_dt, prop_outfile, sep = "\t")
cat("[INFO] Saved cell-type proportions: ", prop_outfile, "\n", sep = "")

## optional long format
prop_long_dt <- melt(
  prop_dt,
  id.vars = "spot_id",
  variable.name = "cell_type",
  value.name = "proportion"
)
prop_long_outfile <- file.path(outdir, "tables", paste0(sample_ST_ID, ".SpaCET.celltype_prop.long.tsv"))
data.table::fwrite(prop_long_dt, prop_long_outfile, sep = "\t")
cat("[INFO] Saved cell-type proportions (long): ", prop_long_outfile, "\n", sep = "")

## save signature genes
sigGenes <- NULL
if ("Ref" %in% names(deconv) && "sigGenes" %in% names(deconv$Ref)) {
  sigGenes <- deconv$Ref$sigGenes
}

if (!is.null(sigGenes)) {
  stopifnot(is.list(sigGenes), !is.null(names(sigGenes)))

  sig_long_dt <- rbindlist(
    lapply(names(sigGenes), function(ct) {
      data.table(cell_type = ct, gene = as.character(sigGenes[[ct]]))
    }),
    use.names = TRUE,
    fill = TRUE
  )

  sig_outfile <- file.path(outdir, "tables", paste0(sample_ST_ID, ".SpaCET.signature_genes.tsv"))
  data.table::fwrite(sig_long_dt, sig_outfile, sep = "\t")
  cat("[INFO] Saved signature genes: ", sig_outfile, "\n", sep = "")
} else {
  cat("[WARN] Signature genes not found at SpaCET_obj@results$deconvolution$Ref$sigGenes\n")
}

############################################################
## Step 6) Cell-cell interaction colocalization
############################################################
SpaCET_obj <- SpaCET.CCI.colocalization(SpaCET_obj)
cat("[INFO] SpaCET CCI colocalization finished.\n")

stopifnot("results" %in% slotNames(SpaCET_obj))
stopifnot("CCI" %in% names(SpaCET_obj@results))
stopifnot("colocalization" %in% names(SpaCET_obj@results$CCI))

coloc_df <- SpaCET_obj@results$CCI$colocalization
stopifnot(is.data.frame(coloc_df))

coloc_dt <- data.table::as.data.table(coloc_df)

print(coloc_dt[fraction_rho > 0 & fraction_pv < 0.05])

expected_cols <- c(
  "cell_type_1", "cell_type_2",
  "fraction_product", "fraction_rho", "fraction_pv",
  "reference_rho", "reference_pv"
)
missing_cols <- setdiff(expected_cols, colnames(coloc_dt))
if (length(missing_cols) > 0) {
  cat("[WARN] Missing expected columns in colocalization table: ",
      paste(missing_cols, collapse = ","), "\n", sep = "")
}

coloc_outfile <- file.path(outdir, "tables", paste0(sample_ST_ID, ".SpaCET.CCI.colocalization.tsv"))
data.table::fwrite(coloc_dt, coloc_outfile, sep = "\t")
cat("[INFO] Saved CCI colocalization table: ", coloc_outfile, "\n", sep = "")

############################################################
## Step 7) Ligand-receptor network score per spot
############################################################
SpaCET_obj <- SpaCET.CCI.LRNetworkScore(SpaCET_obj, coreNo = 8)
cat("[INFO] SpaCET CCI LRNetworkScore finished.\n")

stopifnot("results" %in% slotNames(SpaCET_obj))
stopifnot("CCI" %in% names(SpaCET_obj@results))
stopifnot("LRNetworkScore" %in% names(SpaCET_obj@results$CCI))

lr_mat <- SpaCET_obj@results$CCI$LRNetworkScore
stopifnot(is.matrix(lr_mat))
stopifnot(!is.null(rownames(lr_mat)))
stopifnot(!is.null(colnames(lr_mat)))

cat("[INFO] LRNetworkScore dim (metrics x spots): ",
    nrow(lr_mat), " x ", ncol(lr_mat), "\n", sep = "")
cat("[INFO] LRNetworkScore rows: ", paste(rownames(lr_mat), collapse = ","), "\n", sep = "")

lr_dt <- data.table(
  spot_id = colnames(lr_mat),
  Raw_expr = as.numeric(lr_mat["Raw_expr", ]),
  Network_Score = as.numeric(lr_mat["Network_Score", ]),
  Network_Score_pv = as.numeric(lr_mat["Network_Score_pv", ])
)

lr_outfile <- file.path(outdir, "tables", paste0(sample_ST_ID, ".SpaCET.CCI.LRNetworkScore.tsv"))
data.table::fwrite(lr_dt, lr_outfile, sep = "\t")
cat("[INFO] Saved LRNetworkScore table: ", lr_outfile, "\n", sep = "")

############################################################
## Step 8) Tumor-stroma interface calling
############################################################
SpaCET_obj <- SpaCET.identify.interface(SpaCET_obj)
cat("[INFO] SpaCET identify.interface finished.\n")

stopifnot("results" %in% slotNames(SpaCET_obj))
stopifnot("CCI" %in% names(SpaCET_obj@results))
stopifnot("interface" %in% names(SpaCET_obj@results$CCI))

iface_mat <- SpaCET_obj@results$CCI$interface
stopifnot(is.matrix(iface_mat))
stopifnot(nrow(iface_mat) == 1)

iface_dt <- data.table(
  spot_id = colnames(iface_mat),
  interface_label = as.character(iface_mat[1, ])
)

iface_tab <- iface_dt[, .N, by = interface_label][order(-N)]
print(iface_tab)

iface_outfile <- file.path(outdir, "tables", paste0(sample_ST_ID, ".SpaCET.interface_labels.tsv"))
data.table::fwrite(iface_dt, iface_outfile, sep = "\t")
cat("[INFO] Saved interface labels: ", iface_outfile, "\n", sep = "")

############################################################
## Step 9) Gene set scoring per spot
############################################################
save_geneset_score <- function(SpaCET_obj, geneset_name, outdir, sample_ST_ID) {
  stopifnot("results" %in% slotNames(SpaCET_obj))
  stopifnot("GeneSetScore" %in% names(SpaCET_obj@results))
  stopifnot(is.matrix(SpaCET_obj@results$GeneSetScore))

  gs_mat <- SpaCET_obj@results$GeneSetScore
  stopifnot(!is.null(rownames(gs_mat)), !is.null(colnames(gs_mat)))

  gs_dt <- data.table::as.data.table(t(gs_mat), keep.rownames = "spot_id")

  outfile <- file.path(
    outdir,
    "tables",
    paste0(sample_ST_ID, ".SpaCET.GeneSetScore.", geneset_name, ".tsv")
  )
  data.table::fwrite(gs_dt, outfile, sep = "\t")

  cat("[INFO] Saved GeneSetScore (", geneset_name, "): ", outfile, "\n", sep = "")
  invisible(outfile)
}

## 9.1 Hallmark
SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets = "Hallmark")
cat("[INFO] GeneSetScore finished: Hallmark\n")
save_geneset_score(SpaCET_obj, geneset_name = "Hallmark", outdir = outdir, sample_ST_ID = sample_ST_ID)

## 9.2 CancerCellState
SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets = "CancerCellState")
cat("[INFO] GeneSetScore finished: CancerCellState\n")
save_geneset_score(SpaCET_obj, geneset_name = "CancerCellState", outdir = outdir, sample_ST_ID = sample_ST_ID)

## 9.3 TLS
SpaCET_obj <- SpaCET.GeneSetScore(SpaCET_obj, GeneSets = "TLS")
cat("[INFO] GeneSetScore finished: TLS\n")
save_geneset_score(SpaCET_obj, geneset_name = "TLS", outdir = outdir, sample_ST_ID = sample_ST_ID)

############################################################
## Step 10) Compute spatial weight matrix
############################################################
W <- calWeights(SpaCET_obj, radius = 200, sigma = 100, diagAsZero = TRUE)

############################################################
## Step 11) Spatially variable genes
############################################################
suppressPackageStartupMessages(library(future))
options(future.globals.maxSize = 20 * 1024^3)

SpaCET_obj <- SpaCET.SpatialCorrelation(
  SpaCET_obj,
  mode = "univariate",
  item = NULL,
  W = W,
  nPermutation = 1000
)

stopifnot("results" %in% slotNames(SpaCET_obj))
stopifnot("SpatialCorrelation" %in% names(SpaCET_obj@results))
stopifnot("univariate" %in% names(SpaCET_obj@results$SpatialCorrelation))

svg_df <- SpaCET_obj@results$SpatialCorrelation$univariate
stopifnot(is.data.frame(svg_df))
stopifnot(all(c("p.Moran_I", "p.Moran_Z", "p.Moran_P", "p.Moran_Padj") %in% colnames(svg_df)))

if (!is.null(rownames(svg_df))) {
  svg_dt <- data.table::as.data.table(svg_df, keep.rownames = "gene")
} else if ("gene" %in% colnames(svg_df)) {
  svg_dt <- data.table::as.data.table(svg_df)
} else {
  stop("[ERROR] Cannot find gene identifiers.")
}

setorder(svg_dt, p.Moran_Padj, -p.Moran_I)
svg_dt[, is_SVG := (p.Moran_Padj < 0.05)]

cat("[INFO] SVG table rows: ", nrow(svg_dt), "\n", sep = "")
cat("[INFO] Significant SVGs (FDR < 0.05): ", svg_dt[is_SVG == TRUE, .N], "\n", sep = "")

print(svg_dt[1:min(20, .N)])

svg_outfile <- file.path(outdir, "tables", paste0(sample_ST_ID, ".SpaCET.SVG.univariate.Moran.tsv"))
data.table::fwrite(svg_dt, svg_outfile, sep = "\t")
cat("[INFO] Saved SVG table: ", svg_outfile, "\n", sep = "")

############################################################
## Step 12) Spatially co-expressed ligand-receptor pairs
############################################################
Sys.time()
SpaCET_obj <- SpaCET.SpatialCorrelation(
  SpaCET_obj,
  mode = "bivariate",
  item = NULL,
  W = W,
  nPermutation = 1000
)
Sys.time()

co_expression <- SpaCET_obj@results$SpatialCorrelation$bivariate
stopifnot(is.data.frame(co_expression))

co_expression <- data.table::as.data.table(co_expression, keep.rownames = "pair_id")
co_expression[, c("ligand", "receptor") := tstrsplit(pair_id, "_", fixed = TRUE)]

setorder(co_expression, p.Moran_Padj, -p.Moran_I)
co_expression[, is_sig := (p.Moran_Padj < 0.05)]

print(co_expression[1:20])

outfile_all <- file.path(
  outdir,
  "tables",
  paste0(sample_ST_ID, ".SpaCET.co_expression.spatial_bivariate.tsv")
)
data.table::fwrite(co_expression, outfile_all, sep = "\t")

############################################################
## Final step) Save SpaCET object
############################################################
dir.create(file.path(outdir, "rds"), showWarnings = FALSE, recursive = TRUE)

saveRDS(
  SpaCET_obj,
  file = file.path(outdir, "rds", paste0(sample_ST_ID, ".SpaCET_obj.final.rds")),
  compress = "xz"
)

cat("[INFO] SpaCET object saved as RDS.\n")