############################################################
## Retrieve, expand, and annotate OmniPath interaction edges
## for melanoma network construction
##
## Purpose:
##   - Download curated molecular interaction data from OmniPath
##   - Build two edge tables:
##       (1) post-translational interactions (PPI-like layer)
##       (2) transcriptional regulatory interactions (TR layer)
##   - Expand protein complexes into pairwise component edges
##   - Map source/target protein identifiers to Ensembl gene IDs
##   - Save harmonized edge tables for downstream network analysis
##
## Data sources:
##   - PPI / signaling / ligand-receptor style interactions:
##       omnipath, kinaseextra, pathwayextra, ligrecextra
##   - Transcriptional regulation:
##       collectri, dorothea, tf_target
##
## Author: Shelley
## Date:   2026-01-20
############################################################

############################################################
## Step 1) Load packages
############################################################
suppressPackageStartupMessages({
  library(OmnipathR)
  library(data.table)
})

cat("[INFO] Packages loaded: OmnipathR, data.table\n")

############################################################
## Step 2) Download interaction datasets from OmniPath
############################################################

############################
## Step 2.1) Post-translational interactions
############################
PPI_dt <- import_post_translational_interactions(
  organism    = 9606,
  datasets    = c("omnipath", "kinaseextra", "pathwayextra", "ligrecextra"),
  genesymbols = TRUE,
  fields      = c("type", "databases", "datasets", "entity_type")
)

cat("[INFO] PPI_dt loaded: ", nrow(PPI_dt), " rows x ", ncol(PPI_dt), " cols\n", sep = "")

############################
## Step 2.2) Transcriptional interactions
############################
TR_dt <- import_transcriptional_interactions(
  organism        = 9606,
  dorothea_levels = c("A", "B", "C", "D"),
  datasets        = c("collectri", "dorothea", "tf_target"),
  genesymbols     = TRUE,
  fields          = c("type", "databases", "datasets", "entity_type")
)

cat("[INFO] TR_dt loaded: ", nrow(TR_dt), " rows x ", ncol(TR_dt), " cols\n", sep = "")

############################################################
## Step 3) Expand protein complexes into component edges
##
## Rule:
##   - If source/target starts with "COMPLEX:", split its members
##   - Expand each complex edge into all pairwise combinations
##   - Optionally remove self-loops
##   - Keep all original metadata columns
############################################################

############################
## Step 3.1) Helper functions
############################
split_complex_col <- function(x) {
  x <- as.character(x)
  is_complex <- !is.na(x) & startsWith(x, "COMPLEX:")

  out <- vector("list", length(x))
  out[!is_complex] <- lapply(x[!is_complex], function(v) v)

  if (any(is_complex)) {
    stripped <- sub("^COMPLEX:", "", x[is_complex])
    out[is_complex] <- strsplit(stripped, "_", fixed = TRUE)
  }

  out
}

expand_complex_edges <- function(dt,
                                 from_col = "source",
                                 to_col = "target",
                                 drop_self_loops = TRUE,
                                 dedup_endpoints_only = FALSE) {
  stopifnot(is.data.table(dt))
  stopifnot(all(c(from_col, to_col) %in% names(dt)))

  dt2 <- copy(dt)

  dt2[, src_parts := split_complex_col(get(from_col))]
  dt2[, tgt_parts := split_complex_col(get(to_col))]

  by_cols <- setdiff(names(dt2), c(from_col, to_col, "src_parts", "tgt_parts"))

  dt_exp <- dt2[
    ,
    {
      srcs <- src_parts[[1]]
      tgts <- tgt_parts[[1]]
      CJ(src = srcs, tgt = tgts, sorted = FALSE)
    },
    by = by_cols
  ]

  setnames(dt_exp, c("src", "tgt"), c(from_col, to_col))

  if (drop_self_loops) {
    dt_exp <- dt_exp[get(from_col) != get(to_col)]
  }

  if (dedup_endpoints_only) {
    key_cols <- c(from_col, to_col)
    score_cols <- intersect(c("n_references", "n_resources", "curation_effort"), names(dt_exp))

    if (length(score_cols) == 0) {
      dt_exp <- unique(dt_exp, by = key_cols)
    } else {
      setorderv(
        dt_exp,
        cols = c(key_cols, score_cols),
        order = c(rep(1L, length(key_cols)), rep(-1L, length(score_cols)))
      )
      dt_exp <- dt_exp[!duplicated(dt_exp, by = key_cols)]
    }
  } else {
    dt_exp <- unique(dt_exp)
  }

  if (anyDuplicated(names(dt_exp))) {
    stop(
      "[ERROR] Duplicated column names remain: ",
      paste(names(dt_exp)[duplicated(names(dt_exp))], collapse = ", ")
    )
  }

  dt_exp[]
}

############################
## Step 3.2) Expand PPI and TR edges
############################
PPI_exp_dt <- expand_complex_edges(
  as.data.table(PPI_dt),
  from_col = "source",
  to_col = "target",
  drop_self_loops = TRUE,
  dedup_endpoints_only = FALSE
)

TR_exp_dt <- expand_complex_edges(
  as.data.table(TR_dt),
  from_col = "source",
  to_col = "target",
  drop_self_loops = TRUE,
  dedup_endpoints_only = FALSE
)

cat("[INFO] PPI_exp_dt: ", nrow(PPI_exp_dt), " rows x ", ncol(PPI_exp_dt), " cols\n", sep = "")
cat("[INFO] TR_exp_dt: ",  nrow(TR_exp_dt),  " rows x ", ncol(TR_exp_dt),  " cols\n", sep = "")

cat(
  "[INFO] Remaining COMPLEX entries in PPI_exp_dt: ",
  sum(grepl("^COMPLEX:", PPI_exp_dt$source) | grepl("^COMPLEX:", PPI_exp_dt$target), na.rm = TRUE),
  "\n", sep = ""
)
cat(
  "[INFO] Remaining COMPLEX entries in TR_exp_dt: ",
  sum(grepl("^COMPLEX:", TR_exp_dt$source) | grepl("^COMPLEX:", TR_exp_dt$target), na.rm = TRUE),
  "\n", sep = ""
)

############################################################
## Step 4) Map source and target identifiers to Ensembl IDs
##
## Input:
##   - source / target columns from expanded OmniPath tables
##
## Output:
##   - source_ensg
##   - target_ensg
##
## Strategy:
##   - Use OmnipathR::ensembl_id_mapping_table("ensg")
##   - Collapse to one UniProt -> one ENSG mapping
##   - If multiple ENSG IDs exist for one UniProt ID,
##     keep the most frequent mapping
############################################################

############################
## Step 4.1) Read and prepare UniProt -> ENSG mapping
############################
ensg_up <- ensembl_id_mapping_table("ensg")
ensg_up_dt <- as.data.table(ensg_up)

setnames(ensg_up_dt, c("From", "To"), c("uniprot", "ensg"))
ensg_up_dt <- ensg_up_dt[
  !is.na(uniprot) & uniprot != "" &
    !is.na(ensg) & ensg != ""
]

up2ensg_1to1 <- ensg_up_dt[
  ,
  .(ensg = names(sort(table(ensg), decreasing = TRUE))[1]),
  by = uniprot
]

cat("[INFO] up2ensg_1to1 rows: ", nrow(up2ensg_1to1), "\n", sep = "")

############################
## Step 4.2) Helper to annotate source and target
############################
annotate_ensg <- function(dt, map_dt) {
  dt <- as.data.table(copy(dt))
  stopifnot(all(c("source", "target") %in% names(dt)))

  ## source
  dt <- merge(
    dt,
    map_dt,
    by.x = "source",
    by.y = "uniprot",
    all.x = TRUE
  )
  setnames(dt, "ensg", "source_ensg")

  ## target
  map2 <- copy(map_dt)
  setnames(map2, c("uniprot", "ensg"), c("uniprot_target", "target_ensg"))

  dt <- merge(
    dt,
    map2,
    by.x = "target",
    by.y = "uniprot_target",
    all.x = TRUE
  )

  dt[]
}

############################
## Step 4.3) Annotate both edge tables
############################
PPI_exp_dt <- annotate_ensg(PPI_exp_dt, up2ensg_1to1)
TR_exp_dt  <- annotate_ensg(TR_exp_dt,  up2ensg_1to1)

cat("[INFO] PPI missing source_ensg: ", sum(is.na(PPI_exp_dt$source_ensg)), "\n", sep = "")
cat("[INFO] PPI missing target_ensg: ", sum(is.na(PPI_exp_dt$target_ensg)), "\n", sep = "")
cat("[INFO] TR missing source_ensg: ",  sum(is.na(TR_exp_dt$source_ensg)),  "\n", sep = "")
cat("[INFO] TR missing target_ensg: ",  sum(is.na(TR_exp_dt$target_ensg)),  "\n", sep = "")

############################################################
## Step 5) Save harmonized edge tables
############################################################
out_dir <- "/Users/shelleyz/Projects/Melanoma-WGS/Omnigenic.Architecture/Edges"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ppi_out <- file.path(out_dir, "OmniPath_PPI_expanded_with_ENSG.txt")
tr_out  <- file.path(out_dir, "OmniPath_TR_expanded_with_ENSG.txt")

fwrite(PPI_exp_dt, ppi_out, sep = "\t", quote = FALSE, na = "NA")
fwrite(TR_exp_dt,  tr_out,  sep = "\t", quote = FALSE, na = "NA")

cat("[INFO] Saved PPI_exp_dt -> ", ppi_out, "\n", sep = "")
cat("[INFO] Saved TR_exp_dt  -> ", tr_out,  "\n", sep = "")
cat("[INFO] Finished.\n")