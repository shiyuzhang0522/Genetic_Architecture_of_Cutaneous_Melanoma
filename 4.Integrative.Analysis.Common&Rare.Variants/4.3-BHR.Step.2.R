############################################################
## BHR - Step.2 Format required input for BHR analysis
##
## Author: Shelley
## Date:  2025-12-16
############################################################

############################################################
## Step 1) Load packages
############################################################
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(biomaRt)
})

############################################################
## Step 2) Define inputs and helper functions
############################################################
mask_dir   <- "/path/WGS/Melanoma_WGS/BHR/Coding.GT"
counts_dir <- "/path/WGS/Melanoma_WGS/BHR/Coding.GT/UKB/"
out_dir    <- "/path/WGS/Melanoma_WGS/BHR/pipeline/input"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

chrs <- 1:22

stop_if_missing <- function(x) {
  miss <- x[!file.exists(x)]
  if (length(miss) > 0) {
    stop("Missing files:\n", paste(miss, collapse = "\n"))
  }
  invisible(TRUE)
}

extract_chr_from_file <- function(x) {
  as.integer(str_match(basename(x), "^chr(\\d+)\\.")[[1, 2]])
}

############################################################
## Step 3) Read mask and allele-count files
############################################################

############################
## Step 3.1) Mask annotation files
############################
mask_files <- file.path(mask_dir, sprintf("chr%d.coding_only.long.csv", chrs))
stop_if_missing(mask_files)

mask_list <- lapply(mask_files, function(f) {
  dt <- fread(f)
  dt[, chr_file := extract_chr_from_file(f)]
  dt
})
names(mask_list) <- paste0("chr", chrs)

cat("[INFO] Loaded mask files: ", length(mask_list), " chromosomes\n", sep = "")
cat("[INFO] Example mask table (chr1): ", nrow(mask_list[[1]]), " rows x ", ncol(mask_list[[1]]), " cols\n", sep = "")

############################
## Step 3.2) Case ALT-count files
############################
case_files <- file.path(counts_dir, sprintf("chr%d.case.ALT_counts.acount", chrs))
stop_if_missing(case_files)

case_counts_list <- lapply(case_files, function(f) {
  dt <- fread(f)
  dt[, chr_file := extract_chr_from_file(f)]
  dt
})
names(case_counts_list) <- paste0("chr", chrs)

cat("[INFO] Loaded case .acount files: ", length(case_counts_list), " chromosomes\n", sep = "")
cat("[INFO] Example case .acount (chr1): ", nrow(case_counts_list[[1]]), " rows x ", ncol(case_counts_list[[1]]), " cols\n", sep = "")

############################
## Step 3.3) Control ALT-count files
############################
ctrl_files <- file.path(counts_dir, sprintf("chr%d.control.ALT_counts.acount", chrs))
stop_if_missing(ctrl_files)

control_counts_list <- lapply(ctrl_files, function(f) {
  dt <- fread(f)
  dt[, chr_file := extract_chr_from_file(f)]
  dt
})
names(control_counts_list) <- paste0("chr", chrs)

cat("[INFO] Loaded control .acount files: ", length(control_counts_list), " chromosomes\n", sep = "")
cat("[INFO] Example control .acount (chr1): ", nrow(control_counts_list[[1]]), " rows x ", ncol(control_counts_list[[1]]), " cols\n", sep = "")

cat("[INFO] Columns in case .acount:\n")
print(names(case_counts_list[[1]]))

############################################################
## Step 4) Merge case and control allele counts per chromosome
############################################################
merge_case_control_counts_one_chr <- function(case_dt, ctrl_dt) {
  case_dt <- as.data.table(case_dt)
  ctrl_dt <- as.data.table(ctrl_dt)

  case_dt2 <- case_dt[, .(
    chr = `#CHROM`,
    ID,
    REF,
    ALT,
    ALT_CTS_case = ALT_CTS,
    OBS_CT_case  = OBS_CT
  )]

  ctrl_dt2 <- ctrl_dt[, .(
    chr = `#CHROM`,
    ID,
    REF,
    ALT,
    ALT_CTS_control = ALT_CTS,
    OBS_CT_control  = OBS_CT
  )]

  merge(
    case_dt2,
    ctrl_dt2,
    by = c("chr", "ID", "REF", "ALT"),
    all = TRUE
  )
}

counts_merged_list <- vector("list", length(chrs))
names(counts_merged_list) <- paste0("chr", chrs)

for (chr in chrs) {
  key <- paste0("chr", chr)
  cat("[INFO] Merging case/control counts for ", key, " ...\n", sep = "")

  counts_merged_list[[key]] <- merge_case_control_counts_one_chr(
    case_dt = case_counts_list[[key]],
    ctrl_dt = control_counts_list[[key]]
  )

  cat("[INFO] ", key, ": merged rows = ", nrow(counts_merged_list[[key]]), "\n", sep = "")
}

############################################################
## Step 5) Merge counts with mask annotations per chromosome
############################################################
merge_mask_with_counts_one_chr <- function(mask_dt, counts_dt) {
  mask_dt   <- as.data.table(mask_dt)
  counts_dt <- as.data.table(counts_dt)

  mask_dt2 <- mask_dt[, .(
    SNP_ID_final,
    SNP_ID_merge,
    GENE,
    LOF,
    REVEL_SCORE,
    CADD_PHRED,
    CSQ,
    SpliceAI_max,
    annotation,
    chr_file
  )]

  counts_dt2 <- counts_dt[, .(
    chr,
    ID,
    REF,
    ALT,
    ALT_CTS_case,
    OBS_CT_case,
    ALT_CTS_control,
    OBS_CT_control
  )]

  merged <- merge(
    mask_dt2,
    counts_dt2,
    by.x = "SNP_ID_final",
    by.y = "ID",
    all.x = TRUE
  )

  merged[, missing_counts := is.na(ALT_CTS_case) & is.na(ALT_CTS_control)]
  n_miss <- merged[, sum(missing_counts)]
  if (n_miss > 0) {
    cat("[WARN] Missing counts for ", n_miss, " variants\n", sep = "")
  }

  merged
}

mask_counts_merged_list <- vector("list", length(chrs))
names(mask_counts_merged_list) <- paste0("chr", chrs)

for (chr in chrs) {
  key <- paste0("chr", chr)
  cat("[INFO] Merging mask and counts for ", key, " ...\n", sep = "")

  mask_counts_merged_list[[key]] <- merge_mask_with_counts_one_chr(
    mask_dt   = mask_list[[key]],
    counts_dt = counts_merged_list[[key]]
  )

  cat("[INFO] ", key, ": rows = ", nrow(mask_counts_merged_list[[key]]),
      "; cols = ", ncol(mask_counts_merged_list[[key]]), "\n", sep = "")
}

############################################################
## Step 6) Retrieve gene TSS from Ensembl via biomaRt
############################################################

############################
## Step 6.1) Collect unique genes
############################
all_genes <- unique(unlist(lapply(mask_counts_merged_list, function(dt) unique(dt$GENE))))
all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
all_genes_clean <- unique(gsub("\\..*$", "", all_genes))

cat("[INFO] Unique genes (raw): ", length(all_genes), "\n", sep = "")
cat("[INFO] Unique genes (clean): ", length(all_genes_clean), "\n", sep = "")

############################
## Step 6.2) Connect to Ensembl
############################
connect_ensembl <- function() {
  mirrors <- c("www", "useast", "uswest", "asia")
  for (m in mirrors) {
    cat("[INFO] Trying Ensembl mirror: ", m, "\n", sep = "")
    mart <- tryCatch(
      useEnsembl(
        biomart = "genes",
        dataset = "hsapiens_gene_ensembl",
        mirror = m
      ),
      error = function(e) NULL
    )
    if (!is.null(mart)) return(mart)
  }
  stop("Failed to connect to Ensembl via all mirrors.")
}

mart <- connect_ensembl()

############################
## Step 6.3) Query gene coordinates
############################
attrs <- c(
  "ensembl_gene_id",
  "chromosome_name",
  "start_position",
  "end_position",
  "strand"
)

gene_pos <- getBM(
  attributes = attrs,
  filters    = "ensembl_gene_id",
  values     = all_genes_clean,
  mart       = mart
)

gene_pos_dt <- as.data.table(gene_pos)

cat("[INFO] Rows returned from biomaRt: ", nrow(gene_pos_dt), "\n", sep = "")
cat("[INFO] Unique genes returned: ", uniqueN(gene_pos_dt$ensembl_gene_id), "\n", sep = "")

############################
## Step 6.4) Compute TSS
############################
gene_pos_dt[, gene_tss := fifelse(strand == 1, start_position, end_position)]
gene_pos_dt[, gene_chr := chromosome_name]

gene_tss_dt <- unique(gene_pos_dt[, .(
  GENE = ensembl_gene_id,
  gene_chr,
  gene_start = start_position,
  gene_end   = end_position,
  gene_strand = strand,
  gene_tss
)])

not_found <- setdiff(all_genes_clean, gene_tss_dt$GENE)
cat("[INFO] Genes not found in biomaRt: ", length(not_found), "\n", sep = "")
if (length(not_found) > 0) {
  print(head(not_found, 20))
}

############################################################
## Step 7) Merge TSS back into per-chromosome tables
############################################################
mask_counts_with_tss_list <- vector("list", length(chrs))
names(mask_counts_with_tss_list) <- paste0("chr", chrs)

for (chr in chrs) {
  key <- paste0("chr", chr)
  dt <- as.data.table(mask_counts_merged_list[[key]])

  dt[, GENE_clean := gsub("\\..*$", "", GENE)]

  dt2 <- merge(
    dt,
    gene_tss_dt,
    by.x = "GENE_clean",
    by.y = "GENE",
    all.x = TRUE
  )

  n_miss_tss <- dt2[, sum(is.na(gene_tss))]
  cat("[INFO] ", key, ": variants = ", nrow(dt2),
      "; missing gene_tss = ", n_miss_tss, "\n", sep = "")

  mask_counts_with_tss_list[[key]] <- dt2
}

############################################################
## Step 8) Combine all chromosomes into one table
############################################################
mask_counts_with_tss_list <- lapply(mask_counts_with_tss_list, as.data.table)

mask_counts_with_tss_dt <- rbindlist(
  mask_counts_with_tss_list,
  use.names = TRUE,
  fill = TRUE,
  idcol = "chr_key"
)

cat("[INFO] Combined rows: ", nrow(mask_counts_with_tss_dt), "\n", sep = "")
cat("[INFO] Combined cols: ", ncol(mask_counts_with_tss_dt), "\n", sep = "")

save(mask_counts_with_tss_dt, file = file.path(out_dir, "mask_counts_with_tss_dt.RData"))
cat("[INFO] Saved RData to: ", file.path(out_dir, "mask_counts_with_tss_dt.RData"), "\n", sep = "")

############################################################
## Step 9) Compute ALT allele frequencies
############################################################
dt <- as.data.table(mask_counts_with_tss_dt)

n_na_case    <- dt[, sum(is.na(ALT_CTS_case))]
n_na_control <- dt[, sum(is.na(ALT_CTS_control))]
n_na_any     <- dt[, sum(is.na(ALT_CTS_case) | is.na(ALT_CTS_control))]

cat("[INFO] NA ALT_CTS_case       : ", n_na_case, "\n", sep = "")
cat("[INFO] NA ALT_CTS_control    : ", n_na_control, "\n", sep = "")
cat("[INFO] NA in case OR control : ", n_na_any, "\n", sep = "")

## NOTE:
## Missing ALT counts are not forcibly converted to zero here.
## Uncomment below only if absent variants should be treated as zero-count.
# dt[is.na(ALT_CTS_case),    ALT_CTS_case := 0L]
# dt[is.na(ALT_CTS_control), ALT_CTS_control := 0L]

dt[, ALT_AF_case    := fifelse(OBS_CT_case    > 0, ALT_CTS_case    / OBS_CT_case,    NA_real_)]
dt[, ALT_AF_control := fifelse(OBS_CT_control > 0, ALT_CTS_control / OBS_CT_control, NA_real_)]

dt[, ALT_CTS_total := ALT_CTS_case + ALT_CTS_control]
dt[, OBS_CT_total  := OBS_CT_case  + OBS_CT_control]
dt[, ALT_AF_total  := fifelse(OBS_CT_total > 0, ALT_CTS_total / OBS_CT_total, NA_real_)]

dt_filtered <- dt[!(ALT_CTS_case == 0L & ALT_CTS_control == 0L)]

cat("[INFO] Original variants: ", nrow(dt), "\n", sep = "")
cat("[INFO] Kept variants: ", nrow(dt_filtered), "\n", sep = "")
cat("[INFO] Dropped variants: ", nrow(dt) - nrow(dt_filtered), "\n", sep = "")

############################################################
## Step 10) Compute per-SD and per-allele beta
############################################################
n_cases    <- 3816
n_controls <- 339903
prevalence <- n_cases / (n_cases + n_controls)

cat("[INFO] n_cases = ", n_cases, "\n", sep = "")
cat("[INFO] n_controls = ", n_controls, "\n", sep = "")
cat("[INFO] prevalence = ", prevalence, "\n", sep = "")

table <- as.data.table(dt_filtered)

table[, AF_case := ALT_AF_case]
table[, AF      := ALT_AF_total]

table[, beta := (2 * (AF_case - AF) * prevalence) /
              sqrt(2 * AF * (1 - AF) * prevalence * (1 - prevalence))]

table[, twopq := 2 * AF * (1 - AF)]
table[, beta_perallele := beta / sqrt(twopq)]

bad <- !is.finite(table$beta) | !is.finite(table$beta_perallele)
cat("[INFO] Non-finite beta rows: ", sum(bad), "\n", sep = "")

table[bad, c("beta", "twopq", "beta_perallele") := .(NA_real_, NA_real_, NA_real_)]

############################################################
## Step 11) Build and write the full BHR input table
############################################################
N_total <- n_cases + n_controls

bhr_dt <- table[, .(
  gene            = GENE_clean,
  chromosome      = gene_chr,
  variant_id      = SNP_ID_final,
  consequence     = CSQ,
  mask            = annotation,
  gene_position   = gene_tss,
  twopq           = twopq,
  beta            = beta_perallele,
  AF              = AF,
  AF_cases        = AF_case,
  ALT_CTS_case    = ALT_CTS_case,
  OBS_CT_case     = OBS_CT_case,
  ALT_CTS_control = ALT_CTS_control,
  OBS_CT_control  = OBS_CT_control,
  ALT_AF_case     = ALT_AF_case,
  ALT_AF_control  = ALT_AF_control,
  ALT_CTS_total   = ALT_CTS_total,
  OBS_CT_total    = OBS_CT_total,
  ALT_AF_total    = ALT_AF_total,
  N               = N_total,
  phenotype_key   = "melanoma"
)]

bhr_out_file <- file.path(out_dir, "bhr_input_melanoma.txt")

fwrite(
  bhr_dt,
  file = bhr_out_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("[INFO] BHR input written to:\n", bhr_out_file, "\n", sep = "")

############################################################
## Step 12) Restrict to selected mask groups and filter by call rate
############################################################
dt <- as.data.table(bhr_dt)

keep_masks <- c("pLoF", "damaging_missense_or_protein_altering", "synonymous")
dt <- dt[mask %in% keep_masks]

dt[, MAC_total := ALT_CTS_case + ALT_CTS_control]
dt[, ALT_AF := AF]

obs_max <- 2 * N_total
dt[, callrate_total := OBS_CT_total / obs_max]

n_low_call <- dt[callrate_total < 0.95, .N]
cat("[INFO] Variants with callrate_total < 0.95: ", n_low_call, "\n", sep = "")

cat("[INFO] callrate_total summary:\n")
print(dt[, .(
  min = min(callrate_total, na.rm = TRUE),
  mean = mean(callrate_total, na.rm = TRUE),
  median = median(callrate_total, na.rm = TRUE),
  max = max(callrate_total, na.rm = TRUE)
)])

dt_call95 <- dt[callrate_total >= 0.95]

cat("[INFO] Variants kept after callrate >= 0.95: ", nrow(dt_call95), "\n", sep = "")
cat("[INFO] Variants removed: ", nrow(dt) - nrow(dt_call95), "\n", sep = "")

############################################################
## Step 13) Define frequency bins
############################################################
dt_call95[, MAC_total := ALT_CTS_case + ALT_CTS_control]
dt_call95[, ALT_AF := AF]

dt_call95[, freq_bin := fifelse(
  MAC_total <= 10, "ultra_rare_MAC_le_10",
  fifelse(
    ALT_AF <= 1e-4, "MAC_gt_10_ALT_AF_le_1e-4",
    fifelse(
      ALT_AF > 1e-4 & ALT_AF <= 1e-3, "ALT_AF_gt_1e-4_to_le_1e-3",
      NA_character_
    )
  )
)]

cat("[INFO] freq_bin summary after callrate >= 0.95:\n")
print(
  dt_call95[!is.na(freq_bin),
            .(
              n = .N,
              ALT_AF_min  = min(ALT_AF, na.rm = TRUE),
              ALT_AF_mean = mean(ALT_AF, na.rm = TRUE),
              ALT_AF_max  = max(ALT_AF, na.rm = TRUE),
              MAC_min     = min(MAC_total, na.rm = TRUE),
              MAC_mean    = mean(MAC_total, na.rm = TRUE),
              MAC_max     = max(MAC_total, na.rm = TRUE)
            ),
            by = freq_bin][order(freq_bin)]
)

dt_call95_file <- file.path(out_dir, "dt_call95.txt")
fwrite(
  dt_call95,
  file = dt_call95_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("[INFO] dt_call95 written to:\n", dt_call95_file, "\n", sep = "")

############################################################
## Step 14) Write per-frequency-bin folders split by mask
############################################################
dt_runs <- dt_call95[!is.na(freq_bin) & mask %in% keep_masks]

base_out <- file.path(out_dir, "runs_by_freqbin")
dir.create(base_out, recursive = TRUE, showWarnings = FALSE)

cat("[INFO] Variants per freq_bin × mask:\n")
print(dt_runs[, .N, by = .(freq_bin, mask)][order(freq_bin, mask)])

for (fb in sort(unique(dt_runs$freq_bin))) {
  fb_dir <- file.path(base_out, fb)
  dir.create(fb_dir, recursive = TRUE, showWarnings = FALSE)

  for (mk in keep_masks) {
    sub_dt <- dt_runs[freq_bin == fb & mask == mk]

    out_file <- file.path(fb_dir, paste0("bhr_input_", fb, "_", mk, ".txt"))
    fwrite(sub_dt, out_file, sep = "\t", quote = FALSE, na = "NA")

    cat("[DONE] ", fb, " | ", mk, " : ", nrow(sub_dt),
        " variants -> ", out_file, "\n", sep = "")
  }
}