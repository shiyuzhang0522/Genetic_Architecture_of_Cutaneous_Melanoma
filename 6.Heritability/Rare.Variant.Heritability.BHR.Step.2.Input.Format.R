################## Rscript to format the input for BHR  ##################
################## Author: Shelley ##################
################## Date: 2025-12-16 ##################

######################### 1. load packages ###################################
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

######################### 2. read data into R ################################

# paths
mask_dir    <- "Coding.GT"
counts_dir  <- "UKB_RAP/counts"

# chr index
chrs <- 1:22

# small helper: assert file exists
stop_if_missing <- function(x) {
  miss <- x[!file.exists(x)]
  if (length(miss) > 0) stop("Missing files:\n", paste(miss, collapse = "\n"))
  invisible(TRUE)
}

# 2.1 mask info: chr*.coding_only.long.csv (as a list)
mask_files <- file.path(mask_dir, sprintf("chr%d.coding_only.long.csv", chrs))
stop_if_missing(mask_files)

mask_list <- lapply(mask_files, function(f) {
  dt <- fread(f)
  # optional: keep track of chr from filename
  dt[, chr_file := as.integer(str_match(basename(f), "^chr(\\d+)\\.")[[1, 2]])]
  dt
})
names(mask_list) <- paste0("chr", chrs)

cat("[INFO] Loaded mask files:", length(mask_list), "chromosomes\n")
cat("[INFO] Example mask dt (chr1):", nrow(mask_list[[1]]), "rows x", ncol(mask_list[[1]]), "cols\n")

# 2.2 counts info in the cases: chr*.case.ALT_counts.acount (as a list)
case_files <- file.path(counts_dir, sprintf("chr%d.case.ALT_counts.acount", chrs))
stop_if_missing(case_files)

case_counts_list <- lapply(case_files, function(f) {
  dt <- fread(f)
  dt[, chr_file := as.integer(str_match(basename(f), "^chr(\\d+)\\.")[[1, 2]])]
  dt
})
names(case_counts_list) <- paste0("chr", chrs)

cat("[INFO] Loaded case .acount files:", length(case_counts_list), "chromosomes\n")
cat("[INFO] Example case acount (chr1):", nrow(case_counts_list[[1]]), "rows x", ncol(case_counts_list[[1]]), "cols\n")

# 2.3 counts info in the controls: chr*.control.ALT_counts.acount (as a list)
ctrl_files <- file.path(counts_dir, sprintf("chr%d.control.ALT_counts.acount", chrs))
stop_if_missing(ctrl_files)

control_counts_list <- lapply(ctrl_files, function(f) {
  dt <- fread(f)
  dt[, chr_file := as.integer(str_match(basename(f), "^chr(\\d+)\\.")[[1, 2]])]
  dt
})
names(control_counts_list) <- paste0("chr", chrs)

cat("[INFO] Loaded control .acount files:", length(control_counts_list), "chromosomes\n")
cat("[INFO] Example control acount (chr1):", nrow(control_counts_list[[1]]), "rows x", ncol(control_counts_list[[1]]), "cols\n")

# quick peek at the column names in the plink2 .acount (helps next steps)
cat("[INFO] Columns in case acount:\n")
print(names(case_counts_list[[1]]))

######################### 3. merge case + control counts per variant ##########
# Goal: for each chr, merge .acount tables by variant ID (and keep REF/ALT)

merge_case_control_counts_one_chr <- function(case_dt, ctrl_dt) {
  case_dt <- as.data.table(case_dt)
  ctrl_dt <- as.data.table(ctrl_dt)

  # keep only needed columns and rename
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

  # merge on variant identity
  merged <- merge(
    case_dt2,
    ctrl_dt2,
    by = c("chr", "ID", "REF", "ALT"),
    all = TRUE
  )
  merged
}

# build a list: per chr merged counts
counts_merged_list <- vector("list", length(chrs))
names(counts_merged_list) <- paste0("chr", chrs)

for (chr in chrs) {
  key <- paste0("chr", chr)
  cat("[INFO] merging case/control counts for", key, "...\n")

  counts_merged_list[[key]] <- merge_case_control_counts_one_chr(
    case_dt = case_counts_list[[key]],
    ctrl_dt = control_counts_list[[key]]
  )

  cat("[INFO] ", key, ": merged rows =", nrow(counts_merged_list[[key]]), "\n")
}

# quick look
str(counts_merged_list$chr22)
head(counts_merged_list$chr22)

######################### 4. merge variant counts with mask info (per chr) ####
merge_mask_with_counts_one_chr <- function(mask_dt, counts_dt) {
  mask_dt   <- as.data.table(mask_dt)
  counts_dt <- as.data.table(counts_dt)

  # keep/rename the join key in mask table
  mask_dt2 <- mask_dt %>%
    as.data.table() %>%
    .[, .(
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

  # counts table already has: chr, ID, REF, ALT, ALT_CTS_case, OBS_CT_case, ALT_CTS_control, OBS_CT_control
  counts_dt2 <- counts_dt %>%
    as.data.table() %>%
    .[, .(
      chr,
      ID,
      REF,
      ALT,
      ALT_CTS_case,
      OBS_CT_case,
      ALT_CTS_control,
      OBS_CT_control
    )]

  # merge
  merged <- merge(
    mask_dt2,
    counts_dt2,
    by.x = "SNP_ID_final",
    by.y = "ID",
    all.x = TRUE
  )

  # sanity checks
  merged[, missing_counts := is.na(ALT_CTS_case) & is.na(ALT_CTS_control)]
  n_miss <- merged[, sum(missing_counts)]
  if (n_miss > 0) {
    cat("[WARN] Missing counts for", n_miss, "variants (ALT_CTS_case/control both NA)\n")
  }
  merged
}

mask_counts_merged_list <- vector("list", length(chrs))
names(mask_counts_merged_list) <- paste0("chr", chrs)

for (chr in chrs) {
  key <- paste0("chr", chr)
  cat("[INFO] merging mask + counts for", key, "...\n")

  mask_counts_merged_list[[key]] <- merge_mask_with_counts_one_chr(
    mask_dt   = mask_list[[key]],
    counts_dt = counts_merged_list[[key]]
  )

  cat("[INFO] ", key, ": rows =", nrow(mask_counts_merged_list[[key]]),
      "; cols =", ncol(mask_counts_merged_list[[key]]), "\n")
}

# quick peek
str(mask_counts_merged_list$chr22)
head(mask_counts_merged_list$chr22)

######################### 5. retrieve gene_start (TSS) using biomaRt ##########
# Notes:
# - We query Ensembl for gene coordinates and strand.
# - TSS definition:
#     strand ==  1  -> TSS = start_position
#     strand == -1  -> TSS = end_position
# - If your GENE IDs look like "ENSG... .12", we strip the version suffix.

######################### 5.1 load biomaRt ###################################
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(biomaRt)
})

######################### 5.2 collect unique genes ############################
all_genes <- unique(unlist(lapply(mask_counts_merged_list, \(dt) unique(dt$GENE))))
all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]

# strip Ensembl version suffix if present (ENSG... .xx)
all_genes_clean <- unique(gsub("\\..*$", "", all_genes))

cat("[INFO] Unique genes (raw):   ", length(all_genes), "\n")
cat("[INFO] Unique genes (clean): ", length(all_genes_clean), "\n")
# [INFO] Unique genes (raw):    18459 
# [INFO] Unique genes (clean):  18459

######################### 5.3 connect to Ensembl (GRCh38) #####################
# try a few mirrors in case one is down
connect_ensembl <- function() {
  mirrors <- c("www", "useast", "uswest", "asia")
  for (m in mirrors) {
    cat("[INFO] Trying Ensembl mirror:", m, "\n")
    mart <- tryCatch(
      useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = m),
      error = function(e) NULL
    )
    if (!is.null(mart)) return(mart)
  }
  stop("Failed to connect to Ensembl via all mirrors.")
}

mart <- connect_ensembl()

######################### 5.4 query gene coordinates #########################
# attributes: chr, start, end, strand for each gene
attrs <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand")

gene_pos <- getBM(
  attributes = attrs,
  filters    = "ensembl_gene_id",
  values     = all_genes_clean,
  mart       = mart
)

gene_pos_dt <- as.data.table(gene_pos)

# keep standard chromosomes (optional)
# gene_pos_dt <- gene_pos_dt[chromosome_name %in% c(1:22, "X", "Y")]

cat("[INFO] Rows returned from biomaRt:", nrow(gene_pos_dt), "\n")
cat("[INFO] Unique genes returned:     ", uniqueN(gene_pos_dt$ensembl_gene_id), "\n")
# [INFO] Rows returned from biomaRt: 18459 
# [INFO] Unique genes returned:      18459

######################### 5.5 compute TSS ####################################
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

# sanity: genes not found
not_found <- setdiff(all_genes_clean, gene_tss_dt$GENE)
cat("[INFO] Genes not found in biomaRt:", length(not_found), "\n")
if (length(not_found) > 0) {
  print(head(not_found, 20))
}

######################### 5.6 merge TSS back to per-chr variant tables ########
mask_counts_with_tss_list <- vector("list", length(chrs))
names(mask_counts_with_tss_list) <- paste0("chr", chrs)

for (chr in chrs) {
  key <- paste0("chr", chr)
  dt <- as.data.table(mask_counts_merged_list[[key]])

  # create a clean gene id column for joining
  dt[, GENE_clean := gsub("\\..*$", "", GENE)]

  dt2 <- merge(
    dt,
    gene_tss_dt,
    by.x = "GENE_clean",
    by.y = "GENE",
    all.x = TRUE
  )

  # restore original GENE column name expectations
  setnames(dt2, "GENE_clean", "GENE_clean")  # keep for sanity

  n_miss_tss <- dt2[, sum(is.na(gene_tss))]
  cat("[INFO] ", key, ": variants=", nrow(dt2), "; missing gene_tss=", n_miss_tss, "\n")

  mask_counts_with_tss_list[[key]] <- dt2
}

# quick check
head(mask_counts_with_tss_list$chr22[, .(SNP_ID_final, GENE, GENE_clean, gene_chr, gene_tss)])

######################### 6. merge 22 chromosomes into one data.table #########

# ensure each element is a data.table
mask_counts_with_tss_list <- lapply(mask_counts_with_tss_list, as.data.table)

# row-bind all chromosomes (fast)
mask_counts_with_tss_dt <- rbindlist(mask_counts_with_tss_list, use.names = TRUE, fill = TRUE, idcol = "chr_key")

cat("[INFO] Combined rows:", nrow(mask_counts_with_tss_dt), "\n")
cat("[INFO] Combined cols:", ncol(mask_counts_with_tss_dt), "\n")
print(head(mask_counts_with_tss_dt))

######################### 7. save combined data.table as RData ################
out_dir <- "/public/home/hpc8301200407/WGS/Melanoma_WGS/BHR/pipeline/input"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "mask_counts_with_tss_dt.RData")

mask_counts_with_tss_dt <- as.data.table(mask_counts_with_tss_dt)

save(mask_counts_with_tss_dt, file = out_file)

cat("[INFO] Saved RData to:", out_file, "\n")


########## 8. ALT allele frequency (ALT_AF) in cases/controls/overall ##########
dt <- as.data.table(mask_counts_with_tss_dt)

## 8.0 sanity check: NA counts before filling
n_na_case    <- dt[, sum(is.na(ALT_CTS_case))]
n_na_control <- dt[, sum(is.na(ALT_CTS_control))]
n_na_any     <- dt[, sum(is.na(ALT_CTS_case) | is.na(ALT_CTS_control))]

cat("[INFO] NA ALT_CTS_case           :", n_na_case, "\n")
cat("[INFO] NA ALT_CTS_control        :", n_na_control, "\n")
cat("[INFO] NA in case OR control     :", n_na_any, "\n")

## 8.2 ALT allele frequency in cases / controls / pooled
# NOTE: For plink2 --freq counts (.acount), OBS_CT is number of allele observations,
dt[, ALT_AF_case    := fifelse(OBS_CT_case    > 0, ALT_CTS_case    / OBS_CT_case,    NA_real_)]
dt[, ALT_AF_control := fifelse(OBS_CT_control > 0, ALT_CTS_control / OBS_CT_control, NA_real_)]

dt[, ALT_CTS_total := ALT_CTS_case + ALT_CTS_control]
dt[, OBS_CT_total  := OBS_CT_case  + OBS_CT_control]
dt[, ALT_AF_total  := fifelse(OBS_CT_total > 0, ALT_CTS_total / OBS_CT_total, NA_real_)]

## 8.3 quick sanity peek
dt[, .(
  SNP_ID_final, REF, ALT,
  ALT_CTS_case, OBS_CT_case, ALT_AF_case,
  ALT_CTS_control, OBS_CT_control, ALT_AF_control,
  ALT_CTS_total, OBS_CT_total, ALT_AF_total
)][1:10]

## 8.4 drop variants with ALT dosage = 0 in BOTH cases and controls
# assign to a new object: dt_filtered

dt_filtered <- dt[!(ALT_CTS_case == 0L & ALT_CTS_control == 0L)]

cat("[INFO] Original variants:", nrow(dt), "\n")
cat("[INFO] Kept variants:", nrow(dt_filtered), "\n")
cat("[INFO] Dropped variants:", nrow(dt) - nrow(dt_filtered), "\n")

# quick sanity peek
dt_filtered[, .(
  SNP_ID_final, REF, ALT,
  ALT_CTS_case, ALT_CTS_control,
  ALT_AF_case, ALT_AF_control, ALT_AF_total
)][1:10]

range(dt_filtered$ALT_AF_total)
# [1] 1.454677e-06 1.000000e+00

## 8.5 compute sample prevalence (for per-SD beta)
n_cases    <- 3816
n_controls <- 339903

prevalence <- n_cases / (n_cases + n_controls)

cat("[INFO] n_cases     =", n_cases, "\n")
cat("[INFO] n_controls  =", n_controls, "\n")
cat("[INFO] prevalence  =", prevalence, "\n")

# [INFO] prevalence  = 0.01110209 

########## 9. per-SD and per-allele beta (adapted to dt_filtered) ##########
table <- as.data.table(dt_filtered)

# constants
n_cases    <- 3816
n_controls <- 339903
prevalence <- n_cases / (n_cases + n_controls)

cat("[INFO] prevalence =", prevalence, "\n")

############################
# 9.1 allele frequencies
############################
# already correctly defined earlier:
#   ALT_AF_case   : cases
#   ALT_AF_total  : cases + controls (pooled)

table[, AF_case := ALT_AF_case]
table[, AF      := ALT_AF_total]

############################
# 9.2 per-SD beta (BHR)
############################
# beta = (2 * (AF_case - AF) * prevalence) /
#        sqrt(2 * AF * (1 - AF) * prevalence * (1 - prevalence))

table[, beta := (2 * (AF_case - AF) * prevalence) /
              sqrt(2 * AF * (1 - AF) * prevalence * (1 - prevalence))]

############################
# 9.3 variant variance (2pq)
############################
table[, twopq := 2 * AF * (1 - AF)]

############################
# 9.4 per-allele beta
############################
# convert per-SD beta → per-allele beta
table[, beta_perallele := beta / sqrt(twopq)]

############################
# 9.5 handle numerical edge cases
############################
bad <- !is.finite(table$beta) | !is.finite(table$beta_perallele)
cat("[INFO] Non-finite beta rows:", sum(bad), "\n")

table[bad, c("beta", "twopq", "beta_perallele") :=
        .(NA_real_, NA_real_, NA_real_)]

############################
# 9.6 sanity peek
############################
table[, .(
  SNP_ID_final, REF, ALT,
  ALT_CTS_case, ALT_CTS_control,
  AF_case, AF,
  beta, twopq, beta_perallele
)][1:10]

############################
# 10. subset + rename columns for BHR input
############################
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

# quick sanity check
str(bhr_dt)
bhr_dt[1:10]

############################
# 12. write BHR input table to txt
############################
out_dir <- "/public/home/hpc8301200407/WGS/Melanoma_WGS/BHR/pipeline/input"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "bhr_input_melanoma.txt")

fwrite(
  bhr_dt,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("[INFO] BHR input written to:\n", out_file, "\n")

########### 12. split into runs by MAC/MAF bins × mask types ##################
dt <- as.data.table(bhr_dt)

# keep only 3 mask groups of interest
keep_masks <- c("pLoF", "damaging_missense_or_protein_altering", "synonymous")
dt <- dt[mask %in% keep_masks]

# total ALT allele count across cases+controls
dt[, MAC_total := ALT_CTS_case + ALT_CTS_control]

# use ALT frequency (pooled) directly (your AF column is ALT_AF_total by construction)
dt[, ALT_AF := AF]

########### make a filter based on call rates #################
##### call rates (based on allele observations) #####
N_total <- n_cases + n_controls
obs_max <- 2 * N_total

# call rate
dt[, callrate_total := OBS_CT_total / obs_max]

# count variants failing callrate < 95%
n_low_call <- dt[callrate_total < 0.95, .N]
cat("[INFO] Variants with callrate_total < 0.95:", n_low_call, "\n")

# (optional) also show the distribution of callrate_total for context
cat("[INFO] callrate_total summary:\n")
print(dt[, .(
  min = min(callrate_total, na.rm = TRUE),
  mean = mean(callrate_total, na.rm = TRUE),
  median = median(callrate_total, na.rm = TRUE),
  max = max(callrate_total, na.rm = TRUE)
)])

# remove low-call variants (keep in new object)
dt_call95 <- dt[callrate_total >= 0.95]

cat("[INFO] Variants kept after callrate>=0.95:", nrow(dt_call95), "\n")
cat("[INFO] Variants removed:", nrow(dt) - nrow(dt_call95), "\n")

############################
# define bins AFTER callrate filter (no ALT_AF constraint for MAC<=10)
############################

# use dt_call95 from previous step
dt_call95[, MAC_total := ALT_CTS_case + ALT_CTS_control]
dt_call95[, ALT_AF := AF]   # pooled ALT frequency

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

# sanity check: range/mean ALT_AF + counts per bin
cat("[INFO] freq_bin summary after callrate>=0.95:\n")
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
            by = freq_bin
  ][order(freq_bin)]
)

############################
# save dt_call95 to disk
############################
library(data.table)

out_file <- "/public/home/hpc8301200407/WGS/Melanoma_WGS/BHR/pipeline/input/dt_call95.txt"

fwrite(
  dt_call95,
  file = out_file,
  sep = "\t",
  quote = FALSE,
  na = "NA"
)

cat("[DONE] dt_call95 written to:\n", out_file, "\n")

# quick sanity check
cat("[INFO] rows:", nrow(dt_call95), " cols:", ncol(dt_call95), "\n")
print(head(dt_call95))

############################
# 15. write per-freq_bin folders, split by mask
############################
library(data.table)

keep_masks <- c("pLoF", "damaging_missense_or_protein_altering", "synonymous")

dt_runs <- dt_call95[!is.na(freq_bin) & mask %in% keep_masks]

base_out <- "/public/home/hpc8301200407/WGS/Melanoma_WGS/BHR/pipeline/input/runs_by_freqbin"
dir.create(base_out, recursive = TRUE, showWarnings = FALSE)

# sanity: counts per bin × mask
cat("[INFO] variants per freq_bin × mask:\n")
print(dt_runs[, .N, by = .(freq_bin, mask)][order(freq_bin, mask)])

for (fb in sort(unique(dt_runs$freq_bin))) {

  fb_dir <- file.path(base_out, fb)
  dir.create(fb_dir, recursive = TRUE, showWarnings = FALSE)

  for (mk in keep_masks) {
    sub_dt <- dt_runs[freq_bin == fb & mask == mk]

    out_file <- file.path(fb_dir, paste0("bhr_input_", fb, "_", mk, ".txt"))
    fwrite(sub_dt, out_file, sep = "\t", quote = FALSE, na = "NA")

    cat("[DONE] ", fb, " | ", mk, " : ", nrow(sub_dt), " variants -> ", out_file, "\n", sep = "")
  }
}