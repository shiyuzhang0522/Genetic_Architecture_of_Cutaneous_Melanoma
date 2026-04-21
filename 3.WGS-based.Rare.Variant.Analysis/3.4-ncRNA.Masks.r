############################################################
## Build rare-variant masks for ncRNA aggregates
## Author: Shelley
## Date:  2025-06-24
############################################################

############################################################
## Step 1) Parse command-line arguments
############################################################
args <- commandArgs(trailingOnly = TRUE)
chr <- args[1]
subset_id <- args[2]

cat("[INFO] Processing chromosome: ", chr, "\n", sep = "")
cat("[INFO] Processing ncRNA exon aggregate subset: ", subset_id, "\n", sep = "")

############################################################
## Step 2) Load packages
############################################################
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

############################################################
## Step 3) Read score files and QC-passed variant list
############################################################

## score files
file_CADD_SNV <- paste0(
  "/hpc/home/lijc/zhangsy/WGS_melanoma/WGS_annotation/CADD_annotation/CADD-scripts/",
  "CADD.cutoff.20.SNVs.chr", chr, ".tsv"
)

file_CADD_INDEL <- paste0(
  "/hpc/home/lijc/zhangsy/WGS_melanoma/WGS_annotation/CADD_annotation/CADD-scripts/",
  "CADD.cutoff.20.indels.chr", chr, ".tsv"
)

file_GERP <- paste0(
  "/hpc/home/lijc/zhangsy/WGS_melanoma/WGS_annotation/GERP_annotation/",
  "GERP.cutoff.2.chr", chr, ".tsv"
)

file_JARVIS <- paste0(
  "/hpc/home/lijc/zhangsy/WGS_melanoma/WGS_annotation/JARVIS/",
  "jarvis_ranked_deduplicated.", chr, ".cutoff.0.99.hg38.tsv"
)

CADD_SNV_data   <- fread(file_CADD_SNV, header = FALSE)
CADD_INDEL_data <- fread(file_CADD_INDEL, header = FALSE)
GERP_data       <- fread(file_GERP, header = FALSE)
JARVIS_data     <- fread(file_JARVIS, header = FALSE)

## QC-passed WGS variant list
file_WGSQCpassed <- paste0(
  "/hpc/home/lijc/zhangsy/WGS_melanoma/WGS_annotation/WGS_pvar/download_WGS_pvar/",
  "WGS.chr", chr, ".QCpassed.varlist.vcf"
)

WGS_QCpassed <- fread(
  file_WGSQCpassed,
  skip = "#CHROM",
  select = 1:5,
  header = TRUE,
  sep = "\t"
)

colnames(WGS_QCpassed)[1] <- "CHR"

############################################################
## Step 4) Merge CADD, GERP, and JARVIS scores
############################################################

## build variant IDs for CADD
CADD_SNV_data[, ID := paste0("DRAGEN:chr", V1, ":", V2, ":", V3, ":", V4)]
CADD_INDEL_data[, ID := paste0("DRAGEN:chr", V1, ":", V2, ":", V3, ":", V4)]

setnames(CADD_SNV_data, "V6", "CADD_PHRED")
setnames(CADD_INDEL_data, "V6", "CADD_PHRED")

CADD_data <- rbind(
  CADD_SNV_data[, .(ID, CADD_PHRED)],
  CADD_INDEL_data[, .(ID, CADD_PHRED)]
)

WGS_QCpassed <- WGS_QCpassed %>%
  left_join(as.data.frame(CADD_data), by = "ID")

CADD_percent_above_20 <- WGS_QCpassed %>%
  summarise(
    total_variants = n(),
    count_above_20 = sum(CADD_PHRED >= 20, na.rm = TRUE)
  ) %>%
  mutate(percentage = count_above_20 / total_variants)

message("CADD percentage of variants with PHRED score >= 20:")
print(CADD_percent_above_20)

## GERP
setnames(GERP_data, old = "V3", new = "ID")
setnames(GERP_data, old = "V6", new = "GERP")

GERP_keep <- GERP_data[, .(ID, GERP)]

WGS_QCpassed <- WGS_QCpassed %>%
  left_join(as.data.frame(GERP_keep), by = "ID")

GERP_percent_above_2 <- WGS_QCpassed %>%
  summarise(
    total_variants = n(),
    count_above_2 = sum(GERP >= 2, na.rm = TRUE)
  ) %>%
  mutate(percentage = count_above_2 / total_variants)

message("GERP percentage of variants with score >= 2:")
print(GERP_percent_above_2)

## JARVIS
setnames(JARVIS_data, old = c("V1", "V2", "V3"), new = c("CHR", "POS", "JARVIS"))

jarvis_df <- as.data.frame(JARVIS_data[, .(CHR, POS, JARVIS)]) %>%
  rename(JARVIS_score = JARVIS)

WGS_QCpassed <- WGS_QCpassed %>%
  left_join(jarvis_df, by = c("CHR", "POS"))

JARVIS_percent_above_0.99 <- WGS_QCpassed %>%
  summarise(
    total_variants = n(),
    count_above_0.99 = sum(JARVIS_score >= 0.99, na.rm = TRUE)
  ) %>%
  mutate(percentage = count_above_0.99 / total_variants)

message("Percentage of variants with JARVIS_score >= 0.99:")
print(JARVIS_percent_above_0.99)

############################################################
## Step 5) Assign a unique mask to each variant
############################################################

get_mask <- function(cadd, gerp, jarvis) {
  if (is.na(cadd) & is.na(gerp) & is.na(jarvis)) {
    return("Null")
  }

  mask_parts <- c()

  if (!is.na(cadd) && cadd >= 20) {
    mask_parts <- c(mask_parts, "CADD")
  }

  if (!is.na(gerp) && gerp >= 2) {
    mask_parts <- c(mask_parts, "GERP")
  }

  if (!is.na(jarvis) && jarvis >= 0.99) {
    mask_parts <- c(mask_parts, "JARVIS")
  }

  if (length(mask_parts) == 0) {
    return("None")
  }

  paste(mask_parts, collapse = "_")
}

WGS_QCpassed <- WGS_QCpassed %>%
  mutate(mask = mapply(get_mask, CADD_PHRED, GERP, JARVIS_score))

mask_counts <- WGS_QCpassed %>%
  group_by(mask) %>%
  summarise(count = n())

print(mask_counts)

total_assigned <- sum(mask_counts$count)
total_variants <- nrow(WGS_QCpassed)

if (total_assigned == total_variants) {
  message(
    "The sum of all mask counts (", total_assigned,
    ") equals the total number of variants (", total_variants, ")."
  )
} else {
  message(
    "Mismatch: The sum of mask counts (", total_assigned,
    ") does not equal the total number of variants (", total_variants, ")."
  )
}

## validate score-based counts against mask counts
num_cadd_variants   <- sum(WGS_QCpassed$CADD_PHRED >= 20, na.rm = TRUE)
num_gerp_variants   <- sum(WGS_QCpassed$GERP >= 2, na.rm = TRUE)
num_jarvis_variants <- sum(WGS_QCpassed$JARVIS_score >= 0.99, na.rm = TRUE)

num_mask_CADD   <- sum(str_detect(WGS_QCpassed$mask, "CADD"))
num_mask_GERP   <- sum(str_detect(WGS_QCpassed$mask, "GERP"))
num_mask_JARVIS <- sum(str_detect(WGS_QCpassed$mask, "JARVIS"))

cat("Comparison for CADD:\n")
cat("  Count based on score (CADD_PHRED >= 20): ", num_cadd_variants, "\n", sep = "")
cat("  Count based on mask (contains 'CADD'):   ", num_mask_CADD, "\n\n", sep = "")

cat("Comparison for GERP:\n")
cat("  Count based on score (GERP >= 2):        ", num_gerp_variants, "\n", sep = "")
cat("  Count based on mask (contains 'GERP'):   ", num_mask_GERP, "\n\n", sep = "")

cat("Comparison for JARVIS:\n")
cat("  Count based on score (JARVIS_score >= 0.99): ", num_jarvis_variants, "\n", sep = "")
cat("  Count based on mask (contains 'JARVIS'):      ", num_mask_JARVIS, "\n", sep = "")

## optional save of variant-level mask table
output_file <- paste0(
  "/hpc/home/lijc/zhangsy/WGS_melanoma/masks/",
  "WGS_QCpassed.chr", chr, ".three.score.mask.txt"
)

# fwrite(WGS_QCpassed, file = output_file, sep = "\t")

############################################################
## Step 6) Read ncRNA aggregate definitions
############################################################

ncRNA_file <- paste0(
  "/hpc/home/lijc/zhangsy/WGS_melanoma/RV_aggregates/ncRNA_submask/chr",
  chr, "/ncRNA.exons_chr", chr, "_part_", subset_id, ".txt"
)

ncRNA_aggregates <- fread(ncRNA_file, select = 1:2, header = TRUE)

message("Processing chr ", chr, ", subset ", subset_id, " ncRNA exon aggregates")

############################################################
## Step 7) Expand variants and attach masks
############################################################

## target format:
## Aggregate_name var  ID1   ID2   ID3
## Aggregate_name anno Mask1 Mask2 Mask3

ncRNA_expanded <- ncRNA_aggregates %>%
  separate_rows(variants, sep = "\\|") %>%
  rename(var = variants) %>%
  left_join(WGS_QCpassed %>% select(ID, mask), by = c("var" = "ID")) %>%
  group_by(ncRNA_uniqueID) %>%
  mutate(v_order = row_number()) %>%
  ungroup()

############################################################
## Step 8) Convert to wide format
############################################################

var_wide <- ncRNA_expanded %>%
  select(ncRNA_uniqueID, v_order, var) %>%
  pivot_wider(
    names_from = v_order,
    values_from = var,
    names_prefix = "ID"
  ) %>%
  mutate(row_type = "var") %>%
  select(ncRNA_uniqueID, row_type, everything())

mask_wide <- ncRNA_expanded %>%
  select(ncRNA_uniqueID, v_order, mask) %>%
  pivot_wider(
    names_from = v_order,
    values_from = mask,
    names_prefix = "Mask"
  ) %>%
  mutate(row_type = "anno") %>%
  select(ncRNA_uniqueID, row_type, everything())

############################################################
## Step 9) Save outputs
############################################################

output_dir <- paste0("/hpc/home/lijc/zhangsy/WGS_melanoma/ncRNA_mask/chr", chr)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

var_wide_file <- paste0(
  output_dir,
  "/var_wide.ncRNA.exons.chr", chr, ".subset", subset_id, ".txt"
)

mask_wide_file <- paste0(
  output_dir,
  "/mask_wide.ncRNA.exons.chr", chr, ".subset", subset_id, ".txt"
)

fwrite(var_wide, file = var_wide_file, sep = " ", quote = FALSE, col.names = FALSE)
fwrite(mask_wide, file = mask_wide_file, sep = " ", quote = FALSE, col.names = FALSE)

message(
  "ncRNA exon-only mask has been successfully created for chr",
  chr, ", subset ", subset_id
)