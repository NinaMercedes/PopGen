################################################################################
# FILTER HIGH-FST SNPs BY GENE WINDOWS AND WATCHLIST
# Please note this code has been cleaned using Claude
################################################################################

library(data.table)
library(stringr)
library(dplyr)
library(GenomicRanges)

# ==============================================================================
# USER-DEFINED INPUTS
# ==============================================================================

promoter_size  <- 1000   # bp upstream/downstream of gene
fst_threshold  <- 0.2    # Fst threshold for high-Fst SNPs
fst_dir        <- "path/to/downstream_analysis_directory"  # Directory containing *_high_fst_snps_detailed.tsv files
gff_file       <- "path/to/Pfalciparum.genome.modified.new.gff3"
watchlist_file <- "path/to/watchlist.txt"
output_dir     <- "path/to/output_directory"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# LOAD GFF
# ==============================================================================

cat("Loading gene annotations...\n")
gff <- read.delim(gff_file, header = FALSE, comment.char = "#")
setnames(gff, c("seqname","source","type","start","end","score","strand","phase","attributes"))

genes <- gff[gff$type == "gene", ]

# Convert to data.table and clean chromosome names once
setDT(genes)
genes[, seqname := gsub("Pf3D7_0", "", seqname)]
genes[, seqname := gsub("Pf3D7_",  "", seqname)]
genes[, seqname := gsub("_v3",     "", seqname)]

# Extract gene_ID and gene_name
genes[, gene_ID   := sub(".*ID=([^;]+).*",   "\\1", attributes)]
genes[, gene_name := sub(".*Name=([^;]+).*", "\\1", attributes)]

# Convert to GRanges
genes_gr <- GRanges(
  seqnames  = genes$seqname,
  ranges    = IRanges(start = genes$start, end = genes$end),
  gene_ID   = genes$gene_ID,
  gene_name = genes$gene_name
)

# Prepare genes as data.table for joining
genes_dt <- data.table(
  gene_ID   = genes_gr$gene_ID,
  gene_name = genes_gr$gene_name,
  seqname   = as.character(seqnames(genes_gr)),
  gene_start = start(genes_gr),
  gene_end   = end(genes_gr)
)

cat("  Loaded", nrow(genes_dt), "genes\n\n")

# ==============================================================================
# LOAD WATCHLIST
# ==============================================================================

cat("Loading watchlist...\n")
watchlist       <- fread(watchlist_file, header = FALSE)$V1
watchlist_clean <- paste0("gene:", watchlist)
cat("  Watchlist genes:", length(watchlist_clean), "\n\n")

# ==============================================================================
# FIND INPUT FILES
# ==============================================================================

fst_files <- list.files(fst_dir, pattern = "_high_fst_snps_detailed.tsv",
                        full.names = TRUE, recursive = TRUE)

if (length(fst_files) == 0) {
  stop("No *_high_fst_snps_detailed.tsv files found in: ", fst_dir)
}

cat("Found", length(fst_files), "region files to process\n\n")

# ==============================================================================
# PROCESS EACH REGION
# ==============================================================================

for (i in seq_along(fst_files)) {
  fname  <- basename(fst_files[i])
  region <- gsub("_high_fst_snps_detailed.tsv", "", fname)
  cat("Processing region:", region, "\n")

  fst_dt <- fread(fst_files[i])

  if (nrow(fst_dt) == 0) {
    cat("  Empty file, skipping.\n")
    next
  }

  # Join with gene coordinates using explicit, unambiguous column names
  fst_dt <- merge(fst_dt, genes_dt, by = "gene_ID", all.x = TRUE)

  # Remove SNPs with no gene coordinates
  fst_dt <- fst_dt[!is.na(gene_start) & !is.na(gene_end)]

  if (nrow(fst_dt) == 0) {
    cat("  No SNPs matched to gene coordinates, skipping.\n")
    next
  }

  # Keep SNPs within gene +/- promoter window
  fst_dt <- fst_dt[
    POS >= (gene_start - promoter_size) &
    POS <= (gene_end   + promoter_size)
  ]

  # Apply Fst threshold
  fst_dt <- fst_dt[Fst >= fst_threshold]

  # Deduplicate
  fst_dt <- unique(fst_dt)

  cat("  SNPs passing filters:", nrow(fst_dt), "\n")

  # Save all passing SNPs
  out_all <- file.path(output_dir, paste0(region, "_fst_filtered_all.tsv"))
  fwrite(fst_dt, out_all, sep = "\t")
  cat("  Saved all filtered SNPs to:", out_all, "\n")

  # Filter for watchlist genes and save separately
  fst_watchlist <- fst_dt[gene_ID %in% watchlist_clean]
  cat("  Watchlist SNPs:", nrow(fst_watchlist), "\n")

  out_watchlist <- file.path(output_dir, paste0(region, "_fst_filtered_watchlist.tsv"))
  fwrite(fst_watchlist, out_watchlist, sep = "\t")
  cat("  Saved watchlist SNPs to:", out_watchlist, "\n")
}

cat("\nAll regions processed.\n")
