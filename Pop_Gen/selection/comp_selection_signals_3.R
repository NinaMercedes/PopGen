################################################################################
# FILTER HIGH-FST SNPs BY GENE WINDOWS AND WATCHLIST
################################################################################

library(data.table)
library(stringr)
library(GenomicRanges)

# -----------------------------
# PARAMETERS
# -----------------------------
promoter_size <- 1000  # bp upstream/downstream of gene
fst_dir <- "C:/Users/ninam/OneDrive - London School of Hygiene and Tropical Medicine/Pf_09_24 - Copy/pop_gen/Selection/downstream_analysis"
gff_file <- "C:/Users/ninam/OneDrive - London School of Hygiene and Tropical Medicine/Pf_09_24/pop_gen/Fst/Pfalciparum.genome.modified.new.gff3"
watchlist_file <- "C:/Users/ninam/OneDrive - London School of Hygiene and Tropical Medicine/Pf_09_24 - Copy/pop_gen/Selection/new/watchlist.txt"
output_dir <- file.path("C:/Users/ninam/OneDrive - London School of Hygiene and Tropical Medicine/Pf_09_24 - Copy/pop_gen/Selection/filtered_fst_watchlist")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# LOAD GFF
# -----------------------------
cat("Loading gene annotations...\n")
gff <- read.delim(gff_file, header=FALSE, comment.char="#")
setnames(gff, c("seqname","source","type","start","end","score","strand","phase","attributes"))
genes <- gff[gff$type == "gene",]
genes$seqname <- gsub("Pf3D7_0","",genes$seqname)
genes$seqname <- gsub("Pf3D7_","",genes$seqname)
genes$seqname <- gsub("_v3","",genes$seqname)

# Convert to data.table first
setDT(genes)

# Clean chromosome names
genes[, seqname := gsub("Pf3D7_0","",seqname)]
genes[, seqname := gsub("Pf3D7_","",seqname)]
genes[, seqname := gsub("_v3","",seqname)]

# Extract gene_ID and gene_name
genes[, gene_ID := sub(".*ID=([^;]+).*", "\\1", attributes)]
genes[, gene_name := sub(".*Name=([^;]+).*", "\\1", attributes)]

# Convert to GRanges
genes_gr <- GRanges(
  seqnames = genes$seqname,
  ranges = IRanges(start=genes$start, end=genes$end),
  gene_ID = genes$gene_ID,
  gene_name = genes$gene_name
)

# -----------------------------
# LOAD WATCHLIST
# -----------------------------
watchlist <- fread(watchlist_file, header=FALSE)$V1

# -----------------------------
# PROCESS EACH REGION
# -----------------------------
# Clean watchlist
watchlist_clean <- paste0("gene:",watchlist)

for (i in 1:length(fst_files)) {
  fname <- basename(fst_files[[i]])
  region <- gsub("_high_fst_snps_detailed.tsv","",fname)
  cat("\nProcessing region:", region, "\n")
  
  # Read detailed Fst
  fst_dt <- fread(fst_files[[i]])
  
  
  # Merge with GFF to get coordinates
  genes_dt <- data.table(
    gene_ID = genes_gr$gene_ID,
    gene_name = genes_gr$gene_name,
    seqname = as.character(seqnames(genes_gr)),
    start = start(genes_gr),
    end = end(genes_gr)
  )
  
  fst_dt <- left_join(fst_dt, genes_dt, by=c("gene_ID"))
  
  
  # Remove SNPs with no gene coordinates
  fst_dt <- fst_dt[!is.na(start.y) & !is.na(end.y)]
  
  # Keep SNPs within gene ± promoter
  fst_dt <- fst_dt[
    POS >= (start.y - promoter_size) &
      POS <= (end.y + promoter_size)
  ]
  
  # Keep high Fst SNPs
  fst_dt <- fst_dt[Fst >= 0.2]
  fst_dt <- fst_dt %>% unique()
  out_file <- file.path(output_dir, paste0(region, "_fst_filtered_all.tsv"))
  fwrite(fst_dt, out_file, sep="\t")
  # Filter for watchlist genes
  fst_dt <- fst_dt[gene_ID %in% watchlist_clean]
  # Save
  out_file <- file.path(output_dir, paste0(region, "_fst_filtered_watchlist.tsv"))
  fwrite(fst_dt, out_file, sep="\t")
  cat("  Saved filtered Fst for watchlist genes to:", out_file, "\n")
}

cat("\nAll regions processed.\n")
