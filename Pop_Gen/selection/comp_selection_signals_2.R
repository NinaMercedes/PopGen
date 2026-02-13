################################################################################
# DOWNSTREAM ANALYSIS LOOP - ALL GEOGRAPHIC REGIONS
# Creates UpSet plots, extracts genes, and analyzes high Fst SNPs
################################################################################

library(dplyr)
library(stringr)
library(data.table)
library(GenomicRanges)
library(UpSetR)
library(ggplot2)

# ==============================================================================
# SETUP: Define parameters and load reference data
# ==============================================================================

promoter_size <- 1000  # bp upstream/downstream of TD windows

# Define color palette
upset_colors <- c("#77AADD", "#EE8866", 
                  "#EEDD88", "#FFAABB",  "#44BB99",
                  "#AA66CC")

# Load GFF (genes)
cat("Loading gene annotations...\n")
gff <- read.delim("C:/Users/ninam/OneDrive - London School of Hygiene and Tropical Medicine/Pf_09_24/pop_gen/Fst/Pfalciparum.genome.modified.new.gff3", 
                  header=FALSE, comment.char="#")
setnames(gff, c("seqname","source","type","start","end","score","strand","phase","attributes"))

# Keep genes only and clean chromosome names
genes <- gff[gff$type == "gene",]
genes$seqname <- gsub("Pf3D7_0","",genes$seqname)
genes$seqname <- gsub("Pf3D7_","",genes$seqname)
genes$seqname <- gsub("_v3","",genes$seqname)

setDT(genes)
genes[, gene_ID := sub(".*ID=([^;]+).*", "\\1", attributes)]
genes[, gene_name := sub(".*Name=([^;]+).*", "\\1", attributes)]

# Convert to GRanges
genes_gr <- GRanges(
  seqnames = genes$seqname,
  ranges = IRanges(start = genes$start, end = genes$end),
  gene_ID = genes$gene_ID,
  gene_name = genes$gene_name
)

# Load watchlist
watchlist <- read.table("C:/Users/ninam/OneDrive - London School of Hygiene and Tropical Medicine/Pf_09_24 - Copy/pop_gen/Selection/new/watchlist.txt")

cat("Setup complete.\n\n")

# ==============================================================================
# IDENTIFY ALL UNIQUE REGIONS
# ==============================================================================

cat("Identifying geographic regions...\n")

# Create presence/absence columns
df <- result %>%
  mutate(
    TD_present   = ifelse(!is.na(TD), 1, 0),
    iHS_present   = ifelse(!is.na(iHS_start), 1, 0),
    IBD_present   = ifelse(!is.na(IBD_start), 1, 0),
    XPEHH_present = ifelse(!is.na(XPEHH_start), 1, 0)
  )

# Get all unique regions
all_regions <- unique(c(
  df$TD_Region[!is.na(df$TD_Region)],
  df$iHS_Region[!is.na(df$iHS_Region)],
  df$IBD_Region[!is.na(df$IBD_Region)]
))

# Also extract regions from XP-EHH comparisons
xp_regions <- df$XPEHH_Comparison[!is.na(df$XPEHH_Comparison)]
xp_split <- unique(unlist(strsplit(xp_regions, "\\|")))
xp_split <- trimws(xp_split)

all_regions <- unique(c(all_regions, xp_split))
all_regions <- all_regions[!is.na(all_regions) & all_regions != ""]

cat("Found", length(all_regions), "geographic regions:\n")
cat(paste(all_regions, collapse="\n"), "\n\n")

# ==============================================================================
# CREATE OUTPUT DIRECTORY
# ==============================================================================

output_dir <- "C:/Users/ninam/OneDrive - London School of Hygiene and Tropical Medicine/Pf_09_24 - Copy/pop_gen/Selection/downstream_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# MAIN LOOP: Process each region
# ==============================================================================

# Store results for summary
summary_results <- list()

for (region in all_regions) {
  # use as test case: region <- all_regions[3]
  cat("\n", rep("=", 80), "\n", sep="")
  cat("PROCESSING REGION:", region, "\n")
  cat(rep("=", 80), "\n\n", sep="")
  
  # Create region-specific subdirectory
  region_dir <- file.path(output_dir, gsub(" ", "_", region))
  dir.create(region_dir, showWarnings = FALSE)
  
  # --------------------------------------------------------------------------
  # STEP 1: Filter data for this region
  # --------------------------------------------------------------------------
  
  cat("Step 1: Filtering data for", region, "...\n")
  
  df_region <- df %>% filter(
    TD_Region == region |
      IBD_Region == region |
      iHS_Region == region |
      (
        str_detect(XPEHH_Comparison, fixed(region)) &
          (TD_Region == region | IBD_Region == region | iHS_Region == region)
      )
  )
  
  cat("  Found", nrow(df_region), "genomic windows\n")
  
  if (nrow(df_region) == 0) {
    cat("  No data for this region. Skipping.\n")
    next
  }
  
  # --------------------------------------------------------------------------
  # STEP 2: Create UpSet plot (before filtering)
  # --------------------------------------------------------------------------
  
  cat("\nStep 2: Creating UpSet plot...\n")
  
  mat <- df_region %>% select(TD_present, iHS_present, IBD_present, XPEHH_present)
  
  cat("  Matrix dimensions:", nrow(mat), "rows x", ncol(mat), "columns\n")
  cat("  Column sums:\n")
  print(colSums(mat))
  
  # Save plot
  plot_file <- file.path(region_dir, paste0(gsub(" ", "_", region), "_upset_plot.png"))
  
  png(plot_file, width=10, height=6, units="in", res=300)
  print(upset(
    as.data.frame(mat),
    sets = c("TD_present", "iHS_present", "IBD_present", "XPEHH_present"),
    order.by = "freq",
    empty.intersections = "on",
    main.bar.color = upset_colors[1],
    sets.bar.color = upset_colors[2:5],
    text.scale = 1.5,
    mainbar.y.label = "Intersection Size",
    sets.x.label = "Total Windows"
  ))
  dev.off()
  
  cat("  UpSet plot saved to:", plot_file, "\n")
  
  # --------------------------------------------------------------------------
  # STEP 3: Filter for regions with ≥2 supporting evidence & classify
  # --------------------------------------------------------------------------
  
  cat("\nStep 3: Filtering for ≥2 supporting evidence...\n")
  
  df_region <- df_region %>%
    mutate(n_support = TD_present + iHS_present + IBD_present + XPEHH_present) %>%
    filter(n_support >= 2) %>%  # CHANGED: Now filtering for ≥2
    mutate(evidence_strength = case_when(
      n_support == 2 ~ "weak",
      n_support == 3 ~ "moderate",
      n_support == 4 ~ "high",
      TRUE ~ NA_character_
    ))
  
  cat("  Windows with ≥2 evidence:", nrow(df_region), "\n")
  cat("    Evidence distribution:\n")
  print(table(df_region$evidence_strength))
  
  if (nrow(df_region) == 0) {
    cat("  No windows meet criteria. Skipping.\n")
    summary_results[[region]] <- list(
      n_windows_total = 0,
      n_windows_weak = 0,
      n_windows_moderate = 0,
      n_windows_high = 0,
      n_genes = 0,
      n_watchlist_genes = 0,
      matching_genes_summary = "",
      n_high_fst_snps = 0
    )
    next
  }
  
  # --------------------------------------------------------------------------
  # STEP 4: Extract TD windows and find overlapping genes
  # --------------------------------------------------------------------------
  
  cat("\nStep 4: Finding genes in selected regions...\n")
  
  # -----------------------------
  # FILL TD_START / TD_END if missing
  # -----------------------------
  df_region$TD_start_filled <- df_region$TD_start
  df_region$TD_end_filled   <- df_region$TD_end
  
  # Replace NAs with iHS_start/iHS_end, then IBD_start/IBD_end, then XPEHH_start/XPEHH_end
  df_region$TD_start_filled[is.na(df_region$TD_start_filled)] <- 
    coalesce(df_region$iHS_start[is.na(df_region$TD_start_filled)],
             df_region$IBD_start[is.na(df_region$TD_start_filled)],
             df_region$XPEHH_start[is.na(df_region$TD_start_filled)])
  
  df_region$TD_end_filled[is.na(df_region$TD_end_filled)] <- 
    coalesce(df_region$iHS_end[is.na(df_region$TD_end_filled)],
             df_region$IBD_end[is.na(df_region$TD_end_filled)],
             df_region$XPEHH_end[is.na(df_region$TD_end_filled)])
  
  # Use filled columns for regions
  regions <- df_region %>% select(Chr, TD_start_filled, TD_end_filled,
                                  TD_present, iHS_present, IBD_present, XPEHH_present,
                                  n_support, evidence_strength)
  
  # Rename to original names
  setnames(regions, c("TD_start_filled","TD_end_filled"), c("TD_start","TD_end"))
  
  # Remove rows where still NA
  regions <- regions[!is.na(regions$TD_start), ]
  
  # Ensure Chr is character for consistency
  regions$Chr <- as.character(regions$Chr)
  
  # Convert to GRanges
  regions_gr <- GRanges(
    seqnames = regions$Chr,
    ranges = IRanges(start = regions$TD_start, end = regions$TD_end),
    TD_present = regions$TD_present,
    iHS_present = regions$iHS_present,
    IBD_present = regions$IBD_present,
    XPEHH_present = regions$XPEHH_present,
    n_support = regions$n_support,
    evidence_strength = regions$evidence_strength
  )
  
  # --------------------------------------------------------------------------
  # STEP 4b: Find overlapping genes
  # --------------------------------------------------------------------------
  hits <- findOverlaps(regions_gr, genes_gr)
  
  overlapping <- data.table(
    Chr = as.character(regions$Chr[queryHits(hits)]),
    TD_start = regions$TD_start[queryHits(hits)],
    TD_end = regions$TD_end[queryHits(hits)],
    TD_present = regions$TD_present[queryHits(hits)],
    iHS_present = regions$iHS_present[queryHits(hits)],
    IBD_present = regions$IBD_present[queryHits(hits)],
    XPEHH_present = regions$XPEHH_present[queryHits(hits)],
    n_support = regions$n_support[queryHits(hits)],
    evidence_strength = regions$evidence_strength[queryHits(hits)],
    gene_ID = genes_gr$gene_ID[subjectHits(hits)],
    gene_name = genes_gr$gene_name[subjectHits(hits)]
  )
  
  cat("  Found", length(unique(overlapping$gene_ID)), "unique genes\n")
  
  # --------------------------------------------------------------------------
  # Matching watchlist genes
  # --------------------------------------------------------------------------
  matching_genes <- watchlist$V1[sapply(watchlist$V1, function(g) any(grepl(g, overlapping$gene_ID)))]
  
  if (length(matching_genes) > 0) {
    cat("  Watchlist genes found:", length(matching_genes), "\n")
    cat("    ", paste(matching_genes, collapse=", "), "\n")
    
    matching_genes_with_evidence <- overlapping[
      gene_ID %in% matching_genes | sapply(gene_ID, function(gid) any(grepl(paste(matching_genes, collapse="|"), gid)))
    ][, .(
      max_evidence_strength = evidence_strength[which.max(n_support)],
      max_n_support = max(n_support)
    ), by = gene_ID]
    
    genes_summary_list <- paste0(
      matching_genes_with_evidence$gene_ID, 
      " (", 
      matching_genes_with_evidence$max_evidence_strength, 
      ")"
    )
    matching_genes_summary <- paste(genes_summary_list, collapse = "; ")
  } else {
    cat("  No watchlist genes found\n")
    matching_genes_summary <- ""
  }
  
  # --------------------------------------------------------------------------
  # Save gene summary
  # --------------------------------------------------------------------------
  gene_summary <- overlapping[, .(
    Chr = Chr[1],
    TD_start = min(TD_start),
    TD_end = max(TD_end),
    n_support = max(n_support),
    evidence_strength = evidence_strength[which.max(n_support)],
    TD_present = max(TD_present),
    iHS_present = max(iHS_present),
    IBD_present = max(IBD_present),
    XPEHH_present = max(XPEHH_present)
  ), by = .(gene_ID, gene_name)]
  
  setorder(gene_summary, -n_support, Chr, TD_start)
  
  fwrite(gene_summary, 
         file.path(region_dir, paste0(gsub(" ", "_", region), "_genes_high_support.tsv")),
         sep="\t")
  
  # --------------------------------------------------------------------------
  # STEP 5: Load and process Fst data
  # --------------------------------------------------------------------------
  
  cat("\nStep 5: Analyzing Fst SNPs...\n")
  
  # Construct Fst filename
  fst_file <- file.path(
    "C:/Users/ninam/OneDrive - London School of Hygiene and Tropical Medicine/Pf_09_24 - Copy/pop_gen/Fst",
    paste0(gsub(" ", "_", region), "_vs_all.single.snp_fst_df.csv")
  )
  
  if (!file.exists(fst_file)) {
    cat("  WARNING: Fst file not found:", fst_file, "\n")
    cat("  Skipping Fst analysis for this region.\n")
    
    summary_results[[region]] <- list(
      n_windows_total = nrow(df_region),
      n_windows_weak = sum(df_region$evidence_strength == "weak"),
      n_windows_moderate = sum(df_region$evidence_strength == "moderate"),
      n_windows_high = sum(df_region$evidence_strength == "high"),
      n_genes = length(unique(overlapping$gene_ID)),
      n_watchlist_genes = length(matching_genes),
      matching_genes_summary = matching_genes_summary,
      n_high_fst_snps = NA
    )
    next
  }
  
  # Load Fst data
  fst_data <- fread(fst_file)
  
  # Filter for high Fst
  fst <- fst_data[Fst >= 0.2]
  fst[, CHR := gsub("Pf3D7_0","", CHR)]
  fst[, CHR := gsub("Pf3D7_","", CHR)]
  fst[, CHR := gsub("_v3","", CHR)]
  fst[, CHR := as.character(CHR)]
  
  cat("  Total high Fst SNPs (≥0.4):", nrow(fst), "\n")
  
  if (nrow(fst) == 0) {
    cat("  No high Fst SNPs found.\n")
    summary_results[[region]] <- list(
      n_windows_total = nrow(df_region),
      n_windows_weak = sum(df_region$evidence_strength == "weak"),
      n_windows_moderate = sum(df_region$evidence_strength == "moderate"),
      n_windows_high = sum(df_region$evidence_strength == "high"),
      n_genes = length(unique(overlapping$gene_ID)),
      n_watchlist_genes = length(matching_genes),
      matching_genes_summary = matching_genes_summary,
      n_high_fst_snps = 0
    )
    next
  }
  
  # Ensure overlapping Chr is character
  overlapping[, Chr := as.character(Chr)]
  
  # Prepare SNP table as 1-bp intervals
  fst_dt <- fst[, .(CHR, POS, Fst)]
  fst_dt[, `:=`(start = POS, end = POS)]
  setkeyv(fst_dt, c("CHR", "start", "end"))
  
  # Prepare TD/promoter windows
  td_dt <- overlapping[, .(
    Chr,
    start = pmax(1, TD_start - promoter_size),
    end   = TD_end + promoter_size,
    TD_start,
    TD_end,
    gene_ID,
    gene_name,
    evidence_strength
  )]
  setcolorder(td_dt, c("Chr", "start", "end", "TD_start", "TD_end", "gene_ID", "gene_name", "evidence_strength"))
  setkeyv(td_dt, c("Chr", "start", "end"))
  
  # Perform overlap
  snp_near_TD <- foverlaps(
    fst_dt,
    td_dt,
    by.x = c("CHR", "start", "end"),
    by.y = c("Chr", "start", "end"),
    nomatch = 0L
  )
  
  # Classify SNP relative to TD window
  snp_near_TD[, region_type := fifelse(
    POS < TD_start, "upstream",
    fifelse(POS > TD_end, "downstream", "TD_window")
  )]
  
  cat("  High Fst SNPs in/near selected regions:", nrow(snp_near_TD), "\n")
  
  if (nrow(snp_near_TD) > 0) {
    cat("    By location:\n")
    print(table(snp_near_TD$region_type))
    
    cat("    By evidence strength:\n")
    print(table(snp_near_TD$evidence_strength))
    
    # Summary by gene
    snp_summary <- snp_near_TD[, .(
      Chr = CHR[1],
      evidence_strength = evidence_strength[1],
      n_snps = .N,
      mean_Fst = mean(Fst),
      max_Fst = max(Fst),
      n_upstream = sum(region_type == "upstream"),
      n_TD = sum(region_type == "TD_window"),
      n_downstream = sum(region_type == "downstream")
    ), by = .(gene_ID, gene_name)]
    
    setorder(snp_summary, -n_snps)
    
    # Save SNP results
    fwrite(snp_near_TD, file.path(region_dir, paste0(gsub(" ", "_", region), "_high_fst_snps_detailed.tsv")), sep="\t")
    fwrite(snp_summary, file.path(region_dir, paste0(gsub(" ", "_", region), "_high_fst_snps_by_gene.tsv")), sep="\t")
    cat("  Results saved.\n")
  }
  
  # --------------------------------------------------------------------------
  # Store summary
  # --------------------------------------------------------------------------
  summary_results[[region]] <- list(
    n_windows_total = nrow(df_region),
    n_windows_weak = sum(df_region$evidence_strength == "weak"),
    n_windows_moderate = sum(df_region$evidence_strength == "moderate"),
    n_windows_high = sum(df_region$evidence_strength == "high"),
    n_genes = length(unique(overlapping$gene_ID)),
    n_watchlist_genes = length(matching_genes),
    matching_genes_summary = matching_genes_summary,
    n_high_fst_snps = nrow(snp_near_TD)
  )
  
  cat("\nRegion", region, "complete.\n")
}

# ==============================================================================
# CREATE OVERALL SUMMARY
# ==============================================================================

cat("\n", rep("=", 80), "\n", sep="")
cat("CREATING OVERALL SUMMARY\n")
cat(rep("=", 80), "\n\n", sep="")

summary_df <- data.frame(
  Region = names(summary_results),
  Total_windows_2plus = sapply(summary_results, function(x) x$n_windows_total),
  Windows_weak = sapply(summary_results, function(x) x$n_windows_weak),
  Windows_moderate = sapply(summary_results, function(x) x$n_windows_moderate),
  Windows_high = sapply(summary_results, function(x) x$n_windows_high),
  Unique_genes = sapply(summary_results, function(x) x$n_genes),
  Watchlist_genes = sapply(summary_results, function(x) x$n_watchlist_genes),
  Matching_genes_with_evidence = sapply(summary_results, function(x) x$matching_genes_summary),
  High_Fst_SNPs = sapply(summary_results, function(x) x$n_high_fst_snps)
)

summary_df <- summary_df %>% arrange(desc(Total_windows_2plus))
fwrite(summary_df, file.path(output_dir, "SUMMARY_all_regions.tsv"), sep="\t")

cat("Summary across all regions:\n")
print(summary_df)
cat("\n", rep("=", 80), "\n", sep="")
cat("ANALYSIS COMPLETE\n")
cat("Results saved to:", output_dir, "\n")
cat(rep("=", 80), "\n", sep="")
