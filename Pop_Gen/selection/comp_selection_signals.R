################################################################################
# GENOMIC SELECTION ANALYSIS INTEGRATION - CLEAN XP-EHH HANDLING
# Combines Tajima's D, iHS, XP-EHH, and IBD analyses
################################################################################

library(dplyr)
library(stringr)
library(data.table)

# ==============================================================================
# PART 1: LOAD AND PREPARE TAJIMA'S D DATA
# ==============================================================================

cat("Loading Tajima's D data...\n")

setwd("C:/Users/ninam/OneDrive - London School of Hygiene and Tropical Medicine/Pf_09_24 - Copy/pop_gen/Tajima")

list_files <- list.files(pattern = ".Tajima.D")
list_f <- list()
for (file in list_files) {
  list_f[[file]] <- read.csv(file, sep = "\t")
}

tajima_d <- rbindlist(list_f, idcol = TRUE)
tajima_d$Region <- gsub("_Pf_Nov_24.txt.genotyped.vcf.gz_tajima.txt.Tajima.D", "", tajima_d$.id)
tajima_d$Region <- gsub("_", " ", tajima_d$Region)
tajima_d <- tajima_d %>% filter(!CHROM %in% c("Pf_M76611", "Pf3D7_API_v3"))

percentiles <- tajima_d %>%
  group_by(Region) %>%
  summarise(p_low = quantile(TajimaD, 0.10, na.rm = TRUE),
            p_high = quantile(TajimaD, 0.90, na.rm = TRUE),
            .groups = "drop")

tajima_d2 <- tajima_d %>% left_join(percentiles, by = "Region")
extremes_df <- tajima_d2 %>% filter(!is.na(TajimaD) & (TajimaD <= p_low))

cat("  Extreme negative Tajima's D windows:", nrow(extremes_df), "\n\n")

# ==============================================================================
# PART 2: LOAD SELECTION DATA
# ==============================================================================

cat("Loading selection data...\n")
cr_ihs <- read.delim("C:/Users/ninam/OneDrive - London School of Hygiene and Tropical Medicine/Pf_09_24 - Copy/pop_gen/Selection/new/cr_ihs_all_categories_annot.tsv")
cr_xpehh <- read.delim("C:/Users/ninam/OneDrive - London School of Hygiene and Tropical Medicine/Pf_09_24 - Copy/pop_gen/Selection/new/cr_xpehh_all_categories_annot.tsv")
cr_xpehh <- cr_xpehh[!(grepl("Horn", cr_xpehh$category_name)),]
cr_ihs <- cr_ihs[!(grepl("Horn", cr_ihs$category_name)),]
cat("  iHS regions:", nrow(cr_ihs), "\n")
cat("  XP-EHH regions:", nrow(cr_xpehh), "\n\n")


# ==============================================================================
# PART 3: LOAD IBD DATA
# ==============================================================================

cat("Loading IBD data...\n")
ibd <- read.delim("C:/Users/ninam/OneDrive - London School of Hygiene and Tropical Medicine/Pf_09_24 - Copy/pop_gen/IBD/18_12_2024_hmmIBD_ibd_win10kb.tsv")

percentiles_ibd <- ibd %>%
  group_by(category) %>%
  summarise(p_low = quantile(fraction, 0.15, na.rm = TRUE),
            p_high = quantile(fraction, 0.85, na.rm = TRUE),
            .groups = "drop")

ibd <- ibd %>%
  left_join(percentiles_ibd, by = "category") %>%
  filter(!is.na(fraction), fraction >= p_high)

cat("  High IBD windows:", nrow(ibd), "\n\n")

# ==============================================================================
# PART 4: STANDARDIZE FORMATS
# ==============================================================================

cat("Standardizing formats...\n")

# Tajima's D
TD <- extremes_df %>% select(CHROM, BIN_START, N_SNPS, TajimaD, Region)
colnames(TD) <- c("Chr", "TD_start", "TD_SNP_n", "TD", "TD_Region")
TD$Chr <- gsub("Pf3D7_0", "", TD$Chr)
TD$Chr <- gsub("Pf3D7_", "", TD$Chr)
TD$Chr <- gsub("_v3", "", TD$Chr)
TD$Chr <- as.numeric(as.character(TD$Chr))
TD$TD_end <- TD$TD_start + 9999

# iHS
cr_ihs <- cr_ihs %>% select(chr, start, end, category_name, genes, products)
colnames(cr_ihs) <- c("Chr", "iHS_start", "iHS_end", "iHS_Region", "iHS_genes", "iHS_products")
cr_ihs$iHS_Region <- gsub("_", " ", cr_ihs$iHS_Region)

# IBD
ibd <- ibd %>% select(chr, win_start, win_end, fraction, category)
colnames(ibd) <- c("Chr", "IBD_start", "IBD_end", "IBD_fraction", "IBD_Region")
ibd$IBD_Region <- gsub("_", " ", ibd$IBD_Region)

# XP-EHH - Split comparison into two populations
cr_xpehh <- cr_xpehh %>% select(chr, start, end, category_name, PERC_EXTR_MRK)
colnames(cr_xpehh) <- c("Chr", "XPEHH_start", "XPEHH_end", "XPEHH_Comparison", "XPEHH_extra_marks_p")
cr_xpehh$XPEHH_Comparison <- gsub("_", " ", cr_xpehh$XPEHH_Comparison)

# Convert to data.table
setDT(TD)
setDT(cr_ihs)
setDT(ibd)
setDT(cr_xpehh)

cat("  Done.\n\n")

# ==============================================================================
# PART 5: SPLIT XP-EHH INTO TWO POPULATIONS
# ==============================================================================

cat("Splitting XP-EHH comparisons into population pairs...\n")

# Split "Pop1|Pop2" into separate columns
cr_xpehh[, c("XPEHH_pop1", "XPEHH_pop2") := tstrsplit(XPEHH_Comparison, "|", fixed = TRUE)]
cr_xpehh[, XPEHH_pop1 := trimws(XPEHH_pop1)]
cr_xpehh[, XPEHH_pop2 := trimws(XPEHH_pop2)]

cat("  Split", nrow(cr_xpehh), "XP-EHH regions into pop1 and pop2\n\n")

# ==============================================================================
# PART 6: MERGE STEP-BY-STEP WITH FULL OUTER JOINS
# ==============================================================================

cat("Merging datasets...\n\n")

all_chrs <- unique(c(TD$Chr, cr_ihs$Chr, ibd$Chr, cr_xpehh$Chr))
all_chrs <- sort(all_chrs[!is.na(all_chrs)])

results_list <- list()

for (chr in all_chrs) {
  cat("  Chr", chr, "...\n")
  
  # Subset to chromosome
  TD_chr <- TD[Chr == chr]
  ihs_chr <- cr_ihs[Chr == chr]
  ibd_chr <- ibd[Chr == chr]
  xp_chr <- cr_xpehh[Chr == chr]
  
  # Track rows with IDs
  TD_chr[, .td_id := .I]
  ihs_chr[, .ihs_id := .I]
  ibd_chr[, .ibd_id := .I]
  xp_chr[, .xp_id := .I]
  
  cat("    TD:", nrow(TD_chr), "| iHS:", nrow(ihs_chr), 
      "| IBD:", nrow(ibd_chr), "| XP-EHH:", nrow(xp_chr), "\n")
  
  # === STEP 1: Merge TD + iHS (same region, overlapping intervals) ===
  result <- data.table()
  td_matched <- integer(0)
  ihs_matched <- integer(0)
  
  if (nrow(TD_chr) > 0 && nrow(ihs_chr) > 0) {
    setkeyv(TD_chr, c("Chr", "TD_start", "TD_end"))
    setkeyv(ihs_chr, c("Chr", "iHS_start", "iHS_end"))
    
    ov <- foverlaps(TD_chr, ihs_chr, type = "any", nomatch = NULL)
    ov <- ov[TD_Region == iHS_Region]  # Only keep matching regions
    
    if (nrow(ov) > 0) {
      result <- rbind(result, ov, fill = TRUE)
      td_matched <- ov$.td_id
      ihs_matched <- ov$.ihs_id
    }
  }
  
  # Add TD-only rows
  if (nrow(TD_chr) > 0) {
    td_only <- TD_chr[!.td_id %in% td_matched]
    if (nrow(td_only) > 0) result <- rbind(result, td_only, fill = TRUE)
  }
  
  # Add iHS-only rows
  if (nrow(ihs_chr) > 0) {
    ihs_only <- ihs_chr[!.ihs_id %in% ihs_matched]
    if (nrow(ihs_only) > 0) result <- rbind(result, ihs_only, fill = TRUE)
  }
  
  cat("    After TD+iHS:", nrow(result), "rows\n")
  
  # === STEP 2: Merge with IBD (same region, overlapping intervals) ===
  if (nrow(result) > 0 && nrow(ibd_chr) > 0) {
    # Create unified interval for result
    result[, unified_start := pmin(TD_start, iHS_start, na.rm = TRUE)]
    result[, unified_end := pmax(TD_end, iHS_end, na.rm = TRUE)]
    result[, .res_id := .I]
    
    # Separate valid and NA intervals
    res_valid <- result[!is.na(unified_start) & !is.na(unified_end)]
    res_na <- result[is.na(unified_start) | is.na(unified_end)]
    
    if (nrow(res_valid) > 0) {
      setkeyv(res_valid, c("Chr", "unified_start", "unified_end"))
      setkeyv(ibd_chr, c("Chr", "IBD_start", "IBD_end"))
      
      ov_ibd <- foverlaps(res_valid, ibd_chr, type = "any", nomatch = NULL)
      # Match if any existing region matches IBD
      ov_ibd <- ov_ibd[
        (!is.na(TD_Region) & TD_Region == IBD_Region) |
          (!is.na(iHS_Region) & iHS_Region == IBD_Region)
      ]
      
      res_matched <- ov_ibd[!is.na(.res_id), unique(.res_id)]
      ibd_matched <- ov_ibd[!is.na(.ibd_id), unique(.ibd_id)]
      
      res_unmatched <- res_valid[!.res_id %in% res_matched]
      
      result <- rbind(ov_ibd, res_unmatched, res_na, fill = TRUE)
      
      # Add IBD-only
      ibd_only <- ibd_chr[!.ibd_id %in% ibd_matched]
      if (nrow(ibd_only) > 0) result <- rbind(result, ibd_only, fill = TRUE)
    } else {
      result <- rbind(result, ibd_chr, fill = TRUE)
    }
    
    result[, c("unified_start", "unified_end", ".res_id") := NULL]
  } else if (nrow(ibd_chr) > 0) {
    result <- rbind(result, ibd_chr, fill = TRUE)
  }
  
  cat("    After +IBD:", nrow(result), "rows\n")
  
  # === STEP 3: Merge with XP-EHH (region matches pop1 OR pop2, overlapping) ===
  if (nrow(result) > 0 && nrow(xp_chr) > 0) {
    result[, unified_start := pmin(TD_start, iHS_start, IBD_start, na.rm = TRUE)]
    result[, unified_end := pmax(TD_end, iHS_end, IBD_end, na.rm = TRUE)]
    result[, .res_id := .I]
    
    res_valid <- result[!is.na(unified_start) & !is.na(unified_end)]
    res_na <- result[is.na(unified_start) | is.na(unified_end)]
    
    if (nrow(res_valid) > 0) {
      setkeyv(res_valid, c("Chr", "unified_start", "unified_end"))
      setkeyv(xp_chr, c("Chr", "XPEHH_start", "XPEHH_end"))
      
      ov_xp <- foverlaps(res_valid, xp_chr, type = "any", nomatch = NULL)
      
      # Match if any existing region matches EITHER pop1 OR pop2
      ov_xp <- ov_xp[
        (!is.na(TD_Region) & (TD_Region == XPEHH_pop1 | TD_Region == XPEHH_pop2)) |
          (!is.na(iHS_Region) & (iHS_Region == XPEHH_pop1 | iHS_Region == XPEHH_pop2)) |
          (!is.na(IBD_Region) & (IBD_Region == XPEHH_pop1 | IBD_Region == XPEHH_pop2))
      ]
      
      res_matched <- ov_xp[!is.na(.res_id), unique(.res_id)]
      xp_matched <- ov_xp[!is.na(.xp_id), unique(.xp_id)]
      
      res_unmatched <- res_valid[!.res_id %in% res_matched]
      
      result <- rbind(ov_xp, res_unmatched, res_na, fill = TRUE)
      
      # Add XP-only
      xp_only <- xp_chr[!.xp_id %in% xp_matched]
      if (nrow(xp_only) > 0) result <- rbind(result, xp_only, fill = TRUE)
    } else {
      result <- rbind(result, xp_chr, fill = TRUE)
    }
    
    result[, c("unified_start", "unified_end", ".res_id") := NULL]
  } else if (nrow(xp_chr) > 0) {
    result <- rbind(result, xp_chr, fill = TRUE)
  }
  
  cat("    After +XP-EHH:", nrow(result), "rows\n")
  
  # Clean up temp columns
  temp_cols <- grep("^\\.|^i\\.", names(result), value = TRUE)
  if (length(temp_cols) > 0) result[, (temp_cols) := NULL]
  
  results_list[[as.character(chr)]] <- result
  
  rm(TD_chr, ihs_chr, ibd_chr, xp_chr, result)
  gc(verbose = FALSE)
}

# ==============================================================================
# PART 7: COMBINE AND FINALIZE
# ==============================================================================

cat("\nCombining all chromosomes...\n")
result <- rbindlist(results_list, fill = TRUE, use.names = TRUE)

# Create combined intervals
result[, `:=`(
  combined_start = pmin(TD_start, iHS_start, IBD_start, XPEHH_start, na.rm = TRUE),
  combined_end = pmax(TD_end, iHS_end, IBD_end, XPEHH_end, na.rm = TRUE)
)]

setorder(result, Chr, combined_start, na.last = TRUE)

# ==============================================================================
# PART 8: SUMMARY
# ==============================================================================

cat("\n", rep("=", 80), "\n", sep = "")
cat("INTEGRATION COMPLETE\n")
cat(rep("=", 80), "\n\n", sep = "")

cat("Final dataset: 'result'\n")
cat("  Total rows:", nrow(result), "\n")
cat("  Chromosomes:", paste(sort(unique(result$Chr)), collapse = ", "), "\n\n")

# Count dataset contributions
result[, n_datasets := (!is.na(TD_Region)) + (!is.na(iHS_Region)) + 
         (!is.na(IBD_Region)) + (!is.na(XPEHH_Comparison))]

cat("Dataset overlap distribution:\n")
print(table(result$n_datasets))
cat("\n")

# Show which combinations exist
cat("Dataset combinations (top 10):\n")
result[, combo := paste(
  ifelse(!is.na(TD_Region), "TD", ""),
  ifelse(!is.na(iHS_Region), "iHS", ""),
  ifelse(!is.na(IBD_Region), "IBD", ""),
  ifelse(!is.na(XPEHH_Comparison), "XP", ""),
  sep = "+"
)]
result[, combo := gsub("\\+{2,}", "+", combo)]
result[, combo := gsub("^\\+|\\+$", "", combo)]
print(head(sort(table(result$combo), decreasing = TRUE), 10))

result[, c("n_datasets", "combo") := NULL]

cat("\nFirst 20 rows:\n")
print(head(result, 20))

cat("\n", rep("=", 80), "\n", sep = "")
cat("Result stored in 'result' data.table\n")
cat(rep("=", 80), "\n", sep = "")

# Optional: Save
# fwrite(result, "integrated_selection_analysis.tsv", sep = "\t")