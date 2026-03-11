################################################################################
# EPISTASIS ANALYSIS - KELCH13 INTERACTION SCREEN
# Please note this code has been cleaned using Claude
################################################################################

library(data.table)
library(dplyr)
library(parallel)
library(stringr)

# ==============================================================================
# USER-DEFINED INPUTS
# ==============================================================================

geno_file      <- "path/to/interaction_candidates.mat.bin"
pcs_file       <- "path/to/cmd_points.csv"
output_file    <- "path/to/kelch13_interaction_results.csv"
n_cores        <- 16
fst_threshold  <- 0.2    # kept for consistency if needed downstream
kelch_row      <- NULL   # Set to row number (e.g. 8958) or NULL to auto-detect
kelch_id       <- NULL   # Set to a SNP identifier to auto-detect kelch row (e.g. "PF3D7_1343700")
                         # Only used if kelch_row is NULL

# ==============================================================================
# LOAD DATA
# ==============================================================================

cat("Loading genotype matrix...\n")
geno <- fread(geno_file, sep = "\t", header = TRUE, data.table = FALSE)
cat("  Dimensions:", nrow(geno), "SNPs x", ncol(geno) - 3, "samples\n\n")

cat("Loading PC covariates...\n")
pcs <- fread(pcs_file)
colnames(pcs) <- c("sample", "PC1", "PC2", "PC3", "PC4", "PC5")

# ==============================================================================
# ALIGN SAMPLES
# ==============================================================================

cat("Aligning samples...\n")
common_samples <- intersect(colnames(geno)[4:ncol(geno)], pcs$sample)
cat("  Common samples:", length(common_samples), "\n\n")

geno <- geno[, c("chr", "pos", "ref", common_samples)]

# Align PC rows to match genotype column order exactly
pcs <- pcs[match(common_samples, pcs$sample), ]

stopifnot(all(pcs$sample == common_samples))  # Safety check

# ==============================================================================
# EXTRACT KELCH13 ROW
# ==============================================================================

if (!is.null(kelch_row)) {
  cat("Using user-specified kelch13 row:", kelch_row, "\n")
} else if (!is.null(kelch_id)) {
  cat("Auto-detecting kelch13 row from identifier:", kelch_id, "\n")
  kelch_row <- which(apply(geno[, c("chr","pos","ref")], 1, function(r) any(grepl(kelch_id, r))))
  if (length(kelch_row) == 0) stop("Could not find kelch_id in genotype matrix. Check identifier.")
  if (length(kelch_row) > 1) {
    warning("Multiple rows match kelch_id — using first match.")
    kelch_row <- kelch_row[1]
  }
  cat("  Found kelch13 at row:", kelch_row, "\n")
} else {
  stop("Must specify either kelch_row (integer) or kelch_id (string) to identify the kelch13 SNP.")
}

kelch_geno <- as.numeric(geno[kelch_row, 4:ncol(geno)])
cat("  Kelch13 genotype — 0s:", sum(kelch_geno == 0, na.rm = TRUE),
    "| 1s:", sum(kelch_geno == 1, na.rm = TRUE),
    "| NAs:", sum(is.na(kelch_geno)), "\n\n")

# ==============================================================================
# FILTER GENOTYPE MATRIX
# ==============================================================================

# Exclude kelch13 row from testing
geno_test <- geno[-kelch_row, ]
cat("Testing", nrow(geno_test), "SNPs\n\n")

# ==============================================================================
# PARALLEL ASSOCIATION TEST
# ==============================================================================

cat("Running parallel association tests on", n_cores, "cores...\n")

results <- mclapply(seq_len(nrow(geno_test)), function(i) {

  row  <- geno_test[i, ]
  snp2 <- as.numeric(row[1, 4:ncol(row)])

  # Only test SNPs that co-occur with kelch13 mutation
  interaction_term <- kelch_geno * snp2
  if (sum(interaction_term, na.rm = TRUE) <= 1) return(NULL)

  chr <- row[["chr"]]
  pos <- row[["pos"]]
  ref <- row[["ref"]]

  df <- data.frame(
    y   = kelch_geno,
    snp2 = snp2,
    PC1  = pcs$PC1,
    PC2  = pcs$PC2,
    PC3  = pcs$PC3,
    PC4  = pcs$PC4,
    PC5  = pcs$PC5
  )
  df <- df[complete.cases(df), ]

  if (nrow(df) <= 10) return(NULL)

  fit <- tryCatch(
    glm(y ~ snp2 + PC1 + PC2 + PC3 + PC4 + PC5, data = df, family = "binomial"),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NULL)

  coef_table <- coef(summary(fit))

  if (!"snp2" %in% rownames(coef_table)) return(NULL)

  return(data.frame(
    chr      = chr,
    pos      = pos,
    ref      = ref,
    beta_snp2 = coef_table["snp2", "Estimate"],
    p_snp2    = coef_table["snp2", "Pr(>|z|)"],
    n         = nrow(df)
  ))

}, mc.cores = n_cores)

# ==============================================================================
# COMBINE AND SAVE
# ==============================================================================

cat("Combining results...\n")
results_df <- do.call(rbind, Filter(Negate(is.null), results))

if (is.null(results_df) || nrow(results_df) == 0) {
  stop("No SNPs passed the co-occurrence filter. Check kelch13 genotype or interaction threshold.")
}

cat("  SNPs passing filters:", nrow(results_df), "\n")

results_df$p_adj <- p.adjust(results_df$p_snp2, method = "BH")
results_df       <- results_df[order(results_df$p_adj), ]

write.csv(results_df, output_file, row.names = FALSE)
cat("Results saved to:", output_file, "\n")
