# Load required libraries
library(data.table)
library(dplyr)
library(parallel)
library(stringr)

# Set working directory and number of cores
setwd("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/epistasis/)
n_cores <- 16

# Load genotype matrix
geno <- fread("interaction_candidates.mat.bin", sep = "\t", header = TRUE, data.table = FALSE)




# Load PC covariates
pcs <- fread("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/GEMMA/cmd_points.csv")
colnames(pcs) <- c("sample", "PC1", "PC2", "PC3", "PC4","PC5")

# Keep only samples present in both genotype and PCs
common_samples <- intersect(colnames(geno)[4:ncol(geno)], pcs$sample)
geno <- geno[, c("chr", "pos", "ref", common_samples)]
pcs <- pcs %>% filter(sample %in% common_samples)

# Extract kelch13 mutation row (e.g., by annotation)
kelch_row <- 8958
kelch_geno <- as.numeric(geno[kelch_row, 4:ncol(geno)])

# Filter genotype matrix to exclude kelch13 itself
geno_test <- geno 

# Parallel association test
geno_rows <- split(geno_test, seq(nrow(geno_test)))



results <- mclapply(geno_rows, function(row) {
  snp2 <- as.numeric(row[1, 4:ncol(row)])
  interaction_term <- kelch_geno * snp2
  
  if (sum(interaction_term, na.rm = TRUE) > 1) {
    print("passed interaction")
    chr <- row["chr"]
    pos <- row["pos"]
    ref <- row["ref"]
    df <- data.frame(
      y = kelch_geno,
      snp2 = snp2,
      PC1 = pcs$PC1,
      PC2 = pcs$PC2,
      PC3 = pcs$PC3,
      PC4 = pcs$PC4, 
      PC5 = pcs$PC5
    )
    df <- df[complete.cases(df), ]
    
    if (nrow(df) > 10) {
      print("exceeds sample number")
      fit <- glm(y ~ snp2 + PC1 + PC2 + PC3 + PC4 + PC5, data = df, family = "binomial")
      coef_table <- coef(summary(fit))
      return(data.frame(
        chr = chr,
        pos = pos,
        ref = ref,
        beta_snp2 = coef_table["snp2", "Estimate"],
        p_snp2 = coef_table["snp2", "Pr(>|z|)"],
        n = nrow(df)
      ))
    }
  }
  return(NULL)
}, mc.cores = n_cores)

# Combine and save results
results_df <- do.call(rbind, results)
results_df$p_adj <- p.adjust(results_df$p_snp2, method="BH")
write.csv(results_df, "kelch13_interaction_results.csv", row.names = FALSE)
