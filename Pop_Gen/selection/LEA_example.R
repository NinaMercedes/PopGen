# =====================
#        SETUP
# =====================
library(LEA)
library(data.table)
library(dplyr)

# ----------- INPUT FILES (edit these paths) -----------
geno_path <- "PATH/TO/YOUR/genotype_matrix.mat.bin"
sample_list_path <- "PATH/TO/YOUR/sample_list.txt"
metadata_path <- "PATH/TO/YOUR/metadata.tsv"

# ----------- OUTPUT FILES -----------
geno_output_path <- "OUTPUT/filtered_genotypes.geno"
snmf_project_prefix <- "OUTPUT/project_prefix"
cross_entropy_pdf <- "OUTPUT/cross_entropy_plot.pdf"
ancestry_pdf <- "OUTPUT/ancestry_barplots.pdf"
lfmm_output_csv <- "OUTPUT/lfmm_association_results.csv"

# =====================
#     LOAD + FILTER
# =====================

# Load genotype matrix
geno_file <- fread(geno_path)

# Create SNP IDs
snpids <- paste0(geno_file$chr, "_", geno_file$pos, "_", geno_file$ref, "_", rownames(geno_file))

# Load samples to keep
samples_of_interest <- read.table(sample_list_path)$V1

# Filter genotype matrix for selected samples
geno_filtered <- geno_file %>% select(any_of(samples_of_interest))
sample_names <- colnames(geno_filtered)

# Transpose to samples x SNPs
geno_t <- t(geno_filtered)
colnames(geno_t) <- snpids
rownames(geno_t) <- sample_names

# Remove monomorphic or missing SNPs (all 0 or NA)
is_all_0_or_na <- apply(geno_t, 2, function(col) all(is.na(col) | col == 0))
geno_t <- geno_t[, !is_all_0_or_na]

# Extract one SNP as environmental variable
env_snp_id <- "Pf3D7_13_v3_1725259_C_747398"  # CHANGE if needed
env_snp <- geno_t[, env_snp_id]
geno_t <- geno_t[, colnames(geno_t) != env_snp_id]

# Convert to matrix and recode for LEA
geno_mat <- t(geno_t)  # SNPs x individuals
geno_mat[geno_mat == 0.5] <- NA  # Treat 0.5 as missing
allele_freqs <- rowMeans(geno_mat, na.rm = TRUE)
maf <- pmin(allele_freqs, 1 - allele_freqs)

# Apply MAF filter (e.g., > 0.01)
maf_threshold <- 0.01
geno_mat <- geno_mat[maf >= maf_threshold, ]
snp_ids_filtered <- rownames(geno_mat)

# Recode for LEA (missing = 9)
geno_mat[is.na(geno_mat)] <- 9
storage.mode(geno_mat) <- "integer"

# Write to .geno format
write.table(geno_mat,
            file = geno_output_path,
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "")

# =====================
#         SNMF
# =====================
project <- snmf(geno_output_path,
                K = 1:12,
                repetitions = 5,
                entropy = TRUE,
                project = "new",
                CPU = 20,
                ploidy = 1,
                iterations = 100)

# Cross-entropy plot
pdf(cross_entropy_pdf, width = 7, height = 5)
plot(project, lwd = 5, col = "blue", pch = 1)
dev.off()

# Ancestry barplots
pdf(ancestry_pdf, width = 10, height = 4)
for (k in 1:10) {
  best_run <- which.min(cross.entropy(project, K = k))
  qmatrix <- Q(project, K = k, run = best_run)

  barplot(t(qmatrix),
          col = rainbow(k), border = NA, space = 0,
          main = paste("Ancestry Coefficients (K =", k, ")"),
          xlab = "Individuals", ylab = "Ancestry")
}
dev.off()

# Select K (adjust as needed)
best_K <- 7
best_run <- which.min(cross.entropy(project, K = best_K))

# =====================
#        METADATA
# =====================

# Build environmental matrix
env <- data.frame(sample_names, env = env_snp)
meta <- read.csv(metadata_path, sep = "\t")
env <- left_join(env, meta, by = c("sample_names" = "wgs_id"))  # Adjust join key if needed

# Define environmental variables
env$time <- ifelse(env$Year >= 2008, 1, 0)  # Example temporal proxy
env$Year_std <- scale(env$Year)
env$Long_std <- scale(env$Admin.level.1.longitude)
env$Lat_std <- scale(as.numeric(env$Admin.level.1.latitude))

# Final environmental matrix
env_mat <- as.matrix(env %>% select(env, Year_std, Long_std, Lat_std))
env_mat[is.na(env_mat[, 1]), 1] <- median(env_mat[, 1], na.rm = TRUE)
env_var_numeric <- as.numeric(env_mat[, 1])

# =====================
#        LFMM2
# =====================
# Impute missing values
impute(project, geno_output_path, method = 'mode', K = best_K, run = best_run)

# Run LFMM2
mod <- lfmm2("output_0.95.lfmm_imputed.lfmm", env_var_numeric, K = best_K)
pv <- lfmm2.test(object = mod, input.file = "output_0.95.lfmm_imputed.lfmm",
                 env = env_var_numeric, linear = FALSE, genomic.control = FALSE)

# Save results
results <- data.frame(snp_id = snp_ids_filtered,
                      p = pv$pvalues,
                      adjust_p = p.adjust(pv$pvalues, method = "BH"))
write.csv(results, lfmm_output_csv, row.names = FALSE)

# Optional: print significant hits
sig_hits <- results %>% filter(adjust_p < 0.05)
cat("Number of significant SNPs (FDR < 0.05):", nrow(sig_hits), "\n")
