library(tess3r)
library(readr)
library(dplyr)
library(optparse)

option_list = list(
    make_option("--wd", type = "character", default = NULL,
              help = "wd",
              metavar = "character"),
    make_option("--metadata", type = "character", default = NULL,
              help = "metadata with coordinates, csv file",
              metavar = "character"),
    make_option("--mat_bin", type = "character", default = NULL,
              help = "metadata with coordinates",
              metavar = "character")
              )

opt_parser = OptionParser(option_list = option_list);
opt = parse_args(opt_parser);

wd <- opt$wd
meta <- opt$metadata
mat_bin <- opt$mat_bin


# METADATA
# Make GPS coordinates uniform across site/location
# Keep columns numeric
metadata <- readr::read_csv(meta, col_names = T) %>%
    select(wgs_id, Site_latitude, Site_longitude, Site, Country, Region) %>%
     mutate(latitude = as.numeric(Site_latitude), longitude = as.numeric(Site_longitude))

# MATRIX
# Binary matrix
# Extract header and save sample order
matbin <- data.table::fread(mat_bin)
matbin_s <- matbin[ , -c(1:3)]
samples <- colnames(matbin_s)


# Check overlap
# Sample save
metadata <- metadata %>% filter(wgs_id %in% samples)
matbin_s <- matbin_s %>% select(metadata$wgs_id)
saveRDS(colnames(matbin_s), file.path(wd, "/tess3r_samples.rds"))

# Matrix save
# Transposed recoded mixed calls to alternative - 0.5 to 1
matbin_s[matbin_s == 0.5] <- 1
all(colnames(matbin_s) == metadata$wgs_id)
matbinT <- t(matbin_s)
saveRDS(matbinT, file.path(wd, "/tess3r_matbinT.mis.rds"))

# Coordinates
# Save as matrix
coords <- metadata %>% select(longitude, latitude) %>%
      mutate(longitude = as.numeric(longitude), latitude = as.numeric(latitude))
colnames(coords) <- c("V1", "V2")
coordsM <- data.matrix(coords, rownames.force=FALSE)
saveRDS(coordsM, file.path(wd, "/tess3r_coordsM.rds"))