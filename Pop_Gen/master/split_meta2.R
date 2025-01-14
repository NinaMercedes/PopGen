# Example how to create sub files in R
library(readr)
library(dplyr)
library(optparse)

option_list <- list( 
  make_option(c("-s", "--suffix"), type="character", default="pf_meta",
    help="Suffix to add to files"),
  make_option(c("-m", "--metadata"), type="character", default="metadata.tsv",
    help="Metadata_file, tsv format"),
  make_option(c("-p", "--population"), type="character", default="country",
    help="Population class column"),
  make_option(c("-w", "--wgs_id"), type="character", default="wgs_id",
    help="wgs_id column")  
    )
# Parse Options
parser <- OptionParser(option_list=option_list)
opt = parse_args(parser)

# Functions


make_files <- function(suff, meta, pop, wgs){
  met <- readr::read_tsv(meta) # Load metadata
  suffix <- suff # Add additional suffix after country name
  met %>%
    group_by(.data[[pop]]) %>%
    dplyr::select(.data[[wgs]], .data[[pop]]) %>%
    group_walk(~ write.table(.x, paste0(.y[[pop]], "_", suff, ".txt"), quote = F, row.names = F, col.names = F))
}



make_files(opt$suffix, opt$metadata, opt$population, opt$wgs_id)
