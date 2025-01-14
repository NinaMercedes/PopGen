library(dplyr)
library(data.table)
library(optparse)

# Arguments ####
option_list = list(
  make_option(c("-p", "--prefix"), type = "character", default = NULL,
              help = "prefix", metavar = "character"),
  make_option(c("-w", "--wgs_id"), type = "character", default = NULL,
              help = "wgs_id", metavar = "character"),      
  make_option(c("-m", "--metadata"), type = "character", default = NULL,
              help = "metadata tsv file", metavar = "character")      
  
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
list_files <- list.files(pattern=paste0("*", opt$prefix, "_moi_fws.tsv"))
final_list <- list()
for (i in 1:length(list_files)){
    final_list[[i]] <-  read.table(paste0(list_files[i]), header = TRUE)
    colnames(final_list[[i]]) <- gsub("sample", opt$wgs_id, colnames(final_list[[i]]))
}
combine <- rbindlist(final_list)
meta <- read.csv(opt$meta, header=TRUE, sep = "\t")
meta_fws <- left_join(meta, combine)
meta_fws_comp <- meta_fws %>% filter(!fws %in% c(NA))
write.table(meta_fws_comp, paste0(opt$prefix, "_metadata_fws.tsv"), quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
