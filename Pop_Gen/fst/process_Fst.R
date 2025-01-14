library(dplyr)
library(data.table)
library(stringr)
library(geosphere)
library(ggplot2)


setwd("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/Pf_Nov24_region_FST_25_11_24/")

## First stage is to read in the files 
file_list1 <- list.files(pattern="*.paired.snp_fst_df.csv")
df_list1_95 <- list()
df_list1_75 <- list()
n_list1 <- list()
n_list2 <- list()
for (i in 1:length(file_list1)){
    file <- read.csv(file_list1[i], header=TRUE)
    file1 <- file %>% filter(Fst>0.95)
    file2 <- file %>% filter(Fst>0.75)
    n1 <- nrow(file1)
    n2 <- nrow(file2)
    if (n1 > 0){
        df_list1_95[[i]] <- file1
        df_list1_95[[i]]$comparison <- file_list1[i]
    }
    if (n2 > 0){
        df_list1_75[[i]] <- file2
        df_list1_75[[i]]$comparison <- file_list1[i]
    }
    n_list1[[i]] <- data.frame(file_list1[i], n1)
    n_list2[[i]] <- data.frame(file_list1[i], n2)
}

df_list1_75_2 <- rbindlist(df_list1_75)
df_list1_95_2 <- rbindlist(df_list1_95)

frequencies <- read.table("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/e4a34f3d-a503-49b4-83af-2717f986bd39.frq", header=TRUE)
annotations <- read.csv("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/annotations_final.tsv", sep="\t", header=TRUE)

colnames(annotations) <- c("CHR", "POS", "REF", "ALT", "GENE")
colnames(frequencies) <- c("CHR", "POS", "REF", "ALT", "MAF", "NCHROBS")

df_list1_75_3 <- left_join(df_list1_75_2, annotations)
df_list1_75_3 <- left_join(df_list1_75_3, frequencies)
df_list1_75_4 <- df_list1_75_3 %>% filter(MAF<1)
df_list1_75_4 <- df_list1_75_4 %>% filter(!GENE %in% "")


df_list1_95_3 <- left_join(df_list1_95_2, annotations)
df_list1_95_3 <- left_join(df_list1_95_3, frequencies)
df_list1_95_4 <- df_list1_95_3 %>% filter(MAF<1)
df_list1_95_4 <- df_list1_95_4 %>% filter(!GENE %in% "")


write.csv(df_list1_75_4, "snps_0.75_region.fst")
write.csv(df_list1_95_4, "snps_0.95_region.fst")


n1_95 <- df_list1_95_4 %>% group_by(comparison) %>% summarise(n=n())
n2_75 <- df_list1_75_4 %>% group_by(comparison) %>% summarise(n=n())



write.csv(n1_95, "num_snps_0.95_region.fst")

write.csv(n2_75, "num_snps_0.75_region.fst")


## Let's get overall and draw a graph with correlation
setwd("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/Pf_Nov24_country_FST_25_11_24/")
file_list2 <- list.files(pattern="*.paired.hudson_fst_df.csv")
df_list2 <- list()
for (i in 1:length(file_list2)){
    file <- read.csv(file_list2[i], header=TRUE)
    df_list2[[i]] <- file
}

fst_average <- rbindlist(df_list2)
head(fst_average)
splits <- str_split_fixed(fst_average$Region, "_vs_", 2)
fst_average$Region1 <- splits[,1]
fst_average$Region2 <- splits[,2]
head(fst_average)

metadata <- read.csv("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/Pf_Nov24_metadata_complete.tsv", header=TRUE, sep="\t")
metadata$Region1 <- metadata$Country
metadata$Lat1 <- metadata$Country.latitude
metadata$Long1 <- metadata$Country.longitude
metadata$Lat2 <- metadata$Country.latitude
metadata$Long2 <- metadata$Country.longitude
metadata$Region2 <- metadata$Country
metadata1 <- metadata %>% select(Region1,  Lat1, Long1)
metadata1 <- unique(metadata1)
metadata2 <- metadata %>% select(Region2,  Lat2, Long2)
metadata2 <- unique(metadata2)
fst_average <- left_join(fst_average, metadata1)
fst_average <- left_join(fst_average, metadata2)
head(fst_average)

for ( i in 1:nrow(fst_average)){
fst_average$dist[i] <- distm(c( fst_average$Long1[i], fst_average$Lat1[i]), c(fst_average$Long2[i],fst_average$Lat2[i]), fun = distHaversine)
}
plot <- ggplot(fst_average, aes(x=dist, y=Hudson_Fst)) + geom_point() + theme_classic()
ggsave("plot_test.png",plot)
corr <- cor.test(as.numeric(as.character(fst_average$dist)), as.numeric(as.character(fst_average$Hudson_Fst)), use="complete.obs")
corr
str(corr)
corr <- data.frame(corr$p.value, corr$estimate, corr$conf.int[1], corr$conf.int[2])
colnames(corr) <- c("P-value", "Estimate", "CI_lower", "CI_Upper")
write.csv(fst_average, "fst_average with distances.csv")
write.csv(corr, "Correlation.csv")



library(dplyr)
library(data.table)
library(stringr)
library(geosphere)
library(ggplot2)


setwd("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/Pf_Nov24_country_FST_25_11_24/")

## First stage is to read in the files 
file_list1 <- list.files(pattern="*.paired.snp_fst_df.csv")
df_list1_95 <- list()
df_list1_75 <- list()
n_list1 <- list()
n_list2 <- list()
for (i in 1:length(file_list1)){
    file <- read.csv(file_list1[i], header=TRUE)
    file1 <- file %>% filter(Fst>0.95)
    file2 <- file %>% filter(Fst>0.75)
    n1 <- nrow(file1)
    n2 <- nrow(file2)
    if (n1 > 0){
        df_list1_95[[i]] <- file1
        df_list1_95[[i]]$comparison <- file_list1[i]
    }
    if (n2 > 0){
        df_list1_75[[i]] <- file2
        df_list1_75[[i]]$comparison <- file_list1[i]
    }
    n_list1[[i]] <- data.frame(file_list1[i], n1)
    n_list2[[i]] <- data.frame(file_list1[i], n2)
}

df_list1_75_2 <- rbindlist(df_list1_75)
df_list1_95_2 <- rbindlist(df_list1_95)

frequencies <- read.table("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/e4a34f3d-a503-49b4-83af-2717f986bd39.frq", header=TRUE)
annotations <- read.csv("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/annotations_final.tsv", sep="\t", header=TRUE)

colnames(annotations) <- c("CHR", "POS", "REF", "ALT", "GENE")
colnames(frequencies) <- c("CHR", "POS", "REF", "ALT", "MAF", "NCHROBS")

df_list1_75_3 <- left_join(df_list1_75_2, annotations)
df_list1_75_3 <- left_join(df_list1_75_3, frequencies)
df_list1_75_4 <- df_list1_75_3 %>% filter(MAF<1)
df_list1_75_4 <- df_list1_75_4 %>% filter(!GENE %in% "")


df_list1_95_3 <- left_join(df_list1_95_2, annotations)
df_list1_95_3 <- left_join(df_list1_95_3, frequencies)
df_list1_95_4 <- df_list1_95_3 %>% filter(MAF<1)
df_list1_95_4 <- df_list1_95_4 %>% filter(!GENE %in% "")


write.csv(df_list1_75_4, "snps_0.75_region.fst")
write.csv(df_list1_95_4, "snps_0.95_region.fst")


n1_95 <- df_list1_95_4 %>% group_by(comparison) %>% summarise(n=n())
n2_75 <- df_list1_75_4 %>% group_by(comparison) %>% summarise(n=n())

write.csv(n1_95, "num_snps_0.95_country.fst")

write.csv(n2_75, "num_snps_0.75_country.fst")


## Let's get overall and draw a graph with correlation
setwd("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/Pf_Nov24_region_FST_25_11_24/")
file_list2 <- list.files(pattern="*.paired.hudson_fst_df.csv")
df_list2 <- list()
for (i in 1:length(file_list2)){
    file <- read.csv(file_list2[i], header=TRUE)
    df_list2[[i]] <- file
}

fst_average <- rbindlist(df_list2)
head(fst_average)
splits <- str_split_fixed(fst_average$Region, "_vs_", 2)
fst_average$Region1 <- splits[,1]
fst_average$Region2 <- splits[,2]
head(fst_average)

metadata <- read.csv("region_coords.csv", header=TRUE)
metadata$Region1 <- metadata$Region
metadata$Lat1 <- metadata$Latitude
metadata$Long1 <- metadata$Longitude
metadata$Lat2 <- metadata$Latitude
metadata$Long2 <- metadata$Longitude
metadata$Region2 <- metadata$Region
metadata1 <- metadata %>% select(Region1,  Lat1, Long1)
metadata1 <- unique(metadata1)
metadata2 <- metadata %>% select(Region2,  Lat2, Long2)
metadata2 <- unique(metadata2)
fst_average <- left_join(fst_average, metadata1)
fst_average <- left_join(fst_average, metadata2)
head(fst_average)

for ( i in 1:nrow(fst_average)){
fst_average$dist[i] <- distm(c( fst_average$Long1[i], fst_average$Lat1[i]), c(fst_average$Long2[i],fst_average$Lat2[i]), fun = distHaversine)
}
plot <- ggplot(fst_average, aes(x=dist, y=Hudson_Fst)) + geom_point() + theme_classic()
ggsave("plot_test.png",plot)
corr <- cor.test(as.numeric(as.character(fst_average$dist)), as.numeric(as.character(fst_average$Hudson_Fst)), use="complete.obs")
corr
str(corr)
corr <- data.frame(corr$p.value, corr$estimate, corr$conf.int[1], corr$conf.int[2])
colnames(corr) <- c("P-value", "Estimate", "CI_lower", "CI_Upper")
write.csv(fst_average, "fst_average with distances.csv")
write.csv(corr, "Correlation.csv")






setwd("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/Pf_Nov24_region_FST_25_11_24/")

## First stage is to read in the files 
file_list1 <- list.files(pattern="*.single.snp_fst_df.csv")
df_list1_95 <- list()
df_list1_75 <- list()
n_list1 <- list()
n_list2 <- list()
for (i in 1:length(file_list1)){
    file <- read.csv(file_list1[i], header=TRUE)
    file1 <- file %>% filter(Fst>0.95)
    file2 <- file %>% filter(Fst>0.75)
    n1 <- nrow(file1)
    n2 <- nrow(file2)
    if (n1 > 0){
        df_list1_95[[i]] <- file1
        df_list1_95[[i]]$comparison <- file_list1[i]
    }
    if (n2 > 0){
        df_list1_75[[i]] <- file2
        df_list1_75[[i]]$comparison <- file_list1[i]
    }
    n_list1[[i]] <- data.frame(file_list1[i], n1)
    n_list2[[i]] <- data.frame(file_list1[i], n2)
}

df_list1_75_2 <- rbindlist(df_list1_75)
df_list1_95_2 <- rbindlist(df_list1_95)

frequencies <- read.table("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/e4a34f3d-a503-49b4-83af-2717f986bd39.frq", header=TRUE)
annotations <- read.csv("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/annotations_final.tsv", sep="\t", header=TRUE)

colnames(annotations) <- c("CHR", "POS", "REF", "ALT", "GENE")
colnames(frequencies) <- c("CHR", "POS", "REF", "ALT", "MAF", "NCHROBS")

df_list1_75_3 <- left_join(df_list1_75_2, annotations)
df_list1_75_3 <- left_join(df_list1_75_3, frequencies)
df_list1_75_4 <- df_list1_75_3 %>% filter(MAF<1)
df_list1_75_4 <- df_list1_75_4 %>% filter(!GENE %in% "")


df_list1_95_3 <- left_join(df_list1_95_2, annotations)
df_list1_95_3 <- left_join(df_list1_95_3, frequencies)
df_list1_95_4 <- df_list1_95_3 %>% filter(MAF<1)
df_list1_95_4 <- df_list1_95_4 %>% filter(!GENE %in% "")


write.csv(df_list1_75_4, "snps_0.75_region_against_all.fst")
write.csv(df_list1_95_4, "snps_0.95_region_against_all.fst")


n1_95 <- df_list1_95_4 %>% group_by(comparison) %>% summarise(n=n())
n2_75 <- df_list1_75_4 %>% group_by(comparison) %>% summarise(n=n())



write.csv(n1_95, "num_snps_0.95_region_against_all.fst")

write.csv(n2_75, "num_snps_0.75_region_against_all.fst")




setwd("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/Pf_Nov24_country_FST_25_11_24/")

## First stage is to read in the files 
file_list1 <- list.files(pattern="*.single.snp_fst_df.csv")
df_list1_95 <- list()
df_list1_75 <- list()
n_list1 <- list()
n_list2 <- list()
for (i in 1:length(file_list1)){
    file <- read.csv(file_list1[i], header=TRUE)
    file1 <- file %>% filter(Fst>0.95)
    file2 <- file %>% filter(Fst>0.75)
    n1 <- nrow(file1)
    n2 <- nrow(file2)
    if (n1 > 0){
        df_list1_95[[i]] <- file1
        df_list1_95[[i]]$comparison <- file_list1[i]
    }
    if (n2 > 0){
        df_list1_75[[i]] <- file2
        df_list1_75[[i]]$comparison <- file_list1[i]
    }
    n_list1[[i]] <- data.frame(file_list1[i], n1)
    n_list2[[i]] <- data.frame(file_list1[i], n2)
}

df_list1_75_2 <- rbindlist(df_list1_75)
df_list1_95_2 <- rbindlist(df_list1_95)

frequencies <- read.table("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/e4a34f3d-a503-49b4-83af-2717f986bd39.frq", header=TRUE)
annotations <- read.csv("/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/annotations_final.tsv", sep="\t", header=TRUE)

colnames(annotations) <- c("CHR", "POS", "REF", "ALT", "GENE")
colnames(frequencies) <- c("CHR", "POS", "REF", "ALT", "MAF", "NCHROBS")

df_list1_75_3 <- left_join(df_list1_75_2, annotations)
df_list1_75_3 <- left_join(df_list1_75_3, frequencies)
df_list1_75_4 <- df_list1_75_3 %>% filter(MAF<1)
df_list1_75_4 <- df_list1_75_4 %>% filter(!GENE %in% "")


df_list1_95_3 <- left_join(df_list1_95_2, annotations)
df_list1_95_3 <- left_join(df_list1_95_3, frequencies)
df_list1_95_4 <- df_list1_95_3 %>% filter(MAF<1)
df_list1_95_4 <- df_list1_95_4 %>% filter(!GENE %in% "")


write.csv(df_list1_75_4, "snps_0.75_region_against_all.fst")
write.csv(df_list1_95_4, "snps_0.95_region_against_all.fst")


n1_95 <- df_list1_95_4 %>% group_by(comparison) %>% summarise(n=n())
n2_75 <- df_list1_75_4 %>% group_by(comparison) %>% summarise(n=n())



write.csv(n1_95, "num_snps_0.95_region_against_all.fst")

write.csv(n2_75, "num_snps_0.75_region_against_all.fst")

