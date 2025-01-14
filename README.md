# PopGen :boom:
This repository contains a master script/ code for performing **Population Genetics Analysis** e.g. for *Plasmodium falciparum*!
Let's start by making a repository containing all of the code, packages and software we will need for the analysis...


```
conda create -n pop_stat
conda activate pop_stat
conda install python=3.8 scikit-allel zarr r-base matplotlib seaborn



conda create -n pop_gen
conda activate pop_gen

## Lot's of installations...
install fastq2matrix
conda install python=3.7 bwa samtools bcftools parallel datamash gatk4=4.1.4.1 delly tqdm trimmomatic minimap2 biopython bedtools r-ggplot2 iqtree plink
cd /path/to/fastq2matrix
git pull
python setup.py install


conda install r-base r-BiocManager r-Remotes r-knitr r-lme4 r-mitml r-nloptr r-mice r-mice r-logistf r-devtools
conda install cmake

## MOI in R
BiocManager::install("SeqVarTools")
BiocManager::install("SeqArray", force=TRUE)
BiocManager::install("bahlolab/moimix") # ran with no updates, quite problematic to install
install.packages("optparse")
install.packages("readr")
install.packages("dplyr")

## update bcftools due to dependency clash
conda update bcftools

## install admixture
conda install bioconda::admixture bioconda::plink2


## selection
conda install r-base r-data.table r-stringr r-ggplot2 r-ggrepel r-tidyr --channel conda-forge
## open R
install.packages("rehh")


##IBD

cd ~/PopGen # or your path
git clone https://github.com/glipsnort/hmmIBD.git
cc -o hmmIBD -O3 -Wall hmmIBD.c -lm


## tess3r
# Open R
devtools::install_github("bcm-uga/TESS3_encho_sen")
```
## Prepare Files
```
conda activate pop_gen

python "/mnt/storage13/nbillows/Pop_Gen/master/prep_files.py" --coding_regions "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pf3D7_R.coding.regions.bed" --prefix "dummy_data.2024_10_29.filt.csq.bi.GT.miss0.4.vqslod.filt.snps" --path "/mnt/storage13/nbillows/Pop_Gen/" --metadata "dummy.tsv" --analysis "country_test" --wgs_id "wgs_id" --population "Country" --maf 0.001 --threads 8
```


## Run ALL Steps PopGen
```
python "/mnt/storage13/nbillows/Pop_Gen/master/PopGen.py" ALL --path "/mnt/storage13/nbillows/Pop_Gen/" --binary_matrix "/mnt/storage13/nbillows/Pf_09_24/dummy/dummy_data.2024_10_29.filt.bi.GT.miss0.4.vqslod.filt.snps_coding_sorted.pop_maf_filt_0.001.mat.bin" --analysis "checking_all"  --parallel 5 --metadata "dummy.tsv" --annotation "/mnt/storage13/nbillows/Pf_09_24/dummy/annotations_final.tsv" --product "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/pf_genome_product_v3.tsv" --wgs_id "wgs_id" --population "Country" --date 2024_11_07 --comparison paired --ref_index "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta.fai" --species P_falciparum --vcf "/mnt/storage13/nbillows/Pf_09_24/dummy/dummy_data.2024_10_29.filt.bi.GT.miss0.4.vqslod.filt.snps_coding_sorted.pop_maf_filt_0.001.vcf.gz"
```
## Just MOI
```
python "/mnt/storage13/nbillows/Pop_Gen/master/PopGen.py" MOI --path "/mnt/storage13/nbillows/Pop_Gen/" --vcf "/mnt/storage13/nbillows/Pf_09_24/dummy/dummy_data.2024_10_29.filt.bi.GT.miss0.4.vqslod.filt.snps_coding_sorted.pop_maf_filt_0.001.vcf.gz" --metadata "dummy.tsv" --wgs_id "wgs_id" --analysis "check" --population "Country" --parallel 5
```
## Just Admixture
```
### Note args for species are P_falciparum or P_vivax. For any other species, please provide a bed file using the option --bed_file- should have chr integers.
python "/mnt/storage13/nbillows/Pop_Gen/master/PopGen.py" ADMIXTURE --path "/mnt/storage13/nbillows/Pop_Gen/" --vcf "/mnt/storage13/nbillows/Pf_09_24/dummy/dummy_data.2024_10_29.filt.bi.GT.miss0.4.vqslod.filt.snps_coding_sorted.pop_maf_filt_0.001.vcf.gz" --analysis "check" --species "P_falciparum" --parallel 5
## handy if using own bed file:
plink --vcf <vcf_file> --set-missing-var-ids @:# --keep-allele-order --const-fid --allow-extra-chr --make-bed --out <prefix>
sed -ie 's/PvP01_//g; s/_v1//g; s/^0//g; s/API/15/g; s/MIT/16/g' <prefix>.bim
```
## Just Selection
```
python "/mnt/storage13/nbillows/Pop_Gen/master/PopGen.py" SELECTION --path "/mnt/storage13/nbillows/Pop_Gen/" --binary_matrix "/mnt/storage13/nbillows/Pf_09_24/dummy/dummy_data.2024_10_29.filt.bi.GT.miss0.4.vqslod.filt.snps_coding_sorted.pop_maf_filt_0.001.mat.bin" --analysis "check"  --parallel 5 --metadata "dummy.tsv" --annotation "/mnt/storage13/nbillows/Pf_09_24/dummy/annotations_final.tsv" --product "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/pf_genome_product_v3.tsv" --wgs_id "wgs_id" --population "Country" --date 2024_11_06 --comparison paired
```
## Just IBD
```
python "/mnt/storage13/nbillows/Pop_Gen/master/PopGen.py" IBD --path "/mnt/storage13/nbillows/Pop_Gen/" --binary_matrix "/mnt/storage13/nbillows/Pf_09_24/dummy/dummy_data.2024_10_29.filt.bi.GT.miss0.4.vqslod.filt.snps_coding_sorted.pop_maf_filt_0.001.mat.bin" --analysis "check"  --parallel 5 --metadata "dummy.tsv" --product "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/pf_genome_product_v3.tsv" --wgs_id "wgs_id" --population "Country" --ref_index "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pfalciparum.genome.fasta.fai" --date 2024_11_06 --comparison paired
```
## Just tess3r
```
python "/mnt/storage13/nbillows/Pop_Gen/master/PopGen.py" TESS3R --path "/mnt/storage13/nbillows/Pop_Gen/" --binary_matrix "/mnt/storage13/nbillows/Pf_09_24/dummy/dummy_data.2024_10_29.filt.bi.GT.miss0.4.vqslod.filt.snps_coding_sorted.pop_maf_filt_0.001.mat.bin" --metadata_csv "dummy_meta.csv" --date 2024_11_06 --analysis check
```

## Population Statistics Analysis
```
conda activate pop_stat
python "/mnt/storage13/nbillows/Pop_Gen/master/PopStat.py" --path "/mnt/storage13/nbillows/Pop_Gen/" --vcf "/m
nt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/Pf_dataset_Oct24_filt.csq.bi.GT.miss0.4.vqslod.filt.snps.vcf.gz" --vcf_coding "/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/Pf_dataset_Oct24_filt.csq.bi.GT.miss0.4.vqslod.filt.snps_codingR.vcf.gz" --comparison both --population Region --make_zarr True --metadata "/mnt/storage13/nbillows/Pf_09_24/Pfalciparum_09_24_v2/analysis_09_24_v2/Pf_Nov24_metadata_complete.csv" --do_pca True --prefix Pf_Nov24 --date 25_11_24
```

