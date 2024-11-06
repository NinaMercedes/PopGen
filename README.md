# PopGen :boom:
This repository contains a master script/ code for performing **Population Genetics Analysis** e.g. for *Plasmodium falciparum*!
Let's start by making a repository containing all of the code, packages and software we will need for the analysis...
```
conda create -n pop_gen
conda activate pop_gen

## Lot's of installations...
install fastq2matrix
conda install python=3.7 bwa samtools bcftools parallel datamash gatk4=4.1.4.1 delly tqdm trimmomatic minimap2 biopython bedtools r-ggplot2 iqtree plink
cd /path/to/fastq2matrix
git pull
python setup.py install


conda install r-BiocManager r-Remotes r-knitr r-lme4 r-mitml r-nloptr r-mice r-mice r-logistf
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
```

```
conda activate pop_gen

python "/mnt/storage13/nbillows/Pop_Gen/master/prep_files.py" --coding_regions "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pf3D7_R.coding.regions.bed" --prefix "dummy_data.2024_10_29.filt.bi.GT.miss0.4.vqslod.filt.snps" --path "/mnt/storage13/nbillows/Pop_Gen/" --metadata "dummy.tsv" --analysis "country_test" --wgs_id "wgs_id" --population "Country" --maf 0.001 --threads 8

python "/mnt/storage13/nbillows/Pop_Gen/master/PopGen.py" MOI --path "/mnt/storage13/nbillows/Pop_Gen/" --vcf "/mnt/storage13/nbillows/Pf_09_24/dummy/dummy_data.2024_10_29.filt.bi.GT.miss0.4.vqslod.filt.snps_coding_sorted.pop_maf_filt_0.001.vcf.gz" --metadata "dummy.tsv" --wgs_id "wgs_id" --analysis "check" --population "Country" --parallel 5

### Note args for species are P_falciparum or P_vivax. For any other species, please provide a bed file using the option --bed_file- should have chr integers.
python "/mnt/storage13/nbillows/Pop_Gen/master/PopGen.py" ADMIXTURE --path "/mnt/storage13/nbillows/Pop_Gen/" --vcf "/mnt/storage13/nbillows/Pf_09_24/dummy/dummy_data.2024_10_29.filt.bi.GT.miss0.4.vqslod.filt.snps_coding_sorted.pop_maf_filt_0.001.vcf.gz" --analysis "check" --species "P_falciparum" --parallel 5

## handy if using own bed file:
plink --vcf <vcf_file> --set-missing-var-ids @:# --keep-allele-order --const-fid --allow-extra-chr --make-bed --out <prefix>
sed -ie 's/PvP01_//g; s/_v1//g; s/^0//g; s/API/15/g; s/MIT/16/g' <prefix>.bim



parser.set_defaults(func=main)
```
"/mnt/storage13/nbillows/Pop_Gen/"  --path  --vcf VCF --metadata METADATA --wgs_id
                 WGS_ID --analysis ANALYSIS --population POPULATION --maf MAF
                 --threads THREADS

Going to use fastq2 matrix to get matrix and also do maf filtering first....


DR haploytpes etcc
