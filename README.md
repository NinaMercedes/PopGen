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


## selection
conda install r-base r-data.table r-stringr r-ggplot2 r-ggrepel r-tidyr --channel conda-forge
## open R
install.packages("rehh")


##IBD

cd ~/PopGen # or your path
git clone https://github.com/glipsnort/hmmIBD.git
cc -o hmmIBD -O3 -Wall hmmIBD.c -lm

```

```
conda activate pop_gen

python "/mnt/storage13/nbillows/Pop_Gen/master/prep_files.py" --coding_regions "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/Pf3D7_R.coding.regions.bed" --prefix "dummy_data.2024_10_29.filt.csq.bi.GT.miss0.4.vqslod.filt.snps" --path "/mnt/storage13/nbillows/Pop_Gen/" --metadata "dummy.tsv" --analysis "country_test" --wgs_id "wgs_id" --population "Country" --maf 0.001 --threads 8

python "/mnt/storage13/nbillows/Pop_Gen/master/PopGen.py" MOI --path "/mnt/storage13/nbillows/Pop_Gen/" --vcf "/mnt/storage13/nbillows/Pf_09_24/dummy/dummy_data.2024_10_29.filt.bi.GT.miss0.4.vqslod.filt.snps_coding_sorted.pop_maf_filt_0.001.vcf.gz" --metadata "dummy.tsv" --wgs_id "wgs_id" --analysis "check" --population "Country" --parallel 5

### Note args for species are P_falciparum or P_vivax. For any other species, please provide a bed file using the option --bed_file- should have chr integers.
python "/mnt/storage13/nbillows/Pop_Gen/master/PopGen.py" ADMIXTURE --path "/mnt/storage13/nbillows/Pop_Gen/" --vcf "/mnt/storage13/nbillows/Pf_09_24/dummy/dummy_data.2024_10_29.filt.bi.GT.miss0.4.vqslod.filt.snps_coding_sorted.pop_maf_filt_0.001.vcf.gz" --analysis "check" --species "P_falciparum" --parallel 5

## handy if using own bed file:
plink --vcf <vcf_file> --set-missing-var-ids @:# --keep-allele-order --const-fid --allow-extra-chr --make-bed --out <prefix>
sed -ie 's/PvP01_//g; s/_v1//g; s/^0//g; s/API/15/g; s/MIT/16/g' <prefix>.bim



python "/mnt/storage13/nbillows/Pop_Gen/master/PopGen.py" SELECTION --path "/mnt/storage13/nbillows/Pop_Gen/" --binary_matrix "/mnt/storage13/nbillows/Pf_09_24/dummy/dummy_data.2024_10_29.filt.bi.GT.miss0.4.vqslod.filt.snps_coding_sorted.pop_maf_filt_0.001.mat.bin" --analysis "check"  --parallel 5 --metadata "dummy.tsv" --annotation "/mnt/storage13/nbillows/Pf_09_24/dummy/annotations_final.tsv" --product "/mnt/storage13/nbillows/Pf_09_24/Pf3D7_v3/pf_genome_product_v3.tsv" --wgs_id "wgs_id" --population "Country"




parser.set_defaults(func=main)
```
"/mnt/storage13/nbillows/Pop_Gen/"  --path  --vcf VCF --metadata METADATA --wgs_id
                 WGS_ID --analysis ANALYSIS --population POPULATION --maf MAF
                 --threads THREADS

Going to use fastq2 matrix to get matrix and also do maf filtering first....





parser_sub = subparsers.add_parser('SELECTION', help='Estimate Fws for MOI', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--path',type=str, help='path to PopGen script directory',required=True)
parser_sub.add_argument('--prior',type=str, default='prior', help='if the pipeline has been run previously ie is there fws metadata file',required=True)
parser_sub.add_argument('--metadata',type=str, help='metadata tsv file',required=True)
parser_sub.add_argument('--wgs_id',type=str, help='name of column with wgs id',required=True)
parser_sub.add_argument('--analysis',type=str, help='analysis name, different from prep files arg',required=True)
parser_sub.add_argument('--population',type=str, help='column name containing populations of interest',required=True)
parser_sub.add_argument('--parallel',default=5, type=int, help='Number of threads for parallel operations', required=True)
parser_sub.add_argument('--rem_chr', type=str, default='Pf3D7_API_v3,Pf_M76611', help='chr to remove', required=True)
parser_sub.add_argument('--annotation',type=str, help='path to annotation file', required=True)
parser_sub.add_argument('--product',type=str, help='path to product file', required=True)
parser_sub.add_argument('--fws_th',type=int, default=5, help='fws threshold', required=True)
parser_sub.add_argument('--ihs_th',type=int, default=4 help='fws threshold', required=True)
parser_sub.add_argument('--rsb_th',type=int, default=5, help='fws threshold', required=True)
parser_sub.add_argument('--xpehh_th',type=int, default=5, help='fws threshold', required=True)
parser_sub.add_argument('--binary_matrix',type=str, help='binary matrix to input', required=True)


DR haploytpes etcc
