#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
from fastq2matrix import run_cmd, nofile
from tqdm import tqdm
import time
import os.path
import pandas as pd


def main_moi(args):
    today = args.date 
    out_dir = args.analysis +"_MOI_" + today
    run_cmd("Rscript %(path)s/master/split_meta2.R --metadata %(metadata)s --suffix %(analysis)s --population %(population)s" % vars(args))
    run_cmd("ls *%(analysis)s.txt > %(analysis)s_population_files_list.txt" % vars(args))   
    run_cmd("cat %(analysis)s_population_files_list.txt | xargs -I {} -P %(parallel)s sh -c 'bcftools view -c1 -S {} --force-samples -Oz --threads 20 -o {}.genotyped.vcf.gz %(vcf)s' " % vars(args))
    run_cmd("echo 'Estimating MOI using Fws'")
    run_cmd("mkdir " + out_dir)
    run_cmd("ls *%(analysis)s.txt.genotyped.vcf.gz | sed 's/.genotyped.vcf.gz//g' > %(analysis)s_population_names_list.txt" % vars(args))     
    run_cmd("cat %(analysis)s_population_names_list.txt | xargs -I {} -P %(parallel)s sh -c 'Rscript %(path)s/moi/calculate_fws.R -p %(population)s -d . -f {}.genotyped.vcf.gz'" % vars(args))
    run_cmd("Rscript %(path)s/moi/merge_fws.R -p %(analysis)s -m %(metadata)s -w %(wgs_id)s" % vars(args))
    run_cmd("mv *%(analysis)s.txt* " % vars(args) + out_dir)
    run_cmd("mv *_moi_fws.tsv " + out_dir)
    run_cmd("mv *metadata_fws.tsv " + out_dir)
    
    

def main_admixture(args):
    today = args.date 
    out_dir = args.analysis +"_ADMIXTURE_" + today
    if "P_vivax" in args.species:
      run_cmd("plink2 --vcf %(vcf)s --set-missing-var-ids @:# --keep-allele-order --const-fid --allow-extra-chr --make-bed --out %(analysis)s_admixture" % vars(args))
      run_cmd("sed -ie 's/PvP01_//g; s/_v1//g; s/^0//g; s/API/15/g; s/MIT/16/g' %(analysis)s_admixture.bim" % vars(args))
      run_cmd("plink --indep-pairwise 50 10 0.1 --bfile  %(analysis)s_admixture --out %(analysis)s_admixture --allow-extra-chr" % vars(args))
      run_cmd("plink --bfile %(analysis)s_admixture --extract %(analysis)s_admixture.prune.in --make-bed --out %(analysis)s_admixture2" % vars(args))
      run_cmd("cat %(path)s/admixture/K_runs.txt | xargs -I {} -P %(parallel)s sh -c 'admixture --cv=10 -j8 --haploid='*' -s 12345 %(analysis)s_admixture2.bed {} | tee log{}.cv.haploid.seed12345.out'" % vars(args))
    elif "P_falciparum" in args.species:
      run_cmd("plink2 --vcf %(vcf)s --set-missing-var-ids @:# --keep-allele-order --const-fid --allow-extra-chr --make-bed --out %(analysis)s_admixture" % vars(args))
      run_cmd("sed -ie 's/Pf3D7_01_v3/1/g; s/Pf3D7_02_v3/2/g; s/Pf3D7_03_v3/3/g; s/Pf3D7_04_v3/4/g; s/Pf3D7_05_v3/5/g; s/Pf3D7_06_v3/6/g; s/Pf3D7_07_v3/7/g; s/Pf3D7_08_v3/8/g; s/Pf3D7_09_v3/9/g; s/Pf3D7_10_v3/10/g; s/Pf3D7_11_v3/11/g; s/Pf3D7_12_v3/12/g; s/Pf3D7_13_v3/13/g; s/Pf3D7_14_v3/14/g; s/Pf_M76611/16/g; s/Pf3D7_API_v3/17/g'  %(analysis)s_admixture.bim" % vars(args)) 
      run_cmd("plink --indep-pairwise 50 10 0.1 --bfile  %(analysis)s_admixture --out %(analysis)s_admixture --allow-extra-chr" % vars(args))
      run_cmd("plink --bfile %(analysis)s_admixture --extract %(analysis)s_admixture.prune.in --make-bed --out %(analysis)s_admixture2" % vars(args))
      run_cmd("cat %(path)s/admixture/K_runs.txt | xargs -I {} -P %(parallel)s sh -c 'admixture --cv=10 -j8 --haploid='*' -s 12345 %(analysis)s_admixture2.bed {} | tee log{}.cv.haploid.seed12345.out'" % vars(args))
    else:
      run_cmd("cat %(path)s/admixture/K_runs.txt | xargs -I {} -P %(parallel)s sh -c 'admixture --cv=10 -j8 --haploid='*' -s 12345 %(bed_file)s {} | tee log{}.cv.haploid.seed12345.out'" % vars(args))
    run_cmd("mkdir " + out_dir)
    run_cmd("mv %(analysis)s_admixture* " % vars(args) + out_dir)
    run_cmd("mv *.haploid.seed12345.out " + out_dir)
    
      
      
def main_selection(args):
  today = args.date  
  out_dir = args.analysis +"_SELECTION_" + today
  run_cmd("mkdir " + out_dir)
  if "prior" in args.prior_run:
    args.metadata = args.analysis +"_MOI_" + today + "/" + args.analysis + "_metadata_fws.tsv"
  if "paired" in args.comparison:
    meta = pd.read_csv(args.metadata, sep="\t")
    values = pd.unique(meta[args.population])  
    with open("pop_list.txt", "w") as txt_file:
      for line in values:
        txt_file.write(line + "\n")
    run_cmd("cat pop_list.txt | xargs -I {} -P %(parallel)s sh -c 'Rscript %(path)s/selection/prepare_input_rehh_per_category.R --path %(path)s " % vars(args) + "  -d ./" + out_dir + "  -b %(binary_matrix)s  --remove_chr %(rem_chr)s  -m %(metadata)s  --annotation %(annotation)s  -c {}  --label_category %(population)s  --label_fws fws  --fws_th %(fws_th)s  --label_id %(wgs_id)s  --forced_mixed --maf %(maf_th)s'" % vars(args))
    run_cmd("Rscript %(path)s/selection/calculate_rehh_metrics.R --path %(path)s " % vars(args) + "  -d ./" + out_dir + "  --prefix scanned_haplotypes  --remove_chr %(rem_chr)s  --list_category pop_list.txt  --annotation %(annotation)s  --gene_product %(product)s  --ihs_th %(ihs_th)s  --rsb_th %(rsb_th)s  --xpehh_th %(xpehh_th)s" % vars(args)) 
  if "single" in args.comparison: 
    values = args.pop_of_interest
    with open("pop_list.txt", "w") as txt_file:
      txt_file.write(values + "\n")
    run_cmd("Rscript %(path)s/selection/prepare_input_rehh_per_category.R --path %(path)s " % vars(args) + "  -d ./" + out_dir + "  -b %(binary_matrix)s  --remove_chr %(rem_chr)s -m %(metadata)s  --annotation %(annotation)s  -c %(pop_of_interest)s  --label_category %(population)s --label_fws fws  --fws_th %(fws_th)s  --label_id %(wgs_id)s  --forced_mixed --maf %(maf_th)s" % vars(args))
    run_cmd("Rscript %(path)s/selection/calculate_rehh_metrics.R --path %(path)s " % vars(args) + "  -d ./" + out_dir + "  --prefix scanned_haplotypes  --remove_chr %(rem_chr)s  --list_category pop_list.txt  --annotation %(annotation)s  --gene_product %(product)s  --ihs_th %(ihs_th)s  --rsb_th %(rsb_th)s  --xpehh_th %(xpehh_th)s" % vars(args)) 




def main_ibd(args):
  today = args.date 
  out_dir = args.analysis +"_IBD_" + today
  run_cmd("mkdir " + out_dir)
  if "prior" in args.prior_run:
    args.metadata = args.analysis +"_MOI_" + today + "/" + args.analysis + "_metadata_fws.tsv"
  if "paired" in args.comparison:
    meta = pd.read_csv(args.metadata, sep="\t")
    values = pd.unique(meta[args.population])  
    with open("pop_list.txt", "w") as txt_file:
      for line in values:
        txt_file.write(line + "\n")
    run_cmd("cat pop_list.txt | xargs -I {} -P %(parallel)s sh -c 'Rscript %(path)s/ibd/run_hmmIBD_per_category.R --path %(path)s " % vars(args) + "  -d . -b %(binary_matrix)s  --remove_chr %(rem_chr)s  -m %(metadata)s  -c {}  --label_category %(population)s  --label_fws fws  --fws_th %(fws_th)s  --label_id %(wgs_id)s  --maf %(maf_th)s'" % vars(args))
    run_cmd("Rscript %(path)s/ibd/summary_hmmIBD_results.R --path %(path)s -d . --list_category pop_list.txt  --gene_product %(product)s --legend ibd_matrix_hap_leg.tsv --ref_index %(ref_index)s --window_size %(window)s --quantile_cutoff %(quantile)s --remove_chr %(rem_chr)s" % vars(args))
  if "single" in args.comparison: 
    values = args.pop_of_interest
    with open("pop_list.txt", "w") as txt_file:
      txt_file.write(values + "\n")
    run_cmd("Rscript %(path)s/ibd/run_hmmIBD_per_category.R --path %(path)s   -d .  -b %(binary_matrix)s  --remove_chr %(rem_chr)s  -m %(metadata)s   -c %(pop_of_interest)s --label_category %(population)s  --label_fws fws  --fws_th %(fws_th)s  --label_id %(wgs_id)s --maf %(maf_th)s" % vars(args))
    run_cmd("Rscript %(path)s/ibd/summary_hmmIBD_results.R --path %(path)s " % vars(args) + "  -d . --list_category pop_list.txt  --gene_product %(product)s --legend ibd_matrix_hap_leg.tsv --ref_index %(ref_index)s --window_size %(window)s --quantile_cutoff %(quantile)s --remove_chr %(rem_chr)s" % vars(args)) 
  run_cmd("mv *hmmIBD* "  + out_dir)
  run_cmd("mv ibd* "  + out_dir)

def main_tess3r(args):
  today = args.date 
  out_dir = args.analysis +"_TESS3R_" + today
  run_cmd("mkdir " + out_dir)
  run_cmd("Rscript %(path)s/tess3r/prepare_tess3r_files.R --wd ./"   % vars(args) +  out_dir + "  --mat_bin %(binary_matrix)s --metadata %(metadata_csv)s" % vars(args))
  run_cmd("cat %(path)s/tess3r/K_runs.txt | xargs -I {} sh -c 'Rscript %(path)s/tess3r/run_tess3r_general.R "  % vars(args) + "  -b ./" + out_dir + "/tess3r_matbinT.mis.rds --coord ./" + out_dir + "/tess3r_coordsM.rds --K {} --rep %(rep)s --threads %(threads)s' " % vars(args))
  run_cmd("mv *_tess3r.rds " % vars(args) + out_dir)
  

def main_all(args):
  if "run" in args.moi:
    main_moi(args)
  if "run" in args.adm:
    main_admixture(args)
  if "run" in args.selection:
    main_selection(args)
  if "run" in args.ibd:
    main_ibd(args)
    

parser = argparse.ArgumentParser(description='PopGen pipeline: Analysis',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
subparsers = parser.add_subparsers(help="Task to perform")

parser_sub = subparsers.add_parser('ALL', help='Estimate Fws for MOI', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--analysis',type=str, help='analysis name, different from prep files arg',required=True)
parser_sub.add_argument('--annotation',type=str, help='path to annotation file', required=True)
parser_sub.add_argument('--bed_file',type=str, help='bed file if species is other',required=False)
parser_sub.add_argument('--binary_matrix',type=str, help='binary matrix to input', required=True)
parser_sub.add_argument('--comparison',type=str, help='paired or single', required=True)
parser_sub.add_argument('--date',default="today", type=str, help='date', required=True)
parser_sub.add_argument('--fws_th',type=float, default=0.95, help='fws threshold', required=False)
parser_sub.add_argument('--ihs_th',type=int, default=4, help='ihs threshold', required=False)
parser_sub.add_argument('--maf_th',type=float, default=0.01, help='maf threshold', required=False)
parser_sub.add_argument('--metadata',type=str, help='metadata tsv file',required=True)
parser_sub.add_argument('--parallel',default=5, type=int, help='Number of threads for parallel operations', required=True)
parser_sub.add_argument('--path',type=str, help='path to PopGen script directory',required=True)
parser_sub.add_argument('--pop_of_interest',type=str, help='single population to input', required=False)
parser_sub.add_argument('--population',type=str, help='column name containing populations of interest',required=True)
parser_sub.add_argument('--prior_run',type=str, default='prior', help='if the pipeline has been run previously ie is there fws metadata file',required=False)
parser_sub.add_argument('--product',type=str, help='path to product file', required=True)
parser_sub.add_argument('--quantile',type=float, default=0.95, help='fws threshold', required=False)
parser_sub.add_argument('--ref_index',type=str, help='path to ref index', required=True)
parser_sub.add_argument('--rem_chr', type=str, default='Pf3D7_API_v3,Pf_M76611', help='chr to remove', required=False)
parser_sub.add_argument('--rsb_th',type=int, default=5, help='rsb threshold', required=False)
parser_sub.add_argument('--species',type=str, help='P_vivax or P_falciparum',required=True)
parser_sub.add_argument('--vcf',type=str, help='input vcf',required=True)
parser_sub.add_argument('--wgs_id',type=str, help='name of column with wgs id',required=True)
parser_sub.add_argument('--window',type=int, default=10000, help='ihs threshold', required=False)
parser_sub.add_argument('--xpehh_th',type=int, default=5, help='xpehh threshold', required=False)
parser_sub.add_argument('--moi',type=str, default='run', help='run moi?', required=False)
parser_sub.add_argument('--adm',type=str, default='run', help='run adm?', required=False)
parser_sub.add_argument('--selection',type=str, default='run', help='run_selection', required=False)
parser_sub.add_argument('--ibd',type=str, default='run', help='run ibd?', required=False)
parser_sub.set_defaults(func=main_all)


parser_sub = subparsers.add_parser('MOI', help='Estimate Fws for MOI', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--path',type=str, help='path to PopGen script directory',required=True)
parser_sub.add_argument('--vcf',type=str, help='input vcf',required=True)
parser_sub.add_argument('--metadata',type=str, help='metadata tsv file',required=True)
parser_sub.add_argument('--wgs_id',type=str, help='name of column with wgs id',required=True)
parser_sub.add_argument('--analysis',type=str, help='analysis name, different from prep files arg',required=True)
parser_sub.add_argument('--population',type=str, help='column name containing populations of interest',required=True)
parser_sub.add_argument('--parallel',default=5, type=int, help='Number of threads for parallel operations', required=True)
parser_sub.add_argument('--date',default="today", type=str, help='date', required=True)
parser_sub.set_defaults(func=main_moi)


parser_sub = subparsers.add_parser('ADMIXTURE', help='Admixture analysis for inferring ancestry', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--path',type=str, help='path to PopGen script directory',required=True)
parser_sub.add_argument('--vcf',type=str, help='input vcf',required=True)
parser_sub.add_argument('--species',type=str, help='P_vivax or P_falciparum',required=True)
parser_sub.add_argument('--bed_file',type=str, help='bed file if species is other',required=False)
parser_sub.add_argument('--analysis',type=str, help='analysis name, different from prep files arg',required=True)
parser_sub.add_argument('--parallel',default=5, type=int, help='Number of threads for parallel operations', required=True)
parser_sub.add_argument('--date',default="today", type=str, help='date', required=True)
parser_sub.set_defaults(func=main_admixture)


parser_sub = subparsers.add_parser('SELECTION', help='Calculates iHS, XPEHH, rsb selection metrics', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--path',type=str, help='path to PopGen script directory',required=True)
parser_sub.add_argument('--prior_run',type=str, default='prior', help='if the pipeline has been run previously ie is there fws metadata file',required=False)
parser_sub.add_argument('--metadata',type=str, help='metadata tsv file',required=True)
parser_sub.add_argument('--wgs_id',type=str, help='name of column with wgs id',required=True)
parser_sub.add_argument('--analysis',type=str, help='analysis name, different from prep files arg',required=True)
parser_sub.add_argument('--population',type=str, help='column name containing populations of interest',required=True)
parser_sub.add_argument('--parallel',default=5, type=int, help='Number of threads for parallel operations', required=True)
parser_sub.add_argument('--rem_chr', type=str, default='Pf3D7_API_v3,Pf_M76611', help='chr to remove', required=False)
parser_sub.add_argument('--annotation',type=str, help='path to annotation file', required=True)
parser_sub.add_argument('--product',type=str, help='path to product file', required=True)
parser_sub.add_argument('--fws_th',type=float, default=0.95, help='fws threshold', required=False)
parser_sub.add_argument('--ihs_th',type=int, default=4, help='ihs threshold', required=False)
parser_sub.add_argument('--rsb_th',type=int, default=5, help='rsb threshold', required=False)
parser_sub.add_argument('--xpehh_th',type=int, default=5, help='xpehh threshold', required=False)
parser_sub.add_argument('--maf_th',type=float, default=0.01, help='maf threshold', required=False)
parser_sub.add_argument('--binary_matrix',type=str, help='binary matrix to input', required=True)
parser_sub.add_argument('--comparison',type=str, help='paired or single', required=True)
parser_sub.add_argument('--pop_of_interest',type=str, help='single population to input', required=False)
parser_sub.add_argument('--date',default="today", type=str, help='date', required=True)
parser_sub.set_defaults(func=main_selection)


parser_sub = subparsers.add_parser('IBD', help='Identity by descent analysis', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--path',type=str, help='path to PopGen script directory',required=True)
parser_sub.add_argument('--prior_run',type=str, default='prior', help='if the pipeline has been run previously ie is there fws metadata file',required=False)
parser_sub.add_argument('--metadata',type=str, help='metadata tsv file',required=True)
parser_sub.add_argument('--wgs_id',type=str, help='name of column with wgs id',required=True)
parser_sub.add_argument('--analysis',type=str, help='analysis name, different from prep files arg',required=True)
parser_sub.add_argument('--population',type=str, help='column name containing populations of interest',required=True)
parser_sub.add_argument('--parallel',default=5, type=int, help='Number of threads for parallel operations', required=True)
parser_sub.add_argument('--rem_chr', type=str, default='Pf3D7_API_v3,Pf_M76611', help='chr to remove', required=False)
parser_sub.add_argument('--product',type=str, help='path to product file', required=True)
parser_sub.add_argument('--fws_th',type=float, default=0.95, help='fws threshold', required=False)
parser_sub.add_argument('--quantile',type=float, default=0.95, help='fws threshold', required=False)
parser_sub.add_argument('--window',type=int, default=10000, help='ihs threshold', required=False)
parser_sub.add_argument('--maf_th',type=float, default=0.01, help='maf threshold', required=False)
parser_sub.add_argument('--binary_matrix',type=str, help='binary matrix to input', required=True)
parser_sub.add_argument('--ref_index',type=str, help='path to ref index', required=True)
parser_sub.add_argument('--comparison',type=str, help='paired or single', required=True)
parser_sub.add_argument('--pop_of_interest',type=str, help='single population to input', required=False)
parser_sub.add_argument('--date',default="today", type=str, help='date', required=True)
parser_sub.set_defaults(func=main_ibd)


parser_sub = subparsers.add_parser('TESS3R', help='Geographic patterns in admixture', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser_sub.add_argument('--analysis',type=str, help='analysis name, different from prep files arg',required=True)
parser_sub.add_argument('--path',type=str, help='path to PopGen script directory',required=True)
parser_sub.add_argument('--metadata_csv',type=str, help='metadata csv file, with colnames: wgs_id, Site_latitude, Site_longitude, Site, Country, Region',required=True)
parser_sub.add_argument('--binary_matrix',type=str, help='binary matrix to input', required=True)
parser_sub.add_argument('--rep',type=float, default=50, help='reps', required=False)
parser_sub.add_argument('--threads',type=float, default=10, help='number of threads', required=False)
parser_sub.add_argument('--date',default="today", type=str, help='date', required=True)
parser_sub.set_defaults(func=main_tess3r)



args = parser.parse_args()
if vars(args) == {}:
    parser.print_help(sys.stderr)
else:
    args.func(args)
