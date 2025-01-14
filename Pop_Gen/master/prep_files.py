#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
from fastq2matrix import run_cmd
from tqdm import tqdm
import time
import pandas as pd

def make_anno(annotations):
  anno = pd.read_csv(annotations, sep="\t", header=None)
  gene_name1 = anno[4].str.split('|').str[1]
  gene_name2 = anno[4].str.split('|').str[2]
  gene_name = gene_name1 + " (" + gene_name2 + ")"
  anno_final = anno.iloc[:,0:4]
  anno_final['gene_name'] = gene_name
  anno_final.columns =  ['chr', 'pos', 'ref', 'alt', 'Gene_name']
  anno_final.to_csv("annotations_final.tsv", index=None, sep="\t")

def main(args):
    run_cmd("bcftools view -R %(coding_regions)s -Oz -o %(prefix)s_codingR.vcf.gz %(prefix)s.vcf.gz" % vars(args))
    run_cmd("Rscript %(path)s/master/split_meta.R --metadata %(metadata)s --suffix %(analysis)s --population %(population)s --wgs_id %(wgs_id)s" % vars(args))
    run_cmd("bcftools sort -Oz -o %(prefix)s_coding_sorted.vcf.gz --temp-dir . %(prefix)s_coding.vcf.gz" % vars(args))
    run_cmd("vcf_population_maf_filter.py --vcf %(prefix)s_codingR.vcf.gz --pop-file %(analysis)s.txt --maf %(maf)s --threads %(threads)s" % vars(args))
    run_cmd("vcf2matrix.py --vcf %(prefix)s_codingR.pop_maf_filt_%(maf)s.vcf.gz --na NA --threads %(threads)s " % vars(args))
    run_cmd("bcftools query -f '%%CHROM\\t%%POS\\t%%REF\\t%%ALT\\t%%BCSQ\\n' %(prefix)s_codingR.pop_maf_filt_%(maf)s.vcf.gz > annotations.tsv" % vars(args))
    make_anno("annotations.tsv")
    
    


parser = argparse.ArgumentParser(description='PopGen pipeline: prepare files',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--coding_regions',help='path to coding regions',required=True)
parser.add_argument('--prefix',help='vcf file prefix',required=True)
parser.add_argument('--path',type=str, help='path to PopGen script directory',required=True)
parser.add_argument('--metadata',type=str, help='metadata tsv file',required=True)
parser.add_argument('--wgs_id',type=str, help='name of column with wgs id',required=True)
parser.add_argument('--analysis',type=str, help='analysis name',required=True)
parser.add_argument('--population',type=str, help='column name containing populations of interest',required=True)
parser.add_argument('--maf',type=float,default=0.001, help='Minor allele freq cutoff', required=True)
parser.add_argument('--threads',default=4, type=int, help='Number of threads for parallel operations', required=True)


parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)

