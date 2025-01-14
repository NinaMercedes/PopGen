#! /usr/bin/env python

import sys
import os
import argparse
import subprocess
import itertools
import pandas as pd
import numpy as np


### Functions
def clean_data(metadata):
  meta = pd.read_csv(metadata, sep='\t')
  return(meta)


def region_fst_pair_analysis(meta, pop, out_dir, thresh):
  geo_r = meta[pop].unique()
  geo_r  = [x for x in geo_r if str(x) != 'nan'] 
  high_fst_snps =[]
  fst_region =[]
  for comb in itertools.combinations(geo_r,2): #check combs vs permutations
    pop1 = comb[0]
    pop2 = comb[1]
    filename = out_dir + "/" + comb[0] + "_vs_" + comb[1] + ".paired.snp_fst_df.csv"
    region_fst = pd.read_csv(filename,header=0)
    region_fst[region_fst['Fst']<0] = 0
    mean_val = np.mean(region_fst['Fst'])
    fst_80 = region_fst[region_fst['Fst']>thresh]
    fst_80["Population 1"] = pop1
    fst_80["Population 2"] = pop2
    num_80 = len(fst_80)
    fst_region.append(pd.DataFrame({'Population 1': pop1, 'Population 2': pop2, 'Fst>threshold':num_80}, index=[0]))
    high_fst_snps.append(fst_80)
  fst_region_df = pd.concat(fst_region, ignore_index=True)
  high_fst_snps_df = pd.concat(high_fst_snps, ignore_index=True)
  return(fst_region_df, high_fst_snps_df)
  
def region_fst_one_analysis(meta, pop, out_dir, thresh):
  geo_r = meta[pop].unique()
  geo_r  = [x for x in geo_r if str(x) != 'nan']
  high_fst_snps =[]
  fst_region =[]  
  for region in geo_r:
    filename = out_dir + "/" + region + "_vs_all" + "single.paired.snp_fst_df.csv"
    region_fst = pd.read_csv(filename,header=0)
    fst_80 = region_fst[region_fst['Fst']>thresh]
    fst_80["Region"] = region
    num_80 = len(fst_80)
    fst_region.append(pd.DataFrame({'Region': region,  'Fst>thresh':num_80}, index=[0]))
    high_fst_snps.append(fst_80)
  fst_region_df = pd.concat(fst_region, ignore_index=True)
  high_fst_snps_df = pd.concat(high_fst_snps, ignore_index=True)
  return(fst_region_df, high_fst_snps_df)


#UP TO THIS BIT
def main(args):
  #params = {"metadata": args.metadata}
  meta_clean = clean_data(args.metadata)
  if "paired" or "both" in args.comparison:
    fst_region_df1, high_fst_snps_df1 = region_fst_pair_analysis(meta_clean, args.pop, args.out_dir, args.thresh)
    fst_region_df1.to_csv(args.out_dir + "/" + "fst_region_pairs.csv")
    high_fst_snps_df1.to_csv(args.out_dir + "/" + "high_fst_snps_region_pairs.csv")
  if "one_against_all" or "both" in args.comparison:
    fst_region_df2, high_fst_snps_df2 = region_fst_one_analysis(meta_clean, args.pop, args.out_dir, args.thresh)
    fst_region_df2.to_csv(args.out_dir + "/" + "fst_region_single.csv")
    high_fst_snps_df2.to_csv(args.out_dir + "/" + "high_fst_snps_region_single.csv")

parser = argparse.ArgumentParser(description='Fst statistics analysis',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--metadata',help='csv file containing metadata',required=True)
parser.add_argument('--comparison',help='paired or one_against_all or both',required=True)
parser.add_argument('--thresh', default=80, help='paired or one_against_all or both',required=True)
parser.add_argument('--out_dir',help='out_dir',required=True)
parser.add_argument('--pop',help='population column',required=True)
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)
