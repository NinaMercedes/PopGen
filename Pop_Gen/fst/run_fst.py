#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
from tqdm import tqdm
import time
import os.path
import pandas as pd
import numpy as np
import scipy
import h5py
import allel
from allel import chunked
import pandas as pd
import zarr
import itertools
print('scikit-allel', allel.__version__)


def main(args):
  #today = args.date 
  #out_dir = args.prefix +"_FST_" + args.comparison + "_" + today
  #subprocess.run("mkdir " + out_dir, shell=True)
  #out_zarr = out_dir + "/" + args.prefix + ".zarr"
  #allel.vcf_to_zarr(args.vcf, out_zarr) #'/mnt/storage12/nbillows/Projects/Malaria/20k_finale_QC_pass.bi.GT.miss0.4.vqslod.filt.snps.vcf.gz'
  callset = zarr.open(args.zarr, mode='r')
  sample_names = callset['samples'][:]
#testing
  pos_index = allel.ChromPosIndex(callset['variants/CHROM'][:], callset['variants/POS'][:])
  genotype_all = allel.GenotypeChunkedArray(callset['calldata']['GT'])
  genotype_all
  df_samples = pd.read_csv(args.meta, sep="\t")
  df_samples['domain'] = df_samples[args.wgs_id]
  df_samples =pd.DataFrame({'domain':sample_names}).merge(df_samples, on='domain', how='left')
  #df_samples = df_samples.dropna(subset=[args.population])
  geo_r = df_samples[args.population].unique() 
  geo_r  = [x for x in geo_r if str(x) != 'nan']  
  if "paired" or "both" in args.comparison:
    for comb in itertools.combinations(geo_r,2): #check combs vs permutations
      pop1 = comb[0]
      pop2 = comb[1]
      n_samples_pop1 = np.count_nonzero(df_samples[args.population] == pop1)
      n_samples_pop2 = np.count_nonzero(df_samples[args.population] == pop2)
      print(pop1, n_samples_pop1, pop2, n_samples_pop2)     
      
      subpops = {
          pop1: df_samples[df_samples[args.population] == pop1].index,
          pop2: df_samples[df_samples[args.population] == pop2].index,
      }
      # allele counts
      acs = genotype_all.count_alleles_subpops(subpops)
      acs
      
      acu = allel.AlleleCountsArray(acs[pop1][:] + acs[pop2][:])
      flt = acu.is_segregating() & (acu.max_allele() == 1)
      print('retaining', np.count_nonzero(flt), 'SNPs')
      
      ac1 = allel.AlleleCountsArray(acs[pop1].compress(flt, axis=0)[:, :2])
      ac2 = allel.AlleleCountsArray(acs[pop2].compress(flt, axis=0)[:, :2])
      genotype = genotype_all.compress(flt, axis=0)
      genotype
      # sample indices for population 1
      pop1_idx = subpops[pop1]
      # sample indices for population 2
      pop2_idx = subpops[pop2]
      num, den = allel.hudson_fst(ac1, ac2)
      snp_fst_hudson = num / den
      snp_fst_hudson
      fst_hudson, se_hudson, vb_hudson, _ = allel.average_hudson_fst(ac1, ac2, blen=10000)
      print('%.04f +/- %.04f (Hudson)' % (fst_hudson, se_hudson))
    
      ## now make some dfs!!
      regions = pop1 + "_vs_" + pop2
      f_name1=  args.out_dir + "/"+ regions + ".paired.snp_fst_df.csv"
      f_name2= args.out_dir + "/" + regions + ".paired.hudson_fst_df.csv"
      snp_fst_df = pd.DataFrame({'CHR': callset['variants/CHROM'][:], 'POS': callset['variants/POS'][:] , 'retained':flt})
      snp_fst_df_flt = snp_fst_df[(snp_fst_df['retained']==True)]
      snp_fst_df_flt['Fst'] = snp_fst_hudson
      snp_fst_df_flt.to_csv(f_name1)
      fst_hudson_df = pd.DataFrame({'Region':regions, 'Hudson_Fst': fst_hudson, 'SE': se_hudson, 'n_retained': [np.count_nonzero(flt)]})
      fst_hudson_df.to_csv(f_name2)
  
  if "one_against_all" or "both" in args.comparison:
    for pop in geo_r: #check combs vs permutations
      pop1 = pop
      pop2 = "Other"
      n_samples_pop1 = np.count_nonzero(df_samples[args.population] == pop)
      n_samples_pop2 = np.count_nonzero(df_samples[args.population] != pop)
      print(pop1, n_samples_pop1, "other", n_samples_pop2)
      
      subpops = {
      pop1: df_samples[df_samples[args.population] == pop].index,
      pop2: df_samples[df_samples[args.population] != pop].index,
      }
      # allele counts
      acs = genotype_all.count_alleles_subpops(subpops)
      acs
      
      acu = allel.AlleleCountsArray(acs[pop1][:] + acs[pop2][:])
      flt = acu.is_segregating() & (acu.max_allele() == 1)
      print('retaining', np.count_nonzero(flt), 'SNPs')
      
      ac1 = allel.AlleleCountsArray(acs[pop1].compress(flt, axis=0)[:, :2])
      ac2 = allel.AlleleCountsArray(acs[pop2].compress(flt, axis=0)[:, :2])
      genotype = genotype_all.compress(flt, axis=0)
      genotype
      # sample indices for population 1
      pop1_idx = subpops[pop1]
      # sample indices for population 2
      pop2_idx = subpops[pop2]
      num, den = allel.hudson_fst(ac1, ac2)
      snp_fst_hudson = num / den
      snp_fst_hudson
      fst_hudson, se_hudson, vb_hudson, _ = allel.average_hudson_fst(ac1, ac2, blen=10000)
      print('%.04f +/- %.04f (Hudson)' % (fst_hudson, se_hudson))
      
      ## now make some dfs!!
      regions = pop1 + "_vs_all"
      f_name1=  args.out_dir + "/"+ regions +  ".single.snp_fst_df.csv"
      f_name2= args.out_dir + "/" + regions + ".single.hudson_fst_df.csv"
      snp_fst_df = pd.DataFrame({'CHR': callset['variants/CHROM'][:], 'POS': callset['variants/POS'][:] , 'retained':flt})
      snp_fst_df_flt = snp_fst_df[(snp_fst_df['retained']==True)]
      snp_fst_df_flt['Fst'] = snp_fst_hudson
      snp_fst_df_flt.to_csv(f_name1)
      fst_hudson_df = pd.DataFrame({'Region':regions, 'Hudson_Fst': fst_hudson, 'SE': se_hudson, 'n_retained': [np.count_nonzero(flt)]})
      fst_hudson_df.to_csv(f_name2)


parser = argparse.ArgumentParser(description='run Fst analysis: genome-wide biallelic',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--zarr',help='Zarr file path',required=True)
parser.add_argument('--comparison', help='paired or one_against_all or both')
parser.add_argument('--population', help='Column name containing populations')
parser.add_argument('--out_dir', help='out_dir')
parser.add_argument('--meta', help='metadata tsv'),
parser.add_argument('--wgs_id', default="wgs_id",  help='column name with wgs_id')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)