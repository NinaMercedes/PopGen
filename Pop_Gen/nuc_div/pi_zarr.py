import numpy as np
import scipy
import h5py
import allel
from allel import chunked
import pandas as pd
import zarr
import itertools
import sys
import os
import argparse
import subprocess
print('scikit-allel', allel.__version__)

def nucleotide_diversity(zarrfile, metadata, out_dir, population, wgs_id):
  callset = zarr.open(zarrfile, mode='r')
  sample_names = callset['samples'][:]
  pos_index = allel.ChromPosIndex(callset['variants/CHROM'][:], callset['variants/POS'][:])
  genotype_all = allel.GenotypeChunkedArray(callset['calldata']['GT'])
  df_samples = pd.read_csv(metadata, sep="\t")
  df_samples['domain'] = df_samples[wgs_id]
  df_samples =pd.DataFrame({'domain':sample_names}).merge(df_samples, on='domain', how='left')
  geo_r = df_samples[population].unique() 
  geo_r  = [x for x in geo_r if str(x) != 'nan'] 
  #geo_r = ["Indonesia", "Eritrea", "Brazil", "Ecuador", "French_Guiana", "Guyana", "Zambia", "Burkina_Faso", "Guinea_Bissau", "ivory_coast"] 
  for pop in geo_r: #check combs vs permutations
    print(pop)
    subpops = {
    pop: df_samples[df_samples[population] == pop].index
    }
    # allele counts
    acs = genotype_all.count_alleles_subpops(subpops)
    acu = allel.AlleleCountsArray(acs[pop])
    pi = allel.mean_pairwise_difference(acu)
    pi_df = pd.DataFrame({'pi':pi})
    f_name = out_dir + "/" + pop + "_pi.csv"
    pi_df.to_csv(f_name)   


def main(args):
  nucleotide_diversity(args.zarr, args.metadata, args.out_dir, args.population, args.wgs_id)

parser = argparse.ArgumentParser(description='Pairwise diversity statistics analysis',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--metadata',help='tsv file containing metadata',required=True)
parser.add_argument('--zarr',help='zarr file',required=True)
parser.add_argument('--population', help='Column name containing populations')
parser.add_argument('--out_dir', help='out_dir')
parser.add_argument('--wgs_id', default="wgs_id",  help='column name with wgs_id')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)