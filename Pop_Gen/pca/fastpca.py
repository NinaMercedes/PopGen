import random
import time
import numpy as np
import h5py
import matplotlib.pyplot as plt
import seaborn as sns
#import bcolz
import pandas
import allel
import zarr 
import argparse
import pandas as pd

print('scikit-allel', allel.__version__)
random.seed(35)
sns.set_style('white')
sns.set_style('ticks')

def ld_prune(gn, size, step, threshold=.1, n_iter=1):
    for i in range(n_iter):
        loc_unlinked = allel.locate_unlinked(gn, size=size, step=step, threshold=threshold)
        n = np.count_nonzero(loc_unlinked)
        n_remove = gn.shape[0] - n
        print('iteration', i+1, 'retaining', n, 'removing', n_remove, 'variants')
        gn = gn.compress(loc_unlinked, axis=0)
    return gn

def main(args):
  callset = zarr.open(args.zarr, mode='r')
  g = allel.GenotypeChunkedArray(callset['calldata']['GT'])
  ac = g.count_alleles()[:]
  ac
  np.count_nonzero(ac.max_allele() > 1)
  np.count_nonzero((ac.max_allele() == 1) & ac.is_singleton(1))
  flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
  gf = g.compress(flt, axis=0)
  gf
  gn = gf.to_n_alt()
  gn
  gnu = ld_prune(gn, size=50, step=10, threshold=.1, n_iter=1)
  coords1, model1 = allel.pca(gnu, n_components=10, scaler='patterson')
  
  samples = callset['samples'][:]
  pca_df =  pd.DataFrame(np.vstack(coords1))
  pca_df['samples'] = samples
  pca_df2 =  pd.DataFrame(np.vstack(model1.explained_variance_ratio_[:]))
  num_snps = pd.DataFrame({'num_snps': gnu.shape[0]}, index=[0])
  
  pca_df.to_csv(args.out_dir +"/" + args.prefix + "_pca_coordinates.csv")
  pca_df2.to_csv(args.out_dir +"/" + args.prefix +  "_pca_var_explained.csv")
  pca_df2.to_csv(args.out_dir +"/" + args.prefix +  "_pca_num_snps.csv")
  
  
  
parser = argparse.ArgumentParser(description='run Fst PCA: genome-wide biallelic LD removed',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--zarr',help='Zarr file path',required=True)
parser.add_argument('--prefix', default="pf",help = "prefix")
parser.add_argument('--out_dir', help='out_dir')

parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)

  