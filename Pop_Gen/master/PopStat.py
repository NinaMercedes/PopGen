#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
from tqdm import tqdm
import time
import zarr
import allel
import os.path
import pandas as pd


def main_zarr(args):
  today = args.date 
  out_dir = args.prefix + "_ZARR_" +  today
  subprocess.run("mkdir " + out_dir, shell=True)
  out_zarr = out_dir + "/" + args.prefix + ".zarr"
  allel.vcf_to_zarr(args.vcf, out_zarr) #'/mnt/storage12/nbillows/Projects/Malaria/20k_finale_QC_pass.bi.GT.miss0.4.vqslod.filt.snps.vcf.gz'
  out_zarr_coding = out_dir + "/" + args.prefix + "_coding.zarr"
  allel.vcf_to_zarr(args.vcf_coding, out_zarr_coding)
  
def main_pca(args):
  today = args.date 
  out_dir = args.prefix + "_PCA_" + today
  subprocess.run("mkdir " + out_dir, shell=True)
  if "run" in args.make_zarr:
    out_dir_zarr = args.prefix + "_ZARR_" + today
    zarr = "./" + out_dir_zarr + "/" + args.prefix + ".zarr"
  else:
    zarr = args.zarr
  cmd = "python " + args.path + "pca/fastpca.py " + " --zarr " + zarr + " --out_dir " + out_dir + " --prefix " + args.prefix
  subprocess.run(cmd, shell=True)
  
def main_fst(args):
  today = args.date 
  out_dir = args.prefix + "_FST_" + today
  subprocess.run("mkdir " + out_dir, shell=True)
  if "run" in args.make_zarr:
    out_dir_zarr = args.prefix + "_ZARR_" + today
    zarr = out_dir_zarr + "/" + args.prefix + ".zarr"
  else:
    zarr = args.zarr
  cmd = "python " + args.path + "/fst/run_fst.py " + " --zarr " + zarr + " --out_dir " + out_dir + " --population " + args.population + " --meta " + args.metadata + " --wgs_id " + args.wgs_id + " --comparison " + args.comparison
  subprocess.run(cmd, shell=True)

def main_pi(args):
  today = args.date 
  out_dir = args.prefix + "_PI_" + today
  subprocess.run("mkdir " + out_dir, shell=True)
  if "run" in args.make_zarr:
    out_dir_zarr = args.prefix + "_ZARR_" + today
    args.zarr_coding = out_dir_zarr + "/" + args.prefix + "_coding.zarr"
  cmd = "python " + args.path + "/nuc_div/pi_zarr.py " + " --zarr " + args.zarr_coding + " --out_dir " + out_dir + " --population " + args.population + " --metadata " + args.metadata + " --wgs_id " + args.wgs_id
  subprocess.run(cmd, shell=True)
  
def main(args):
  if "run" in args.make_zarr:
    main_zarr(args)
  if "run" in args.do_pca:  
    main_pca(args)
  if "run" in args.do_fst:  
    main_fst(args)
  if "run" in args.do_pi:  
    main_pi(args)
  
  
  
parser = argparse.ArgumentParser(description='run Pop stat analysis',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--path',help='path to PopGen Pipeline scripts',required=True)
parser.add_argument('--vcf',help='vcf file path',required=False)
parser.add_argument('--vcf_coding',help='vcf file path',required=False)
parser.add_argument('--zarr',help='Zarr file path',required=False)
parser.add_argument('--zarr_coding',help='Zarr file path fro coding regions',required=False)
parser.add_argument('--comparison', help='paired or one_against_all or both')
parser.add_argument('--population', help='Column name containing populations')
parser.add_argument('--make_zarr', default='False', help='run')
parser.add_argument('--metadata', help='metadata tsv'),
parser.add_argument('--wgs_id', default="wgs_id",  help='column name with wgs_id')
parser.add_argument('--do_pca', default='False', help='run')
parser.add_argument('--do_fst', default='False', help='run')
parser.add_argument('--do_pi', default='False', help='run')
parser.add_argument('--date',  help='2024_11_09')
parser.add_argument('--prefix',  help='check')
parser.set_defaults(func=main)

args = parser.parse_args()
args.func(args)