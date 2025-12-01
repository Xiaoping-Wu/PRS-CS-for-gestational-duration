import pandas as pd
import scipy
import numpy as np
import h5py
import os

if ((snakemake.wildcards.pheno == "BW") & (snakemake.wildcards.sample == "fetal")):
	sample_size= 298142
if ((snakemake.wildcards.pheno == "BW") & (snakemake.wildcards.sample == "maternal")):
	sample_size= 210267
if ((snakemake.wildcards.pheno == "PW") & (snakemake.wildcards.sample == "fetal")):
	sample_size= 65405
if ((snakemake.wildcards.pheno == "PW") & (snakemake.wildcards.sample == "maternal")):
	sample_size= 61228
if ((snakemake.wildcards.pheno == "GD") & (snakemake.wildcards.sample == "fetal")):
	sample_size= 84689
if ((snakemake.wildcards.pheno == "GD") & (snakemake.wildcards.sample == "maternal")):
	sample_size= 151987

print(sample_size)

os.system("N_THREADS=1; export MKL_NUM_THREADS=$N_THREADS; export NUMEXPR_NUM_THREADS=$N_THREADS; export OMP_NUM_THREADS=$N_THREADS; python --version; python PRScs/PRScs.py --ref_dir={snakemake.params[1]} --bim_prefix={snakemake.input[2]} --sst_file={snakemake.input[0]} --n_gwas={sample_size} --out_dir={snakemake.params[0]} --chrom=22")

