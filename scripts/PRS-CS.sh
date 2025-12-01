#!/bin/bash

ref_dire=$(realpath $1)
sumstat=$(realpath $2)
tag=$(echo $3)
out_prefix=$(realpath $4)
bimfile=$(realpath $5)
chr=$(echo $6)
code=$(realpath $7)


if [[ $tag == "EUR_Mother_GA" ]]; then
	sample_size=59537
elif [[ $tag == "EUR_Child_GA" ]]; then
	sample_size=47369
else
	echo "No matching condition found"
	exit 1
fi

echo "Sample size is $sample_size"


## prepare bim file
## Full path and the prefix of the bim file for the target (validation/testing) dataset. This file is used to provide a list of SNPs that are available in the target dataset.



N_THREADS=1
export MKL_NUM_THREADS=$N_THREADS
export NUMEXPR_NUM_THREADS=$N_THREADS
export OMP_NUM_THREADS=$N_THREADS


python --version
python $code --ref_dir=${ref_dire} --bim_prefix=${bimfile} --sst_file=${sumstat} --n_gwas=$sample_size --out_dir=${out_prefix} --chrom=${chr}

