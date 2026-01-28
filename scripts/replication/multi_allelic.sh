## SNP 1:154948327 has multiple allele in ALSPAC


#===========for bw, child parity0
plink2 --bfile  /mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr1  --snp 1:154948327 --make-bed --out tmp/chr1_s
awk -v OFS="\t" '{print $1,$2":"$5":"$6,$3,$4,$5,$6}' tmp/chr1_s.bim >tmp/updatesnp.txt

rm tmp/chr1_s.bim
mv tmp/updatesnp.txt tmp/chr1_s.bim

plink2 --bfile tmp/chr1_s --snps 1:154948327:G:T --make-bed --out tmp/chr1_ss

regenie \
    --step 2 \
    --bed tmp/chr1_ss \
    --covarFile results/replication/bw_zscore_Child_parity0/cov.tsv \
    --phenoFile results/replication/bw_zscore_Child_parity0/phe.tsv \
    --covarColList maternal_age_at_delivery,GA,PC{1:10} \
    --catCovarList  sex_assigned_at_birth \
    --bsize 1000 \
    --qt \
    --ignore-pred \
    --lowmem \
    --lowmem-prefix results/replication/bw_zscore_Child_parity0/chr1_s \
    --threads 5 \
    --no-condtl \
    --out results/replication/bw_zscore_Child_parity0/chr1_s

cat results/replication/bw_zscore_Child_parity0/chr1_s_phe.regenie 

#===========for bw, child parityall
plink2 --bfile  /mnt/scratch/xiaoping/gestational_duration/results/ALSPAC_PGS/geno/chr1  --snp 1:154948327 --make-bed --out tmp/chr1_s
awk -v OFS="\t" '{print $1,$2":"$5":"$6,$3,$4,$5,$6}' tmp/chr1_s.bim >tmp/updatesnp.txt

rm tmp/chr1_s.bim
mv tmp/updatesnp.txt tmp/chr1_s.bim



regenie \
    --step 2 \
    --bed tmp/chr1_ss \
    --covarFile results/replication/bw_zscore_Child_parityall/cov.tsv \
    --phenoFile results/replication/bw_zscore_Child_parityall/phe.tsv \
    --covarColList maternal_age_at_delivery,GA,PC{1:10} \
    --catCovarList  sex_assigned_at_birth \
    --bsize 1000 \
    --qt \
    --ignore-pred \
    --lowmem \
    --lowmem-prefix results/replication/bw_zscore_Child_parityall/chr1_s \
    --threads 5 \
    --no-condtl \
    --out results/replication/bw_zscore_Child_parityall/chr1_s

cat  results/replication/bw_zscore_Child_parityall/chr1_s_phe.regenie
