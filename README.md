## Estimate polygenic score for gestational duration---------------------------
## @Xiaoping WU
## 2025-11-24


## Data

GWAS summary statistics: gestataional duration, European ancestry,Mother
LD reference: 1000 Genomes Project phase 3 samples
Targeted sample: gestataional duration, European ancestry, Child

## Code
/mnt/scratch/xiaoping/GxE_parity_PRS-CS/workflow/Snakefile
/mnt/scratch/xiaoping/GxE_parity_PRS-CS/workflow/scripts/*

## Results

### PRS-CS posterior SNP effect size estimates for each chromosome 
/mnt/scratch/xiaoping/GxE_parity_PRS-CS/results/PRS-CS_EUR_Mother_GA/

### PGS for targeted sample
/mnt/scratch/xiaoping/GxE_parity_PRS-CS/results/PGS_EUR_Mother_GA/PGS.txt

### PGS results check

#### variance
mod0 <- lm(phe ~ mor_age + mor_age2 + PARITET_5 + KJONN + PC1 + PC2 + PC3 + PC4 + PC5, data=df)
mod1 <- lm(phe ~ mor_age + mor_age2 + PARITET_5 + KJONN + PC1 + PC2 + PC3 + PC4 + PC5 + PRS_z, data=df)

ΔR² = 0.1605 → adding the PRS to covariates explains ~16% more variance of the phenotype.

#### quantile plot
/mnt/scratch/xiaoping/GxE_parity_PRS-CS/results/PGS_EUR_Mother_GA/quantile_plot.png