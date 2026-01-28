## Estimate polygenic score for gestational duration
## @Xiaoping Wu


## 2026-01-28 Polygenic score for Maternal Gestational duration(GD)
/mnt/scratch/xiaoping/GxE_parity_PRS-CS/workflow/rules/replication.smk

GWAS summary statistics: Moba, Maternal, GD, first pregnancy (parity0)/later pregnancy (parity1)/all pregnancy (parityall)

target sample: ALSPAC cohorts

### Results
1. interaction effect p-value
2. Polygenic score effect on first pregnancy VS later pregnancy
  model: GD ~ INT + PRS + parity + prs*parity  
  if parity=0, GD ~ INT + PRS, so plot a line with intercept=INT, slope=effect of PRS  
  if parity=1, GD ~ INT + PRS + parity + prs*parity, so plot a line with intercept = INT + effect of parity, slope = effect of PRS + effect of interaction

  
  mod <- lm(phe ~ othercov+PARITET + PRS_z + PRS_z*PARITET,data=df)  
  s <- coef(summary(mod))  
  int0 <- s["(Intercept)","Estimate"]  
 slope0 <- s["PRS_z","Estimate"]  
  
  int1 <- s["(Intercept)","Estimate"] + s["PARITET","Estimate"]  
  slope1 <- s["PRS_z","Estimate"] + s["PARITET:PRS_z","Estimate"]  
 
  interactionpvalue = formatC(s["PARITET:PRS_z","Pr(>|t|)"], format = "e", digits = 2)  


## 2026-01-28 Replication of parity-stratified GWASs in ALSPAC

/mnt/scratch/xiaoping/GxE_parity_PRS-CS/workflow/rules/PGS.smk



## 2025-11-24 test for PRScs


## Data

GWAS summary statistics: gestataional duration, European ancestry,Mother

LD reference: 1000 Genomes Project phase 3 samples

Targeted sample: gestataional duration, European ancestry, Child

## Code
/mnt/scratch/xiaoping/GxE_parity_PRS-CS/workflow/rules/PGS_test.smk


### PGS results check

#### variance
mod0 <- lm(phe ~ mor_age + mor_age2 + PARITET_5 + KJONN + PC1 + PC2 + PC3 + PC4 + PC5, data=df)

mod1 <- lm(phe ~ mor_age + mor_age2 + PARITET_5 + KJONN + PC1 + PC2 + PC3 + PC4 + PC5 + PRS_z, data=df)

ΔR² = 0.1605 → adding the PRS to covariates explains ~16% more variance of the phenotype.

#### quantile plot
/mnt/scratch/xiaoping/GxE_parity_PRS-CS/results/PGS_EUR_Mother_GA/quantile_plot.png
