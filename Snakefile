## Estimate polygenic score for gestational duration---------------------------
## @Xiaoping WU
## 2025-11-24


# --------------------------- Dictionaries ---------------------------
CHR=list(range(1, 23))

COHORTS=["Mother"]
TRAITS=["GA"]
ORIGIN=["EUR"]



##--------------------------- Target rules ---------------------------
 
rule run_PGS:
	input:
		expand('results/PGS_{origin}_{cohorts}_{traits}/quantile_plot.png',cohorts=COHORTS,traits=TRAITS,origin=ORIGIN),
		expand('results/PGS_{origin}_{cohorts}_{traits}/PGS.txt',cohorts=COHORTS,traits=TRAITS,origin=ORIGIN),
		expand('results/PGS_{origin}_{cohorts}_{traits}/chr{ichr}-PGS.sscore',ichr=CHR,cohorts=COHORTS,traits=TRAITS,origin=ORIGIN)
	   
rule run_PRS_CS:
  input:
    expand('results/PRS-CS_{origin}_{cohorts}_{traits}/beta_pst_eff_a1_b0.5_phiauto_chr{ichr}.txt',
    ichr=CHR,cohorts=COHORTS,traits=TRAITS,origin=ORIGIN)
      
     
 
##--------------------------- analysis rules -------------------------

rule format_sumstats:
	'Format summary statistics according to the PRS-CS.'
	input:
		'data/gwas_{origin}_{cohorts}_{traits}.tsv.gz'
	output:
		'results/prepare/gwas_{origin}_{cohorts}_{traits}.txt'
	log:
	  'log/format_sumstats_{origin}_{cohorts}_{traits}.log'
	conda:
		'./envs/Rbase.yml'
	shell:
		'''
		chmod +x ./workflow/scripts/format_sumstats.R
		Rscript --vanilla ./workflow/scripts/format_sumstats.R \
      --gwas {input[0]}  \
      --gwas_format {output[0]} \
      >& {log[0]}
		'''

rule download_LD_reference:
	'Download reference panel.'
	output:
		'data/LD-reference/ldblk_1kg_eur.tar.gz'
	params:
	   'data/LD-reference/'
	shell:
		'''
		wget https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz?dl=0 -O {output[0]}
		cd {params[0]}
		tar -zxvf {output[0]}
		'''

rule bim:
  'prepare bim file'
  input:
    '/mnt/scratch/HDGB-MoBaGenetics/2025.09.25/delivery/no_phase_moba_genotypes_2025.09.25_common.pvar',
  output:
    'results/prepare/no_phase_moba_genotypes_2025.09.25_common.bim'
  shell:
    '''
    grep -v "^#" {input[0]} |awk 'BEGIN{{OFS="\t"}} ($1!="PAR1" && $1!="PAR2" && $1!="X") {{print $1, $3,0,$2,$4,$5}}' >{output[0]}
    '''

rule PRS_CS:
  'Run PRS CS'
  input:
    sumstats='results/prepare/gwas_{origin}_{cohorts}_{traits}.txt',
    ls_ref='data/LD-reference/ldblk_1kg_eur.tar.gz',
    bim='results/prepare/no_phase_moba_genotypes_2025.09.25_common.bim',
    code="workflow/envs/PRScs/PRScs.py"
  output:
    'results/PRS-CS_{origin}_{cohorts}_{traits}/beta_pst_eff_a1_b0.5_phiauto_chr{ichr}.txt'
  params:
    ld_ref='data/LD-reference/ldblk_1kg_eur',
    out_prefix='results/PRS-CS_{origin}_{cohorts}_{traits}/beta',
    ichr='{ichr}',
    tag='{origin}_{cohorts}_{traits}',
    genotype='results/prepare/no_phase_moba_genotypes_2025.09.25_common'
  log:
    'log/PRScs_chr{ichr}_{origin}_{cohorts}_{traits}.log'
  conda:
    'envs/PRS-CS.yml'
  shell:
    """
    chmod +x ./workflow/scripts/PRS-CS.sh
    ./workflow/scripts/PRS-CS.sh {params.ld_ref}  {input.sumstats}  {params.tag}  {params.out_prefix} {params.genotype} {params.ichr} {input.code} &> {log} 
    """		

rule target_sample:
  'Use child data as target'
  input:
    multiext("/mnt/scratch/HDGB-MoBaGenetics/2025.09.25/delivery/no_phase_moba_genotypes_2025.09.25_common", ".pgen", ".psam", ".pvar"),
    'data/EUR_Child_GA_phe.tsv',
    'data/LD-reference/ldblk_1kg_eur/snpinfo_1kg_hm3'
  output:
    'results/prepare/EUR_Child_GA.pgen',
    'results/prepare/EUR_Child_GA.psam',
    'results/prepare/EUR_Child_GA.pvar'
  params:
    '/mnt/scratch/HDGB-MoBaGenetics/2025.09.25/delivery/no_phase_moba_genotypes_2025.09.25_common',
    'results/prepare/EUR_Child_GA'
  threads: 20
  resources:
    mem_mb=20000
  conda:
    'envs/plink2.yml'
  shell:
	  '''
		jobid=$$
    mkdir -p "tmp/job{jobid}"
    awk 'NR>1 ' {input[3]} >tmp/job{jobid}/keepid.txt
    awk 'NR>1 {{print $2}}' {input[4]}  >tmp/job{jobid}/snp.txt
		plink2 --pfile {params[0]} --keep tmp/job{jobid}/keepid.txt --chr 1-22 --extract tmp/job{jobid}/snp.txt --make-pgen --threads {threads} --memory {resources.mem_mb} --out {params[1]}
		rm -r tmp/job{jobid}
		'''

rule PGS:
	'Create PGS for all individuals'
	input:
		multiext("results/prepare/EUR_Child_GA", ".pgen", ".psam", ".pvar"),
		'results/PRS-CS_{origin}_{cohorts}_{traits}/beta_pst_eff_a1_b0.5_phiauto_chr{ichr}.txt',
	output:
		'results/PGS_{origin}_{cohorts}_{traits}/chr{ichr}-PGS.sscore'
	params:
		'results/prepare/EUR_Child_GA',
		'results/PGS_{origin}_{cohorts}_{traits}/chr{ichr}-PGS'
	threads:1
	resources: mem_mb=4000
	conda:
		'envs/plink2.yml'
	shell:
		'''
		 plink2 --pfile {params[0]} --score {input[3]} 2 4 6  --threads {threads} --memory {resources.mem_mb} --out {params[1]}
		'''

rule concat_PGS:
	'Concat per-CHR PGS into a single PGS.'
	input:
		expand('results/PGS_{{origin}}_{{cohorts}}_{{traits}}/chr{ichr}-PGS.sscore', ichr=CHR)
	output:
		'results/PGS_{origin}_{cohorts}_{traits}/PGS.txt'
	params:
	  'results/PGS_{origin}_{cohorts}_{traits}'
	log:
	  'log/concat_PGS_{origin}_{cohorts}_{traits}.log'
	conda:
		'envs/Rbase.yml'
	shell:
	  '''
	  chmod +x ./workflow/scripts/PRS_sum.R
		Rscript --vanilla ./workflow/scripts/PRS_sum.R \
      --PRS_chr {params[0]}  \
      --PRS {output[0]} \
      > {log} 2>&1
	  '''

rule PGS_check:
	'PGS quantile plot'
	input:
		'results/PGS_{origin}_{cohorts}_{traits}/PGS.txt',
		'data/EUR_Child_GA_phe.tsv',
		'data/EUR_Child_GA_cov.tsv'
	output:
		'results/PGS_{origin}_{cohorts}_{traits}/quantile_plot.png'
	log:
	  'log/quantile_plot_{origin}_{cohorts}_{traits}.log'
	conda:
		'envs/Rbase.yml'
	shell:
	  '''
	  chmod +x ./workflow/scripts/PRS_accuracy.R
		Rscript --vanilla ./workflow/scripts/PRS_accuracy.R \
      --PRS {input[0]}  \
      --phe {input[1]}  \
      --cov {input[2]}  \
      --quantile_plot {output[0]} \
      > {log} 2>&1
	  '''
