library(optparse)
library(data.table)
library(dplyr)

parsed_opts <- list(
  make_option("--gwas",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--gwas_format",
              action = "store",
              default = NA,
              type = 'character')
  
)
if (interactive()) {
  opt <- list(gwas="data/gwas_EUR_Mother_GA.tsv.gz",
              gwas_format="results/sumstats/gwas_EUR_Mother_GA.txt"
  )
} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
  
}

d= fread(opt$gwas) %>%
  select(ID,ALLELE1,ALLELE0,BETA,SE)
d$ALLELE1= toupper(d$ALLELE1)
d$ALLELE0= toupper(d$ALLELE0)
names(d)
names(d)= c('SNP', 'A1', 'A2', 'BETA', 'SE')

print("any missing value")
print(colSums(is.na(d)))

print("SNP without rsid")
print(sum(grepl(':', d$SNP)))

print("BETA")
print(summary(d$BETA))

print("SE")
print(summary(d$SE))


d= na.omit(d)
gwas <- d %>%
  filter(!grepl(":", SNP), !is.na(SNP)) %>%
  filter(SE>0)




fwrite(gwas, opt$gwas_format, sep= '\t')

