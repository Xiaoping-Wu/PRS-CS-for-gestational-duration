library(optparse)
library(data.table)
library(dplyr)
library(ggplot2)

parsed_opts <- list(
  make_option("--PRS",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--phe",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--cov",
              action = "store",
              default = NA,
              type = 'character'),
  make_option("--quantile_plot",
              action = "store",
              default = NA,
              type = 'character') 
)


if (interactive()) {
  opt <- list(PRS="results/PGS_EUR_Mother_GA/PGS.txt",
              phe="data/EUR_Child_GA_phe.tsv",
              cov="data/EUR_Child_GA_cov.tsv",
              quantile_plot="results/PGS_EUR_Mother_GA/quantile_plot.png"
  )
} else {
  opt <- parse_args(OptionParser(option_list = parsed_opts))
  
}
#incremental R² and β per SD

phe <- fread(opt$phe)
cov <- fread(opt$cov) %>%
  inner_join(phe,by=c("FID","IID"))

df <- fread(opt$PRS) %>%
  inner_join(cov,by=c("FID","IID"))

df$PRS_z <- scale(df$PRS)

# baseline model with covariates
mod0 <- lm(phe ~ mor_age + mor_age2 + PARITET_5 + KJONN + PC1 + PC2 + PC3 + PC4 + PC5, data=df)
mod1 <- lm(phe ~ mor_age + mor_age2 + PARITET_5 + KJONN + PC1 + PC2 + PC3 + PC4 + PC5 + PRS_z, data=df)

summary(mod1)$coef["PRS_z",]         
r2_0 <- summary(mod0)$r.squared
r2_1 <- summary(mod1)$r.squared
delta_r2 <- r2_1 - r2_0
cat("Incremental R2:", delta_r2, "\n")
#	ΔR² = 0.1605 → adding the PRS to your covariates explains ~16% more variance of the phenotype.


##=====================PRS quantile
get_quantile <- function(x, num.quant, quant.ref){
  quant <- as.numeric(cut(x,
                          breaks = unique(quantile(
                            x, probs = seq(0, 1, 1 / num.quant)
                          )),
                          include.lowest = T))
  if(is.null(quant.ref)){
    quant.ref <- ceiling(num.quant / 2)       #define reference factor
  }
  quant <- factor(quant, levels = c(quant.ref, seq(min(quant), max(quant), 1)[-quant.ref]))
  return(quant)
}




lr <- lm(phe ~ mor_age + mor_age2 + KJONN + PC1 + PC2 + PC3 + PC4 + PC5, data=df)
lr_res <- resid(lr)
lr_res_zscore <- (lr_res-mean(lr_res))/sd(lr_res)
df$lr_res_zscore <- lr_res_zscore

df$quantile_group <- get_quantile(df$PRS_z,num.quant=10, quant.ref=1)
q975 <- qnorm(0.975)

lr_case <- summary(lm(lr_res_zscore ~ -1 + quantile_group, data=subset(df, PARITET_5 == 1)))
lr_ctl <- summary(lm(lr_res_zscore ~ -1 + quantile_group, data=subset(df, PARITET_5 == 0)))


lr_case_m <- lr_case$coefficients[,1]
lr_case_sd <- lr_case$coefficients[,2]
lr_case_counts <- table(df[df[,PARITET_5]==1,"quantile_group"])
clr_case <- cbind(c(1:10),lr_case_m,lr_case_m+q975*lr_case_sd,lr_case_m-q975*lr_case_sd,lr_case_counts,"Parity=1")
colnames(clr_case) <- c("Quantile","Coef","CI.U","CI.L","N","sample")

lr_ctl_m <- lr_ctl$coefficients[,1]
lr_ctl_sd <- lr_ctl$coefficients[,2]
lr_ctl_counts <- table(df[df[,PARITET_5]==0,"quantile_group"])
clr_ctl <- cbind(c(1:10),lr_ctl_m,lr_ctl_m+q975*lr_ctl_sd,lr_ctl_m-q975*lr_ctl_sd,lr_ctl_counts,"Parity=0")
colnames(clr_ctl) <- c("Quantile","Coef","CI.U","CI.L","N","sample")

effect <- rbind(clr_case,clr_ctl)
effect <- as.data.frame(effect)
effect$Quantile <- as.numeric(effect$Quantile)
effect$Coef <- as.numeric(effect$Coef)
effect$CI.U <- as.numeric(effect$CI.U)
effect$CI.L <- as.numeric(effect$CI.L)

myplot <- ggplot(effect,aes(x=Quantile,y=Coef,ymin=CI.L,ymax=CI.U,color=sample))+
  geom_point(position=position_dodge(0.5))+
  scale_color_manual(values=c("#F8766D","#A3A500"))+
  geom_errorbar(width=0.1,position=position_dodge(0.5)) +
  xlab("PGS Quantiles")+
  ylab("Effect (95% CI)")+
  scale_x_continuous(breaks=seq(1,10,1)) +
  theme_bw() +
  theme(legend.title = element_blank())
myplot
ggsave(opt$quantile_plot,myplot,height=4,width=4)


