## effect of interaction effect on GD

phe='results/replication/gd_Mother_parity0/phe.tsv'
cov='results/replication/gd_Mother_parity0/cov.tsv'
prs="results/PGS/gd_Mother_parity0/PGS.txt"

phe <- fread(phe)
cov <- fread(cov)
prs <- fread(prs)
df <- prs %>%
  mutate(PRS_z=scale(as.numeric(PRS), scale = F)) %>%
  inner_join(phe,by=c("FID","IID"))%>%
  inner_join(cov,by=c("FID","IID"))
mod <- lm(phe ~ maternal_age_at_delivery + sex_assigned_at_birth +  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+PRS_z,data=df)
con_mod <- confint(mod)

out <- as.data.frame(matrix(2,2,9))
names(out) <- c("Parity","Estimate",  "SE","tvalue","Pval","CIL","CIU","rsq","rsq_adj")
out[1,1] <- "First pregancy"
out[1,2:5] <- coef(summary(mod))["PRS_z",]
out[1,6:7] <- con_mod["PRS_z",]
out[1,8] <- summary(mod)$r.squared
out[1,9] <- summary(mod)$adj.r.squared

##
phe='results/replication/gd_Mother_parity1/phe.tsv'
cov='results/replication/gd_Mother_parity1/cov.tsv'
prs="results/PGS/gd_Mother_parity1/PGS.txt"

phe <- fread(phe)
cov <- fread(cov)
prs <- fread(prs)
df <- prs %>%
  mutate(PRS_z=scale(as.numeric(PRS), scale = F)) %>%
  inner_join(phe,by=c("FID","IID"))%>%
  inner_join(cov,by=c("FID","IID"))
mod <- lm(phe ~ maternal_age_at_delivery + sex_assigned_at_birth +  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+PRS_z,data=df)
con_mod <- confint(mod)

out[2,1] <- "Later pregancy"
out[2,2:5] <- coef(summary(mod))["PRS_z",]
out[2,6:7] <- con_mod["PRS_z",]
out[2,8] <- summary(mod)$r.squared
out[2,9] <- summary(mod)$adj.r.squared


print(out)
# Parity Estimate        SE   tvalue         Pval      CIL      CIU         rsq      rsq_adj
# 1 First pregancy 4.246204 0.7806338 5.439432 5.823109e-08 2.715505 5.776903 0.021575286 0.0168747975
# 2 Later pregancy 1.768221 0.6558957 2.695887 7.054802e-03 0.482230 3.054212 0.004558525 0.0007434481
