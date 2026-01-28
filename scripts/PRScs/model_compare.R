prs0 <- fread("results/PGS/gd_Mother_parity0/PGS.txt") %>%
  mutate(PRS_z=scale(PRS))
prs1 <- fread("results/PGS/gd_Mother_parity1/PGS.txt")%>%
  mutate(PRS_z=scale(PRS))
prsall <- fread("results/PGS/gd_Mother_effectSNPnoInteraction_parityall/PGS.txt")%>%
  mutate(PRS_z=scale(PRS))


df0 <- prs0 %>%
  inner_join(cov,by=c("FID","IID"))
df1 <- prs1 %>%
  inner_join(cov,by=c("FID","IID"))
dfall <- prsall %>%
  inner_join(cov,by=c("FID","IID"))

mod_p0 <- lm(phe ~ maternal_age_at_delivery + sex_assigned_at_birth +  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=df0)
mod0 <- lm(phe ~ maternal_age_at_delivery + sex_assigned_at_birth +  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PRS_z, data=df0)
summary(mod0)$r.squared - summary(mod_p0)$r.squared

mod_p1 <- lm(phe ~ maternal_age_at_delivery + sex_assigned_at_birth +  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=df1)
mod1 <- lm(phe ~ maternal_age_at_delivery + sex_assigned_at_birth +  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PRS_z, data=df1)
summary(mod1)$r.squared - summary(mod_p1)$r.squared

mod_pall <- lm(phe ~ maternal_age_at_delivery + sex_assigned_at_birth + PARITET +  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=dfall)
mod_pall_pgs <- lm(phe ~ maternal_age_at_delivery + sex_assigned_at_birth + PARITET +  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+ PRS_z, data=dfall)
summary(mod_pall_pgs)$r.squared - summary(mod_pall)$r.squared

modall_int <- lm(phe ~ maternal_age_at_delivery + sex_assigned_at_birth + PARITET + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + PRS_z + PRS_z*PARITET, data=dfall)
summary(modall_int)$r.squared - summary(mod_pall_pgs)$r.squared
