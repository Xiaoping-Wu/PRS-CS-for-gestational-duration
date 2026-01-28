m_p0 <- fread("results/replication/gd_Mother_parity0/gwas_replication.tsv") %>%
  mutate(gwas="Mother_parity0")
m_p1 <- fread("results/replication/gd_Mother_parity1/gwas_replication.tsv") %>%
  mutate(gwas="Mother_parity1")
m_p <- fread("results/replication/gd_Mother_parityall/gwas_replication.tsv") %>%
  mutate(gwas="Mother_parityall")

df <- bind_rows(m_p0,m_p1,m_p)

write.table(df,"tmp/df.tsv",col.names=T,row.names=F,quote=F,sep="\t")

c1 <- cor(m_p0[same_sign==T,]$BETA.d,m_p0[same_sign==T,]$BETA.r.harmonized,use="complete.obs")
c2 <- cor(m_p1[same_sign==T,]$BETA.d,m_p1[same_sign==T,]$BETA.r.harmonized,use="complete.obs")
c3 <- cor(m_p[same_sign==T,]$BETA.d,m_p[same_sign==T,]$BETA.r.harmonized,use="complete.obs")
cor(df[same_sign==T,]$BETA.d,df[same_sign==T,]$BETA.r.harmonized,use="complete.obs") #0.8568011
cat(c1,c2,c3)
#0.8180256 -1 0.894897
##==============
m_p0 <- fread("results/replication/gd_Child_parity0/gwas_replication.tsv") %>%
  mutate(gwas="Child_parity0")
m_p <- fread("results/replication/gd_Child_parityall/gwas_replication.tsv") %>%
  mutate(gwas="Child_parityall")

df <- bind_rows(m_p0,m_p)

write.table(df,"tmp/df.tsv",col.names=T,row.names=F,quote=F,sep="\t")
c1 <- cor(m_p0[same_sign==T,]$BETA.d,m_p0[same_sign==T,]$BETA.r.harmonized,use="complete.obs")
c3 <- cor(m_p[same_sign==T,]$BETA.d,m_p[same_sign==T,]$BETA.r.harmonized,use="complete.obs")
cor(df[same_sign==T,]$BETA.d,df[same_sign==T,]$BETA.r.harmonized,use="complete.obs") #0.8732692
cat(c1,c3) 
#NA 0.9999532
##================
m_p0 <- fread("results/replication/bw_zscore_Mother_parity0/gwas_replication.tsv") %>%
  mutate(gwas="Mother_parity0")
m_p1 <- fread("results/replication/bw_zscore_Mother_parity1/gwas_replication.tsv") %>%
  mutate(gwas="Mother_parity1")
m_p <- fread("results/replication/bw_zscore_Mother_parityall/gwas_replication.tsv") %>%
  mutate(gwas="Mother_parityall")

df <- bind_rows(m_p0,m_p1,m_p)

write.table(df,"tmp/df.tsv",col.names=T,row.names=F,quote=F,sep="\t")

c1 <- cor(m_p0[same_sign==T,]$BETA.d,m_p0[same_sign==T,]$BETA.r.harmonized,use="complete.obs")
c2 <- cor(m_p1[same_sign==T,]$BETA.d,m_p1[same_sign==T,]$BETA.r.harmonized,use="complete.obs")
c3 <- cor(m_p[same_sign==T,]$BETA.d,m_p[same_sign==T,]$BETA.r.harmonized,use="complete.obs")
cor(df[same_sign==T,]$BETA.d,df[same_sign==T,]$BETA.r.harmonized,use="complete.obs") ##0.9029054
cat(c1,c2,c3) 
#0.8848123 0.9241725 0.905004

##==============
m_p0 <- fread("results/replication/bw_zscore_Child_parity0/gwas_replication.tsv") %>%
  mutate(gwas="Child_parity0")
m_p1 <- fread("results/replication/bw_zscore_Child_parity1/gwas_replication.tsv") %>%
  mutate(gwas="Child_parity1")
m_p <- fread("results/replication/bw_zscore_Child_parityall/gwas_replication.tsv") %>%
  mutate(gwas="Child_parityall")

df <- bind_rows(m_p0,m_p1,m_p)

write.table(df,"tmp/df.tsv",col.names=T,row.names=F,quote=F,sep="\t")
c1 <- cor(m_p0[same_sign==T,]$BETA.d,m_p0[same_sign==T,]$BETA.r.harmonized,use="complete.obs")
c2 <- cor(m_p1[same_sign==T,]$BETA.d,m_p1[same_sign==T,]$BETA.r.harmonized,use="complete.obs")
c3 <- cor(m_p[same_sign==T,]$BETA.d,m_p[same_sign==T,]$BETA.r.harmonized,use="complete.obs")
cor(df[same_sign==T,]$BETA.d,df[same_sign==T,]$BETA.r.harmonized,use="complete.obs") #0.8885855
cat(c1,c2,c3)
#0.8248598 0.9235442 0.9207352
