## interaction effect plot

## model: GD ~ INT + PRS + parity + prs*parity
## if parity=0, GD ~ INT + PRS, so plot a line with intercept=INT, slope=effect of PRS
## if parity=1, GD ~ INT + PRS + parity + prs*parity, so  plot a line with intercept = INT + effect of parity, slope = effect of PRS + effect of interaction

library(data.table)
library(dplyr)
library(ggplot2)


phe='results/replication/gd_Mother_parityall/phe.tsv'
cov='results/replication/gd_Mother_parityall/cov.tsv'
prs="results/PGS/gd_Mother_effectSNPnoInteraction_parityall/PGS.txt"
genome="Maternal"
traits="GD"
outfile="results/PGS/gd_Mother_effectSNPnoInteraction_parityall/interaction_plot.png"

  phe <- fread(phe)
  cov <- fread(cov)
  prs <- fread(prs)
  df <- prs %>%
    mutate(PRS_z=scale(as.numeric(PRS), scale = F)) %>%
    inner_join(phe,by=c("FID","IID"))%>%
    inner_join(cov,by=c("FID","IID"))
  
#  mod0 <- lm(phe ~ maternal_age_at_delivery + sex_assigned_at_birth +  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10,data=df)
#  res <- resid(mod0)
#  res_zscore <- (res-mean(res))/sd(res)
#  df$res_zscore <- res_zscore
#  df$GD <- mean(df$phe)+resid(mod0)
  
#  mod <- lm(GD ~ PARITET + PRS_z + PRS_z*PARITET, data=df)
  mod <- lm(phe ~ maternal_age_at_delivery + sex_assigned_at_birth +  PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10+PARITET + PRS_z + PRS_z*PARITET,data=df)
    print(summary(mod))
    # Call:
    #   lm(formula = phe ~ maternal_age_at_delivery + sex_assigned_at_birth + 
    #        PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + 
    #        PARITET + PRS_z + PRS_z * PARITET, data = df)
    # 
    # Residuals:
    #   Min      1Q  Median      3Q     Max 
    # -99.001  -5.025   1.907   8.159  25.064 
    # 
    # Coefficients:
    #   Estimate Std. Error t value Pr(>|t|)    
    # (Intercept)               2.777e+02  8.547e-01 324.919  < 2e-16 ***
    #   maternal_age_at_delivery -3.797e-02  3.097e-02  -1.226  0.22019    
    # sex_assigned_at_birth     8.033e-04  3.417e-04   2.351  0.01875 *  
    #   PC1                       1.785e+01  1.086e+01   1.643  0.10036    
    # PC2                      -9.101e+00  1.084e+01  -0.840  0.40102    
    # PC3                      -3.488e-01  1.084e+01  -0.032  0.97433    
    # PC4                      -7.235e+00  1.084e+01  -0.668  0.50433    
    # PC5                       1.197e+01  1.084e+01   1.105  0.26922    
    # PC6                      -1.922e+01  1.083e+01  -1.774  0.07613 .  
    # PC7                      -1.693e+01  1.082e+01  -1.565  0.11775    
    # PC8                       1.588e+00  1.083e+01   0.147  0.88340    
    # PC9                      -1.500e+01  1.084e+01  -1.383  0.16667    
    # PC10                      3.151e+00  1.085e+01   0.290  0.77147    
    # PARITET                   1.239e+00  2.917e-01   4.247 2.20e-05 ***
    #   PRS_z                     7.247e+00  9.859e-01   7.350 2.24e-13 ***
    #   PARITET:PRS_z            -3.752e+00  1.312e+00  -2.860  0.00426 ** 
    #   ---
    #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
    # 
    # Residual standard error: 10.79 on 6073 degrees of freedom
    # Multiple R-squared:  0.01738,	Adjusted R-squared:  0.01496 
    # F-statistic: 7.163 on 15 and 6073 DF,  p-value: 7.189e-16
    # 
  s <- coef(summary(mod))
  int0 <- s["(Intercept)","Estimate"]
  slope0 <- s["PRS_z","Estimate"]
  
  int1 <- s["(Intercept)","Estimate"] + s["PARITET","Estimate"]
  slope1 <- s["PRS_z","Estimate"] + s["PARITET:PRS_z","Estimate"]
  
  interactionpvalue = formatC(s["PARITET:PRS_z","Pr(>|t|)"], format = "e", digits = 2)
  
  plotdata = c()
  plotdata = as.data.frame(rbind(plotdata, c(int0,int1,slope0,slope1,interactionpvalue,traits,genome)))
  names(plotdata) <- c("int0","int1","slope0","slope1","interactionpvalue","traits","genome")
  
  legend_df <- data.frame(
    Parity = c("= 0", "> 0")
  )
  p <- plotdata %>% ggplot() +
    coord_cartesian(
      xlim = c(round(min(df$PRS_z),0), round(max(df$PRS_z),0)),
      ylim = c(259, 300)
    )+
    geom_hline(yintercept = 0, alpha = .5, color = "black") +
 #   geom_vline(xintercept = 0, alpha = .5, color = "black") +
    xlab(paste("Standardized polygenic genetic score", paste(genome, "genome",sep =" "), sep="\n")) +
    ylab(paste(traits)) +
    geom_abline(intercept = as.numeric(plotdata$int0), slope = as.numeric(plotdata$slope0), color="#008837", linewidth=1.5) +
    geom_abline(intercept = as.numeric(plotdata$int1), slope = as.numeric(plotdata$slope1), color="#7b3294", linewidth=1.5) +
#    geom_abline(aes(intercept = as.numeric(plotdata$int0), slope = as.numeric(plotdata$slope0), color="= 0"), linewidth=1.5) +
 #   geom_abline(aes(intercept = as.numeric(plotdata$int1), slope = as.numeric(plotdata$slope1), color="> 0"), linewidth=1.5) +
#    scale_colour_manual(name = "Parity", values=c("= 0"="#008837","> 0"="#7b3294"))+
    annotate("text",x = max(df$PRS_z, na.rm = TRUE),y = 300,hjust = 1,label = paste0("PGS × Parity p = ", interactionpvalue), size = 3) +
    theme_bw()+
    theme(
      legend.position = "right",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size=11, color="black"),
      axis.title = element_text(size=11, color="black"),
      strip.background = element_blank(),
      strip.placement = "outside",
      strip.text.x = element_text(size = 11, color="black"),
      panel.border = element_rect(color ="black", fill=NA, size=0.5),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA)
    ) +
    # Add unrelated legend
    geom_blank(
      data = legend_df,
      aes(color = Parity)
    ) +
    scale_color_manual(
      name = "Parity", 
      values = c("= 0" = "#008837", "> 0" = "#7b3294")
    ) 
    
  ggsave(filename=outfile,plot=p)

 
# ggplot(df,aes(x=PRS_z,y=phe,color=as.factor(PARITET))) +
#   # geom_point(size=0.5) +
#   geom_smooth(method=lm, se=FALSE, fullrange=TRUE,linewidth=0.5)+
#   scale_color_manual(values=c('#008837','#7b3294'))+
#   theme(legend.position="top")