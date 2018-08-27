#PAWS HGHEDLV2 hghedlv2-3cat Regressions

#Load Data
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(lme4)

setwd("/big_data/redgar/PAWS/PAWS_final")
load("combat_PAWS_Beta_norep.RData")
load("PAWS_meta_sentrix_genetic_clusters.RData")
meta<-meta[which(meta$Meth_ID%in%colnames(PAWS_Beta)),]
meta<-meta[match(colnames(PAWS_Beta), meta$Meth_ID),]

#M Values 
Mval<-function(beta) log2(beta/(1-beta))
PAWS_Mval = apply(PAWS_Beta, 1, Mval) # need mvalues for combat
PAWS_Mval = as.data.frame(PAWS_Mval)
PAWS_Mval = t(PAWS_Mval)


## Variable Distribution
hist_plot<-meta[,c(1,20,21)]
hist_plot<-melt(hist_plot, id="Meth_ID")
ggplot(hist_plot, aes(value, fill=variable))+geom_histogram(position="identity", color="black", alpha=0.5, binwidth = 0.25)+theme_bw()+
  scale_fill_manual(values=c("lightblue","blue","darkblue"))+facet_wrap(~variable, ncol=1)



# Liner Mixed Effects Model for HGHEDLV2 MVAL
# LME P value calculation
meta$TWIN_Pair<-as.factor(meta$TWIN_Pair)
Likelihood_Ratio_Test<-lapply(1:nrow(PAWS_Mval), function(CpG){
  metaex<-meta
  metaex$Beta<-PAWS_Mval[CpG,]
  mod_edu<-lmer(Beta ~ HGHEDLV2 + as.factor(minor_child)+as.factor(Genetic_cluster)+(1|TWIN_Pair) , data=metaex)
  mod_covar<-lmer(Beta ~ as.factor(minor_child)+as.factor(Genetic_cluster)+(1|TWIN_Pair) , data=metaex)
  Likelihood_Ratio_Test<-anova(mod_edu, mod_covar)
  Likelihood_Ratio_Test})
save(Likelihood_Ratio_Test, file="Likelihood_Ratio_Test_lmer_HGHEDLV2_PAWS_Mval_genetic_cluster_and_twin_fixed_allFactors.RData")

Likelihood_Ratio_Test_Mval<-sapply(1:length(Likelihood_Ratio_Test), function(x) Likelihood_Ratio_Test[[x]][["Pr(>Chisq)"]][2])

pvalue_dist_hghedlv2<-data.frame(CpG=rownames(PAWS_Mval), Nominal_P=Likelihood_Ratio_Test_Mval)
ggplot(pvalue_dist_hghedlv2, aes(Nominal_P))+geom_histogram(fill="grey90", color="black")+theme_bw()+xlab("Nominal P Value")+
  ylim(0,50000)

length(which(Likelihood_Ratio_Test_Mval<0.001))# 4172 Mvalue (3961 when used beta)
#multiple test correction
Multi_test_corr_relaxed<-p.adjust(Likelihood_Ratio_Test_Mval, method = "fdr", n = length(Likelihood_Ratio_Test_Mval))
stat_hits<-as.data.frame(PAWS_Beta)[which(Multi_test_corr_relaxed<=0.1),] # 830 at fdr 0.05; 7046 at fdr 0.1; 16224 at fdr 0.15



# Delta beta cutoff (Mvalue cutoff doesnt really translate)
delbeta<-sapply(1:nrow(PAWS_Beta), function(x) {
  z<-lm(unlist(PAWS_Beta[x,]) ~ meta$HGHEDLV2 + meta$minor_child+meta$Genetic_cluster)
  intercept=z$coefficients[1]
  slope=z$coefficients[2]
  as.numeric((slope*6+intercept)-intercept)})

bio_hits<-PAWS_Beta[which(abs(delbeta)>=0.1),] # 1450 CpGs
save(Multi_test_corr_relaxed, stat_hits, delbeta, bio_hits, file="HGHEDLV2_PAWS_hits_Mval.RData")# thresholded on Mval but data in stahits and biohits if betas for plotting


load("HGHEDLV2_PAWS_hits_Mval.RData")
sta_bio_hits<-stat_hits[which(rownames(stat_hits)%in%rownames(bio_hits)),]#13 at fdr 0.05; 70 at fdr 0.1; 110 at fdr 0.15



## VOLCANO
volcano<-data.frame(Adjusted_Pvalue=Multi_test_corr_relaxed, Delta_Beta=delbeta)

dB<-0.1 #delta beta cutoff
Pv<-0.1 #Pvalue cutoff

## positive delta beta more methylated with higher family education
coloredu<-sapply(1:nrow(volcano), function(x) if(volcano$Adjusted_Pvalue[x]<=Pv){
  if(abs(volcano$Delta_Beta[x])>dB){
    if(volcano$Delta_Beta[x]>dB){"Hypermethylated\n(with Potential Biological Impact)"}else{"Hypomethylated\n (with Potential Biological Impact)"}
  }else{if(volcano$Delta_Beta[x]>0){"Hypermethylated"}else{"Hypomethylated"}}}else{"Not Significantly Different"})

volcano$Interesting_CpG3<-coloredu
volcano<-volcano[which(volcano$Adjusted_Pvalue<1),]

#omg
library(scales)
ggplot(volcano, aes(Delta_Beta, -log10(Adjusted_Pvalue), color=Interesting_CpG3))+
  geom_point(shape=19, size=1)+theme_bw()+scale_y_log10()+
  scale_color_manual(values=c(muted("red", l=80, c=30),"red",muted("blue", l=70, c=40),"blue","grey"),name="Methylation Change \nWith Highest Household Education")+
  geom_vline(xintercept=c(-dB,dB), color="grey60")+geom_hline(yintercept=-log10(Pv), color="grey60")+
  ylab("Multiple Test Corrected P Value (-log10)")+xlab("Delta Beta")+xlim(-1, 1)+
  theme(axis.text = element_text(size =14, color="black"))+ guides(color = guide_legend(override.aes = list(size = 4)))



