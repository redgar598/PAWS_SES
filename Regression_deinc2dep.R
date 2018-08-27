library(reshape)
library(ggplot2)
library(RColorBrewer)
library(lme4)

# Linear regression deinc2dep (did total household income as well but not shown here as did not make the paper)
setwd("/big_data/redgar/PAWS/PAWS_final")
load("combat_PAWS_Beta_norep.RData")
load("PAWS_meta_sentrix_genetic_clusters.RData")

meta<-meta[which(meta$Meth_ID%in%colnames(PAWS_Beta)),]
meta<-meta[match(colnames(PAWS_Beta), meta$Meth_ID),]

Mval<-function(beta) log2(beta/(1-beta))
PAWS_Mval = apply(PAWS_Beta, 1, Mval) # need mvalues for combat
PAWS_Mval = as.data.frame(PAWS_Mval)
PAWS_Mval = t(PAWS_Mval)

## Variable Distribution
hist_plot<-meta[,c(1,22,24,32)]
hist_plot<-melt(hist_plot, id="Meth_ID")
ggplot(hist_plot, aes(value, fill=variable))+geom_histogram(position="identity", color="black", alpha=0.5, binwidth = 0.25)+theme_bw()+
  scale_fill_manual(values=c("lightblue","blue","darkblue"))+facet_wrap(~variable, ncol=1)





# Liner Mixed Effects Model for INCOME per dependent (deinc2dep) MVAL
# LME P value calculation
meta$TWIN_Pair<-as.factor(meta$TWIN_Pair)
Likelihood_Ratio_Test<-lapply(1:nrow(PAWS_Mval), function(CpG){
  metaex<-meta
  metaex$Beta<-PAWS_Mval[CpG,]
  mod_income<-lmer(Beta ~ deinc2dep + Age_Genetic_Collection + as.factor(minor_child)+as.factor(Genetic_cluster)+(1|TWIN_Pair) , data=metaex)
  mod_covar<-lmer(Beta ~ Age_Genetic_Collection + as.factor(minor_child)+as.factor(Genetic_cluster)+(1|TWIN_Pair) , data=metaex)
  Likelihood_Ratio_Test<-anova(mod_income, mod_covar)
  Likelihood_Ratio_Test})
save(Likelihood_Ratio_Test, file="Likelihood_Ratio_Test_lmer_deinc2dep_linear_all_CpGs.RData")

Likelihood_Ratio_Test_Mval<-sapply(1:length(Likelihood_Ratio_Test), function(x) Likelihood_Ratio_Test[[x]][["Pr(>Chisq)"]][2])

## P value distribution
pvalue_dist_deinc2dep<-data.frame(CpG=rownames(PAWS_Mval), Nominal_P=Likelihood_Ratio_Test_Mval)
ggplot(pvalue_dist_deinc2dep, aes(Nominal_P))+geom_histogram(fill="grey90", color="black")+theme_bw()+xlab("Nominal P Value")+
  ylim(0,50000)

#multiple test correction
Multi_test_corr_relaxed<-p.adjust(Likelihood_Ratio_Test_Mval, method = "fdr", n = length(Likelihood_Ratio_Test_Mval))
stat_hits<-as.data.frame(PAWS_Beta)[which(Multi_test_corr_relaxed<=0.1),] # 65 at fdr 0.05; 1445 at fdr 0.1; 4823 at fdr 0.15

# Delta beta cutoff (Mvalue cutoff doesnt really translate)

delbeta<-sapply(1:nrow(PAWS_Beta), function(x) {
  z<-lm(unlist(PAWS_Beta[x,]) ~ meta$deinc2dep + meta$Age_Genetic_Collection + meta$minor_child + meta$Genetic_cluster)
  intercept=z$coefficients[1]
  slope=z$coefficients[2]
  as.numeric((slope*max(meta$deinc2dep, na.rm=T)+intercept)-intercept)
})

bio_hits<-PAWS_Beta[which(abs(delbeta)>=0.1),] # 1540 CpGs
save(Multi_test_corr_relaxed, stat_hits, delbeta, bio_hits, file="deinc2dep_PAWS_hits_Mval.RData")# thresholded on Mval but data in stahits and biohits if betas for plotting



############# VOLCANO

load("deinc2dep_PAWS_hits_Mval.RData")
sta_bio_hits<-stat_hits[which(rownames(stat_hits)%in%rownames(bio_hits)),]# at fdr 0.05; 16 at fdr 0.1;  at fdr 0.15

# psot hoc SNP look
load("SNPCpG.RData")
SnpatCpG<-SNPCpG[which(SNPCpG$SNPCpG!=""),]


## VOLCANO
volcano<-data.frame(Adjusted_Pvalue=Multi_test_corr_relaxed, Delta_Beta=delbeta)

dB<-0.05 #delta beta cutoff
Pv<-0.25 #Pvalue cutoff

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
  scale_color_manual(values=c(muted("red", l=80, c=30),"red",muted("blue", l=70, c=40),"blue","grey"),name="Methylation Change \nWith Income (detotin2)  \nper Dependent")+
  geom_vline(xintercept=c(-dB,dB), color="grey60")+geom_hline(yintercept=-log10(Pv), color="grey60")+
  ylab("Multiple Test Corrected P Value (-log10)")+xlab("Delta Beta")+xlim(-1, 1)+
  theme(axis.text = element_text(size =14, color="black"))+ guides(color = guide_legend(override.aes = list(size = 4)))
