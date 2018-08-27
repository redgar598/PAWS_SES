library(ggplot2)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(reshape)
library(lme4)
library(plyr)


setwd("/big_data/redgar/PAWS/PAWS_final/reviews_response/")
load("/big_data/redgar/PAWS/PAWS_final/combat_PAWS_Beta_norep.RData")
load("/big_data/redgar/PAWS/PAWS_final/PAWS_meta_sentrix_genetic_clusters.RData")

meta<-meta[which(meta$Meth_ID%in%colnames(PAWS_Beta)),]
meta<-meta[match(colnames(PAWS_Beta), meta$Meth_ID),]

Mval<-function(beta) log2(beta/(1-beta))
PAWS_Mval = apply(PAWS_Beta, 1, Mval) # need mvalues for combat
PAWS_Mval = as.data.frame(PAWS_Mval)
PAWS_Mval = t(PAWS_Mval)



# going to do a 90% winsorization of DNAm data and hits

######################################
## Income 
######################################
load(file="/big_data/redgar/PAWS/PAWS_final/Likelihood_Ratio_Test_lmer_deinc2dep_linear_all_CpGs.RData")
load("/big_data/redgar/PAWS/PAWS_final/deinc2dep_PAWS_hits_Mval.RData")
paws_pval<-sapply(1:length(Likelihood_Ratio_Test), function(x) Likelihood_Ratio_Test[[x]][["Pr(>Chisq)"]][2])
paws_data<-data.frame(CpG=rownames(PAWS_Beta), pval=paws_pval, fdr=Multi_test_corr_relaxed, db=delbeta)


## only hits for windorizing
paws_income<-read.csv("/big_data/redgar/PAWS/PAWS_final/CpG_hits_published/Hitsdeinc2dep.csv")
PAWS_Mval_incomehit<-PAWS_Mval[which(rownames(PAWS_Mval)%in%paws_income$x),] #488 CpGs
PAWS_beta_incomehit<-PAWS_Beta[which(rownames(PAWS_Beta)%in%paws_income$x),] #488 CpGs

identical(as.character(meta$Meth_ID),colnames(PAWS_Mval_incomehit))

paws_data_hits<-paws_data[which(paws_data$CpG%in%paws_income$x),]


## winsorize Mvalues
minmval<-sapply(1:nrow(PAWS_Mval_incomehit), function(x) quantile(PAWS_Mval_incomehit[x,], 0.05, na.rm=T))
maxmval<-sapply(1:nrow(PAWS_Mval_incomehit), function(x) quantile(PAWS_Mval_incomehit[x,], 0.95, na.rm=T))

PAWS_Mval_incomehit_windsorized<-lapply(1:nrow(PAWS_Mval_incomehit), function(x){
  y<-PAWS_Mval_incomehit[x,]
  
  y[which(y<minmval[x])]<-minmval[x]
  y[which(y>maxmval[x])]<-maxmval[x]
  
  y
})

PAWS_Mval_incomehit_windsorized<-do.call(rbind, PAWS_Mval_incomehit_windsorized)
rownames(PAWS_Mval_incomehit_windsorized)<-rownames(PAWS_Mval_incomehit)

# Liner Mixed Effects Model for INCOME per dependent (deinc2dep) MVAL
# LME P value calculation
meta$TWIN_Pair<-as.factor(meta$TWIN_Pair)
Likelihood_Ratio_Test<-lapply(1:nrow(PAWS_Mval_incomehit_windsorized), function(CpG){
  metaex<-meta
  metaex$Beta<-PAWS_Mval_incomehit_windsorized[CpG,]
  mod_income<-lmer(Beta ~ deinc2dep + Age_Genetic_Collection + as.factor(minor_child)+as.factor(Genetic_cluster)+(1|TWIN_Pair) , data=metaex)
  mod_covar<-lmer(Beta ~ Age_Genetic_Collection + as.factor(minor_child)+as.factor(Genetic_cluster)+(1|TWIN_Pair) , data=metaex)
  Likelihood_Ratio_Test<-anova(mod_income, mod_covar)
  Likelihood_Ratio_Test})




#multiple test correction Does not actually make sense to do FDR on candidates as it is dependent on p value distribution
pval<-sapply(1:length(Likelihood_Ratio_Test), function(x) Likelihood_Ratio_Test[[x]][["Pr(>Chisq)"]][2])
Multi_test_corr_relaxed<-p.adjust(pval, method = "fdr", n = nrow(PAWS_Mval))
stat_hits<-as.data.frame(PAWS_beta_incomehit)[which(Multi_test_corr_relaxed<=0.2),] # 15/488



# Delta beta cutoff (Mvalue cutoff doesnt really translate)
# winsroize betas

minbeta<-sapply(1:nrow(PAWS_beta_incomehit), function(x) quantile(PAWS_beta_incomehit[x,], 0.05, na.rm=T))
maxbeta<-sapply(1:nrow(PAWS_beta_incomehit), function(x) quantile(PAWS_beta_incomehit[x,], 0.95, na.rm=T))

PAWS_beta_incomehit_windsorized<-lapply(1:nrow(PAWS_beta_incomehit), function(x){
  y<-PAWS_beta_incomehit[x,]
  
  y[which(y<minbeta[x])]<-minbeta[x]
  y[which(y>maxbeta[x])]<-maxbeta[x]
  
  y
})

PAWS_beta_incomehit_windsorized<-do.call(rbind, PAWS_beta_incomehit_windsorized)
rownames(PAWS_beta_incomehit_windsorized)<-rownames(PAWS_beta_incomehit)

delbeta<-sapply(1:nrow(PAWS_beta_incomehit_windsorized), function(x) {
  z<-lm(unlist(PAWS_beta_incomehit_windsorized[x,]) ~ meta$deinc2dep + meta$Age_Genetic_Collection + meta$minor_child + meta$Genetic_cluster)
  intercept=z$coefficients[1]
  slope=z$coefficients[2]
  as.numeric((slope*max(meta$deinc2dep, na.rm=T)+intercept)-intercept)
})

bio_hits<-PAWS_beta_incomehit[which(abs(delbeta)>=0.05),] # 348/488 CpGs

paws_winsor_data<-data.frame(CpG=rownames(PAWS_beta_incomehit), pval_windsor=pval, fdr_windsor=Multi_test_corr_relaxed, db_windsor=delbeta)

cor(paws_data_hits$pval, paws_winsor_data$pval)
cor(paws_data_hits$db, paws_winsor_data$db)


## number passing winsor
nrow(paws_winsor_data[which(paws_winsor_data$fdr<=0.2 & abs(paws_winsor_data$db)>0.05),]) #12
nrow(paws_winsor_data[which(paws_winsor_data$pval<=0.005 & abs(paws_winsor_data$db)>0.05),]) #215
nrow(paws_winsor_data[which(paws_winsor_data$pval<=0.014 & abs(paws_winsor_data$db)>0.05),]) #335 same nom pval as fdr=0.2 in og data


paws_data_plt<-merge(paws_data_hits, paws_winsor_data, by="CpG")

ggplot(paws_data_plt, aes(db,db_windsor))+geom_hline(yintercept=0,linetype="longdash", color="lightgrey")+
  geom_vline(xintercept=0, linetype="longdash",  color="lightgrey")+
  geom_point(shape=21, color="black", fill="grey",size=2)+
  theme_bw()+xlab("Delta Beta Original DNAm Values")+ylab("Delta Beta with 90% winsorization of DNAm")+
  ylim(-0.25, 0.25)+xlim(-0.25,0.25)


paws_data_plt$color<-sapply(1:nrow(paws_data_plt), function(x){
  if(paws_data_plt$pval_windsor[x]<=0.014){"sig"}else{"not sig"}
})

ggplot(paws_data_plt, aes(-log10(pval),-log10(pval_windsor), fill=color))+
  geom_hline(yintercept=-log10(max(paws_data_plt$pval)),linetype="longdash", color="lightgrey")+
  geom_vline(xintercept=-log10(max(paws_data_plt$pval)), linetype="longdash",  color="lightgrey")+
  geom_point(shape=21, color="black",size=2)+scale_fill_manual(values=c("grey","cornflowerblue"))+
  theme_bw()+xlab("-log10 P Value Original DNAm Values ")+ylab("-log10 P Value with 90% winsorization of DNAm")+
  ylim(0,10)+xlim(0,10)













######################################
## education 
######################################
load(file="/big_data/redgar/PAWS/PAWS_final/Likelihood_Ratio_Test_lmer_HGHEDLV2_PAWS_Mval_genetic_cluster_and_twin_fixed_allFactors.RData")
load("/big_data/redgar/PAWS/PAWS_final/HGHEDLV2_PAWS_hits_Mval.RData")
paws_pval<-sapply(1:length(Likelihood_Ratio_Test), function(x) Likelihood_Ratio_Test[[x]][["Pr(>Chisq)"]][2])
paws_data<-data.frame(CpG=rownames(PAWS_Beta), pval=paws_pval, fdr=Multi_test_corr_relaxed, db=delbeta)


## only hits for windorizing
paws_edu<-read.csv("/big_data/redgar/PAWS/PAWS_final/CpG_hits_published/HitsHGHEDLV2.csv")
PAWS_Mval_eduhit<-PAWS_Mval[which(rownames(PAWS_Mval)%in%paws_edu$x),] #354 CpGs
PAWS_beta_eduhit<-PAWS_Beta[which(rownames(PAWS_Beta)%in%paws_edu$x),] #354 CpGs

identical(as.character(meta$Meth_ID),colnames(PAWS_Mval_eduhit))

paws_data_hits<-paws_data[which(paws_data$CpG%in%paws_edu$x),]


## winsorize Mvalues
minmval<-sapply(1:nrow(PAWS_Mval_eduhit), function(x) quantile(PAWS_Mval_eduhit[x,], 0.05, na.rm=T))
maxmval<-sapply(1:nrow(PAWS_Mval_eduhit), function(x) quantile(PAWS_Mval_eduhit[x,], 0.95, na.rm=T))

PAWS_Mval_eduhit_windsorized<-lapply(1:nrow(PAWS_Mval_eduhit), function(x){
  y<-PAWS_Mval_eduhit[x,]
  
  y[which(y<minmval[x])]<-minmval[x]
  y[which(y>maxmval[x])]<-maxmval[x]
  
  y
})

PAWS_Mval_eduhit_windsorized<-do.call(rbind, PAWS_Mval_eduhit_windsorized)
rownames(PAWS_Mval_eduhit_windsorized)<-rownames(PAWS_Mval_eduhit)

# Liner Mixed Effects Model for INCOME per dependent (deinc2dep) MVAL
# LME P value calculation
meta$TWIN_Pair<-as.factor(meta$TWIN_Pair)
Likelihood_Ratio_Test<-lapply(1:nrow(PAWS_Mval_eduhit_windsorized), function(CpG){
  metaex<-meta
  metaex$Beta<-PAWS_Mval_eduhit_windsorized[CpG,]
  mod_edu<-lmer(Beta ~ HGHEDLV2 + as.factor(minor_child)+as.factor(Genetic_cluster)+(1|TWIN_Pair) , data=metaex)
  mod_covar<-lmer(Beta ~ as.factor(minor_child)+as.factor(Genetic_cluster)+(1|TWIN_Pair) , data=metaex)
  Likelihood_Ratio_Test<-anova(mod_edu, mod_covar)
  Likelihood_Ratio_Test})



#multiple test correction Does not actually make sense to do FDR on candidates as it is dependent on p value distribution
pval<-sapply(1:length(Likelihood_Ratio_Test), function(x) Likelihood_Ratio_Test[[x]][["Pr(>Chisq)"]][2])
Multi_test_corr_relaxed<-p.adjust(pval, method = "fdr", n = nrow(PAWS_Mval))
stat_hits<-as.data.frame(PAWS_beta_eduhit)[which(Multi_test_corr_relaxed<=0.2),] # 0/354



# Delta beta cutoff (Mvalue cutoff doesnt really translate)
# winsroize betas

minbeta<-sapply(1:nrow(PAWS_beta_eduhit), function(x) quantile(PAWS_beta_eduhit[x,], 0.05, na.rm=T))
maxbeta<-sapply(1:nrow(PAWS_beta_eduhit), function(x) quantile(PAWS_beta_eduhit[x,], 0.95, na.rm=T))

PAWS_beta_eduhit_windsorized<-lapply(1:nrow(PAWS_beta_eduhit), function(x){
  y<-PAWS_beta_eduhit[x,]
  
  y[which(y<minbeta[x])]<-minbeta[x]
  y[which(y>maxbeta[x])]<-maxbeta[x]
  
  y
})

PAWS_beta_eduhit_windsorized<-do.call(rbind, PAWS_beta_eduhit_windsorized)
rownames(PAWS_beta_eduhit_windsorized)<-rownames(PAWS_beta_eduhit)

delbeta<-sapply(1:nrow(PAWS_beta_eduhit_windsorized), function(x) {
  z<-lm(unlist(PAWS_beta_eduhit_windsorized[x,]) ~ meta$HGHEDLV2 + meta$minor_child+meta$Genetic_cluster)
  intercept=z$coefficients[1]
  slope=z$coefficients[2]
  as.numeric((slope*6+intercept)-intercept)})


bio_hits<-PAWS_beta_eduhit[which(abs(delbeta)>=0.05),] # 251/354 CpGs

paws_winsor_data<-data.frame(CpG=rownames(PAWS_beta_eduhit), pval_windsor=pval, fdr_windsor=Multi_test_corr_relaxed, db_windsor=delbeta)

cor(paws_data_hits$pval, paws_winsor_data$pval)
cor(paws_data_hits$db, paws_winsor_data$db)


## number passing winsor
nrow(paws_winsor_data[which(paws_winsor_data$fdr<=0.2 & abs(paws_winsor_data$db)>0.05),]) #0
nrow(paws_winsor_data[which(paws_winsor_data$pval<=0.005 & abs(paws_winsor_data$db)>0.05),]) #151
nrow(paws_winsor_data[which(paws_winsor_data$pval<=0.013 & abs(paws_winsor_data$db)>0.05),]) #238 same nom pval as fdr=0.2 in og data


paws_data_plt<-merge(paws_data_hits, paws_winsor_data, by="CpG")

ggplot(paws_data_plt, aes(db,db_windsor))+geom_hline(yintercept=0,linetype="longdash", color="lightgrey")+
  geom_vline(xintercept=0, linetype="longdash",  color="lightgrey")+
  geom_point(shape=21, color="black", fill="grey",size=2)+
  theme_bw()+xlab("Delta Beta Original DNAm Values")+ylab("Delta Beta with 90% winsorization of DNAm")+
  ylim(-0.25, 0.25)+xlim(-0.25,0.25)


paws_data_plt$color<-sapply(1:nrow(paws_data_plt), function(x){
  if(paws_data_plt$pval_windsor[x]<=0.013){"sig"}else{"not sig"}
})

ggplot(paws_data_plt, aes(-log10(pval),-log10(pval_windsor), fill=color))+
  geom_hline(yintercept=-log10(max(paws_data_plt$pval)),linetype="longdash", color="lightgrey")+
  geom_vline(xintercept=-log10(max(paws_data_plt$pval)), linetype="longdash",  color="lightgrey")+
  geom_point(shape=21, color="black",size=2)+scale_fill_manual(values=c("grey","cornflowerblue"))+
  theme_bw()+xlab("-log10 P Value Original DNAm Values ")+ylab("-log10 P Value with 90% winsorization of DNAm")+
  ylim(0,10)+xlim(0,10)
















######################################
## adversity 
######################################
load(file="/big_data/redgar/PAWS/PAWS_final/Likelihood_Ratio_Test_lmer_FADV_6_PAWS_Mval_genetic_cluster_and_twin_fixed_allFactors.RData")
load("/big_data/redgar/PAWS/PAWS_final/FADV_6_PAWS_hits_Mval.RData")
paws_pval<-sapply(1:length(Likelihood_Ratio_Test), function(x) Likelihood_Ratio_Test[[x]][["Pr(>Chisq)"]][2])
paws_data<-data.frame(CpG=rownames(PAWS_Beta), pval=paws_pval, fdr=Multi_test_corr_relaxed, db=delbeta)


## only hits for windorizing
paws_adv<-read.csv("/big_data/redgar/PAWS/PAWS_final/CpG_hits_published/HitsFADV_6.csv")
PAWS_Mval_advhit<-PAWS_Mval[which(rownames(PAWS_Mval)%in%paws_adv$x),] #102 CpGs
PAWS_beta_advhit<-PAWS_Beta[which(rownames(PAWS_Beta)%in%paws_adv$x),] #102 CpGs

identical(as.character(meta$Meth_ID),colnames(PAWS_Mval_advhit))

paws_data_hits<-paws_data[which(paws_data$CpG%in%paws_adv$x),]


## winsorize Mvalues
minmval<-sapply(1:nrow(PAWS_Mval_advhit), function(x) quantile(PAWS_Mval_advhit[x,], 0.05, na.rm=T))
maxmval<-sapply(1:nrow(PAWS_Mval_advhit), function(x) quantile(PAWS_Mval_advhit[x,], 0.95, na.rm=T))

PAWS_Mval_advhit_windsorized<-lapply(1:nrow(PAWS_Mval_advhit), function(x){
  y<-PAWS_Mval_advhit[x,]
  
  y[which(y<minmval[x])]<-minmval[x]
  y[which(y>maxmval[x])]<-maxmval[x]
  
  y
})

PAWS_Mval_advhit_windsorized<-do.call(rbind, PAWS_Mval_advhit_windsorized)
rownames(PAWS_Mval_advhit_windsorized)<-rownames(PAWS_Mval_advhit)

# Liner Mixed Effects Model for INCOME per dependent (deinc2dep) MVAL
# LME P value calculation
meta$TWIN_Pair<-as.factor(meta$TWIN_Pair)
Likelihood_Ratio_Test<-lapply(1:nrow(PAWS_Mval_advhit_windsorized), function(CpG){
  metaex<-meta
  metaex$Beta<-PAWS_Mval_advhit_windsorized[CpG,]
  mod_adversity<-lmer(Beta ~ FADV_6 +as.factor(Genetic_cluster)+(1|TWIN_Pair) , data=metaex)
  mod_covar<-lmer(Beta ~ as.factor(Genetic_cluster)+(1|TWIN_Pair) , data=metaex)
  Likelihood_Ratio_Test<-anova(mod_adversity, mod_covar)
  Likelihood_Ratio_Test})



#multiple test correction Does not actually make sense to do FDR on candidates as it is dependent on p value distribution
pval<-sapply(1:length(Likelihood_Ratio_Test), function(x) Likelihood_Ratio_Test[[x]][["Pr(>Chisq)"]][2])
Multi_test_corr_relaxed<-p.adjust(pval, method = "fdr", n = nrow(PAWS_Mval))
stat_hits<-as.data.frame(PAWS_beta_advhit)[which(Multi_test_corr_relaxed<=0.2),] # 1/102



# Delta beta cutoff (Mvalue cutoff doesnt really translate)
# winsroize betas

minbeta<-sapply(1:nrow(PAWS_beta_advhit), function(x) quantile(PAWS_beta_advhit[x,], 0.05, na.rm=T))
maxbeta<-sapply(1:nrow(PAWS_beta_advhit), function(x) quantile(PAWS_beta_advhit[x,], 0.95, na.rm=T))

PAWS_beta_advhit_windsorized<-lapply(1:nrow(PAWS_beta_advhit), function(x){
  y<-PAWS_beta_advhit[x,]
  
  y[which(y<minbeta[x])]<-minbeta[x]
  y[which(y>maxbeta[x])]<-maxbeta[x]
  
  y
})

PAWS_beta_advhit_windsorized<-do.call(rbind, PAWS_beta_advhit_windsorized)
rownames(PAWS_beta_advhit_windsorized)<-rownames(PAWS_beta_advhit)

delbeta<-sapply(1:nrow(PAWS_beta_advhit_windsorized), function(x) {
  z<-lm(unlist(PAWS_beta_advhit_windsorized[x,]) ~ meta$FADV_6 +meta$Genetic_cluster)
  intercept=z$coefficients[1]
  slope=z$coefficients[2]
  as.numeric((slope*max(meta$FADV_6)+intercept)-intercept)})


bio_hits<-PAWS_beta_advhit[which(abs(delbeta)>=0.05),] # 13/102 CpGs

paws_winsor_data<-data.frame(CpG=rownames(PAWS_beta_advhit), pval_windsor=pval, fdr_windsor=Multi_test_corr_relaxed, db_windsor=delbeta)

cor(paws_data_hits$pval, paws_winsor_data$pval)
cor(paws_data_hits$db, paws_winsor_data$db)


## number passing winsor
nrow(paws_winsor_data[which(paws_winsor_data$fdr<=0.2 & abs(paws_winsor_data$db)>0.05),]) #0
nrow(paws_winsor_data[which(paws_winsor_data$pval<=0.005 & abs(paws_winsor_data$db)>0.05),]) #13
nrow(paws_winsor_data[which(paws_winsor_data$pval<=0.002 & abs(paws_winsor_data$db)>0.05),]) #11 same nom pval as fdr=0.2 in og data


paws_data_plt<-merge(paws_data_hits, paws_winsor_data, by="CpG")

ggplot(paws_data_plt, aes(db,db_windsor))+geom_hline(yintercept=0,linetype="longdash", color="lightgrey")+
  geom_vline(xintercept=0, linetype="longdash",  color="lightgrey")+
  geom_point(shape=21, color="black", fill="grey",size=2)+
  theme_bw()+xlab("Delta Beta Original DNAm Values")+ylab("Delta Beta with 90% winsorization of DNAm")+
  ylim(-0.25, 0.25)+xlim(-0.25,0.25)


paws_data_plt$color<-sapply(1:nrow(paws_data_plt), function(x){
  if(paws_data_plt$pval_windsor[x]<=0.002){"sig"}else{"not sig"}
})

ggplot(paws_data_plt, aes(-log10(pval),-log10(pval_windsor), fill=color))+
  geom_hline(yintercept=-log10(max(paws_data_plt$pval)),linetype="longdash", color="lightgrey")+
  geom_vline(xintercept=-log10(max(paws_data_plt$pval)), linetype="longdash",  color="lightgrey")+
  geom_point(shape=21, color="black",size=2)+scale_fill_manual(values=c("grey","cornflowerblue"))+
  theme_bw()+xlab("-log10 P Value Original DNAm Values ")+ylab("-log10 P Value with 90% winsorization of DNAm")+
  ylim(0,10)+xlim(0,10)






#######################################
### plot the representative CpGs
#######################################
# FADV_6_genes<-read.csv("/big_data/redgar/PAWS/PAWS_final/CpG_hits_published/FADV_6_genes_FDR02_DB05.csv")
# HGHEDLV2_genes<-read.csv("/big_data/redgar/PAWS/PAWS_final/CpG_hits_published/HGHEDLV2_genes_FDR02_DB05.csv")
# deinc2dep_genes<-read.csv("/big_data/redgar/PAWS/PAWS_final/CpG_hits_published/deinc2dep_genes_FDR02_DB05.csv")


### rep CpG list for gene list

adv_repcpg<-c("cg06588556","cg21676617","cg26053697","cg10581375")
edu_rep_cpg<-c("cg02918577", "cg14245471", "cg20944157", "cg26511075")
inc_repcpg<-c("cg08200582", "cg27149937", "cg17632375", "cg21502834")


## Gene CpG number adjustment
load("/big_data/redgar/PAWS/PAWS_final/Gene_CpG_Relations_updatejune2015.RData")
Gene_CpG_Relations_update$gene<-as.character(Gene_CpG_Relations_update$gene)
Overrep<-as.data.frame(tapply(Gene_CpG_Relations_update$Probe_ID, Gene_CpG_Relations_update$gene, length))
Overrep$Gene<-rownames(Overrep)
colnames(Overrep)<-c("CpG_number", "Gene")
Overrep<-Overrep[which(Overrep$Gene!="None"),]
mean(Overrep$CpG_number, na.rm=T)# 25
Overrep$Enrichment_fromAverage<-Overrep$CpG_number/mean(Overrep$CpG_number, na.rm=T)

## Gene summaries
load("/big_data/redgar/PAWS/PAWS_final/Price_annotation.RData")
annotation$CpG<-rownames(annotation)

Format_gene_table<-function(Gene_CpG_Relations_update_subset){
  print(paste("CpGs Associated: ", length(unique(Gene_CpG_Relations_update_subset$Probe_ID)), sep=""))
  print(paste("Genes Associated: ", length(unique(Gene_CpG_Relations_update_subset$gene)), sep=""))
  Overrep_subset<-as.data.frame(tapply(Gene_CpG_Relations_update_subset$Probe_ID, Gene_CpG_Relations_update_subset$gene, length))
  Overrep_subset$Gene<-rownames(Overrep_subset)
  colnames(Overrep_subset)<-c("CpG_number", "Gene")
  Overrep_subset<-Overrep_subset[which(Overrep_subset$Gene!="None"),]
  Overrep_subset_merge<-merge(Overrep_subset, Overrep, by="Gene")
  colnames(Overrep_subset_merge)<-c("Gene","CpG_Associated","CpG_in_Gene", "Enrichment_fromAverage")
  Overrep_subset_merge$Suprise<-Overrep_subset_merge$CpG_Associated/Overrep_subset_merge$Enrichment_fromAverage
  Gene_table<-merge(Gene_CpG_Relations_update_subset, Overrep_subset_merge, by.x="gene", by.y="Gene")
  Gene_table<-merge(Gene_table, annotation[,c(49,50,58)], by.x="Probe_ID", by.y="CpG")
  
  Gene_table<-Gene_table[,c(2,7,8,9,10,1,4,3,6,5,11,12)]
  Gene_table<-Gene_table[order(-Gene_table$Suprise, Gene_table$gene),]
  Gene_table}



#Genes 
FADV_6<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%adv_repcpg),]
FADV_6<-FADV_6[!duplicated(FADV_6),]
FADV_6<-FADV_6[!duplicated(FADV_6[,c(1,4)]),]#remove duplicate CpG to gene associations
FADV_6_genes<-Format_gene_table(FADV_6)

deinc2dep<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%inc_repcpg),]
deinc2dep<-deinc2dep[!duplicated(deinc2dep),]
deinc2dep<-deinc2dep[!duplicated(deinc2dep[,c(1,4)]),]#remove duplicate CpG to gene associations
deinc2dep_genes<-Format_gene_table(deinc2dep)

HGHEDLV2<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%edu_rep_cpg),]
HGHEDLV2<-HGHEDLV2[!duplicated(HGHEDLV2),]
HGHEDLV2<-HGHEDLV2[!duplicated(HGHEDLV2[,c(1,4)]),]#remove duplicate CpG to gene associations
HGHEDLV2_genes<-Format_gene_table(HGHEDLV2)



## FUNCTIONS
CpG_plot_inc<-function(CpG_list){
  CpGs<-CpG_list
  plot_beta<-as.data.frame(PAWS_Beta[which(rownames(PAWS_Beta)%in%CpGs),])
  minbeta<-sapply(1:nrow(plot_beta), function(x) quantile(plot_beta[x,], 0.05, na.rm=T))
  maxbeta<-sapply(1:nrow(plot_beta), function(x) quantile(plot_beta[x,], 0.95, na.rm=T))
  
  plot_beta_win<-lapply(1:nrow(plot_beta), function(x){
    y<-plot_beta[x,]
    y[which(y<minbeta[x])]<-minbeta[x]
    y[which(y>maxbeta[x])]<-maxbeta[x]
    y  })
  
  plot_beta_win<-do.call(rbind, plot_beta_win)
  rownames(plot_beta_win)<-rownames(plot_beta)
  plot_beta_win$CpG<-rownames(plot_beta_win)
  plot<-melt(plot_beta_win, id="CpG")
  
  plot_CpG<-merge(plot, meta, by.x="variable", by.y="Meth_ID")
  plot_CpG<-merge(plot_CpG, deinc2dep_genes, by.x="CpG",by.y="Probe_ID")
  plot_CpG$Label<-paste(plot_CpG$gene," (", plot_CpG$CpG,")", sep="")
  
  
  lm_eqn_income <- function(df){
    df$Genetic_cluster<-as.factor(df$Genetic_cluster)
    m <- lm(value ~ deinc2dep + Age_Genetic_Collection + minor_child + Genetic_cluster, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2), 
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
  }
  
  eq <- ddply(plot_CpG,.(Label),lm_eqn_income)
  
  
  ggplot(plot_CpG, aes(deinc2dep, value))+
    geom_point(aes(color=as.factor(Genetic_cluster)),shape=19)+theme_bw()+
    stat_smooth(method = "lm", size = 0.75, se=F, color="grey20")+
    facet_wrap(~Label)+ylab("Beta Value")+ylim(0,1)+
    scale_color_manual(values=colorRampPalette(brewer.pal(11,"RdYlGn")[1:8])(4), name="PLINK Genetic Cluster")+
    ylab("Methylation Level (Beta Value)")+xlab("Income-Per-Dependent")+
    geom_text(data=eq,aes(x = 2, y = 0.98,label=V1), parse = TRUE, inherit.aes=FALSE, size=3, color="grey40")}

CpG_plot_adv<-function(CpG_list){
  CpGs<-CpG_list
  plot_beta<-as.data.frame(PAWS_Beta[which(rownames(PAWS_Beta)%in%CpGs),])
  
  minbeta<-sapply(1:nrow(plot_beta), function(x) quantile(plot_beta[x,], 0.05, na.rm=T))
  maxbeta<-sapply(1:nrow(plot_beta), function(x) quantile(plot_beta[x,], 0.95, na.rm=T))
  
  plot_beta_win<-lapply(1:nrow(plot_beta), function(x){
    y<-plot_beta[x,]
    y[which(y<minbeta[x])]<-minbeta[x]
    y[which(y>maxbeta[x])]<-maxbeta[x]
    y  })
  
  plot_beta_win<-do.call(rbind, plot_beta_win)
  rownames(plot_beta_win)<-rownames(plot_beta)
  
  plot_beta_win$CpG<-rownames(plot_beta_win)
  plot<-melt(plot_beta_win, id="CpG")
  plot_CpG<-merge(plot, meta, by.x="variable", by.y="Meth_ID")
  plot_CpG<-merge(plot_CpG, FADV_6_genes, by.x="CpG",by.y="Probe_ID")
  plot_CpG$Label<-paste(plot_CpG$gene," (", plot_CpG$CpG,")", sep="")
  
  
  lm_eqn_adv <- function(df){
    df$Genetic_cluster<-as.factor(df$Genetic_cluster)
    m <- lm(value ~ FADV_6 + Genetic_cluster, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2), 
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
  }
  
  eq <- ddply(plot_CpG,.(Label),lm_eqn_adv)
  
  
  ggplot(plot_CpG, aes(FADV_6, value))+
    geom_point(aes(color=as.factor(Genetic_cluster)),shape=19)+theme_bw()+
    stat_smooth(method = "lm", size = 0.75, se=F, color="grey20")+
    facet_wrap(~Label)+ylab("Beta Value")+ylim(0,1)+
    scale_color_manual(values=colorRampPalette(brewer.pal(11,"RdYlGn")[1:8])(4), name="PLINK Genetic Cluster")+
    ylab("Methylation Level (Beta Value)")+xlab("Family Adversity")+
    geom_text(data=eq,aes(x = 0, y = 0.98,label=V1), parse = TRUE, inherit.aes=FALSE, size=3, color="grey40")}

CpG_plot_edu<-function(CpG_list){
  CpGs<-CpG_list
  plot_beta<-as.data.frame(PAWS_Beta[which(rownames(PAWS_Beta)%in%CpGs),])
  
  minbeta<-sapply(1:nrow(plot_beta), function(x) quantile(plot_beta[x,], 0.05, na.rm=T))
  maxbeta<-sapply(1:nrow(plot_beta), function(x) quantile(plot_beta[x,], 0.95, na.rm=T))
  
  plot_beta_win<-lapply(1:nrow(plot_beta), function(x){
    y<-plot_beta[x,]
    y[which(y<minbeta[x])]<-minbeta[x]
    y[which(y>maxbeta[x])]<-maxbeta[x]
    y  })
  
  plot_beta_win<-do.call(rbind, plot_beta_win)
  rownames(plot_beta_win)<-rownames(plot_beta)
  
  
  plot_beta_win$CpG<-rownames(plot_beta_win)
  plot<-melt(plot_beta_win, id="CpG")
  plot_CpG<-merge(plot, meta, by.x="variable", by.y="Meth_ID")
  plot_CpG<-merge(plot_CpG, HGHEDLV2_genes, by.x="CpG",by.y="Probe_ID")
  plot_CpG$Label<-paste(plot_CpG$gene," (", plot_CpG$CpG,")", sep="")
  
  lm_eqn_edu <- function(df){
    df$Genetic_cluster<-as.factor(df$Genetic_cluster)
    m <- lm(value ~ HGHEDLV2 +minor_child + Genetic_cluster, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                     list(a = format(coef(m)[1], digits = 2), 
                          b = format(coef(m)[2], digits = 2), 
                          r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));                 
  }
  
  eq <- ddply(plot_CpG,.(Label),lm_eqn_edu)
  
  ggplot(plot_CpG, aes(HGHEDLV2, value))+
    geom_point(aes(color=as.factor(Genetic_cluster)),shape=19)+theme_bw()+
    stat_smooth(method = "lm", size = 0.75, se=F, color="grey20")+
    facet_wrap(~Label)+ylab("Beta Value")+ylim(0,1)+
    scale_color_manual(values=colorRampPalette(brewer.pal(11,"RdYlGn")[1:8])(4), name="PLINK Genetic Cluster")+
    ylab("Methylation Level (Beta Value)")+xlab("Parental Education")+
    geom_text(data=eq,aes(x = 3, y = 0.98,label=V1), parse = TRUE, inherit.aes=FALSE, size=3, color="grey40")}



adv<-CpG_plot_adv(adv_repcpg)
edu<-CpG_plot_edu(edu_rep_cpg)
inc<-CpG_plot_inc(inc_repcpg)

grid.arrange(inc,edu,adv, ncol=3)

