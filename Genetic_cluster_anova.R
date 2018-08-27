#PAWS Genetic Cluster ANOVA

#Using normalized combatted data regressions can now be applied

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



# ANOVA Genetic Cluster
Genetic_cluster_aov<-sapply(1:nrow(PAWS_Mval), function(CpG){
  x<-aov(PAWS_Mval[CpG,]~as.factor(meta$Genetic_cluster)+meta$FADV_6+meta$deinc2dep+meta$HGHEDLV2)
  summary(x)[[1]][["Pr(>F)"]][1]})

pvalue_dist<-data.frame(CpG=rownames(PAWS_Mval), Nominal_P=Genetic_cluster_aov)
ggplot(pvalue_dist, aes(Nominal_P))+geom_histogram(fill="grey90", color="black")+theme_bw()+xlab("Nominal P Value")


#multiple test correction
Multi_test_corr<-p.adjust(Genetic_cluster_aov, method = "fdr", n = length(Genetic_cluster_aov))
stat_hits<-as.data.frame(PAWS_Beta)[which(Multi_test_corr<=0.1),] # 830 at fdr 0.05; 7046 at fdr 0.1; 16224 at fdr 0.15



# Delta beta cutoff (Mvalue cutoff doesnt really translate)
delbeta<-sapply(1:nrow(PAWS_Beta), function(x) {
  GC_means<-tapply(PAWS_Beta[x,], meta$Genetic_cluster, mean, na.rm=T)
  max(GC_means)-min(GC_means)})

bio_hits<-PAWS_Beta[which(abs(delbeta)>=0.1),] # 1450 CpGs
save(Multi_test_corr, stat_hits, delbeta, bio_hits, file="Genetic_cluster_hits.RData")# thresholded on Mval but data in stahits and biohits if betas for plotting


# PLot hits
fdr<-0.0001
dB<-0.25
stat_hits<-as.data.frame(PAWS_Beta)[which(Multi_test_corr<=fdr),] 
bio_hits<-PAWS_Beta[which(abs(delbeta)>=dB),] 
sta_bio_hits<-stat_hits[which(rownames(stat_hits)%in%rownames(bio_hits)),]
sta_bio_hits<-sta_bio_hits[which(!(rownames(sta_bio_hits)%in%rownames(SnpatCpG))),]

sta_bio_hits<-sta_bio_hits[sample(1:nrow(sta_bio_hits), 4),]
sta_bio_hits$CpG<-rownames(sta_bio_hits)
sta_bio_hits_melt<-melt(sta_bio_hits, id="CpG")

sta_bio_hits_melt_merge<-merge(sta_bio_hits_melt, meta, by.x="variable", by.y="Meth_ID")

ggplot(sta_bio_hits_melt_merge, aes(as.factor(Genetic_cluster), value, fill=as.factor(Genetic_cluster)))+geom_boxplot()+theme_bw()+
  scale_fill_manual(values=colorRampPalette(brewer.pal(11,"RdYlGn")[1:8])(4), name="Genetic Cluster")+facet_wrap(~CpG)

