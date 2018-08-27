library(ggplot2)
library(grid)
library(gridExtra)
library(reshape)
library(RColorBrewer)


setwd("/big_data/redgar/PAWS/PAWS_final")

load("combat_PAWS_Beta_norep.RData")
load("PAWS_meta_sentrix_genetic_clusters.RData")
meta<-meta[which(meta$Meth_ID%in%colnames(PAWS_Beta)),]
meta<-meta[match(colnames(PAWS_Beta), meta$Meth_ID),]

load("SNPCpG.RData")
SnpatCpG<-SNPCpG[which(SNPCpG$SNPCpG!=""),]

#### Get gene hits
load("Genetic_cluster_hits.RData")
fdr<-0.25
dB<-0.05
stat_hits<-as.data.frame(PAWS_Beta)[which(Multi_test_corr<=fdr),] 
bio_hits<-PAWS_Beta[which(abs(delbeta)>=dB),] 
sta_bio_hits<-stat_hits[which(rownames(stat_hits)%in%rownames(bio_hits)),]
sta_bio_hits_Genetic_cluster<-sta_bio_hits[which(!(rownames(sta_bio_hits)%in%rownames(SnpatCpG))),] ## 8445 CpGs

## gene table stuff
load("Gene_CpG_Relations_updatejune2015.RData")
Gene_CpG_Relations_update$gene<-as.character(Gene_CpG_Relations_update$gene)
Overrep<-as.data.frame(tapply(Gene_CpG_Relations_update$Probe_ID, Gene_CpG_Relations_update$gene, length))
Overrep$Gene<-rownames(Overrep)
colnames(Overrep)<-c("CpG_number", "Gene")
Overrep<-Overrep[which(Overrep$Gene!="None"),]
mean(Overrep$CpG_number, na.rm=T)# 25
Overrep$Enrichment_fromAverage<-Overrep$CpG_number/mean(Overrep$CpG_number, na.rm=T)

## Gene summaries
load("Price_annotation.RData")
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
  pval<-data.frame(CpG=rownames(PAWS_Beta)[which(rownames(PAWS_Beta)%in%Gene_table$Probe_ID)],
                   corr_pval=Multi_test_corr_relaxed[which(rownames(PAWS_Beta)%in%Gene_table$Probe_ID)],
                   db=delbeta[which(rownames(PAWS_Beta)%in%Gene_table$Probe_ID)])
  Gene_table<-merge(Gene_table, pval, by.x="Probe_ID", by.y="CpG")
  Gene_table<-Gene_table[,c(2,7,8,9,10,1,4,3,6,5,11,12,13,14)]
  Gene_table<-Gene_table[order(-Gene_table$Suprise, Gene_table$gene),]
  Gene_table}



### Get gene hits (did the pyro selection at the relaxed FDR level)
fdr<-0.25
dB<-0.05



load("FADV_6_PAWS_hits_Mval.RData")
stat_hits<-as.data.frame(PAWS_Beta)[which(Multi_test_corr_relaxed<=fdr),] 
bio_hits<-PAWS_Beta[which(abs(delbeta)>=dB),] 
sta_bio_hits<-stat_hits[which(rownames(stat_hits)%in%rownames(bio_hits)),]
FADV_6_sta_bio_hits<-sta_bio_hits[which(!(rownames(sta_bio_hits)%in%rownames(SnpatCpG))),]
FADV_6_sta_bio_hits<-FADV_6_sta_bio_hits[which(!(rownames(FADV_6_sta_bio_hits)%in%rownames(sta_bio_hits_Genetic_cluster))),]
#Genes 
FADV_6<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%rownames(FADV_6_sta_bio_hits)),]
FADV_6<-FADV_6[!duplicated(FADV_6),]
FADV_6<-FADV_6[!duplicated(FADV_6[,c(1,4)]),]#remove duplicate CpG to gene associations

FADV_6_genes<-Format_gene_table(FADV_6)

load("HGHEDLV2_PAWS_hits_Mval.RData")
stat_hits<-as.data.frame(PAWS_Beta)[which(Multi_test_corr_relaxed<=fdr),] 
bio_hits<-PAWS_Beta[which(abs(delbeta)>=dB),] 
sta_bio_hits<-stat_hits[which(rownames(stat_hits)%in%rownames(bio_hits)),]
HGHEDLV2_sta_bio_hits<-sta_bio_hits[which(!(rownames(sta_bio_hits)%in%rownames(SnpatCpG))),]
HGHEDLV2_sta_bio_hits<-HGHEDLV2_sta_bio_hits[which(!(rownames(HGHEDLV2_sta_bio_hits)%in%rownames(sta_bio_hits_Genetic_cluster))),]
#Genes
HGHEDLV2<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%rownames(HGHEDLV2_sta_bio_hits)),]
HGHEDLV2<-HGHEDLV2[!duplicated(HGHEDLV2),]
HGHEDLV2<-HGHEDLV2[!duplicated(HGHEDLV2[,c(1,4)]),]#remove duplicate CpG to gene associations

HGHEDLV2_genes<-Format_gene_table(HGHEDLV2)

load("deinc2dep_PAWS_hits_Mval.RData")
stat_hits<-as.data.frame(PAWS_Beta)[which(Multi_test_corr_relaxed<=fdr),] 
bio_hits<-PAWS_Beta[which(abs(delbeta)>=dB),] 
sta_bio_hits<-stat_hits[which(rownames(stat_hits)%in%rownames(bio_hits)),]
deinc2dep_sta_bio_hits<-sta_bio_hits[which(!(rownames(sta_bio_hits)%in%rownames(SnpatCpG))),]
deinc2dep_sta_bio_hits<-deinc2dep_sta_bio_hits[which(!(rownames(deinc2dep_sta_bio_hits)%in%rownames(sta_bio_hits_Genetic_cluster))),]
#Genes
deinc2dep<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%rownames(deinc2dep_sta_bio_hits)),]
deinc2dep<-deinc2dep[!duplicated(deinc2dep),]
deinc2dep<-deinc2dep[!duplicated(deinc2dep[,c(1,4)]),]#remove duplicate CpG to gene associations

deinc2dep_genes<-Format_gene_table(deinc2dep)


##################
## CpG PLot
##################

## FUNCTIONS: define em
CpG_plot_inc<-function(CpG_list){
  CpGs<-CpG_list
  plot_beta<-as.data.frame(PAWS_Beta[which(rownames(PAWS_Beta)%in%CpGs),])
  plot_beta$CpG<-rownames(plot_beta)
  plot<-melt(plot_beta, id="CpG")
  plot_CpG<-merge(plot, meta, by.x="variable", by.y="Meth_ID")
  plot_CpG<-merge(plot_CpG, deinc2dep_genes, by.x="CpG",by.y="Probe_ID")
  plot_CpG$Label<-paste(plot_CpG$gene," (", plot_CpG$CpG,")", sep="")
  
  ggplot(plot_CpG, aes(deinc2dep, value))+
    geom_point(aes(color=as.factor(Genetic_cluster)),shape=19)+theme_bw()+
    stat_smooth(method = "lm", size = 1)+
    facet_wrap(~Label)+ylab("Beta Value")+ylim(0,1)+
    scale_color_manual(values=colorRampPalette(brewer.pal(11,"RdYlGn")[1:8])(4), name="Genetic Cluster")+
    ylab("Methylation Level (Beta Value)")+xlab("Income Per Dependent")}

CpG_plot_adv<-function(CpG_list){
  CpGs<-CpG_list
  plot_beta<-as.data.frame(PAWS_Beta[which(rownames(PAWS_Beta)%in%CpGs),])
  plot_beta$CpG<-rownames(plot_beta)
  plot<-melt(plot_beta, id="CpG")
  plot_CpG<-merge(plot, meta, by.x="variable", by.y="Meth_ID")
  plot_CpG<-merge(plot_CpG, FADV_6_genes, by.x="CpG",by.y="Probe_ID")
  plot_CpG$Label<-paste(plot_CpG$gene," (", plot_CpG$CpG,")", sep="")
  
  ggplot(plot_CpG, aes(FADV_6, value))+
    geom_point(aes(color=as.factor(Genetic_cluster)),shape=19)+theme_bw()+
    stat_smooth(method = "lm", size = 1)+
    facet_wrap(~Label)+ylab("Beta Value")+ylim(0,1)+
    scale_color_manual(values=colorRampPalette(brewer.pal(11,"RdYlGn")[1:8])(4), name="Genetic Cluster")+
    ylab("Methylation Level (Beta Value)")+xlab("Fall Adversity Composite")}

CpG_plot_edu<-function(CpG_list){
  CpGs<-CpG_list
  plot_beta<-as.data.frame(PAWS_Beta[which(rownames(PAWS_Beta)%in%CpGs),])
  plot_beta$CpG<-rownames(plot_beta)
  plot<-melt(plot_beta, id="CpG")
  plot_CpG<-merge(plot, meta, by.x="variable", by.y="Meth_ID")
  plot_CpG<-merge(plot_CpG, HGHEDLV2_genes, by.x="CpG",by.y="Probe_ID")
  plot_CpG$Label<-paste(plot_CpG$gene," (", plot_CpG$CpG,")", sep="")
  
  ggplot(plot_CpG, aes(HGHEDLV2, value))+
    geom_point(aes(color=as.factor(Genetic_cluster)),shape=19)+theme_bw()+
    stat_smooth(method = "lm", size = 1)+
    facet_wrap(~Label)+ylab("Beta Value")+ylim(0,1)+
    scale_color_manual(values=colorRampPalette(brewer.pal(11,"RdYlGn")[1:8])(4), name="Genetic Cluster")+
    ylab("Methylation Level (Beta Value)")+xlab("Highest Level Family Education")}



############### 
## Pyro site selection
############### 

CpGsadv<-FADV_6_genes[which(abs(FADV_6_genes$db)>0.12 & FADV_6_genes$corr_pval<0.25),]
CpGsedu<-HGHEDLV2_genes[which(abs(HGHEDLV2_genes$db)>0.15 & HGHEDLV2_genes$corr_pval<0.15),]
CpGsinc<-deinc2dep_genes[which(abs(deinc2dep_genes$db)>0.15 & deinc2dep_genes$corr_pval<0.25),]

CpG_plot_adv(CpGsadv$Probe_ID)
CpG_plot_edu(CpGsedu$Probe_ID)
CpG_plot_inc(CpGsinc$Probe_ID)


pyro<-c("cg10581375","cg14245471","cg26479374")# db 0.13, 0.21, 19

cg02872767 # SNP in region
cg26479374 # Biased pyro
cg17369694 # SNP in region

cg21502834## as alternative (good alternative) income
cg10581375# working great for adversity

cg26479374# needs replacement pyro bad education

cg26511075# SNP kinda close
cg09549987# weird distribution
cg09885502# SNP right beside
cg20944157## seems good

pyro<-c("cg26511075","cg21502834","cg10581375")# db 0.13, 0.17, 0.16


##################################################################################################################################### 
# Confirm Pyro 1

library(hydroGOF)

cg10581375_all<-read.csv("PAWS_cg10581375_Assay 1.csv")
cg10581375_all<-cg10581375_all[which(cg10581375_all$Quality=="Passed"),]
cg10581375<-cg10581375_all[which(cg10581375_all$CpG=="cg10581375"),]

cg10581375_450K<-PAWS_Beta[which(rownames(PAWS_Beta)=="cg10581375"),]
cg10581375_450K<-data.frame(Sample=meta$PAWSG_ID, Methylation450=unlist(cg10581375_450K))
cg10581375_450K_pyro<-merge(cg10581375_450K, cg10581375, by.x="Sample", by.y="Sample.ID")
cor(cg10581375_450K_pyro$Methylation450, cg10581375_450K_pyro$Methylation, method="spearman", use="complete.obs")

rmse(cg10581375_450K_pyro$Methylation450, (cg10581375_450K_pyro$Methylation)/100)

##################################################################################################################################### 
# Confirm Pyro 4

cg21502834_all<-read.csv("PAWS_cg21502834_Assay 4.csv")
cg21502834_all<-cg21502834_all[which(cg21502834_all$Quality=="Passed"),]
cg21502834<-cg21502834_all[which(cg21502834_all$CpG=="cg21502834"),]

cg21502834_450K<-PAWS_Beta[which(rownames(PAWS_Beta)=="cg21502834"),]
cg21502834_450K<-data.frame(Sample=meta$PAWSG_ID, Methylation450=unlist(cg21502834_450K))
cg21502834_450K_pyro<-merge(cg21502834_450K, cg21502834, by.x="Sample", by.y="Sample.ID")
cor(cg21502834_450K_pyro$Methylation450, cg21502834_450K_pyro$Methylation,method="spearman", use="complete.obs")

rmse(cg21502834_450K_pyro$Methylation450, (cg21502834_450K_pyro$Methylation)/100)


##################################################################################################################################### 
# Confirm Pyro 6

cg26511075_all<-read.csv("PAWS_cg26511075_Assay 6.csv")
cg26511075_all<-cg26511075_all[which(cg26511075_all$Quality=="Passed"),]
cg26511075<-cg26511075_all[which(cg26511075_all$CpG=="cg26511075"),]

cg26511075_450K<-PAWS_Beta[which(rownames(PAWS_Beta)=="cg26511075"),]
cg26511075_450K<-data.frame(Sample=meta$PAWSG_ID, Methylation450=unlist(cg26511075_450K))
cg26511075_450K_pyro<-merge(cg26511075_450K, cg26511075, by.x="Sample", by.y="Sample.ID")
cor(cg26511075_450K_pyro$Methylation450, cg26511075_450K_pyro$Methylation, method="spearman", use="complete.obs")

rmse(cg26511075_450K_pyro$Methylation450, (cg26511075_450K_pyro$Methylation)/100)


####### Plot all CpGs
arrya_pyro_plot<-rbind(cg21502834_450K_pyro, cg10581375_450K_pyro, cg26511075_450K_pyro)

arrya_pyro_plot<-merge(arrya_pyro_plot, meta, by.x="Sample", by.y="PAWSG_ID")

arrya_pyro_plot$Methylation<-(arrya_pyro_plot$Methylation/100)
arrya_pyro_plot$Mean<-sapply(1:nrow(arrya_pyro_plot), function(x) mean(arrya_pyro_plot$Methylation450[x], arrya_pyro_plot$Methylation[x]))
arrya_pyro_plot$difference<-sapply(1:nrow(arrya_pyro_plot), function(x) arrya_pyro_plot$Methylation450[x]-arrya_pyro_plot$Methylation[x])
mn<-tapply(arrya_pyro_plot$difference, as.character(arrya_pyro_plot$CpG), mean)
SD<-tapply(arrya_pyro_plot$difference, as.character(arrya_pyro_plot$CpG), sd)
sd_data<-data.frame(CpG=names(mn),
                    mn=mn,
                    sd=SD)

##### adjacent sites show same trends?
cg26511075_450K_pyroall<-merge(cg26511075_450K, cg26511075_all, by.x="Sample", by.y="Sample.ID")#edu
cg26511075_450K_pyroall<-merge(cg26511075_450K_pyroall, meta, by.x="Sample", by.y="PAWSG_ID")

cg21502834_450K_pyroall<-merge(cg21502834_450K, cg21502834_all, by.x="Sample", by.y="Sample.ID")#inc
cg21502834_450K_pyroall<-merge(cg21502834_450K_pyroall, meta, by.x="Sample", by.y="PAWSG_ID")

cg10581375_450K_pyroall<-merge(cg10581375_450K, cg10581375_all, by.x="Sample", by.y="Sample.ID")#adv
cg10581375_450K_pyroall<-merge(cg10581375_450K_pyroall, meta, by.x="Sample", by.y="PAWSG_ID")



## Correlation and bland altman
cor_plot<-ggplot(arrya_pyro_plot)+
  geom_line(aes(Methylation450, Methylation450),arrya_pyro_plot, color="black")+
  geom_point(aes(Methylation450, Methylation),arrya_pyro_plot, color="cornflowerblue",shape=19)+facet_wrap(~CpG, ncol=1)+theme_bw()+
  xlab("Methylation Level (450K)")+  ylab("Methylation Level (Pyro)")+
  theme(axis.text = element_text(size =12, color="black"),
        axis.title = element_text(size =15),
        legend.text = element_text(size =14),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12))

bland_altman<-ggplot(arrya_pyro_plot, aes(Mean, difference))+
  facet_wrap(~CpG, ncol=1)+theme_bw()+
  geom_hline(aes(yintercept=(mn+(sd*2))), sd_data,linetype="dashed")+
  geom_hline(aes(yintercept=(mn-(sd*2))), sd_data, linetype="dashed")+
  geom_hline(aes(yintercept=(mn)), sd_data)+
  geom_point(shape=19, color="cornflowerblue")+
  xlab("Mean Methylation")+ylab("Methylation Difference")+
  theme(axis.text = element_text(size =12, color="black"),
        axis.title = element_text(size =15),
        legend.text = element_text(size =14),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12))

grid.arrange(cor_plot, bland_altman, ncol=2)





## plot adjacent sites

levels(cg26511075_450K_pyroall$CpG)<-c("cg26511075","cg26511075 Adjacent")
edu<-ggplot(cg26511075_450K_pyroall, aes(HGHEDLV2,Methylation))+
  geom_point(aes(color=as.factor(Genetic_cluster)),shape=19)+facet_wrap(~CpG)+theme_bw()+
  stat_smooth(method="lm")+xlab("Parental Education")+
  scale_color_manual(values=colorRampPalette(brewer.pal(11,"RdYlGn")[1:8])(4), name="Genetic Cluster")+ylim(0,100)+
  theme(axis.text = element_text(size =12, color="black"),
        axis.title = element_text(size =15),
        legend.text = element_text(size =14),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12))

levels(cg21502834_450K_pyroall$CpG)<-c("cg21502834","cg21502834 Adjacent")
inc<-ggplot(cg21502834_450K_pyroall, aes(deinc2dep,Methylation))+
  geom_point(aes(color=as.factor(Genetic_cluster)),shape=19)+facet_wrap(~CpG)+theme_bw()+
  stat_smooth(method="lm")+xlab("Income-per-dependent")+
  scale_color_manual(values=colorRampPalette(brewer.pal(11,"RdYlGn")[1:8])(4), name="Genetic Cluster")+ylim(0,100)+
  theme(axis.text = element_text(size =12, color="black"),
        axis.title = element_text(size =15),
        legend.text = element_text(size =14),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12))

levels(cg10581375_450K_pyroall$CpG)<-c("cg10581375","cg10581375 Adjacent")
adv<-ggplot(cg10581375_450K_pyroall, aes(FADV_6,Methylation))+
  geom_point(aes(color=as.factor(Genetic_cluster)),shape=19)+facet_wrap(~CpG)+theme_bw()+
  stat_smooth(method="lm")+xlab("Family Adversity")+
  scale_color_manual(values=colorRampPalette(brewer.pal(11,"RdYlGn")[1:8])(4), name="Genetic Cluster")+ylim(0,100)+
  theme(axis.text = element_text(size =12, color="black"),
        axis.title = element_text(size =15),
        legend.text = element_text(size =14),
        legend.title = element_text(size =12),
        strip.text.x = element_text(size = 12))

grid.arrange(edu, inc, adv, ncol=1)

