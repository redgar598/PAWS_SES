library(ggplot2)
library(RColorBrewer)
library(scales)
library(gridExtra)
library(reshape)
library(lme4)


setwd("/big_data/redgar/PAWS/PAWS_final/reviews_response/")
load("/big_data/redgar/PAWS/PAWS_final/combat_PAWS_Beta_norep.RData")
load("/big_data/redgar/PAWS/PAWS_final/PAWS_meta_sentrix_genetic_clusters.RData")
load("/big_data/redgar/PAWS/PAWS_final/Price_annotation.RData")


meta<-meta[which(meta$Meth_ID%in%colnames(PAWS_Beta)),]
meta<-meta[match(colnames(PAWS_Beta), meta$Meth_ID),]

Mval<-function(beta) log2(beta/(1-beta))
PAWS_Mval = apply(PAWS_Beta, 1, Mval) # need mvalues for combat
PAWS_Mval = as.data.frame(PAWS_Mval)
PAWS_Mval = t(PAWS_Mval)



# going to look at p values within 1000bp of hits

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


adj_CpG<-function(cpg){
  coor<-as.numeric(as.character(annotation$MAPINFO[which(annotation$TargetID==cpg)]))
  chr<-as.numeric(as.character(annotation$CHR[which(annotation$TargetID==cpg)]))
  
  chr_match<-annotation[which(annotation$CHR==chr),]
  region<-chr_match[which(chr_match$MAPINFO>(coor-1000) & chr_match$MAPINFO<(coor+1000)),c("TargetID","CHR","MAPINFO")]
  
  adj_pval<-merge(region, paws_data, by.x="TargetID", by.y="CpG")
  
  adj_pval$distance<-coor-adj_pval$MAPINFO
  adj_pval$og_cpg_db<-paws_data$db[which(paws_data$CpG==cpg)]
  adj_pval}
  

adj_pvalues_list<-lapply(1:nrow(paws_data_hits), function(x) adj_CpG(as.character(paws_data_hits$CpG[x])))
sapply(1:nrow(paws_data_hits), function(x) nrow(adj_pvalues_list[[x]]))

##145/488 are singletons

adj_pvalues<-do.call(rbind, adj_pvalues_list)
# 
# ggplot(adj_pvalues, aes(distance, pval))+geom_point()+stat_smooth()
# ggplot(adj_pvalues[which(adj_pvalues$distance!=0),], aes(distance, pval))+geom_point()+stat_smooth()


### DB plot
adj_pvalues_adjonly_inc<-adj_pvalues[which(adj_pvalues$distance!=0),]
cor(adj_pvalues_adjonly_inc$og_cpg_db,adj_pvalues_adjonly_inc$db)

inc<-ggplot(adj_pvalues_adjonly_inc, aes(og_cpg_db,db))+geom_hline(yintercept=0,linetype="longdash", color="lightgrey")+
  geom_vline(xintercept=0, linetype="longdash",  color="lightgrey")+
  geom_point(shape=21, color="black", fill="grey",size=2)+
  theme_bw()+xlab("Delta Beta\nSignificant Income-per-Dependent Associated CpGs")+ylab("Delta Beta\nCpGs Within 1kb of Significant EWAS CpGs")+
  ylim(-0.35, 0.35)+xlim(-0.35,0.35)+stat_smooth(se=F, color="black", method="lm")+
  annotate("text",x=0.25,y=-0.2, label = as.character(expression(italic("r")[s]~"=")), parse = TRUE)+
  annotate("text",x=0.3,y=-0.2, label = round(cor(adj_pvalues_adjonly_inc$og_cpg_db,adj_pvalues_adjonly_inc$db, method="spearman"), 2))



  
# 
# ### rep DMR
# # tapply(adj_pvalues_adjonly$pval, adj_pvalues_adjonly$og_cpg_db, mean)[order(tapply(adj_pvalues_adjonly$pval, adj_pvalues_adjonly$og_cpg_db, mean))]
# # table(adj_pvalues_adjonly$og_cpg_db)[order(names(table(adj_pvalues_adjonly$og_cpg_db)))]
# # adj_pvalues[which(adj_pvalues$og_cpg_db>0.05819 & adj_pvalues$og_cpg_db<0.0582),]
# # 
# # 
# # rep_DMR<-adj_CpG("cg01760090") # intergenic
# # rep_DMR<-adj_CpG("cg05041795") # intergenic
# # rep_DMR<-adj_CpG("cg10456990") # SIM2
# 
# rep_DMR<-adj_CpG("cg01881182") # ABAT
# # rep_DMR<-adj_CpG("cg08202494") # OSR2
# # 
# # rep_DMR<-adj_CpG("cg09577511") # AQP2
# 
# 
# 
# 
# ### rep CpG list for gene list
# 
# # adv_repcpg<-c("cg06588556","cg21676617","cg26053697","cg10581375")
# # edu_rep_cpg<-c("cg02918577", "cg14245471", "cg20944157", "cg26511075")
# inc_repcpg<-rep_DMR$TargetID
# 
# 
# ## Gene CpG number adjustment
# load("/big_data/redgar/PAWS/PAWS_final/Gene_CpG_Relations_updatejune2015.RData")
# Gene_CpG_Relations_update$gene<-as.character(Gene_CpG_Relations_update$gene)
# Overrep<-as.data.frame(tapply(Gene_CpG_Relations_update$Probe_ID, Gene_CpG_Relations_update$gene, length))
# Overrep$Gene<-rownames(Overrep)
# colnames(Overrep)<-c("CpG_number", "Gene")
# Overrep<-Overrep[which(Overrep$Gene!="None"),]
# mean(Overrep$CpG_number, na.rm=T)# 25
# Overrep$Enrichment_fromAverage<-Overrep$CpG_number/mean(Overrep$CpG_number, na.rm=T)
# 
# ## Gene summaries
# load("/big_data/redgar/PAWS/PAWS_final/Price_annotation.RData")
# annotation$CpG<-rownames(annotation)
# 
# Format_gene_table<-function(Gene_CpG_Relations_update_subset){
#   print(paste("CpGs Associated: ", length(unique(Gene_CpG_Relations_update_subset$Probe_ID)), sep=""))
#   print(paste("Genes Associated: ", length(unique(Gene_CpG_Relations_update_subset$gene)), sep=""))
#   Overrep_subset<-as.data.frame(tapply(Gene_CpG_Relations_update_subset$Probe_ID, Gene_CpG_Relations_update_subset$gene, length))
#   Overrep_subset$Gene<-rownames(Overrep_subset)
#   colnames(Overrep_subset)<-c("CpG_number", "Gene")
#   Overrep_subset<-Overrep_subset[which(Overrep_subset$Gene!="None"),]
#   Overrep_subset_merge<-merge(Overrep_subset, Overrep, by="Gene")
#   colnames(Overrep_subset_merge)<-c("Gene","CpG_Associated","CpG_in_Gene", "Enrichment_fromAverage")
#   Overrep_subset_merge$Suprise<-Overrep_subset_merge$CpG_Associated/Overrep_subset_merge$Enrichment_fromAverage
#   Gene_table<-merge(Gene_CpG_Relations_update_subset, Overrep_subset_merge, by.x="gene", by.y="Gene")
#   Gene_table<-merge(Gene_table, annotation[,c(49,50,58)], by.x="Probe_ID", by.y="CpG")
#   pval<-data.frame(CpG=rownames(PAWS_Beta)[which(rownames(PAWS_Beta)%in%Gene_table$Probe_ID)],
#                    corr_pval=Multi_test_corr_relaxed[which(rownames(PAWS_Beta)%in%Gene_table$Probe_ID)],
#                    db=delbeta[which(rownames(PAWS_Beta)%in%Gene_table$Probe_ID)])
#   Gene_table<-merge(Gene_table, pval, by.x="Probe_ID", by.y="CpG")
#   Gene_table<-Gene_table[,c(2,7,8,9,10,1,4,3,6,5,11,12,13,14)]
#   Gene_table<-Gene_table[order(-Gene_table$Suprise, Gene_table$gene),]
#   Gene_table}
# 
# 
# 
# #Genes 
# # FADV_6<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%adv_repcpg),]
# # FADV_6<-FADV_6[!duplicated(FADV_6),]
# # FADV_6<-FADV_6[!duplicated(FADV_6[,c(1,4)]),]#remove duplicate CpG to gene associations
# # FADV_6_genes<-Format_gene_table(FADV_6)
# 
# deinc2dep<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%inc_repcpg),]
# deinc2dep<-deinc2dep[!duplicated(deinc2dep),]
# deinc2dep<-deinc2dep[!duplicated(deinc2dep[,c(1,4)]),]#remove duplicate CpG to gene associations
# deinc2dep_genes<-Format_gene_table(deinc2dep)
# 
# # HGHEDLV2<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%edu_rep_cpg),]
# # HGHEDLV2<-HGHEDLV2[!duplicated(HGHEDLV2),]
# # HGHEDLV2<-HGHEDLV2[!duplicated(HGHEDLV2[,c(1,4)]),]#remove duplicate CpG to gene associations
# # HGHEDLV2_genes<-Format_gene_table(HGHEDLV2)
# 
# ## income
# CpGs<-inc_repcpg
# plot_beta<-as.data.frame(PAWS_Beta[which(rownames(PAWS_Beta)%in%CpGs),])
# minbeta<-sapply(1:nrow(plot_beta), function(x) quantile(plot_beta[x,], 0.05, na.rm=T))
# maxbeta<-sapply(1:nrow(plot_beta), function(x) quantile(plot_beta[x,], 0.95, na.rm=T))
# 
# plot_beta$CpG<-rownames(plot_beta)
# plot<-melt(plot_beta, id="CpG")
# 
# plot_CpG<-merge(plot, meta, by.x="variable", by.y="Meth_ID")
# plot_CpG<-merge(plot_CpG, deinc2dep_genes, by.x="CpG",by.y="Probe_ID")
# plot_CpG$Label<-paste(plot_CpG$gene," (", plot_CpG$CpG,")", sep="")
# 
# plot_CpG$Label<-as.factor( plot_CpG$Label)
# plot_CpG$Label = factor(plot_CpG$Label, levels=paste(unique(plot_CpG$gene)," (", as.character(rep_DMR$TargetID[order(rep_DMR$MAPINFO)]),")", sep=""))
# 
# 
# ggplot(plot_CpG, aes(deinc2dep, value))+
#   geom_point(aes(color=as.factor(Genetic_cluster)),shape=19)+theme_bw()+
#   stat_smooth(method = "lm", size = 0.75, se=F, color="grey20")+
#   facet_wrap(~Label, nrow=1)+ylab("Beta Value")+ylim(0,1)+
#   scale_color_manual(values=colorRampPalette(brewer.pal(11,"RdYlGn")[1:8])(4), name="PLINK Genetic Cluster")+
#   ylab("Methylation Level (Beta Value)")+xlab("Income-Per-Dependent")
# 
# 
# # 
# # 
# # CpG_plot_adv<-function(CpG_list){
# #   CpGs<-CpG_list
# #   plot_beta<-as.data.frame(PAWS_Beta[which(rownames(PAWS_Beta)%in%CpGs),])
# #   
# #   minbeta<-sapply(1:nrow(plot_beta), function(x) quantile(plot_beta[x,], 0.05, na.rm=T))
# #   maxbeta<-sapply(1:nrow(plot_beta), function(x) quantile(plot_beta[x,], 0.95, na.rm=T))
# #   
# #   plot_beta_win<-lapply(1:nrow(plot_beta), function(x){
# #     y<-plot_beta[x,]
# #     y[which(y<minbeta[x])]<-minbeta[x]
# #     y[which(y>maxbeta[x])]<-maxbeta[x]
# #     y  })
# #   
# #   plot_beta_win<-do.call(rbind, plot_beta_win)
# #   rownames(plot_beta_win)<-rownames(plot_beta)
# #   
# #   plot_beta_win$CpG<-rownames(plot_beta_win)
# #   plot<-melt(plot_beta_win, id="CpG")
# #   plot_CpG<-merge(plot, meta, by.x="variable", by.y="Meth_ID")
# #   plot_CpG<-merge(plot_CpG, FADV_6_genes, by.x="CpG",by.y="Probe_ID")
# #   plot_CpG$Label<-paste(plot_CpG$gene," (", plot_CpG$CpG,")", sep="")
# #   
# #   ggplot(plot_CpG, aes(FADV_6, value))+
# #     geom_point(aes(color=as.factor(Genetic_cluster)),shape=19)+theme_bw()+
# #     stat_smooth(method = "lm", size = 0.75, se=F, color="grey20")+
# #     facet_wrap(~Label)+ylab("Beta Value")+ylim(0,1)+
# #     scale_color_manual(values=colorRampPalette(brewer.pal(11,"RdYlGn")[1:8])(4), name="PLINK Genetic Cluster")+
# #     ylab("Methylation Level (Beta Value)")+xlab("Family Adversity")}
# # 
# # CpG_plot_edu<-function(CpG_list){
# #   CpGs<-CpG_list
# #   plot_beta<-as.data.frame(PAWS_Beta[which(rownames(PAWS_Beta)%in%CpGs),])
# #   
# #   minbeta<-sapply(1:nrow(plot_beta), function(x) quantile(plot_beta[x,], 0.05, na.rm=T))
# #   maxbeta<-sapply(1:nrow(plot_beta), function(x) quantile(plot_beta[x,], 0.95, na.rm=T))
# #   
# #   plot_beta_win<-lapply(1:nrow(plot_beta), function(x){
# #     y<-plot_beta[x,]
# #     y[which(y<minbeta[x])]<-minbeta[x]
# #     y[which(y>maxbeta[x])]<-maxbeta[x]
# #     y  })
# #   
# #   plot_beta_win<-do.call(rbind, plot_beta_win)
# #   rownames(plot_beta_win)<-rownames(plot_beta)
# #   
# #   
# #   plot_beta_win$CpG<-rownames(plot_beta_win)
# #   plot<-melt(plot_beta_win, id="CpG")
# #   plot_CpG<-merge(plot, meta, by.x="variable", by.y="Meth_ID")
# #   plot_CpG<-merge(plot_CpG, HGHEDLV2_genes, by.x="CpG",by.y="Probe_ID")
# #   plot_CpG$Label<-paste(plot_CpG$gene," (", plot_CpG$CpG,")", sep="")
# #   
# #   ggplot(plot_CpG, aes(HGHEDLV2, value))+
# #     geom_point(aes(color=as.factor(Genetic_cluster)),shape=19)+theme_bw()+
# #     stat_smooth(method = "lm", size = 0.75, se=F, color="grey20")+
# #     facet_wrap(~Label)+ylab("Beta Value")+ylim(0,1)+
# #     scale_color_manual(values=colorRampPalette(brewer.pal(11,"RdYlGn")[1:8])(4), name="PLINK Genetic Cluster")+
# #     ylab("Methylation Level (Beta Value)")+xlab("Parental Education")}
# # 
# # 
# # 
# # adv<-CpG_plot_adv(adv_repcpg)
# # edu<-CpG_plot_edu(edu_rep_cpg)
# # 
# # 
# # 
# # ## FUNCTIONS
# # 
# # 
# 
# 
# 


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



adj_CpG<-function(cpg){
  coor<-as.numeric(as.character(annotation$MAPINFO[which(annotation$TargetID==cpg)]))
  chr<-as.numeric(as.character(annotation$CHR[which(annotation$TargetID==cpg)]))
  
  chr_match<-annotation[which(annotation$CHR==chr),]
  region<-chr_match[which(chr_match$MAPINFO>(coor-1000) & chr_match$MAPINFO<(coor+1000)),c("TargetID","CHR","MAPINFO")]
  
  adj_pval<-merge(region, paws_data, by.x="TargetID", by.y="CpG")
  
  adj_pval$distance<-coor-adj_pval$MAPINFO
  adj_pval$og_cpg_db<-paws_data$db[which(paws_data$CpG==cpg)]
  adj_pval}


adj_pvalues_list<-lapply(1:nrow(paws_data_hits), function(x) adj_CpG(as.character(paws_data_hits$CpG[x])))
sapply(1:nrow(paws_data_hits), function(x) nrow(adj_pvalues_list[[x]]))

adj_pvalues<-do.call(rbind, adj_pvalues_list)
##99/354 are singletons

### DB plot
adj_pvalues_adjonlyedu<-adj_pvalues[which(adj_pvalues$distance!=0),]
cor(adj_pvalues_adjonlyedu$og_cpg_db,adj_pvalues_adjonlyedu$db)

edu<-ggplot(adj_pvalues_adjonlyedu, aes(og_cpg_db,db))+geom_hline(yintercept=0,linetype="longdash", color="lightgrey")+
  geom_vline(xintercept=0, linetype="longdash",  color="lightgrey")+
  geom_point(shape=21, color="black", fill="grey",size=2)+
  theme_bw()+xlab("Delta Beta\nSignificant Parental Education Associated CpGs")+ylab("Delta Beta\nCpGs Within 1kb of Significant EWAS CpGs")+
  ylim(-0.35, 0.35)+xlim(-0.35,0.35)+stat_smooth(se=F, color="black", method="lm")+
  annotate("text",x=0.25,y=-0.2, label = as.character(expression(italic("r")[s]~"=")), parse = TRUE)+
  annotate("text",x=0.3,y=-0.2, label = round(cor(adj_pvalues_adjonlyedu$og_cpg_db,adj_pvalues_adjonlyedu$db, method="spearman"), 2))


















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


adj_CpG<-function(cpg){
  coor<-as.numeric(as.character(annotation$MAPINFO[which(annotation$TargetID==cpg)]))
  chr<-as.numeric(as.character(annotation$CHR[which(annotation$TargetID==cpg)]))
  
  chr_match<-annotation[which(annotation$CHR==chr),]
  region<-chr_match[which(chr_match$MAPINFO>(coor-1000) & chr_match$MAPINFO<(coor+1000)),c("TargetID","CHR","MAPINFO")]
  
  adj_pval<-merge(region, paws_data, by.x="TargetID", by.y="CpG")
  
  adj_pval$distance<-coor-adj_pval$MAPINFO
  adj_pval$og_cpg_db<-paws_data$db[which(paws_data$CpG==cpg)]
  adj_pval}


adj_pvalues_list<-lapply(1:nrow(paws_data_hits), function(x) adj_CpG(as.character(paws_data_hits$CpG[x])))
sapply(1:nrow(paws_data_hits), function(x) nrow(adj_pvalues_list[[x]]))

adj_pvalues<-do.call(rbind, adj_pvalues_list)
##31/102 are singletons

### DB plot
adj_pvalues_adjonlyadv<-adj_pvalues[which(adj_pvalues$distance!=0),]
cor(adj_pvalues_adjonlyadv$og_cpg_db,adj_pvalues_adjonlyadv$db)

adv<-ggplot(adj_pvalues_adjonlyadv, aes(og_cpg_db,db))+geom_hline(yintercept=0,linetype="longdash", color="lightgrey")+
  geom_vline(xintercept=0, linetype="longdash",  color="lightgrey")+
  geom_point(shape=21, color="black", fill="grey",size=2)+
  theme_bw()+xlab("Delta Beta\nSignificant Family Adversity Associated CpGs")+ylab("Delta Beta\nCpGs Within 1kb of Significant EWAS CpGs")+
  ylim(-0.35, 0.35)+xlim(-0.35,0.35)+stat_smooth(se=F, color="black", method="lm")+
  annotate("text",x=0.25,y=-0.2, label = as.character(expression(italic("r")[s]~"=")), parse = TRUE)+
  annotate("text",x=0.3,y=-0.2, label = round(cor(adj_pvalues_adjonlyadv$og_cpg_db,adj_pvalues_adjonlyadv$db, method="spearman"), 2))





grid.arrange(inc, edu, adv)
