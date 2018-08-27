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



## Variable comparison
fdr<-0.2
dB<-0.05


load("FADV_6_PAWS_hits_Mval.RData")
stat_hits<-as.data.frame(PAWS_Beta)[which(Multi_test_corr_relaxed<=fdr),] 
bio_hits<-PAWS_Beta[which(abs(delbeta)>=dB),] 
sta_bio_hits<-stat_hits[which(rownames(stat_hits)%in%rownames(bio_hits)),]
FADV_6_sta_bio_hits<-sta_bio_hits[which(!(rownames(sta_bio_hits)%in%rownames(SnpatCpG))),]
FADV_6_sta_bio_hits<-FADV_6_sta_bio_hits[which(!(rownames(FADV_6_sta_bio_hits)%in%rownames(sta_bio_hits_Genetic_cluster))),]


load("HGHEDLV2_PAWS_hits_Mval.RData")
stat_hits<-as.data.frame(PAWS_Beta)[which(Multi_test_corr_relaxed<=fdr),] 
bio_hits<-PAWS_Beta[which(abs(delbeta)>=dB),] 
sta_bio_hits<-stat_hits[which(rownames(stat_hits)%in%rownames(bio_hits)),]
HGHEDLV2_sta_bio_hits<-sta_bio_hits[which(!(rownames(sta_bio_hits)%in%rownames(SnpatCpG))),]
HGHEDLV2_sta_bio_hits<-HGHEDLV2_sta_bio_hits[which(!(rownames(HGHEDLV2_sta_bio_hits)%in%rownames(sta_bio_hits_Genetic_cluster))),]


load("deinc2dep_PAWS_hits_Mval.RData")
stat_hits<-as.data.frame(PAWS_Beta)[which(Multi_test_corr_relaxed<=fdr),] 
bio_hits<-PAWS_Beta[which(abs(delbeta)>=dB),] 
sta_bio_hits<-stat_hits[which(rownames(stat_hits)%in%rownames(bio_hits)),]
deinc2dep_sta_bio_hits<-sta_bio_hits[which(!(rownames(sta_bio_hits)%in%rownames(SnpatCpG))),]
deinc2dep_sta_bio_hits<-deinc2dep_sta_bio_hits[which(!(rownames(deinc2dep_sta_bio_hits)%in%rownames(sta_bio_hits_Genetic_cluster))),]


## written to hit lists to CpG_hits_published 
## Now updated with stricter genetic cluster removal
      write.csv(rownames(FADV_6_sta_bio_hits), file="CpG_hits_published/HitsFADV_6.csv")
      write.csv(rownames(HGHEDLV2_sta_bio_hits), file="CpG_hits_published/HitsHGHEDLV2.csv")
      write.csv(rownames(deinc2dep_sta_bio_hits), file="CpG_hits_published/Hitsdeinc2dep.csv")


      
      
## FDR Table load the Multi_test_corr_relaxed and delbeta from above for variable of interest
Hit_count<-function(dB, fdr){
  stat_hits<-as.data.frame(PAWS_Beta)[which(Multi_test_corr_relaxed<=fdr),] 
  bio_hits<-PAWS_Beta[which(abs(delbeta)>=dB),] 
  sta_bio_hits<-stat_hits[which(rownames(stat_hits)%in%rownames(bio_hits)),]
  nrow(sta_bio_hits[which(!(rownames(sta_bio_hits)%in%rownames(SnpatCpG))),])# remove SNP CpGs
}
Hit_count(0.05, 0.25)

## FDR Table (with conservative Genetic Cluster correction)
Hit_count_GC<-function(dB, fdr){
  stat_hits<-as.data.frame(PAWS_Beta)[which(Multi_test_corr_relaxed<=fdr),] 
  bio_hits<-PAWS_Beta[which(abs(delbeta)>=dB),] 
  sta_bio_hits<-stat_hits[which(rownames(stat_hits)%in%rownames(bio_hits)),]
  sta_bio_hits<-sta_bio_hits[which(!(rownames(sta_bio_hits)%in%rownames(SnpatCpG))),]
  sta_bio_hits<-sta_bio_hits[which(!(rownames(sta_bio_hits)%in%rownames(sta_bio_hits_Genetic_cluster))),]
  nrow(sta_bio_hits)
}

Hit_count_GC(0.05, 0.2)


## hypo or hyper (with conservative GC filtering)
Hit_count_direction<-function(dB, fdr){
  stat_hits<-as.data.frame(PAWS_Beta)[which(Multi_test_corr_relaxed<=fdr),] 
  bio_hits<-PAWS_Beta[which(delbeta>=dB),] 
  sta_bio_hitshypo<-stat_hits[which(rownames(stat_hits)%in%rownames(bio_hits)),]
  sta_bio_hitshypo<-sta_bio_hitshypo[which(!(rownames(sta_bio_hitshypo)%in%rownames(SnpatCpG))),]
  sta_bio_hitshypo<-sta_bio_hitshypo[which(!(rownames(sta_bio_hitshypo)%in%rownames(sta_bio_hits_Genetic_cluster))),]
  print(paste(nrow(sta_bio_hitshypo), " hits with a delta beta >= ", dB, sep="" ))
  
  bio_hits<-PAWS_Beta[which(delbeta<=-dB),] 
  sta_bio_hitshyper<-stat_hits[which(rownames(stat_hits)%in%rownames(bio_hits)),]
  sta_bio_hitshyper<-sta_bio_hitshyper[which(!(rownames(sta_bio_hitshyper)%in%rownames(SnpatCpG))),]
  sta_bio_hitshyper<-sta_bio_hitshyper[which(!(rownames(sta_bio_hitshyper)%in%rownames(sta_bio_hits_Genetic_cluster))),]
  print(paste(nrow(sta_bio_hitshyper), " hits with a delta beta <= ", -dB, sep="" ))
}

Hit_count_direction(0.05, 0.2)




#Genes assocaited with each variable
load("Gene_CpG_Relations_updatejune2015.RData")

FADV_6<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%rownames(FADV_6_sta_bio_hits)),] # rownames(FADV_6_sta_bio_hits) CpG hit list
FADV_6<-FADV_6[!duplicated(FADV_6),]
FADV_6<-FADV_6[!duplicated(FADV_6[,c(1,4)]),]#remove duplicate CpG to gene associations (ie associated wiht two isoforms of a gene)

HGHEDLV2<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%rownames(HGHEDLV2_sta_bio_hits)),]
HGHEDLV2<-HGHEDLV2[!duplicated(HGHEDLV2),]
HGHEDLV2<-HGHEDLV2[!duplicated(HGHEDLV2[,c(1,4)]),]#remove duplicate CpG to gene associations

deinc2dep<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%rownames(deinc2dep_sta_bio_hits)),]
deinc2dep<-deinc2dep[!duplicated(deinc2dep),]
deinc2dep<-deinc2dep[!duplicated(deinc2dep[,c(1,4)]),]#remove duplicate CpG to gene associations

## written to hit lists to CpG_hits_published 
    write.csv(unique(FADV_6$gene), file="CpG_hits_published/geneHitsFADV_6.csv")
    write.csv(unique(HGHEDLV2$gene), file="CpG_hits_published/geneHitsHGHEDLV2.csv")
    write.csv(unique(deinc2dep$gene), file="CpG_hits_published/geneHitsdeinc2dep.csv")


    
###############
### Is the overlap suprising?
###############
Permutate_overlap<-function(probe_list1, probe_list2){
  len1<-nrow(probe_list1)
  len2<-nrow(probe_list2)
  rnd1<-rownames(PAWS_Beta)[sample(1:nrow(PAWS_Beta), len1)]
  rnd2<-rownames(PAWS_Beta)[sample(1:nrow(PAWS_Beta), len2)]
  
  Gene1<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%rnd1),]
  Gene1<-Gene1[!duplicated(Gene1),]
  Gene1<-Gene1[!duplicated(Gene1[,c(1,4)]),]
  
  Gene2<-Gene_CpG_Relations_update[which(Gene_CpG_Relations_update$Probe_ID%in%rnd2),]
  Gene2<-Gene2[!duplicated(Gene2),]
  Gene2<-Gene2[!duplicated(Gene2[,c(1,4)]),]
  length(intersect(unique(Gene1$gene), unique(Gene2$gene))) 
}


# FADV_6, HGHEDLV2
FADV_6_HGHEDLV2_expected_overlap<-sapply(1:10000, function(seed){
  set.seed(seed)
  Permutate_overlap(FADV_6_sta_bio_hits, HGHEDLV2_sta_bio_hits)
})
length(intersect(unique(FADV_6$gene), unique(HGHEDLV2$gene))) #15
mean(FADV_6_HGHEDLV2_expected_overlap)# 18
sd(FADV_6_HGHEDLV2_expected_overlap)

# FADV_6, deinc2dep
FADV_6_deinc2dep_expected_overlap<-sapply(1:10000, function(seed){
  set.seed(seed)
  Permutate_overlap(FADV_6_sta_bio_hits, deinc2dep_sta_bio_hits)
})
length(intersect(unique(FADV_6$gene), unique(deinc2dep$gene))) #42
mean(FADV_6_deinc2dep_expected_overlap)# 20
sd(FADV_6_deinc2dep_expected_overlap)

#deinc2dep, HGHEDLV2
deinc2dep_HGHEDLV2_expected_overlap<-sapply(1:10000, function(seed){
  set.seed(seed)
  Permutate_overlap(deinc2dep_sta_bio_hits, HGHEDLV2_sta_bio_hits)
})
length(intersect(unique(deinc2dep$gene), unique(HGHEDLV2$gene))) #35
mean(deinc2dep_HGHEDLV2_expected_overlap)#56
sd(deinc2dep_HGHEDLV2_expected_overlap)

barplot<-data.frame(Comparison=c("FADV_6:HGHEDLV2","FADV_6:deinc2dep","deinc2dep:HGHEDLV2",
                        "FADV_6:HGHEDLV2","FADV_6:deinc2dep","deinc2dep:HGHEDLV2"),
           Overlap=c(length(intersect(unique(FADV_6$gene), unique(HGHEDLV2$gene))), 
                     length(intersect(unique(FADV_6$gene), unique(deinc2dep$gene))),
                     length(intersect(unique(deinc2dep$gene), unique(HGHEDLV2$gene))),
                     mean(FADV_6_HGHEDLV2_expected_overlap),
                     mean(FADV_6_deinc2dep_expected_overlap),
                     mean(deinc2dep_HGHEDLV2_expected_overlap)),
           Sd=c(0,0,0,
                sd(FADV_6_HGHEDLV2_expected_overlap), sd(FADV_6_deinc2dep_expected_overlap), sd(deinc2dep_HGHEDLV2_expected_overlap)),
           CpG_list=c("Significantly Associated","Significantly Associated","Significantly Associated",
                      "Random", "Random", "Random"))

levels(barplot$Comparison)<-c("Family Adversity \nParental Education",
                              "Family Adversity \nIncome-per-dependent",
                              "Income-per-dependent \nParental Education")

levels(barplot$CpG_list)<-c("Random Gene Lists \nPermutations", "Genes with \nDifferential DNAm")

ggplot(barplot, aes(Comparison, Overlap, group=CpG_list, fill=CpG_list))+
  geom_bar(stat="identity",position="dodge", color="black")+theme_bw()+
  geom_errorbar(aes(ymin=Overlap-Sd, ymax=Overlap+Sd),width=.1, position=position_dodge(width=0.9))+
  scale_fill_manual(values=c("grey", "cornflowerblue"), name="Gene List")+ylab("Gene Overlap Percent")




## permutation Pvalue
real_overlap<-length(intersect(unique(FADV_6$gene), unique(HGHEDLV2$gene)))
length(which(FADV_6_HGHEDLV2_expected_overlap>real_overlap))/length(FADV_6_HGHEDLV2_expected_overlap)

real_overlap<-length(intersect(unique(deinc2dep$gene), unique(HGHEDLV2$gene))) 
length(which(deinc2dep_HGHEDLV2_expected_overlap>real_overlap))/length(deinc2dep_HGHEDLV2_expected_overlap)

real_overlap<-length(intersect(unique(deinc2dep$gene), unique(FADV_6$gene))) 
length(which(FADV_6_deinc2dep_expected_overlap>=real_overlap))/length(FADV_6_deinc2dep_expected_overlap)
p.adjust(length(which(FADV_6_deinc2dep_expected_overlap<=real_overlap))/length(FADV_6_deinc2dep_expected_overlap), method="fdr", n=3)


Overlap_pvalue_function<-function(genelist1, genelist2, permutated_overlaps, multiple_test_correction_number){
  real_overlap<-length(intersect(unique(genelist1), unique(genelist2))) 
  print(paste("Corrected P value for the question are the lists more overlapping than by chance?",
        p.adjust(length(which(permutated_overlaps>=real_overlap))/length(permutated_overlaps), method="fdr", n=multiple_test_correction_number), sep=" "))
  print(paste("Corrected P value for the question are the lists more distinct than by chance?",
        p.adjust(length(which(permutated_overlaps<=real_overlap))/length(permutated_overlaps), method="fdr", n=multiple_test_correction_number), sep=" "))}

## I had three gene lists to compare
Overlap_pvalue_function(unique(deinc2dep$gene),unique(FADV_6$gene), FADV_6_deinc2dep_expected_overlap,3)
Overlap_pvalue_function(unique(FADV_6$gene),unique(HGHEDLV2$gene), FADV_6_HGHEDLV2_expected_overlap,3)
Overlap_pvalue_function(unique(deinc2dep$gene),unique(HGHEDLV2$gene), deinc2dep_HGHEDLV2_expected_overlap,3)






######### intergenic CpG overlap (as these CpGs were not compared in the gene level overlap)
deinc2dep_inter<-deinc2dep[which(deinc2dep$region=="intergenic"),]
FADV_6_inter<-FADV_6[which(FADV_6$region=="intergenic"),]
HGHEDLV2_inter<-HGHEDLV2[which(HGHEDLV2$region=="intergenic"),]


### Is the overlap suprising?
Permutate_overlap_CpGs<-function(probe_list1, probe_list2){
  len1<-length(probe_list1)
  len2<-length(probe_list2)
  rnd1<-rownames(PAWS_Beta)[sample(1:nrow(PAWS_Beta), len1)]
  rnd2<-rownames(PAWS_Beta)[sample(1:nrow(PAWS_Beta), len2)]
  length(intersect(rnd1, rnd2)) 
}


# FADV_6, HGHEDLV2
FADV_6_HGHEDLV2_inter_expected_overlap<-sapply(1:10000, function(seed){
  set.seed(seed)
  print(seed)
  Permutate_overlap_CpGs(FADV_6_inter, HGHEDLV2_inter)
})

length(intersect(FADV_6_inter$Probe_ID, HGHEDLV2_inter$Probe_ID)) #0
mean(FADV_6_HGHEDLV2_inter_expected_overlap)# 18
sd(FADV_6_HGHEDLV2_inter_expected_overlap)


# deinc2dep, HGHEDLV2
deinc2dep_HGHEDLV2_inter_expected_overlap<-sapply(1:10000, function(seed){
  set.seed(seed)
  print(seed)
  Permutate_overlap_CpGs(deinc2dep_inter, HGHEDLV2_inter)
})

length(intersect(deinc2dep_inter$Probe_ID, HGHEDLV2_inter$Probe_ID)) #0
mean(deinc2dep_HGHEDLV2_inter_expected_overlap)# 18
sd(deinc2dep_HGHEDLV2_inter_expected_overlap)

# deinc2dep, FADV_6
deinc2dep_FADV_6_inter_expected_overlap<-sapply(1:10000, function(seed){
  set.seed(seed)
  print(seed)
  Permutate_overlap_CpGs(deinc2dep_inter, FADV_6_inter)
})

length(intersect(deinc2dep_inter$Probe_ID, FADV_6_inter$Probe_ID)) #0
mean(deinc2dep_FADV_6_inter_expected_overlap)#0
sd(deinc2dep_FADV_6_inter_expected_overlap)



## permutation Pvalue
real_overlap<-length(intersect(FADV_6_inter$Probe_ID, HGHEDLV2_inter$Probe_ID))
length(which(FADV_6_HGHEDLV2_inter_expected_overlap>real_overlap))/length(FADV_6_HGHEDLV2_inter_expected_overlap)

real_overlap<-length(intersect(deinc2dep_inter$Probe_ID, HGHEDLV2_inter$Probe_ID))
length(which(deinc2dep_HGHEDLV2_inter_expected_overlap>real_overlap))/length(deinc2dep_HGHEDLV2_inter_expected_overlap)

real_overlap<-length(intersect(deinc2dep_inter$Probe_ID, FADV_6_inter$Probe_ID))
length(which(deinc2dep_FADV_6_inter_expected_overlap>=real_overlap))/length(deinc2dep_FADV_6_inter_expected_overlap)

p.adjust(length(which(deinc2dep_FADV_6_inter_expected_overlap<=real_overlap))/length(deinc2dep_FADV_6_inter_expected_overlap), method="fdr", n=3)


Overlap_pvalue_function<-function(genelist1, genelist2, permutated_overlaps, multiple_test_correction_number){
  real_overlap<-length(intersect(unique(genelist1), unique(genelist2))) 
  print(paste("Corrected P value for the question are the lists more overlapping than by chance?",
              p.adjust(length(which(permutated_overlaps>=real_overlap))/length(permutated_overlaps), method="fdr", n=multiple_test_correction_number), sep=" "))
  print(paste("Corrected P value for the question are the lists more distinct than by chance?",
              p.adjust(length(which(permutated_overlaps<=real_overlap))/length(permutated_overlaps), method="fdr", n=multiple_test_correction_number), sep=" "))}






###############
## GENE TABLE
###############

## Gene CpG number adjustment
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
write.csv(FADV_6_genes, file="CpG_hits_published/FADV_6_genes_FDR02_DB05.csv")


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
write.csv(HGHEDLV2_genes, file="CpG_hits_published/HGHEDLV2_genes_FDR02_DB05.csv")


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
write.csv(deinc2dep_genes, file="CpG_hits_published/deinc2dep_genes_FDR02_DB05.csv")







################
## CpG PLot
################
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

# Cherry pick some CpGs to plot
CpGsadv<-FADV_6_genes$Probe_ID[which(abs(FADV_6_genes$db)>0.09 & FADV_6_genes$corr_pval<0.2)]
CpG_plot_adv(CpGsadv)

CpGsedu<-HGHEDLV2_genes$Probe_ID[which(abs(HGHEDLV2_genes$db)>0.11 & HGHEDLV2_genes$corr_pval<0.075)]
CpG_plot_edu(CpGsedu)

CpGsinc<-deinc2dep_genes$Probe_ID[which(abs(deinc2dep_genes$db)>0.1 & deinc2dep_genes$corr_pval<0.06325)]
CpG_plot_inc(CpGsinc)




#############
## Enrichment
#############
load("Gene_CpG_Relations_updatejune2015.RData")
load("Price_annotation.RData") #annotations based on Price et al. 2013 Epigenetics & Chromatin
annotation$CpG<-rownames(annotation)

source("Gene_enrichment_functions.R")


# define background of CpGs and the CpG list of interest (hits)
background.PAWS<-rownames(PAWS_Beta) #426893
deinc2dep_CpGs<-rownames(deinc2dep_sta_bio_hits) #488
HGHEDLV2_CpGs<-rownames(HGHEDLV2_sta_bio_hits) #354
FADV_6_CpGs<-rownames(FADV_6_sta_bio_hits) #102


## plot fold enrichment and the permutation pvalues
CGI_Gene_permutation_enrichment(deinc2dep_CpGs,background.PAWS, 1000)
CGI_Gene_permutation_enrichment(HGHEDLV2_CpGs,background.PAWS, 1000)
CGI_Gene_permutation_enrichment(FADV_6_CpGs,background.PAWS, 1000)
0.05
