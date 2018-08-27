# PAWS_SES
Analysis is DNAm associations with SES in the PAWS cohort

DNAm data is located at GSE94734


| Script                                 | Action                            | Input | Output |
|----------------------------------------|-----------------------------------|-------|--------|
| **SNP_final_report_Loading_annotated.rmd** |Pulls the files from genome studio;  Runs quality control to remove bad samples; Runs quality controls to remove bad probes;  Does and initial PCA to pull out genetic clusters (but PLINK is the more comprehensive method); Builds ped and map files for PLINK analysis|    PAWS_psychchip_DNAReport.csv PAWS_psychchip_FinalReport_B_allele_freq.txt   |     PAWS_snps_GCScore.RData Quality_control_removed_SNPs.RData   **PAWS_snps_filtered.RData PAWS_snps_Ballelefreq.RData**|
|   **Normalization.Rmd**                        |    Pull files from genome studio; Match PsychChip and 450K based on 15 SNPs;  Probe Filter (forgot polymorphic CpGs so did it later); BMIQ normalization         |   PAWS.all-alldata.txt PAWS-GdatasetforKoborlabJuly2014-ID_ordered.csv   PAWS.all-samplefile.txt |    **PAWS_sampFiltered_prbFilter_BMIQ.RData**  PAWS_sample_info_with_replicates.csv  |
| **PLINK_Data_Prep_annotated.Rmd**          | Make ped and map files for PLINK  |PAWS_snps_filtered.RData PAWS-GdatasetforKoborlabJuly2014-ID_ordered.csv|**PAWS_snps_filtered_MAP.txt PAWS_snps_filtered_ped.txt**|
|**GEO_Data_blood_w_age.Rmd**                |Pulls blood form GEO with similar age range|GEO|GEO_blood_age_samples.RData **GEO_PAWS_aged_matched_whole_blood.RData** GSE53191.RData GSE50222.RData|
|**Post_Normalization_QC.Rmd**|Sorts our the mets data; T-cell and blood contamination testing; Heat scree plots; Meta data correlation; covariate testing for model building; ComBat; Replicate RMSE|PAWS_Beta_norm_filtered.RData PAWS_sample_info_with_replicates.csv PAWS-GdatasetforKoborlabJan2015-.csv 2014-09-03_pawsSampleSheet_all.csv GSE53191.RData GSE50222.RData GEO_PAWS_aged_matched_whole_blood.RData|**PAWS_meta_sentrix_genetic_clusters.RData combat_PAWS_Beta_norep.RData**|
|**Regression_deinc2dep.R**|Runs LME on mvalues for income per dependent|combat_PAWS_Beta_norep.RData PAWS_meta_sentrix_genetic_clusters.RData SNPCpG.RData|deinc2dep_PAWS_hits_Mval.RData deinc2dep_PAWS_hits_Mval.RData|
|**Regression_FADV_6.R**|Runs LME on mvalues for adversity|combat_PAWS_Beta_norep.RData PAWS_meta_sentrix_genetic_clusters.RData SNPCpG.RData|Likelihood_Ratio_Test_lmer_FADV_6_PAWS_Mval_genetic_cluster_and_twin_fixed_allFactors.RData FADV_6_PAWS_hits_Mval.RData|
|**Regression_HGHEDLV2.R**|Runs LME on mvalues for highest education level|combat_PAWS_Beta_norep.RData PAWS_meta_sentrix_genetic_clusters.RData SNPCpG.RData|Likelihood_Ratio_Test_lmer_HGHEDLV2_PAWS_Mval_genetic_cluster_and_twin_fixed_allFactors.RData HGHEDLV2_PAWS_hits_Mval.RData|
|**Genetic_cluster_anova.R**|Collect CpG at which methylation associates with genetic ancestry for removal from main hit lists|combat_PAWS_Beta_norep.RData  PAWS_meta_sentrix_genetic_clusters.RData|**Genetic_cluster_hits.RData**|
|**Overlap_regression interpretation.R**|Generates numbers for hit count tables; Makes CpG-Gene association tables; Tests overlap between income, adversity and education at gene level; Plot individual CpGs; Island Gene body enrichment|combat_PAWS_Beta_norep.RData PAWS_meta_sentrix_genetic_clusters.RData SNPCpG.RData Genetic_cluster_hits.RData **FADV_6_PAWS_hits_Mval.RData HGHEDLV2_PAWS_hits_Mval.RData deinc2dep_PAWS_hits_Mval.RData** Gene_enrichment_functions.R|Gene/CpG tables and figures no major R objects|
|**Pyro_site_selection_and_confirmation.R**|Selectionof CpGs with big delta beta for Pyro|All the things above; csv files from the pyro for each sample|Plots|


# Revisions
| Script                                 | Action                            | Input | Output |
|----------------------------------------|-----------------------------------|-------|--------|
| **sensitivity_analysis_blood_gender.R** |Run original linear models but with coraviates for gender and blood contmination|    combat_PAWS_Beta_norep.RData PAWS_meta_sentrix_genetic_clusters.RData SNPCpG.RData  |     plots and p values   
