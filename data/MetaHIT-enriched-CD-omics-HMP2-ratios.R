#################################
#Enriched clusters (CD_C) in MetaHIT
#HMP2-IBD profiler
#################################
library(tidyverse)
library(ggsci)
library(ggpubr)

#########################
#download data 
quantifier_results_DNA <- read_delim("combined-HMP2-IBD-metagenomes-spanish-clusters-quantifier-results.txt", col_names = F, delim = "\t")
names(quantifier_results_DNA) <- c("bgcName", "Coverage", "RPKM", "Sample") # add column names

quantifier_results_RNA <- read_delim("combined-HMP2-IBD-metatranscriptomes-spanish-clusters-quantifier-results.txt", col_names = F, delim = "\t")
names(quantifier_results_RNA) <- c("bgcName", "Coverage", "RPKM", "Sample") # add column names

enriched_metahit_cd <- read_csv("supplementary-tables/lefse_CD_C_supp_table_1.csv",col_names = T) %>% filter(class_discriminitive == "CD" ) #43
################################################################################
# HMP2 IBD metadata 
sample_metadata<-read_csv("hmp2_project_metadata_2016-10-15.csv", col_names = TRUE) %>% select(., `External ID`, `Participant ID`,`data_type`, `visit_num`, sex, hispa, race,specify_race, diagnosis )
colnames(sample_metadata)<- c("ExternalID", "ParticipantID", "data_type","visit_num","sex","hispa","race","specify_race", "GroupAttribute")

################### 
#split metadata by omic's data 
DNA.metadata<-sample_metadata %>% filter(data_type == "metagenomics")
RNA.metadata<-sample_metadata %>% filter(data_type == "metatranscriptomics")

#################
#filter by Coverage 50% in DNA 
DNA_quantifier_df <-quantifier_results_DNA %>% inner_join(., DNA.metadata, by =c("Sample" = "ExternalID")) %>% filter(Coverage >=50)

RNA_quantifier_df <-quantifier_results_RNA %>% inner_join(., RNA.metadata, by =c("Sample" = "ExternalID"))
##
#how many samples/participants are matched for DNA and RNA analyses?
matched_samples <-DNA.metadata %>% inner_join(.,RNA.metadata, by=c("ExternalID", "ParticipantID","visit_num","sex","hispa","race","specify_race", "GroupAttribute" ))
paste0( "Total number of Samples:", unique(matched_samples$ExternalID) %>% length())
paste0( "Total number of Participants:", unique(matched_samples$ParticipantID) %>% length())
quantifier_df_matched <- rbind(DNA_quantifier_df,RNA_quantifier_df ) %>% filter(Sample %in% matched_samples$ExternalID )
#################
#filter data by enriched CD in HMP2 only
quantifier_df_matched_enriched <- quantifier_df_matched %>% filter(bgcName %in% enriched_metahit_cd$BGC)
#plot RNA and DNA by BGC 
quantifier_df_matched_enriched %>% ggplot(., aes(Sample, Coverage,fill=data_type))+  geom_bar(stat = "identity",position = "dodge")+
  facet_wrap(~bgcName)+ theme(axis.text.x = element_text(angle = 90, hjust = .5, size =6)) + geom_hline(aes(yintercept=10))+theme(
  strip.text.x = element_text(size = 6))

missing_bgcs <- enriched_metahit_cd %>% filter(!BGC %in% quantifier_df_matched_enriched$bgcName) # missing 6 BGCs that are not in HMP2 
#############
#RNA /DNA log ratios 
quantifier_df_matched_enriched_wide <- quantifier_df_matched_enriched %>% select(-c(Coverage, Sample)) %>% group_by(ParticipantID, bgcName, data_type) %>% mutate(meanRPKM = mean(RPKM)) %>% 
  select(-c(visit_num, sex, hispa, race,specify_race,RPKM))%>% distinct() %>% spread(., data_type, meanRPKM)


#405rows for DNA is 0s 
#5 rows for RNA is 0s 

####################################
#filter for data that has both DNA and RNA 

enriched_match_dna_rna<-quantifier_df_matched_enriched_wide %>% filter(!is.na(metagenomics) & !is.na(metatranscriptomics)) %>% mutate(logRatio = log2(metatranscriptomics/metagenomics))

enriched_match_dna_rna %>% ggplot(., aes(ParticipantID, logRatio, fill=GroupAttribute))+  geom_col()+facet_wrap(~bgcName)+ theme(plot.title = element_text(hjust = 0.5,face = "bold",size = "24"))+
  theme(axis.text.x = element_text(angle = 90, hjust = .5, size =6)) + scale_fill_npg()


