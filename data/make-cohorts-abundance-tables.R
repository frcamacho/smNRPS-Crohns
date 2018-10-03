require(tidyverse)
source("format_data_for_lefse.R")
prepare_data<-function(df, metadata_df){
  data_wide <- dcast(df, bgcName ~ Sample, value.var="RPKM") 
  data_wide[is.na(data_wide)] <- 0 
  data_melt <- melt(data_wide) %>%  mutate_if(is.factor, as.character) %>% inner_join(.,metadata_df, by = c("variable"="Sample") )
  names(data_melt)[2]<-"Sample"
  names(data_melt)[3]<-"RPKM"
  return(data_melt)
  
}
prepare_data_hmp2<-function(df, metadata_df){
  data_wide <- dcast(df, bgcName ~ Sample, value.var="RPKM") 
  data_wide[is.na(data_wide)] <- 0 
  data_melt <- melt(data_wide) %>%  mutate_if(is.factor, as.character) %>% inner_join(.,metadata_df, by = c("variable"="ExternalID") )
  names(data_melt)[2]<-"Sample"
  names(data_melt)[3]<-"RPKM"
  return(data_melt)
  
}
#############################################################
#load quantifier data 
metahit_quantifier_df <- read_delim("combined-MetaHIT_Gut-spanish-quantifier-results.txt", col_names = F, delim = "\t")
names(metahit_quantifier_df) <- c("bgcName", "Coverage", "RPKM", "Sample") # add column names

filter_metahit_quantifier_df <- metahit_quantifier_df %>% filter(Coverage >= 50)  # filter by coverage 
metahit_sample_metadata <- read_csv("sample_metadata.csv", col_names = T) # load metadata

metahit_results<- prepare_data(filter_metahit_quantifier_df, metahit_sample_metadata) 
metahit_results_raw_table<-format_data_for_lefse(metahit_results) #filtered at 50% coverage 

write_csv(metahit_results_raw_table, "MetaHIT-Gut_abundance-table.csv",col_names = T)

##########################################################
hmp2_metagenomes_quantifier_df <- read_delim("HMP2-IBD/combined-HMP2-IBD-metagenomes-quantifier-results.txt", col_names = F, delim = "\t")
names(hmp2_metagenomes_quantifier_df) <- c("bgcName", "Coverage", "RPKM", "Sample") # add column names

filter_hmp2_metagenomes_quantifier_df <- hmp2_metagenomes_quantifier_df %>% filter(Coverage >= 50)  # filter by coverage

hmp2_sample_metadata<-read_csv("hmp2_project_metadata_2016-10-15.csv", col_names = TRUE) %>% select(., `External ID`, `Participant ID`,`data_type`, diagnosis )
colnames(hmp2_sample_metadata)<- c("ExternalID", "ParticipantID", "data_type","GroupAttribute")
hmp2_sample_metadata_metagenomes <- hmp2_sample_metadata %>% filter(data_type == "metagenomics") 

hmp2_results<- prepare_data_hmp2(filter_hmp2_metagenomes_quantifier_df, hmp2_sample_metadata_metagenomes) %>% filter(GroupAttribute != "Other")
hmp2_results$GroupAttribute[hmp2_results$GroupAttribute == "Crohn's Disease"]<-"CD"
hmp2_results$GroupAttribute[hmp2_results$GroupAttribute == "Ulcerative colitis"]<-"UC"
hmp2_results$GroupAttribute[hmp2_results$GroupAttribute == "Healthy control"]<-"C"

hmp2_results_raw_table<-format_data_for_lefse(hmp2_results)
write_csv(hmp2_results_raw_table, "HMP2-IBD_abundance-table.csv",col_names = T)

