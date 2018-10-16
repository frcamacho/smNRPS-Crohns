library(tidyverse)
library(ggsci)
library(ggpubr)
library(reshape2)

setwd("/Users/francinecamacho/Documents/Donia_analysis/smNRPS-manuscript/smNRPS-Crohns/data/HMP2-IBD/")

sample_metadata<-read_csv("hmp2_project_metadata_2016-10-15.csv", col_names = TRUE) %>% select(., `External ID`, `Participant ID`,`data_type`, `visit_num`, sex, hispa, race,specify_race, diagnosis )
colnames(sample_metadata)<- c("ExternalID", "ParticipantID", "data_type","visit_num","sex","hispa","race","specify_race", "GroupAttribute")
################### 
#split metadata by omic's data 
DNA.metadata<-sample_metadata %>% filter(data_type == "metagenomics")
RNA.metadata<-sample_metadata %>% filter(data_type == "metatranscriptomics")
#################
#how many samples/participants are matched for DNA and RNA analyses?
matched_samples <-DNA.metadata %>% inner_join(.,RNA.metadata, by=c("ExternalID", "ParticipantID","visit_num","sex","hispa","race","specify_race", "GroupAttribute" ))
paste0( "Total number of Samples:", unique(matched_samples$ExternalID) %>% length())
paste0( "Total number of Participants:", unique(matched_samples$ParticipantID) %>% length())


rainin.exact_results_DNA <- read_delim("smNRPS-exact-unique-reads-HMP2-IBD/combined-HMP2-IBD-metagenomes-rainin-exact-unique-reads-quantifier-results.txt", col_names = F, delim = "\t")
names(rainin.exact_results_DNA) <- c("bgcName", "Coverage", "RPKM", "Sample") # add column names
rainin.exact_results_DNA$data_type<-"metagenomics"
# no cutoff for DNA we have an assumption that if it is detected in the RNA than it will be detected in the DNA 


rainin.exact_results_RNA <- read_delim("smNRPS-exact-unique-reads-HMP2-IBD/combined-HMP2-IBD-metatranscriptomes-rainin-exact-unique-reads-quantifier-results.txt", col_names = F, delim = "\t")
names(rainin.exact_results_RNA) <- c("bgcName", "Coverage", "RPKM", "Sample") # add column names

filter_rainin.exact_results_RNA <- rainin.exact_results_RNA %>% filter(Coverage >= 10)  # filter by coverage 
filter_rainin.exact_results_RNA$data_type<-"metatranscriptomics"



rainin.exact_matched <- rbind(rainin.exact_results_DNA,filter_rainin.exact_results_RNA ) %>% filter(Sample %in% matched_samples$ExternalID )%>% inner_join(.,matched_samples%>% select(-c(3,10)), by = c("Sample"= "ExternalID"))
############################################

prepare_rainin_data <- function(df, metadata_df){
  missingsamples<- metadata_df %>% filter(!ExternalID %in% df$Sample)
  df_wide <- dcast(df, bgcName ~ Sample, value.var="RPKM") 
  i<-1
  for (i in 1:nrow(missingsamples)){
    sample_col <- missingsamples[i,]$ExternalID
    df_wide[, sample_col]<-NA
    i<- i+1 
    
  }
  df_wide[is.na(df_wide)] <- 0 
  df_melt <- melt(df_wide) %>%  mutate_if(is.factor, as.character) %>% inner_join(.,metadata_df, by = c("variable"="ExternalID") )
  names(df_melt)[2]<-"Sample"
  names(df_melt)[3]<-"RPKM"
  return(df_melt)
  
}


calculate_hitRatio <-function(data, metadata){
  metadata_counts<- metadata %>% group_by(GroupAttribute) %>% count() %>% ungroup()
  bgc_counts <- data %>% group_by(bgcName,GroupAttribute) %>% count() %>% ungroup()
  hitRatio_df <- data.frame() # intialize results df 
  for (i in 1:nrow(metadata_counts)){
    diseaseStatus <- metadata_counts[i,]$GroupAttribute
    diseaseStatus_count <- metadata_counts[i,]$n
    disease_hitratio <- bgc_counts %>% filter(GroupAttribute == diseaseStatus) %>%  mutate(., hitRatio = case_when(GroupAttribute == diseaseStatus ~ (n/diseaseStatus_count)*100))
    hitRatio_df<-rbind(hitRatio_df,disease_hitratio)
    
  }
  return(hitRatio_df)
}
calculate_hitRatio_byPart <-function(data, metadata){
  metadata_counts<- metadata %>% group_by(ParticipantID,GroupAttribute) %>% count() %>% ungroup() %>% group_by(GroupAttribute) %>% count() %>% ungroup()
  bgc_counts <- data %>% group_by(bgcName,GroupAttribute,ParticipantID) %>% count() %>% ungroup() %>% group_by(bgcName,GroupAttribute) %>% count() %>% ungroup()
  hitRatio_df <- data.frame() # intialize results df 
  for (i in 1:nrow(metadata_counts)){
    diseaseStatus <- metadata_counts[i,]$GroupAttribute
    diseaseStatus_count <- metadata_counts[i,]$nn
    disease_hitratio <- bgc_counts %>% filter(GroupAttribute == diseaseStatus) %>%  mutate(., hitRatio = case_when(GroupAttribute == diseaseStatus ~ (nn/diseaseStatus_count)*100))
    hitRatio_df<-rbind(hitRatio_df,disease_hitratio)
    
  }
  return(hitRatio_df)
}


rainin.exact_results_DNA_filter<- rainin.exact_results_DNA %>% filter(Coverage >=50)%>% inner_join(.,DNA.metadata%>% select(-c(3,10)), by = c("Sample"= "ExternalID"))


rainin_df_hitratio <- calculate_hitRatio(rainin.exact_results_DNA_filter, DNA.metadata)
rainin_df_hitratiobyPart <- calculate_hitRatio_byPart(rainin.exact_results_DNA_filter, DNA.metadata)
rainin.exact_results_DNA_filter_plots<-prepare_rainin_data(rainin.exact_results_DNA_filter,DNA.metadata)
rainin.exact_results_DNA_filter_plots$GroupAttribute <-factor(rainin.exact_results_DNA_filter_plots$GroupAttribute, levels=c("Healthy control","Crohn's Disease","Ulcerative colitis"))

rainin.exact_results_DNA_filter_plots %>% filter(bgcName !="smNRPS-BP") %>%filter(GroupAttribute!="Other") %>%group_by(GroupAttribute,bgcName) %>% mutate(ymean= mean(RPKM)) %>% mutate(ymedian= median(RPKM)) %>%  ggplot(., aes(Sample,RPKM)) + geom_col()+geom_hline(aes(yintercept=ymedian), color="red", linetype="dashed")+geom_hline(aes(yintercept=ymean), color="red", linetype="solid")+theme_pubclean()+ facet_grid(bgcName~GroupAttribute, switch= "x",scales = "free_x", space = "free_x") + xlab("Disease")+  theme(
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank()) + theme(plot.title = element_text(hjust = 0.5,face = "bold",size = "24"))

rainin.exact_results_DNA_filter_plots %>% filter(bgcName !="smNRPS-BP") %>% filter(GroupAttribute != "Other") %>%  group_by(GroupAttribute) %>% mutate(ymean= mean(RPKM)) %>% mutate(ymedian= median(RPKM)) %>% ggplot(., aes(Sample,RPKM)) + geom_bar(stat = "identity", aes(fill=bgcName))+ theme_pubclean()+geom_hline(aes(yintercept=ymedian), color="black", linetype="dashed")+geom_hline(aes(yintercept=ymean), color="black", linetype="solid")+ facet_grid(~GroupAttribute, switch= "x",scales = "free_x", space = "free_x") + xlab("Disease")+  theme(
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank()) + theme(plot.title = element_text(hjust = 0.5,face = "bold",size = "24")) + scale_fill_npg()
ggsave("combined-smNRPS-exact-unique-reads-HMP2-barplot.eps", width = 11, height = 5, units = "in", device = "eps") 

rainin.exact_results_DNA_filter_plots_new<- rainin.exact_results_DNA_filter_plots%>% filter(bgcName !="smNRPS-BP") %>% filter(GroupAttribute != "Other") 
rainin.exact_results_DNA_filter_plots_new %>%  group_by(GroupAttribute) %>% mutate(ymean= mean(RPKM)) %>% mutate(ymedian= median(RPKM)) %>% ggplot(., aes(Sample,RPKM)) + 
  geom_bar(stat = "identity", aes(fill=bgcName))+ theme_pubclean()+geom_hline(aes(yintercept=ymedian), color="black", linetype="dashed")+geom_hline(aes(yintercept=ymean), color="black", linetype="solid")+ 
  facet_grid(~GroupAttribute, switch= "x",scales = "free_x", space = "free_x") + xlab("Disease")+  theme(
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank()) + coord_cartesian(ylim=c(0, 20))+
  theme(plot.title = element_text(hjust = 0.5,face = "bold",size = "24")) + scale_fill_npg()
ggsave("combined-smNRPS-exact-unique-reads-HMP2-barplot-20-max.eps", width = 11, height = 5, units = "in", device = "eps") 

###############################
#RNA analysis with 10% with DNA matched samples 


rainin_df_hitratio_RNA <- calculate_hitRatio(rainin.exact_matched %>% filter(data_type == "metatranscriptomics"), RNA.metadata)
rainin.exact_results_RNA_filter_plots<-prepare_rainin_data(rainin.exact_matched %>% filter(data_type == "metatranscriptomics"), RNA.metadata)

rainin.exact_results_RNA_filter_plots %>% filter(bgcName !="smNRPS-BP") %>%filter(GroupAttribute!="Other") %>%group_by(GroupAttribute,bgcName) %>% mutate(ymean= mean(RPKM)) %>% mutate(ymedian= median(RPKM)) %>%  ggplot(., aes(Sample,RPKM)) + geom_col()+geom_hline(aes(yintercept=ymedian), color="red", linetype="dashed")+geom_hline(aes(yintercept=ymean), color="red", linetype="solid")+theme_pubclean()+ facet_grid(bgcName~GroupAttribute, switch= "x",scales = "free_x", space = "free_x") + xlab("Disease")+  theme(
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank()) + theme(plot.title = element_text(hjust = 0.5,face = "bold",size = "24"))

