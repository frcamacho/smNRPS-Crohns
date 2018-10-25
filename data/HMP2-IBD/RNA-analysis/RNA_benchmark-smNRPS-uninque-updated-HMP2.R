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

########################
rainin.exact_results_DNA <- read_delim("smNRPS-exact-unique-reads-HMP2-IBD-updated/combined-HMP2-IBD-metagenomes-rainin-exact-unique-reads-updated-quantifier-results.txt", col_names = F, delim = "\t")
names(rainin.exact_results_DNA) <- c("bgcName", "Coverage", "RPKM", "Sample") # add column names
rainin.exact_results_DNA$data_type<-"metagenomics"
rainin.exact_results_DNA$Coverage[rainin.exact_results_DNA$Coverage>100] <-100
# add cutoff for median 19% 
# add cutodd for 50% 

rainin_DNA_median_cutoff<-rainin.exact_results_DNA %>% filter(Coverage >=19)
rainin_DNA_50_cutoff<-rainin.exact_results_DNA %>% filter(Coverage >=50)
########################
rainin.exact_results_RNA <- read_delim("smNRPS-exact-unique-reads-HMP2-IBD-updated/combined-HMP2-IBD-metatranscriptomes-rainin-exact-unique-reads-updated-quantifier-results.txt", col_names = F, delim = "\t")
names(rainin.exact_results_RNA) <- c("bgcName", "Coverage", "RPKM", "Sample") # add column names
rainin.exact_results_RNA$data_type<-"metatranscriptomics"

# filter RNA data with samples that passed DNA cutoff and find the matched samples 
rainin.exact_matched_50_cutoff <- rainin.exact_results_RNA %>% filter(Sample%in% rainin_DNA_50_cutoff$Sample) %>% rbind(.,rainin_DNA_50_cutoff) %>% filter(Sample %in% matched_samples$ExternalID )%>%
  inner_join(.,matched_samples %>% select(-c(3,10)), by = c("Sample"= "ExternalID")) %>%select(-c(2,7:11)) %>% spread(., data_type, RPKM) %>% filter(!is.na(metagenomics))

rainin.exact_matched_median_cutoff <- rainin.exact_results_RNA %>% filter(Sample%in% rainin_DNA_median_cutoff$Sample) %>% rbind(.,rainin_DNA_median_cutoff) %>% filter(Sample %in% matched_samples$ExternalID )%>%
  inner_join(.,matched_samples %>% select(-c(3,10)), by = c("Sample"= "ExternalID")) %>%  select(-c(2,7:11)) %>% spread(., data_type, RPKM) %>% filter(!is.na(metagenomics))


rainin.exact_matched_50_cutoff  %>% filter(!is.na(metatranscriptomics)) %>% mutate(logRatio = log2(metatranscriptomics/metagenomics)) %>% ggplot(., aes(GroupAttribute, logRatio)) +
  geom_boxplot() + geom_jitter(width=.1)

rainin.exact_matched_median_cutoff %>% filter(!is.na(metatranscriptomics)) %>% mutate(logRatio = log2(metatranscriptomics/metagenomics)) %>% ggplot(., aes(GroupAttribute, logRatio)) + 
  geom_boxplot()+geom_jitter(width=.1)

#########################
#try 12% 

rainin_DNA_12_cutoff<-rainin.exact_results_DNA %>% filter(Coverage >=12)
rainin.exact_matched_12_cutoff <- rainin.exact_results_RNA %>% filter(Sample%in% rainin_DNA_12_cutoff$Sample) %>% rbind(.,rainin_DNA_12_cutoff) %>% filter(Sample %in% matched_samples$ExternalID )%>%
  inner_join(.,matched_samples %>% select(-c(3,10)), by = c("Sample"= "ExternalID")) %>%select(-c(2,7:11)) %>% spread(., data_type, RPKM) %>% filter(!is.na(metagenomics))
rainin.exact_matched_12_cutoff %>% filter(!is.na(metatranscriptomics)) %>% mutate(logRatio = log2(metatranscriptomics/metagenomics)) %>% ggplot(., aes(GroupAttribute, logRatio)) + 
  geom_boxplot()+geom_jitter(width=.11)
