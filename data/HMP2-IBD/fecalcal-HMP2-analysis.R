require(tidyverse)


hmp2_metadata_dna<- read_csv("HMP2-IBD/hmp2_metadata_downloaded_09_14_2018.csv",col_names = T) %>% select(., `External ID`, `Participant ID`,`data_type`, diagnosis, fecalcal) %>%
  filter(data_type == "metagenomics") 
names(hmp2_metadata_dna)[1:2] <- c("ExternalID", "ParticipantID")


hmp2_metadata_dna <- hmp2_metadata_dna%>% separate(., ExternalID, into = c("ExternalID", "extra"), sep = "_") %>% select(-c(extra))



analysis_hmp2_dna<- read_csv("HMP2-IBD/hmp2_project_metadata_2016-10-15.csv",col_names = T) %>% select(., `External ID`, `Participant ID`,`data_type`, diagnosis) %>%
  filter(data_type == "metagenomics") 
names(analysis_hmp2_dna)[1:2] <- c("ExternalID", "ParticipantID")


rainin.exact_results_DNA <- read_delim("HMP2-IBD/smNRPS-exact-unique-reads-HMP2-IBD/combined-HMP2-IBD-metagenomes-rainin-exact-unique-reads-quantifier-results.txt", col_names = F, delim = "\t")
names(rainin.exact_results_DNA) <- c("bgcName", "Coverage", "RPKM", "Sample") # add column names

rainin.exact_results_DNA_filter <- rainin.exact_results_DNA %>% filter(Coverage >=50)


smnrps_c_cd<-hmp2_metadata_dna %>% full_join(.,rainin.exact_results_DNA_filter, by = c("ExternalID" = "Sample"))
#%>% filter(diagnosis == "CD" | diagnosis == "nonIBD")


smNRPS_present <- smnrps_c_cd %>% filter(!is.na(fecalcal) & !is.na(bgcName)) #37 samples; 42
smNRPS_present$smNRPS <- "Present"

smNRPS_absent<-smnrps_c_cd %>% filter(!ExternalID %in%smNRPS_present$ExternalID) %>% filter(!is.na(fecalcal)) #286; 404
smNRPS_absent$smNRPS<- "Absent"


smNRPS_data <- rbind(smNRPS_present, smNRPS_absent) %>% filter(ExternalID %in% analysis_hmp2_dna$ExternalID)

#smnrps_c_cd %>% filter(is.na(fecalcal) & !is.na(bgcName)) %>% filter(!ExternalID %in% smNRPS_data$ExternalID)
#smnrps_c_cd %>% filter(is.na(fecalcal) & is.na(bgcName)) %>% filter(!ExternalID %in% smNRPS_data$ExternalID)
#how many samples are missing metadata 

#Total samples unique == 1151
#Total samples with any smNRPS and fecal data== 37
#Total samples with any smNRPS and no fecal data == 11
#Total samples with no smNRPS and fecal data == 286
#Total samples with no smNRPS and no fecal data == 817
#####################################################################################
#These are did not have any BLAST hits 
#HSM5MD4P
#HSM5MD8H
#MSM5LLD6
#MSM5LLHE
#MSM5LLIS
#PSM6XBRK
ggplot(smNRPS_data, aes(x = smNRPS, y = fecalcal)) + geom_boxplot(outlier.colour = "red") + geom_jitter(width = 0.1, alpha = 0.25)
ggplot(smNRPS_data, aes(x = diagnosis, y = fecalcal)) + geom_boxplot() + geom_jitter(width = 0.1, alpha = 0.25) + facet_wrap(~smNRPS)
#ggplot(smNRPS_data, aes(x = diagnosis, y = fecalcal)) + geom_violin(adjust=1.2)+ geom_jitter(width=0.2, alpha=0.5)+ facet_wrap(~smNRPS)
library("ggpubr")
smNRPS_data$RPKM[is.na(smNRPS_data$RPKM)] <-0
ggscatter(smNRPS_data, x = "fecalcal", y = "RPKM", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "fecalcal", ylab = "RPKM")

res2 <-cor.test(smNRPS_data$fecalcal, smNRPS_data$RPKM,  method = "kendall")
#Just non-IBD
smNRPS_data %>% filter(diagnosis == "nonIBD") %>% ggscatter(., x = "fecalcal", y = "RPKM", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "fecalcal", ylab = "RPKM")

smNRPS_data %>% filter(diagnosis == "nonIBD")  %>% ggplot(., aes(x = smNRPS, y = fecalcal)) + geom_boxplot(outlier.colour = "red") + geom_jitter(width = 0.1, alpha = 0.25)


wilcox.test(fecalcal ~ smNRPS, data=smNRPS_data, paired=T,subset = diagnosis %in% "nonIBD" )

kruskal.test(fecalcal ~ smNRPS, data = smNRPS_data %>% filter(diagnosis == "nonIBD"))
