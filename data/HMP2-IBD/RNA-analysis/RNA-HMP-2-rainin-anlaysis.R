library(tidyverse)
library(ggsci)
library(ggpubr)
library(reshape2)

setwd("/Users/francinecamacho/Documents/Donia_analysis/smNRPS-manuscript/smNRPS-Crohns/data/HMP2-IBD/")

## Download BLAST data 
blast_data<- read_delim("smNRPS-HMP2-IBD-mapped-metatranscriptomes-filter-reads.txt", col_names = T, delim = "\t") #1579

### Get annotation files 
cb_annotate <- read_csv("Annotation-table-CB.csv", col_names = T)
cc_annotate <- read_csv("Annotation-table-CC.csv", col_names = T)
##################################
#FUnction to map reads to Gene annotation tables 
mapBlastReads<-function(metadataDF, blastDF){
  results<-data.frame()
  for(f in 1:nrow(metadataDF)){
    domainRow <-metadataDF[f,]
    geneResults<-subset(blastDF, blastDF$qstart %in% domainRow$Minimum:domainRow$Maximum | blastDF$qend %in% domainRow$Minimum:domainRow$Maximum )
    if (nrow(geneResults)!= 0){
      geneResults$Gene<-domainRow$Name
      results<-rbind(results,geneResults)
    }
  }
  return(results)
  
}
##################################
cb_map_data <- mapBlastReads(cb_annotate, blast_data %>% filter(qseqid== "smNRPS-CB"))
cb_map_data_gene_count <-cb_map_data %>% select(Sample, Gene) %>% distinct() %>% count(Sample)
names(cb_map_data_gene_count)[2] <- c("totalGenes")
cb_map_count_data <- cb_map_data %>% inner_join(.,cb_map_data_gene_count )%>% filter(totalGenes >=2) %>% ungroup()

cb_final<-cb_map_count_data  %>% group_by(qseqid, Sample,GroupAttribute) %>% count() %>% arrange(desc(n)) %>% as.data.frame() %>% 
  inner_join(.,cb_map_data_gene_count )
names(cb_final)<-c("smNRPS_type", "Sample", "disease_type", "read_count", "gene_count")

###########################################################################################################

cc_map_data <- mapBlastReads(cc_annotate, blast_data %>% filter(qseqid== "smNRPS-CC"))
cc_map_data_gene_count <-cc_map_data %>% select(Sample, Gene) %>% distinct() %>% count(Sample)
names(cc_map_data_gene_count)[2] <- c("totalGenes")
cc_map_count_data <- cc_map_data %>% inner_join(.,cc_map_data_gene_count )%>% filter(totalGenes >=2) %>% ungroup()

cc_final<-cc_map_count_data  %>% group_by(qseqid, Sample,GroupAttribute) %>% count() %>% arrange(desc(n)) %>% as.data.frame() %>% 
  inner_join(.,cc_map_data_gene_count )
names(cc_final)<-c("smNRPS_type", "Sample", "disease_type", "read_count", "gene_count")

final_df <- rbind(cb_final,cc_final)

write_csv(final_df, "HMP2-smNRPS-RNA-read_counts.csv", col_names = T)
