require(tidyverse)
#compare BLAST reads that mapped to our smNRPS-clusters  

blast_df <- read_delim("smNRPS-exact-unique-reads-HMP2-IBD/rainin-combined-HMP-2-IBD-BLAST-data.txt", col_names = T, delim= "\t")
filterBlastReads<-function(blastDF){
  
  filterDF<-data.frame()
  for (i in 1:nrow(blastDF)){
    internal <- NULL
    row <- blastDF[i,]
    smin <- min(row$sstart, row$send)
    smax <- max(row$sstart, row$send)
    qmin <- min(row$qstart, row$qend)
    qmax <- max(row$qstart, row$qend)
    percentReadCoveredByCluster <- ((smax - smin + 1) / (row$slen) )* 100
    if(qmin - 1 >= smin - 1  && row$qlen - qmax >= row$slen - smax){
      internal <- TRUE
    }else{
      internal <- FALSE}
    
    if (internal == TRUE){
      if ((row$pident >= 95.0) && (percentReadCoveredByCluster>= 90)){
        row$readPos <- "internal"
        row$readCov <- percentReadCoveredByCluster
        filterDF<-rbind(row,filterDF )
      }} else if (internal == FALSE){
        if ((row$pident >= 95.0) && (percentReadCoveredByCluster>= 50)){
          row$readPos <- "edge"
          row$readCov <- percentReadCoveredByCluster
          filterDF<-rbind(row,filterDF )}}
  }
  return(filterDF)
}


blast_df_filter<-filterBlastReads( blast_df %>% filter(qseqid != "smNRPS-BP"))

blast_df_filter_DNA<- blast_df_filter %>% filter(data_type == "HMP2-IBD-metagenomes" )
blast_df_filter_RNA<- blast_df_filter %>% filter(data_type == "HMP2-IBD-metatranscriptomes" )
#blast_df_filter_count<-blast_df_filter %>% group_by(sseqid) %>% count()%>% filter(n==2) %>% ungroup()#12110/48972 ~25%
# Take the read with the max pident first, then rank them by min evalue for the ties of pident

blast_df_filter_DNA_dups_max <- blast_df_filter_DNA%>% group_by(sseqid) %>% top_n(1,pident) %>% mutate(rank = rank(evalue, ties.method = "min"))%>% filter(rank == 1)%>%
  select(c(-rank))%>%ungroup()

blast_df_filter_DNA_dups_max_ties<- blast_df_filter_DNA_dups_max %>% group_by(sseqid) %>% mutate(rank = rank(evalue, ties.method = "random"))%>% filter(rank == 1)

blast_df_filter_RNA_dups_max <- blast_df_filter_RNA%>% group_by(sseqid) %>% top_n(1,pident) %>% mutate(rank = rank(evalue, ties.method = "min"))%>% filter(rank == 1)%>%
  select(c(-rank))%>%ungroup()

blast_df_filter_RNA_dups_max_ties<- blast_df_filter_RNA_dups_max %>% group_by(sseqid) %>% mutate(rank = rank(evalue, ties.method = "random"))%>% filter(rank == 1)


write_delim(blast_df_filter_DNA_dups_max_ties %>% select(c(-rank, -readPos,-readCov)), "smNRPS-HMP2-IBD-mapped-metagenomes-filter-reads.txt", delim = "\t")
write_delim(blast_df_filter_RNA_dups_max_ties %>% select(c(-rank, -readPos,-readCov)), "smNRPS-HMP2-IBD-mapped-metatranscriptomes-filter-reads.txt", delim = "\t")

#####################################################
#Parse reads #Rainin-clusters-exact-vs-O2.UC43.2-tabular-blast

sample_smnrps_df_dna <- read_delim("smNRPS-HMP2-IBD-mapped-metagenomes-filter-reads.txt", col_names = T, delim = "\t")
sample_smnrps_df_rna <- read_delim("smNRPS-HMP2-IBD-mapped-metatranscriptomes-filter-reads.txt", col_names = T, delim = "\t")
#Combine the smNRPS reads for each sample
parseReadID<- function(df, dir_path){
  setwd(dir_path)
  samples<-unique(df$Sample)
  
  for (s in 1:length(samples)){
    currentSample<-samples[s]
    currentSampleResults<- df %>% filter(Sample == currentSample)
    #fileName<-paste0("smNRPS-exact-unique-reads-vs-",currentSample,"-tabular-blast", ".txt", sep="")
    fileName<-paste0("smNRPS-exact-unique-reads-updated-vs-",currentSample,"-tabular-blast", ".txt", sep="")
    write_delim(currentSampleResults %>% select(-c(Sample, GroupAttribute,ParticipantID,data_type)),fileName, delim="\t", col_names = F )
  }}


#parseReadID(sample_smnrps_df_dna, "/Users/francinecamacho/Documents/Donia_analysis/smNRPS-manuscript/smNRPS-Crohns/data/HMP2-IBD/smNRPS-exact-unique-reads-HMP2-IBD/smNRPS-exact-unique-HMP2-metagenomes")
#parseReadID(sample_smnrps_df_rna, "/Users/francinecamacho/Documents/Donia_analysis/smNRPS-manuscript/smNRPS-Crohns/data/HMP2-IBD/smNRPS-exact-unique-reads-HMP2-IBD/smNRPS-exact-unique-HMP2-metatranscriptomes" )

########################
#combined reads to one clusters 
sample_smnrps_df_dna$qseqid <-"smNRPS-1/2"
sample_smnrps_df_dna$qlen<-5910

sample_smnrps_df_rna$qseqid <-"smNRPS-1/2"
sample_smnrps_df_rna$qlen<-5910
parseReadID(sample_smnrps_df_dna, "/Users/francinecamacho/Documents/Donia_analysis/smNRPS-manuscript/smNRPS-Crohns/data/HMP2-IBD/smNRPS-exact-unique-reads-HMP2-IBD-updated/smNRPS-exact-unique-HMP2-metagenomes/" )
parseReadID(sample_smnrps_df_rna, "/Users/francinecamacho/Documents/Donia_analysis/smNRPS-manuscript/smNRPS-Crohns/data/HMP2-IBD/smNRPS-exact-unique-reads-HMP2-IBD-updated/smNRPS-exact-unique-HMP2-metatranscriptomes/" )
