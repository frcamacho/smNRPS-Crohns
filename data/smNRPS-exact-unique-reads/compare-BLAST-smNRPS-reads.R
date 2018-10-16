require(tidyverse)
#compare BLAST reads that mapped to our smNRPS-clusters  

blast_df <- read_delim("smNRPS-exact-unique-reads/rainin-combined-BLAST-data.txt", col_names = T, delim= "\t")
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

#blast_df_filter_count<-blast_df_filter %>% group_by(sseqid) %>% count()%>% filter(n==2) %>% ungroup()#11363/34211 ~33%

#blast_df_filter_dups <- blast_df_filter %>% filter(sseqid %in% blast_df_filter_count$sseqid)

# Take the read with the max pident first, then rank them by min evalue for the ties of pident
blast_df_filter_dups_max <-blast_df_filter %>% group_by(sseqid) %>% top_n(1,pident) %>% mutate(rank = rank(evalue, ties.method = "min"))%>% filter(rank == 1)
#write_delim(blast_df_filter_dups_max %>% select(c(-rank, -readPos,-readCov)), "smNRPS-MetaHIT-mapped-filter-reads.txt", delim = "\t")
#####################################################
#Parse reads #Rainin-clusters-exact-vs-O2.UC43.2-tabular-blast

sample_smnrps_df <- read_delim("smNRPS-exact-unique-reads/smNRPS-MetaHIT-mapped-filter-reads.txt", col_names = T, delim = "\t")
#Combine the smNRPS reads for each sample
parseReadID<- function(df){
  #setwd("/Users/francinecamacho/Documents/Donia_analysis/smNRPS-manuscript/smNRPS-Crohns/data/smNRPS-exact-unique-reads/BLAST-unique-sample-data")
  setwd("/Users/francinecamacho/Documents/Donia_analysis/smNRPS-manuscript/smNRPS-Crohns/data/smNRPS-exact-unique-reads-updated")
  samples<-unique(df$Sample)
  
  for (s in 1:length(samples)){
    currentSample<-samples[s]
    currentSampleResults<- df %>% filter(Sample == currentSample)
    #fileName<-paste0("smNRPS-exact-unique-reads-vs-",currentSample,"-tabular-blast", ".txt", sep="")
    fileName<-paste0("smNRPS-exact-unique-reads-updated-vs-",currentSample,"-tabular-blast", ".txt", sep="")
    write_delim(currentSampleResults %>% select(-c(Sample)),fileName, delim="\t", col_names = F )
  }}

#parseReadID(sample_smnrps_df)


sample_smnrps_df$qseqid <- "smNRPS-1/2"
sample_smnrps_df$qlen<-5910
parseReadID(sample_smnrps_df)
