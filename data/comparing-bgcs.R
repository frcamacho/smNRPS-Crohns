

require(tidyverse)
############################################################################################################################################################################
compare_C_clusters <- read_delim("compare-MetaHIT-HMP2-enriched-bgcs/compare-enriched-C_against-CD",
                                 col_names = F, delim = "\t")

names(compare_C_clusters)<- c("sseqid", "qseqid", "slen","qlen", "qcovs", "pident", "evalue", "bitscore", "qstart", "qend", "sstart", "send")

compare_C_clusters %>% filter(qlen<slen) %>% filter(qcovs >=95 )
compare_C_clusters %>% filter(slen<qlen) %>% mutate(sqcovs=((pmax(sstart,send)-pmin(sstart,send))/slen )*100) %>% filter( sqcovs >=95)

#C2075__NODE_793_length_16630_cov_377.864__16630__other__ANTISMASH__0_16630 42532 16630  & O2.UC23.2__NODE_123_length_53602_cov_5.17215__42532__other__ANTISMASH__8083_50615

compare_CD_clusters <- read_delim("compare-MetaHIT-HMP2-enriched-bgcs/compare-enriched-CD_against-C",
                                 col_names = F, delim = "\t")
names(compare_CD_clusters)<- c("sseqid", "qseqid", "slen","qlen", "qcovs", "pident", "evalue", "bitscore", "qstart", "qend", "sstart", "send")
compare_CD_clusters %>% filter(qlen<slen) %>% filter(qcovs >=95 )
compare_CD_clusters %>% filter(slen<qlen) %>% mutate(sqcovs=((pmax(sstart,send)-pmin(sstart,send))/slen )*100) %>% filter( sqcovs >=95)
#C4004__NODE_591_length_10305_cov_3.58927__10305__nrps__ANTISMASH__0_10305 ::: V1.CD18.3__NODE_78_length_114577_cov_12.9704__43190__nrps__ANTISMASH__49316_92506 smNRPS: smNRPS-CC 
#C3008__NODE_58_length_80320_cov_9.89467__38512__nrps__ANTISMASH__41808_80320 ::: V1.UC61.0__NODE_930_length_26002_cov_3.77986__26002__nrps__ANTISMASH__0_26002: smNRPS: smNRPS-CB 
#71.82909 sqcovs 
#C3017__NODE_742_length_38583_cov_142.623__25546__nrps__ANTISMASH__13037_38583 ::V1.CD18.3__NODE_78_length_114577_cov_12.9704__43190__nrps__ANTISMASH__49316_92506 


compare_UC_clusters <- read_delim("compare-MetaHIT-HMP2-enriched-bgcs/compare-enriched-UC_against-C",
                                  col_names = F, delim = "\t")
compare_C_UC_clusters <- read_delim("compare-MetaHIT-HMP2-enriched-bgcs/compare-enriched-C_against-UC",
                                  col_names = F, delim = "\t")
names(compare_C_UC_clusters)<- c("sseqid", "qseqid", "slen","qlen", "qcovs", "pident", "evalue", "bitscore", "qstart", "qend", "sstart", "send")
compare_C_UC_clusters %>% filter(qlen<slen) %>% filter(qcovs >=95 )

compare_C_UC_clusters %>% filter(slen<qlen) %>% mutate(sqcovs=((pmax(sstart,send)-pmin(sstart,send))/slen )*100) %>% filter( sqcovs >=95)

#C4006__NODE_15_length_222149_cov_19.3876__20894__terpene__ANTISMASH__77199_98093 ::: O2.UC50.0__NODE_43_length_109678_cov_2.66805__26723__terpene__ANTISMASH__82955_109678
########################################################
#compare BLAST reads that mapped to our smNRPS-clusters  

blast_df <- read_delim("compare-MetaHIT-HMP2-enriched-bgcs/rainin-combined-BLAST-data.txt", col_names = T, delim= "\t")
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

blast_df_filter_count<-blast_df_filter %>% group_by(sseqid) %>% count()%>% filter(n==2) %>% ungroup()#11363/34211 ~33%

blast_df_filter_dups <- blast_df_filter %>% filter(sseqid %in% blast_df_filter_count$sseqid)

# Take the read with the max pident first, then rank them by min evalue for the ties of pident
blast_df_filter_dups_max <-blast_df_filter_dups %>% group_by(sseqid) %>% top_n(1,pident) %>% mutate(rank = rank(evalue, ties.method = "min"))%>% filter(rank == 1)

#write_delim(blast_df_filter_dups_max %>% select(c(-rank, -readPos,-readCov)), "smNRPS-MetaHIT-mapped-filter-reads.txt", delim = "\t")
#####################################################
#Parse reads 
parseReadID<- function(HMMdf){
  samples<-unique(HMMdf$Sample)
  
  for (s in 1:length(samples)){
    currentSample<-samples[s]
    currentSampleResults<- HMMdf %>% filter(Sample == currentSample)
    sampleType <- unique(currentSampleResults$sampleType)
    currentSampleReads<- unique(currentSampleResults$readID)
    fileName<-paste0(currentSample,paste("-",sampleType,sep =""),"-hmm-unique-reads", ".txt", sep ="")
    metaDataFile<-paste0(currentSample,paste("-",sampleType,sep =""),"-hmm-unique-reads-MetaData", ".txt", sep ="")
    write.table(currentSampleReads,fileName, quote = F, row.names = F, col.names = F )
    write.table(currentSampleResults,metaDataFile, quote = F, row.names = F, col.names = T, sep = "\t" )
  }

