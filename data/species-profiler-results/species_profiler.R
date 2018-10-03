#Species profiler 
require(tidyverse)

species_profiler<-function(blastdf, coverage_cutoff, bgc_list){
  
  blastdf_cutoff <- blastdf %>% filter(qcovs >=coverage_cutoff)
  bgc_list_df <- as.data.frame(bgc_list)
  names(bgc_list_df) <- c("BGC_NAME")
  filter_species_df<-blastdf_cutoff %>% group_by(qseqid) %>% filter(qcovs == max(qcovs)) %>% filter(bitscore == max(bitscore)) %>% mutate(the_rank  = rank(-qcovs, ties.method = "first")) %>%
    filter(the_rank == 1) %>% select(-the_rank)%>% ungroup() # first filter multiple hits by qcovs then by bitscore and then if there are ties than take the first one 
  results_df <-filter_species_df %>% select(c(qseqid,stitle,sacc,pident,qcovs))
  names(results_df) <- c("BGC_NAME", "TAXA_NAME", "ACC_ID", "PERC_IDENT", "COVERAGE")
  full_results <- results_df %>% full_join(.,bgc_list_df,by = "BGC_NAME")
  full_results$TAXA_NAME[is.na( full_results$TAXA_NAME)]<-"N/A"
  full_results$ACC_ID[is.na( full_results$ACC_ID)]<-"N/A"
  full_results$PERC_IDENT[is.na( full_results$PERC_IDENT)]<-0.0
  full_results$COVERAGE[is.na( full_results$COVERAGE)]<-0.0
  return(full_results)
}

species_df<- read_delim("/Users/francinecamacho/Documents/Donia_analysis/smNRPS-manuscript/duplication-test/species-profiler/refseq_against_MetaHIT-SPAdes_uniqueBGCs_tabular-file-tophit_500_hits-bitscore",
                        delim = "\t", col_names = F) 
names(species_df)<-c( "sseqid","stitle", "sacc", "qseqid", "qlen", "qcovs", "pident", "evalue", "bitscore", "qstart", "qend")

cluster_ids <- read_lines("MetaHIT_Gut-spanish_uniqueBGCs-list.txt")
data_df<- species_profiler(species_df, 50, cluster_ids)
write_delim(data_df, "species_results-MetaHIT-spanish.txt", delim="\t", col_names = T)