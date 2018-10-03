require("tidyverse")

# Retrieve Accession IDs for Uniprot 
ncbi_df<-read_delim("species_results-HMP2-IBD.txt", 
                    delim = "\t", col_names = TRUE)
names(ncbi_df)<-c("BGC_NAME", "TAXA_NAME", "ACC_ID", "PERC_IDENT", "COVERAGE")

filter_ncbi_df <- ncbi_df %>% filter(ACC_ID != "N/A") # 163 mapped
write_lines(unique(filter_ncbi_df$ACC_ID), "refseq_against_HMP2-IBD-SPAdes_uniqueBGCs-accession_ids-500_hits-bitscore.txt")

uniprot_df <- read_delim("uniprot-results-accesion-id-mapped.txt",
                         col_names = T, delim = "\t") %>% select( Organism, `Taxonomic lineage (PHYLUM)`, `yourlist:M20180930A7434721E10EE6586998A056CCD0537EB941995` ) 

colnames(uniprot_df) <- c( "ORGANISM", "PHYLUM", "ACCID")
uniprot_df_unique <- uniprot_df %>%  group_by(ORGANISM,PHYLUM,ACCID ) %>% distinct() 

#Output from species_profiler to combined results from NCBI and UniProt
species_profile_df<-read_delim("species_results-HMP2-IBD.txt", 
                               delim = "\t", col_names = TRUE)
combined_results <- species_profile_df %>% left_join(., uniprot_df_unique, by = c("ACC_ID" = "ACCID"))
combined_results[is.na(combined_results)] <- "N/A"
write_delim(combined_results, "species_results-HMP2-IBD_UniProt.txt", delim = "\t")
