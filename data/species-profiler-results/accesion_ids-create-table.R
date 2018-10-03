require("tidyverse")

# Retrieve Accession IDs for Uniprot 
ncbi_df<-read_delim("species_results-MetaHIT-spanish.txt", 
                    delim = "\t", col_names = TRUE)
names(ncbi_df)<-c("BGC_NAME", "TAXA_NAME", "ACC_ID", "PERC_IDENT", "COVERAGE")

filter_ncbi_df <- ncbi_df %>% filter(ACC_ID != "N/A") # 400 mapped
write_lines(unique(filter_ncbi_df$ACC_ID), "refseq_against_MetaHIT-SPAdes_uniqueBGCs-accession_ids-500_hits-bitscore.txt")

uniprot_df <- read_delim("uniprot-results-accesion-id-mapped.txt",
                         col_names = T, delim = "\t") %>% select( Organism, `Taxonomic lineage (PHYLUM)`, `yourlist:M20180928AAFB7E4D2F1D05654627429E83DA5CCE0EA8F1R` ) 

colnames(uniprot_df) <- c( "ORGANISM", "PHYLUM", "ACCID")
uniprot_df_unique <- uniprot_df %>%  group_by(ORGANISM,PHYLUM,ACCID ) %>% distinct() 

#Output from species_profiler to combined results from NCBI and UniProt
species_profile_df<-read_delim("species_results-MetaHIT-spanish.txt", 
                    delim = "\t", col_names = TRUE)
combined_results <- species_profile_df %>% left_join(., uniprot_df_unique, by = c("ACC_ID" = "ACCID"))
combined_results[is.na(combined_results)] <- "N/A"
write_delim(combined_results, "species_results-MetaHIT-spanish_UniProt.txt", delim = "\t")
