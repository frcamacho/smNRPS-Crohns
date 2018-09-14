require("tidyverse")

# Retrieve Accession IDs for Uniprot 
ncbi_df<-read_delim("refseq_against_MetaHIT-SPAdes_uniqueBGCs_tabular-file-tophit_90-50", 
           delim = "\t", col_names = FALSE)
colnames(ncbi_df) <- c("sseqid", "stitle", "sacc", "qseqid", "qlen", "qcovs", "pident", "evalue","qstart", "qend")

filter_ncbi_df <- ncbi_df %>% filter(qcovs >=50 & pident >=90) # 278 mapped
write_lines(unique(filter_ncbi_df$sacc), "refseq_against_MetaHIT-SPAdes_uniqueBGCs-accession_ids-tophit_90-50.txt")

uniprot_df <- read_delim("uniprot-results-accesion-id-mapped.txt",
                         col_names = T, delim = "\t") %>% select( Organism, `Taxonomic lineage (PHYLUM)`, `yourlist:M20180914AAFB7E4D2F1D05654627429E83DA5CCE0CF9567` ) 

colnames(uniprot_df) <- c( "ORGANISM", "PHYLUM", "ACCID")
uniprot_df_unique <- uniprot_df %>%  group_by(ORGANISM,PHYLUM,ACCID ) %>% distinct()  %>% separate(., ACCID, into = c("OtherACCID", "ACCID"), 
                                                                                                   sep = ","  )

#Output from species_profiler to combined results from NCBI and UniPTO
species_profile_df<-read_delim("species_results-MetaHIT-spanish.txt", 
                    delim = "\t", col_names = TRUE)
combined_results <- species_profile_df %>% left_join(., uniprot_df_unique, by = c("ACC_ID" = "ACCID"))
combined_results[is.na(combined_results)] <- "N/A"
write_delim(combined_results, "species_results-MetaHIT-spanish_UniProt.txt", delim = "\t")
