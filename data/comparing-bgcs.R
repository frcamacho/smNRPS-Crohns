

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