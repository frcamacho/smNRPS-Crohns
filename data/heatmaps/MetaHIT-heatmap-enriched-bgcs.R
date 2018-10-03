#load quantifier data 
require(tidyverse)
require(pheatmap)
setwd("/Users/francinecamacho/Documents/Donia_analysis/smNRPS-manuscript/smNRPS-Crohns/data")
quantifier_df <- read_delim("combined-MetaHIT_Gut-spanish-quantifier-results.txt", col_names = F, delim = "\t")
names(quantifier_df) <- c("bgcName", "Coverage", "RPKM", "Sample") # add column names

filter_quantifier_df <- quantifier_df %>% filter(Coverage >= 50)  # filter by coverage 
sample_metadata <- read_csv("sample_metadata.csv", col_names = T) # load metadata
filter_quantifier_df <- filter_quantifier_df %>% inner_join(., sample_metadata, by = "Sample")

lefse_CD_C <- read_csv("supplementary-tables/lefse_CD_C_supp_table_1.csv",col_names = T)
lefse_CD_C_discriminative_features <- lefse_CD_C %>% filter(class_discriminitive!= "-") 

annotations<-sample_metadata %>% as.data.frame()
rownames(annotations)<-annotations$Sample
annotations_col<-annotations[-c(1)]

annotations_row <-lefse_CD_C %>% filter(class_discriminitive!="-") %>% select(c(BGC,class_discriminitive)) %>% as.data.frame()
rownames(annotations_row)<-annotations_row$BGC
annotations_row<-annotations_row[-c(1)]
#######################################################
filter_quantifier_df_enriched_CD_C <- filter_quantifier_df %>% filter(bgcName %in% lefse_CD_C_discriminative_features$BGC) %>% filter(GroupAttribute!="UC")

quantile(filter_quantifier_df_enriched_CD_C$RPKM, probs = c(.05, .90)) #    5% = 0.3330851  90% =12.3655056 


heatmap_RS<-filter_quantifier_df_enriched_CD_C  %>% select(c(-Coverage))
heatmap_RS$RPKM[heatmap_RS$RPKM > 12.3655056 ] <- 12.3655056  #recode 
#heatmap_RS$RPKM[heatmap_RS$RPKM > 20 ] <- 20 

data_wide <- dcast(heatmap_RS, bgcName ~ Sample, value.var="RPKM")
row.names(data_wide)<- data_wide$bgcName

data_wide<-data_wide[,-1]
data_wide[is.na(data_wide)] <- 0
data_wide<-as.matrix(data_wide)
#############################################
#Colors for heatmap 
bk = unique(c(seq(0,.05, length=100), seq(.06,12.3655056 ,length=1000)))
#bk = unique(c(seq(0,.05, length=100), seq(.06,90 ,length=1000)))
#bk = unique(c(seq(0,.05, length=100), seq(.06,20 ,length=1000)))

hmcols<- colorRampPalette(c("white","blue","green","yellow","orange","red"))(length(bk)-1)


# pheatmap::pheatmap(data_wide, scale = "none",breaks=bk, color=hmcols,annotation_col= annotations_col, annotation_row = annotations_row, show_rownames=TRUE, show_colnames = FALSE,fontsize_row = 6,
#                    clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", cellwidth = 5,cellheight = 5,
#                    clustering_method = "ward.D",main="Ward.D Clustering with Correlation")

pheatmap::pheatmap(data_wide, scale = "none",breaks=bk, color=hmcols,annotation_col= annotations_col, annotation_row = annotations_row, show_rownames=TRUE, show_colnames = FALSE,fontsize_row = 6,
                   clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", cellwidth = 5,cellheight = 5,
                   clustering_method = "ward.D2",main="Ward.D2 Clustering with Correlation")
pheatmap::pheatmap(data_wide, scale = "none",breaks=bk, color=hmcols,annotation_col= annotations_col, annotation_row = annotations_row, show_rownames=TRUE, show_colnames = FALSE,fontsize_row = 6,
                   clustering_distance_rows = "binary", clustering_distance_cols = "binary", cellwidth = 5,cellheight = 5,
                   clustering_method = "ward.D",main="Ward.D Clustering with Binary")


###########################################
#With UC 

filter_quantifier_df_enriched_all <- filter_quantifier_df %>% filter(bgcName %in% lefse_CD_C_discriminative_features$BGC) 


heatmap_RS_all<-filter_quantifier_df_enriched_all  %>% select(c(-Coverage))

quantile(filter_quantifier_df_enriched_all$RPKM, probs = c(.05, .90)) # 0.3323036 12.1994997 
#heatmap_RS_all$RPKM[heatmap_RS_all$RPKM > 12.1994997 ] <- 12.1994997  #recode 
heatmap_RS_all$RPKM[heatmap_RS_all$RPKM > 20 ] <- 20 
data_wide_all <- dcast(heatmap_RS_all, bgcName ~ Sample, value.var="RPKM")
row.names(data_wide_all)<- data_wide_all$bgcName



data_wide_all<-data_wide_all[,-1]
data_wide_all[is.na(data_wide_all)] <- 0
data_wide_all<-as.matrix(data_wide_all)

bk_2 = unique(c(seq(0,.05, length=100), seq(.06,20 ,length=1000)))
#bk_2 = unique(c(seq(0,.05, length=100), seq(.06,12.1994997 ,length=1000)))
#bk = unique(c(seq(0,.05, length=100), seq(.06,90 ,length=1000)))
hmcols_2<- colorRampPalette(c("white","blue","green","yellow","orange","red"))(length(bk_2)-1)

# pheatmap::pheatmap(data_wide_all, scale = "none",breaks=bk_2, color=hmcols,annotation_col= annotations_col, annotation_row = annotations_row, show_rownames=TRUE, show_colnames = FALSE,fontsize_row = 6,
#                    clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", cellwidth = 3,cellheight = 3,
#                    clustering_method = "ward.D",main="Ward.D Clustering with Correlation")
# 
# pheatmap::pheatmap(data_wide_all, scale = "none",breaks=bk_2, color=hmcols,annotation_col= annotations_col, annotation_row = annotations_row, show_rownames=TRUE, show_colnames = FALSE,fontsize_row = 6,
#                    clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", cellwidth = 3,cellheight = 3,
#                    clustering_method = "ward.D2",main="Ward.D2 Clustering with Correlation")
pheatmap::pheatmap(data_wide_all, scale = "none",breaks=bk_2, color=hmcols_2,annotation_col= annotations_col, annotation_row = annotations_row, show_rownames=TRUE, show_colnames = FALSE,fontsize_row = 6,
                   clustering_distance_rows = "binary", clustering_distance_cols = "binary", cellwidth = 3,cellheight = 3,
                   clustering_method = "ward.D",main="Ward.D Clustering with Binary")

############################################
#using top 30 enriched BGCs in each cohort (Crohns and healthy)

lefse_CD_C_features_top30 <- lefse_CD_C_discriminative_features %>% group_by(class_discriminitive) %>% arrange(desc(LDA_score)) %>% group_by(class_discriminitive) %>% top_n(30, LDA_score)
df_enriched_CD_C_top30 <- filter_quantifier_df %>% filter(bgcName %in% lefse_CD_C_features_top30$BGC) %>% filter(GroupAttribute!="UC")

quantile(df_enriched_CD_C_top30$RPKM, probs = c(.05, .90)) #    5% = 0.4223711  90% =17.7351081 


heatmap_RS_top30<-df_enriched_CD_C_top30  %>% select(c(-Coverage))

heatmap_RS_top30$RPKM[heatmap_RS_top30$RPKM > 20 ] <- 20 
data_wide_top30 <- dcast(heatmap_RS_top30, bgcName ~ Sample, value.var="RPKM")
row.names(data_wide_top30)<- data_wide_top30$bgcName



data_wide_top30<-data_wide_top30[,-1]
data_wide_top30[is.na(data_wide_top30)] <- 0
data_wide_top30<-as.matrix(data_wide_top30)

bk_3 = unique(c(seq(0,.05, length=100), seq(.06,20 ,length=1000)))
hmcols_3<- colorRampPalette(c("white","blue","green","yellow","orange","red"))(length(bk_3)-1)

pheatmap::pheatmap(data_wide_top30, scale = "none",breaks=bk_3, color=hmcols_3,annotation_col= annotations_col, annotation_row = annotations_row, show_rownames=TRUE, show_colnames = FALSE,fontsize_row = 6,
                   clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", cellwidth = 3,cellheight = 3,
                   clustering_method = "ward.D",main="Ward.D Clustering with correlation")
pheatmap::pheatmap(data_wide_top30, scale = "none",breaks=bk_3, color=hmcols_3,annotation_col= annotations_col, annotation_row = annotations_row, show_rownames=TRUE, show_colnames = FALSE,fontsize_row = 6,
                   clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", cellwidth = 3,cellheight = 3,
                   clustering_method = "ward.D2",main="Ward.D2 Clustering with correlation")

pheatmap::pheatmap(data_wide_top30, scale = "none",breaks=bk_3, color=hmcols_3,annotation_col= annotations_col, annotation_row = annotations_row, show_rownames=TRUE, show_colnames = FALSE,fontsize_row = 6,
                   clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", cellwidth = 3,cellheight = 3,
                   clustering_method = "mcquitty",main="Mcquitty with correlation")
pheatmap::pheatmap(data_wide_top30, scale = "none",breaks=bk_3, color=hmcols_3,annotation_col= annotations_col, annotation_row = annotations_row, show_rownames=TRUE, show_colnames = FALSE,fontsize_row = 6,
                   clustering_distance_rows = "correlation", clustering_distance_cols = "binary", cellwidth = 3,cellheight = 3,
                   clustering_method = "complete",main="Complete with correlation")
