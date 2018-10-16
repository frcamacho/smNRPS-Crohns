require(tidyverse)

setwd("/tigress/DONIA/fcamacho/Crohns-manuscript/MetaHIT_Gut/HMP2-IBD-spanish-quanification")
raininBlastfiles<- list.files(path= "/tigress/DONIA/fcamacho/Crohns-manuscript/MetaHIT_Gut/HMP2-IBD-spanish-quanification",
                              pattern="Rainin-clusters-exact-vs-*", full.names=T, recursive=TRUE)


parseBlastFiles<-function(files){
  
  df<-data.frame()
  for (x in 1:length(files)) {
    if (!file.size(files[x]) == 0) {
      t <- read.delim(files[x], header=F, stringsAsFactors = FALSE) # load file
      filename<-strsplit(files[x], "/")[[1]][12] # "SRS017916-stool-against-bt1fas-PARSED.txt"
      names(t)<-(c("sseqid","slen", "sstart", "send", "qseqid", "qlen", "qstart", "qend", "pident", "evalue"))
      t$GroupAttribute<-strsplit(files[x], "/")[[1]][9] 
      t$ParticipantID<-strsplit(files[x], "/")[[1]][10] 
      t$data_type<-strsplit(files[x], "/")[[1]][8] 
      t$Sample<-strsplit(files[x], "/")[[1]][12] 
      df<-rbind(t, df)
    } else {
      paste("Error file is empty",files[x], sep = ": ")
    }
  }
  return (df)
}

raininBlastDF<-parseBlastFiles(raininBlastfiles)
write_delim(raininBlastDF, "rainin-combined-HMP-2-IBD-BLAST-data.txt", col_names = T, delim ="\t")
