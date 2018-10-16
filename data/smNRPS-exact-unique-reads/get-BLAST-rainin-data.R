require(tidyverse)

setwd("/tigress/DONIA/fcamacho/Crohns-manuscript/compare-MetaHIT-HMP2-enriched-bgcs")
raininBlastfiles<- list.files(path= "/tigress/DONIA/fcamacho/Crohns-manuscript/MetaHIT_Gut/quantification-analysis",
                             pattern="Rainin-clusters-exact-vs-*", full.names=T, recursive=TRUE)


parseBlastFiles<-function(files){
  
  df<-data.frame()
  for (x in 1:length(files)) {
    if (!file.size(files[x]) == 0) {
      t <- read.delim(files[x], header=F, stringsAsFactors = FALSE) # load file
      filename<-strsplit(files[x], "/")[[1]][10] # "SRS017916-stool-against-bt1fas-PARSED.txt"
      names(t)<-(c("sseqid","slen", "sstart", "send", "qseqid", "qlen", "qstart", "qend", "pident", "evalue"))
      t$Sample<-strsplit(files[x], "/")[[1]][9] 
      df<-rbind(t, df)
    } else {
      paste("Error file is empty",files[x], sep = ": ")
    }
  }
  return (df)
}

raininBlastDF<-parseBlastFiles(raininBlastfiles)
write_delim(raininBlastDF, "rainin-combined-BLAST-data.txt", col_names = T, delim ="\t")
