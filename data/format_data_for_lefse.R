library(reshape2)
library(dplyr)
#Function to convert data structure back to lefse format 
format_data_for_lefse<-function(df){
  dcast_df<-dcast(df, bgcName ~ Sample + GroupAttribute,value.var = "RPKM")
  
  
  temprow <- matrix(c(rep.int(df$GroupAttribute,length(dcast_df))),nrow=1,ncol=length(dcast_df))
  
  
  newrow <- data.frame(temprow)
  colnames(newrow) <- colnames(dcast_df)
  i<-1
  for (i in 1:ncol(newrow)){
    c.name<-colnames(newrow)[i]
    if(c.name == "bgcName"){
      newrow[,i]<-"Disease"
    }
    else{
      diseaseStatus<-strsplit(c.name,split = "_")[[1]][2]
      
      if (diseaseStatus !=newrow[1,i]){
        newrow[,i]<-diseaseStatus
      }}
    i<-i+1
  }
  
  data <- rbind(dcast_df,newrow)
  
  arrange.data<-dplyr::arrange(data, -row_number())
  
  arrange.data[is.na(arrange.data)]<-0
  names(arrange.data) <-sapply(strsplit(names(arrange.data),"\\_"), function(x) paste0(head(x,-1),collapse="_") )
  
  names(arrange.data)[1]<-"Sample" # the sapply function removed the ID cols so I am putting it back
  return(arrange.data)
}