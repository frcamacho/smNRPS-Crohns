
require(tidyverse)
require(RColorBrewer)

setwd("/Users/francinecamacho/Documents/Donia_analysis/smNRPS-manuscript/smNRPS-Crohns/data")
metaphlanDF<- read_delim("metaphlan/all-combined-metaphlan_profiled.txt", col_names = T, delim = "\t")

spanish_metadata <- read_csv("sample_metadata.csv", col_names = T) # load metadata


metaphlanDF<-metaphlanDF[-c(1),]

metaphlanDF.long<-metaphlanDF%>% gather(sampleID, relativeAbundance, MH0001_profiled:V1.UC58.0_profiled) %>% separate(.,col=sampleID, into = c("Sample", "test"), sep = "_", remove = T ) %>% select(-test) 
metaphlanDF.long$relativeAbundance<-as.numeric(metaphlanDF.long$relativeAbundance)

metaphlanDF.long.metadata<-metaphlanDF.long %>% inner_join(.,spanish_metadata)

# filter taxa to only include Phylum 
metaphlan.phylumnDF<- metaphlanDF.long.metadata %>% filter(str_detect(ID, "p__") & str_detect(ID, "k__") ) %>% filter(., !grepl('c__', ID)) %>% filter(relativeAbundance != 0)

#colors for plot 
colourCount = length(unique(metaphlan.phylumnDF$ID))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

metaphlan.phylumnDF %>% ggplot(., aes(x = Sample, y = relativeAbundance, fill = ID)) + 
  geom_bar(stat = "identity") + facet_grid(~GroupAttribute, switch = "x", scales = "free_x", space = "free_x") +
  theme(panel.spacing = unit(0, "lines"), 
        strip.background = element_blank()) + xlab("Disease")+ theme_classic() + theme(
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank())+scale_fill_manual(values = getPalette(colourCount))+ ggtitle("Bacterial Phyla Relative Abundance in Spanish-MetaHit Cohort")

ggsave("MetaPhlAn-Phyla-barplot-MetaHit.eps", width = 7.3, height = 6.31, units = "in", device = "eps") 


############################## 
#Just focus on Clostridia and Ecoli genera 

metaphlan.generaDF<- metaphlanDF.long.metadata %>% filter(str_detect(ID, "p__") & str_detect(ID, "k__") & str_detect(ID, "c__") & str_detect(ID, "o__") &str_detect(ID, "f__") & str_detect(ID, "g__"))  %>% 
  filter(., !grepl('s__', ID))%>% filter(relativeAbundance != 0) #175 genera

#Clostrida/E.coli Distribution 

selected_genera <- metaphlan.generaDF %>% filter(str_detect(ID,"g__Clostridium") | str_detect(ID, "g__Escherichia"))


selected_genera %>% filter(str_detect(ID,"g__Clostridium")) %>% group_by(GroupAttribute) %>% mutate(ymean= mean(relativeAbundance)) %>% mutate(ymedian= median(relativeAbundance)) %>%  ggplot(., aes(Sample,relativeAbundance)) + geom_col()+geom_hline(aes(yintercept=ymedian), color="red", linetype="dashed")+geom_hline(aes(yintercept=ymean), color="red", linetype="solid")+theme_pubclean()+ facet_grid(~GroupAttribute, switch= "x",scales = "free_x", space = "free_x") + xlab("Disease")+  theme(
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank()) + ggtitle("g__Clostridium") + theme(plot.title = element_text(hjust = 0.5,face = "bold",size = "24"))

ggsave("MetaPhlAn-g__Clostridium-barplot-MetaHit.eps", width = 7.3, height = 6.31, units = "in", device = "eps") 


selected_genera %>% filter(str_detect(ID,"g__Escherichia")) %>% group_by(GroupAttribute) %>% mutate(ymean= mean(relativeAbundance)) %>% mutate(ymedian= median(relativeAbundance)) %>%  ggplot(., aes(Sample,relativeAbundance)) + geom_col()+geom_hline(aes(yintercept=ymedian), color="red", linetype="dashed")+geom_hline(aes(yintercept=ymean), color="red", linetype="solid")+theme_pubclean()+ facet_grid(~GroupAttribute, switch= "x",scales = "free_x", space = "free_x") + xlab("Disease")+  theme(
  axis.text.x=element_blank(),
  axis.ticks.x=element_blank()) + ggtitle("g__Escherichia") + theme(plot.title = element_text(hjust = 0.5,face = "bold",size = "24"))

ggsave("MetaPhlAn-g__Escherichia-barplot-MetaHit.eps", width = 7.3, height = 6.31, units = "in", device = "eps") 

