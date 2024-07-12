# load and install required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace('readr', quietly = TRUE)) install.packages('readr'); library('readr')
if (!requireNamespace('GenomicRanges', quietly = TRUE)) BiocManager::install('GenomicRanges'); library('GenomicRanges')
if (!requireNamespace('dplyr', quietly = TRUE)) install.packages('dplyr'); library('dplyr')
if (!requireNamespace('rtracklayer', quietly = TRUE)) BiocManager::install('rtracklayer'); library('dplyr')
if (!requireNamespace('stringr', quietly = TRUE)) install.packages('stringr'); library('stringr')


featCountFolder<-commandArgs(TRUE)[1]
DNXMatOut<-commandArgs(TRUE)[2]
featCountFiles<-paste0(featCountFolder,"/*.txt")

# collect files in "~/dnx_projects/EXP23003664_aged_HDF/peakCountBySample_broadPeaks_featureCounts"
DNX_exp_data<-system(paste0("ls -d ",featCountFiles),intern = TRUE)
names(DNX_exp_data)<-strsplit(DNX_exp_data,"featureCounts/") %>% sapply("[[",2) %>% strsplit(".markdu") %>% sapply("[[",1)
DNX_exp_data

# loop through each DNX data and read in
DNX_data<-list()
for (i in 1:length(DNX_exp_data)){
  print(i)
  DNX_data[[names(DNX_exp_data)[i]]]<-read.delim(DNX_exp_data[i],comment.char='#')
  
  colnames(DNX_data[[i]])[ncol(DNX_data[[i]])]<-"readCount"
}
# normalize RNA expression in each DNA_data
for (i in 1:length(DNX_data)){
  print(i)
  DNX_data[[i]]$norm<-log2(DNX_data[[i]]$readCount/sum(DNX_data[[i]]$readCount)*1e6+1)
}
# keep only standard chromosomes in DNX_data
for (i in 1:length(DNX_data)){
  print(i)
  DNX_data[[i]]<-DNX_data[[i]][DNX_data[[i]]$Chr %in% c(1:22,"X","Y"),]
}

# load chromatin annotations
library(readr)
hg38lift_genome_100_browser <- data.frame(read_table("/weka/mwang/dnx_projects/full_stack_ChromHMM_annotations/hg38lift_genome_100_browser.bed.gz"))
chromAnnoMeta <- utils::read.csv("/weka/mwang/dnx_projects/full_stack_ChromHMM_annotations/state_annotations_processed.csv.gz")
rownames(chromAnnoMeta) <- chromAnnoMeta$mneumonics

hg38lift_genome_100_browser$length<-unlist(hg38lift_genome_100_browser[,3]-hg38lift_genome_100_browser[,2])
hg38lift_genome_100_browser$group<-chromAnnoMeta[as.character(hg38lift_genome_100_browser[,4]),"Group"]
nrow(hg38lift_genome_100_browser)
# remove "chr" in chromosome names
hg38lift_genome_100_browser$track<-gsub("chr","",hg38lift_genome_100_browser$track)
colnames(hg38lift_genome_100_browser)<-c("chrom","start","end","name","strand","length","group")
hg38lift_genome_100_browser$strand<-"."



# for each unique hg38lift_genome_100_browser$name, sum DNX_data$readCount
# initialize a list with length of DNX_data
DNX_data_count_anno<-data.frame(matrix(ncol=length(DNX_data),nrow=length(unique(hg38lift_genome_100_browser$name))+1))
colnames(DNX_data_count_anno)<-names(DNX_data)
rownames(DNX_data_count_anno)<-c("uannotated",unique(hg38lift_genome_100_browser$name))
for(i in rownames(DNX_data_count_anno)){
  for(j in colnames(DNX_data_count_anno)){
    DNX_data_count_anno[i,j]<-sum(DNX_data[[j]]$readCount[grepl(i,DNX_data[[j]]$Geneid)])
  }
}

# save as tsv
write.table(DNX_data_count_anno,file=DNXMatOut,quote=F,row.names=T,col.names=T,sep="\t")
