# load and install required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace('readr', quietly = TRUE)) install.packages('readr'); library('readr')
if (!requireNamespace('GenomicRanges', quietly = TRUE)) BiocManager::install('GenomicRanges'); library('GenomicRanges')
if (!requireNamespace('dplyr', quietly = TRUE)) install.packages('dplyr'); library('dplyr')
if (!requireNamespace('rtracklayer', quietly = TRUE)) BiocManager::install('rtracklayer'); library('dplyr')
if (!requireNamespace('stringr', quietly = TRUE)) install.packages('stringr'); library('stringr')

#library(polyester);library(ggrepel)

# NA_peaks_neg <- read.table("~/dnx_projects/EXP23003664_aged_HDF/peakCountBySample_broadPeaks/AG06291_rep1_S13.markdup.sorted.neg/MACS_RNA/NA_peaks.xls", header=TRUE, quote="\"")
# NA_peaks_pos <- read.table("~/dnx_projects/EXP23003664_aged_HDF/peakCountBySample_broadPeaks/AG06291_rep1_S13.markdup.sorted.plus/MACS_RNA/NA_peaks.xls", header=TRUE, quote="\"")

# read in NA_peaks_neg and NA_peaks_pos as 2 inputs from command line
NA_peaks_neg <- read.table(commandArgs(TRUE)[1], header=TRUE, quote="\"")
NA_peaks_pos <- read.table(commandArgs(TRUE)[2], header=TRUE, quote="\"")


NA_peaks_neg$strand <- "-"
NA_peaks_pos$strand <- "+"

# filter for qValue above 2
NA_peaks_neg <- NA_peaks_neg[NA_peaks_neg$X.log10.qvalue. > 2,]
NA_peaks_pos <- NA_peaks_pos[NA_peaks_pos$X.log10.qvalue. > 2,]

# zoom in on x axis without removing data
#ggplot() +
#  geom_histogram(data=NA_peaks_neg, aes(x=length), bins=500, fill="blue", alpha=0.5) +
#  geom_histogram(data=NA_peaks_pos, aes(x=length), bins=500, fill="red", alpha=0.5) +
#  theme_minimal() +
#  theme(legend.position="none") +
#  labs(title="Distribution of peak widths in NA peaks", x="Peak width", y="Count") +
#  coord_cartesian(xlim=c(0, 2000))

# load metadata
chromAnnoMeta <- read.csv("/weka/mwang/dnx_projects/full_stack_ChromHMM_annotations/state_annotations_processed.csv.gz")
rownames(chromAnnoMeta) <- chromAnnoMeta$mneumonics

library(readr)
hg38lift_genome_100_browser <- data.frame(read_table("/weka/mwang/dnx_projects/full_stack_ChromHMM_annotations/hg38lift_genome_100_browser.bed.gz"))

hg38lift_genome_100_browser$length<-unlist(hg38lift_genome_100_browser[,3]-hg38lift_genome_100_browser[,2])
hg38lift_genome_100_browser$group<-chromAnnoMeta[as.character(hg38lift_genome_100_browser[,4]),"Group"]
nrow(hg38lift_genome_100_browser)
# remove "chr" in chromosome names
hg38lift_genome_100_browser$track<-gsub("chr","",hg38lift_genome_100_browser$track)
colnames(hg38lift_genome_100_browser)<-c("chrom","start","end","name","strand","length","group")
hg38lift_genome_100_browser$strand<-"."

overlapInfo<-findOverlaps(query = GRanges(NA_peaks_neg), subject = GRanges(hg38lift_genome_100_browser))
overlapInfo<-as.data.frame(overlapInfo)
overlapInfo$name<-hg38lift_genome_100_browser$name[overlapInfo$subjectHits]
overlapInfo$group<-hg38lift_genome_100_browser$group[overlapInfo$subjectHits]
# for each queryHits, collect the set of names and groups
overlapInfo2<-overlapInfo %>% group_by(queryHits) %>% summarize(name=paste(name,collapse=","),group=paste(group,collapse=",")) %>% data.frame()
rownames(overlapInfo2)<-overlapInfo2$queryHits
# make overlapInfo2 name and group unique
overlapInfo2$name<-sapply(overlapInfo2$name,function(x) paste(unique(unlist(strsplit(x,","))),collapse=","))
overlapInfo2$group<-sapply(overlapInfo2$group,function(x) paste(unique(unlist(strsplit(x,","))),collapse=","))
NA_peaks_neg$chromAnno_name<-"uannotated"
NA_peaks_neg$chromAnno_name[overlapInfo2$queryHits]<-overlapInfo2$name
NA_peaks_neg$chromAnno_group<-"uannotated"
NA_peaks_neg$chromAnno_group[overlapInfo2$queryHits]<-overlapInfo2$group

# repeat for NA_peaks_pos
overlapInfo<-findOverlaps(query = GRanges(NA_peaks_pos), subject = GRanges(hg38lift_genome_100_browser))
overlapInfo<-as.data.frame(overlapInfo)
overlapInfo$name<-hg38lift_genome_100_browser$name[overlapInfo$subjectHits]
overlapInfo$group<-hg38lift_genome_100_browser$group[overlapInfo$subjectHits]

# for each queryHits, collect the set of names and groups
overlapInfo2<-overlapInfo %>% group_by(queryHits) %>% summarize(name=paste(name,collapse=","),group=paste(group,collapse=",")) %>% data.frame()
rownames(overlapInfo2)<-overlapInfo2$queryHits
# make overlapInfo2 name and group unique
overlapInfo2$name<-sapply(overlapInfo2$name,function(x) paste(unique(unlist(strsplit(x,","))),collapse=","))
overlapInfo2$group<-sapply(overlapInfo2$group,function(x) paste(unique(unlist(strsplit(x,","))),collapse=","))
NA_peaks_pos$chromAnno_name<-"uannotated"
NA_peaks_pos$chromAnno_name[overlapInfo2$queryHits]<-overlapInfo2$name
NA_peaks_pos$chromAnno_group<-"uannotated"
NA_peaks_pos$chromAnno_group[overlapInfo2$queryHits]<-overlapInfo2$group

# combine pos and neg
NA_peaks<-rbind(NA_peaks_neg,NA_peaks_pos)
NA_peaks_GR<-GRanges(NA_peaks)
NA_peaks_GR$gene_id<-paste(NA_peaks_GR$name,NA_peaks$strand,NA_peaks_GR$chromAnno_name,NA_peaks_GR$chromAnno_group,sep="_")
NA_peaks_GR$gene_name<-NA_peaks_GR$gene_id

# export as gtf
# rtracklayer::export(NA_peaks_GR, "~/dnx_projects/EXP23003664_aged_HDF/peakCountBySample_broadPeaks_gtf/AG06291_rep1_S13.gtf", format = "gtf")

# export as gtf based on command line input 3, only gene_name and gene_id as attributes
# 



# only export attributes gene_name and gene_id
# rtracklayer::export(NA_peaks_GR, "~/dnx_projects/EXP23003664_aged_HDF/peakCountBySample_broadPeaks_gtf/AG06291_rep1_S13.gtf", format = "gtf", attributes = c("gene_name", "gene_id"))

# subset NA_peaks_GR for only gene_name and gene_id attributes
NA_peaks_GR_subset <- NA_peaks_GR[,c("gene_name", "gene_id")]

rtracklayer::export(NA_peaks_GR_subset, commandArgs(TRUE)[3], format = "gtf")
