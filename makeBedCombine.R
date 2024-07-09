# read in input from command line, should be 2 input bed.plus and bed.neg
bed.plus<-read.table(commandArgs(TRUE)[1],header=F,stringsAsFactors=F)
bed.neg<-read.table(commandArgs(TRUE)[2],header=F,stringsAsFactors=F)

colnames(bed.plus)<-c("chr","start","end","score")
bed.plus$name<-"."
bed.plus$strand<-"+"
bed.plus<-bed.plus[,c("chr","start","end","name","score","strand")]

# do the same for bed.neg
colnames(bed.neg)<-c("chr","start","end","score")
bed.neg$name<-"."
bed.neg$strand<-"-"
bed.neg<-bed.neg[,c("chr","start","end","name","score","strand")]

#write out as bed file
write.table(bed.plus,file=commandArgs(TRUE)[3],quote=F,row.names=F,col.names=F,sep="\t")
write.table(bed.neg,file=commandArgs(TRUE)[4],quote=F,row.names=F,col.names=F,sep="\t")
