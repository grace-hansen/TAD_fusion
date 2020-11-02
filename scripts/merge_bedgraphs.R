#!/usr/bin/Rscript
args=commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("USAGE: merge_bedgraphs.R <sample_list> <dir> \n", call.=FALSE)
}
library(data.table)

samples<-scan(args[1],what='character',sep='\n')
dir=args[2]
setwd(dir)
############ Author: Grace Hansen ##############
# This script merges different bedgraphs and writes out a bedgraph representing the merged value from the bedgraphs.

input<-paste(paste(samples,".DI.bedGraph",sep=''),collapse=' ')
all_DI<-system(paste("bedtools unionbedg -i ",input,"",sep=''),intern=TRUE)
all_DI<-as.data.frame(do.call('rbind',strsplit(as.data.frame(all_DI,stringsAsFactors=FALSE)$all_DI,'\t',fixed=TRUE)),stringsAsFactors=FALSE)
all_DI[,4:ncol(all_DI)] <- sapply(all_DI[,4:ncol(all_DI)], as.numeric)


means<-rowMeans(all_DI[,4:ncol(all_DI)])
mean_DI<-cbind(all_DI[,1:3],means)

write.table(mean_DI,"merged.DI.bedGraph",col.names = FALSE,row.names=FALSE,quote=FALSE,sep='\t')