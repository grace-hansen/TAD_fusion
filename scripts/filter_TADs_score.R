#!/usr/bin/Rscript
library(optparse)
library(tidyverse)
library(data.table)
arguments <- parse_args(OptionParser(usage = "filter_TADs_score.R <tad_bed> <N>",option_list=list()),
                        positional_arguments = 2)
opt<-arguments$opt
tads<-arguments$args[1]
outdir=paste(strsplit(tads,'/')[[1]][1:length(strsplit(tads,'/')[[1]])-1],collapse='/')
N<-arguments$args[2]
setwd("~/midway/TAD_fusion")

############ Author: Grace Hansen ##############
# This script takes as input the merged TAD bed file and a number between 0 and 1
# It then keeps the N tads with highest mean score
tads<-as.data.frame(fread(tads),stringsAsFactors=FALSE)
tads$score<-rowMeans(tads[,8:ncol(tads)])
tads<-tads[order(tads$score,decreasing = TRUE),]
tads<-tads[1:N,]
write.table(tads,paste(outdir,"/tads_filtered_",N,".bed",sep=''),quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
