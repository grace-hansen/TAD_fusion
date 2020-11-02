#!/usr/bin/Rscript
library(optparse)
library(tidyverse)
library(data.table)
arguments <- parse_args(OptionParser(usage = "filter_TADs_DIs.R <tad_bed> <DI bedgraph> <RE> <N> <dir>",option_list=list()),
                        positional_arguments = 5)
opt<-arguments$opt
tads<-arguments$args[1]
DI<-arguments$args[2]
RE<-arguments$args[3]
N<-arguments$args[4]
setwd(arguments$args[5])

############ Author: Grace Hansen ##############
# This script takes as input the filtered merged.tad.2d.bed file
# It filters this file: 
    # Remove tads that overlap with other tads (shouldn't be any, but is a few (need to look up why this is on HOMER))
    # Remove tads < 100kb (average of 200 frags in a TAD of this size, so 25 tads/side is reasonable)
    # Keep N tads with highest score (can't go too much over 2000, only start with 3500 TADs)
# It merges these TADs with DI information from merged.DI.bedGraph data
# It ranks the TADs by abs(DI), and outputs the top N TADs by DI

#Filter TADs:
## Remove TADs that overlap other TADs
sorted_tads=gsub("2D","sorted",tads)
merged_tads<-system(paste("bedtools merge -i ",sorted_tads,".gz",sep=''),intern=TRUE)
merged_tads<-data.frame(do.call('rbind',strsplit(as.data.frame(merged_tads,stringsAsFactors=FALSE)$merged_tads,'\t',fixed=TRUE)),stringsAsFactors = FALSE)
colnames(merged_tads)<-c("chr","start","stop")
merged_tads$start<-as.integer(merged_tads$start)
merged_tads$stop<-as.integer(merged_tads$stop)
merged_tads$chr<-sapply(strsplit(merged_tads$chr,'r'),'[[',2)

tads<-as.data.frame(fread(tads),stringsAsFactors=FALSE)
colnames(tads)<-c("chr","start","stop","chr1","start1","stop1","color","info1","info2")
tads$score<-rowMeans(tads[,8:ncol(tads)])
tads<-merge(tads,merged_tads,by=c("chr","start","stop")) 

## Remove TADs < 100kb
tads<-tads[(tads$stop-tads$start)>=100000,]

## Keep top N TADs by score 
tads<-tads[order(tads$score,decreasing = TRUE),]
tads<-tads[1:N,]
write.table(tads,paste("tads_filtered_",N,".bed",sep=''),quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)

#Get DI for fragments that overlap TADs
print("Getting average DI per restriction enzyme in TADs...")
system(paste("sed '1d' tads_filtered_",N,".bed | awk '$1=\"chr\"$1' | sort -k1,1 -k2,2n | tr ' ' '\t' > tads_",N,".sorted.bed",sep=''))
system(paste("bedtools intersect -wa -a ",RE,"_digest/",RE,"_frags_sizefilt.bed -b tads_",N,".sorted.bed | sort -u > ",RE,"_frags_TADs_",N,".bed",sep=''))#subset RE frags to only those overapping filtered TADs
cmd=paste("multiBigwigSummary BED-file -b ",DI," --BED ",RE,"_frags_TADs_",N,".bed --outRawCounts ",RE,"_frags_TADs_",N,"_avgDI.txt -o ",RE,"_frags_filtTADs_avgDI.npz",sep='')
print(cmd)
system(cmd)
system(paste("rm tads_",N,".sorted.bed ",RE,"_frags_filtTADs_avgDI.npz",sep=''))

#Grab top 25 pos and top 25 neg DI fragments within each of the top N TADs
print(paste("Writing out fragments with highest abs(DI) in top ",N," TADs with highest abs(DI)",sep=''))
frags<-fread(paste(RE,"_frags_TADs_",N,"_avgDI.txt",sep=''))
colnames(frags)<-c('chr','start','end','DI')

RE_frags<-matrix(ncol=5)
for (i in 1:nrow(tads)) {
  tad_frags<-frags[frags$chr==paste("chr",tads$chr[i],sep='') & frags$start>tads$start[i] & frags$end< tads$stop[i]] #Grab fragments within TAD
  top_neg_frags<-tad_frags[order(tad_frags$DI)][1:25,]
  top_neg_frags<-top_neg_frags[top_neg_frags$DI<0,] #Don't keep if not negative
  top_pos_frags<-tad_frags[order(tad_frags$DI,decreasing=TRUE)][1:25,]
  top_pos_frags<-top_pos_frags[top_pos_frags$DI>0,] #Don't keep if not positive
  top_frags<-rbind(top_pos_frags,top_neg_frags)
  top_frags$index<-rep(paste(tads$chr[i],tads$start[i],tads$stop[i],sep='-'),nrow(top_frags))
  RE_frags<-rbind(RE_frags,top_frags,use.names=FALSE)
}
RE_frags<-RE_frags[2:nrow(RE_frags),]
ind_counts<-as.data.frame(table(RE_frags$V5))
full_TADs<-ind_counts$Var1[ind_counts$Freq==50]
RE_frags<-RE_frags[RE_frags$V5 %in% full_TADs,] #Remove probes from TADs with fewer than 50 probes
write.table(RE_frags,paste("top_frags_selected_tads.bed",sep=''),row.names=FALSE,col.names=FALSE,sep='\t',quote=FALSE)

#Write out TADs to match probes after filtering
tads$index<-paste(tads$chr,tads$start,tads$stop,sep='-')
tads<-tads[tads$index %in% RE_frags$V5,]
write.table(tads,"selected_tads.bed",quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
