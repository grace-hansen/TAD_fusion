#!/usr/bin/Rscript
library(optparse)
library(tidyverse)
library(data.table)
arguments <- parse_args(OptionParser(usage = "get_probes_from_frags.R <frags> <probe length> <fasta> <out>",option_list=list()),
                        positional_arguments = 4)
opt<-arguments$opt
frags<-arguments$args[1]
len<-as.numeric(arguments$args[2])
fasta<-arguments$args[3]
setwd(arguments$args[4])

############ Author: Grace Hansen ##############
# This script takes as input top_frags_selected_tads.bed, and selects probes based on it.
# It grabs sequences of X base pairs as potential probes, and performs the below filtering steps:
  # Remove probes with <25% or >65% GC content
  # Remove probes that don't map uniquely to genome
##############################################

#Get sequences for frags
system(paste("bedtools getfasta -fi ",fasta," -bed ",frags," -tab -fo top_frags_seqs",sep=''))
seqs<-fread("top_frags_seqs",header=FALSE)
colnames(seqs)<-c("fragloc","seq")

frags<-as.data.frame(fread(frags),stringAsFactors=FALSE)
colnames(frags)<-c("chr","start","stop","score","tad")
frags$fragloc<-paste(frags$chr,paste(frags$start,frags$stop,sep='-'),sep=':')
frags<-merge(frags,seqs,by="fragloc")
seqs<-NULL
system("rm top_frags_seqs")

#For each frag sequence, create 3 candidate probes from beginning, middle, and end
frags$probe1<-str_sub(frags$seq,1,len)
frags$probe2<-str_sub(frags$seq,-(len),-1)

frags$loc1<-paste(frags$chr,frags$start+1,frags$start+len+1,sep='-')
frags$loc2<-paste(frags$chr,frags$stop-len+1,frags$stop+1,sep='-')#Make sure these match

frags$GC1<-lengths(regmatches(frags$probe1,gregexpr("G|C",frags$probe1)))/len
frags$GC2<-lengths(regmatches(frags$probe2,gregexpr("G|C",frags$probe2)))/len

#Overlap probes from "" with rmsk, returning probes with <60% overlap.
#bed1<-separate(frags,col=loc1, into = c("chr","start","stop"), sep = "\\-") %>% select(chr,start,stop,loc)
#bed2<-separate(frags,col=loc2, into = c("chr","start","stop"), sep = "\\-") %>% select(chr,start,stop,loc)
#write.table(bed1,"probe_frags_1.bed",sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
#write.table(bed2,"probe_frags_2.bed",sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)
#system(paste("bedtools intersect -a probe_frags_1.bed -b ../rmsk_hg19.bed -f 0.6 -v > probe_frags_norepeats_1.bed",sep=''))
#system(paste("bedtools intersect -a probe_frags_2.bed -b ../rmsk_hg19.bed -f 0.6 -v > probe_frags_norepeats_2.bed",sep=''))

#For each frag sequence, save out as fasta. Map to genome with no multimapping allowed. Load mapped reads back in, and merge with frags dataframe
get_aligned_reads<-function(probes,i) {
  colnames(probes)<-c("tad","seq","loc")
  fa<-as.data.table(rep(">",2*nrow(probes)))
  fa[seq(2,nrow(fa),2)]<-probes$seq
  write.table(fa,paste("probes",i,".fa",sep=''),quote=FALSE,row.names=FALSE,col.names=FALSE,sep='\t')
  system(paste(" /home/gthansen/STAR-2.7.3a/source/STAR --readFilesIn probes",i,".fa --outFileNamePrefix probes",i," --parametersFiles ../STAR_params",sep=''))
  aligned=fread(paste("grep -v '@' probes",i,"Aligned.out.sam | cut -f3,4,6,10",sep=''))
  colnames(aligned)<-c("chr","start","mapped","seq")
  return(aligned)
}

probes1<-frags[,c("tad","probe1","loc1")]
probes2<-frags[,c("tad","probe2","loc2")]
aligned1<-get_aligned_reads(probes1,1)
aligned1$mapped1<-aligned1$mapped
aligned1$probe1<-aligned1$seq
aligned2<-get_aligned_reads(probes2,2)
aligned2$mapped2<-aligned2$mapped
aligned2$probe2<-aligned2$seq

frags<-merge(frags,aligned1[,c("mapped1","probe1")],all.x=TRUE)
frags<-merge(frags,aligned2[,c("mapped2","probe2")],all.x=TRUE)
frags$mapped1[is.na(frags$mapped1)]<-"0"
frags$mapped2[is.na(frags$mapped2)]<-"0"

#For each fragment, find a probe that works, if present
frags$probeseq=NA
frags$probeloc=NA
frags$probefragloc=NA
frags$probeGC=NA
for (i in 1:nrow(frags)) {
  if (frags$GC1[i]>0.25 & frags$GC1[i]<0.6) {
    if (frags$mapped1[i]=="120M") {
      frags$probefragloc[i]="beginning"
      frags$probeseq[i]=frags$probe1[i]
      frags$probeloc[i]=frags$loc1[i]
      frags$probeGC[i]=frags$GC1[i]
    }
  } else if (frags$GC2[i]>0.25 & frags$GC2[i]<0.6) {
    if (frags$mapped2[i]=="120M") {
      frags$probefragloc[i]="end"
      frags$probeseq[i]=frags$probe2[i]
      frags$probeloc[i]=frags$loc2[i]
      frags$probeGC[i]=frags$GC2[i]
    }
  }
}

####### Keep only probes in TADs where there are >5 probes on each side of TAD ##########
frags$probechr<-sapply(strsplit(frags$probeloc,'-'),'[[',1)
frags$probestart[!(is.na(frags$probeloc))]<-sapply(strsplit(frags$probeloc[!(is.na(frags$probeloc))],'-'),'[[',2)
frags$probestop[!(is.na(frags$probeloc))]<-sapply(strsplit(frags$probeloc[!(is.na(frags$probeloc))],'-'),'[[',3)
probes<-frags[!(is.na(frags$probeloc)),]
write.table(probes[,c("probechr","probestart","probestop")],"probes_temp.bed",row.names=FALSE,col.names=FALSE,sep='\t',quote=FALSE)

#Get avg DI per probe
system("sort -k1,1 -k2,2 probes_temp.bed > temp && mv temp probes_temp.bed")
cmd=paste("multiBigwigSummary BED-file -b merged_HIC.DI.bw --BED probes_temp.bed --outRawCounts probes_avgDI.txt -o probes_avgDI.npz",sep='')
print(cmd)
system(cmd)
probeDI<-fread("probes_avgDI.txt")
colnames(probeDI)<-c("chr","start","end","DI")
probeDI$probeloc<-paste(probeDI$chr,probeDI$start,probeDI$end,sep='-')
frags<-merge(frags,probeDI[,c("probeloc","DI")],all.x=TRUE)
write.table(frags,"candidate_probes.txt",quote=FALSE,sep='\t',row.names=FALSE)
system(paste("rm tads_sorted.bed probes_avgDI.txt probes_avgDI.npz probes_temp.bed",sep=''))

#Get TADs with > 10 probes per side (i.e. 10 probes with pos DI, 10 probes with neg DI)
N=5
probes<-frags[!(is.na(frags$probeloc)),c("probeloc","probeseq","probefragloc","probeGC","tad","DI")]
probes<-as.data.table(probes)
filtered_probes<-matrix(ncol=6)
for (i in unique(probes$tad)) {
  neg_probes<-probes[probes$tad==i & probes$DI<0,]
  pos_probes<-probes[probes$tad==i & probes$DI<0,]
  if (nrow(neg_probes)>=N && nrow(pos_probes)>=N) {
    filtered_probes<-rbind(filtered_probes,neg_probes,pos_probes,use.names=FALSE)
  }
}

#Re-subset probes and tads
filtered_probes<-filtered_probes[2:nrow(filtered_probes),]
colnames(filtered_probes)<-c("loc","seq","probefragloc","GC","tad","DI")
write.table(filtered_probes,paste("selected_probes.txt",sep=''),row.names=FALSE,col.names=FALSE,sep='\t',quote=FALSE)

#Write out TADs to match frags after filtering
tads<-fread("filtered_tads.bed")
selected_tads<-tads[tads$V11 %in% filtered_probes$tad]
write.table(selected_tads,"selected_tads.bed",quote=FALSE,sep='\t',row.names=FALSE,col.names=FALSE)
system("rm probes1* probes2*")
