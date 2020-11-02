#!/bin/bash

########### Author: Grace Hansen #########3
# This script is a wrapper for RestrictionDigest
RE=$1
fasta=$2
outdir=$3
logdir=$4

mkdir ${outdir}${RE}_digest
perl scripts/digest_genome.pl -f $fasta -r $RE -o ${outdir}${RE}_digest
cd ${outdir}${RE}_digest

fastaname=$(echo $fasta | rev | cut -f1 -d'/' | rev)
mv digestion_summary_${fastaname}_by_${RE} coverage_every_chromosome_${fastaname}_by_${RE} $logdir
rm seq_frags_*${fastaname}_by_${RE}

cut -f1 -d' ' position_frags_${fastaname}_by_${RE}  | sed "s|>||g" > chr
cut -f2 -d' ' position_frags_${fastaname}_by_${RE} | cut -f1 > ind
cut -f2-3 position_frags_${fastaname}_by_${RE} > pos
paste chr pos ind > ${RE}_frags.bed

cut -f1 -d' ' position_frags_in_range_${fastaname}_by_${RE} | sed "s|>||g" > chr
cut -f2 -d' ' position_frags_in_range_${fastaname}_by_${RE} | cut -f1 > ind
cut -f2-3 position_frags_in_range_${fastaname}_by_${RE} > pos
paste chr pos ind > ${RE}_frags_sizefilt.bed

rm chr pos ind
rm position_frags_*${fastaname}_by_${RE}