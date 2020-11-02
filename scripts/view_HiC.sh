#!/bin/bash
sample=$1
dir=$2

cd $dir
awk '$1="chr"$1' $sample.tad.2D.bed | sort -k1,1 -k2,2n | cut -f1-3 -d' ' | tr ' ' '\t' > $sample.tad.sorted.bed
bgzip $sample.tad.sorted.bed
tabix -p bed $sample.tad.sorted.bed.gz

awk '$1="chr"$1' $sample.loop.2D.bed | sort -k1,1 -k2,2n | cut -f1-3 -d' ' | tr ' ' '\t' > $sample.loop.sorted.bed
bgzip $sample.loop.sorted.bed
tabix -p bed $sample.loop.sorted.bed.gz

awk '$1="chr"$1' <(sed '1d' $sample.DI.bedGraph) | sort -k1,1 -k2,2n | tr ' ' '\t' > $sample.DI.sorted.bedGraph
/project2/nobrega/grace/bin/bedGraphToBigWig $sample.DI.sorted.bedGraph /project2/nobrega/grace/genos/hg19/hg19.chrom.sizes $sample.DI.bw
rm $sample.DI.sorted.bedGraph

sed '1d' $sample.Insulation.bedGraph | awk '$1="chr"$1' | sort -k1,1 -k2,2n | tr ' ' '\t' > $sample.Insulation.sorted.bedGraph
/project2/nobrega/grace/bin/bedGraphToBigWig $sample.Insulation.sorted.bedGraph /project2/nobrega/grace/genos/hg19/hg19.chrom.sizes $sample.Insulation.bw
rm $sample.Insulation.sorted.bedGraph