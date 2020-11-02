#!/bin/bash

########### Author: Grace Hansen #########
# This script merges per-chromosome output from HOMER findTADsandloops.pl
#Files to merge:
  # *.DI.bedGraph
  # *.interactions.txt
  # *.tad.2D.bed
  # *.ucsc.bed
  # *.Insulation.bedGraph
  # *.loop.2D.bed
  # *.tad.peaks.txt

prefix=$1
outdir=$2

cd $outdir

echo -e "track name=\"${outdir}/${prefix} DI (/scratch/midway2/gthansen/${prefix})\" type=bedGraph visibility=2" > ${outdir}/${prefix}.DI.bedGraph
echo -e "#InteractionID\tPeakID(1)\tchr(1)\tstart(1)\tend(1)\tstrand(1)\tTotal Reads(1)\tPeakID(2)\tchr(2)\tstart(2)\tend(2)\tstrand(2)\tTotal Reads(2)\tDistance\tInteraction Reads\tExpected Reads\tZ-score\tLogP\tFDR" > ${outdir}/${prefix}.interactions.txt
> ${outdir}/${prefix}.tad.2D.bed
> ${outdir}/${prefix}.ucsc.bed
echo -e "track name=\"${outdir}/${prefix} Insulation Ratio (/scratch/midway2/gthansen/${prefix})\" type=bedGraph visibility=2" > ${outdir}/${prefix}.Insulation.bedGraph
> ${outdir}/${prefix}.loop.2D.bed
> ${outdir}/${prefix}.tad.peaks.txt

for i in $(seq 1 22) X Y; do
	sed '1D' ${outdir}/${prefix}_$i.DI.bedGraph >> ${prefix}.DI.bedGraph
	sed '1D' ${outdir}/${prefix}_$i.interactions.txt >> ${prefix}.interactions.txt
	cat ${outdir}/${prefix}_$i.tad.2D.bed >> ${outdir}/${prefix}.tad.2D.bed
	sed '1D' ${outdir}/${prefix}_$i.ucsc.bed >> ${prefix}.ucsc.bed
	sed '1D' ${outdir}/${prefix}_$i.Insulation.bedGraph >> ${outdir}/${prefix}.Insulation.bedGraph
	cat ${outdir}/${prefix}_$i.loop.2D.bed >> ${prefix}.loop.2D.bed
	cat ${outdir}/${prefix}_$i.tad.peaks.txt >> ${outdir}/${prefix}.tad.peaks.txt
done

# rm ${outdir}/${prefix}_*.DI.bedGraph
# rm ${outdir}/${prefix}_*.interactions.txt
# rm ${outdir}/${prefix}_*.tad.2D.bed
# rm ${outdir}/${prefix}_*.ucsc.bed
# rm ${outdir}/${prefix}_*.Insulation.bedGraph
# rm ${outdir}/${prefix}_*.loop.2D.bed
# rm ${outdir}/${prefix}_*.tad.peaks.txt