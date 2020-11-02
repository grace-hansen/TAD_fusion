#!/bin/sh

file=$1
logs=$2
prefix=$(echo $file | cut -f1 -d '.')
sample=$(echo $prefix | rev | cut -f1 -d'/' | rev)

/home/gthansen/bin/homer/bin/makeTagDirectory $prefix -format HiCsummary $file > ${logs}/${sample}_makeTAGdir.log

#Calculate DI and identify TADs
/home/gthansen/bin/homer/bin/findTADsAndLoops.pl find $prefix -cpu 16 -res 10000 -window 15000 -genome hg19 > ${logs}/${sample}_calcDI.log