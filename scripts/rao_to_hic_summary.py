#! /usr/bin/python

import argparse, io, gzip
parser = argparse.ArgumentParser()
parser.add_argument("file", help="path of file being converted")
parser.add_argument("out", help="directory to put output files")
args = parser.parse_args()

######### Author: Grace Hansen #########
# This  script takes data downloaded from the Rao et al paper () and converts it into a format that HOMER can process.

###################################################################
file=args.file
out=args.out
prefix=file.split('/')[-1].split('.')[0]

#Make list of all variants in GWAS summary statistics to filter by
file=io.TextIOWrapper(gzip.open(file,'r'))
out=gzip.open("%s/%s.hicsummary.gz"%(out,prefix),'wt')

for ln in file:
	read,strand1,chr1,pos1,index1,strand2,chr2,pos2,index2,QC1,QC2=ln.strip().split()
	if strand1=='0':
		newstrand1='+'
	elif strand1=='16':
		newstrand1='-'
	if strand2=='16':
		newstrand2='+'
	elif strand2=='0':
		newstrand2='-'
	out.write('\t'.join([read,chr1,pos1,newstrand1,chr2,pos2,newstrand2])+'\n')
