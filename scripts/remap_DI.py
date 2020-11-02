#!/usr/bin/python
import subprocess, pandas as pd, argparse, numpy as np
parser = argparse.ArgumentParser()
parser.add_argument("DI", help="bedGraph representing mean DI (directionality index)")
parser.add_argument("frags", help="bed file with restriction fragments from whole genome")
parser.add_argument("chr", help="chr to run, in format '1")
parser.add_argument("dir", help="directory in which files are located, and which to write to")

args = parser.parse_args()

################## Author: Grace Hansen ########################
# This  script takes directionality index (DI) values from a bedGraph file and re-maps them to
# a bed file consisting of restriction enzyme fragments.
# the value per fragment represents the mean value over that area. 
# To do this, it calculates the per-base DI, then averages this over the bases in the fragment.
# As such, it makes big files and takes a lot of time.
###################################################################
DI=args.DI
frags=args.frags
chr=args.chr
dir=args.dir
RE=DI.split('/')[-1].split('_')[0]
buffer=100

def extract_conservation(frags,DI,chr):
	with open(frags,'r') as frags:
			for ln in frags:
				chrom,start,stop,ind=ln.strip().split('\t')[0:4]
				if chrom=='chr'+chr:
					print(chr)
					perbase_out=open("tempfrag_%s.bed"%chr,'w')
					for bp in range(int(start),int(stop)):
						perbase_out.write(chr+'\t'+str(bp)+'\t'+str(bp+1)+'\t'+str(ind)+'\n')
					perbase_out.close()
					winstart=int(start)-buffer
					winend=int(stop)+buffer
					subprocess.check_call("bedtools intersect -wb -a %s -b temp_%s_bed.bed > DItemp.bed"%(DI,chr),shell=True)
					subprocess.check_call("""bedmap --count --echo --echo-map-score --delim '\t' tempfrag_%s.bed DItemp.bed | awk '{ if ($1=="0") { print $0"0.0"; } else { print $0; } }' | cut -f2- >> %s_%s.bed"""%(chr,RE,chr),shell=True)
					subprocess.check_call("rm tempfrag_%s.bed DItemp.bed"%chr,shell=True)

extract_conservation(frags,DI,chr)