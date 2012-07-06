#!/bin/sh

# 
# Created by:  on 
# Purpose of this script: Makes chr9-only bam files from the bam files in the current directory
# Input(s): 
# Output(s):


for F in `ls *.bam | grep -P '_\d\.bam' - | tr '\n' ' '`; do
	F_PREFIX=`echo $F | tr -d '.bam'`
	CHR9=$F_PREFIX\_chr9.bam
	samtools index $F
	samtools view -b $F chr9 > $CHR9
	samtools view -H $F | \
		samtools reheader - $CHR9 > $F_PREFIX\_chr9_withheader.bam
done