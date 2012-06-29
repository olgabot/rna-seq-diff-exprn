#!/bin/sh

####################### !!! ##########################
# Make sure you add these to your path before running:
# export PYTHONPATH=/share/apps/RSeQC/usr/local/lib/python2.7/site-packages:$PYTHONPATH
# export PATH=/share/apps/RSeQC/usr/local/bin:$PATH
####################### !!! ##########################

# Example run:
# rseqc.sh my_favorite_sample.bam rseqc /share/reference/hg19/genes/hg19_ucsc-genes.bed /share/reference/hg19/genes/hg19_ucsc-genes.gtf
##
# Example using bash scripting:
# rseqc.sh $BAM $OUT_DIR $BED

# BAM file that you want to do the analysis on
BAM=$1

# Output directory where you want all the results to go
# Files will be prepended with "$OUT_DIR/"
OUT_DIR=$2

# BED files of which you want to see the gene body coverage
BED=$3

# Make plots via PREFIX
# BAM=tophat_6n6_trimmed_34_merged_rmdup.bam ; REF_BED=/home/obot/bed_and_gtf/hg19_ucsc-genes.bed ; SINO_BED=/home/obot/bed_and_gtf/hg19_ucsc-genes.bed ; PREFIX=PREFIX

PREFIX=$OUT_DIR/

echo 'clipping_profile.py' ; clipping_profile.py --input-file $BAM --out-prefix $PREFIX
echo 'geneBody_coverage.py' ; geneBody_coverage.py --input-file $BAM --refgene $BED --out-prefix $PREFIX
echo 'inner_distance.py' ; inner_distance.py --input-file $BAM --refgene $BED --out-prefix $PREFIX
echo 'junction_annotation.py' ; junction_annotation.py --input-file $BAM --refgene $BED --out-prefix $PREFIX
echo 'junction_saturation.py' ; junction_saturation.py --input-file $BAM --refgene $BED --out-prefix $PREFIX
echo 'read_duplication.py' ; read_duplication.py --input-file $BAM --out-prefix $PREFIX
echo 'read_distribution.py' ; read_distribution.py --input-file $BAM --refgene $BED >$PREFIX\_read_distribution.txt
echo 'read_GC.py' ; read_GC.py --input-file $BAM --out-prefix $PREFIX
echo 'read_NVC.py' ; read_NVC.py --input-file $BAM --out-prefix $PREFIX
echo 'read_quality.py' ; read_quality.py --input-file $BAM --out-prefix $PREFIX
echo 'RPKM_saturation.py' ; RPKM_saturation.py --input-file $BAM --out-prefix $PREFIX --refgene $BED 
echo 'RPKM_count.py' ; RPKM_count.py --input-file $BAM --out-prefix $PREFIX --refgene $BED