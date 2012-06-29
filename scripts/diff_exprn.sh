#!/usr/bin/sh

BAM=$1
GFF=$2
BAM_PREFIX=$3

PATH=$PATH:/share/apps/samtools-0.1.18

# Gene Expression Estimation via HTSeq
samtools view | python2.7 htseq-count --stranded=no - $GFF 

# Gene Expression Estimation via bedtools coverage
GENE_COUNTS=$BAM_PREFIX\_gene_counts_bedtools.txt
/share/apps/BEDTools-Version-2.15.0/bin/bedtools coverage -abam $SINO_BAM -b /home/obot/hg19_ensembl-genes_s.bed > $GENE_COUNTS
