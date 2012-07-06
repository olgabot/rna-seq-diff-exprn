#!/bin/sh -x

# Author: Olga Botvinnik (olga.botvinnik@gmail.com)
# Date: 21 June 2012
# Purpose of this script: Finds gene- and genome-wide counts
# of a BED file using both bedtools and HTSeq
#
# This gene_counts.sh script does:
# 1. Estimate gene counts from BAM via bedtools and HTSeq
# 2. Calls another script to:
#    a. Estimate genome-wide coverage via bedtools and HTSeq
#    b. Make circos plot of coverage for this sample
# 3. Estimates differential exon usage counts using
#    dexseq_counts.py (does not perform DEXSeq analysis)

BAM_PREFIX=$1
# OUT_DIR=$2
GENDER=$2

# Need the ID for circos file creation/plotting
ID=$3
COMMON_VARS=$4

# Initialize common variables
source $COMMON_VARS

BAM=$BAM_PREFIX.bam
echo "finding gene counts for" $BAM

######## Begin gene count estimation via Bedtools #########
# Do the intersect
BAM_INTERSECT=$BAM_PREFIX\_intersect.bam
samtools view -b $BAM | \
    bedtools intersect -abam - \
    -b $BED > $BAM_INTERSECT

# Gene expression estimation via bedtools
# -counts: Only report the count of overlaps, 
#          don't compute fraction, etc.
# -hist: Report a histogram of coverage for each feature in B
#        as well as a summary histogram for _all_ 
#        features in B.
#        Output (tab delimited) after each feature in B:
#           1) depth
#           2) # bases at depth
#           3) size of B
#           4) % of B at depth
# -d: Report the depth at each position in each B feature.
#     Positions reported are one based.  Each position
#     and depth follow the complete B feature.


COUNTS_BED=bedtools_gene_counts.txt
echo 'COUNTS_BED=$COUNTS_BED' | cat - >> $COMMON_VARS
THIS_COUNTS_BED=$OUT_DIR/$COUNTS_BED
coverageBed_options='-counts -hist -d'

bedtools coverage \
    $coverageBed_options \
    -abam $BAM_INTERSECT \
    -b $BED >$THIS_COUNTS_BED
######## END gene count estimation via Bedtools #########


######## BEGIN gene count estimation via HTSeq #########
############### For some reason this does NOT work when
############### executed from this script but works fine 
############### from the command line

COUNTS_HTSEQ_PREFIX=htseq_counts
COUNTS_HTSEQ=$COUNTS_HTSEQ_PREFIX.txt
echo 'COUNTS_HTSEQ=$COUNTS_HTSEQ' | cat - >> $COMMON_VARS

THIS_COUNTS_HTSEQ_PREFIX=$OUT_DIR/$COUNTS_HTSEQ_PREFIX

# 'sort -s -k 1,1': Sort SAM file by read name before HTSeq
SAM_SORTED=$BAM_PREFIX\_sorted_read_name.sam
if [ ! -e $SAM_SORTED ] ; then
    samtools view $BAM | \
	sort -s -k 1,1 >$SAM_SORTED
fi

HTSeqCount_options='--stranded=no'
/usr/local/bin/htseq-count \
    $HTSeqCount_options $SAM_SORTED $GFF \
    >$THIS_COUNTS_HTSEQ_PREFIX.txt\
    2>$THIS_COUNTS_HTSEQ_PREFIX.err
######## END gene count estimation via HTSeq #########


$SCRIPTS_DIR/rna-seq_diff_exprn_pipeline_circos.sh \
    $OUT_DIR $BAM $SAM_SORTED $GENDER $ID $COMMON_VARS


############# BEGIN DEXSeq counts ##################
DEXSEQ_OUT=$OUT_DIR/dexseq_counts.txt
if [ ! -e $DEXSEQ_OUT ]; then
    samtools view $BAM | \
	dexseq_count.py --stranded=no \
	$DEXSEQ_GFF - $DEXSEQ_OUT
fi
########################