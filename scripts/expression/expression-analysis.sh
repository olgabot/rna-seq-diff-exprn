#!/bin/sh

IN_DIR=$1
OUT_DIR=$2
SAMPLE_ID=$3

GENE_SYMBOLS=/home/obot/single-cell/results/expression/hg19_ucsc-genes_symbol.tab

# Take only columns 4 and 13 of the bedtools counts
# file so we get a file that looks like the HTSEQ file:
# uc022cox.1      0
# uc022coy.1      0
# uc022coz.1      0
# uc022cpa.1      0
# uc022cpb.1      0
# uc022cpc.1      0
# uc022cpd.1      0
# uc022cpe.1      0
# uc022cpf.1      0
# uc022cpg.1      0
# Also need to sort the file by gene name because that's
# what HTSEQ goes by, and this way we can compare counts
# between the counting methods
BEDTOOLS_COUNTS_ORIG=$IN_DIR/ref_bedtools_counts.txt
BEDTOOLS_COUNTS=$IN_DIR/ref_bedtools_counts_cleaned.txt
cut -f 4,13 $BEDTOOLS_COUNTS_ORIG \
    | sort -k 1,1 > $BEDTOOLS_COUNTS

HTSEQ_COUNTS_ORIG=$IN_DIR/ref_htseq_counts.txt
HTSEQ_COUNTS=$IN_DIR/ref_htseq_counts_cleaned.txt
LINES=`wc -l $HTSEQ_COUNTS_ORIG | awk ' { print $1 } '` 
# Exclude the last 5 lines of the HTSEQ_COUNTS file,
# which looks like this:
#
# no_feature      12400494
# ambiguous       4939662
# too_low_aQual   0
# not_aligned     0
# alignment_not_unique    6588235

END=`expr $LINES - 5`
head -n $END $HTSEQ_COUNTS_ORIG > $HTSEQ_COUNTS

COMBINED_COUNTS_TMP1=$IN_DIR/ref_combined_counts.tmp1
COMBINED_COUNTS_TMP2=$IN_DIR/ref_combined_counts.tmp2
HEADER='gene\tbedtools\thtseq'
cut -f 2 $HTSEQ_COUNTS \
    | paste $BEDTOOLS_COUNTS - > $COMBINED_COUNTS_TMP1
echo -e $HEADER | cat - $COMBINED_COUNTS_TMP1 \
    > $COMBINED_COUNTS_TMP2

COMBINED_COUNTS=$IN_DIR/ref_combined_counts.txt
cut -f 2,3 $COMBINED_COUNTS_TMP2 | \
    paste $GENE_SYMBOLS - > $COMBINED_COUNTS
rm $COMBINED_COUNTS_TMP1 ; rm $COMBINED_COUNTS_TMP2