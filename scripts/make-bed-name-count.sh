#!/bin/sh

#DIR=$1
COUNTS_BED=$1 #$DIR/ref_bedtools_counts.txt
COUNTS_BED_NAME_COUNT=$2 #$DIR/ref_bedtools_counts_name-count.txt
cut -f 4,13 $COUNTS_BED > $COUNTS_BED_NAME_COUNT
sort -k 1,1 $COUNTS_BED_NAME_COUNT > tmp
mv tmp $COUNTS_BED_NAME_COUNT