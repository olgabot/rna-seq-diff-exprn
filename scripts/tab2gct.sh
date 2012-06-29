#!/bin/sh

# tab2gct.sh
#
# Created by Olga Botvinnik on 2012-05-18.

IN_TAB=$1

# Add "required" line of .gct file
echo "#1.2"

# Get number of lines in the tab file, subtract one to get the lines of
# data without the header
echo `wc -l counts_per-symbol-max.tab | awk ' { print ($1-1) } '`"\c"
echo "\t\c"

# Get number of columns in tab file, subtract one to get the actual
# columns of data
echo `head -n 1 counts_per-symbol-max.tab | awk ' { print NF-1 } '`

# Copy the first column as the "description" for gct format
# Rename the first column to "Name" so ssGSEA knows to use that column
# when finding the genes
ORIG_GENE_COL_NAME=`cut -f 1 $IN_TAB | head -n 1`
cut -f 1 $IN_TAB | sed 's/'"$ORIG_GENE_COL_NAME"'/Name/' | paste - $IN_TAB