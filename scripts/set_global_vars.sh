#!/bin/sh

# do this:
# source set_global_vars.sh
# so that all these variables are available

SINGLE_CELL_DIR=/home/obot/single-cell
SCRIPTS_DIR=$SINGLE_CELL_DIR/scripts
DATA_DIR=$SINGLE_CELL_DIR/data
RESULTS_DIR=$SINGLE_CELL_DIR/results
EXPR_DIR=$RESULTS_DIR/expression
SV_DIR=$RESULTS_DIR/structural-variation
CIRCOS_DIR=$RESULTS_DIR/circos

BED=/share/references/hg19/hg19_ucsc-genes.bed
GFF=/share/references/hg19/hg19_ucsc-genes.gtf

# Structure/location of BAM files
# Make a sed-able template
# can replace #ID# with whatever the id number is via:
# SAMPLE_DIR=`echo $SAMPLE_DIR_TEMPLATE | sed 's/#ID#/'"$ID"'/'`
# BAM=`echo $BAM_PREFIX_TEMPLATE | sed 's/#ID#/'"$ID"'/'`
# where $ID is the identification number of the sample
# in your script
SAMPLE_DIR_TEMPLATE=$DATA_DIR/Sample_#ID#
BAM_PREFIX_TEMPLATE=$OUT_DIR_TEMPLATE/tophat_6n6_trimmed_#ID#_merged_rmdup

# File that links the sample numbers to their treatment groups
# via a tab-delimited file:
# [obot@grettir expression]$ head ../../data/samplenumbers_to_treatmentgroups.txt
# 3        Untreated
# 4        Untreated
# 6        Untreated
# 7        Untreated
# 8        Untreated
# 9        Untreated
# 34       HighDose
# 35       HighDose
# 36       HighDose
# 38       HighDose
METADATA_FILE=$DATA_DIR/samplenumbers_to_treatmentgroups.txt
SAMPLE_IDS=`sed 1d $METADATA_FILE | cut -f1 - | tr "\n" " "`
TREATMENT_GROUPS=`cut -f2 $METADATA_FILE | uniq | tr "\n" " "`