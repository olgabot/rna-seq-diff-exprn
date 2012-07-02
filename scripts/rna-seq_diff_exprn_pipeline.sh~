#!/bin/sh

SCRIPTS_DIR=/share/apps/bin
#SCRIPTS_DIR=/home/obot/single-cell/scripts

##### Error-handling #####
PROGNAME=$(basename $0)
function error_exit
{
#   ---------------------------
#   Function for exit due to fatal program error
#       Accepts 1 argument:
#               string containing descriptive error message
#   ---------------------------
    echo "${PROGNAME}: ${1:-"Unknown Error"}" 1>&2
    exit 1
}
##### End error-handling ####



############## Sample-specific variables ##############
# The bam files you want to use (comma-separated)
BAM_PREFIXES=$1

# Where you want the output of all this craziness to go
# (comma-separated, one for each in $BAM_PREFIXES)
OUT_DIRS=$2

# Sample identification for this sample. 
# Must be the same as in the first column of the conditions file
# I chose to do it this way so that you don't have to
# specify the BAM files in the same order that you want them
# displayed - you can specify that via the conditions file
IDS=$3

# To determine whether ChrY should be included or omitted
# from the analyses, specify the gender of each sample
# comma separated, e.g.: male,female,female,male
GENDERS=$4
############## END Sample-specific variables ##############






############## Globally-used variables ##############

# tab-delimited file of the conditions of each sample, e.g.:
# sample1 untreated
# sample2 untreated
# sample3 untreated
# sample4 treated
# sample5 treated
# sample6 treated
COND=$5

# GFF files that you want to use to estimate gene counts. 
# used by: HTSEQ and DEXSeq (differential exon usage)
GFF=$6

# BED files that you want to use to estimate gene counts.
# Used by: bedtools coverage
BED=$7

# Expression results output location
# Will create coverageBed and HTSeq folders in this location,
# which will have the combined gene expression tables
# in the order specified by the conditions file (COND)
EXPR_DIR=$8

# Circos results output location
# Will create a folder for each sample in this location


if ( ${#*} > 8 );
do
########## Only used if specified ##############
# number of groups to create for each condition group
# e.g. in the example above, if you specify 2, you will
# get a counts file that has all the individual samples and
# the grouped ones, with a header like this:
# (with spaces instead of tabs)
# sample1 sample2 sample3 sample4 sample5 sample6 untreatedGroup1 untreatedGroup2 treatedGroup1 treatedGroup2
## Where untreatedGroup1 has sample1 and sample2, while
## untreatedGroup2 contains only sample3.
# To change the ordering of samples into groupings, change
# the conditions file, as that is what all the sample
# typing goes off of.
NUM_GROUPS=$9
else
    NUM_GROUPS=0
fi 
############## END Globally-used variables ##############


############## Begin scripting! ##############

# Get number of files to iterate over
# Subtract one because bash uses 0-based indexing,
# while awk is 1-based
END=`awk -F"," '{print NF-1}' <$BAM_PREFIXES`

# Make arrays for easy access to corresponding
# bam files / out directories / ids
declare -a BAM_PREFIX_ARRAY=( `echo $BAM_PREFIXES | tr "," " "`)
declare -a OUT_DIR_ARRAY=( `echo $OUT_DIRS | tr "," " "`)
declare -a ID_ARRAY=( `echo $IDS | tr "," " "`)

if( ${#BAM_PREFIX_ARRAY[@]} ne ${#OUT_DIR_ARRAY[@]} ); do
   error_exit "number of bam prefixes does not match number of out directories specified. stopping."
fi
if( ${#ID_ARRAY[@]} ne ${#OUT_DIR_ARRAY[@]} ); do
   error_exit "number of out directories does not match number of ids specified. stopping."
fi

for ((i=0;i<=END;++i)); do
    BAM_PREFIX=${#BAM_PREFIX_ARRAY[$i]}
    OUT_DIR=${OUT_DIR_ARRAY[$i]}
    # Do quality control on this via 
    # RNA-Seq-QualityControl, aka RSeqQC:
    rseqc.sh $BAM_PREFIX.bam $OUT_DIR $BED

# samtools version in this is: samtools-0.1.18
# 
# if you want a later/different version of these, change it
# in ~line 20 of:
#  /home/obot/single-cell/scripts/gene-counts-for-pipeline.sh
    $SCRIPTS_DIR/rna-seq_diff_exprn_pipeline_gene_counts.sh \
	$BAM_PREFIX $OUT_DIR \
	$GFF $BED

# Make circos plot of coverage for this sample
    
done

# Make circos plot of coverage for all samples



# Merge the gene counts



# Do differential expression analysis