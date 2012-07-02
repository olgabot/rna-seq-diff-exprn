#!/bin/sh

############# BEGIN Initialization ###############
PATH=/usr/local/bin/python2.7:$PATH 
PATH=$PATH:/share/apps/samtools-0.1.18:/share/apps/BEDTools-Version-2.15.0/
export PATH
export LC_ALL=POSIX # Forces char encoding to be POSIX

# Save all the variables into a common variable .sh file
# Variables that are specific to a single sample (ie a BAM
# file or an ID) are NOT stored, but the comma-separated
# list of *all* BAM files or IDs *are* stored
COMMON_VARS=$1
if [ -e $COMMON_VARS ] ; then
    # remove the common-variables file to make sure 
    # we're starting anew
    rm $COMMON_VARS
fi
echo '#!/bin/sh\n' | cat - >> $COMMON_VARS
SCRIPTS_DIR=/share/apps/bin
echo 'SCRIPTS_DIR=$SCRIPTS_DIR' | cat - >> $COMMON_VARS
CIRCOS_BIN=/share/apps/circos-0.60/bin/circos
echo 'CIRCOS_BIN=$CIRCOS_BIN' | cat - >> $COMMON_VARS
echo 'GENOME=/share/apps/BEDTools-Version-2.15.0/genomes/human.hg19.genome\n' \
    | cat - >> $COMMON_VARS
echo '# File that has UCSC transcript IDs in column 1\n# and gene symbol in column 2\nUCSC_SYMBOL=/share/reference/hg19/genes/hg19_ucsc-genes_symbol_no-header.tab' \
    | cat - >> $COMMON_VARS
echo '# For DEXSeq\nDEXSEQ_GFF=/home/obot/bed_and_gtf/hg19_ucsc-genes_dexseq.gff'
############# END Initialization ###############

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


############## Globally-used variables ##############
echo '# Globally used variables' | cat - >> $COMMON_VARS

# Where you want the output of all this craziness to go
# (This is a single directory: rseqc, gene counts, and genome 
# coverage output go to $BASE_OUT_DIR/$ID for each sample)
# e.g. ~/data
BASE_OUT_DIR=$2
echo 'BASE_OUT_DIR=$BASE_OUT_DIR' | cat - >> $COMMON_VARS

# tab-delimited file of the conditions of each sample, e.g.:
# bam_prefix                   id      group     gender read_type
# (no header in actual file, this is an example)
# ~/data/sample1/accepted_hits sample1 untreated female paired_end
# ~/data/sample2/accepted_hits sample2 untreated female paired_end
# ~/data/sample3/accepted_hits sample3 untreated male paired_end
# ~/data/sample4/accepted_hits sample4 treated female paired_end
# ~/data/sample5/accepted_hits sample5 treated female paired_end
# ~/data/sample6/accepted_hits sample6 treated male paired_end
COND=$3
echo 'COND=$COND' | cat - >> $COMMON_VARS

COND_WITHOUT_COMMENTS=`sed -e 's/#.*//' -e 's/[ ^I]*$//' -e '/^$/ d' < $COND`

# The bam prefixes of the files you want to use (comma-separated)
BAM_PREFIXES=`echo $COND_WITHOUT_COMMENTS | cut -f1 | tr "\n" ","`
echo 'BAM_PREFIXES=$BAM_PREFIXES' | cat - >> $COMMON_VARS

# Sample identification for each sample (comma-separated)
IDS=`echo $COND_WITHOUT_COMMENTS | cut -f2 | tr "\n" ","`
echo 'IDS=$IDS' | cat - >> $COMMON_VARS

# Treatment groups (space-separated)
GROUPS=`echo $COND_WITHOUT_COMMENTS | cut -f3 | uniq | tr "\n" " "`
echo 'GROUPS=$GROUPS' | cat - >> $COMMON_VARS

# To determine whether ChrY should be included or omitted
# from the analyses, specify the gender of each sample
# (comma-separated)
GENDERS=`echo $COND_WITHOUT_COMMENTS | cut -f3 | tr "\n" ","` 
echo 'GENDERS=$GENDERS\n' | cat - >> $COMMON_VARS

# GFF files that you want to use to estimate gene counts. 
# used by: HTSEQ and DEXSeq (differential exon usage)
GFF=$4
echo 'GFF=$GFF' | cat - >> $COMMON_VARS

# BED files that you want to use to estimate gene counts.
# Used by: bedtools coverage
BED=$5
echo 'BED=$BED' | cat - >> $COMMON_VARS

# Expression results output location
# Will create coverageBed and HTSeq folders in this location,
# which will have the combined gene expression tables
# in the order specified by the conditions file (COND)
EXPR_DIR=$6
echo 'EXPR_DIR=$EXPR_DIR' | cat - >> $COMMON_VARS

# Circos results output location
# Will create a folder for each sample in this location
# (each folder's name determined by sample id)
CIRCOS_OUT_DIR=$7
echo 'CIRCOS_OUT_DIR=$CIRCOS_OUT_DIR' | cat - >> $COMMON_VARS

# Create a script of commonly used variables to avoid
# sending them everywhere

if ( ${#*} > 7 );
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
NUM_GROUPS=$8
else
    NUM_GROUPS=0
fi 
echo 'NUM_GROUPS=$NUM_GROUPS\n' | cat - >> $COMMON_VARS
############## END Globally-used variables ##############


############## Begin scripting! ##############

# Get number of files to iterate over
# Subtract one because bash uses 0-based indexing,
# while awk is 1-based
echo '# Array variables' | cat - >> $COMMON_VARS
END=`awk -F"," '{print NF-1}' <$BAM_PREFIXES`
echo 'END=$END' | cat - >> $COMMON_VARS

# Make arrays for easy access to corresponding
# bam files / out directories / ids
declare -a BAM_PREFIX_ARRAY=( `echo $BAM_PREFIXES | tr "," " "`)
echo 'declare -a BAM_PREFIX_ARRAY=$BAM_PREFIX_ARRAY' \
    | cat - >> $COMMON_VARS

#declare -a OUT_DIR_ARRAY=( `echo $OUT_DIRS | tr "," " "`)
#echo 'declare -a OUT_DIR_ARRAY=$OUT_DIR_ARRAY' \
#    | cat - >> $COMMON_VARS

declare -a ID_ARRAY=( `echo $IDS | tr "," " "`)
echo 'declare -a ID_ARRAY=$ID_ARRAY' \
    | cat - >> $COMMON_VARS

declare -a GROUPS_ARRAY=( `echo $GROUPS` )
echo 'declare -a GROUPS_ARRAY=$GROUPS_ARRAY' \
    | cat - >> $COMMON_VARS

declare -a GENDER_ARRAY=( `echo $GENDERS | tr "," " "`)
echo 'declare -a GENDER_ARRAY=$GENDER_ARRAY' \
    | cat - >> $COMMON_VARS

if( ${#BAM_PREFIX_ARRAY[@]} ne ${#ID_ARRAY[@]} ); do
   error_exit "number of bam prefixes does not match number of out directories specified. stopping."
fi
#if( ${#ID_ARRAY[@]} ne ${#OUT_DIR_ARRAY[@]} ); do
#   error_exit "number of out directories does not match number of ids specified. stopping."
#fi

for ((i=0;i<=$END;++i)); do
    BAM_PREFIX=${BAM_PREFIX_ARRAY[$i]}
    ID=${ID_ARRAY[$i]}
    OUT_DIR=$BASE_OUT_DIR/$ID

    # Do quality control on this via 
    # RNA-Seq-QualityControl, aka RSeqQC:
    rseqc.sh $BAM_PREFIX.bam $OUT_DIR $BED

# This ...gene_counts.sh script does:
# 1. Estimate gene counts from BAM via bedtools and HTSeq
# 2. Calls another script to:
#    a. Estimate genome-wide coverage via bedtools and HTSeq
#    b. Make circos plot of coverage for this sample
# 3. Estimates differential exon usage counts using
#    dexseq_counts.py
    $SCRIPTS_DIR/rna-seq_diff_exprn_pipeline_gene_counts.sh \
	$BAM_PREFIX $OUT_DIR \
	$GENDER $ID $COMMON_VARS
done

# Merge the gene counts
$SCRIPTS_DIR/make_expression_counts_table.sh \
    $COMMON_VARS

# Do differential expression analysis