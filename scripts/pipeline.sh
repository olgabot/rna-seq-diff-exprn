#!/bin/sh -vx

# Example run:
# scripts/pipeline.sh test-results test-data/conditions.tab test-data/hg19_ucsc_genes.gtf test-data/hg19_ucsc_genes_dexseq.gtf test-data/hg19_ucsc_genes.bed test-data/hg19_id_symbol.txt test-data/human.hg19.genome 1

############# BEGIN Initialization ###############
PATH=/usr/local/bin/python2.7:$PATH 
PATH=$PATH:/share/apps/samtools-0.1.18:/share/apps/BEDTools-Version-2.15.0/
export PATH
export LC_ALL=POSIX # Forces char encoding to be POSIX

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
# Where you want the output of all this craziness to go
# (This is a single directory: rseqc, gene counts, and genome 
# coverage output go to $BASE_OUT_DIR/$ID for each sample)
# e.g. ~/data
BASE_OUT_DIR=`echo $1 | sed 's:/$::'`

if [[ ! -d $BASE_OUT_DIR ]]; then
    # If this directory does not yet exist, create it
    mkdir $BASE_OUT_DIR
fi

# Save all the variables into a common variable .sh file
# Variables that are specific to a single sample (ie a BAM
# file or an ID) are NOT stored, but the comma-separated
# list of *all* BAM files or IDs *are* stored
COMMON_VARS=$BASE_OUT_DIR/common_variables.sh
if [ -e $COMMON_VARS ] ; then
    # remove the common-variables file to make sure 
    # we're starting anew
    rm $COMMON_VARS
fi
echo '#!/bin/sh\n' | cat - >> $COMMON_VARS
echo '# Globally used variables' | cat - >> $COMMON_VARS
SCRIPTS_DIR=/share/apps/bin
echo 'SCRIPTS_DIR=$SCRIPTS_DIR' | cat - >> $COMMON_VARS
CIRCOS_BIN=/share/apps/circos-0.60/bin/circos
echo 'CIRCOS_BIN=$CIRCOS_BIN' | cat - >> $COMMON_VARS
echo 'GENOME=/share/apps/BEDTools-Version-2.15.0/genomes/human.hg19.genome\n' \
    | cat - >> $COMMON_VARS
echo '# File that has UCSC transcript IDs in column 1\n# and gene symbol in column 2\nUCSC_SYMBOL=/share/reference/hg19/genes/hg19_ucsc-genes_symbol_no-header.tab' \
    | cat - >> $COMMON_VARS
echo '# For DEXSeq\nDEXSEQ_GFF=/home/obot/bed_and_gtf/hg19_ucsc-genes_dexseq.gff' | cat - >> $COMMON_VARS
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

COND_WITHOUT_COMMENTS=$COND.without_comments
sed -e 's/#.*//' -e 's/[ ^I]*$//' -e '/^$/ d' \
    < $COND > $COND_WITHOUT_COMMENTS

# The bam prefixes of the files you want to use (comma-separated)
BAM_PREFIXES=`cut -f1 $COND_WITHOUT_COMMENTS | tr "\n" ","`
echo 'BAM_PREFIXES=$BAM_PREFIXES' | cat - >> $COMMON_VARS

# Sample identification for each sample (comma-separated)
IDS=`cut -f2 $COND_WITHOUT_COMMENTS | tr "\n" ","`
echo 'IDS=$IDS' | cat - >> $COMMON_VARS

# Treatment groups (comma-separated)
TREATMENT_GROUPS=`cut -f3 $COND_WITHOUT_COMMENTS | uniq | tr "\n" ","`
echo 'TREATMENT_GROUPS=$TREATMENT_GROUPS' | cat - >> $COMMON_VARS

# To determine whether ChrY should be included or omitted
# from the analyses, specify the gender of each sample
# (comma-separated)
GENDERS=`cut -f4 $COND_WITHOUT_COMMENTS | tr "\n" ","` 
echo 'GENDERS=$GENDERS\n' | cat - >> $COMMON_VARS

# GTF (Gene Transfer Format) files that you want to use to estimate gene 
# counts. 
# used by: HTSEQ and DEXSeq (differential exon usage)
# GTF files are a subset of GFF (General Feature Format) files
# To get one of these files, do the following steps:
# (there is probably a similar method to use the ENSEMBL website,
# but I am not familiar with it so I am giving these instructions that
# I myself have followed many times)
# 1. Go to http://genome.ucsc.edu/cgi-bin/hgTables
# 2. Choose your clade and organism of interest
# 3. Choose these settings:
#   group: "Genes and Gene Prediction Tracks"
#   track: (whatever you want, but these instructions are built on using
#   "UCSC Genes." You can use Ensembl or other transcript IDs but then
#   you will need to choose different columns from the kgXref file for
#   the )
#   table: "knownGene"
#   region: "genome"
#   output format: "GTF - gene transfer format"
# 4. output file: (whatever you want, but I suggest something informative
#    like hg19_ucsc_genes.gtf)
#    Make sure to include the file extension (.gtf) in the filename
# 5. Press "get output"
#    --> A file will be downloaded to your "Downloads" folder
GTF=$4
echo 'GFF=$GFF' | cat - >> $COMMON_VARS

# GTF (Gene Transfer Format) file specially formatted for use with
# DEXSeq, which measures differential exon usage. To create this file,
# You need a GTF file (which can be obtained as described in the GTF
# section), and then use the script included in rna-seq-diff-exprn/scripts/external/dexseq_prepare_annotation.py:
#  $ python2.7 scripts/external/dexseq_prepare_annotation.py test-data/hg19_ucsc_genes.gtf test-data/hg19_ucsc_genes_dexseq.gtf
# [The dollar sign `$' indicates a bash shell and shows that we are
# using a command-line interface, as opposed to a command embedded
# in source code such as this document.]
DEXSEQ_GTF=$5

# 
# Tab-delimited file
#   $ cut -f 1,5 < hg19_kgXref.txt | sed 1d > hg19_id_symbol.txt
TXPTID_SYMBOL=$6

# BED (Browser Extensible Data) files that you want to use to estimate 
# gene counts.
# Used by: bedtools coverage
# To get one of these files, do the following steps:
# (there is probably a similar method to use the ENSEMBL website, but I 
# am not familiar with it so I am giving these instructions that I myself
# have followed many times and can help you out with if you are stuck)
# 1. Go to http://genome.ucsc.edu/cgi-bin/hgTables
# 2. Choose your clade and organism of interest
# 3. Choose these settings:
#   group: "Genes and Gene Prediction Tracks"
#   track: (whatever you want, but these instructions are built on using
#   "UCSC Genes." You can use Ensembl or other transcript IDs but then
#   you will need to choose different columns from the kgXref file for
#   the )
#   table: "knownGene"
#   region: "genome"
#   output format: "BED - browser extensible data"
# 4. output file: (whatever you want, but I suggest something informative
#    like hg19_ucsc_genes.gtf)
#    Make sure to include the file extension (.gtf) in the filename
# 5. Press "get output"
#    --> A file will be downloaded to your "Downloads" folder
BED=$7
echo 'BED=$BED' | cat - >> $COMMON_VARS

# "Genome" file which really just says how long each chromosome is.
# An example file is included in the test-data/human.hg19.genome, 
# but you can also use ones shipped with BEDTools. On my machine, these 
# files are located in ~/packages/BEDTools/genomes:
#  $ ls ~/packages/BEDTools-Version-2.16.2/genomes/
#   human.hg18.genome human.hg19.genome mouse.mm8.genome  mouse.mm9.genome
# As to how to create these files for non-mouse or human organisms,
# my suggestion (while rather unwieldy) is to:
# 1. Go to http://genome.ucsc.edu/cgi-bin/hgTables
# 2. Choose your clade and organism of interest
# 3. Choose these settings:
#   group: "All Tables"
#   table: "chromInfo"
#   output format: "all fields from selected table"
#   output file: (anything you want, but preferably something informative
#   like platypus.ornAna1.genome)
# 4. Press "get output"
# 5. Remove the first line and the third column of the file, which you
#    could do in Microsoft Excel (since this file will be comparatively
#    small), or by using this shell command:
#      $ cut -f 1,2 < platypus.ornAna1.genome | sed 1d > platypus.ornAna1.genome.fixed
GENOME=$7

# Expression results output location
# Will create coverageBed and HTSeq folders in this location,
# which will have the combined gene expression tables
# in the order specified by the conditions file (COND)
EXPR_DIR=$BASE_OUT_DIR/expression
echo 'EXPR_DIR=$EXPR_DIR' | cat - >> $COMMON_VARS

# Circos results output location
# Will create a folder for each sample in this location
# (each folder's name determined by sample id)
CIRCOS_OUT_DIR=$BASE_OUT_DIR/circos
echo 'CIRCOS_OUT_DIR=$CIRCOS_OUT_DIR' | cat - >> $COMMON_VARS

# Create a script of commonly used variables to avoid
# sending them everywhere

if [[${#*} > 7]]; then
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
# Make arrays for easy access to corresponding
# bam files / out directories / ids
declare -a BAM_PREFIX_ARRAY=( `echo $BAM_PREFIXES | tr "," " "`)

echo '\n# Array variables' | cat - >> $COMMON_VARS

# Get number of files to iterate over
END={#BAM_PREFIX_ARRAY[@]} #`awk -F"," '{print NF-1}' <$BAM_PREFIXES`
echo 'END=$END' | cat - >> $COMMON_VARS

echo 'declare -a BAM_PREFIX_ARRAY=$BAM_PREFIX_ARRAY' \
    | cat - >> $COMMON_VARS

declare -a ID_ARRAY=( `echo $IDS | tr "," " "`)
echo 'declare -a ID_ARRAY=$ID_ARRAY' \
    | cat - >> $COMMON_VARS

declare -a GROUPS_ARRAY=( `echo $GROUPS` )
echo 'declare -a GROUPS_ARRAY=$GROUPS_ARRAY' \
    | cat - >> $COMMON_VARS

declare -a GENDER_ARRAY=( `echo $GENDERS | tr "," " "`)
echo 'declare -a GENDER_ARRAY=$GENDER_ARRAY' \
    | cat - >> $COMMON_VARS

if [[${#BAM_PREFIX_ARRAY[@]} -ne ${#ID_ARRAY[@]}]] ; then
   error_exit "number of bam prefixes does not match number of out directories specified. stopping."
fi

exit

for ((i=0;i<=$END;++i)); do
    BAM_PREFIX=${BAM_PREFIX_ARRAY[$i]}
    ID=${ID_ARRAY[$i]}
    GENDER=${GENDER_ARRAY[$i]}
    RSEQC_OUT_DIR=$BASE_OUT_DIR/rseqc/$ID

    # Do quality control on this via 
    # RNA-Seq-QualityControl, aka RSeqQC:
    $SCRIPTS_DIR/rseqc.sh $BAM_PREFIX.bam $RSEQC_OUT_DIR $BED

# This ...gene_counts.sh script does:
# 1. Estimate gene counts from BAM via bedtools and HTSeq
# 2. Calls another script to:
#    a. Estimate genome-wide coverage via bedtools and HTSeq
#    b. Make circos plot of coverage for this sample
# 3. Estimates differential exon usage counts using
#    dexseq_counts.py
    $SCRIPTS_DIR/gene_counts.sh \
	   $BAM_PREFIX $OUT_DIR \
	   $GENDER $ID $COMMON_VARS
done

# Merge the gene counts
$SCRIPTS_DIR/make_expression_counts_table.sh \
    $COMMON_VARS

# Do differential expression analysis