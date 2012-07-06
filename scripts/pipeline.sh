#!/bin/sh -x

# Example run:
# scripts/pipeline.sh test-results test-data/conditions.tab test-data/hg19_ucsc_genes.gtf test-data/hg19_ucsc_genes_dexseq.gtf test-data/hg19_ucsc_genes.bed test-data/hg19_id_symbol.txt test-data/human.hg19.genome 1
# A bash version of above:
# scripts/pipeline.sh $OUT_DIR $CONDITIONS $GTF $DEXSEQ_GTF $BED $TXPTID_SYMBOL $GENOME $NUM_GROUPS

############# BEGIN Initialization ###############
# PATH=/usr/local/bin/python2.7:$PATH 
# PATH=$PATH:/share/apps/samtools-0.1.18:/share/apps/BEDTools-Version-2.15.0/
# export PATH

# Force char encoding to be POSIX, necessary for HTSeq
export LC_ALL=POSIX
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
# Use sed to remove trailing forward slashes
BASE_OUT_DIR=`echo $1 | sed 's:/$::'`

if [[ ! -d $BASE_OUT_DIR ]]; then
    # If this directory does not yet exist, create it
    mkdir -p $BASE_OUT_DIR
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

# Check which directory this file is in - this is the 'scripts'
# directory that all the other scripts are in
SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "SCRIPTS_DIR='$SCRIPTS_DIR'" | cat - >> $COMMON_VARS

# 
CIRCOS_BIN="`which circos`"
echo "CIRCOS_BIN='$CIRCOS_BIN'" | cat - >> $COMMON_VARS

# Finally, write the base out directory
echo "BASE_OUT_DIR='$BASE_OUT_DIR'" | cat - >> $COMMON_VARS


# tab-delimited file of the conditions of each sample, e.g.:
# bam_prefix                   id      group     gender read_type
# (no header in actual file, this is an example)
# ~/data/sample1/accepted_hits sample1 untreated female paired_end
# ~/data/sample2/accepted_hits sample2 untreated female paired_end
# ~/data/sample3/accepted_hits sample3 untreated male paired_end
# ~/data/sample4/accepted_hits sample4 treated female paired_end
# ~/data/sample5/accepted_hits sample5 treated female paired_end
# ~/data/sample6/accepted_hits sample6 treated male paired_end
COND="$2"
echo "COND='$COND'" | cat - >> $COMMON_VARS

COND_WITHOUT_COMMENTS=$COND.without_comments
sed -e 's/#.*//' -e 's/[ ^I]*$//' -e '/^$/ d' \
    < $COND > $COND_WITHOUT_COMMENTS

# The bam prefixes of the files you want to use (comma-separated)
BAM_PREFIXES=`cut -f1 $COND_WITHOUT_COMMENTS | tr "\n" "," | sed 's:,$::'`
echo "BAM_PREFIXES='$BAM_PREFIXES'" | cat - >> $COMMON_VARS

# Sample identification for each sample (comma-separated)
IDS=`cut -f2 $COND_WITHOUT_COMMENTS | tr "\n" "," | sed 's:,$::' `
echo "IDS='$IDS'" | cat - >> $COMMON_VARS

# Treatment groups (comma-separated)
TREATMENT_GROUPS=`cut -f3 $COND_WITHOUT_COMMENTS | uniq | tr "\n" "," | sed 's:,$::'`
echo "TREATMENT_GROUPS='$TREATMENT_GROUPS'" | cat - >> $COMMON_VARS

# To determine whether ChrY should be included or omitted
# from the analyses, specify the gender of each sample
# (comma-separated)
GENDERS=`cut -f4 $COND_WITHOUT_COMMENTS | tr "\n" "," | sed 's:,$::'` 
echo "GENDERS='$GENDERS'" | cat - >> $COMMON_VARS

echo "\n# Gene and species-specific variables" | cat - >> $COMMON_VARS
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
#   the TXPTID_SYMBOL file, so I recommend sticking with UCSC IDs for now)
#   table: "knownGene"
#   region: "genome"
#   output format: "GTF - gene transfer format"
# 4. output file: (whatever you want, but I suggest something informative
#    like hg19_ucsc_genes.gtf)
#    Make sure to include the file extension (.gtf) in the filename
# 5. Press "get output"
#    --> A file will be downloaded to your "Downloads" folder
GTF="$3"
echo "GTF='$GTF'" | cat - >> $COMMON_VARS

# GTF (Gene Transfer Format) file specially formatted for use with
# DEXSeq, which measures differential exon usage. To create this file,
# You need a GTF file (which can be obtained as described in the GTF
# section), and then use the script included in rna-seq-diff-exprn/scripts/external/dexseq_prepare_annotation.py:
#  $ python2.7 scripts/external/dexseq_prepare_annotation.py test-data/hg19_ucsc_genes.gtf test-data/hg19_ucsc_genes_dexseq.gtf
# [The dollar sign `$' indicates a bash shell and shows that we are
# using a command-line interface, as opposed to a command embedded
# in source code such as this document.]
DEXSEQ_GTF="$4"
echo "DEXSEQ_GTF='$DEXSEQ_GTF'" | cat - >> $COMMON_VARS

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
#   the TXPTID_SYMBOL file, so I recommend sticking with UCSC IDs for now)
#   table: "knownGene"
#   region: "genome"
#   output format: "BED - browser extensible data"
# 4. output file: (whatever you want, but I suggest something informative
#    like hg19_ucsc_genes.gtf)
#    Make sure to include the file extension (.gtf) in the filename
# 5. Press "get output"
#    --> A file will be downloaded to your "Downloads" folder
BED="$5"
echo "BED='$BED'" | cat - >> $COMMON_VARS

# Tab-delimited file with the transcript ID in column 1 and the gene
# symbols you'd like to use in column 2, without a header line. You can
# get one of these files by:
# 1. Go to http://genome.ucsc.edu/cgi-bin/hgTables
# 2. Choose your clade and organism of interest
# 3. Choose these settings:
#   group: "Genes and Gene Prediction Tracks"
#   track: (whatever you want, but these instructions are built on using
#   "UCSC Genes." You can use Ensembl or other transcript IDs but then
#   you will need to choose different columns from the kgXref file for
#   the TXPTID_SYMBOL file, so I recommend sticking with UCSC IDs for now)
#   table: "kgXref"
#   region: "genome"
#   output format: "GTF - gene transfer format"
# 4. output file: (whatever you want, but I suggest something informative
#    like hg19_kgXref.txt)
#    Make sure to include the file extension (.gtf) in the filename
# 5. Press "get output"
# 6. Now you need to take an extra step to get just the UCSC IDs, 
#    e.g. uc002gig.1 (column 1 in the kgXref file) and the gene symbols,
#    e.g. TP53  (column 5 in the kgXref file), a known tumor suppressor 
#    gene.
#    You could do this in Microsoft Excel, but the human file
#    (for example) has 80,923 lines in it and will crash Excel. For 
#    organisms with fewer documented genes, using Excel to push columns 
#    around may be just fine.
#    The Linux/UNIX (lovingly called *NIX) commands to take columns is 
#    called "cut" (there is also "paste" to put together columns from 
#    different files but that's out of the scope of what we're doing here)
#    We want columns 1 and 5 (the UCSC ID and the gene symbol - 
#    take a peek at the file by typing "head hg19_kgXref.txt" on the 
#    command line in the directory - this will show the first 10 lines of 
#    the file), so we'll say "cut -f 1,5" where the "-f" indicates the 
#    "fields" we want to "cut." Then we use "sed 1d" to skip the first 
#    line (skipping more than one line has a slighly different command,
#    my favorite Sed tutorial is: http://www.grymoire.com/Unix/Sed.html
#    if you're interested in learning more). And as before, the "<" 
#    indicates the input file, the "|" indicates that the output of
#    the previous command is treated as input to the next command,
#    and the ">" indicates the output file. Note that we created a *new*
#    file and did not overwrite the old one. In general, it is best 
#    practices to create a new file rather than overwrite the old one.
#    Also, if you try to make your input and output files the same, the 
#    commands may get confused and you could lose your original data. :(
#      $ cut -f 1,5 < hg19_kgXref.txt | sed 1d > hg19_id_symbol.txt
TXPTID_SYMBOL="$6"
echo "TXPTID_SYMBOL='$TXPTID_SYMBOL'" | cat - >> $COMMON_VARS

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
GENOME="$7"
echo "GENOME='$GENOME'" | cat - >> $COMMON_VARS

echo "\n# Number of groups to make from the TREATMENT_GROUPS" \
    | cat - >> $COMMON_VARS
if [[ ${#*} > 7 ]]; then
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
# 
# If you don't want ANY groupings, omit this variable from the command 
# line.
NUM_GROUPS=$8
else
    NUM_GROUPS=0
fi 
echo "NUM_GROUPS='$NUM_GROUPS'" | cat - >> $COMMON_VARS

echo "\n# More output directories" | cat - >> $COMMON_VARS
# Expression results output location
# Will create coverageBed and HTSeq folders in this location,
# which will have the combined gene expression tables
# in the order specified by the conditions file (COND)
EXPRN_DIR="$BASE_OUT_DIR/expression"
echo "EXPRN_DIR='$EXPRN_DIR'" | cat - >> $COMMON_VARS

# Circos results output location
# Will create a folder for each sample in this location
# (each folder's name determined by sample id)
CIRCOS_OUT_DIR="$BASE_OUT_DIR/circos"
echo "CIRCOS_OUT_DIR='$CIRCOS_OUT_DIR'" | cat - >> $COMMON_VARS
############## END Globally-used variables ##############


############## Begin scripting! ##############
# Make arrays for easy access to corresponding
# bam files / out directories / ids
# declare -a 
BAM_PREFIX_ARRAY=( `echo $BAM_PREFIXES | tr "," " "`)

echo '\n# Array variables' | cat - >> $COMMON_VARS

# Get number of files to iterate over
END=${#BAM_PREFIX_ARRAY[@]} #`awk -F"," '{print NF-1}' <$BAM_PREFIXES`
echo "END=$END" | cat - >> $COMMON_VARS

echo "BAM_PREFIX_ARRAY=( ${BAM_PREFIX_ARRAY[@]} )" \
    | cat - >> $COMMON_VARS

declare -a ID_ARRAY=( `echo $IDS | tr "," " "`)
echo "ID_ARRAY=( ${ID_ARRAY[@]} )" \
    | cat - >> $COMMON_VARS

declare -a GROUPS_ARRAY=( `echo $TREATMENT_GROUPS | tr "," " "` )
echo "GROUPS_ARRAY=( ${GROUPS_ARRAY[@]} )" \
    | cat - >> $COMMON_VARS

declare -a GENDER_ARRAY=( `echo $GENDERS | tr "," " "`)
echo "GENDER_ARRAY=( ${GENDER_ARRAY[@]} )" \
    | cat - >> $COMMON_VARS

if [[ ${#BAM_PREFIX_ARRAY[@]} -ne ${#ID_ARRAY[@]} ]] ; then
   error_exit "number of bam prefixes does not match number of out directories specified. stopping."
fi

for ((i=0;i<=$END;++i)); do
    BAM_PREFIX=${BAM_PREFIX_ARRAY[$i]}
    ID=${ID_ARRAY[$i]}
    GENDER=${GENDER_ARRAY[$i]}
    RSEQC_OUT_DIR=$BASE_OUT_DIR/rseqc/$ID

    if [[ ! -d $RSEQC_OUT_DIR ]]; then
        # If this directory doesn't yet exist, make this directory
        # recursively (`-p') to make all the in-between directories
        # on this path
        mkdir -p $RSEQC_OUT_DIR
    fi

    # Do quality control on this via 
    # RNA-Seq-QualityControl, aka RSeqQC:
    $SCRIPTS_DIR/rseqc.sh $BAM_PREFIX.bam $RSEQC_OUT_DIR $BED

    # Detect structural variants via SVDetect for this sample

# This ...gene_counts.sh script does:
# 1. Estimate gene counts from BAM via bedtools and HTSeq
# 2. Calls another script to:
#    a. Estimate genome-wide coverage via bedtools and HTSeq
#    b. Make circos plot of coverage for this sample
# 3. Estimates differential exon usage counts using
#    dexseq_counts.py
    $SCRIPTS_DIR/gene_counts.sh $BAM_PREFIX $EXPRN_OUT_DIR \
	   $GENDER $ID $COMMON_VARS
done

# Merge the gene counts
$SCRIPTS_DIR/make_expression_counts_table.sh \
    $COMMON_VARS

# Do differential expression analysis