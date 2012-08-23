#!/bin/sh

# Example run (without groups):
# scripts/pipeline.sh test-results test-data/conditions_chr9.tab test-data/hg19_ucsc_genes.gtf test-data/hg19_ucsc_genes_chr9_dexseq.gtf test-data/hg19_ucsc_genes.bed test-data/hg19_id_symbol.txt test-data/human.hg19.genome test-data/karyotype/karyotype.human.hg19.txt test-data/hg19_gene_density_1e5bins.txt test-data/hg19_gc_content_circos_chr9.txt

# Example run (with groups):
# scripts/pipeline.sh test-results test-data/conditions_chr9.tab test-data/hg19_ucsc_genes.gtf test-data/hg19_ucsc_genes_chr9_dexseq.gtf test-data/hg19_ucsc_genes.bed test-data/hg19_id_symbol.txt test-data/human.hg19.genome test-data/karyotype/karyotype.human.hg19.txt test-data/hg19_gene_density_1e5bins.txt test-data/hg19_gc_content_circos_chr9.txt 2

# Test for making sure there are at least two samples per
# treatment type (This *should* break)
# scripts/pipeline.sh test-results test-data/conditions_chr9_noPrEC_2.tab test-data/hg19_ucsc_genes.gtf test-data/hg19_ucsc_genes_chr9_dexseq.gtf test-data/hg19_ucsc_genes.bed test-data/hg19_id_symbol.txt test-data/human.hg19.genome test-data/karyotype/karyotype.human.hg19.txt test-data/hg19_gene_density_1e5bins.txt test-data/hg19_gc_content_circos_chr9.txt

# scripts/pipeline.sh test-results test-data/conditions_chr9.tab test-data/hg19_ucsc_genes_chr9.gtf test-data/hg19_ucsc_genes_chr9_dexseq.gtf test-data/hg19_ucsc_genes_chr9.bed test-data/hg19_id_symbol.txt test-data/human.hg19.genome test-data/karyotype/karyotype.human.hg19.txt 1

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

# --- This function isn't working and I can't figure out why --- #
# function make_absolute_path
# {
#     # Checks if this path name is an absolute path,
#     # if it's not an absolute path, make it an absolute path
#     if [[ $1 != /* ]]; then
#         echo "`pwd`/$1" | sed"s:/$::"
#     else
#         echo $1
#     fi
# }
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

# If BASE_OUT_DIR is not an absolute path (starts with `/'), make it so
if [[ $BASE_OUT_DIR != /* ]]; then
    BASE_OUT_DIR=`pwd`/$BASE_OUT_DIR
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
echo "COMMON_VARS='$COMMON_VARS'" | cat - >> $COMMON_VARS

# Find the index of a value in an array
IndexOf()    {
    local i=0 S=$1; shift
    while [ $S != $1 ]
    do    ((i++)); shift
        [ -z "$1" ] && { i=-1; break; }
    done
    echo $i
}
echo '\nIndexOf()    {
    local i=0 S=$1; shift
    while [ $S != $1 ]
    do    ((i++)); shift
        [ -z "$1" ] && { i=-1; break; }
    done
    echo $i
}\n' | cat - >> $COMMON_VARS

# Check which directory this file is in - this is the 'scripts'
# directory that all the other scripts are in
SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
echo "SCRIPTS_DIR='$SCRIPTS_DIR'" | cat - >> $COMMON_VARS

# Find where the circos binary is located
# This location needs to be added to your path to be sensed
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
if [[ $COND != /* ]]; then
    # if $COND is not an absolute path (starts with `/'), 
    # then make it one
    COND=`pwd`/$COND
fi
echo "COND='$COND'" | cat - >> $COMMON_VARS

COND_WITHOUT_COMMENTS=$COND.without_comments
sed -e 's/#.*//' -e 's/[ ^I]*$//' -e '/^$/ d' \
    < $COND > $COND_WITHOUT_COMMENTS
echo "COND_WITHOUT_COMMENTS='$COND_WITHOUT_COMMENTS'" \
    | cat - >> $COMMON_VARS

# The bam prefixes of the files you want to use (comma-separated)
BAM_PREFIXES=`cut -f1 $COND_WITHOUT_COMMENTS | \\
    tr "\n" "," | sed 's:,$::'`
BAM_PREFIXES_FULL_PATH=''
for B in `echo $BAM_PREFIXES | tr , ' '` ; do
    if [[ $B != /* ]]; then
        # if $B is not an absolute path (starts with `/'), 
        # then make it one
        B=`pwd`/$B
    fi
    BAM_PREFIXES_FULL_PATH=$BAM_PREFIXES_FULL_PATH,$B
done
BAM_PREFIXES=`echo $BAM_PREFIXES_FULL_PATH | sed 's/^,//'`
echo "BAM_PREFIXES='$BAM_PREFIXES'" | cat - >> $COMMON_VARS

# Sample identification for each sample (comma-separated)
IDS=`cut -f2 $COND_WITHOUT_COMMENTS | \\
    tr "\n" "," | sed 's:,$::' `
echo "IDS='$IDS'" | cat - >> $COMMON_VARS

# Treatment groups (unique list, comma-separated)
TREATMENT_GROUPS=`cut -f3 $COND_WITHOUT_COMMENTS | \\
    uniq | tr "\n" "," | sed 's:,$::'`
echo "TREATMENT_GROUPS='$TREATMENT_GROUPS'" | cat - >> $COMMON_VARS
TREATMENT_GROUPS_ALL=( `cut -f3 $COND_WITHOUT_COMMENTS` ) #| \\
    # tr "\n" ' '` )
# echo `cut -f3 $COND_WITHOUT_COMMENTS`
# echo ${TREATMENT_GROUPS_ALL[@]}
echo "TREATMENT_GROUPS_ALL=( ${TREATMENT_GROUPS_ALL[@]} )" | \
    cat - >> $COMMON_VARS

GROUPS_ARRAY=( `echo $TREATMENT_GROUPS | tr "," " "` )
echo "GROUPS_ARRAY=( ${GROUPS_ARRAY[@]} )" \
    | cat - >> $COMMON_VARS

# Check that there are at least two samples in each treatment group,
# preventing DESeq from getting angry when it sees that there are
# no replicates for a particular treatment group.
# --- Initial test: make sure the unique treatment groups --- #
# --- aren't the same length as ALL the treatment groups  --- #
if [[ $TREATMENT_GROUPS == `${TREATMENT_GROUPS_ALL[@]} | tr ' ' ,` ]]; then
    error_exit "The number of treatment groups is the same as the
number of samples, i.e. each sample is a member of a treatment
group different from every other sample. This means that this experiment
has no replicates, but DESeq, the differential expression package we use,
requires at least one replicate per treatment type.
Please re-run this analysis with at least one replicate per treatment
type."
fi
# Declare OBSERVED_TREATMENT_GROUPS as an array.
# This will have the same length as TREATMENT_GROUPS but
# the values at an index will be the number of times a treatment
# group is observed. This is a workaround since Bash 3 (what I have
# on my machine and what most people have since Bash 4 came out in
# 2009) does not have hash tables/associative arrays as in Bash 4.
declare -a OBSERVED_TREATMENT_GROUPS

# Iterate over all treatment groups to count them up
for G in ${TREATMENT_GROUPS_ALL[@]}; do
    i=`IndexOf $G ${GROUPS_ARRAY[@]}`
    if [[ ${OBSERVED_TREATMENT_GROUPS[$i]} == '' ]]; then
        # echo 'true'
        OBSERVED_TREATMENT_GROUPS[$i]=1
    else
        # echo 'false'
        # echo ${OBSERVED_TREATMENT_GROUPS[$i]}
        OBSERVED_TREATMENT_GROUPS[$i]=`
            echo "${OBSERVED_TREATMENT_GROUPS[$i]}+1" | bc`
    fi
done

# echo ${OBSERVED_TREATMENT_GROUPS[@]}

# Check if any treatment group has been observed only
# once, i.e. if there is any value in OBSERVED_TREATMENT_GROUPS
# equal to 1
i=`IndexOf 1 ${OBSERVED_TREATMENT_GROUPS[@]}`
if [[ $i -gt -1 ]]; then
    error_exit "One of the treatment groups, ${GROUPS_ARRAY[$i]},
has only one sample. The differential expression software used in this
package, DESeq, requires at least two samples (one replicate) per
group. Please add samples to this treatment group."
fi

# To determine whether ChrY should be included or omitted
# from the analyses, specify the gender of each sample
# (comma-separated)
GENDERS=`cut -f4 $COND_WITHOUT_COMMENTS | tr "\n" "," | sed 's:,$::'` 
echo "GENDERS='$GENDERS'" | cat - >> $COMMON_VARS

STRAND_SPECIFICITIES=`cut -f6 $COND_WITHOUT_COMMENTS | \\
    tr "\n" "," | sed 's:,$::'` 
echo "STRAND_SPECIFICITIES='$STRAND_SPECIFICITIES'" | \
    cat - >> $COMMON_VARS

echo "\n# Gene and species-specific variables" | cat - >> $COMMON_VARS
# GTF (Gene Transfer Format) files that you want to use to estimate 
# gene counts. 
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
if [[ $GTF != /* ]]; then
    # if $GTF is not an absolute path (starts with `/'), then make it one
    GTF=`pwd`/$GTF
fi
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
if [[ $DEXSEQ_GTF != /* ]]; then
    # if $DEXSEQ_GTF is not an absolute path (starts with `/'), 
    # then make it one
    DEXSEQ_GTF=`pwd`/$DEXSEQ_GTF
fi
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

# Need to the first 'track name' line in a BED file or else rseqc gets mad
BED=$5.without_track_name
sed -e 's/^track name.*//' -e '/^$/ d' \
    <$5 >$BED
if [[ $BED != /* ]]; then
    # if $BED is not an absolute path (starts with `/'), then make it one
    BED=`pwd`/$BED
fi
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
#    Also you need to remove the first line, the header of the file
#    that explains what is in which column
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
#    Or, if you want to create a chromosome-specific file like I did,
#    use your .bed file to search through your kgXref file:
#      $ cut -f 4 hg19_ucsc_genes_chr9.bed | grep --fixed-strings - hg19_kgXref.txt >hg19_kgXref_chr9.txt
#    Then do the same as above, but with your chr9 file:
#      $ cut -f 1,5 < hg19_kgXref_chr9.txt | sed 1d > hg19_id_symbol_chr9.txt
TXPTID_SYMBOL=$6.sorted

# Sort the transcript ID file so we have an established order for
# the gene counts tables.
sort -k 1,1 $6 > $TXPTID_SYMBOL
if [[ $TXPTID_SYMBOL != /* ]]; then
    # if $GENOME is not an absolute path (starts with `/'), 
    # then make it one
    TXPTID_SYMBOL=`pwd`/$TXPTID_SYMBOL
fi
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
if [[ $GENOME != /* ]]; then
    # if $GENOME is not an absolute path (starts with `/'), 
    # then make it one
    GENOME=`pwd`/$GENOME
fi
echo "GENOME='$GENOME'" | cat - >> $COMMON_VARS

# Karyotype file used by Circos, which specifies the chromosome lengths
# The third column, the chromosome name, MUST use "chr1"-type notation,
# and not the typical "hs1" notation for Homo sapiens chromosome 1
# This is because bedtools and friends use chr1 notation, but I didn't
# want to require the organism name and then lookup the conversion.
# Presumably, you would have samples from all the same organism since
# you are comparing gene expression and coverage across different
# treatments, so I felt this was a safe assumption. I also didn't want
# to lock you into ONLY using human data, because there are plenty
# of interesting organisms out there
KARYOTYPE="$8"
if [[ $KARYOTYPE != /* ]]; then
    # if $KARYOTYPE is not an absolute path (starts with `/'), 
    # then make it one
    KARYOTYPE=`pwd`/$KARYOTYPE
fi
echo "KARYOTYPE='$KARYOTYPE'" | cat - >> $COMMON_VARS

# Gene density file, can be created with 
# 1. Go to http://genome.ucsc.edu/cgi-bin/hgTables
# 2. Choose your clade and organism of interest
# 3. Choose these settings:
#   group: "All Tables"
#   table: "knownCanonical"
#   - "knownCanonical" ignores multiple isoforms that are present in
#      other knownGene files, so it's a more accurate representation
#      of gene density. If you did this with knownGene, you'd get an
#      overly optimistic view of gene density
#   region: "genome"
#   output format: "all fields from selected table"
# 4. output file: (whatever you want, but I suggest something informative
#    like hg19_ucsc_knownCanonical.tab)
#    Make sure to include a file extension (.tab) in the filename
#    .tab is used to show that the output is tab-delimited
# 5. Press "get output"
#    --> A file will be downloaded to your "Downloads" folder
# Then run scripts/get_gene_density.R on the file to get the
# gene density, like this:
#   $ cd rna-seq-diff-exprn/test-data
#   $ ../scripts/get_gene_density.R hg19_ucsc_knownCanonical.tab hg19_gene_density.tab
GENE_DENSITY_FILE=$9
if [[ $GENE_DENSITY_FILE != /* ]]; then
    # if $GENE_DENSITY_FILE is not an absolute path (starts with `/'), 
    # then make it one
    GENE_DENSITY_FILE=`pwd`/$GENE_DENSITY_FILE
fi
echo "GENE_DENSITY_FILE='$GENE_DENSITY_FILE'" | cat - >> $COMMON_VARS

MAX_GENE_DENSITY=`sed 1d $GENE_DENSITY_FILE | \\
awk '{if ($4>max) max=$4} END{print max}'`
echo "MAX_GENE_DENSITY='$MAX_GENE_DENSITY'" | cat - >> $COMMON_VARS

GENE_DENSITY_MIN_VALUE_CHANGE=`echo "scale=10;\\
    $MAX_GENE_DENSITY/100" | bc`
echo "GENE_DENSITY_MIN_VALUE_CHANGE='$GENE_DENSITY_MIN_VALUE_CHANGE'" \
    | cat - >> $COMMON_VARS


# GC content file, can be created by converting a .wig (wiggle)
# format file that's used for a genome browser into a circos format
# using:
# ../scripts/wig_to_circos.R hg19_gc1000Base.txt hg19_gc_content.circos
GC_CONTENT_FILE="${10}"
if [[ $GC_CONTENT_FILE != /* ]]; then
    # if $GC_CONTENT_FILE is not an absolute path (starts with `/'), 
    # then make it one
    GC_CONTENT_FILE=`pwd`/$GC_CONTENT_FILE
fi
echo "GC_CONTENT_FILE='$GC_CONTENT_FILE'" | cat - >> $COMMON_VARS

echo "\n# Number of groups to make from the TREATMENT_GROUPS" \
    | cat - >> $COMMON_VARS
if [[ ${#*} > 10 ]]; then
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
    NUM_GROUPS=${11}
else
    NUM_GROUPS=0
fi 

if [[ $NUM_GROUPS == 1 ]]; then
    error_exit "Cannot specify a single group for each treatment type. 
This is because DESeq (used in differential expression analysis)
breaks if there are no biological replicates of an experiment."
fi
echo "NUM_GROUPS='$NUM_GROUPS'" | cat - >> $COMMON_VARS

echo "\n# More output directories" | cat - >> $COMMON_VARS
# Expression results output location
# Will create coverageBed and HTSeq folders in this location,
# which will have the combined gene expression tables
# in the order specified by the conditions file (COND)
EXPRN_DIR="$BASE_OUT_DIR/expression"
echo "EXPRN_DIR='$EXPRN_DIR'" | cat - >> $COMMON_VARS

# BEDTools and HTSeq are the two main programs used to estimate
# expression on a per-gene level
BEDTOOLS_DIR="$EXPRN_DIR/bedtools"
echo "BEDTOOLS_DIR='$BEDTOOLS_DIR'" | cat - >> $COMMON_VARS

HTSEQ_DIR="$EXPRN_DIR/htseq"
echo "HTSEQ_DIR='$HTSEQ_DIR'" | cat - >> $COMMON_VARS

DEXSEQ_DIR="$EXPRN_DIR/dexseq"
echo "DEXSEQ_DIR='$DEXSEQ_DIR'" | cat - >> $COMMON_VARS

DEXSEQ_OUT=dexseq_counts.txt
echo "DEXSEQ_OUT='$DEXSEQ_OUT'" | cat - >> $COMMON_VARS

# Circos results output location
# Will create a folder for each sample in this location
# (each folder's name determined by sample id)
CIRCOS_OUT_DIR="$BASE_OUT_DIR/circos"
echo "CIRCOS_OUT_DIR='$CIRCOS_OUT_DIR'" | cat - >> $COMMON_VARS

CIRCOS_TEMPLATE_DIR="$SCRIPTS_DIR/circos-templates"
echo "CIRCOS_TEMPLATE_DIR='$CIRCOS_TEMPLATE_DIR'" | \
    cat - >> $COMMON_VARS



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

ID_ARRAY=( `echo $IDS | tr "," " "`)
echo "ID_ARRAY=( ${ID_ARRAY[@]} )" \
    | cat - >> $COMMON_VARS

GENDER_ARRAY=( `echo $GENDERS | tr "," " "`)
echo "GENDER_ARRAY=( ${GENDER_ARRAY[@]} )" \
    | cat - >> $COMMON_VARS

STRAND_ARRAY=( `echo $STRAND_SPECIFICITIES | tr "," " "`)
echo "STRAND_ARRAY=( ${STRAND_ARRAY[@]} )" \
    | cat - >> $COMMON_VARS


if [[ ${#BAM_PREFIX_ARRAY[@]} -ne ${#ID_ARRAY[@]} ]] ; then
   error_exit "number of bam prefixes does not match number of out directories specified. stopping."
fi

###### Some more variables for gene counts files
echo "\n# Variables for gene counts and genome-wide coverage files"\
    | cat - >> $COMMON_VARS
COUNTS_BED=bedtools_gene_counts.txt
echo "COUNTS_BED='$COUNTS_BED'" | cat - >> $COMMON_VARS

COUNTS_HTSEQ_PREFIX=htseq_gene_counts
COUNTS_HTSEQ=$COUNTS_HTSEQ_PREFIX.txt
echo "COUNTS_HTSEQ='$COUNTS_HTSEQ'" | cat - >> $COMMON_VARS
echo "COUNTS_HTSEQ_PREFIX='$COUNTS_HTSEQ_PREFIX'" | cat - >> $COMMON_VARS

# --- Circos-specific variables --- #
CIRCOS_ALL_COLORS_ARRAY=( dark2-8-qual-{1,2,3,4,5,6,7,8} )
echo "CIRCOS_ALL_COLORS_ARRAY=( ${CIRCOS_ALL_COLORS_ARRAY[@]} )" \
    | cat - >> $COMMON_VARS
# TODO: cannot have more than 6 groups if gene density and GC
#       content are provided because we need different colors for
#       plotting on circos. If gene density and GC aren't provided,
#       can have up to 8 groups.

# Make array of colors for each sample, based on its group
for (( i = 0; i < ${#TREATMENT_GROUPS_ALL[@]}; i++ )); do
    # Get the name of the group for this sample index
    GROUP=${TREATMENT_GROUPS_ALL[$i]}

    # Find the index of this group in the array of unique group names
    IND=`IndexOf $GROUP ${GROUPS_ARRAY[@]}`

    # Assign this ID the IND'th color in the list of possible colors
    CIRCOS_ID_COLOR_ARRAY[$i]=${CIRCOS_ALL_COLORS_ARRAY[$IND]}
done

echo "CIRCOS_ID_COLOR_ARRAY=( ${CIRCOS_ID_COLOR_ARRAY[@]} )" \
    | cat - >> $COMMON_VARS

GC_CONTENT_COLOR=${CIRCOS_ALL_COLORS_ARRAY[6]}
echo "GC_CONTENT_COLOR='$GC_CONTENT_COLOR'" | \
    cat - >> $COMMON_VARS

GENE_DENSITY_COLOR=${CIRCOS_ALL_COLORS_ARRAY[7]}
echo "GENE_DENSITY_COLOR='$GENE_DENSITY_COLOR'" | \
    cat - >> $COMMON_VARS

UPPER_LIMITS=''
echo "UPPER_LIMITS='$UPPER_LIMITS'" | cat - >> $COMMON_VARS

MIN_VALUE_CHANGES=''
echo "MIN_VALUE_CHANGES='$MIN_VALUE_CHANGES'" | \
    cat - >> $COMMON_VARS

# echo ${CIRCOS_ID_COLOR_ARRAY[@]}
# --- End Circos-specific variables --- #
# exit

## More variables for genome-wide counts
COVERAGE_HTSEQ_PREFIX=htseq_genome_coverage
COVERAGE_HTSEQ=$COVERAGE_HTSEQ_PREFIX.wig
echo "COVERAGE_HTSEQ_PREFIX='$COVERAGE_HTSEQ_PREFIX'" \
    | cat - >> $COMMON_VARS
echo "COVERAGE_HTSEQ='$COVERAGE_HTSEQ'" | cat - >> $COMMON_VARS

COVERAGE_BEDTOOLS_PREFIX=bedtools_genome_coverage
COVERAGE_BEDTOOLS=$COVERAGE_BEDTOOLS_PREFIX.txt
echo "COVERAGE_BEDTOOLS_PREFIX='$COVERAGE_BEDTOOLS_PREFIX'" \
    | cat - >> $COMMON_VARS
echo "COVERAGE_BEDTOOLS='$COVERAGE_BEDTOOLS'" | cat - >> $COMMON_VARS

# --- Find the htseq-counts file --- #
echo 'looking for your htseq-counts file ... hopefully you installed HTSeq'

# On Olga's machine for speed:
# HTSEQ_BIN=/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/htseq-count
# --- BEGIN On anyone else's machine --- #
HTSEQ_BIN=`find /opt/local -name 'htseq-count' -exec \\
    echo - grep 'htseq-count' {} \; -quit 2>find.err | awk -F' ' '{ print $4 }' `
# --- END On anyone else's machine --- #

# The real slim shady:
# HTSEQ_BIN=`find / -name 'htseq-count' -exec echo - grep 'htseq-count' {} \; -quit 2>find.err | awk -F' ' '{ print $4 }' `

if [[ $HTSEQ_BIN == '' ]]; then
    # If the search found no file
    error_exit 'You do not have HTSeq installed.\nPlease visit http://www-huber.embl.de/users/anders/HTSeq/doc/overview.html and follow the instructions there on how to install the software package.'
else
    echo 'Way to go! You installed HTSeq. htseq-count is located:' \
        $HTSEQ_BIN
fi
echo "HTSEQ_BIN='$HTSEQ_BIN'" | cat - >> $COMMON_VARS

#for (( i = 0 ; i < $END ; ++i )); do
for (( i = 0 ; i < $END ; ++i )); do
    BAM_PREFIX=${BAM_PREFIX_ARRAY[$i]}
    ID=${ID_ARRAY[$i]}
    GENDER=${GENDER_ARRAY[$i]}
    STRAND=${STRAND_ARRAY[$i]}
    RSEQC_OUT_DIR=$BASE_OUT_DIR/rseqc/$ID

    if [[ $BAM_PREFIX != /* ]]; then
        # If the $BAM_PREFIX does not start with a `/' then it is not
        # an absolute path and we need to prepend the current directory
        # to the $BAM_PREFIX
        BAM_PREFIX=`pwd`/$BAM_PREFIX
        BAM_PREFIX_ARRAY[$i]=$BAM_PREFIX
        # echo BAM_PREFIX: $BAM_PREFIX
    fi  

    # Debugging for now, remove the rseqc files before executing
    # so that only the new ones appear in this folder:
    # rm $RSEQC_OUT_DIR/*  ------------ Debugging

    if [[ ! -d $RSEQC_OUT_DIR ]]; then
        # If this directory doesn't yet exist, make this directory
        # recursively (`-p') to make all the in-between directories
        # on this path
        mkdir -p $RSEQC_OUT_DIR
    fi

    # Do quality control on this via 
    # RNA-Seq-QualityControl, aka RSeqQC:
    # ------------ BEGIN Debugging
    pushd $RSEQC_OUT_DIR 
    $SCRIPTS_DIR/rseqc.sh $BAM_PREFIX.bam $RSEQC_OUT_DIR $BED
    popd 
    # ------------ END Debugging (uncomment this for the real thing)

    # Detect structural variants via SVDetect for this sample

# This gene_counts.sh script does:
# 1. Estimate gene counts from BAM via bedtools and HTSeq
# 2. Calls another script to:
#    a. Estimate genome-wide coverage via bedtools and HTSeq
#    b. Make circos plot of coverage for this sample
# 3. Estimates differential exon usage counts using
#    dexseq_counts.py
    # ------------ BEGIN Debugging
    $SCRIPTS_DIR/gene_counts.sh $BAM_PREFIX $EXPRN_OUT_DIR \
	   $GENDER $ID $STRAND $COMMON_VARS $i
    # ------------ END Debugging (uncomment this for the real thing)
done

# Do group gene counts
# ------------ BEGIN Debugging
if [[ $NUM_GROUPS > 0 ]]; then
    TREATMENT_GROUPS_DIR=$BASE_OUT_DIR/merged_groups_bam
    echo "TREATMENT_GROUPS_DIR='$TREATMENT_GROUPS_DIR'" | \
        cat - >> $COMMON_VARS

    $SCRIPTS_DIR/group_gene_counts.sh $COMMON_VARS \
        $TREATMENT_GROUPS_DIR
fi
# ------------ END Debugging (uncomment this for the real thing)


# Regardless of whether the number of groups was specified
# or not, plot all the samples on the same circos plot
# --- BEGIN Debugging comments --- #
$SCRIPTS_DIR/circos_all_samples.sh $COMMON_VARS
# --- END Debugging comments --- #



# Merge the gene counts onto one table
# This script also takes care of the header in the bedtools
# counts file and the 5 extra lines in the HTSEQ_COUNT file.
# --- BEGIN Debugging comments --- #
$SCRIPTS_DIR/make_gene_counts_table.sh \
    $COMMON_VARS
# --- END Debugging comments --- #

# Make a table with the gene symbols and all possible transcripts,
# and one with only the max for each gene symbol. 
# 
# For example, for
# the "all" gene counts, if we start with these counts for the
# gene TTLL11:
# uc004blr.3    TTLL11  397 0   456 450 298 207 1744    838
# uc004blt.1    TTLL11  78  0   92  98  76  53  356 166
# uc004blu.1    TTLL11  78  0   92  98  76  53  356 166
# uc011lyl.2    TTLL11  397 0   456 450 298 207 1744    838
# uc011lym.1    TTLL11  71  0   88  87  61  37  322 147
# Then we end with these counts:
# TTLL11 (uc004blr.3)   397 0   456 450 298 207 1744    838
# TTLL11 (uc004blt.1)   78  0   92  98  76  53  356 166
# TTLL11 (uc004blu.1)   78  0   92  98  76  53  356 166
# TTLL11 (uc011lyl.2)   397 0   456 450 298 207 1744    838
# TTLL11 (uc011lym.1)   71  0   88  87  61  37  322 147
#
# But for the "max" counts, if we start with the same counts for
# TTLL11 as above, then we get:
# TTLL11    397 0   456 450 298 207 1744    838
# Which was the maximum for each sample, so we can't associate it
# with a particular transcript ID.
EXPRN_SCRIPTS_DIR=$SCRIPTS_DIR/expression
echo "EXPRN_SCRIPTS_DIR='$EXPRN_SCRIPTS_DIR'" | \
    cat - >> $COMMON_VARS

BEDTOOLS_ALL=$BEDTOOLS_DIR/bedtools_counts_table_all.tab
echo "BEDTOOLS_ALL='$BEDTOOLS_ALL'" | \
    cat - >> $COMMON_VARS
BEDTOOLS_MAX=$BEDTOOLS_DIR/bedtools_counts_table_max.tab
echo "BEDTOOLS_MAX='$BEDTOOLS_MAX'" | \
    cat - >> $COMMON_VARS
R --slave --args $BEDTOOLS_DIR/bedtools_counts_table.tab \
    $BEDTOOLS_ALL $BEDTOOLS_MAX \
    < $EXPRN_SCRIPTS_DIR/make_max_and_all_txpt_counts.R

HTSEQ_ALL=$HTSEQ_DIR/htseq_counts_table_all.tab
echo "HTSEQ_ALL='$HTSEQ_ALL'" | \
    cat - >> $COMMON_VARS
HTSEQ_MAX=$HTSEQ_DIR/htseq_counts_table_max.tab
echo "HTSEQ_MAX='$HTSEQ_MAX'" | \
    cat - >> $COMMON_VARS
R --slave --args $HTSEQ_DIR/htseq_counts_table.tab \
    $HTSEQ_ALL $HTSEQ_MAX \
    < $EXPRN_SCRIPTS_DIR/make_max_and_all_txpt_counts.R

# Do differential expression analysis

# Testing:
# cd ~/workspace/rna-seq-diff-exprn
# scripts/expression/make_expression_common.sh test-results/common_variables.sh test-results/expression/expression_common.R

# --- Make tab-delimited file that has all the compared groups --- #
# --- for differential expression analysis                     --- #
# Need to source COMMON_VARS to get access to GROUP_IDS, a variable
# containing the grouping information and IDs (e.g. PrEC_group1of2) that 
# was created in group_gene_counts.sh
source $COMMON_VARS
DIFF_EXPRN_GROUPS=$EXPRN_DIR/diff_exprn_groups.tab
echo ${TREATMENT_GROUPS_ALL[@]} `echo $GROUP_IDS | \\
    sed -E 's/(_group)[[:digit:]]of[[:digit:]]/\1/g' | tr , ' '` \
    | tr ' ' "\t" > $DIFF_EXPRN_GROUPS
echo "DIFF_EXPRN_GROUPS='$DIFF_EXPRN_GROUPS'" | \
    cat - >> $COMMON_VARS

# --- Make BEDTools and HTSeq folders for DESeq figures --- #
BEDTOOLS_FIGS=$BEDTOOLS_DIR/figures
BEDTOOLS_DESEQ_DIR=$BEDTOOLS_FIGS/DESeq
if [[ ! -d $BEDTOOLS_DESEQ ]]; then
    mkdir -p $BEDTOOLS_DESEQ_DIR
fi
BEDTOOLS_DESEQ_PREFIX=$BEDTOOLS_DESEQ/deseq
echo '
# --- Variables for DESEQ analysis --- #' | cat - >> $COMMON_VARS
echo "BEDTOOLS_FIGS='$BEDTOOLS_FIGS'
BEDTOOLS_DESEQ_DIR='$BEDTOOLS_DESEQ_DIR'
BEDTOOLS_DESEQ_PREFIX='$BEDTOOLS_DESEQ_DIR/'" | \
    cat - >> $COMMON_VARS

HTSEQ_FIGS=$HTSEQ_DIR/figures
HTSEQ_DESEQ_DIR=$HTSEQ_FIGS/DESeq
if [[ ! -d $HTSEQ_DESEQ ]]; then
    mkdir -p $HTSEQ_DESEQ_DIR
fi
HTSEQ_DESEQ_PREFIX=$HTSEQ_DESEQ/deseq
echo "HTSEQ_FIGS='$HTSEQ_FIGS'
HTSEQ_DESEQ_DIR='$HTSEQ_DESEQ_DIR'
HTSEQ_DESEQ_PREFIX='$HTSEQ_DESEQ_DIR/'" | \
    cat - >> $COMMON_VARS

# --- Make an R script that contains common variables to all   --- #
# --- differential expression scripts (colors, filenames, etc) --- #
# This script also makes the (1) "all" and (2) "max" gene counts 
# tables, which contain (1) all the transcripts for a gene symbol,
# and (2) only the count for the maximum transcript for a particular
# sample and gene symbol. This is described in greater detail
# (and with an example!) in:
# rna-seq-diff-exprn/scripts/expression/make_max_and_all_txpt_counts.R
EXPRESSION_COMMON=$EXPRN_DIR/expression_common.R
echo "EXPRESSION_COMMON='$EXPRESSION_COMMON'" | \
    cat - >> $COMMON_VARS
DESEQ_SCRIPTS_DIR=$EXPRN_SCRIPTS_DIR/DESeq
echo "DESEQ_SCRIPTS_DIR='$DESEQ_SCRIPTS_DIR'" | \
    cat - >> $COMMON_VARS
$EXPRN_SCRIPTS_DIR/make_expression_common.sh $COMMON_VARS

# --- Run DESeq analysis on BEDTools and HTSeq counts --- #
R --slave --args $EXPRESSION_COMMON \
    < $DESEQ_SCRIPTS_DIR/DESeqAnalysis.R

echo "Number of seconds since starting this script: $SECONDS"
