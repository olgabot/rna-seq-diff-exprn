#!/bin/sh -xv

# Author: Olga Botvinnik (olga.botvinnik@gmail.com)
# Date: 21 June 2012
# Purpose of this script: Finds genome-wide counts
# of a BED file using both bedtools and HTSeq
# (PS: in my experience, genome-wide coverage is EXACTLY
# the same regardless of whether you use bedtools or HTSeq)


# In this file, COVERAGE = GENOME COVERAGE
# TODO: Need to consider male/female about whether to \
# include ChrY or not

# OUT_DIR=$1
BAM=$1
SAM_SORTED=$2
# BED=$4
GENDER=$3  # male or female
ID=$4
COMMON_VARS=$5

source $COMMON_VARS

echo "doing circos for" $BAM

# Don't need to make these directories because
# we already made them in the gene_counts.sh file
BEDTOOLS_DIR=$EXPRN_DIR/bedtools/$ID
HTSEQ_DIR=$EXPRN_DIR/htseq/$ID

THIS_COVERAGE_BEDTOOLS=$BEDTOOLS_DIR/$COVERAGE_BEDTOOLS
THIS_COVERAGE_HTSEQ=$HTSEQ_DIR/$COVERAGE_HTSEQ
if [[ ! -e $THIS_COVERAGE_BEDTOOLS ]] ; then
    genomeCoverageBed -bg -ibam $BAM -g $GENOME \
	   >$THIS_COVERAGE_BEDTOOLS
fi

# Genome-wide coverage via HTSeq
if [[ ! -e $THIS_COVERAGE_HTSEQ ]] ; then
    $SCRIPTS_DIR/rna-seq_diff_exprn_pipeline_htseq_coverage.py \
        -asam $SAM_SORTED -o $THIS_COVERAGE_HTSEQ
fi

COVERAGE_BEDTOOLS_CIRCOS=$BEDTOOLS_DIR/$COVERAGE_BEDTOOLS_PREFIX.circos
COVERAGE_HTSEQ_CIRCOS=$HTSEQ_DIR/$COVERAGE_HTSEQ_PREFIX.circos

echo "making circos data out of coverage data (remove chrM and add stats at top)"
# Make circos data out of coverage data (remove chrM and add stats at top)
COVERAGE_BEDTOOLS_NO_CHR_M=$BEDTOOLS_DIR/$COVERAGE_BEDTOOLS_PREFIX.no_chrM
COVERAGE_HTSEQ_NO_CHR_M=$HTSEQ_DIR/$COVERAGE_HTSEQ_PREFIX.no_chrM

# First line of HTSEQ is track type=BedGraph so need to remove first line,
# but only in HTSeq genome coverage data. Also remove
# all zero entries for faster circos processing.
# Need to remove all mitochondrial (chromosome M, chrM) 
# DNA mapping for both
grep -v chrM $THIS_COVERAGE_HTSEQ | sed 1d  \
    | awk ' { if ( $4 > 0 ) print $0 } ' \
    > $COVERAGE_HTSEQ_NO_CHR_M
grep -v chrM $THIS_COVERAGE_BEDTOOLS \
    > $COVERAGE_BEDTOOLS_NO_CHR_M

if [[ $GENDER == 'female' ]] ; then
    # Gender is female --> remove ChrY from all analyses
    grep -v chrY $COVERAGE_HTSEQ_NO_CHR_M \
        | awk ' { if ( $4 > 0 ) print $0 } ' \
        > $COVERAGE_HTSEQ_NO_CHR_Y
    grep -v chrY $COVERAGE_BEDTOOLS_NO_CHR_M \
	   > $COVERAGE_BEDTOOLS_NO_CHR_Y
    mv $COVERAGE_HTSEQ_NO_CHR_Y $COVERAGE_HTSEQ_NO_CHR_M
    mv $COVERAGE_BEDTOOLS_NO_CHR_Y $COVERAGE_BEDTOOLS_NO_CHR_M
fi

# exit

HTSEQ_STATS=`awk '{sumsq+=$4*$4; sum+=$4<0?-$4:$4} END {print int(sqrt(sumsq/NR - (sum/NR)**2)); print int(sum/NR)}' $COVERAGE_HTSEQ_NO_CHR_M | tr "\n" " " | awk ' { print "# coverage stats: mean="$1"  std dev="$2 } ' `
# test:
# awk '{sumsq+=$4*$4; sum+=$4<0?-$4:$4} END {print int(sqrt(sumsq/NR - (sum/NR)**2)); print int(sum/NR)}' genome_coverage_bedtools.txt | tr "\n" " " | awk ' { print "# coverage stats: mean="$1"  std dev="$2 } '

BEDTOOLS_STATS=`awk '{sumsq+=$4*$4; sum+=$4<0?-$4:$4} END {print int(sqrt(sumsq/NR - (sum/NR)**2)); print int(sum/NR)}' $COVERAGE_BEDTOOLS_NO_CHR_M | tr "\n" " " | awk ' { print "# coverage stats: mean="$1"  std dev="$2 } ' `
echo "finished making HTSEQ_STATS and BEDTOOLS_STATS"
echo "BEDTOOLS_STATS" $BEDTOOLS_STATS
echo "HTSEQ_STATS" $HTSEQ_STATS

echo $HTSEQ_STATS | cat - $COVERAGE_HTSEQ_NO_CHR_M \
    > $COVERAGE_HTSEQ_CIRCOS
echo $BEDTOOLS_STATS | cat - $COVERAGE_BEDTOOLS_NO_CHR_M \
    > $COVERAGE_BEDTOOLS_CIRCOS
rm $COVERAGE_BEDTOOLS_NO_CHR_M ; rm $COVERAGE_HTSEQ_NO_CHR_M
echo 'finished' $COVERAGE_HTSEQ_CIRCOS 'and' \
    $COVERAGE_BEDTOOLS_CIRCOS


######### BEGIN Average over every megabase (1Mb) #########
######### Makes circos plotting ~10x faster #########
INTERVAL=1000000  # Average over every megabase
AVG_METHOD=median # take the median over every interval

# Average the circos data over each 1Mb ($INTERVAL)
BEDTOOLS_CIRCOS=$COVERAGE_BEDTOOLS_CIRCOS.$INTERVAL
R --slave --vanilla --quiet --no-save --args \
    $COVERAGE_BEDTOOLS_CIRCOS\
    $INTERVAL \
    $AVG_METHOD \
    < $SCRIPTS_DIR/average_coverage.R \
    > $BEDTOOLS_CIRCOS

HTSEQ_CIRCOS=$COVERAGE_HTSEQ_CIRCOS.$INTERVAL
R --slave --vanilla --quiet --no-save --args \
    $COVERAGE_HTSEQ_CIRCOS\
    $INTERVAL \
    $AVG_METHOD \
    < $SCRIPTS_DIR/average_coverage.R \
    > $HTSEQ_CIRCOS
######### END Average over every megabase (1Mb) #########

####### BEGIN create circos files & execute circos #######
####### create circos files #######
THIS_CIRCOS_OUT=$CIRCOS_OUT_DIR/$ID
CIRCOS_TEMPLATE_DIR=/share/apps/circos_templates
cp -r $CIRCOS_TEMPLATE_DIR/template_individual \
    $THIS_CIRCOS_OUT

if [[ $GENDER -eq 'female' ]] ;
# Gender is female --> remove ChrY from circos plot
# by uncommenting the 'chromosomes = -chrY line
then
    sed -i -e's/^#chromosomes/chromosomes/' \
	$THIS_CIRCOS_OUT/etc/circos.conf
fi

# Make the upper limit of the plot be the mean + (2 * std dev)
UPPER_LIMIT=`echo $BEDTOOLS_STATS | tr '=' " " | awk -F' ' '{ mu=$5 ; sigma=$8 ; print mu+2*sigma}'`

# Divide UPPER_LIMIT by 100 to get minimum value change
# Last awk command is a floor function
MIN_VALUE_CHANGE=`echo $UPPER_LIMIT | awk '{print $1/100}' | awk -F. '{print $1}'`

sed -i '' -e's:UPPER_LIMIT:"$UPPER_LIMIT":' \
    $THIS_CIRCOS_OUT/etc/histogram.conf

sed -i '' -e's:MIN_VALUE_CHANGE:"$MIN_VALUE_CHANGE":' \
    $THIS_CIRCOS_OUT/etc/histogram.conf

sed -i '' -e's:BEDTOOLS_COVERAGE:"$BEDTOOLS_CIRCOS":' \
    $THIS_CIRCOS_OUT/etc/histogram.conf

sed -i '' -e's:HTSEQ_HTSEQ:"$HTSEQ_CIRCOS":' \
    $THIS_CIRCOS_OUT/etc/histogram.conf

####### execute circos #######
$CIRCOS_BIN -conf $THIS_CIRCOS_OUT/etc/circos.conf \
    -outputdir $THIS_CIRCOS_OUT
####### end create circos files & execute circos #######
