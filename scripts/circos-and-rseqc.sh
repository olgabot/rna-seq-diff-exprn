#!/bin/sh

# In this file, COVERAGE = GENOME COVERAGE
# TODO: Need to consider male/female about whether to include ChrY or not

OUT_DIR=$1
BAM=$2
BED=$3
GENDER=$4  # male or female
echo "doing circos for" $BAM

# Rename to actual file name (before it was just a prefix)

HTSEQ_PREFIX=$OUT_DIR/genome_coverage_htseq
BEDTOOLS_PREFIX=$OUT_DIR/genome_coverage_bedtools
COVERAGE_HTSEQ=$HTSEQ_PREFIX.wig
COVERAGE_BEDTOOLS=$BEDTOOLS_PREFIX.txt
if [ ! -e $COVERAGE_BEDTOOLS ] ; then
    GENOME=/share/apps/BEDTools-Version-2.15.0/genomes/human.hg19.genome
    genomeCoverageBed -bg -ibam $RMDUP_BAM -g $GENOME \
	>$COVERAGE_BEDTOOLS
fi

COVERAGE_HTSEQ_CIRCOS=$OUT_DIR/genome_coverage_htseq.circos
COVERAGE_BEDTOOLS_CIRCOS=$OUT_DIR/genome_coverage_bedtools.circos

echo "making circos data out of coverage data (remove chrM and add stats at top)"
# Make circos data out of coverage data (remove chrM and add stats at top)
COVERAGE_BEDTOOLS_NO_CHR_M=$OUT_DIR/genome_coverage_bedtools.tmp
COVERAGE_HTSEQ_NO_CHR_M=$OUT_DIR/genome_coverage_htseq.tmp

# First line of HTSEQ is track type=BedGraph so need to remove first line,
# but only in HTSeq genome coverage data. Also remove
# all zero entries for faster circos processing.
# Need to remove all mitochondrial DNA mapping for both
grep -v chrM $COVERAGE_HTSEQ | sed 1d  \
    | awk ' { if ( $4 > 0 ) print $0 } ' \
    > $COVERAGE_HTSEQ_NO_CHR_M
grep -v chrM $COVERAGE_BEDTOOLS \
    > $COVERAGE_BEDTOOLS_NO_CHR_M

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
echo 'finished' $COVERAGE_HTSEQ_CIRCOS 'and' \
    $COVERAGE_BEDTOOLS_CIRCOS

INTERVAL=1000000  # Average over every megabase
AVG_METHOD=median # take the median over every interval
# Average the circos data over each 1Mb ($INTERVAL)
R --slave --vanilla --quiet --no-save --args \
    $COVERAGE_BEDTOOLS_CIRCOS\
    $INTERVAL \
    $AVG_METHOD \
    < $SCRIPTS_DIR/average_coverage.R \
    >$COVERAGE_BEDTOOLS_CIRCOS.$INTERVAL


rm $COVERAGE_BEDTOOLS_NO_CHR_M ; rm $COVERAGE_HTSEQ_NO_CHR_M

#SCRIPTS_DIR=/home/obot/single-cell/scripts
#$SCRIPTS_DIR/rseqc.sh $RMDUP_BAM $OUT_DIR $REF_BED