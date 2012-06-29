#!/bin/sh

SINGLE_CELL_DIR=/home/obot/single-cell
EXPR_DIR=$SINGLE_CELL_DIR/results/expression
MAPPED_DIR=$SINGLE_CELL_DIR/data
SCRIPTS_DIR=$SINGLE_CELL_DIR/scripts

# File that has UCSC transcript IDs in column 1
# and gene symbol in column 2
UCSC_SYMBOL=$EXPR_DIR/hg19_ucsc-genes_symbol_no-header.tab

IN_FILES_BED=$UCSC_SYMBOL
COL_NAMES="#ucsc_id,geneSymbol"

COUNTS_BED=ref_bedtools_counts.txt
COUNTS_BED_TAB=ref_bedtools_counts_name-count.txt

# Need to skip the last 5 lines of an htseq-counts file
COUNTS_HTSEQ=ref_htseq_counts.txt
HTSEQ_TAB=ref_htseq_counts_notail.txt
IN_FILES_HTSEQ=$UCSC_SYMBOL

# For DEXSeq
DEXSEQ_GFF=/home/obot/bed_and_gtf/hg19_ucsc-genes_dexseq.gff

for i in `cut -f1 \
$SINGLE_CELL_DIR/untreated_highdose_survivors.txt  \
| tr "\n" " "`
do
    echo "making name-count file for Sample_$i"
    SAMPLE_DIR=$MAPPED_DIR/Sample_$i
    BAM=$SAMPLE_DIR/tophat_6n6_trimmed_$i\_merged_rmdup.bam
    THIS_COUNTS_BED=$SAMPLE_DIR/$COUNTS_BED
    THIS_COUNTS_TAB=$SAMPLE_DIR/$COUNTS_TAB
    if [ ! -e $THIS_COUNTS_TAB ]; then
	$SCRIPTS_DIR/make-bed-name-count.sh $THIS_COUNTS_BED \
	    $THIS_COUNTS_TAB
    fi
    IN_FILES_BED=$IN_FILES_BED,$SAMPLE_DIR/$COUNTS_TAB
    COL_NAMES=$COL_NAMES,Sample_$i

    THIS_COUNTS_HTSEQ=$SAMPLE_DIR/$COUNTS_HTSEQ
    THIS_HTSEQ_TAB=$SAMPLE_DIR/$HTSEQ_TAB
    if [ ! -e $THIS_COUNTS_HTSEQ ]; then
	REF_GFF=/home/obot/bed_and_gtf/hg19_ucsc-genes.gtf
	REF_HTSEQ_PREFIX=$SAMPLE_DIR/ref_htseq_counts
	samtools view $BAM | sort -k 1,1 |
	/usr/local/bin/htseq-count \
	    --stranded=no - $REF_GFF \
	    >$REF_HTSEQ_PREFIX.txt\
	    2>$REF_HTSEQ_PREFIX.err
    fi
    if [ ! -e $THIS_HTSEQ_TAB ]; then
    # Take everything except the last 5 lines
	head -n \
	    `wc -l $THIS_COUNTS_HTSEQ | awk '{ print $1-5 }'` \
	    $THIS_COUNTS_HTSEQ > $THIS_HTSEQ_TAB
    fi
    IN_FILES_HTSEQ=$IN_FILES_HTSEQ,$THIS_HTSEQ_TAB
    DEXSEQ_OUT=$SAMPLE_DIR/dexseq_counts.txt
    if [ ! -e $DEXSEQ_OUT ]; then
	samtools view $BAM | \
	    dexseq_count.py --stranded=no \
	    $DEXSEQ_GFF - $DEXSEQ_OUT
    fi
done

for GROUP in `cut -f2 /home/obot/single-cell/untreated_\
highdose_survivors.txt | uniq | tr "\n" " "`
do
  THIS_GROUP_DIR=$MAPPED_DIR/treatment_groups/$GROUP
  echo "making results/expression tab files for" $GROUP
  for i in 1 2
  do
      G_DIR=$THIS_GROUP_DIR/$GROUP$i
      THIS_COUNTS_BED=$G_DIR/$COUNTS_BED
      THIS_COUNTS_TAB=$G_DIR/$COUNTS_TAB
      if [ ! -e $THIS_COUNTS_TAB ]; then
      $SCRIPTS_DIR/make-bed-name-count.sh $THIS_COUNTS_BED \
	  $THIS_COUNTS_TAB
      fi
      IN_FILES_BED=$IN_FILES_BED,$THIS_COUNTS_TAB
      COL_NAMES=$COL_NAMES,$GROUP$i

      # .bam comes before _ref.bam and _sino.bam so this
      # will get the real bam file
#      echo G_DIR: $G_DIR
      BAM=`ls $G_DIR/*.bam  | head -n 1`
#      echo after making bam file

      THIS_COUNTS_HTSEQ=$G_DIR/$COUNTS_HTSEQ
      THIS_HTSEQ_TAB=$G_DIR/$HTSEQ_TAB
    echo THIS_HTSEQ_TAB $THIS_HTSEQ_TAB
#      if [ ! -e $THIS_COUNTS_HTSEQ ]; then
	  REF_GFF=/home/obot/bed_and_gtf/hg19_ucsc-genes.gtf
	  REF_HTSEQ_PREFIX=$G_DIR/ref_htseq_counts
	  samtools view $BAM | sort -k 1,1 |
	  /usr/local/bin/htseq-count \
	      --stranded=no - $REF_GFF \
	      >$REF_HTSEQ_PREFIX.txt\
	    2>$REF_HTSEQ_PREFIX.err
#      fi


#      if [ ! -e $THIS_HTSEQ_TAB ]; then
    # Take everything except the last 5 lines
	  head -n \
	      `wc -l $THIS_COUNTS_HTSEQ | awk '{print $1-5}'`\
	      $THIS_COUNTS_HTSEQ > $THIS_HTSEQ_TAB
#      fi
      IN_FILES_HTSEQ=$IN_FILES_HTSEQ,$THIS_HTSEQ_TAB
      DEXSEQ_OUT=$G_DIR/dexseq_counts.txt
      if [ ! -e $DEXSEQ_OUT ]; then
	  samtools view $BAM | \
	      dexseq_count.py --stranded=no \
	      $DEXSEQ_GFF - $DEXSEQ_OUT
      fi

  done
done

BED_COUNTS_TABLE=$EXPR_DIR/counts_bedtools.tab
HTSEQ_COUNTS_TABLE=$EXPR_DIR/counts_htseq.tab

$SCRIPTS_DIR/merge-results/expression-counts.sh $IN_FILES_BED \
    $COL_NAMES $BED_COUNTS_TABLE

$SCRIPTS_DIR/merge-results/expression-counts.sh $IN_FILES_HTSEQ \
    $COL_NAMES $HTSEQ_COUNTS_TABLE