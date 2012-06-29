#!/bin/sh

SINGLE_CELL_DIR=/home/obot/single-cell
EXPR_DIR=$SINGLE_CELL_DIR/results/expression
MAPPED_DIR=$SINGLE_CELL_DIR/data
SCRIPTS_DIR=$SINGLE_CELL_DIR/scripts

AVG_INTERVAL=1000000
AVG_METHOD=median

COVERAGE=genome_coverage_bedtools.circos

for i in `cut -f1 $SINGLE_CELL_DIR/untreated_highdose_survivors.txt | tr "\n" " "`
do
    echo "averaging coverage for Sample_$i"
    SAMPLE_DIR=$MAPPED_DIR/Sample_$i
    THIS_COVERAGE=$SAMPLE_DIR/$COVERAGE
    THIS_COVERAGE_AVG=$THIS_COVERAGE.$AVG_INTERVAL
#    if [ ! -e $THIS_COVERAGE_AVG ]; then
#	qsub -e $THIS_COVERAGE_AVG.qsub.err \
#	    -o $THIS_COVERAGE_AVG.qsub.out \
	/share/apps/bin/R --slave --args \
	    $THIS_COVERAGE $AVG_INTERVAL \
	    $THIS_COVERAGE_AVG $AVG_METHOD \
	    <$SCRIPTS_DIR/average_coverage.R 	    
#    fi
done

for GROUP in `cut -f2 /home/obot/single-cell/untreated_\
highdose_survivors.txt | uniq | tr "\n" " "`
do
  THIS_GROUP_DIR=$MAPPED_DIR/treatment_groups/$GROUP
  for i in 1 2
  do
      SAMPLE_DIR=$THIS_GROUP_DIR/$GROUP$i
    THIS_COVERAGE=$SAMPLE_DIR/$COVERAGE
    THIS_COVERAGE_AVG=$THIS_COVERAGE.$AVG_INTERVAL
#    if [ ! -e $THIS_COVERAGE_AVG ]; then
#	qsub -e $THIS_COVERAGE_AVG.qsub.err \
#	    -o $THIS_COVERAGE_AVG.qsub.out \
	/share/apps/bin/R --slave --args \
	    $THIS_COVERAGE $AVG_INTERVAL \
	    $THIS_COVERAGE_AVG $AVG_METHOD \
	    <$SCRIPTS_DIR/average_coverage.R 

#    fi      
  done
done