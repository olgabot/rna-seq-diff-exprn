#!/bin/sh


SINGLE_CELL_DIR=/home/obot/single-cell
IN_DIR=/home/wlee/RNASeq/single-cell/Salk_Taxol_BamFiles

for i in `sed 1,12d \
$SINGLE_CELL_DIR/untreated_highdose_survivors.txt | \
cut -f1 - | tr "\n" " "`
do
#  for RUN in `ls $HOME_DIR/unmapped`
#    do
  OUT_DIR=$SINGLE_CELL_DIR/data/Sample_$i
  if [ ! -d $OUT_DIR ]; then
      mkdir -p $OUT_DIR
  fi  
  IN_BAM=$IN_DIR/tophat_6n6_trimmed_$i\_merged.bam
#  if [ -e $OUT_DIR/accepted_hits.bam ]; then
  RMDUP_PREFIX=$OUT_DIR/tophat_6n6_trimmed_$i\_merged_rmdup
  if [ -e $OUT_DIR/qsub.err ]; then
      rm $OUT_DIR/qsub.err; rm $OUT_DIR/qsub.out   # Start with fresh error/outputs every time
  fi
#  qsub -e $OUT_DIR/qsub.err -o $OUT_DIR/qsub.out \
#      -l mppmem=8G \
      $SINGLE_CELL_DIR/scripts/rmdup_gene-counts_wlee.sh \
	  $OUT_DIR $IN_BAM $RMDUP_PREFIX
done