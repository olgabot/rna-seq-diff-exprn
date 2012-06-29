#/bin/sh
# To run:
# /home/obot/single-cell/scripts/remove-pcr-duplicates.sh
SINGLE_CELL_DIR=/home/obot/single-cell
PATH=$PATH:/share/apps/samtools-0.1.18
BAM_DIR=/home/wlee/RNASeq/single-cell/Salk_Taxol_BamFiles

for i in `cat $HOME_DIR/untreated_highdose_survivors.txt | cut -f1 - | tr "\n" " "`
do
#  for RUN in `ls $HOME_DIR/unmapped`
#    do
  OUT_DIR=$SINGLE_CELL_DIR/data/Sample_$i
  IN_BAM=$BAM_DIR/tophat_6n6_trimmed_$i\_merged.bam
#  if [ -e $OUT_DIR/accepted_hits.bam ]; then
  RMDUP_BAM=$OUT_DIR/tophat_6n6_trimmed_$i\_merged_rmdup.bam
  if [ ! -e $RMDUP_BAM ]; then
      echo 'processing:' $IN_BAM
      samtools rmdup $IN_BAM $RMDUP_BAM >$OUT_DIR/rmdup.log
      $RMDUP_SORTED=$OUT_DIR/accepted_hits_rmdup_sorted
      samtools sort $RMDUP_BAM $RMDUP_SORTED
      samtools index $RMDUP_SORTED.bam >$OUT_DIR/index.log
      STATS=$OUT_DIR/idxstats.tab
      samtools idxstats $RMDUP_SORTED >$STATS
#	MAPPED=$OUT_DIR/idxstats.mapped
	    # Remove the mitochondrial DNA because its circular chromosome gets amplified like whoa
      MAPPED=`cat $STATS | sed '/^chrM/d' | awk 'BEGIN { sum=0 } { sum+=$3 } END { print sum }'`
#	    echo 'mapped:' $MAPPED
      STATS_FRACS=$OUT_DIR/idxstats_fractions.tab
      cat $STATS | awk '{print $3/$2}' > $STATS_FRACS
#	    QSUB_ERR=$OUT_DIR/qsub.err
      RIGHT_INFO=$OUT_DIR/right_kept_reads.info
      LEFT_INFO=$OUT_DIR/left_kept_reads.info
      RIGHT_UNMAPPED=`grep 'reads_in =' $RIGHT_INFO | awk '{print(substr($2,index($2, "=")+1, length($2)))}'`
      LEFT_UNMAPPED=`grep 'reads_in =' $LEFT_INFO | awk '{print(substr($2,index($2, "=")+1, length($2)))}'`
      echo 'right unmapped:' $RIGHT_UNMAPPED
      echo 'left unmapped:' $LEFT_UNMAPPED
      MAPPED_STATS=$OUT_DIR/mapping_stats.tab
      echo $RIGHT_UNMAPPED $LEFT_UNMAPPED $MAPPED | awk '{print "#mapped\t#unmapped\t%mapped"} \
{ UNMAPPED+=$1+$2 }  { print $3"\t"UNMAPPED"\t"$3/UNMAPPED }' >$MAPPED_STATS
  fi
	# Estimate gene expression counts
  $HOME_DIR/scripts/gene-counts.sh $OUT_DIR
 # fi
 # done
done