#!/bin/sh
#$ -q all.q
#$ -cwd
#$ -u obot
#$ -pe mpi 4
#$ -S /bin/sh

# Example run: rm test/qsub.out ; rm test/qsub.err ; qsub -o test/qsub.out -e test/qsub.err /home/obot/single-cell/scripts/rmdup_gene-counts_wlee.sh test /home/wlee/RNASeq/single-cell/Salk_Taxol_BamFiles/tophat_6n6_trimmed_6_merged.bam test/test_rmdup


SINGLE_CELL_DIR=/home/obot/single-cell
PATH=$PATH:/share/apps/samtools-0.1.18

# Original for loop
#for i in `cat $SINGLE_CELL_DIR/untreated_highdose_survivors.txt | cut -f1 - | tr "\n" " "`
# Skip the first line:
#for i in `sed 1d $SINGLE_CELL_DIR/untreated_highdose_survivors.txt | cut -f1 - | tr "\n" " "`
#do
#  for RUN in `ls $HOME_DIR/unmapped`
#    do
  OUT_DIR=$1 #$SINGLE_CELL_DIR/data/Sample_$i
  if [ ! -d $OUT_DIR ]; then
      mkdir -p $OUT_DIR
  fi  
  IN_BAM=$2 #$IN_DIR/tophat_6n6_trimmed_$i\_merged.bam
  RMDUP_PREFIX=$3 #$OUT_DIR/tophat_6n6_trimmed_$i\_merged_rmdup
  RMDUP_BAM=$RMDUP_PREFIX.bam
#  if [ ! -e $RMDUP_BAM ]; then
      echo 'processing:' $IN_BAM
      SORTED_PREFIX=$3\_sorted #$OUT_DIR/tophat_6n6_trimmed_$i\_merged_sorted
      samtools sort $IN_BAM $SORTED_PREFIX
      # 2>$OUT_DIR/index.log redirects stderr ("2") to index.log
      samtools rmdup $SORTED_PREFIX.bam \
	  $RMDUP_BAM 2>$OUT_DIR/rmdup.log
      # remove any reads mapping to the chromosome "*"
      NOSTAR_BAM=$RMDUP_PREFIX\_nostar.bam
      samtools view -h $RMDUP_BAM | \
	  awk -F" " '{ if ( $1 ~ "^@" ) print $0; \
else if ( $3 !~ "\*") print $0 } ' | \
	  samtools view -S -b - >$NOSTAR_BAM
      mv $NOSTAR_BAM $RMDUP_BAM
      samtools index $RMDUP_BAM 2>$OUT_DIR/index.log
      rm $SORTED_PREFIX.bam
 # fi
	# Estimate gene expression counts
  echo "estimating gene expression for:" $RMDUP_BAM
  $SINGLE_CELL_DIR/scripts/gene-counts.sh \
      $RMDUP_PREFIX $OUT_DIR 2>$OUT_DIR/gene-counts.sh.err
 # fi
 # done
#done