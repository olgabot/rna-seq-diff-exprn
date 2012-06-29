#!/bin/sh

GFF=/home/obot/bed_and_gtf/hg19_ucsc-genes_dexseq.gff
MAPPED_DIR=/home/obot/single-cell/data
SINGLE_CELL_DIR=/home/obot/single-cell

for i in `sed 1d $SINGLE_CELL_DIR/untreated_highdose_survivors.txt | cut -f1 - | tr "\n" " "`
do
    OUT_DIR=$MAPPED_DIR/Sample_$i
    BAM=$OUT_DIR/tophat_6n6_trimmed_$i\_merged_rmdup.bam
    OUT_COUNTS=$OUT_DIR/dexseq_counts.txt
    echo "making DEXSeq files for" $BAM
    samtools view $BAM | sort -k 1,1 | dexseq_count.py \
	--paired=yes --stranded=no $GFF - $OUT_COUNTS
done

for GROUP in `cut -f2 /home/obot/single-cell/untreated_\
highdose_survivors.txt | uniq | tr "\n" " "`
do
  THIS_GROUP_DIR=$MAPPED_DIR/treatment_groups/$GROUP

  GROUP1_DIR=$THIS_GROUP_DIR/$GROUP\1
  GROUP2_DIR=$THIS_GROUP_DIR/$GROUP\2

  for i in 1 2
  do
      OUT_DIR=$THIS_GROUP_DIR/$GROUP$i

      # .bam comes before _ref.bam and _sino.bam so this
      # will get the real bam file
      BAM=`ls *.bam $OUT_DIR | head -n 1`
      echo "making DEXSeq files for" $BAM      
      OUT_COUNTS=$OUT_DIR/dexseq_counts.txt
      samtools view $BAM | sort -k 1,1 | dexseq_count.py \
	  --paired=yes --stranded=no $GFF - $OUT_COUNTS
  done
done