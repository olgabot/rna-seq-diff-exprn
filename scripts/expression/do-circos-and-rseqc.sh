#!/bin/sh

SCRIPTS_DIR=/home/obot/single-cell/scripts
BASE_DIR=/home/obot/single-cell/data

for i in 3 6 7 8 9 #62 63 64 67 68
do
    OUT_DIR=$BASE_DIR/Sample_$i
    IN_BAM=$OUT_DIR/tophat_6n6_trimmed_$i\_merged_rmdup.bam
    $SCRIPTS_DIR/circos-and-rseqc.sh $OUT_DIR $IN_BAM \
	2>$OUT_DIR/circos-and-rseqc.sh.err
done

#GROUPS_DIR=$BASE_DIR/treatment_groups

#for GROUP in `cut -f2 /home/obot/single-cell/untreated_\
#highdose_survivors.txt | uniq | tr "\n" " "`
#do
#  echo "finding gene counts for:" $GROUP
#  THIS_GROUP_DIR=$GROUPS_DIR/$GROUP
#  echo "THIS_GROUP_DIR:" $THIS_GROUP_DIR

#  GROUP1_DIR=$THIS_GROUP_DIR/$GROUP\1
#  GROUP2_DIR=$THIS_GROUP_DIR/$GROUP\2
#  if [ ! -d $GROUP1_DIR ] ; then
#      mkdir -p $GROUP1_DIR
#      mkdir -p $GROUP2_DIR
#  fi
#  echo "finding samples...."
#  SAMPLES=`grep $GROUP $SINGLE_CELL_DIR/untreated_highdose_survivors.txt | cut -f1 | tr "\n" " "`
#  echo "samples:" $SAMPLES
#  GROUP1=`echo $SAMPLES | awk -F" " ' { print $1" "$2" "$3 } '`
  # for GROUP2: if there's only 5, $6 will be empty
#  GROUP2=`echo $SAMPLES | awk -F" " ' { print $4" "$5" "$6 } '`

#  BAM1_OUT_PREFIX=`echo $GROUP1 | tr " " "-" | awk -F" " -v DIR=$GROUP1_DIR -v GROUP=$GROUP ' { print DIR"/"GROUP"1_"$1 } '`
#  BAM2_OUT_PREFIX=`echo $GROUP2 | tr " " "-" | awk -F" " -v DIR=$GROUP2_DIR -v GROUP=$GROUP ' { print DIR"/"GROUP"2_"$1 } '`

#  BAM1_OUT=$BAM1_OUT_PREFIX.bam
#  BAM2_OUT=$BAM2_OUT_PREFIX.bam

#  $SCRIPTS_DIR/circos-and-rseqc.sh $GROUP1_DIR $BAM1_OUT \
#      2>$GROUP1_DIR/circos-and-rseq.sh.err
#  $SCRIPTS_DIR/circos-and-rseqc.sh $GROUP2_DIR $BAM2_OUT \
#      2>$GROUP2_DIR/circos-and-rseq.sh.err
#done