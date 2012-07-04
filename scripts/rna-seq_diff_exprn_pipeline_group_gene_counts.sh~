#!/bin/sh
SINGLE_CELL_DIR=/home/obot/single-cell
MAPPED_DIR=$SINGLE_CELL_DIR/data
GROUPS_DIR=$MAPPED_DIR/treatment_groups
SCRIPTS_DIR=$SINGLE_CELL_DIR/scripts

echo "MAPPED_DIR:" $MAPPED_DIR

for GROUP in `cut -f2 /home/obot/single-cell/untreated_\
highdose_survivors.txt | uniq | sed 1d | tr "\n" " "`
do
  echo "finding gene counts for:" $GROUP
  THIS_GROUP_DIR=$GROUPS_DIR/$GROUP
  echo "THIS_GROUP_DIR:" $THIS_GROUP_DIR

  GROUP1_DIR=$THIS_GROUP_DIR/$GROUP\1
  GROUP2_DIR=$THIS_GROUP_DIR/$GROUP\2
  if [ ! -d $GROUP1_DIR ] ; then
      mkdir -p $GROUP1_DIR
      mkdir -p $GROUP2_DIR
  fi
  echo "finding samples...."
  SAMPLES=`grep $GROUP $SINGLE_CELL_DIR/untreated_highdose_survivors.txt | cut -f1 | tr "\n" " "`
  echo "samples:" $SAMPLES
  GROUP1=`echo $SAMPLES | awk -F" " ' { print $1" "$2" "$3 } '`
  # for GROUP2: if there's only 5, $6 will be empty
  GROUP2=`echo $SAMPLES | awk -F" " ' { print $4" "$5" "$6 } '`

#  echo "MAPPED_DIR:" $MAPPED_DIR
  export MAPPED_DIR

  BAM1_IN=`echo $GROUP1 | tr " " "\n" | awk -F" " -v DIR=$MAPPED_DIR ' { print DIR"/Sample_"$1"/tophat_6n6_trimmed_"$1"_merged_rmdup.bam" } END { } ' | tr "\n" " "`
  BAM2_IN=`echo $GROUP2 | tr " " "\n" | awk -F" " -v DIR=$MAPPED_DIR ' { print DIR"/Sample_"$1"/tophat_6n6_trimmed_"$1"_merged_rmdup.bam" } END { } ' | tr "\n" " "`

  BAM1_OUT_PREFIX=`echo $GROUP1 | tr " " "-" | awk -F" " -v DIR=$GROUP1_DIR -v GROUP=$GROUP ' { print DIR"/"GROUP"1_"$1 } '`
  BAM2_OUT_PREFIX=`echo $GROUP2 | tr " " "-" | awk -F" " -v DIR=$GROUP2_DIR -v GROUP=$GROUP ' { print DIR"/"GROUP"2_"$1 } '`

  BAM1_OUT=$BAM1_OUT_PREFIX.bam
  BAM2_OUT=$BAM2_OUT_PREFIX.bam

  if [ ! -e $BAM1_OUT ] ; then
      echo "merging:" $BAM1_IN " >" $BAM1_OUT
      samtools merge $BAM1_OUT $BAM1_IN
  fi
  if [ ! -e $BAM2_OUT ] ; then
      echo "merging:" $BAM2_IN " >" $BAM2_OUT
      samtools merge $BAM2_OUT $BAM2_IN
  fi

  $SCRIPTS_DIR/gene-counts.sh $BAM1_OUT_PREFIX $GROUP1_DIR\
      2>$GROUP1_DIR/gene-counts.sh.err
  $SCRIPTS_DIR/gene-counts.sh $BAM2_OUT_PREFIX $GROUP2_DIR\
      2>$GROUP2_DIR/gene-counts.sh.err
done