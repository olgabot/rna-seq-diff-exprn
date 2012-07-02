#!/bin/sh

COMMON_VARS=$1
source $COMMON_VARS

for GROUP in $GROUPS ; do
  echo "finding gene counts for:" $GROUP
  THIS_GROUP_DIR=$GROUPS_DIR/$GROUP
  echo "THIS_GROUP_DIR:" $THIS_GROUP_DIR

  echo "finding samples...."
  SAMPLES=`grep $GROUP $COND | cut -f3 | tr "\n" " "`
  echo "samples:" $SAMPLES
  declare -a SAMPLE_ARRAY=( `echo $SAMPLES` )

  for (( i = 0 ; i < ${#GROUPS_ARRAY[@]} ; i++ ));  do
      
  done

  MIN_SAMPLES_PER_GROUP=`echo ${#SAMPLE_ARRAY[@]} $NUM_GROUPS | awk -F' ' '{ print $1/$2 }' | awk -F. '{print $1}'`

  for ((i=1;i<=$NUM_GROUPS;++i)); do
      GROUPi_DIR=$THIS_GROUP_DIR/$GROUP$i
      if [ ! -d $GROUPi_DIR ] ; then
	  mkdir -p $GROUPi_DIR
      fi
      
      ind1=`echo $i $MIN_SAMPLES_PER_GROUP | awk -F' ' '{ print 0+($1-1)*$2'}`
      if [ $i -lt $NUM_GROUPS ] ; then
	  ind2=`echo $i $MIN_SAMPLES_PER_GROUP | awk -F' ' '{ print 0+($1)*$2-1'}`
      else
	  ind2=`echo ${#SAMPLE_ARRAY[@]}`
      fi
      GROUPi_SAMPLES=`echo ${SAMPLE_ARRAY[@]ind1:ind2}`

      BAM_IN=`echo $GROUPi_SAMPLES | tr " " "\n" | awk -F" " -v DIR=$MAPPED_DIR ' { print DIR"/Sample_"$1"/tophat_6n6_trimmed_"$1"_merged_rmdup.bam" } END { } ' | tr "\n" " "`

      BAM_OUT_PREFIX=`echo $GROUP1 | tr " " "-" | awk -F" " -v DIR=$GROUP1_DIR -v GROUP=$GROUP ' { print DIR"/"GROUP"1_"$1 } '`
      BAM_OUT=$BAM_OUT_PREFIX.bam
  
      if [ ! -e $BAM1_OUT ] ; then
	  echo "merging:" $BAM1_IN " >" $BAM1_OUT
	  samtools merge $BAM1_OUT $BAM1_IN
      fi

      $SCRIPTS_DIR/gene-counts.sh \
	  $BAM1_OUT_PREFIX $GROUP1_DIR \
	  2>$GROUP1_DIR/gene-counts.sh.err
done