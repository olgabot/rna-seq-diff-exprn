#!/bin/sh


SINGLE_CELL_DIR=/home/obot/single-cell
MAPPED_DIR=$SINGLE_CELL_DIR/mapped
#awk '{ printf("\trun2\trun2\trun2\trun3\trun3\trun3\trun4\trun4\trun4\trun5\trun5\trun5") }'
for RUN in `ls $SINGLE_CELL_DIR/unmapped`; do echo -ne "\t"$RUN"_date\t"$RUN"_size\t"$RUN"_%mapped"; done; echo

while read LINE
  do 
  i=`echo $LINE | awk -F" " '{ print $1 }'`
  STATUS=`echo $LINE | awk -F" " '{ print $2 }'`
  echo -ne $i\_$STATUS
  for RUN in run2 run3 run4 run5
    do 
    DIR=$MAPPED_DIR/Sample_$i/$RUN
    if [ -e $DIR/accepted_hits.bam ]; then
	LS=`ls -lh $DIR/accepted_hits.bam`
	DATE=`echo $LS | awk -F' ' '{ print $6 $7 }'`
	SIZE=`echo $LS | awk -F' ' '{ print $5 }'`
	if [ -e $DIR/mapping_stats.tab ]; then
	    MAPPING_PERCENTAGE=`cat $DIR/mapping_stats.tab | awk -F'\t' '{ if ($3 != "%mapped") { print $3} }'`
	else
	    MAPPING_PERCENTAGE=""
	fi
	echo -ne "\t"$DATE"\t"$SIZE"\t"$MAPPING_PERCENTAGE
    else
	echo -ne "\t\t\t"
    fi
  done
  echo
# -rw-r--r-- 1 obot obot 17M Apr 22 20:21 /home/obot/single-cell/mapped/Sample_3/run2/accepted_hits.bam
done <$SINGLE_CELL_DIR/untreated_highdose_survivors.txt