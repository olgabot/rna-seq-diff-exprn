#!/bin/sh

# To run:
# /home/obot/single-cell/scripts/do-tophat.sh

HOME_DIR=/home/obot/single-cell

# Map reads from individual cells
for i in `head -n 1 $HOME_DIR/untreated_highdose_survivors.txt | cut -f1 - | tr "\n" " "`
do
  for RUN in `ls $HOME_DIR/unmapped`
    do
    IN_PREFIX=$HOME_DIR/unmapped/$RUN/Sample_$i/trimmed_$i
    OUT_DIR=$HOME_DIR/mapped/Sample_$i/$RUN
    OUT_PREFIX=$OUT_DIR/Sample_$i/$RUN/$RUN
    rm $OUT_DIR/qsub.err; rm $OUT_DIR/qsub.out   # Start with fresh error/outputs every time

  # READ_ONE is a filename of trimmed reads of the first pair for that sample
  # Same goes for READ_TWO, the corresponding filename of trimmed reads of the 2nd pair
  # --- These will be used as input into the ---
  # --- Tophat mapping script                ---
    READ_ONE=`ls $IN_PREFIX\_*_R1\_*trim9.fastq | tr "\n" "," | sed "s/,$//"`
#  do
    READ_TWO=`echo $READ_ONE | sed "s/_R1_/_R2_/g"`
# Make the output directory if it's not there already
    if [ ! -d $OUT_DIR ]; then
	mkdir $OUT_DIR
    fi
    qsub -e $OUT_DIR/qsub.err -o $OUT_DIR/qsub.out \
	-v "TMPDIR=$OUT_DIR/tmp-qsub" \
	/home/obot/single-cell/scripts/tophat.sh \
	$READ_ONE $READ_TWO \
	$OUT_PREFIX $OUT_DIR $OUT_PREFIX\_read_stats.txt
  done
done

# To merge reads: merge bowtie/tophat files

# Map reads from merged cell conditions
#for condition in `cut -f2 README.txt | uniq | tr "\n" " "`
#do
#  if [-d /home/obot/single-cell/mapped/$condition/ ]
#      mkdir /home/obot/single-cell/mapped/$condition/
#  fi
#  for i in `grep $condition README.txt | cut -f1 | tr "\n" " "`
#  do    
#  done
#done