#!/bin/sh
# To run:
# /home/obot/single-cell/scripts/do-map-reads.sh
HOME_DIR=/home/obot/single-cell

# Map reads from individual cells
# Do all but the first line in the README file
for i in `sed 1d $HOME_DIR/README.txt | cut -f1 - | tr "\n" " "`
do
  for RUN in `ls $HOME_DIR/unmapped`
    do
    IN_PREFIX=$HOME_DIR/unmapped/$RUN/Sample_$i/trimmed_$i
    OUT_DIR=$HOME_DIR/mapped/Sample_$i/$RUN
    OUT_PREFIX=$OUT_DIR/$RUN

    READ_ONE=`ls $IN_PREFIX\_*_R1\_*trim9.fastq | tr "\n" "," | sed "s/,$//"`
    READ_TWO=`echo $READ_ONE | sed "s/_R1_/_R2_/g"`
#  echo 'Sample ' $i
  # READS_ONE is all the filenames of the reads of the first pair for that sample, separated by spaces
  # Same goes for READS_TWO, all the filenames  reads of the 2nd pair, separated by spaces
  # --- These will be used as input into trimming ---
  # --- off the first 9 bp from the 5' end        ---
#  READS_ONE=`ls --ignore=*trim9* /home/obot/single-cell/unmapped/run*/Sample_$i/trimmed_$i\_*_R1\_*.fastq | tr "\n" " "`
#  READS_TWO=`echo $READS_ONE | sed "s/_R1_/_R2_/g"`

  # TRIM_READS_ONE is all the filenames of the reads of the
  # first pair for that sample, separated by a COMMA: ","
  # Same goes for TRIM_READS_TWO, all the filenames of the 2nd
  # pair of reads, separated by COMMA
  # --- These will be used as input ---
  # --- for the map-reads script    ---
#  TRIM_READS_ONE=`echo $READS_ONE | sed 's:\.fastq:_trim9\.fastq:' | tr " " "," | sed "s/,$//"`
#  TRIM_READS_TWO=`echo $READS_TWO | sed 's:\.fastq:_trim9\.fastq:' | tr " " "," | sed "s/,$//"`
    if [ ! -d $OUT_DIR ]; then
	mkdir $OUT_DIR
    fi
#  for R1 in $READS_ONE
#    do
	# TRIM_R1 and _R1 are the filenames used to save 
	# the output of trimming the first 9 base pairs
	# from each of READ_ONE and READ_TWO
#    R2=`echo $R1 | sed "s/_R1_/_R2_/g"`
#    TRIM_R1=`echo $R1 | sed 's:\.fastq:_trim9\.fastq:'`
#    TRIM_R2=`echo $R2 | sed 's:\.fastq:_trim9\.fastq:'`
    # Wendy's magic script to cut off the first 9 base pairs
#    cat $R1 | sed 'n;s/^.\{9\}//g' > $TRIM_R1
#    cat $R2 | sed 'n;s/^.\{9\}//g' > $TRIM_R2
#    rm $R1; rm $R2
#  done
#  echo 'TRIM_READS_ONE: ' $TRIM_READS_ONE
#  echo 'TRIM_READS_TWO: ' $TRIM_READS_TWO
#  OUTPUT_DIR=/home/obot/single-cell/mapped/Sample_$i
    if [ -e $OUT_DIR/qsub.err ]; then
        rm $OUT_DIR/qsub.err; rm $OUT_DIR/qsub.out   # Start with fresh error/outputs every time
    fi
    qsub -e $OUT_DIR/qsub.err -o $OUT_DIR/qsub.out \
	-v "TMPDIR=$OUT_DIR/tmp-qsub" \
	$HOME_DIR/scripts/map-reads.sh \
	$READ_ONE $READ_TWO \
	$OUT_PREFIX $OUT_DIR
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