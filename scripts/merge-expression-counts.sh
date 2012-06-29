#!/bin/sh

# 1st argument is comma-delimited files 
# (col1 is gene, col2 is count, e.g.)
IN_FILES=$1
# 2nd argument is comma-delimited column names
COL_NAMES=$2
# 3rd argument is output file name
OUT_FILE=$3

# initialize the output file
FILE1=`echo $IN_FILES | tr "," "\n" | head -n 1`
cp $FILE1 $OUT_FILE

# skip the first line via sed 1d
for F in `echo $IN_FILES | \
tr "," "\n" | sed 1d | tr "\n" " " `
do
#    echo "in loop, getting column 2 from file:" $F
    cut -f 2 $F | paste $OUT_FILE - > tmp
    mv tmp $OUT_FILE
done

# Add column names to top
echo $COL_NAMES | tr "," "\t" | \
    cat - $OUT_FILE > tmp
mv tmp $OUT_FILE
