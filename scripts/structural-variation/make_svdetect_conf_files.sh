#!/bin/sh

source /home/obot/single-cell/scripts/set_global_vars.sh

OLD_ID=#ID#
OUT_DIR=$SV_DIR/SVDetect
TEMPLATE=$OUT_DIR/template.sv.conf

#echo "SAMPLE_IDS:" $SAMPLE_IDS

for i in $SAMPLE_IDS
do
#    echo Sample_$i
    sed 's/'"$OLD_ID"'/'"$i"'/g' $TEMPLATE > $OUT_DIR/Sample_$i/sample_$i.sv.conf
done