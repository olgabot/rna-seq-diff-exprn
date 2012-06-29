#!/bin/sh
for f in `ls /home/wlee/RNASeq/single-cell/Taxol_Fastq/*.tgz`
do
  tar -xvzf $f
done