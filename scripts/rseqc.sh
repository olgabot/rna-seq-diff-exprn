#!/bin/sh

####################### !!! ##########################
# Make sure you add these to your path before running:
# export PYTHONPATH=/share/apps/RSeQC/usr/local/lib/python2.7/site-packages:$PYTHONPATH
# export PATH=/share/apps/RSeQC/usr/local/bin:$PATH
####################### !!! ##########################

# Example run:
# rseqc.sh my_favorite_sample.bam rseqc /share/reference/hg19/genes/hg19_ucsc-genes.bed /share/reference/hg19/genes/hg19_ucsc-genes.gtf
##
# Example using bash scripting:
# rseqc.sh $BAM $OUT_DIR $BED

# things needed by SGE
#$ -cwd
#$ -q big.q
#$ -S /bin/bash
export PYTHONPATH=/share/apps/RSeQC/usr/local/lib/python2.7/site-packages:$PYTHONPATH
export PATH=/share/apps/RSeQC/usr/local/bin:/share/apps/python-2.7/bin/:/share/apps/bin:$PATH
RSEQC_BIN=/share/apps/RSeQC/usr/local/bin


# BAM file that you want to do the analysis on
BAM=$1

# Output directory where you want all the results to go
# Files will be prepended with "$OUT_DIR/"
OUT_DIR=$2

# BED files of which you want to see the gene body coverage
BED=$3

# Make plots via PREFIX
PREFIX=$OUT_DIR/rseqc

echo 'executing: clipping_profile.py' ; python2.7 $RSEQC_BIN/clipping_profile.py --input-file $BAM --out-prefix $PREFIX
echo 'executing: geneBody_coverage.py' ; python2.7 $RSEQC_BIN/geneBody_coverage.py --input-file $BAM --refgene $BED --out-prefix $PREFIX
echo 'executing: inner_distance.py' ; python2.7 $RSEQC_BIN/inner_distance.py --input-file $BAM --refgene $BED --out-prefix $PREFIX
echo 'executing: junction_annotation.py' ; python2.7 $RSEQC_BIN/junction_annotation.py --input-file $BAM --refgene $BED --out-prefix $PREFIX
echo 'executing: junction_saturation.py' ; python2.7 $RSEQC_BIN/junction_saturation.py --input-file $BAM --refgene $BED --out-prefix $PREFIX
echo 'executing: read_distribution.py' ; python2.7 $RSEQC_BIN/read_distribution.py - $BAM -r $BED  > $PREFIX/read_distribution.txt
echo 'executing: read_duplication.py' ; python2.7 $RSEQC_BIN/read_duplication.py --input-file $BAM --out-prefix $PREFIX
echo 'executing: read_GC.py' ; python2.7 $RSEQC_BIN/read_GC.py --input-file $BAM --out-prefix $PREFIX
echo 'executing: read_NVC.py' ; python2.7 $RSEQC_BIN/read_NVC.py --input-file $BAM --out-prefix $PREFIX
echo 'executing: read_quality.py' ; python2.7 $RSEQC_BIN/read_quality.py --input-file $BAM --out-prefix $PREFIX
echo 'executing: RPKM_saturation.py' ; python2.7 $RSEQC_BIN/RPKM_saturation.py --input-file $BAM --out-prefix $PREFIX --refgene $BED 
echo 'executing: RPKM_count.py' ; python2.7 $RSEQC_BIN/RPKM_count.py --input-file $BAM --out-prefix $PREFIX --refgene $BED
