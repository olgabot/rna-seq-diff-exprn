#!/bin/sh

# Globally used variables
SCRIPTS_DIR=$SCRIPTS_DIR
CIRCOS_BIN=$CIRCOS_BIN
GENOME=/share/apps/BEDTools-Version-2.15.0/genomes/human.hg19.genome

# File that has UCSC transcript IDs in column 1
# and gene symbol in column 2
UCSC_SYMBOL=/share/reference/hg19/genes/hg19_ucsc-genes_symbol_no-header.tab
# For DEXSeq
DEXSEQ_GFF=/home/obot/bed_and_gtf/hg19_ucsc-genes_dexseq.gff
BASE_OUT_DIR=$BASE_OUT_DIR
COND=$COND
