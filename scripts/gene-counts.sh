#!/bin/sh

# Example run:
### Old version! do not use! ./gene-counts.sh /home/obot/single-cell/mapped/Sample_3/run2
# /home/obot/single-cell/scripts/gene-counts.sh /home/obot/single-cell/data/Sample_3/tophat_6n6_trimmed_3_merged_rmdup /home/obot/single-cell/data/Sample_3 3


#echo `ls $1/*.bam`

BAM_PREFIX=$1 #$OUT_DIR/tophat_6n6_trimmed_$i\_merged_rmdup
OUT_DIR=$2

echo "finding gene counts for" $BAM_PREFIX.bam

#echo 'OUT_DIR:' $OUT_DIR
RMDUP_BAM=$BAM_PREFIX.bam

PATH=$PATH:/share/apps/samtools-0.1.18:/share/apps/BEDTools-Version-2.15.0/
export PATH

SCRIPTS_DIR=/home/obot/single-cell/scripts

SINO_BED=/home/obot/bed_and_gtf/hg19_ucsc-genes_single-exon.bed
REF_BED=/home/obot/bed_and_gtf/hg19_ucsc-genes.bed

# Do the intersect
SINO_BAM=$BAM_PREFIX\_sino.bam
REF_BAM=$BAM_PREFIX\_reference.bam
if [ ! -e $SINO_BAM ]; then
# -F 256 : filters out reads that are "not primary alignment"
    samtools view -b $RMDUP_BAM | \
	bedtools intersect -abam - \
	-b $SINO_BED >$SINO_BAM
fi

if [ ! -e $REF_BAM ]; then
    samtools view -b $RMDUP_BAM | \
	bedtools intersect -abam - \
	-b $REF_BED >$REF_BAM   
fi

# Gene expression estimation via bedtools
# -counts: Only report the count of overlaps, 
#          don't compute fraction, etc.
# -hist: Report a histogram of coverage for each feature in B
#        as well as a summary histogram for _all_ 
#        features in B.
#        Output (tab delimited) after each feature in B:
#           1) depth
#           2) # bases at depth
#           3) size of B
#           4) % of B at depth
# -d: Report the depth at each position in each B feature.
#     Positions reported are one based.  Each position
#     and depth follow the complete B feature.
coverageBed_options='-counts -hist -d'
SINO_BEDTOOLS_PREFIX=$OUT_DIR/sino_bedtools_counts
bedtools coverage \
    $coverageBed_options \
    -abam $SINO_BAM \
    -b $SINO_BED >$SINO_BEDTOOLS_PREFIX.txt
REF_BEDTOOLS_PREFIX=$OUT_DIR/ref_bedtools_counts
bedtools coverage \
    $coverageBed_options \
    -abam $REF_BAM \
    -b $REF_BED >$REF_BEDTOOLS_PREFIX.txt

# Turn bedtools counts files into circos data
awk '{ print $1"\t"$2"\t"$3"\t"$13"\tgene="$4}' \
    $SINO_BEDTOOLS_PREFIX.txt > $SINO_BEDTOOLS_PREFIX.circos
awk '{ print $1"\t"$2"\t"$3"\t"$13"\tgene="$4}' \
    $REF_BEDTOOLS_PREFIX.txt > $REF_BEDTOOLS_PREFIX.circos

# Gene Expression estimation via HTSeq
############### For some reason this does NOT work when
############### executed from this script but works fine 
############### from the command line
PATH=$PATH:/usr/local/bin/python2.7 ; export PATH

SINO_GFF=/home/obot/bed_and_gtf/hg19_ucsc-genes_single-exon.gtf
REF_GFF=/home/obot/bed_and_gtf/hg19_ucsc-genes.gtf

HTSeqCount_options='--stranded=no'

SINO_HTSEQ_PREFIX=$OUT_DIR/sino_htseq_counts
export LC_ALL=POSIX # Forces char encoding to be POSIX (maybe?)
# 'sort -s -k 1,1': Sort SAM file by read name before HTSeq
RMDUP_READ_SORTED=$BAM_PREFIX\_sorted_read_name.sam
samtools view $RMDUP_BAM | \
    sort -s -k 1,1 >$RMDUP_READ_SORTED

#samtools view $SINO_BAM | \
#    sort -s -k 1,1 | \
    /usr/local/bin/htseq-count \
    $HTSeqCount_options $RMDUP_READ_SORTED $SINO_GFF \
    >$SINO_HTSEQ_PREFIX.txt\
    2>$SINO_HTSEQ_PREFIX.err


# one line:
# SINO_BAM=tophat_6n6_trimmed_3_merged_rmdup_sino.bam; HTSeqCount_options='--stranded=no'; SINO_HTSEQ_PREFIX=sino_htseq_counts; SINO_GFF=/home/obot/bed_and_gtf/hg19_ucsc-genes_single-exon.gtf; samtools view $SINO_BAM | sort -s -k 1,1 | /usr/local/bin/htseq-count $HTSeqCount_options - $SINO_GFF >$SINO_HTSEQ_PREFIX.txt 2>$SINO_HTSEQ_PREFIX.err

#REF_HTSEQ=$OUT_DIR/ref_htseq_counts
# 'sort -s -k 1,1': Sort SAM file by read name before HTSeq
#samtools view $REF_BAM | \
#    sort -s -k 1,1 | \
    /usr/local/bin/htseq-count \
    $HTSeqCount_options $RMDUP_READ_SORTED $REF_GFF \
    >$REF_HTSEQ.txt \
    2>$REF_HTSEQ.err
# one line:
# REF_GFF=/home/obot/bed_and_gtf/hg19_ucsc-genes.gtf ; OUT_DIR=/home/obot/single-cell/data/Sample_4 ; REF_HTSEQ=$OUT_DIR/ref_htseq_counts ; RMDUP_READ_SORTED=$OUT_DIR/tophat_6n6_trimmed_4_merged_rmdup_sorted_read_name.sam ; /usr/local/bin/htseq-count --stranded=no $RMDUP_READ_SORTED $REF_GFF >$REF_HTSEQ.txt 2>$REF_HTSEQ.err


# Make circos data from HTSeq output
SINO_SORTED_BED=/home/obot/bed_and_gtf/hg19_ucsc-genes_single-exon_sorted_featurename.bed
REF_SORTED_BED=/home/obot/bed_and_gtf/hg19_ucsc-genes_sorted_featurename.bed
TMP=tmp.htseq
awk '{ print $2"\tgene_htseq="$1 }' $SINO_HTSEQ.txt > $TMP
paste $SINO_SORTED_BED $TMP | \
    awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4","$6 }' \
    >$SINO_HTSEQ.circos
awk '{ print $2"\tgene_htseq="$1 }' $REF_HTSEQ.txt > $TMP
paste $REF_SORTED_BED $TMP | \
    awk '{ print $1"\t"$2"\t"$3"\t"$5"\t"$4","$6 }' \
    >$REF_HTSEQ.circos

# Genome-wide coverage via bedtools
GENOME=/share/apps/BEDTools-Version-2.15.0/genomes/human.hg19.genome
GENOME_COVERAGE_BED=$OUT_DIR/genome_coverage_bedtools
# -d:depth value -bg:bedGraph format
genomeCoverageBed -bg \
    -ibam $RMDUP_BAM -g $GENOME \
    >$GENOME_COVERAGE_BED.txt
# Rename to actual file name (before it was just a prefix)
GENOME_COVERAGE_BED=$GENOME_COVERAGE_BED.txt


# Genome-wide coverage via HTSeq
PATH=$PATH:/usr/local/bin/python2.7 ; export PATH
GENOME_COVERAGE_HTSEQ=$OUT_DIR/genome_coverage_htseq
# 'sort -s -k 1,1': Sort SAM file by read name before HTSeq
$SCRIPTS_DIR/get_transcript_coverage.py \
    -asam $RMDUP_READ_SORTED -o $GENOME_COVERAGE_HTSEQ
rm $RMDUP_READ_SORTED

$SCRIPTS_DIR/circos-and-rseqc.sh $OUT_DIR $RMDUP_BAM

