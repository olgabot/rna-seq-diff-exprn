#!/bin/sh

#$ -q all.q
#$ -cwd
#$ -u obot
#$ -pe mpi 4
#$ -S /bin/sh

# Example run: qsub -o test/test.out -e test/test.err tophat.sh run2/Sample_3/trimmed_3_TAGCTT_L006_R1_001.fastq run2/Sample_3/trimmed_3_TAGCTT_L006_R2_001.fastq test/test test test/test_read_stats.txt
echo
echo "Job PID:" $job_pid
echo "arguments:" $*

READ_ONE=$1
READ_TWO=$2
BAM_PREFIX=$3
#SAM=$BAM_PREFIX.sam
#BAM=$BAM_PREFIX.bam
#TRIM_BAM=$BAM_PREFIX\_trim9.bam
#UNALIGNED=$BAM_PREFIX\_unaligned.fastq
#SINO_BAM=$BAM_PREFIX\_sino.bam
#SINO_TRIM_BAM=$BAM_PREFIX\_sino_trim9.bam
OUT_DIR=$4
READ_STATS=$5

# Gene models file
GTF=/home/pvcastro/reference_known_genes.gtf

#echo 'GTF:' $GTF

export BOWTIE_INDEXES=~hyjkim/reference/hg19/bowtie/

#echo 'Bowtie indexes:' $BOWTIE_INDEXES

echo "READ_ONE: " $READ_ONE
echo "READ_TWO: " $READ_TWO
PATH=$PATH:/share/apps/bowtie-0.12.7:/share/apps/samtools-0.1.18

INSERT_SIZE_MEAN=`sed -n '2p;d' $READ_STATS`
INSERT_SIZE_SD=`sed -n '1p;d' $READ_STATS`

echo 'read stats - mean:' $INSERT_SIZE_MEAN '  std dev:' $INSERT_SIZE_SD

# Run TopHat to determine splice sites de novo
/share/apps/tophat-1.4.1.Linux_x86_64/tophat --num-threads 4 \
    -G $GTF --mate-inner-dist $INSERT_SIZE_MEAN \
    --mate-std-dev $INSERT_SIZE_SD \
    --segment-length 20 \
    --library-type fr-unstranded -o $OUT_DIR \
    hg19 $READ_ONE $READ_TWO

#echo 'unaligned:' $UNALIGNED
# Take the unmapped reads, trim the first 9 base pairs and map
#/share/apps/bowtie-0.12.7/bowtie -u 100000 --chunkmbs 256 -S -m 1 --maxins 500 -p4 --phred64-quals -ul $UNALIGNED -1 $READ_ONE -2 $READ_TWO hg19 | /share/apps/samtools-0.1.18/samtools view -F 4 -S -b - > $BAM

#/share/apps/bowtie-0.12.7/bowtie -5 9 --chunkmbs 256 -S -m 1 --maxins 500 -p4 --phred64-quals -1 $READ_ONE -2 $READ_TWO hg19 -ul | /share/apps/samtools-0.1.18/samtools view -F 4 -S -b - >$TRIM_BAM

# Do the intersect
#/share/apps/BEDTools-Version-2.15.0/bin/bedtools intersect -abam $BAM -b /home/obot/hg19_ensembl-genes_single-exon.bed > $SINO_BAM
#/share/apps/BEDTools-Version-2.15.0/bin/bedtools intersect -abam $TRIM_BAM -b ~/hg19_ensembl-genes_single-exon.bed > $SINO_TRIM_BAM



# Gene expression estimation
#GENE_COUNTS=$BAM_PREFIX\_gene_counts_sino.txt
#/share/apps/BEDTools-Version-2.15.0/bin/bedtools coverage -abam $SINO_BAM -b /home/obot/hg19_ensembl-genes_single-exon.bed > $GENE_COUNTS

# Randomly sample reads
# your script
#READ_DIST_R=read_distribution.R
#READ_STATS=$BAM_PREFIX\_readlengths_stats.txt
#/share/apps/bin/R < ~/get_read_distribution.R --vanilla --args $READ_LENGTHS $READ_STATS
#/home/obot/get_read_distribution.R $READ_LENGTHS $READ_STATS
#echo "read_lengths = read.delim("$READ_LENGTHS", header=FALSE)" >>$READ_DIST_R
#echo "mu = mean(abs(read_lengths[,1]))" >>$READ_DIST_R
#echo "sigma = sd(abs(read_lengths[,1]))" >>$READ_DIST_R
#echo "pdf(paste("$READ_LENGTHS", '_read_lengths.pdf'))" >>$READ_DIST_R
#echo "hist(abs(as.numeric(read_lengths[,1]), breaks=dim(read_lengths)[1]/10," >>$READ_DIST_R
#echo "main=paste("$READ_LENGTHS", 'read length distribution\nmu:'," >>$READ_DIST_R
#echo "signif(mu), ' sigma:', sigma), xlab='Read length')" >>$READ_DIST_R
#echo "dev.off()" >>$READ_DIST_R
#echo "write(paste('mu:', mu,'\nsigma:',sigma),file="$READ_STATS")" >>$READ_DIST_R

#RANDOM_READS=$BAM_PREFIX\_simulated_
#/home/obot/simulate_read_coverage.py --fasta /home/hyjkim/reference/mm9/genome/mm9.validated.fa --bed mm9_ensembl-genes_single-exon.bed --read-stats $READ_STATS -o $RANDOM_READS