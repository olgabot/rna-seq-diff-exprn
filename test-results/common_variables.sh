#!/bin/sh

# Globally used variables
SCRIPTS_DIR='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/scripts'
CIRCOS_BIN='/usr/bin/circos/bin/circos'
BASE_OUT_DIR='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results'
COND='test-data/conditions_chr9.tab'
COND_WITHOUT_COMMENTS='test-data/conditions_chr9.tab.without_comments'
BAM_PREFIXES='test-data/GSM721117_mctp_20F0GAAXX_1_chr9_withheader,test-data/GSM721119_mctp_20F0GAAXX_2_chr9_withheader,test-data/GSM721118_mctp_20F0GAAXX_3_chr9_withheader,test-data/GSM721116_mctp_20F0GAAXX_4_chr9_withheader,test-data/GSM721123_mctp_30CYNAAXX_5_chr9_withheader,test-data/GSM721124_mctp_209ENAAXX_8_chr9_withheader'
IDS='LNCaP_1,LNCaP_2,LNCaP_3,LNCaP_4,PrEC_1,PrEC_2'
TREATMENT_GROUPS='LNCaP,PrEC'
GENDERS='male,male,male,male,male,male'
STRAND_SPECIFICITIES='not_strand_specific,not_strand_specific,not_strand_specific,not_strand_specific,not_strand_specific,not_strand_specific'

# Gene and species-specific variables
GTF='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-data/hg19_ucsc_genes.gtf'
DEXSEQ_GTF='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-data/hg19_ucsc_genes_chr9_dexseq.gtf'
BED='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-data/hg19_ucsc_genes.bed.without_track_name'
TXPTID_SYMBOL='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-data/hg19_id_symbol.txt.sorted'
GENOME='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-data/human.hg19.genome'
KARYOTYPE='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-data/karyotype/karyotype.human.hg19.txt'
GENE_DENSITY_FILE='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-data/hg19_gene_density_1e5bins.txt'
MAX_GENE_DENSITY='0.001770017700177'
GENE_DENSITY_MIN_VALUE_CHANGE='.0000177001'
GC_CONTENT_FILE='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-data/hg19_gc_content_circos_chr9.txt'

# Number of groups to make from the TREATMENT_GROUPS
NUM_GROUPS='1'

# More output directories
EXPRN_DIR='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression'
BEDTOOLS_DIR='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/bedtools'
HTSEQ_DIR='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/htseq'
DEXSEQ_DIR='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/dexseq'
DEXSEQ_OUT='dexseq_counts.txt'
CIRCOS_OUT_DIR='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/circos'
CIRCOS_TEMPLATE_DIR='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/scripts/circos-templates'

# Array variables
END=6
BAM_PREFIX_ARRAY=( test-data/GSM721117_mctp_20F0GAAXX_1_chr9_withheader test-data/GSM721119_mctp_20F0GAAXX_2_chr9_withheader test-data/GSM721118_mctp_20F0GAAXX_3_chr9_withheader test-data/GSM721116_mctp_20F0GAAXX_4_chr9_withheader test-data/GSM721123_mctp_30CYNAAXX_5_chr9_withheader test-data/GSM721124_mctp_209ENAAXX_8_chr9_withheader )
ID_ARRAY=( LNCaP_1 LNCaP_2 LNCaP_3 LNCaP_4 PrEC_1 PrEC_2 )
GROUPS_ARRAY=( LNCaP PrEC )
GENDER_ARRAY=( male male male male male male )
STRAND_ARRAY=( not_strand_specific not_strand_specific not_strand_specific not_strand_specific not_strand_specific not_strand_specific )

# Variables for gene counts and genome-wide coverage files
COUNTS_BED='bedtools_gene_counts.txt'
COUNTS_HTSEQ='htseq_gene_counts.txt'
COUNTS_HTSEQ_PREFIX='htseq_gene_counts'
CIRCOS_ALL_COLORS_ARRAY=( dark2-8-qual-1 dark2-8-qual-2 dark2-8-qual-3 dark2-8-qual-4 dark2-8-qual-5 dark2-8-qual-6 dark2-8-qual-7 dark2-8-qual-8 )
CIRCOS_ID_COLOR_ARRAY=( dark2-8-qual-1 dark2-8-qual-1 dark2-8-qual-1 dark2-8-qual-1 dark2-8-qual-2 dark2-8-qual-2 )
GC_CONTENT_COLOR='dark2-8-qual-7'
GENE_DENSITY_COLOR='dark2-8-qual-8'
UPPER_LIMITS=''
MIN_VALUE_CHANGES=''
COVERAGE_HTSEQ_PREFIX='htseq_genome_coverage'
COVERAGE_HTSEQ='htseq_genome_coverage.wig'
COVERAGE_BEDTOOLS_PREFIX='bedtools_genome_coverage'
COVERAGE_BEDTOOLS='bedtools_genome_coverage.txt'
HTSEQ_BIN='/opt/local/Library/Frameworks/Python.framework/Versions/2.7/bin/htseq-count'
TREATMENT_GROUPS_DIR='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/merged_groups'
GROUP_IDS='LNCaP_group1of1,PrEC_group1of1'
COUNTS_BED_TAB='bedtools_txptID_count.txt'
COUNTS_HTSEQ_TAB='htseq_txptID_count.txt'

# Variables for expression counts tables
BED_COUNTS_TABLE='bedtools_counts_table.tab'
HTSEQ_COUNTS_TABLE='htseq_counts_table.tab'
