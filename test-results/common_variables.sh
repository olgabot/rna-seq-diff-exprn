#!/bin/sh

# Globally used variables
SCRIPTS_DIR='/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/scripts'
CIRCOS_BIN='/usr/bin/circos/bin/circos'
BASE_OUT_DIR='test-results'
COND='test-data/conditions.tab'
BAM_PREFIXES='test-data/GSM721117_mctp_20F0GAAXX_1,test-data/GSM721119_mctp_20F0GAAXX_2,test-data/GSM721118_mctp_20F0GAAXX_3,test-data/GSM721116_mctp_20F0GAAXX_4,test-data/GSM721123_mctp_30CYNAAXX_5,test-data/GSM721124_mctp_209ENAAXX_8'
IDS='LNCaP_1,LNCaP_2,LNCaP_3,LNCaP_4,PrEC_1,PrEC_2'
TREATMENT_GROUPS='LNCaP,PrEC'
GENDERS='male,male,male,male,male,male'

# Gene and species-specific variables
GTF='test-data/hg19_ucsc_genes.gtf'
DEXSEQ_GTF='test-data/hg19_ucsc_genes_dexseq.gtf'
BED='test-data/hg19_ucsc_genes.bed'
TXPTID_SYMBOL='test-data/hg19_id_symbol.txt'
GENOME='test-data/human.hg19.genome'

# Number of groups to make from the TREATMENT_GROUPS
NUM_GROUPS='1'

# More output directories
EXPRN_DIR='test-results/expression'
CIRCOS_OUT_DIR='test-results/circos'

# Array variables
END=6
BAM_PREFIX_ARRAY=( test-data/GSM721117_mctp_20F0GAAXX_1 test-data/GSM721119_mctp_20F0GAAXX_2 test-data/GSM721118_mctp_20F0GAAXX_3 test-data/GSM721116_mctp_20F0GAAXX_4 test-data/GSM721123_mctp_30CYNAAXX_5 test-data/GSM721124_mctp_209ENAAXX_8 )
ID_ARRAY=( LNCaP_1 LNCaP_2 LNCaP_3 LNCaP_4 PrEC_1 PrEC_2 )
GROUPS_ARRAY=( LNCaP PrEC )
GENDER_ARRAY=( male male male male male male )
