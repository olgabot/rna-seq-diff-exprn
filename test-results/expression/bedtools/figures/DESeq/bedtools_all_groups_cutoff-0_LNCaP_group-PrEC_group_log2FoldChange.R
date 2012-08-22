#/usr/bin/Rscript
# This script plots the normalized mean (x-axis) versus log2 fold change (y-axis). This plot is also called the 'MA-plot.'

# Load the data object which contains the results from the n-binomial test we used to test for differential expression.
# This includes the log2FoldChange, baseMean, and adjusted and unadjusted p-values.
res = read.delim('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/bedtools/figures/DESeq/bedtools_all_groups_cutoff-0_LNCaP_group-PrEC_group_DEseq.txt')

# The expression_common.R file sources the DESeq_utils.R file that has the plotting functions we need
source('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/expression_common.R')

# Open a pdf file to put the graphics we're about to create
pdf('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/bedtools/figures/DESeq/bedtools_all_groups_cutoff-0_LNCaP_group-PrEC_group_log2FoldChange.pdf')

# Plot of normalized mean versus log2FoldChange
plotDE(res)

# Turn off the graphics device so the pdf is open-able.
# If you do not do this, your pdf reader will claim this file is corrupt and will not open it.
dev.off()
