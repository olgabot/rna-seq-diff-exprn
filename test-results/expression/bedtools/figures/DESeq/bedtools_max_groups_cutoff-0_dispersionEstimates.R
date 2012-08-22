#/usr/bin/Rscript
# This script plots the counts per gene (normalized across samples) on the x-axis versus the empirical dispersion values per gene on the y-axis, found by estimateDispersions in DESeq.

# Load the DESeq library
library(DESeq)

# Load the data object which contains the dispersion and counts data we need to plot
load('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/bedtools/figures/DESeq/bedtools_max_groups_cutoff-0_dispersionEstimates.Robject')

# The expression_common.R file sources the DESeq_utils.R file that has the plotting functions we need
source('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/expression_common.R')

# Open a pdf file to put the graphics we're about to create
pdf('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/bedtools/figures/DESeq/bedtools_max_groups_cutoff-0_dispersionEstimates.pdf')

# Plot the counts versus the dispersion estimates
plotDispEsts(cds)

# Turn off the graphics device so the pdf is open-able.
# If you do not do this, your pdf reader will claim this file is corrupt and will not open it.
dev.off()
