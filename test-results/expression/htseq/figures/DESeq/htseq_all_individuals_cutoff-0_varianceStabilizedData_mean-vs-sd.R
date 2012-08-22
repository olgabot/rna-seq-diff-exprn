#/usr/bin/Rscript
# This script plots the per-gene rank of the mean counts (normalized across samples) on the x-axis versus the per-gene standard deviation of the transformed data on the y-axis, found by getVarianceStabilizedData in DESeq.

# Load library containing "meanSdPlot" function
library("vsn")

# Load the data object which contains the standard deviation and counts data we need to plot
load('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/htseq/figures/DESeq/htseq_all_individuals_cutoff-0_varianceStabilizedData_mean-vs-sd.Robject')

# The expression_common.R file sources the DESeq_utils.R file that has the plotting functions we need
source('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/expression_common.R')

# Open a pdf file to put the graphics we're about to create
pdf('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/htseq/figures/DESeq/htseq_all_individuals_cutoff-0_dispersionEstimates.pdf')

# Plot the counts versus the dispersion estimates
meanSdPlot(vsd)

# Turn off the graphics device so the pdf is open-able.
# If you do not do this, your pdf reader will claim this file is corrupt and will not open it.
dev.off()
