#/usr/bin/Rscript
# This script plots the counts per gene (normalized across samples) on the x-axis versus the empirical dispersion values per gene on the y-axis, found by estimateDispersions in DESeq.
# Estimate the dispersions 'blindly,' i.e. without knowledge
# of the experimental design, such as which samples are
# treated and which are untreated.
# Necessary to create the variance stabilized data file

# Load the DESeq library
library(DESeq)

# Load the data object which contains the dispersion and counts data we need to plot
load('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/htseq/figures/DESeq/htseq_all_individuals_cutoff-0_dispersionEstimates_blind.Robject')

# The expression_common.R file sources the DESeq_utils.R file that has the plotting functions we need
source('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/expression_common.R')

# Open a pdf file to put the graphics we're about to create
pdf('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/htseq/figures/DESeq/htseq_all_individuals_cutoff-0_dispersionEstimates_blind.pdf')

# Plot the counts versus the dispersion estimates
plotDispEsts(cds.blind)

# Turn off the graphics device so the pdf is open-able.
# If you do not do this, your pdf reader will claim this file is corrupt and will not open it.
dev.off()
