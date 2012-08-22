#/usr/bin/Rscript
# /Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/bedtools/figures/DESeq/bedtools_max_individuals_cutoff-0_varianceStabilizedData.R
# This script plots, per-sample, the per-gene counts  on the x-axis versus the per-gene size factor in the transformed data on the y-axis, found by getVarianceStabilizedData in DESeq.
# This transformation is used because the high variability between genes can be difficult to visualize and cluster for many algorithms, so we moderate the variability.

# Load the data object which contains the standard deviation and counts data we need to plot
load('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/bedtools/figures/DESeq/bedtools_max_individuals_cutoff-0_varianceStabilizedData_mean-vs-sd.Robject')
load('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/bedtools/figures/DESeq/bedtools_max_individuals_cutoff-0_dispersionEstimates.Robject')
filename.base = '/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/bedtools/figures/DESeq/bedtools_max_individuals_cutoff-0'

# The expression_common.R file sources the DESeq_utils.R file that has the plotting functions we need
source('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/expression_common.R')

# Plot the counts versus the dispersion estimates
# includes pdf and dev.off functions
plotVSD(cds, vsd, filename.base)
