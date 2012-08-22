#/usr/bin/Rscript
# This script plots a histogram showing the distribution of p-values adjusted by the Benjamini-Hochberg correction for controlling False Discovery Rate (FDR). 
#See http://en.wikipedia.org/wiki/False_discovery_rate for more information.
#Ideally, a p-value histogram would be flat (uniform distribution), with equal probability of every value.

# Vector of adjusted p-values
padj = c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.542953783051901,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.149149980887576,1,0.447788249814043,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.313552285593236,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.62951711741908,0.447788249814043,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.80119072735728,1,0.80119072735728,1,1,1,1,1,0.80119072735728,1,1,1,0.951900203233576,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0.542953783051901,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

# The expression_common.R file sources the DESeq_utils.R file that has the plotting functions we need
source('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/expression_common.R')

# Open a pdf file to put the graphics we're about to create
pdf('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/htseq/figures/DESeq/htseq_all_groups_cutoff-0_LNCaP_group-PrEC_group_padj.pdf')

# Plot adjusted p-values on a histogram
pValHist(padj, "LNCaP_group vs PrEC_group 
adjusted p-values")

# Turn off the graphics device so the pdf is open-able.
# If you do not do this, your pdf reader will claim this file is corrupt and will not open it.
dev.off()
