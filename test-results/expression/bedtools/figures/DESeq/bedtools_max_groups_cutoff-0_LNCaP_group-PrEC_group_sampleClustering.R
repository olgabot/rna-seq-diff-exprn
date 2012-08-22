#/usr/bin/Rscript
# This script plots a heatmap of the Euclidean distances between all the samples, and clusters the samples according to similarity.

# Variance-stabilized data used for plotting
vsd = read.delim('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/bedtools/figures/DESeq/bedtools_max_groups_cutoff-0_vsd.txt', row.names=1)

# Get the Euclidean distance between samples
dists = dist(t(vsd))

# Labels for the rows based on the sample conditions, e.g. their treatment groups
labRow = c('LNCaP_group','LNCaP_group','PrEC_group','PrEC_group')

# Open a pdf file to put the graphics we're about to create
pdf('/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/bedtools/figures/DESeq/bedtools_max_groups_cutoff-0_LNCaP_group-PrEC_group_sampleClustering.pdf')

# Plot the heatmap, clustering by both rows and columns.
heatmap(as.matrix(dists), symm=TRUE, scale="none", margins=c(10,10), col=colorRampPalette(c("darkblue","white"))(100),labRow=labRow

# Turn off the graphics device so the pdf is open-able.
# If you do not do this, your pdf reader will claim this file is corrupt and will not open it.
dev.off()
