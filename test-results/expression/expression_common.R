#!/usr/bin/Rscript
options(error=utils::recover)
##
## expression_common
##
## Created by Olga Botvinnik on 2012-05-30.
## Copyright (c) 2012 __MyCompanyName__. All rights reserved.
##
## Collaborators: 
## Usage: ./expression_common.R
## Example run: ./expression_common.R

###### Get initial conditions ######
conditionsFile = "/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/diff_exprn_groups.tab"
make.pinkogram.R = "/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/scripts/expression/make.pinkogram.R"

# row.names=1 indicates that the row names are in the first column
# otherwise, no row names are assigned.
bedtools.counts.max = read.delim("/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/bedtools/bedtools_counts_table_max.tab", 
	row.names=1)
bedtools.counts.all = read.delim("/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/bedtools/bedtools_counts_table_all.tab", 
	row.names=1)
bedtools.figures = "/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/bedtools/figures"
bedtools.deseq.dir = "/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/bedtools/figures/DESeq"
bedtools.deseq.prefix = "/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/bedtools/figures/DESeq/"

htseq.counts.max = read.delim("/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/htseq/htseq_counts_table_max.tab", row.names=1)
htseq.counts.all = read.delim("/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/htseq/htseq_counts_table_all.tab", row.names=1)
htseq.figures = "/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/htseq/figures"
htseq.deseq.dir = "/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/htseq/figures/DESeq"
htseq.deseq.prefix = "/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/test-results/expression/htseq/figures/DESeq/"
####################################
# print(ls.str())

# --- Make list type to access name and dataset for each of --- #
# --- max and all types                                     --- #
bedtools.counts = list(list(name="all", ds=bedtools.counts.all), 
	list(name="max", ds=bedtools.counts.max))
htseq.counts = list(list(name="all", ds=htseq.counts.all), 
	list(name="max", ds=htseq.counts.max))

conds = factor(unlist(read.delim(conditionsFile, header=FALSE))) #t(as.matrix(read.delim(conditionsFile, header=FALSE)))
indGroup = grep("group", as.character(conds), ignore.case=TRUE, value=FALSE)
indNotGroup = grep("group", as.character(conds), ignore.case=TRUE, value=FALSE, 
	invert = TRUE)

# print(ls.str())

countTypes = list(list(name="bedtools", counts=bedtools.counts,
	figures=bedtools.figures, deseq.dir=bedtools.deseq.dir,
	deseq.prefix=bedtools.deseq.prefix),
	list(name="htseq", counts=htseq.counts,
	figures=htseq.figures, deseq.dir=htseq.deseq.dir,
	deseq.prefix=htseq.deseq.prefix))

# Red-Blue "pinkogram" for heatmap colors
source(make.pinkogram.R)
pinkogram = make.pinkogram()

# Labels for treatment groups
library(RColorBrewer)
treatmentColorsAllAvailable = brewer.pal(8, "Dark2")
# > library(RColorBrewer)
# > brewer.pal(8, "Dark2")
# [1] "#1B9E77" "#D95F02" "#7570B3" "#E7298A" "#66A61E" "#E6AB02" "#A6761D"
# [8] "#666666"

treatmentColorsSelected = c("#D95F02", "#7570B3")

treatmentColorsPerIndividual = c("#D95F02", "#D95F02", "#D95F02", "#D95F02", "#7570B3", "#7570B3")
treatmentColorsPerGroup = c("#D95F02", "#D95F02", "#7570B3", "#7570B3")

individualTreatmentNames = c("LNCaP_1", "LNCaP_2", "LNCaP_3", "LNCaP_4", "PrEC_1", "PrEC_2")
groupTreatmentNames = c("LNCaP_group1of2", "LNCaP_group1of2", "LNCaP_group2of2", "LNCaP_group2of2")

# --- DESeq-specific --- #

# Get the plotting functions
source("/Users/olgabotvinnik/workspace/rna-seq-diff-exprn/scripts/expression/DESeq/DESeq_utils.R")
