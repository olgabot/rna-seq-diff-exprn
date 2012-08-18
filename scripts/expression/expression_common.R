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
conditionsFile = "CONDITIONS_FILE"
make.pinkogram.R = "MAKE_PINKOGRAM_R"

# row.names=1 indicates that the row names are in the first column
# otherwise, no row names are assigned.
bedtools.counts.max = read.delim("BEDTOOLS_COUNTS_MAX_FILE", 
	row.names=1)
bedtools.counts.all = read.delim("BEDTOOLS_COUNTS_ALL_FILE", 
	row.names=1)
bedtools.figures = "BEDTOOLS_FIGS"
bedtools.deseq.dir = "BEDTOOLS_DESEQ_DIR"
bedtools.deseq.prefix = "BEDTOOLS_DESEQ_PREFIX"

htseq.counts.max = read.delim("HTSEQ_COUNTS_MAX_FILE", row.names=1)
htseq.counts.all = read.delim("HTSEQ_COUNTS_ALL_FILE", row.names=1)
htseq.figures = "HTSEQ_FIGS"
htseq.deseq.dir = "HTSEQ_DESEQ_DIR"
htseq.deseq.prefix = "HTSEQ_DESEQ_PREFIX"
####################################
# print(ls.str())

# --- Make list type to access name and dataset for each of --- #
# --- max and all types                                     --- #
bedtools.counts = list(list(name="all", ds=bedtools.counts.all), 
	list(name="max",ds=bedtools.counts.max))
htseq.counts = list(list(name="all", ds=htseq.counts.all), 
	list(name="max",ds=htseq.counts.max))

conds = factor(unlist(read.delim(conditionsFile, header=FALSE))) #t(as.matrix(read.delim(conditionsFile, header=FALSE)))
indGroup = grep("group", conds, ignore.case=TRUE, value=FALSE)
indNotGroup = grep("group", conds, ignore.case=TRUE, value=FALSE, 
	invert = TRUE)

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

treatmentColorsSelected = c(TREATMENT_COLORS_SELECTED)

treatmentColorsPerIndividual = c(INDIVIDUAL_TREATMENT_COLORS)
treatmentColorsPerGroup = c("GROUP_TREATMENT_COLORS")

individualTreatmentNames = c(INDIVIDUAL_TREATMENT_NAMES)
groupTreatmentNames = c("GROUP_TREATMENT_NAMES")

# --- DESeq-specific --- #

# Get the plotting functions
source("DESEQ_SCRIPTS_DIR/DESeq_utils.R")
