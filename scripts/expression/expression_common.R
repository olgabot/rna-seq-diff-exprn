#!/usr/bin/Rscript
options(error=utils::recover)
##
## expression_common
##
## Created by Olga Botvinnik on 2012-05-30.
## Copyright (c) 2012 __MyCompanyName__. All rights reserved.
##
## Collaborators: 
## Usage: ./expression_common.R <arg1> <arg2> <arg3>
## Example run: ./expression_common.R

###### Get initial conditions ######
conditionsFile = "/home/obot/single-cell/results/expression/conditions.tab"
pinkogram.R = "/home/obot/single-cell/scripts/expression/make.pinkogram.R"
####################################
# print(ls.str())

conds = factor(unlist(read.delim(conditionsFile, header=FALSE))) #t(as.matrix(read.delim(conditionsFile, header=FALSE)))
indGroup = grep("group", conds, ignore.case=TRUE, value=FALSE)
indNotGroup = grep("group", conds, ignore.case=TRUE, value=FALSE, 
	invert = TRUE)


# Red-Blue "pinkogram" for heatmap colors
source(pinkogram.R)
pinkogram = make.pinkogram()

# Labels for treatment groups
library(RColorBrewer)
treatmentColors = brewer.pal(3, "Set2")

individualTreatmentColors = c(rep(treatmentColors[3],6),rep(treatmentColors[2],6),
	rep(treatmentColors[1],5))
groupTreatmentColors = c(rep(treatmentColors[3],2),rep(treatmentColors[2],2),
	rep(treatmentColors[1],2))
individualTreatmentNames = c(paste("untreated", 1:6), paste("treated",1:6),
		paste("survivor", 1:5))
groupTreatmentNames = 	c(paste("untreated group", 1:2), 
		paste("treated group",1:2), paste("survivor group", 1:2))