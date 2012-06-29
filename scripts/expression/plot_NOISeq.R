#!/usr/bin/Rscript
options(error=utils::recover)
library(gplots) # for heatmap.2
##
## plot_NOISeq
##
## Created by Olga Botvinnik on 2012-06-07.
## Copyright (c) 2012 __MyCompanyName__. All rights reserved.
##
## Collaborators: 
## Usage: ./plot_NOISeq.R <arg1> <arg2> <arg3>
## Example run: ./plot_NOISeq.R $EXPR_DIR/coverageBed/DESeq/counts_per-symbol "cutoff-10_vsd.txt" $EXPR_DIR/coverageBed/NOISeq/NOISeq "topDE.txt"

# on grettir:
# R --slave --args $EXPR_DIR/coverageBed/DESeq/counts_per-symbol "_cutoff-10_vsd.txt" $EXPR_DIR/coverageBed/NOISeq/NOISeq "topDE.txt" < ./plot_NOISeq.R

###### Get initial conditions ######
args = commandArgs(TRUE)

# sample filename:
# counts_per-symbol_max_groups_cutoff-10_vsd.txt

# Before "max" or "uniq" and "groups" or "individuals"
vsdPrefix = args[1]
#"/Users/olgabotvinnik/workspace/UCSC/pourmand/2012-Spring/results/coverageBed/expression/DESeq/counts_per-symbol_"

# After "groups" or "individuals"
vsdSuffix = args[2] #"_cutoff-10_vsd.txt"

# Example file:
# NOISeq_uniq_individuals_uqua_Untreated-Survivors_topDE.txt

###### On Olga's personal computer ######
# noiseqPrefix = "/Users/olgabotvinnik/workspace/UCSC/pourmand/2012-Spring/results/coverageBed/expression/NOISeq/NOISeq_"

###### On grettir ######
noiseqPrefix = args[3]

noiseqSuffix = args[4] #"_topDE.txt"

# Get red-blue "pinkogram" for plotting and other useful utilities
# source("/Users/olgabotvinnik/workspace/UCSC/pourmand/2012-Spring/expression/DESeq_utils.R")
source("~/single-cell/scripts/expression/expression_common.R")
####################################
# print(ls.str())

countTypes = c("max", "uniq")
groupType = c("groups", "individuals")
normalizationTypes = c(
	"rpkm",  # reads per kilobase of exon mapped, default of NOISeq
	"uqua", # Upper-quartile
	"tmm",  # Trimmed mean of m
	"n")    # No normalization
individualComparisons = c("Untreated-HighDose", "Untreated-Survivors",
	"HighDose-Survivors")
groupComparisons = 	c("UntreatedGroup-HighDoseGroup", 
	"UntreatedGroup-SurvivorsGroup",
	"HighDoseGroup-SurvivorsGroup")

for( countType in countTypes ){
	for( groupType in groupType ){
			
			vsdFile = paste(vsdPrefix, countType, groupType, vsdSuffix,
			sep="_")
			vsd = read.delim(vsdFile, header=FALSE,stringsAsFactors=FALSE) 
			
		for( normType in normalizationTypes){
			
			if( groupType == "groups" ){
				comparisons = groupComparisons
				treatmentColors = groupTreatmentColors
				treatmentNames = groupTreatmentNames
			} else{ 
				comparisons = individualComparisons
				treatmentColors = individualTreatmentColors
				treatmentNames = individualTreatmentNames
			}
			
			for( comparison in comparisons ){
					topDEPrefix = paste(noiseqPrefix, countType,
						groupType, normType, comparison, sep="_")
					topDEFile = paste(topDEPrefix, noiseqSuffix, sep="_")
					topDE = read.delim(topDEFile, stringsAsFactors=FALSE)
					
					# Remove olfactory receptors
					olfactoryReceptors = grep("^OR", topDE[,1])
					if (length(olfactoryReceptors) > 0){
						topDE = topDE[-olfactoryReceptors,]
					}
					
					# print("removed olfactory receptors")
					
					# browser()
					
					pdf(paste(topDEPrefix, ".pdf", sep=""), width=8,
						height=17)
					# vsd[,1] %in% topDE[,1] : get the top differentially 
					# expressed genes
					# [...,-1] : remove the first column because it has gene 
					# names and it forces everything to be a character
					# heatm = as.matrix(vsd[vsd[,1] %in% topDE[,1],-1])
					##### For some reason, labeling via rownames and colnames 
					##### is very slow
					# Get row names from first column
					# rownames(heatm) = vsd[vsd[,1] %in% topDE_noOR[,1],1]
					# colnames(heatm) = treatmentNames
					topDEind = (vsd[,1] %in% topDE[,1])
					if( length(topDEind) > 100){
						topDEind = topDEind[1:100]
					}
					print("found top DE indices in vsd")
					heatmap.2(as.matrix(vsd[topDEind,-1]), 
						col=pinkogram,
						main=paste("Top DE genes by NOISeq:", countType, 
							groupType, normType, comparison),
						margins=c(12,15),
						trace="none",
						labRow = vsd[topDEind,1],
						labCol = treatmentNames,
						ColSideColors=treatmentColors)
					print("plotted heatmap")
					dev.off()
					print(paste("Finished:", countType, groupType, normType, 
						comparison))
				}
			}
		}
	}