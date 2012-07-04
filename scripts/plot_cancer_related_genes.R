#!/usr/bin/Rscript
options(error=utils::recover)
##
## plot_cancer_related_genes
##
## Created by Olga Botvinnik on 2012-05-23.
## Copyright (c) 2012 __MyCompanyName__. All rights reserved.
##
## Collaborators: 
## Usage: ./plot_cancer_related_genes.R <arg1> <arg2> <arg3>
## Example run: ./plot_cancer_related_genes.R

###### Get initial conditions ######
args = commandArgs(TRUE)
####################################
# print(ls.str())
library(gplots)
library(DESeq)
source("make.pinkogram.R")
pinkogram = make.pinkogram()

min.cutoff = 10

under.min.cutoff = apply(counts.symbols.max, 1, 
	# Want all features that have at least one sample 
	# with a value above min.cutoff
	function(x){ sum(x > min.cutoff) == 0 } )
counts.max.over10 = counts.symbols.max[!under.min.cutoff,]

under.min.cutoff = apply(counts.symbols.uniq, 1, 
	# Want all features that have at least one sample 
	# with a value above min.cutoff
	function(x){ sum(x > min.cutoff) == 0 } )
counts.uniq.over10 = counts.symbols.uniq[!under.min.cutoff,]


cds.max.individuals = newCountDataSet(counts.max.over10[,1:17], conds[1:17])
cds.max.groups = newCountDataSet(counts.max.over10[,18:23], conds[18:23])

cds.uniq.individuals = newCountDataSet(counts.uniq.over10[,1:17], conds[1:17])
cds.uniq.groups = newCountDataSet(counts.uniq.over10[,18:23], conds[18:23])

cds.max.individuals = estimateSizeFactors(cds.max.individuals)
cds.max.groups = estimateSizeFactors(cds.max.groups)
cds.uniq.individuals = estimateSizeFactors(cds.uniq.individuals)
cds.uniq.groups = estimateSizeFactors(cds.uniq.groups)

cds.max.individuals.blind = estimateDispersions(cds.max.individuals, method="blind", fitType="local")
cds.max.groups.blind = estimateDispersions(cds.max.groups, method="blind", fitType="local")
cds.uniq.individuals.blind = estimateDispersions(cds.uniq.individuals, method="blind", fitType="local")
cds.uniq.groups.blind = estimateDispersions(cds.uniq.groups, method="blind", fitType="local")

vsd.max.individuals = getVarianceStabilizedData(cds.max.individuals.blind)
vsd.max.groups = getVarianceStabilizedData(cds.max.groups.blind)
vsd.uniq.individuals = getVarianceStabilizedData(cds.uniq.individuals.blind)
vsd.uniq.groups = getVarianceStabilizedData(cds.uniq.groups.blind)

genesets = read.delim("fernando_genesets_noparen.gmt", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
cellCycleInd = grep("CellCycle", genesets[,1])
cellCycle = sapply(cellCycleInd, function(x){ genes = genesets[x,which(genesets[x,] != "")][-1:-2] ; which(rownames(vsd.max.individuals) %in% genes) })

individualTreatmentNames = c(paste("untreated", 1:6), paste("treated",1:6),
	paste("survivor", 1:5))
groupTreatmentNames = 	c(paste("untreated group", 1:2), 
	paste("treated group",1:2), paste("survivor group", 1:2))


colnames(vsd.max.individuals) = individualTreatmentNames
colnames(vsd.max.groups) = groupTreatmentNames
colnames(vsd.uniq.individuals) = individualTreatmentNames
colnames(vsd.uniq.groups) = groupTreatmentNames

library(RColorBrewer)
treatmentColors = brewer.pal(3, "Set2")

individualTreatmentColors = c(rep(treatmentColors[3],6),rep(treatmentColors[2],6),
	rep(treatmentColors[1],5))
groupTreatmentColors = c(rep(treatmentColors[3],2),rep(treatmentColors[2],2),
	rep(treatmentColors[1],2))

unplotted.genes = NULL
apply(genesets, 1, function(x){
	filename.base = "figures/max"
	
	geneset = x[1]
	# print(geneset)
	genes = x[which(x != "")][-1:-2]
	indMax = which(rownames(vsd.max.individuals) %in% genes)
	if( length(indMax) > 1){
		max.individ = vsd.max.individuals[indMax,]
		max.groups = vsd.max.groups[indMax,]
	
	pdf.file = paste(filename.base, "_individuals_unclusteredSamples_", geneset, ".pdf", sep="")
	pdf(pdf.file, height=11, width=8.5)
	heatmap.2(max.individ, col=pinkogram, Colv=FALSE,
		main=geneset, margins=c(7,12),
		ColSideColors=individualTreatmentColors)
	dev.off()
	
	pdf.file = paste(filename.base, "_groups_unclusteredSamples_", geneset, ".pdf", sep="")
	pdf(pdf.file, height=11, width=8.5)
	# par(mar=rep(5,4))
	heatmap.2(max.groups, col=pinkogram, Colv=FALSE,
		main=geneset, margins=c(12,15),
		ColSideColors=groupTreatmentColors)
	# heatmap(vsd[indPval,], col=pinkogram, 
	# 	main=paste(groupOne, "vs", groupTwo),
	# 	sub="plotting method: heatmap()")
	dev.off()
	} else{ unplotted.genes = c(unplotted.genes, genes)}
})
#indPval = order(res$pval)[1:40]
