#!/usr/bin/Rscript
options(error=utils::recover)
##
## DESeq_utils
##
## Created by Olga Botvinnik on 2012-05-20.
## Copyright (c) 2012 __MyCompanyName__. All rights reserved.
##
## Collaborators: 
## Usage: ./DESeq_utils.R <arg1> <arg2> <arg3>
## Example run: ./DESeq_utils.R

###### Get initial conditions ######
# args = commandArgs(TRUE)
####################################

# Helper plotting functions
pValHist = function(pVals, plot.title){
	hist(pVals, breaks=100, col="skyblue", 
		border="slateblue", 
		main=plot.title)
}

plotDispEsts = function(cds){ 
 plot(rowMeans(counts(cds, normalized=TRUE)), fitInfo(cds)$perGeneDispEsts, 
	pch = ".", log="xy")
 xg = 10^seq(-.5, 5, length.out=300)
 lines(xg, fitInfo(cds)$dispFun(xg), col="red") }

plotDE = function(res){
plot(res$baseMean, res$log2FoldChange, log="x", pch=20,cex=.3, 
	col=ifelse(res$padj < .1, "red", "black"))
}

plotVSD = function(cds, vsd, filename.base){
	# Plot variance stabilizing transformation 
	sampleNamesCDS = sampleNames(phenoData(cds))
	# sampleNamesVSD = sampleNames(phenoData(vsd))
	sampleNames = sampleNamesCDS
	nSamples = dim(vsd)[2]
	
	# print("sampleNamesCDS")
	# print(sampleNamesCDS)
	# print("cds")
	# print(head(counts(cds)))
	# print("vsd")
	# print(head(vsd))
	
	pdf(paste(filename.base, "_varianceStabilizingTransformation.pdf", sep=""))
	for( sampleInd in 1:nSamples){
		px     = counts(cds)[,sampleInd] / sizeFactors(cds)[sampleInd]
		ord    = order(px)
		ord    = ord[px[ord] < 150]
		ord    = ord[seq(1, length(ord), length=50)]
		last   = ord[length(ord)]
		colors = c("blue", "black")
		matplot(px[ord],
		        cbind(vsd[, sampleInd], log2(px))[ord, ],
		        type="l", lty=1, col=colors, xlab="n", ylab="f(n)",
				main=paste(sampleNames[sampleInd]))
		legend("bottomright",
		       legend = c(
		        expression("variance stabilizing transformation"),
		        expression(log[2](n/s))),
		       fill=colors)
	}
	dev.off()
}

make.pinkogram = function(){
	# Makes a high-resolution blue-red (blue=low, red=high, can think of it 
	# as cold and hot) color map to represent values in heatmaps. 
	# An alternative to red-green heatmap colors, which are unreadable for
	# those that are color-blind.
	# Author: Olga Botvinnik
	mycol = vector(length=512, mode = "numeric")
	for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
	for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, 
		maxColorValue=255)
	mycol = rev(mycol)
	return(mycol)
}
pinkogram = make.pinkogram()