#!/usr/bin/Rscript

options(error=utils::recover)
library(gplots)  # Plotting with heatmap.2
library(DESeq)   # Differential results/expression
library(vsn)     # Variance stabilization
# source("ssGSEA_analysis.R")
# For heatmap.olga (no clustering) and mutual.inf.2


write_deseq = function(filename){
	write("\n# Load the DESeq library", file=filename, append=TRUE)
	write("library(DESeq)", file=filename, append=TRUE)
}
write_expression_common = function(filename, expression.common.R){
	write("\n# The expression_common.R file sources the DESeq_utils.R file that has the plotting functions we need", file=filename, append=TRUE)
	write(paste(c("source('", expression.common.R, "')"),
		collapse=""), file=filename, append=TRUE)
}

##
## DESeqAnalysis
##
## Created by Olga Botvinnik on 2012-05-18.
## Copyright (c) 2012 __MyCompanyName__. All rights reserved.
##
## Collaborators: 
## Usage: ./DESeqAnalysis.R <arg1> <arg2> <arg3>
## Example run: 
# ./DESeqAnalysis.R counts_per-symbol-max.tab counts_per-symbol-all.tab counts_per-symbol-max.conditions figures/counts_per-symbol 10

## To run on Grettir which doesn't like my /usr/bin/Rscript and won't
## comply with #!/share/apps/bin/R:
## R --slave --args counts_per-symbol-max.tab counts_per-symbol-all.tab counts_per-symbol-max.conditions figures/counts_per-symbol 10 < DESeqAnalysis.R

###### Get initial conditions ######
args = commandArgs(TRUE)

# Gene-level counts, took the max at each gene
# countsMaxFile = args[1]

# Transcript-level counts, have gene symbol but
# looks like: TP53, TP53.1, TP53.2, ... for all the 
# possible transcripts
# countsAllFile = args[2]

# conditionsFile = args[3]

expression.common.R = args[1]

source(expression.common.R)

if( length(args) == 2){
min.cutoff = as.numeric(args[2])   	# Only use genes that have values greater than this
} else{ 				# minimum cutoff for Differential Expression analysis.
	min.cutoff = 0		# If not provided, ignore all genes which are zero
}						# for all samples.
						

# Arguments didn't work, so I'm hard-coding the variables in and testing them
# countsMaxFile = "~/workspace/UCSC/pourmand/2012-spring/results/expression/counts_per-symbol-max.tab"
# countsAllFile = "~/workspace/UCSC/pourmand/2012-spring/results/expression/counts_per-symbol-all.tab"
# conditionsFile = "~/workspace/UCSC/pourmand/2012-spring/results/expression/counts_per-symbol-max.conditions"
# filePrefix = "~/workspace/UCSC/pourmand/2012-spring/results/expression/figures/counts_per-symbol"
# min.cutoff = 10
####################################

# conds = factor(unlist(read.delim(conditionsFile, header=FALSE,
# 	comment.char="#"))) #t(as.matrix(read.delim(conditionsFile, header=FALSE)))
# indGroup = grep("group", conds, ignore.case=TRUE, value=FALSE)
# indNotGroup = grep("group", conds, ignore.case=TRUE, value=FALSE, 
# 	invert = TRUE)

print(ls.str())

# counts.all = list(list(name="ssGSEA", ds=ssGSEA.no.desc))


# print(plotVSD)

# Take out the groups and compare only individual samples
for( counts in countTypes ){
	countTypeName = counts$name
	print(paste("count type name:", countTypeName))
	figures = counts$figures
	deseq.dir =counts$deseq.dir
	deseq.prefix = counts$deseq.prefix

	for( countData in counts$counts){
		# Print out the data type
		print(paste("data type:", countData$name))
		
		counts.ds = countData$ds
		countsName = countData$name
		
		# print("conditions:")
		# print(conds)
		
		# Iterate over comparing 5 or 6 single-cell samples vs comparing 2 groups
		# of 2 or 3 samples in each data type. 
		# e.g. Untreated1 = samples 3, 4, 6 & Untreated2 = samples 7, 8, 9
		# compared vs HighDose1 = 34, 35, 36 & HighDose2 = samples 38, 39, 40
		for( groupType in c("individuals", "groups")){
			
			# Print out the group type
			print(paste("group type:", groupType))

			# Create string of this counts data type and 
			# group for easy access
			min.cutoff.str = paste("cutoff", min.cutoff, sep="-")
			filename.base = paste(deseq.prefix, countTypeName, 
				"_", countsName, "_", groupType, "_", 
				min.cutoff.str, sep="")
			print(paste("filename.base:", filename.base))

			if( groupType == "groups"){
				print("indGroup:")
				print(indGroup)
				if( length(indGroup) == 0){
					print("There are no grouped treatment types to evaluate")
					next
				}
				theseConds = conds[indGroup]
				theseCounts = counts.ds[,indGroup]
				# cds = newCountDataSet(, theseConds)
			} else{
				print("indNotGroup:")
				print(indNotGroup)
				theseConds = conds[indNotGroup]
				theseCounts = counts.ds[,indNotGroup]
				# cds = newCountDataSet(counts.ds.over.min.cutoff[,indNotGroup], theseConds)
			}
			print(head(theseCounts))
			print(theseConds)

			# For debugging, "TTLL11" is a gene on chromosome 9,
			# so this will always give some counts for the example
			# data.
			print("TTLL11 data:")
			print(theseCounts[grep("TTLL11", rownames(theseCounts), value=FALSE),])

			print(paste("number of features before cutoff of", 
				min.cutoff, ":", dim(theseCounts)[1]))
			under.min.cutoff = apply(theseCounts, 1, 
				# Want all features that have at least one sample 
				# with a value above min.cutoff
				function(x){ sum(x > min.cutoff) == 0 } )
			theseCounts = theseCounts[!under.min.cutoff,]
			print(paste("number of features after cutoff of", 
				min.cutoff, ":", dim(theseCounts)[1]))

			# Write the counts that passed the cutoff to a file.
			counts.filename = paste(filename.base, "counts.txt", sep="_")
			write(paste(c("Name", colnames(theseCounts)), collapse="\t"),
				file=counts.filename)
			write.table(theseCounts, counts.filename, 
				quote=FALSE, row.names=TRUE, col.names=FALSE, sep="\t",
				append=TRUE)
			print(paste("wrote counts to", counts.filename))
			
			cds = newCountDataSet(theseCounts, theseConds)
			
			# print("theseConds")
			# print(theseConds)
			possibleComparisons = unique(theseConds)

			# Estimate effective library size
			cds = estimateSizeFactors(cds)
			
			# Estimate gene dispersions
			# Simon Anders' description of why "pooled" from 
			# SEQAnswers forum:
			# "Sometimes, both condition have unequal variance (for 
			# example, knock-down samples might differ strongly from each 
			# other than untreated control samples because knock-down 
			# efficiency is so hard to keep constant), and then, 
			# "per-condition" can give more power. This is why this was 
			# the default. However, I realized recently that our way of 
			# avoiding outliers (see the discussion of 
			# 'sharingMode="maximum"' in the vignette) does not work as 
			# reliably as I hoped when using "per-condition" estimation. 
			# This is why I changed the default to "pooled" and added a 
			# note about this fact to the help page. I have some ideas on 
			# how to improve this matter but pending that I recommend 
			# "pooled"."
			cds = estimateDispersions(cds, fitType="local")
			pdf(paste(filename.base, "_dispersionEstimates.pdf", sep=""))
			save(file=paste(filename.base, 
				"_dispersionEstimates.Robject", sep=""), 
				list="cds")
			plotDispEsts(cds); dev.off()
			cdsFilenameR = paste(filename.base, 
				"_dispersionEstimates.R", sep="")
			write("#/usr/bin/Rscript", file=cdsFilenameR)
			write("# This script plots the counts per gene (normalized across samples) on the x-axis versus the empirical dispersion values per gene on the y-axis, found by estimateDispersions in DESeq.",
				file=cdsFilenameR, append=TRUE)
			write_deseq(cdsFilenameR)
			write("\n# Load the data object which contains the dispersion and counts data we need to plot", file=cdsFilenameR, append=TRUE)
			write(paste(c("load('", paste(filename.base, 
				"_dispersionEstimates.Robject", sep=""), "')"), 
				collapse=""), file=cdsFilenameR, append=TRUE)
			write_expression_common(cdsFilenameR, expression.common.R)
			write("\n# Open a pdf file to put the graphics we're about to create",
						file=cdsFilenameR, append=TRUE)
			write(paste(c("pdf('", paste(filename.base, "_dispersionEstimates.pdf", sep=""), "')"), collapse=""),
				file=cdsFilenameR, append=TRUE)
			write("\n# Plot the counts versus the dispersion estimates",
				file=cdsFilenameR, append=TRUE)
			write("plotDispEsts(cds)", file=cdsFilenameR, append=TRUE)
			write("\n# Turn off the graphics device so the pdf is open-able.\n# If you do not do this, your pdf reader will claim this file is corrupt and will not open it.",
				file=cdsFilenameR, append=TRUE)
			write("dev.off()", file=cdsFilenameR, append=TRUE)

			# Estimate the dispersions "blindly," i.e. without knowledge
			# of the experimental design, such as which samples are 
			# treated and which are untreated. 
			# Necessary to create the variance stabilized data file
			cds.blind = estimateDispersions(cds, method="blind", fitType="local")
			pdf(paste(filename.base, "_dispersionEstimates_blind.pdf", sep=""))
			save(file=paste(filename.base, 
				"_dispersionEstimates_blind.Robject", sep=""), 
				list="cds.blind")
			plotDispEsts(cds.blind); dev.off()
			cds.blindFilenameR = paste(filename.base, 
				"_dispersionEstimates_blind.R", sep="")
			write("#/usr/bin/Rscript", file=cds.blindFilenameR)
			write("# This script plots the counts per gene (normalized across samples) on the x-axis versus the empirical dispersion values per gene on the y-axis, found by estimateDispersions in DESeq.",
				file=cds.blindFilenameR, append=TRUE)
			write("# Estimate the dispersions 'blindly,' i.e. without knowledge\n# of the experimental design, such as which samples are\n# treated and which are untreated.\n# Necessary to create the variance stabilized data file",
				file=cds.blindFilenameR, append=TRUE)
			write_deseq(cds.blindFilenameR)
			write("\n# Load the data object which contains the dispersion and counts data we need to plot", file=cds.blindFilenameR, append=TRUE)
			write(paste(c("load('", paste(filename.base, 
				"_dispersionEstimates_blind.Robject", sep=""), "')"), 
				collapse=""), file=cds.blindFilenameR, append=TRUE)
			write_expression_common(cds.blindFilenameR, 
				expression.common.R)
			write("\n# Open a pdf file to put the graphics we're about to create",
						file=cds.blindFilenameR, append=TRUE)
			write(paste(c("pdf('", paste(filename.base, "_dispersionEstimates_blind.pdf", sep=""), "')"), collapse=""),
				file=cds.blindFilenameR, append=TRUE)
			write("\n# Plot the counts versus the dispersion estimates",
				file=cds.blindFilenameR, append=TRUE)
			write("plotDispEsts(cds.blind)", file=cds.blindFilenameR, append=TRUE)
			write("\n# Turn off the graphics device so the pdf is open-able.\n# If you do not do this, your pdf reader will claim this file is corrupt and will not open it.",
				file=cds.blindFilenameR, append=TRUE)
			write("dev.off()", file=cds.blindFilenameR, append=TRUE)

			
			vsd = getVarianceStabilizedData(cds.blind)
			pdf(paste(filename.base, "_varianceStabilizedData_mean-vs-sd.pdf", 
			sep=""))
			meanSdPlot(vsd)
			dev.off()
			save(list="vsd", file=paste(filename.base, "_varianceStabilizedData_mean-vs-sd.Robject", 
				sep=""))
			meanSdPlotR = paste(filename.base, "_varianceStabilizedData_mean-vs-sd.R", sep="")
			write("#/usr/bin/Rscript", file=meanSdPlotR)
			write("# This script plots the per-gene rank of the mean counts (normalized across samples) on the x-axis versus the per-gene standard deviation of the transformed data on the y-axis, found by getVarianceStabilizedData in DESeq.",
				file=meanSdPlotR, append=TRUE)
			write('\n# Load library containing "meanSdPlot" function', 
				file=meanSdPlotR, append=TRUE)
			write('library("vsn")', file=meanSdPlotR, append=TRUE)
			write("\n# Load the data object which contains the standard deviation and counts data we need to plot", 
				file=meanSdPlotR, append=TRUE)
			write(paste(c("load('", paste(filename.base, "_varianceStabilizedData_mean-vs-sd.Robject", 
				sep=""), "')"), 
				collapse=""), file=meanSdPlotR, append=TRUE)
			write_expression_common(meanSdPlotR, expression.common.R)
			write("\n# Open a pdf file to put the graphics we're about to create",
						file=meanSdPlotR, append=TRUE)
			write(paste(c("pdf('", paste(filename.base, "_dispersionEstimates.pdf", sep=""), "')"), collapse=""),
				file=meanSdPlotR, append=TRUE)
			write("\n# Plot the counts versus the dispersion estimates",
				file=meanSdPlotR, append=TRUE)
			write("meanSdPlot(vsd)", file=meanSdPlotR, append=TRUE)
			write("\n# Turn off the graphics device so the pdf is open-able.\n# If you do not do this, your pdf reader will claim this file is corrupt and will not open it.",
				file=meanSdPlotR, append=TRUE)
			write("dev.off()", file=meanSdPlotR, append=TRUE)
			
			# includes pdf and dev.off functions
			plotVSD(cds, vsd, filename.base)
			# Write the code to re-plot the plotVSD function, if necessary
			vsdR = paste(filename.base, "_varianceStabilizedData.R", sep="")
			write("#/usr/bin/Rscript", file=vsdR)
			write(paste("#", vsdR, collapse=""), file=vsdR, append=TRUE)
			write("# This script plots, per-sample, the per-gene counts  on the x-axis versus the per-gene size factor in the transformed data on the y-axis, found by getVarianceStabilizedData in DESeq.\n# This transformation is used because the high variability between genes can be difficult to visualize and cluster for many algorithms, so we moderate the variability.",
				file=vsdR, append=TRUE)
			write("\n# Load the data object which contains the standard deviation and counts data we need to plot", 
				file=vsdR, append=TRUE)
			write(paste(c("load('", paste(filename.base, "_varianceStabilizedData_mean-vs-sd.Robject", 
				sep=""), "')"), 
				collapse=""), file=vsdR, append=TRUE)
			write(paste(c("load('", paste(filename.base, 
				"_dispersionEstimates.Robject", sep=""), "')"), 
				collapse=""), file=vsdR, append=TRUE)
			# write("print(ls.str())", file=vsdR, append=TRUE)
			write(paste(c("filename.base = '", filename.base, "'"), 
				collapse=""), file=vsdR, append=TRUE)
			write_expression_common(vsdR, expression.common.R)
			write("\n# Plot the counts versus the dispersion estimates\n# includes pdf and dev.off functions",
				file=vsdR, append=TRUE)
			write("plotVSD(cds, vsd, filename.base)", file=vsdR, append=TRUE)

			compared = NULL
			for( groupOne in possibleComparisons){
				for( groupTwo in possibleComparisons){
					# don't compare a group to itself
					if( groupOne == groupTwo) next

					# don't repeat comparisons
					bothGroupsStr.1 = paste(groupOne, groupTwo, sep="-")
					bothGroupsStr.2 = paste(groupTwo, groupOne, sep="-")
					bothGroupsStr = bothGroupsStr.1
					if( length(compared) > 0){
						print(compared)
						print(bothGroupsStr)
						if(sum(bothGroupsStr == compared) > 0) next
					}
					
					# Print out the comparison
					print(paste("comparing", groupOne, "vs", groupTwo))	
					
					compared = c(compared, bothGroupsStr.1, bothGroupsStr.2)
					res = nbinomTest(cds, groupOne, groupTwo)

					# Plot the histogram of unadjusted p-values
					pdf(paste(filename.base, bothGroupsStr,"pval.pdf", 
						sep="_"))
					pValHist(res$pval, paste(groupOne, "vs", groupTwo, "\nunadjusted p-values")); 
					dev.off()
					# write the p-value histogram plot instructions
					# to an R file
					pValHistR = paste(filename.base, bothGroupsStr,
						"pval.R", sep="_")
					write("#/usr/bin/Rscript", file=pValHistR)
					write("# This script plots a histogram showing the distribution of p-values.\n# Ideally, a p-value histogram would be flat (uniform distribution), with equal probability of every value.",
						file=pValHistR, append=TRUE)
					write("\n# Vector of p-values", 
						file=pValHistR, append=TRUE)
					write(paste(c("pval = c(", paste(res$pval,
						collapse=","), ")"), 
						collapse=""), file=pValHistR, append=TRUE)
					write_expression_common(pValHistR, 
						expression.common.R)
					write("\n# Open a pdf file to put the graphics we're about to create",
								file=pValHistR, append=TRUE)
					write(paste(c("pdf('", paste(filename.base, bothGroupsStr, "pval.pdf", sep="_"), "')"), collapse=""),
						file=pValHistR, append=TRUE)
					write("\n# Plot a histogram of the unadjusted p-values",
						file=pValHistR, append=TRUE)
					write(paste(c('pValHist(pval, "', paste(groupOne, "vs", groupTwo, "\nunadjusted p-values"), '")'), 
						collapse=""), file=pValHistR, append=TRUE)
					write("\n# Turn off the graphics device so the pdf is open-able.\n# If you do not do this, your pdf reader will claim this file is corrupt and will not open it.",
						file=pValHistR, append=TRUE)
					write("dev.off()", file=pValHistR, append=TRUE)

					# Plot the histogram of adjusted p-values
					# (adjusted by the Benjamini-Hochberg correction
					# for False Discovery Rates)
					pdf(paste(filename.base, bothGroupsStr,"padj.pdf", sep="_"))
					pValHist(res$padj, paste(groupOne, "vs", groupTwo, "\nadjusted p-values")); 
					dev.off()
					# Write the padjusted histogram plot to an R file
					pAdjHistR = paste(filename.base, bothGroupsStr,
						"padj.R", sep="_")
					write("#/usr/bin/Rscript", file=pAdjHistR)
					write("# This script plots a histogram showing the distribution of p-values adjusted by the Benjamini-Hochberg correction for controlling False Discovery Rate (FDR). \n#See http://en.wikipedia.org/wiki/False_discovery_rate for more information.\n#Ideally, a p-value histogram would be flat (uniform distribution), with equal probability of every value.",
						file=pAdjHistR, append=TRUE)
					write("\n# Vector of adjusted p-values", 
						file=pAdjHistR, append=TRUE)
					write(paste(c("padj = c(", paste(res$padj,
						collapse=","), ")"), 
						collapse=""), file=pAdjHistR, append=TRUE)
					write_expression_common(pAdjHistR, expression.common.R)
					write("\n# Open a pdf file to put the graphics we're about to create",
								file=pAdjHistR, append=TRUE)
					write(paste(c("pdf('", paste(filename.base, bothGroupsStr, "padj.pdf", sep="_"), "')"), collapse=""),
						file=pAdjHistR, append=TRUE)
					write("\n# Plot adjusted p-values on a histogram",
						file=pAdjHistR, append=TRUE)
					write(paste(c('pValHist(padj, "', paste(groupOne, "vs", groupTwo, "\nadjusted p-values"), '")'), 
						collapse=""), file=pAdjHistR, append=TRUE)
					write("\n# Turn off the graphics device so the pdf is open-able.\n# If you do not do this, your pdf reader will claim this file is corrupt and will not open it.",
						file=pAdjHistR, append=TRUE)
					write("dev.off()", file=pAdjHistR, append=TRUE)
					
					# Plot the log2FoldChange
					pdf(paste(filename.base, bothGroupsStr,
						"log2FoldChange.pdf", sep="_"))
					plotDE(res) ; dev.off()
					# Write the log2FoldChange plot instructions to a file
					log2FoldChangeR = paste(filename.base, bothGroupsStr,
						"log2FoldChange.R", sep="_")
					write("#/usr/bin/Rscript", file=log2FoldChangeR)
					write("# This script plots the normalized mean (x-axis) versus log2 fold change (y-axis). This plot is also called the 'MA-plot.'",
						file=log2FoldChangeR, append=TRUE)
					write("\n# Load the data object which contains the results from the n-binomial test we used to test for differential expression.\n# This includes the log2FoldChange, baseMean, and adjusted and unadjusted p-values.", 
						file=log2FoldChangeR, append=TRUE)
					write(paste(c("res = read.delim('", 
						paste(filename.base, bothGroupsStr, 
							"DEseq.txt", sep="_"), "')"), 
						collapse=""), file=log2FoldChangeR, append=TRUE)
					write_expression_common(log2FoldChangeR, expression.common.R)
					write("\n# Open a pdf file to put the graphics we're about to create",
								file=log2FoldChangeR, append=TRUE)
					write(paste(c("pdf('", paste(filename.base, bothGroupsStr, "log2FoldChange.pdf", sep="_"), "')"), collapse=""),
						file=log2FoldChangeR, append=TRUE)
					write("\n# Plot of normalized mean versus log2FoldChange",
						file=log2FoldChangeR, append=TRUE)
					write("plotDE(res)", file=log2FoldChangeR, append=TRUE)
					write("\n# Turn off the graphics device so the pdf is open-able.\n# If you do not do this, your pdf reader will claim this file is corrupt and will not open it.",
						file=log2FoldChangeR, append=TRUE)
					write("dev.off()", file=log2FoldChangeR, append=TRUE)

					
					# Moderated log fold change estimates via variance 
					# stabilizing transformation
					mod_lfc = (rowMeans(vsd[,conditions(cds)==groupOne, 
							drop=FALSE]) -
						rowMeans(vsd[,conditions(cds)==groupTwo, drop=FALSE]))
					lfc = res$log2FoldChange
					finite = is.finite(lfc)
					# print(table(as.character(lfc[!finite]), 
					# 	useNA="always"))
					
					largeNumber = 7
					lfc = ifelse(finite, lfc, sign(lfc)*largeNumber)
					log.decade = 1 + round(log10(1+rowMeans(counts(cds.blind, 
						normalized=TRUE))))
					colors = colorRampPalette(	c("gray", "blue"))(6)[log.decade]
					pdf(paste(filename.base, bothGroupsStr,
						"log2FoldChange-vs-ModeratedLogRatio.pdf", sep="_"))
					# Plot the moderated log fold change
					plot(lfc, mod_lfc, pch=20, cex=.4, asp=1, 
						col=ifelse(finite,colors,"purple"),
						main="Direct (lfc) vs Moderated log ratios (mod_lfc)",
						sub="purple were infinite in lfc, arbitrarily set for plotting")
					abline(a=0, b=1, col="#40C04040")
					dev.off()
					# Save the instructions for plotting moderated log
					# fold change in an R script
					mod_lfcR = paste(filename.base, bothGroupsStr,
						"log2FoldChange-vs-ModeratedLogRatio.R", sep="_")
					write("#/usr/bin/Rscript", file=mod_lfcR)
					write("# This script plots the direct log fold change (lfc, x-axis) versus the moderated log-ratios (mod_lfc, y-axis).\n# [From DESeq documentation] 'The moderation criterion used is variance stabilization.\n# The points are colored in a scale from gray to blue, representing weakly to strongly expressed genes.\n# The purple points correspond to values that were infinite in lfc and were arbitrarily set to +/- 7 for the purpose of plotting.'",
						file=mod_lfcR, append=TRUE)
					write("\n# Vector of log2 fold change values", 
						file=mod_lfcR, append=TRUE)
					write(paste(c("lfc = c(", paste(lfc, collapse=","), ")"), 
						collapse=""), file=mod_lfcR, append=TRUE)
					write("\n# Vector of log2 fold change values, moderated by variance stabilization", 
						file=mod_lfcR, append=TRUE)
					write(paste(c("mod_lfc = c(", paste(mod_lfc, collapse=","), ")"), 
						collapse=""), file=mod_lfcR, append=TRUE)
					write("\n# Boolean vector indicating whether lfc is finite or not", 
						file=mod_lfcR, append=TRUE)
					write(paste(c("finite = c(", paste(finite, collapse=","), ")"), 
						collapse=""), file=mod_lfcR, append=TRUE)
					write("\n# Vector of colors", 
						file=mod_lfcR, append=TRUE)
					write(paste(c("colors = c('", paste(colors, 
						collapse="','"), "')"), 
						collapse=""), file=mod_lfcR, append=TRUE)
					write("\n# Open a pdf file to put the graphics we're about to create",
								file=mod_lfcR, append=TRUE)
					write(paste(c("pdf('", paste(filename.base, bothGroupsStr, "log2FoldChange-vs-ModeratedLogRatio.pdf", sep="_"), "')"), collapse=""),
						file=mod_lfcR, append=TRUE)
					write("\n# Plot of normalized mean versus log2FoldChange",
						file=mod_lfcR, append=TRUE)
					write('plot(lfc, mod_lfc, pch=20, cex=.4, asp=1, 
						col=ifelse(finite,colors,"purple"),
						main="Direct (lfc) vs Moderated log ratios (mod_lfc)",
						sub="purple were infinite in lfc, arbitrarily set for plotting")
					abline(a=0, b=1, col="#40C04040")', file=mod_lfcR, append=TRUE)
					write("\n# Turn off the graphics device so the pdf is open-able.\n# If you do not do this, your pdf reader will claim this file is corrupt and will not open it.",
						file=mod_lfcR, append=TRUE)
					write("dev.off()", file=mod_lfcR, append=TRUE)
					
					olfactoryReceptors = grep("^OR", res$id)
					
					indPval = order(res$pval)[-olfactoryReceptors][1:40]
					heatmapFilename = paste(filename.base, bothGroupsStr, 
						"top40DE-geneClustering", sep="_")
					pdf(paste(heatmapFilename, ".pdf", sep=""))
					heatmap.2(vsd[indPval,], col=pinkogram, 
						main=paste(groupOne, "vs", groupTwo, 
						"\nmin.cutoff =",min.cutoff))
					# heatmap(vsd[indPval,], col=pinkogram, 
					# 	main=paste(groupOne, "vs", groupTwo),
					# 	sub="plotting method: heatmap()")
					dev.off()

					# Save the heatmap plotting instructions to a file
					# for reproducibility
					heatmapFilenameR = paste(heatmapFilename, ".R", 
						sep="")
					write("#/usr/bin/Rscript", file=heatmapFilenameR,
						append=FALSE)
					write(paste(c("# This script plots a heatmap of the top 40 differentially expressed genes (minus the Olfactory Receptors, because they often turn up as highly differentially expressed even in non-olfactory cells) between", 
						groupOne, "and", groupTwo), collapse=" "),
						file=heatmapFilenameR, append=TRUE)
					write("\n# Indices of the top 40 differentially expressed genes with the best adjusted p-value (padj)", file=heatmapFilenameR, append=TRUE)
					write("# Olfactory Receptors (gene symbols starting with OR) have been removed because they are highly homologous with one another and cause problems in identifying top differentially expressed genes.",
						file=heatmapFilenameR, append=TRUE)
					write("# All the genes can be seen in padj order in the file:",
						file=heatmapFilenameR, append=TRUE)
					write(paste("#", paste(filename.base, bothGroupsStr, 
							"DEseq.txt", sep="_")),
						file=heatmapFilenameR, append=TRUE)
					write(paste(c("indPval = c(", paste(indPval, 
						collapse=","), ")"), collapse=""), 
						file=heatmapFilenameR, append=TRUE)
					write("\n# Load gplots library for heatmap.2",
						file=heatmapFilenameR, append=TRUE)
					write("library(gplots)",
						file=heatmapFilenameR, append=TRUE)
					write("\n# Variance-stabilized data used for plotting", file=heatmapFilenameR, append=TRUE)
					write(paste(c("vsd = read.delim('", 
						paste(filename.base, "vsd.txt", sep="_"),
						"', row.names=1)"), collapse=""), file=heatmapFilenameR,
						append=TRUE)
					write("\n# Red-blue 'pinkogram' colors for plotting, red=hot, high expression ; blue=cold, low expression", file=heatmapFilenameR, append=TRUE)
					write(paste(c("pinkogram = c('", paste(pinkogram,
						collapse="','"), "')"), collapse=""), 
						file=heatmapFilenameR,
						append=TRUE)
					write("print(ls.str())",file=heatmapFilenameR, append=TRUE)
					write("\n# Open a pdf file to put the graphics we're about to create",
						file=heatmapFilenameR, append=TRUE)
					write(paste(c("pdf('", 
						paste(heatmapFilename, ".pdf", sep=""), "')"),
						collapse=""), 
						file=heatmapFilenameR, append=TRUE)
					write("\n# Plot the heatmap. \n# To remove the cyan line through the heatmap, set trace='none'\n# If your sample or gene names are longer than the margins, you can increase the margin size from the default margin=c(5,5) to margin=(10,10) for example.", 
						file=heatmapFilenameR, append=TRUE)
					write(paste(c("heatmap.2(as.matrix(vsd[indPval,]), col=pinkogram, main='", paste(groupOne, 
						"vs", groupTwo, 
						"\nmin.cutoff =", min.cutoff), "')"), collapse=""), file=heatmapFilenameR, append=TRUE)
					write("\n# Turn off the graphics device so the pdf is open-able.\n# If you do not do this, your pdf reader will claim this file is corrupt and will not open it.",
						file=heatmapFilenameR, append=TRUE)
					write("dev.off()", file=heatmapFilenameR, append=TRUE)

					dists = dist(t(vsd))
					# Cluster the samples by Euclidean distance of gene
					# expression
					pdf(paste(filename.base, bothGroupsStr, 
						"sampleClustering.pdf", sep="_"))
					heatmap(as.matrix(dists), symm=TRUE, scale="none",
						margins=c(10,10), col=colorRampPalette(c("darkblue",
							"white"))(100),
							labRow=paste(pData(cds.blind)$condition))
					dev.off()
					# Save the script to cluster the samples by Euclidean
					# distance in an R script
					sampleClusteringR = paste(filename.base, 
						bothGroupsStr, 
						"sampleClustering.R", sep="_")
					write("#/usr/bin/Rscript", file=sampleClusteringR,
						append=FALSE)
					write("# This script plots a heatmap of the Euclidean distances between all the samples, and clusters the samples according to similarity.",
						file=sampleClusteringR, append=TRUE)
					write("\n# Variance-stabilized data used for plotting", file=sampleClusteringR, append=TRUE)
					write(paste(c("vsd = read.delim('", 
						paste(filename.base, "vsd.txt", sep="_"),
						"', row.names=1)"), collapse=""), file=sampleClusteringR,
						append=TRUE)
					write("\n# Get the Euclidean distance between samples",
						file=sampleClusteringR, append=TRUE)
					write("dists = dist(t(vsd))",
						file=sampleClusteringR, append=TRUE)
					write("\n# Labels for the rows based on the sample conditions, e.g. their treatment groups",
						file=sampleClusteringR, append=TRUE)
					write(paste(c("labRow = c('", paste(pData(cds.blind)$condition, collapse="','"), "')"), collapse=""),
						file=sampleClusteringR, append=TRUE)
					write("\n# Open a pdf file to put the graphics we're about to create",
						file=sampleClusteringR, append=TRUE)
					write(paste(c("pdf('", paste(filename.base, 
						bothGroupsStr, "sampleClustering.pdf", 
						sep="_"), "')"), collapse=""), 
						file=sampleClusteringR, append=TRUE)
					write("\n# Plot the heatmap, clustering by both rows and columns.", 
						file=sampleClusteringR, append=TRUE)
					write('heatmap(as.matrix(dists), symm=TRUE, scale="none", margins=c(10,10), col=colorRampPalette(c("darkblue","white"))(100),labRow=labRow)', file=sampleClusteringR, append=TRUE)
					write("\n# Turn off the graphics device so the pdf is open-able.\n# If you do not do this, your pdf reader will claim this file is corrupt and will not open it.",
						file=sampleClusteringR, append=TRUE)
					write("dev.off()", file=sampleClusteringR, append=TRUE)
					
					# Write down the genes!
					 write.table(res[order(res$pval),], 
						paste(filename.base, bothGroupsStr, 
							"DEseq.txt", sep="_"), 
						row.names=FALSE, col.names=TRUE, 
						append=FALSE, quote=FALSE, sep="\t")
					
					# indPadj = order(res$padj)[1:40]
					# pdf
					# heatmap.2()
				}
			}
			vsd.filename = paste(filename.base, "vsd.txt", sep="_")
			write(paste(c("Name", colnames(vsd)), collapse="\t"),
				file=vsd.filename)
			write.table(vsd, vsd.filename, 
				quote=FALSE, row.names=TRUE, col.names=FALSE, sep="\t",
				append=TRUE)
			# using vsd for MI produces better results than with cds 
			# (more obvious divide between classes)
			# binary.vector = grep("Survivors", theseConds)
			# vsd.MI = apply(vsd, 1, mutual.inf.2, ifSurvivors)
		}
	}
}