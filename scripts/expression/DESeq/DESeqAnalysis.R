#!/usr/bin/Rscript

options(error=utils::recover)
library(gplots)  # Plotting with heatmap.2
library(DESeq)   # Differential results/expression
library(vsn)     # Variance stabilization
source("ssGSEA_analysis.R")# For heatmap.olga (no clustering) and mutual.inf.2
##
## DESeqAnalysis
##
## Created by Olga Botvinnik on 2012-05-18.
## Copyright (c) 2012 __MyCompanyName__. All rights reserved.
##
## Collaborators: 
## Usage: ./DESeqAnalysis.R <arg1> <arg2> <arg3>
## Example run: ./DESeqAnalysis.R counts_per-symbol-max.tab counts_per-symbol-uniq.tab counts_per-symbol-max.conditions figures/counts_per-symbol 10

## To run on Grettir which doesn't like my /usr/bin/Rscript and won't
## comply with #!/share/apps/bin/R:
## R --slave --args counts_per-symbol-max.tab counts_per-symbol-uniq.tab counts_per-symbol-max.conditions figures/counts_per-symbol 10 < DESeqAnalysis.R

###### Get initial conditions ######
args = commandArgs(TRUE)
countsMaxFile = args[1]   # Gene-level counts, took the max at each gene
countsUniqFile = args[2]  # Transcript-level counts, have gene symbol but
						# looks like: TP53, TP53.1, TP53.2, ... for all the 
						# possible transcripts
conditionsFile = args[3]
filePrefix = args[4]
if( length(args) == 5){
min.cutoff = as.numeric(args[5])   	# Only use genes that have values greater than this
} else{ 				# minimum cutoff for Differential Expression analysis.
	min.cutoff = 0		# If not provided, ignore all genes which are zero
}						# for all samples.
						

# Arguments didn't work, so I'm hard-coding the variables in and testing them
# countsMaxFile = "~/workspace/UCSC/pourmand/2012-spring/results/expression/counts_per-symbol-max.tab"
# countsUniqFile = "~/workspace/UCSC/pourmand/2012-spring/results/expression/counts_per-symbol-uniq.tab"
# conditionsFile = "~/workspace/UCSC/pourmand/2012-spring/results/expression/counts_per-symbol-max.conditions"
# filePrefix = "~/workspace/UCSC/pourmand/2012-spring/results/expression/figures/counts_per-symbol"
# min.cutoff = 10
####################################
# print(ls.str())

conds = factor(unlist(read.delim(conditionsFile, header=FALSE))) #t(as.matrix(read.delim(conditionsFile, header=FALSE)))
indGroup = grep("group", conds, ignore.case=TRUE, value=FALSE)
indNotGroup = grep("group", conds, ignore.case=TRUE, value=FALSE, 
	invert = TRUE)

counts.max = read.delim(countsMaxFile, row.names=1)
counts.uniq = read.delim(countsUniqFile, row.names=1)


# Make list type to access name and dataset for each of max and uniq types
counts.all = list(list(name="uniq", ds=counts.uniq), 
	list(name="max",ds=counts.max))

# counts.all = list(list(name="ssGSEA", ds=ssGSEA.no.desc))

source("DESeq_utils.R")  # Get the plotting functions

# Take out the groups and compare only individual samples
for( countData in counts.all){
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
		if( groupType == "groups"){
			# print("indGroup:")
			# print(indGroup)
			theseConds = conds[indGroup]
			theseCounts = counts.ds[,indGroup]
			# cds = newCountDataSet(, theseConds)
		} else{
			theseConds = conds[indNotGroup]
			theseCounts = counts.ds[,indNotGroup]
			# cds = newCountDataSet(counts.ds.over.min.cutoff[,indNotGroup], theseConds)
		}
		print(paste("number of features before cutoff of", min.cutoff,
			":", dim(theseCounts)[1]))
		under.min.cutoff = apply(theseCounts, 1, 
			# Want all features that have at least one sample 
			# with a value above min.cutoff
			function(x){ sum(x > min.cutoff) == 0 } )
		theseCounts = theseCounts[!under.min.cutoff,]
		print(paste("number of features after cutoff of", 
			min.cutoff, ":", dim(theseCounts)[1]))
		
		cds = newCountDataSet(theseCounts, theseConds)
		
		# print("theseConds")
		# print(theseConds)
		possibleComparisons = unique(theseConds)
		
		# Create string of this counts data type and group for easy access
		min.cutoff.str = paste("cutoff", min.cutoff, sep="-")
		filename.base = paste(filePrefix, countsName, groupType, 
			min.cutoff.str, sep="_")

		# Estimate effective library size
		cds = estimateSizeFactors(cds)
		
		# Estimate gene dispersions
		# Simon Anders' description of why "pooled" from SEQAnswers forum:
		# "Sometimes, both condition have unequal variance (for example, 
		# knock-down samples might differ strongly from each other than 
		# untreated control samples because knock-down efficiency is so hard 
		# to keep constant), and then, "per-condition" can give more power. 
		# This is why this was the default. However, I realized recently that 
		# our way of avoiding outliers (see the discussion of 
		# 'sharingMode="maximum"' in the vignette) does not work as reliably 
		# as I hoped when using "per-condition" estimation. This is why I 
		# changed the default to "pooled" and added a note about this fact to 
		# the help page. I have some ideas on how to improve this matter but 
		# pending that I recommend "pooled"."
		cds = estimateDispersions(cds, fitType="local")
		pdf(paste(filename.base, "_dispersionEstimates.pdf", sep=""))
		plotDispEsts(cds); dev.off()
		
		cds.blind = estimateDispersions(cds, method="blind", fitType="local")
		pdf(paste(filename.base, "_dispersionEstimates_blind.pdf", sep=""))
		plotDispEsts(cds.blind); dev.off()
		
		vsd = getVarianceStabilizedData(cds.blind)
		pdf(paste(filename.base, "_varianceStabilizedData_mean-vs-sd.pdf", 
		sep=""))
		meanSdPlot(vsd)
		dev.off()
		
		plotVSD(cds, vsd, filename.base)  # includes pdf and dev.off functions
		
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
				pdf(paste(filename.base, bothGroupsStr,"pval.pdf", sep="_"))
				pValHist(res$pval, paste(groupOne, "vs", groupTwo)); dev.off()
				pdf(paste(filename.base, bothGroupsStr,"padj.pdf", sep="_"))
				pValHist(res$padj, paste(groupOne, "vs", groupTwo)); dev.off()
				pdf(paste(filename.base, bothGroupsStr,
					"log2FoldChange.pdf", sep="_"))
				plotDE(res) ; dev.off()
				
				# Moderated log fold change estimates via variance stabilizing
				# transformation
				mod_lfc = (rowMeans(vsd[,conditions(cds)==groupOne, 
						drop=FALSE]) -
					rowMeans(vsd[,conditions(cds)==groupTwo, drop=FALSE]))
				lfc = res$log2FoldChange
				finite = is.finite(lfc)
				print(table(as.character(lfc[!finite]), useNA="always"))
				
				largeNumber = 7
				lfc = ifelse(finite, lfc, sign(lfc)*largeNumber)
				log.decade = 1 + round(log10(1+rowMeans(counts(cds.blind, 
					normalized=TRUE))))
				colors = colorRampPalette(c("gray", "blue"))(6)[log.decade]
				pdf(paste(filename.base, bothGroupsStr,
					"log2FoldChange-vs-ModeratedLogRatio.pdf", sep="_"))
				plot(lfc, mod_lfc, pch=20, cex=.4, asp=1, 
					col=ifelse(finite,colors,"purple"),
					main="Direct (lfc) vs Moderated log ratios (mod_lfc)",
					sub="purple were infinite in lfc, arbitrarily set for plotting")
				abline(a=0, b=1, col="#40C04040")
				dev.off()
				
				olfactoryReceptors = grep("^OR", res$id)
				
				indPval = order(res$pval)[-olfactoryReceptors][1:40]
				pdf(paste(filename.base, bothGroupsStr, 
					"top40DE-geneClustering.pdf", sep="_"))
				heatmap.2(vsd[indPval,], col=pinkogram, 
					main=paste(groupOne, "vs", groupTwo, 
					"\nmin.cutoff =",min.cutoff),
					sub="plotting method: heatmap.2")
				# heatmap(vsd[indPval,], col=pinkogram, 
				# 	main=paste(groupOne, "vs", groupTwo),
				# 	sub="plotting method: heatmap()")
				dev.off()
				
				dists = dist(t(vsd))
				pdf(paste(filename.base, bothGroupsStr, 
					"top40DE-sampleClustering.pdf", sep="_"))
				heatmap(as.matrix(dists), symm=TRUE, scale="none",
					margins=c(10,10), col=colorRampPalette(c("darkblue",
						"white"))(100),
						labRow=paste(pData(cds.blind)$condition))
				dev.off()
				
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
		write.table(vsd, paste(filename.base, "vsd.txt", sep="_"), 
			quote=FALSE, row.names=TRUE, col.names=FALSE, sep="\t")
		# using vsd for MI produces better results than with cds 
		# (more obvious divide between classes)
		# binary.vector = grep("Survivors", theseConds)
		# vsd.MI = apply(vsd, 1, mutual.inf.2, ifSurvivors)
	}
}



################## Analyze ssGSEA data ###############
ssGSEA = read.delim(projected, header=TRUE, skip=2, row.names=1)
# Remove description column
ssGSEA = as.matrix(ssGSEA[,-1])
projected = "counts_per-symbol-max.PROJ.gct"

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


genesetsToPlot = list(groups=NULL, individuals=NULL)
for( groupType in c("groups", "individuals")){
	print(paste("group type:", groupType))
	if( groupType == "groups"){
		# print("indGroup:")
		# print(indGroup)
		theseConds = conds[indGroup]
		thesessGSEA = ssGSEA[,indGroup]
		colnames(thesessGSEA) = groupTreatmentNames
		treatmentColors = groupTreatmentColors
		# cds = newCountDataSet(, theseConds)
	} else{
		theseConds = conds[indNotGroup]
		thesessGSEA = ssGSEA[,indNotGroup]
		colnames(thesessGSEA) = individualTreatmentNames
		treatmentColors = individualTreatmentColors
		# cds = newCountDataSet(counts.ds.over.min.cutoff[,indNotGroup], theseConds)
	}
	genesetsInd = NULL
	possibleComparisons = unique(theseConds)
	# compared = NULL
	for( treatmentType in possibleComparisons){
		print(paste("Finding MI of geneset and", treatmentType))
		binary.vector = as.numeric(theseConds == treatmentType)
		MI = apply(thesessGSEA, 1, mutual.inf.2, binary.vector)
		MI.sorted = sort(MI, decreasing=TRUE)
		MI.order = order(MI)
		first5 = 1:5
		last5 = (length(MI.sorted)-4):length(MI.sorted)
		genesetsInd = c(genesetsInd, MI.order[first5], MI.order[last5])
		
		if( groupType == "groups"){
			genesetsToPlot$groups = c(genesetsToPlot$groups, 
				names(MI.sorted)[first5],
		 		names(MI.sorted)[last5])
		} else{
			genesetsToPlot$individuals = c(genesetsToPlot$individuals,
				names(MI.sorted)[first5],
			 	names(MI.sorted)[last5])
		}

		outFile = paste(filePrefix, groupType, treatmentType, 
			"ssGSEA_MI.txt", sep="_")
		
		write.table(MI.sorted, outFile, quote=FALSE,sep="\t")
	}
	# quartz(height=11,width=8.5)
	pdf(paste(filePrefix, groupType,
		 "ssGSEA_MI_unclusteredSamples.pdf",sep="_"),
		height=11, width=8.5)
	heatmap.2(thesessGSEA[genesetsInd,], col=pinkogram,
		main="Genesets with extreme mutual information to treatment type",
		margins=c(12,15),
		ColSideColors=treatmentColors, Colv=FALSE, #Rowv=FALSE,
		#dendrogram="none"
		)
	dev.off()
}

