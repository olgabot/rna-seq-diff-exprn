#!/usr/bin/Rscript
options(error=utils::recover)
##
## NOISeqAnalysis
##
## Created by Olga Botvinnik on 2012-05-27.
## Copyright (c) 2012 __MyCompanyName__. All rights reserved.
##
## Collaborators: 
## Usage: ./NOISeqAnalysis.R <arg1> <arg2> <arg3>
## Example run: ./NOISeqAnalysis.R
## To run on annoying servers that don't interpret R scripts properly:
# R --slave < NOISeqAnalysis.R

###### Get initial conditions ######
args = commandArgs(TRUE)
filePrefix = "/home/obot/single-cell/results/expression/coverageBed/NOISeq"
uniqLengthFile = "/home/obot/bed_and_gtf/hg19_ucsc-genes_transcript_info_symbolUniq-length.txt"
maxLengthFile = "/home/obot/bed_and_gtf/hg19_ucsc-genes_transcript_info_symbolMax-length.txt"
# took the maximum count per gene symbol, so use symbolLengthFile
countsMaxFile = "/home/obot/single-cell/results/expression/coverageBed/counts_per-symbol-max.tab"  
# Each UCSC transcript ID
countsUniqFile = "/home/obot/single-cell/results/expression/coverageBed/counts_per-symbol-uniq.tab"
min.cutoff = 0
conditionsFile = "/home/obot/single-cell/results/expression/conditions.tab"
####################################
# print(ls.str())

conditions = unlist(read.delim(conditionsFile, header=FALSE))
indGroup = grep("group", conditions, ignore.case=TRUE, value=FALSE)
indNotGroup = grep("group", conditions, ignore.case=TRUE, value=FALSE, 
	invert = TRUE)

setwd("/share/reference/hg19")
lengths = read.delim("hg19_ucsc-genes_transcript-info_id-length-symbol_noheader.txt", header=FALSE)
symbols = unique(lengths[,3])
symbolLengths = sapply(symbols, function(symb){max(lengths[lengths[,3]==symb,2])})

uniqLengths = cbind(make.unique(as.character(lengths[,3])), lengths[,2])
write.table(uniqLengths, uniqLengthFile, 
	row.names=FALSE, col.names=FALSE, 
	quote=FALSE, sep="\t")

maxLengths = cbind(as.character(symbols), as.character(symbolLengths))
write.table(maxLengths, maxLengthFile, row.names=FALSE, col.names=FALSE, 
	quote=FALSE, sep="\t")

setwd("/home/obot/single-cell/results/expression/coverageBed/NOISeq")
source("/home/obot/single-cell/scripts/expression/NOISeq/noiseq.r")

# counts.max = read.delim(countsMaxFile, row.names=1)
# counts.uniq = read.delim(countsUniqFile, row.names=1)

# Make list type to access name and dataset for each of max and uniq types
counts.all = list(list(name="uniq", dataFile=countsUniqFile, lengths=readInfo(uniqLengthFile, header=FALSE)), 
	list(name="max",dataFile=countsMaxFile, lengths=readInfo(maxLengthFile, header=FALSE)))

normalizationTypes = c("rpkm",  # default of NOISeq
						"uqua", # Upper-quartile
						"tmm",  # Trimmed mean of m
						"n")    # No normalization

# counts.all = list(list(name="ssGSEA", ds=ssGSEA.no.desc))

# Take out the groups and compare only individual samples
for( countData in counts.all){
	# Print out the data type
	print(paste("data type:", countData$name))
	
	dataFile = countData$dataFile
	dataName = countData$name
	lengths = countData$lengths
	
	# print("conditions:")
	# print(conds)
	
	# Iterate over comparing 5 or 6 single-cell samples vs comparing 2 groups
	# of 2 or 3 samples in each data type. 
	# e.g. Untreated1 = samples 3, 4, 6 & Untreated2 = samples 7, 8, 9
	# compared vs HighDose1 = 34, 35, 36 & HighDose2 = samples 38, 39, 40
	for( groupType in c("individuals", "groups")){
		# Print out the group type
		thisFilePrefix.orig = paste(filePrefix, dataName, groupType, sep="_")
		print(paste("group type:", groupType))
		if( groupType == "groups"){
			# print("indGroup:")
			# print(indGroup)
			theseConds = conds[indGroup]
			# theseCounts = counts.ds[,indGroup]
			possibleComparisons = unique(theseConds)

			# Create string of this counts data type and group for easy access
			# min.cutoff.str = paste("cutoff", min.cutoff, sep="-")
			# filename.base = paste(filePrefix, countsName, groupType, 
			# 	min.cutoff.str, sep="_")

			compared = NULL
			for( groupOne in possibleComparisons){
				# Need to add 1 to ignore 1st column, which is the
				# gene symbol
				cond1 = which(conditions==groupOne)+1
				for( groupTwo in possibleComparisons){
					# Need to add 1 to ignore 1st column, which is the
					# gene symbol					
					cond2 = which(conditions==groupTwo)+1
					
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
					
					theseData = readData(file=dataFile, cond1=cond1,
						cond2=cond2, header=TRUE)
						for(normType in normalizationTypes){
							thisFilePrefix = paste(thisFilePrefix.orig, 
								normType, bothGroupsStr.1,
								 sep="_")
							noiseqResults = noiseq(theseData[[1]], theseData[[2]],
								repl="bio", norm=normType, long=lengths)
							topResultsTable = cbind(
								noiseqResults$Ms[noiseqResults$deg],
								noiseqResults$Ds[noiseqResults$deg],
								noiseqResults$probab[noiseqResults$deg])
							topHeader = paste("Name", "MvalueInSignal",
							"DvalueInSignal", "ProbDiffExprn", collapse="\t")
							topFile = paste(thisFilePrefix, "topDE.txt", 
								sep="_")
							write(topHeader, topFile)
							write.table(topResultsTable, file=topFile, 
								row.names=TRUE, col.names=FALSE, 
								append=TRUE, quote=FALSE, sep="\t")

							allResultsTable = 	cbind(
									noiseqResults$Ms,
									noiseqResults$Ds,
									noiseqResults$probab)
							allHeader = topHeader
							allFile = paste(thisFilePrefix, "all.txt", 
								sep="_")
							write(allHeader, allFile)
							write.table(allResultsTable, file=allFile, 
								row.names=TRUE, col.names=FALSE, 
								append=TRUE, quote=FALSE, sep="\t")

							noiseTable = cbind(noiseqResults$Mn,
							noiseqResults$Dn)
							noiseFile = paste(thisFilePrefix, "noise.txt", 
								sep="_")
							write(paste("MvaluesInNoise", "DvaluesInNoise", 
								collapse="\t"), noiseFile)
							write.table(noiseTable, file=noiseFile, 
								row.names=FALSE, col.names=FALSE, 
								append=TRUE, quote=FALSE, sep="\t")
						}
				}
			}
		} else{
			theseConds = conditions[indNotGroup]
			# theseCounts = counts.ds[,indNotGroup]
			possibleComparisons = unique(theseConds)

			compared = NULL
			for( groupOne in possibleComparisons){
				# Need to add 1 to ignore 1st column, which is the
				# gene symbol
				cond1 = which(conditions==groupOne)+1
				for( groupTwo in possibleComparisons){
					# Need to add 1 to ignore 1st column, which is the
					# gene symbol					
					cond2 = which(conditions==groupTwo)+1
					
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
					
					theseData = readData(file=dataFile, cond1=cond1,
						cond2=cond2, header=TRUE)
					for(normType in normalizationTypes){
						thisFilePrefix = paste(thisFilePrefix.orig, normType,
							bothGroupsStr.1,
							 sep="_")
						noiseqResults = noiseq(theseData[[1]], theseData[[2]],
							repl="bio", norm=normType, long=lengths)
						topResultsTable = cbind(
							noiseqResults$Ms[noiseqResults$deg],
							noiseqResults$Ds[noiseqResults$deg],
							noiseqResults$probab[noiseqResults$deg])
						topHeader = paste("Name", "MvalueInSignal",
						"DvalueInSignal", "ProbDiffExprn", collapse="\t")
						topFile = paste(thisFilePrefix, "topDE.txt", 
							sep="_")
						write(topHeader, topFile)
						write.table(topResultsTable, file=topFile, 
							row.names=TRUE, col.names=FALSE, 
							append=TRUE, quote=FALSE, sep="\t")
						
						allResultsTable = 	cbind(
								noiseqResults$Ms,
								noiseqResults$Ds,
								noiseqResults$probab)
						allHeader = topHeader
						allFile = paste(thisFilePrefix, "all.txt", 
							sep="_")
						write(allHeader, allFile)
						write.table(allResultsTable, file=allFile, 
							row.names=TRUE, col.names=FALSE, 
							append=TRUE, quote=FALSE, sep="\t")
						
						noiseTable = cbind(noiseqResults$Mn,
						noiseqResults$Dn)
						noiseFile = paste(thisFilePrefix, "noise.txt", 
							sep="_")
						write(paste("MvaluesInNoise", "DvaluesInNoise", 
							collapse="\t"), noiseFile)
						write.table(noiseTable, file=noiseFile, 
							row.names=FALSE, col.names=FALSE, 
							append=TRUE, quote=FALSE, sep="\t")
					}
				}
			}
		}
	}
}


