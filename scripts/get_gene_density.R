#!/usr/bin/Rscript
options(error=utils::recover)
##
## 
##
## Created by Olga Botvinnik on 07-07-2012.
## Copyright (c)  . All rights reserved.
##
## Collaborators: 
## Usage: ./.R <arg1> <arg2> <arg3>
## Example run: ./get_gene_density.R hg19_ucsc_knownCanonical.tab hg19_gene_density.txt

###### Get initial conditions ######
args = commandArgs(TRUE)
knownCanonical = read.delim(args[1])
outFile = args[2]

# Measure gene density per megabase
interval = 1e6 
####################################
# print(ls.str())

chromosomes = unique(knownCanonical[,1])

# Remove mitochondrial chromosome (chrM) and anything with underscores
chromosomes = grep("_", chromosomes, invert=TRUE, value=TRUE)
chromosomes = grep("chrM", chromosomes, invert=TRUE, value=TRUE)

# a = sapply(chromosomes, function(chr){
# 	ind = which(knownCanonical[,1] = chr)

# 	})
write("#chr\tchrStart\tchrEnd\tdensity", file=outFile, append=FALSE)
for( chr in chromosomes ){
	ind = which(knownCanonical[,1] == chr)
	chrMin = min(knownCanonical[ind,2])
	chrMax = max(knownCanonical[ind,3])
	numIntervals = ceiling((chrMax-chrMin)/interval)
	for( i in 1:numIntervals ){

    # Get the value of the start and end positions on the chromosome
    # for this interval
		thisIntervalMin = chrMin + (i-1)*interval
		thisIntervalMax = min(chrMin + i*interval-1, chrMax)

    # Get the indices of genes in this interval
		thisIntervalInd = which((knownCanonical[,2] > thisIntervalMin) &
                  (knownCanonical[,3] < thisIntervalMax))

    # Get the density of this interval
		thisIntervalDensity = length(thisIntervalInd)/(
			thisIntervalMax - thisIntervalMin)

		write.table(matrix(c(chr, thisIntervalMin, thisIntervalMax, 
			thisIntervalDensity), nrow=1), file=outFile, append=TRUE,
      sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
	}
}