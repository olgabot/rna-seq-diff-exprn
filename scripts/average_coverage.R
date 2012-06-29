#!/usr/bin/Rscript
options(error=utils::recover)
##
## average_coverage
## Given an input coverage file in circos form, ie:
# # coverage stats: mean=4748 std dev=313
# chr1    14400   14404   2
# chr1    14404   14406   3
# chr1    14406   14407   8
# chr1    14407   14408   13
# chr1    14408   14409   21
# chr1    14409   14410   39
# chr1    14410   14411   54
# chr1    14411   14412   82
# chr1    14412   14413   101
## Take every chromosome and take the median coverage value per interval.
## 
## Created by Olga Botvinnik on 2012-05-22.
## 
##
## Usage: ./average_coverage.R <coverageFile> <interval> <outFile> [optional: method]
## Example run: ./average_coverage.R genome_coverage_bedtools.txt 100000 genome_coverage_bedtools_100000.txt

###### Get initial conditions ######
args = commandArgs(TRUE)
coverageFile = args[1]
interval = as.numeric(args[2])
outFile = args[3]
if( length(args) == 4){
  avgMethod = args[4]
} else{ avgMethod = "median"}
####################################

# if avgMethod is not recognized, you will see an error like:
# Error in eval(expr, envir, enclos) : object 'average' not found
avg = eval(parse(text=avgMethod))
  
coverage = read.delim(coverageFile, header=FALSE, skip=1)
print(paste('time to read file:', coverageFile))
print(proc.time())
chromosomes = unique(as.character(coverage[,1]))

avgCoverage = NULL
for(chr in chromosomes){
	chrInd = which(coverage[,1] == chr)
	chrMin = min(coverage[chrInd,2])  # take the minimum start site
	chrMax = max(coverage[chrInd,3])  # take the maximum end site
	# calculate the number of intervals
	numIntervals = ceiling((chrMax-chrMin)/interval)
        print(paste("number of", interval,
                    "bp intervals for", chr,
                    numIntervals))
        thisChrCoverage = apply(matrix(1:numIntervals), 1,
          function(i){
		thisIntervalMin = chrMin + (i-1)*interval
		thisIntervalMax = min(chrMin + i*interval-1, chrMax)
		thisIntervalInd = which((coverage[,2] > thisIntervalMin) &
                  (coverage[,3] < thisIntervalMax))
                return(c(chr,
                  thisIntervalMin, thisIntervalMax,
                  avg(coverage[thisIntervalInd,4])))
          }
          )
        avgCoverage = cbind(avgCoverage, thisChrCoverage)
        print(paste("finished chromosome:", chr))
        print(proc.time())
}
avgCoverage = t(avgCoverage)
write.table(avgCoverage, outFile, row.names=FALSE, col.names=FALSE,
            quote=FALSE, sep="\t")
