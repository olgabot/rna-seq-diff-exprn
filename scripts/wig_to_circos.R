#!/usr/bin/Rscript
options(error=utils::recover)

## Created by Olga Botvinnik on 07-08-2012.
##
## Usage: ./wig_to_circos.R <wig file> <circos file>
## Example run: 
## scripts/wig_to_circos.R test-data/hg19_gc1000Base.txt test-data/hg19_gc_content_circos.txt

# NOTE: .wig is 1-based while circos is 0-based

###### Get initial conditions ######
args = commandArgs(TRUE)
# wig = read.delim("hg19_gc_data_points.tab", row.names=NULL, 
# 	stringsAsFactors=FALSE)
wig = read.delim(args[1], row.names=NULL, 
	colClasses=c("character", "numeric"), 
	stringsAsFactors=FALSE, header=FALSE)
circosFilename = args[2]
####################################
# print(ls.str())

write.table(matrix(c("#chr", "chrStart", "chrEnd", "value"), nrow=1), 
	file=circosFilename, sep="\t", row.names=FALSE, 
	col.names=FALSE, quote=FALSE)

# Get the indices of ALL chromosome starts
indChromAll = grep("chrom=chr", wig[,1], perl=TRUE)

# Get the indices of the REAL chromosome starts
# only numbered chromosomes (no Un) and X, Y
# (No mitochondrial DNA, ie chrM)
indChromClean = grep("chrom=chr[\\dXY]+(?= )", wig[,1], perl=TRUE)

# Get the bin size for each chromosome,
# since they're allowed to be different.
binSizes = as.numeric(sapply(indChromClean, function(i){
	unlist(strsplit(wig[i,1], split="span="))[2]
	}))

# Get the chromosome names
chromosomes = sapply(indChromClean, function(i) {
	unlist(strsplit(unlist(strsplit(wig[i,1], 
		split="chrom=", perl=TRUE))[2], split=" span="))[1]
	})

# Have to assign to some variable so it doesn't output a bunch of
# R-formatted lists to standard out.
x = sapply(1:length(indChromAll), function(j){
	if(indChromAll[j] %in% indChromClean){
		i = which(indChromClean == indChromAll[j])
		chr = chromosomes[i]
		binSize = binSizes[i]
	
		# Get the start and end of this chromosome's range
		indChromStart = indChromAll[j] + 1
		indChromEnd = ifelse(j < length(indChromAll), 
			indChromAll[j+1]-1,
			nrow(wig))
	
		# Get the start values as defined by the .wig file
		wigStart = as.numeric(wig[indChromStart:indChromEnd,1])
	
		# Turn the .wig 1-based start values into circos 0-based
		# start values by subtracting 1
		# Note: formatC can force scientific notation into integers,
		#      	and Circos requires integers for the chromosome
		# 		start and end sites
		chrStart = wigStart - 1
	
		# Need to subtract 1 from the bin size so we get nice 
		# ranges like: 100-104, 105-109
		# Note: formatC can force scientific notation into integers,
		#      	and Circos requires integers for the chromosome
		# 		start and end sites
		chrEnd = chrStart + binSize - 1
	
		# as.numeric forces weird characters into NAs so it 
		# helps cut out any strangeness in the file
		value = as.numeric(wig[indChromStart:indChromEnd,2])
	
		# Make a data table that has the chromosome name in the first 
		# column, the start of the range in the second column, the end
		# of the range in the third column, and the value of the range
		# in the fourth column
		dataTable = data.frame(cbind(chr, formatC(chrStart, format="d"), 
			formatC(chrEnd, format="d"), value))
		
		# Get rid of any NA values
		indNA = is.na(dataTable[,4])
		dataTable = dataTable[!indNA,]
	
		write.table(dataTable, 
			file=circosFilename, append=TRUE,
			sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
	}
	})