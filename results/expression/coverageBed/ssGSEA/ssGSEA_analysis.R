#!/share/apps/bin/R

source("make.pinkogram.R")

############# Import arguments
args = commandArgs(TRUE)
# expression.table = read.delim(args[1], header=TRUE)
nFeatures = 20
##############################

library(MASS)

############## This does not run!!!!!!!!!!!!!!!!! ###############

# gene.names = make.unique(expression.table[,1])
# gene.names = rownames(counts(cds))
# expression.counts = expression.table[,-1]
	
heatmap.olga = function(mat, feature.order, pdf.file, 
		data.type=c("positive.only", "positive.and.negative")){
	# When using DESeq-type data, variance stabilized data (vsd) works best
	# You can get vsd via :
	# cds = newCountDataSet(myCounts, conds)
	# cds = estimateSizeFactors(cds)
	# cds.blind = estimateDispersions(cds, method="blind")
	# vsd = getVarianceStabilizedData(cds.blind)
	
	pinkogram = make.pinkogram()
	nColors = length(pinkogram)

	# Normalize counts with a function that normalizes positive-count data
	if(data.type == "positive.only"){
	heatm = ceiling(nColors * (mat - min(mat))/(1.001*(max(mat) - min(mat))))  
}	
	else{
		# only appropriate for data with negative values
		max.c <- max(max(mat), -min(mat))
		heatm <- ceiling(nColors * (mat - (- max.c))/(1.001*(max.c - (- max.c))))
	}
	heatm.top = mat[feature.order,]
	
	nFeatures = length(feature.order)
	nSamples = dim(heatm.top)[2]
	
	pdf(pdf.file)
	image(1:nSamples, 1:nFeatures, t(heatm.top), zlim = c(0, nColors),
	      col=pinkogram, axes=FALSE, main ="", sub="", xlab="", ylab="")
	axis(2, at=1:nFeatures, labels=rev(rownames(heatm.top)), adj= 0.5,
	     tick=FALSE, las = 1, cex.axis=0.50, font.axis=1, line=-1)
	axis(1, at=1:nSample, labels=colnames(heatm.top), adj= 0.5, tick=FALSE, las = 3,
	     cex.axis=0.60, font.axis=1, line=-1)
	dev.off()
}




mutual.inf.2 <- function(x, y, n.grid=20, 
		normalize.by ="HXY", # Whether to normalize by HXY, HX, or HY
		pos.and.neg = T) {
	# x and y can be binary or continuous
	# If there is not sufficient variation in x and y, 
	# will take the standard deviation as the bandwidth
	# (IQR finds the inter-quartile range of the data vector)

	## Add random noise if vectors are constant to prevent errors in bcv 
	## (because the inter-quartile range of a constant vector is 0, and the IQR is 
	## used in calculating the bandwidth)
	y = unlist(y); y = as.vector(y)
	x = unlist(x); x = as.vector(x)
	if( length(unique(y)) == 1 || length(unique(x)) == 1){ return(0) }
	if( sd(y) == 0){
		y = y + runif(n=length(y), min=mean(y)-0.001, max=mean(y)+0.001)
	}
	if( sd(x) == 0){
		x = x + runif(n=length(x), min=mean(x)-0.001, max=mean(x)+0.001)
	}

	## Using suppressWarnings because bcv(.) gets mad if the minimum is 
	## on one side of the vector
	kde2d.xy <- kde2d(x, y, n = n.grid, h = c(suppressWarnings(bcv(x)),
		suppressWarnings(bcv(y))) )
	PXY <- kde2d.xy$z/sum(kde2d.xy$z)

	PX <- rowSums(PXY)#apply(PXY, MARGIN=1, sum)
	PX <- PX/sum(PX)
	HX = -sum(PX * log2(PX))
	PX <- matrix(PX, nrow=n.grid, ncol=n.grid)

	PY <- colSums(PXY)#apply(PXY, MARGIN=2, sum)
	PY <- PY/sum(PY)
	HY = -sum( PY * log2(PY))
	PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)

	## Ignore NaN's that are a result of underflow (experimentally derived)
	HXY <- - sum(PXY * log2(PXY), na.rm=TRUE)

	MI.norm =  ifelse(pos.and.neg, sign(cor(x, y)), 1) * ((HX + HY)/HXY - 1)
	return( MI.norm )#list(MI=MI, HXY=HXY))
}


# May need this:
## Tissue Legend
# ifMultipleTissues = FALSE
# if(ifMultipleTissues){
#   #browser()
#   par(mar = c(0, 0, 0, 0))
#   plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
#   legend(x=0, y= 1, #horiz = T, x.intersp = 0.5, y.intersp=.25,
#          legend=tissue.names, bty="n", xjust=0, yjust= 1,
#          fill = tissue.colors, cex = 0.8, #pt.cex=1.75,
#          ncol=4)
# }
