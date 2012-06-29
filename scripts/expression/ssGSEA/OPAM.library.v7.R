library(CvM2SL2Test)
library(MASS)
library(verification)

options(error=utils::recover)
t1 <<- proc.time()

OPAM.Projection.3 <- function(
		data.array,
		gene.names,
		n.cols,
		n.rows,
		weight = 0,
		statistic = "Kolmogorov-Smirnov",  # "Kolmogorov-Smirnov", # "Kolmogorov-Smirnov", "Cramer-von-Mises",
		# "Anderson-Darling", "Zhang_A", "Zhang_C", "Zhang_K",
		# "area.under.RES", or "Wilcoxon"
		gene.set,
		nperm = 200,
		correl.type  = "rank",
		save.gene.set.rankings = FALSE)             # "rank", "z.score", "symm.rank"
# Runs a 2-3x faster (2-2.5x for ES statistic and 2.5-3x faster for area.under.ES statsitic)
# version of GSEA.EnrichmentScore.5 internally that avoids overhead from the function call.
{
	
	ES.vector <- vector(length=n.cols)
	NES.vector <- vector(length=n.cols)
	p.val.vector <- vector(length=n.cols)
	correl.vector <- vector(length=n.rows, mode="numeric")
	
# Compute ES score for signatures in each sample
#	browser()
#   print("Computing GSEA.....")
	if(save.gene.set.rankings){
		gene.ranks.matrix = matrix(ncol = n.cols,  # n.cols = n.samples
				nrow = length(gene.set), # size of gene set
				dimnames=list(gene.set, colnames(data.array)))
		expression.matrix = matrix(ncol = n.cols, nrow = n.rows,  # used to be n.row = number of total genes
				dimnames=list(NULL, colnames(data.array))) 
	}
	
	phi <- array(0, c(n.cols, nperm))
	gene.set2 <- match(gene.set, gene.names)
	N <- n.rows
	Nh <- length(gene.set2) 
	Nm <-  N - Nh 
	for (sample.index in 1:n.cols) {
		gene.list <- order(data.array[, sample.index], decreasing=T)
		
		## Redundant, moved to outside of loop to save time - Olga 07/07/2011
#		gene.set2 <- match(gene.set, gene.names)
		
		#      print(paste("Computing observed enrichment for UP signature in sample:", sample.index, sep=" ")) 
		
		
		if (weight == 0) {
			correl.vector <- rep(1, n.rows)
		} else if (weight > 0) {
			if (correl.type == "rank") {
				correl.vector <- data.array[gene.list, sample.index]
			} else if (correl.type == "symm.rank") {
				correl.vector <- data.array[gene.list, sample.index]
				correl.vector <- ifelse(correl.vector > correl.vector[ceiling(n.rows/2)], 
						correl.vector,
						correl.vector + correl.vector - correl.vector[ceiling(n.rows/2)]) 
			} else if (correl.type == "z.score") {
				x <- data.array[gene.list, sample.index]
				correl.vector <- (x - mean(x))/sd(x)
			}
		}
		### Olga's Speedups ###
#		ptm.new = proc.time()
		tag.indicator <- sign(match(gene.list, gene.set2, nomatch=0))    # notice that the sign is 0 (no tag) or 1 (tag) 
		no.tag.indicator <- 1 - tag.indicator
		### Redundant, can be done once. Moved these to outside of the loop to save time - Olga 07/07/2011
#		N <- length(gene.list)    
#		Nh <- length(gene.set2) 
#		Nm <-  N - Nh 
		orig.correl.vector <- correl.vector
		if (weight == 0) correl.vector <- rep(1, N)   # unweighted case
		ind = which(tag.indicator==1)
		if(save.gene.set.rankings){
#			browser()
			expression.matrix[,sample.index] = abs(correl.vector)^weight #[match(gene.names[gene.list][ind], gene.set)]
			#data.array[gene.list,sample.index]
#			expression.matrix[,sample.index] = sapply(data.array[,sample.index], function(x) (x - mean(x))/sd(x)) 
			#expression.matrix[,sample.index] = correl.vector #
			
			# index of each gene in the gene set, in order
			gene.ranks.matrix[,sample.index] = ind[match(gene.names[gene.list][ind], gene.set)]
			#gene.ranks.matrix[,sample.index] = gene.list[gene.set2]
		}
		correl.vector <- abs(correl.vector[ind])^weight
#		browser()
		
		
		
		sum.correl = sum(correl.vector)
		
		gaps = (c(ind-1, N) - c(0, ind))  # gaps between ranked pathway genes
		down = gaps/Nm
		if( sum.correl != 0){
			up = correl.vector/sum.correl     # "up" represents the peaks in the mountain plot
			
			RES = cumsum(c(up,up[Nh])-down)
			valleys = RES[1:Nh]-up            # "valleys" is every time the mountain plot goes down
			max.ES = max(RES)
			min.ES = min(valleys)
		} else{
			max.ES <- min.ES <- 0
			RES = vector(length=length(RES))
			valleys <- up <- vector(length=length(valleys))
		}
		
		if( statistic == "Kolmogorov-Smirnov" ){
			if( max.ES > -min.ES ){
				ES <- signif(max.ES, digits=5)
				arg.ES <- which.max(RES)
			} else{
				ES <- signif(min.ES, digits=5)
				arg.ES <- which.min(RES)
			}
		} else if( statistic == "area.under.RES"){
			if( max.ES > -min.ES ){
				arg.ES <- which.max(RES)
			} else{
				arg.ES <- which.min(RES)
			}
			gaps = gaps+1
			RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES[1:Nh]) - c(valleys,0) ) * (gaps)
			ES = sum(RES)
		} else if (statistic == "Wilcoxon") {
			# Wilcoxon test score
			library(exactRankTests)
			seq.index <- seq(1, N)
			gene.set.ranks <- seq.index[tag.indicator == 1]
			gene.set.comp.ranks <- seq.index[tag.indicator == 0]
			W <- wilcox.exact(x=gene.set.ranks, y =gene.set.comp.ranks, alternative = "two.sided", mu = 0,
					paired = FALSE, exact = F, conf.int = T, conf.level = 0.95)
			ES <- log(1/W$p.value)
			arg.ES = RES = NULL
			#return(list(ES = ES, arg.ES = NULL, RES = NULL, indicator = tag.indicator))
		} else if (statistic == "Cramer-von-Mises") {
			# Based on T. W. Anderson Annals of Mathematical Statistics 1962.
			# Modified to separate positive and negative enrichment parts
			RES <- c(valleys,0) * (gaps) + 0.5*( c(0,RES[1:Nh]) - c(valleys,0) ) * (gaps)#Fn - F0
			X <- RES^2
			X_p <- X[RES >= 0]
			X_n <- X[RES < 0]
			ES_p <- ((Nh * Nm)/N) * sum(X_p)
			ES_n <- ((Nh * Nm)/N) * sum(X_n)
			if (ES_p > ES_n) {
				ES <- signif(ES_p, digits = 5)
				arg.ES <- which.min(abs(X - max(X_p)))
			} else {
				ES <- - signif(ES_n, digits=5)
				arg.ES <- which.min(abs(X - max(X_n)))
			}
		}
		GSEA.results = list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator)
#		new.time <<- new.time + (proc.time() - ptm.new)
		### End Olga's Speedups ###
		#GSEA.results <- GSEA.EnrichmentScore5(gene.list=gene.list, gene.set=gene.set2,
		#		statistic = statistic, alpha = weight, correl.vector = correl.vector)
		ES.vector[sample.index] <- GSEA.results$ES
		
		if (nperm == 0) {
			NES.vector[sample.index] <- ES.vector[sample.index]
			p.val.vector[sample.index] <- 1
		} else {
			for (r in 1:nperm) {
				reshuffled.gene.labels <- sample(1:n.rows)
				if (weight == 0) {
					correl.vector <- rep(1, n.rows)
				} else if (weight > 0) {
					correl.vector <- data.array[reshuffled.gene.labels, sample.index]
				} 
#				GSEA.results <- GSEA.EnrichmentScore5(gene.list=reshuffled.gene.labels, gene.set=gene.set2,
#						statistic = statistic, alpha = weight, correl.vector = correl.vector)
				### Olga's Additions ###
				tag.indicator <- sign(match(reshuffled.gene.labels, gene.set2, nomatch=0))    
				# notice that the sign is 0 (no tag) or 1 (tag) 
				no.tag.indicator <- 1 - tag.indicator 
				N <- length(reshuffled.gene.labels) 
				Nh <- length(gene.set2) 
				Nm <-  N - Nh 
#   orig.correl.vector <- correl.vector
				if (weight == 0) correl.vector <- rep(1, N)   # unweighted case
				ind <- which(tag.indicator==1)
				correl.vector <- abs(correl.vector[ind])^weight   
				
				sum.correl <- sum(correl.vector)
				up = correl.vector/sum.correl
				gaps = (c(ind-1, N) - c(0, ind))
				down = gaps/Nm
				
				RES = cumsum(c(up,up[Nh])-down)
				valleys = RES[1:Nh]-up
				
				max.ES = max(RES)
				min.ES = min(valleys)
				
				if( statistic == "Kolmogorov-Smirnov" ){
					if( max.ES > -min.ES ){
						ES <- signif(max.ES, digits=5)
						arg.ES <- which.max(RES)
					} else{
						ES <- signif(min.ES, digits=5)
						arg.ES <- which.min(RES)
					}
				}
				
				if( statistic == "area.under.RES"){
					if( max.ES > -min.ES ){
						arg.ES <- which.max(RES)
					} else{
						arg.ES <- which.min(RES)
					}
					gaps = gaps+1
					RES = c(valleys,0) * (gaps) + 0.5*( c(0,RES[1:Nh]) - c(valleys,0) ) * (gaps)
					ES = sum(RES)
				}
				
				GSEA.results = list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator)
				### End Olga's Additions ###
				phi[sample.index, r] <- GSEA.results$ES
			}
			if (ES.vector[sample.index] >= 0) {
				pos.phi <- phi[sample.index, phi[sample.index, ] >= 0]
				if (length(pos.phi) == 0) pos.phi <- 0.5
				pos.m <- mean(pos.phi)
				NES.vector[sample.index] <- ES.vector[sample.index]/pos.m
				s <- sum(pos.phi >= ES.vector[sample.index])/length(pos.phi)
				p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
			} else {
				neg.phi <-  phi[sample.index, phi[sample.index, ] < 0]
				if (length(neg.phi) == 0) neg.phi <- 0.5 
				neg.m <- mean(neg.phi)
				NES.vector[sample.index] <- ES.vector[sample.index]/abs(neg.m)
				s <- sum(neg.phi <= ES.vector[sample.index])/length(neg.phi)
				p.val.vector[sample.index] <- ifelse(s == 0, 1/nperm, s)
			}
		}
	}
#	if(save.gene.set.rankings){
	##		ls.str()
#		browser()
#	}
	if(save.gene.set.rankings){
		return(list(ES.vector = ES.vector, NES.vector =  NES.vector, p.val.vector = p.val.vector, 
						gene.ranks.matrix = gene.ranks.matrix, expression.matrix = expression.matrix))
	} else{
		return(list(ES.vector = ES.vector, NES.vector =  NES.vector, p.val.vector = p.val.vector))
	}
	
	
} # end of OPAM.Projection.3



OPAM.project.dataset.3 <- function( 
		input.ds,
		output.ds,
		gene.set.databases,
		gene.set.selection  = "ALL",  # "ALL" or list with names of gene sets
		sample.norm.type    = "rank",  # "rank", "log" or "log.rank"
		weight              = 0.75,
		statistic           = "area.under.RES",
		output.score.type   = "ES",  # "ES" or "NES". NES is not recommended for REVEALER as it normalizes
		# along a single sample and occludes resolution across samples.
		nperm               = 200,  # number of random permutations for NES case
		combine.mode        = "combine.off",  # "combine.off" do not combine *_UP and *_DN versions in 
		# a single score. "combine.replace" combine *_UP and 
		# *_DN versions in a single score that replaces the individual
		# *_UP and *_DN versions. "combine.add" combine *_UP and 
		# *_DN versions in a single score and add it but keeping 
		# the individual *_UP and *_DN versions.
		correl.type  = "rank",
		save.gene.set.rankings = FALSE,    ## for use in Evaluate.Results.3
		sample.names.for.ordering = NULL,    ## Wanted to use this in Evaluate.Results.3 but
		rankings.data.dir = "",
		target= NULL
## don't have sample names at that point
)             # "rank", "z.score", "symm.rank"
{ #----------------------------------------------------------------------------------------
	
	# Load libraries
	library(gtools)
	library(verification)
	library(RColorBrewer)
	
	# Read input dataset
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	gene.names <- dataset$row.names
	gene.descs <- dataset$descs
	sample.names <- dataset$names
	Ns <- length(m[1,])
	Ng <- length(m[,1])
	temp <- strsplit(input.ds, split="/") # Extract input file name
	s <- length(temp[[1]])
	input.file.name <- temp[[1]][s]
	temp <- strsplit(input.file.name, split=".gct")
	input.file.prefix <-  temp[[1]][1]
	
	# Sample normalization
	
	if (sample.norm.type == "rank") {
		for (j in 1:Ns) {  # column rank normalization 
			m[,j] <- rank(m[,j], ties.method = "average")
		}
		m <- 10000*m/Ng
	} else if (sample.norm.type == "log.rank") {
		for (j in 1:Ns) {  # column rank normalization 
			m[,j] <- rank(m[,j], ties.method = "average")
		}
		m <- log(10000*m/Ng + exp(1))
	} else if (sample.norm.type == "log") {
		m[m < 1] <- 1
		m <- log(m + exp(1))
	}
	
	# Read gene set databases
	
	max.G <- 0
	max.N <- 0
	for (gsdb in gene.set.databases) {
		GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
		max.G <- max(max.G, max(GSDB$size.G))
		max.N <- max.N +  GSDB$N.gs
	}
	N.gs <- 0
	gs <- matrix("null", nrow=max.N, ncol=max.G)
	gs.names <- vector(length=max.N, mode="character")
	gs.descs <- vector(length=max.N, mode="character")
	size.G <- vector(length=max.N, mode="numeric")
	start <- 1
	for (gsdb in gene.set.databases) {
		GSDB <- Read.GeneSets.db(gsdb, thres.min = 2, thres.max = 2000, gene.names = NULL)
		N.gs <- GSDB$N.gs 
		gs.names[start:(start + N.gs - 1)] <- GSDB$gs.names
		gs.descs[start:(start + N.gs - 1)] <- GSDB$gs.desc
		size.G[start:(start + N.gs - 1)] <- GSDB$size.G
		gs[start:(start + N.gs - 1), 1:max(GSDB$size.G)] <- GSDB$gs[1:N.gs, 1:max(GSDB$size.G)]
		start <- start + N.gs
	}
	N.gs <- max.N
	
	if(save.gene.set.rankings){
		gene.set.rankings = vector("list", length=length(gene.set.selection))
		names(gene.set.rankings) = gene.set.selection
	}
	
	# Select desired gene sets
	
	if (gene.set.selection[1] != "ALL") {
		locs <- match(gene.set.selection, gs.names)
		if(save.gene.set.rankings){
#			browser()
			uncombined.gs.names = c(paste(gene.set.selection[which(is.na(locs))], "_UP", sep=""), 
					paste(gene.set.selection[which(is.na(locs))], "_DN", sep=""))
			locs = na.omit(c(locs, match(uncombined.gs.names, gs.names)))
			gs.names = na.omit(gs.names)
#			N.gs = length(locs)
		} #else 
		N.gs <- sum(!is.na(locs))
		if( N.gs > 1){ gs <- gs[locs,]
		} else{ gs <- t(as.matrix(gs[locs,])) }  # Force vector to matrix if only one gene set specified
		gs.names <- gs.names[locs]
		gs.descs <- gs.descs[locs]
		size.G <- size.G[locs]
	}
#	browser()
	# Loop over gene sets
	
	score.matrix <- matrix(0, nrow=N.gs, ncol=Ns)
	for (gs.i in 1:N.gs) {
		#browser()
		if(save.gene.set.rankings && is.na(gs.names[gs.i])){
			next    ## This geneset is a "combined" version - need to look up both "UP" and "DN" parts
		}
		gene.set <- gs[gs.i, 1:size.G[gs.i]]
		gene.overlap <- intersect(gene.set, gene.names)
		print(paste(gs.i, "gene set:", gs.names[gs.i], " overlap=", length(gene.overlap)))
		if (length(gene.overlap) == 0) { 
			score.matrix[gs.i, ] <- runif(Ns, min=1E-06, max=1.1E-06)
			next
		} else {
#			browser()
			gene.set.locs <- match(gene.overlap, gene.set)
			gene.names.locs <- match(gene.overlap, gene.names)
			msig <- m[gene.names.locs,]
			msig.names <- gene.names[gene.names.locs]
			if (output.score.type == "ES") {
				OPAM <- OPAM.Projection.3(data.array = m, gene.names = gene.names, n.cols = Ns, 
						n.rows = Ng, weight = weight, statistic = statistic,
						gene.set = gene.overlap, nperm = 1, correl.type = correl.type,
						save.gene.set.rankings = save.gene.set.rankings)
				score.matrix[gs.i,] <- OPAM$ES.vector
			} else if (output.score.type == "NES") {
				OPAM <- OPAM.Projection.3(data.array = m, gene.names = gene.names, n.cols = Ns, 
						n.rows = Ng, weight = weight, statistic = statistic,
						gene.set = gene.overlap, nperm = nperm, correl.type = correl.type,
						save.gene.set.rankings = save.gene.set.rankings)
				score.matrix[gs.i,] <- OPAM$NES.vector
			}
			################## --- make gene set rankings plot ---			
			if(save.gene.set.rankings){
#				browser()
				sample.order = match(sample.names.for.ordering, colnames(OPAM$gene.ranks.matrix))
				OPAM$gene.ranks.matrix = OPAM$gene.ranks.matrix[,sample.order]
				OPAM$expression.matrix = OPAM$expression.matrix[,sample.order]
				#browser()
#				expression.with.gene.ranks.heatmap(OPAM$expression.matrix, 
#						OPAM$gene.ranks.matrix,
#						gene.set.name = gs.names[gs.i],
#						graphic.dir = rankings.data.dir, #paste(rankings.data.dir, gs.names[gs.i], ".png", sep=""),
#						graphic.filetype = "png",
#						normalize.function = "zero.one",
#						target = target)
#				V0 = apply(OPAM$expression.matrix, 2, function(x) REVEALER.normalize.zero.one(rev(x))); 
#				V1 = length(mycol)*V0; V2=ifelse(V1==0, 1,V1)
#				quartz(height=8.5,width=11)
#				pdf(height=8.5,width=11)
#				image(1:Ns, 1:Ng, t(V2), col=mycol, axes=FALSE, main=gs.names[gs.i], xlab="", ylab="")
#				# Make lines to delineate samples
#				for(x.ind in 1:(Ns-1)){ 
#					lines(c(x.ind+0.5, x.ind+0.5),c(1,Ng), col="white", lwd=0.1)
#				}
#				# Make lines at the gene set ranks of each
#				for(s.ind in 1:Ns){ 
#					for(rank.ind in Ng+1-OPAM$gene.ranks.matrix[,s.ind]){ 
#						lines(c(s.ind-0.49, s.ind+0.49), c(rank.ind, rank.ind), lwd=1 )
#					}
#				}
#				size.col.char <- 30/(Ns + 15)
#				axis(1, at=1:Ns, labels=sample.names, tick=FALSE, las = 3, 
#						cex.axis=size.col.char, font.axis=2, line=-1)
				
				
#				expression.with.gene.ranks.heatmap(OPAM$expression.matrix, 
#						OPAM$gene.ranks.matrix, gs.names[gs.i], 
#						normalize.function = "z.score", pdf.filename)
				
				write.table(OPAM$gene.ranks.matrix, file=paste(#input.ds, ".", 
								rankings.data.dir,
								gs.names[gs.i], ".gene.ranks.txt", sep=""), quote=FALSE, sep="\t")
				write.table(OPAM$expression.matrix, file=paste(#input.ds, ".",
								rankings.data.dir,
								gs.names[gs.i], ".expression.txt", sep=""), quote=FALSE, sep="\t")
			}
			################ --- end make gene set rankings plot ---			
		}
	}
	
	initial.up.entries <- 0
	final.up.entries <- 0
	initial.dn.entries <- 0
	final.dn.entries <- 0
	combined.entries <- 0
	other.entries <- 0
	
	if (combine.mode == "combine.off") {
		score.matrix.2 <- score.matrix
		gs.names.2 <- gs.names
		gs.descs.2 <- gs.descs
	} else if ((combine.mode == "combine.replace") || (combine.mode == "combine.add")) {
		score.matrix.2 <- NULL
		gs.names.2 <- NULL
		gs.descs.2 <- NULL
		k <- 1
		for (i in 1:N.gs) {
			if(is.na(gs.names[i])) next
			temp <- strsplit(gs.names[i], split="_") 
			body <- paste(temp[[1]][seq(1, length(temp[[1]]) -1)], collapse="_")
			suffix <- tail(temp[[1]], 1)
			print(paste("i:", i, "gene set:", gs.names[i], "body:", body, "suffix:", suffix))
			if (suffix == "UP") {  # This is an "UP" gene set
				initial.up.entries <- initial.up.entries + 1
				target <- paste(body, "DN", sep="_")
				loc <- match(target, gs.names)            
				if (!is.na(loc)) { # found corresponding "DN" gene set: create combined entry
					score <- score.matrix[i,] - score.matrix[loc,]
					score.matrix.2 <- rbind(score.matrix.2, score)
					gs.names.2 <- c(gs.names.2, body)
					gs.descs.2 <- c(gs.descs.2, paste(gs.descs[i], "combined UP & DN"))
					combined.entries <- combined.entries + 1
					if (combine.mode == "combine.add") {  # also add the "UP entry
						score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
						gs.names.2 <- c(gs.names.2, gs.names[i])
						gs.descs.2 <- c(gs.descs.2, gs.descs[i])
						final.up.entries <- final.up.entries + 1
					}
				} else { # did not find corresponding "DN" gene set: create "UP" entry
					score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
					gs.names.2 <- c(gs.names.2, gs.names[i])
					gs.descs.2 <- c(gs.descs.2, gs.descs[i])
					final.up.entries <- final.up.entries + 1
				}
			} else if (suffix == "DN") { # This is a "DN" gene set
				initial.dn.entries <- initial.dn.entries + 1
				target <- paste(body, "UP", sep="_")
				loc <- match(target, gs.names)            
				if (is.na(loc)) { # did not find corresponding "UP" gene set: create "DN" entry
					score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
					gs.names.2 <- c(gs.names.2, gs.names[i])
					gs.descs.2 <- c(gs.descs.2, gs.descs[i])
					final.dn.entries <- final.dn.entries + 1
				} else { # it found corresponding "UP" gene set
					if (combine.mode == "combine.add") { # create "DN" entry
						score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
						gs.names.2 <- c(gs.names.2, gs.names[i])
						gs.descs.2 <- c(gs.descs.2, gs.descs[i])
						final.dn.entries <- final.dn.entries + 1
					}
				}
			} else { # This is neither "UP nor "DN" gene set: create individual entry
				score.matrix.2 <- rbind(score.matrix.2, score.matrix[i,])
				gs.names.2 <- c(gs.names.2, gs.names[i])
				gs.descs.2 <- c(gs.descs.2, gs.descs[i])
				other.entries <- other.entries + 1
			}
			#print(paste("gs.descs.2", gs.descs.2))
		} # end for loop over gene sets
		print(paste("initial.up.entries:", initial.up.entries))
		print(paste("final.up.entries:", final.up.entries))
		print(paste("initial.dn.entries:", initial.dn.entries))
		print(paste("final.dn.entries:", final.dn.entries))
		print(paste("other.entries:", other.entries))
		print(paste("combined.entries:", combined.entries))
		print(paste("total entries:", length(score.matrix.2[,1])))
	}            
	#browser()
	V.GCT <- data.frame(score.matrix.2)
	names(V.GCT) <- sample.names
	row.names(V.GCT) <- make.unique(gs.names.2)
	write.gct(gct.data.frame = V.GCT, descs = make.unique(gs.descs.2), filename = output.ds)  
	
} # end of OPAM.project.dataset.2






OPAM.match.projection.to.phenotypes <-  function(
		input.ds,
		input.cls,
		results.dir,
		normalize.score = T,
		normalization.type = "zero.one",  # "zero.one", "z.score" or "r.z.score"
		markers.num=5,
		user.colors = NA,
		markers.metric = "ROC",   # "ROC" or "T.TEST"
		markers.file = NULL,
		sort.phenotypes = T,
		sort.decreasing = T,    # T = decreasing, F = increasing
		sort.expression = T,
		sort.decreasing.genes = T,
		legend = T,
		char.res = 1,
		only.up = F,
		cmap.type = 3,
		show.desc = T,
		row.norm = T)
{
	
	library(gtools)
	library(verification)
	library(ROCR)
	library(MASS)
	library(RColorBrewer)
	library(heatmap.plus)
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	model.names <- dataset$row.names
	model.descs <- dataset$descs
	Ns <- length(m[1,])
	
	for (i in 1:length(m[,1])) {
		if (sd(m[i,]) == 0) {
			val <- m[i, 1]
			m[i,] <- m[i,] + runif(n=Ns, min= val - 0.001, max=val + 0.001)  # add small noise to flat profiles
		}
	}
	dim(m)
	sample.names <- dataset$names
	
	target.ds$nRow <- length(m[,1])
	temp <- strsplit(input.ds, split="/") # Extract test file name
	s <- length(temp[[1]])
	test.file.name <- temp[[1]][s]
	temp <- strsplit(test.file.name, split=".gct")
	test.file.prefix <-  temp[[1]][1]
#   char.res <-  0.013 * target.ds$nRow + 0.65
	
	# normalize scores
	
	if (normalize.score == T) {
		if (normalization.type == "zero.one") {
			for (i in 1:target.ds$nRow) {
				m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
			}
		} else if (normalization.type == "z.score") {
			for (i in 1:target.ds$nRow) {
				m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
			}
		} else if (normalization.type == "r.z.score") {
			for (i in 1:target.ds$nRow) {
				m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
			}
		}         
	}
	
	CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
	cls.labels <- CLS$class.v
	cls.phen <- CLS$phen
	cls.list <- CLS$class.list 
	
	if (is.vector(cls.labels)) {
		n.phen <- 1
	} else {
		n.phen <- length(cls.labels[,1])
	}
	if (!is.na(user.colors)) {
		c.test <- user.colors
	} else {
		if (!is.null(CLS$col.phen)) {
			c.test <- CLS$col.phen
		} else {
			c.test <- c(brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"))
		}
	}
	
	
	if (!is.null(CLS$phen.names)) {
		phen.names <- CLS$phen.names
#      if (is.vector(cls.list)) {
#         cls.phen <- paste(phen.names, cls.phen, collapse="_")
#      } else {
#         for (i in 1:length(cls.phen)) {
#            for (j in 1:length(cls.phen[[i]])) {
#               cls.phen[[i]][j] <- paste(phen.names[i], cls.phen[[i]][j], collapse="_")
#            }
#         }
#      }
	} else {
		phen.names <- "NA"
	}
	
	cls.phen.index <- unlist(cls.phen)
	cls.phen.colors <- c.test[1:length(cls.phen.index)]
	
	n.classes <- vector(length=n.phen, mode="numeric")
	if (n.phen == 1) {
		max.classes <- length(cls.phen)
		n.classes[1] <- max.classes
	} else {
		max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
		for (i in 1:n.phen) {
			n.classes[i] <- length(cls.phen[[i]])
		}
	}
	
	x <- rbind(sample.names, cls.list, cls.labels)
	print("before loop")
	print(x)
	print(cls.phen)
	print(phen.names)
	
	filename <- paste(results.dir, test.file.prefix, ".PHEN.MARKERS.", markers.metric, ".pdf", sep="")
	pdf(file=filename, height = 10, width = 10)
	
	# Loop over phenotypes
	
	for (k.phen in 1:n.phen) {
		
		
		if (is.vector(cls.labels)) {
			k.phen.labels <- cls.labels
			k.phen.list <- cls.list
		} else {
			k.phen.labels <- as.vector(cls.labels[k.phen,])
			k.phen.list <- as.vector(cls.list[k.phen,])
		}
		
		# Sort according to current phenotype
		
		if(sort.expression == T) {
			phen.index <- order(k.phen.labels, decreasing=sort.decreasing)
		} else {
			phen.index <- seq(1, length(k.phen.labels))
		}
		if (is.vector(cls.labels)) {
			cls.labels2 <- cls.labels[phen.index]
			cls.list2 <- cls.list[phen.index]
		} else {
			cls.labels2 <- cls.labels[, phen.index]
			cls.list2 <- cls.list[, phen.index]
		}
		k.phen.labels <- k.phen.labels[phen.index]
		k.phen.list <- k.phen.list[phen.index]
		sample.names2 <- sample.names[phen.index]
		m2 <- m[, phen.index]
		
		x <- rbind(sample.names2, cls.list2, cls.labels2)
		print(paste("inside loop phen=", k.phen))
		print(x)
		print(cls.phen)
		print(phen.names)
		
		# Markers for each class
		
		if (is.vector(cls.labels2)) {
			classes <- unique(cls.list2)
		} else {
			classes <- unique(cls.list2[k.phen, ])
		}
		if (length(classes) > 2) {
			k.only.up <- T
		} else {
			k.only.up <- only.up
		}
		
		if(length(classes) == 2) classes <- classes[1]
		markers <- NULL
		markers.descs <- NULL
		metric.list <- NULL
		p.val.list <- NULL
		k.class <- NULL
		for (k in classes) {
			if (is.vector(cls.labels2)) {
				bin.class <- ifelse(cls.list2 == k, 0, 1)
			} else {
				bin.class <- ifelse(cls.list2[k.phen, ] == k, 0, 1)
			}
			if (markers.metric == "T.TEST") {
				metric <- vector(length=target.ds$nRow, mode="numeric")
				p.val <- vector(length=target.ds$nRow, mode="numeric")
				for (i in 1:target.ds$nRow) {
					temp <- split(m2[i, ], bin.class)
					x <- temp[[1]]
					y <- temp[[2]]
					metric[i] <- signif(t.test(x=x, y=y)$statistic, digits=3)
					p.val[i] <- signif(t.test(x=x, y=y)$p.value, digits=3)
				}
			} else if (markers.metric == "ROC") {
				bin.class <- ifelse(bin.class == 1, 0, 1)
				metric <- vector(length=target.ds$nRow, mode="numeric")
				p.val <- vector(length=target.ds$nRow, mode="numeric")
				for (i in 1:target.ds$nRow) {
					m.score <- m2[i,]
					m.score.norm <- (m.score - min(m.score))/(max(m.score) - min(m.score))
					if (length(table(bin.class)) > 1) {
						perf.auc <- roc.area(bin.class, m.score.norm)
						metric[i] <- signif(perf.auc$A, digits=3)
						p.val[i] <- signif(perf.auc$p.value, digits=3)
					} else {
						metric[i] <- 1
						p.val[i] <- 1
					}
				}
			} else if (markers.metric == "MEAN.DIFF") {
				bin.class <- ifelse(bin.class == 1, 0, 1)
				metric <- vector(length=target.ds$nRow, mode="numeric")
				p.val <- vector(length=target.ds$nRow, mode="numeric")
				for (i in 1:target.ds$nRow) {
					temp <- split(m2[i, ], bin.class)
					x <- temp[[1]]
					y <- temp[[2]]
					metric[i] <- signif(mean(x) - mean(y), digits=3)
					p.val[i] <- signif(t.test(x=x, y=y)$p.value, digits=3)
				}
			}
			
			if (is.na(sort.decreasing.genes)) {
				metric.order <- seq(1, length(metric))
			} else {
				metric.order <- order(metric, decreasing=sort.decreasing.genes)
			}
			if (only.up == TRUE) {
				k.markers.num <- ifelse(markers.num > target.ds$nRow, target.ds$nRow, markers.num)
				
#            if (length(classes) == 2) {
#               k.markers.num <- ifelse(markers.num > target.ds$nRow, target.ds$nRow, markers.num)
#            } else {
#               k.markers.num <- ifelse(length(classes)*markers.num > target.ds$nRow, 
#                                               floor(target.ds$nRow/length(classes)), markers.num)
#            }
				markers <- c(markers, model.names[metric.order][1:k.markers.num])
				markers.descs <- c(markers.descs, model.descs[metric.order][1:k.markers.num])
				metric.list <- c(metric.list, metric[metric.order][1:k.markers.num])
				p.val.list <- c(p.val.list, p.val[metric.order][1:k.markers.num])
				k.class <- c(k.class, rep(k, k.markers.num))
			} else {
				k.markers.num <- ifelse(length(classes)*markers.num > target.ds$nRow, floor(target.ds$nRow/length(classes)), 
						markers.num)
				markers <- c(markers, model.names[metric.order][1:k.markers.num],
						model.names[metric.order][(length(model.names) - k.markers.num +1):length(model.names)])
				markers.descs <- c(markers.descs, model.descs[metric.order][1:k.markers.num],
						model.descs[metric.order][(length(model.names) - k.markers.num + 1):length(model.names)])
				metric.list <- c(metric.list, metric[metric.order][1:k.markers.num],
						metric[metric.order][(length(model.names) - k.markers.num + 1):length(model.names)])
				p.val.list <- c(p.val.list, p.val[metric.order][1:k.markers.num],
						p.val[metric.order][(length(model.names) - k.markers.num + 1):length(model.names)])
				k.class <- c(k.class, rep(k, k.markers.num), rep(paste("not", k), k.markers.num))
			}
		}
		
		V3 <- m2[markers,]
		print(V3)
		print(markers)
		
		if (show.desc == T) {
			model.descs2 <- paste(metric.list, p.val.list, k.class, markers.descs)
		} else {
			model.descs2 <- paste(metric.list, p.val.list)
		}
		height <- ifelse(length(markers) + n.phen >= 9, 10, (length(markers) + n.phen)*0.44 + 5)
#      char.res <-  0.0085 * length(markers) + 0.65
		
		
		# Sort markers inside each phenotype class
		
		if(sort.expression == T) {
			for (j in unique(k.phen.labels)) {
				V4 <- V3[ , k.phen.labels == j]
				sn <- sample.names2[k.phen.labels == j]
				if (is.vector(cls.labels)) {
					clab <- cls.labels2[k.phen.labels == j]
					clis <- cls.list2[k.phen.labels == j]
				} else {
					clab <- cls.labels2[, k.phen.labels == j]
					clis <- cls.list2[, k.phen.labels == j]
				}
				l.phen <- sum(k.phen.labels == j)
				if (l.phen > 1) {
					dist.matrix <- dist(t(V4))
					HC <- hclust(dist.matrix, method="complete")
					HC.order <- HC$order
					V4 <- V4[ , HC.order]
					sn <- sn[HC.order]
					if (is.vector(cls.labels2)) {
						clab <- clab[HC.order]
						clis <- clis[HC.order]
					} else {
						clab <- clab[, HC.order]
						clis <- clis[, HC.order]
					}
				}
				V3[ , k.phen.labels == j] <- V4
				sample.names2[k.phen.labels == j] <- sn
				if (is.vector(cls.labels2)) {
					cls.labels2[k.phen.labels == j] <- clab
					cls.list2[k.phen.labels == j] <- clis
				} else {
					cls.labels2[, k.phen.labels == j] <- clab
					cls.list2[, k.phen.labels == j] <- clis
				}
			}
		}
		x <- rbind(sample.names2, cls.list2, cls.labels2)
		print(paste("inside loop after in-class sort phen=", k.phen))
		print(x)
		print(cls.phen)
		print(phen.names)
		
		# Recompute cls.phen and cls.labels2 as order may have changed
		
		cls.phen2 <- list(NULL)
		if (is.vector(cls.labels2)) {
			classes <- unique(cls.list2)
			cls.phen2 <- classes
			cls.labels2 <- match(cls.list2, cls.phen2)
		} else {
			for (kk in 1:length(cls.list2[, 1])) {
				classes <- unique(cls.list2[kk,])
				cls.phen2[[kk]] <- classes
				cls.labels2[kk,] <- match(cls.list2[kk,], cls.phen2[[kk]])
			}
		}
		
		x <- rbind(sample.names2, cls.list2, cls.labels2)
		print(paste("inside loop after cls.phen renorm phen=", k.phen))
		print(cls.phen2)
		print(phen.names)
		
		
		library(gmodels)
		if (!is.vector(cls.labels2)) {
			if (sort.phenotypes == T) {
				phen.score <- vector(length=n.phen, mode="numeric")
				for (k.lab in 1:n.phen) {
					tab <- table(as.vector(cls.list2[k.lab,]), k.phen.list)
					print(tab)
#              phen.score[k.lab] <- 1 - chisq.test(tab)$p.value
#              phen.score[k.lab] <- 1 - fisher.test(tab)$p.value
					if ((length(tab[,1]) > 1) && (length(tab[1,]) > 1)) { 
						CT <- CrossTable(tab, chisq=T)
						phen.score[k.lab] <- CT$chisq$p.value
						print(phen.score[k.lab])
					} else {
						phen.score[k.lab] <- 0.50
						print(phen.score[k.lab])
					}
				}
				phen.order <- order(phen.score, decreasing= T)
				print(phen.order)
				cls.labels2 <- cls.labels2[phen.order,]
				cls.phen2 <- cls.phen2[phen.order]
				phen.names2 <- phen.names[phen.order]
				main.string <- paste(test.file.prefix, " - ", phen.names2[n.phen], markers.metric, " order")
			} else {
				phen.names2 <- phen.names
				main.string <- paste(test.file.prefix, " - ", phen.names2[k.phen], markers.metric, " order")
			}
		} else {
			phen.names2 <- phen.names[1]
			main.string <- paste(test.file.prefix, " - ", phen.names2, markers.metric, " order")
		}
		
#     windows(width=15, height=height)
		
		
		x <- rbind(sample.names2, cls.list2, cls.labels2)
		print(paste("inside loop after phen sort before figure phen=", k.phen))
		print(x)
		print(cls.phen2)
		print(phen.names2)
		
		phen.list <- unlist(cls.phen2)
		colors.list <- cls.phen.colors[match(phen.list, cls.phen.index)]
		
		print(rbind(phen.list, colors.list))
		
		if (show.desc == T) {
			markers <- paste(markers, seq(1, length(markers)), sep="_")
		}
		
		MSIG.HeatMapPlot.7(V = V3, row.names = markers,
				row.names2 = model.descs2, col.labels = cls.labels2, 
				col.classes = cls.phen2, phen.cmap = colors.list, phen.names = phen.names2,
				col.names = sample.names2, main = main.string, xlab="  ", ylab="  ", 
				row.norm = row.norm,  
				cmap.type = cmap.type, char.rescale = char.res,  legend=legend)
		
		V3 <- data.frame(V3)
		colnames(V3) <- sample.names2
		row.names(V3) <- markers
		
		if (!is.null(markers.file)) {
			write.gct(gct.data.frame = V3, descs = model.descs2, filename = markers.file)  
		}
		
	} # end loop over phenotypes
	
	dev.off()
	
}



OPAM.sort.projection.by.score.7 <- function(
#		input.ds,
		signatures = "NA",
		input.all.pathways.ds,
		input.cls,
		tissue = "NA",
		results.dir = NULL,
		normalize.score = T,
		normalization.type = "zero.one",
		model = "NA",
		target.phen = NA,
		target.class = NA,
		user.colors = NA,
#		decreasing.order = T,
		output.dataset = NA,
		char.rescale = 1,
		cmap.type = 3,
		row.norm = T,
		u.gene.names.known = "NA",
		add.amp.del = FALSE,
		#n.random.signatures = 10,
		show.multiple.tissues = FALSE,
		add.chrom.locs = FALSE,
		file.suffix = "",
		add.mut = FALSE,
		n.iter = 5,
		pdf.height = 11,
		pdf.width = 17,
		#do.mRMR = FALSE,
		#skip.step1 = FALSE,
		skip.aberration.refinement = FALSE,
		de.novo.search = FALSE,
		align.only.no.discovery = FALSE,
		continuous.data = FALSE,
		match.to.blue = FALSE,
		highlight.tissue.name = NULL,
		ifLegend = FALSE,
		ifScoreBarPlot = FALSE,
		make.inbetween.heatmaps = FALSE,
		exclude.samples = NA
)
# Calls MSIG.HeatMapPlot.9 and makes a plot sorted by the highest-scoring
# signatures and abnormalities (gene mutations or copy number alterations)
# i.e. doesn't require a "model" to score by as OPAM.sort.projection.by.score.2 does.
# However, it *will* use "model" if it cannot calculate p-values on the gene signatures, which
# happens when every cell line has a genomic aberration.
#
# Runs 3 passes on the data:
# 1st pass: looks at the genes and copy number alterations specified by u.gene.names.known
# 2nd pass: looks at only the top abnormalities (using a p-value cutoff) from the 1st pass, and adjusts 
# the PATHWAY.MUT vector accordingly (only according to the genes, not by copy number data)
# 3rd pass: Takes the winning signature from the 2nd pass and then looks all the genes available
#
# Very similar to OPAM.sort.projection.by.score.4, however this version uses mutual.inf instead of
# roc.area to calculate mutual information scores and p-values for PATHWAY.MUT, the vector of total genomic aberrations
# in all samples
#
# Differs from OPAM.sort.projection.by.score.6 by requiring the gct file of expression in 
# all pathways by the input tissue ("input.all.pathways.ds")
{
	
	library(gtools)
	library(verification)
	library(ROCR)
	library(MASS)
	library(RColorBrewer)
	library(heatmap.plus)
	
	
	decreasing.order = ifelse(match.to.blue, FALSE, TRUE)
	if( add.chrom.locs ) file.suffix = paste(file.suffix, "_with.chrom.locs", sep="")
	
	dataset.all <- MSIG.Gct2Frame( filename = input.all.pathways.ds)
	m.all <- data.matrix(dataset.all$ds)
	model.names.all <- dataset.all$row.names
	Ns = length(m.all[1,])
	sample.names = dataset.all$names
	if(!is.na(exclude.samples)){
		#browser()
		exclude.samples.ind = which(sample.names %in% exclude.samples)
		sample.names = sample.names[-exclude.samples.ind]
		m.all = m.all[, -exclude.samples.ind]
		Ns = length(m.all[1,])
	}
	
	if( show.multiple.tissues ){
		tissue.type <- vector(length=Ns, mode="character")
#		temp = strsplit(sample.names, split="_")
		for (k in 1:Ns) {
			temp <- strsplit(sample.names[k], split="_") 
			tissue.type[k] <- paste(temp[[1]][2:length(temp[[1]])], collapse="_")
		}
		tissue.names = unique(tissue.type)
		tissue.labels = match(tissue.type, tissue.names)
	} else{
		tissue.names = tissue
		tissue.labels = rep(1, Ns)
	}
	
	if( is.na(signatures[1]) ){
		stop("Must provide a vector of signature names to evaluate, or specify 'ALL'")
	}
	
	## Remove "Commented out" signatures (with # at beginning of name)
	if( length(grep("^#", signatures)) > 0){
		signatures = signatures[-grep("^#", signatures)]
	}
	if( signatures[1] == "ALL"){
		model.names = model.names.all
		m = m.all
		model.descs = dataset.all$descs
	} else{
		model.names = signatures
		model.ind = match(signatures, model.names.all)
		m = m.all[model.ind,]
		model.descs = dataset.all$descs[model.ind]
		if( length(model.ind) == 1 ){
			m = t(as.matrix(m))
			rownames(m) = model.names
		}
		rm(list=c("m.all", "dataset.all"))
#	browser()
	}
	
	## Remove NA rows
	na.row.index = unique(as.vector(apply(m, MARGIN=2, FUN=function(x) which(is.na(x)))))
	if(length(na.row.index) > 0 ) { m = m[-na.row.index,] } 
	target.ds$nRow <- length(m[,1])
	temp <- strsplit(input.all.pathways.ds, split="/") # Extract test file name
	s <- length(temp[[1]])
	test.file.name <- temp[[1]][s]
	temp <- strsplit(test.file.name, split=".gct")
	test.file.prefix <-  temp[[1]][1]
	char.res <-  0.013 * target.ds$nRow + 0.65
	
	#browser()
	if( !is.null(results.dir) ){
		#browser()
		if( length(grep("./$", results.dir)) == 0 ){
			results.dir = paste(results.dir, "/", sep="")
		}
	}
	
	# normalize scores
	
	if (normalize.score == T) {
		if (normalization.type == "zero.one") {
			#browser()
			for (i in 1:target.ds$nRow) {
				m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
			}
		} else if (normalization.type == "z.score") {
			for (i in 1:target.ds$nRow) {
				m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
			}
		} else if (normalization.type == "r.z.score") {
			for (i in 1:target.ds$nRow) {
				m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
			}
		}         
	}
	
#	if( substr(input.cls, nchar(input.cls)-3, nchar(input.cls)) == ".txt"){
#		## Input a text file with
#		
#	}
	
	CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
	cls.labels <- CLS$class.v
	cls.phen <- CLS$phen
	cls.list <- gsub("\\s+", "", CLS$class.list, perl=TRUE)  # get rid of trailing whitespace   
	
	if(!is.na(exclude.samples)){
		cls.labels = cls.labels[, -exclude.samples.ind]
		cls.list = cls.list[, -exclude.samples.ind]
	}
	
#	browser()
	
	if( continuous.data ){
		if( is.vector(cls.list) ){  # If continuous.data=TRUE and cls.list is a vector, then the 
			# one row in cls.list must be continuous and do not need to
			# iterate and test every row as when cls.list is a matrix
			cls.list = as.double(cls.list)
		} else{                     # cls.list is a matrix and therefore need to find the rows which
			# are continuous. If as.double returns NA, then the elements are
			# characters and are left alone. Otherwise if all the elements are
			# 0's and 1's then they are replaced with "WT" and "MUT"
			# respectively. If the elements are not just 0's and 1's then they
			# are assumed to be continuous numbers and are forced to doubles.
			apply(cls.list, MARGIN=1, FUN=function(x) {
						if( !is.na(as.double(x[1,1])) ) {
							if( sum(as.double(x[1,]) %in% 0:1) == length(x[1,]) ){
								return( ifelse(x=="1", "MUT", "WT"))
							} else { return (as.double(x)) } } else return(x) } )
		}
	} else if( "0" %in% cls.list && "1" %in% cls.list){
		cls.list[which(cls.list == "0")] = "WT"
		cls.list[which(cls.list == "1")] = "MUT"
	}
	
#	browser()
	
	if (is.vector(cls.labels)) {
		n.phen <- 1
	} else {
		n.phen <- length(cls.labels[,1])
	}
	if (!is.na(user.colors[1])) {
		c.test <- user.colors
	} else {
		if (!is.null(CLS$col.phen)) {
			c.test <- CLS$col.phen
		} else {
			c.test <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"))
		}
	}
	
	
	if (!is.null(CLS$phen.names)) {
		phen.names <- CLS$phen.names
	} else {
		phen.names <- "NA"
	}
	
	
	cls.phen.index <- unlist(cls.phen)
	cls.phen.colors <- c.test[1:length(cls.phen.index)]
#	print("cls.phen.colors:")
#	print(cls.phen.colors)
	
	n.classes <- vector(length=n.phen, mode="numeric")
	if (n.phen == 1) {
		max.classes <- length(cls.phen)
		n.classes[1] <- max.classes
	} else {
		max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
		n.classes = unlist(lapply(cls.phen, length))
#		for (i in 1:n.phen) {
#			n.classes[i] <- length(cls.phen[[i]])
#		}
	}
	
	pdf.options(height=pdf.height, width=pdf.width, colormodel="rgb", bg="transparent")
	
	## Get rid of trailing whitespace
#	browser()
	cls.list = rbind(rep("WT", length=ifelse(is.vector(cls.list), length(cls.list), length(cls.list[1,]))), 
			apply(cls.list, MARGIN=2, FUN=function(x) gsub("\\s+$", "", x, perl=TRUE)))
	cls.labels = rbind(rep(1, length=ifelse(is.vector(cls.list), length(cls.list), length(cls.list[1,]))), 
			cls.labels)
#			gsub("\\s+$", "", cls.labels, #unlist(cls.phen2.pass1), 
#					perl=TRUE))
	
#	browser()
	
	## Remove "Commented out" gene names (with # at beginning of name)
	if( length(grep("^#", u.gene.names.known)) > 0){
		u.gene.names.known = u.gene.names.known[-grep("^#", u.gene.names.known)]
	}
	if(add.mut){
		u.gene.names.known.mut = paste(u.gene.names.known, "_MUT", sep="")
		
	} else{ u.gene.names.known.mut = NULL }
	if ( add.amp.del ){
		u.gene.names.known.amp.del = c( u.gene.names.known, paste(u.gene.names.known, "_AMP", sep=""), 
				paste(u.gene.names.known, "_DEL", sep="") )
	} else{u.gene.names.known.amp.del = NULL }
	
	if( add.amp.del || add.mut ){
		phen.pass1 = c(u.gene.names.known.mut, u.gene.names.known.amp.del)
	} else{ phen.pass1 = c( u.gene.names.known ) }
	
	
	## Find chromosomal locations of genes specified
	## See "if( find.chromosomal.locations)" for more transparent code
	if( add.chrom.locs ){
		library(org.Hs.eg.db)
		
#			browser()
		phen.pass1.split = strsplit(phen.pass1, split="_")
		phen.pass1.noampdel = unlist( lapply(phen.pass1.split, function(x) x[1]))
		phen.pass1.egIDs = mget(phen.pass1.noampdel, org.Hs.egALIAS2EG, ifnotfound=NA)
		phen.pass1.no.chrom.loc = which(lapply(lapply(phen.pass1.egIDs, is.na), sum) > 0)
		if( length(phen.pass1.no.chrom.loc) == 0){
			phen.pass1.egIDs.w.locs = phen.pass1.egIDs
		} else{
			phen.pass1.egIDs.w.locs = phen.pass1.egIDs[-phen.pass1.no.chrom.loc]
		}
		phen.pass1.locs.list = lapply(phen.pass1.egIDs.w.locs, mget, org.Hs.egMAP)
		phen.pass1.locs = vector(mode="character", length=length(phen.pass1))
		if( length(phen.pass1.no.chrom.loc) == 0){
			phen.pass1.locs = unlist(lapply(phen.pass1.locs.list, function(x) paste(unlist(x), collapse="_")))
		} else{
			phen.pass1.locs[phen.pass1.no.chrom.loc] = "NA"
			phen.pass1.locs[-phen.pass1.no.chrom.loc] = (
						unlist(lapply(phen.pass1.locs.list, function(x) paste(unlist(x), collapse="_"))) )
		}
		phen.pass1.w.locs = paste(phen.pass1.noampdel, ".", phen.pass1.locs, sep="")
		phen.pass1.ampdel.suffix = unlist(lapply(phen.pass1.split, function(x) x[2]))
		phen.pass1 = paste(phen.pass1.w.locs, "_", phen.pass1.ampdel.suffix, sep="")
		print("Found chromosomal locations for all genes in u.gene.names.known!")
		
		
		phen.names.split = strsplit(phen.names, split="_")
		phen.names.noampdel = unlist( lapply(phen.names.split, function(x) x[1]))
		phen.names.egIDs = mget(phen.names.noampdel, org.Hs.egALIAS2EG, ifnotfound=NA)
		phen.names.no.chrom.loc = which(lapply(lapply(phen.names.egIDs, is.na), sum) > 0)
		if( length(phen.names.no.chrom.loc) == 0){
			phen.names.egIDs.w.locs = phen.names.egIDs
		} else{
			phen.names.egIDs.w.locs = phen.names.egIDs[-phen.names.no.chrom.loc]
		}
		phen.names.locs.list = lapply(phen.names.egIDs.w.locs, FUN=mget, org.Hs.egMAP)
		phen.names.locs = vector(mode="character", length=length(phen.names))
		if( length(phen.names.no.chrom.loc) == 0){
			phen.names.locs = unlist(lapply(phen.names.locs.list, function(x) paste(unlist(x), collapse="_")))
		} else{
			phen.names.locs[phen.names.no.chrom.loc] = "NA"
			phen.names.locs[-phen.names.no.chrom.loc] = (
						unlist(lapply(phen.names.locs.list, function(x) paste(unlist(x), collapse="_"))) )
		}
		phen.names.w.locs = paste(phen.names.noampdel, ".", phen.names.locs, sep="")
		phen.names.ampdel.suffix = unlist(lapply(phen.names.split, function(x) x[2]))
		#phen.names.no.suffix = which(is.na(phen.names.ampdel.suffix))
		#phen.names[phen.names.no.suffix] = phen.names.w.locs[phen.names.no.suffix]
#			phen.names[-phen.names.no.suffix] = paste(phen.names.w.locs[-phen.names.no.suffix], 
#					"_", phen.names.ampdel.suffix[-phen.names.no.suffix], sep="")
		phen.names = paste(phen.names.w.locs, "_", phen.names.ampdel.suffix, sep="")
		print("Found chromosomal locations for all genes in cls file!")
		
		rm(list=c("phen.names.split", "phen.names.noampdel", "phen.names.egIDs", "phen.names.no.chrom.loc", 
						"phen.names.egIDs.w.locs", "phen.names.locs.list", "phen.names.locs", "phen.names.w.locs",
						"phen.names.ampdel.suffix", "phen.pass1.split", "phen.pass1.noampdel", 
						"phen.pass1.egIDs", "phen.pass1.no.chrom.loc", 
						"phen.pass1.egIDs.w.locs", "phen.pass1.locs.list", "phen.pass1.locs", "phen.pass1.w.locs",
						"phen.pass1.ampdel.suffix"))
	}
	
	
	## Was originally immediately after "ind.phen.pass1 = ..." but since now want to find the chromosomal
	## locations of the genes, have to first find the indices of the genes specified at the onset of the
	## program, THEN find all the chromosomal locations
	
	#browser()
	
	#browser()
#		phen.pass1 = c(phen.names[1], phen.pass1)
	if( !de.novo.search ){
		
		#if( !skip.step1 ){
		print("--- Begin Pass 1 ---")
		phen.names = c("SUMMARY", phen.names)
		ind.phen.pass1 = c(1, which( phen.names %in% phen.pass1 ))
		## I tried to do this as a tryCatch but couldn't get it to work, so I kept it simple
		if(length(ind.phen.pass1)==1){
			stop(simpleError(paste("The gene names provided do not match with",
									"the cls file. Some common mistakes:",
									"the genes must either have '_MUT', '_AMP',",
									"or '_DEL' suffixed to them, or", 
									"'add.mut' and/or 'add.amp.del' set to TRUE")))
		}
		phen.pass1 = phen.names[ind.phen.pass1]
		n.phen.pass1 = length(phen.pass1)
		MI.list.pass1 = vector( length=target.ds$nRow, mode="numeric" )
		
		
		
		cls.list.pass1 = cls.list[ind.phen.pass1,]
		if( !continuous.data ){ cls.list.pass1.2 = ifelse(cls.list.pass1 == "WT", 0, 1)
		} else{ cls.list.pass1.2 = cls.list.pass1 }
		cls.labels.pass1 = cls.labels[ind.phen.pass1,]
		
		if (!is.na(target.phen)) {
			if( length( phen.pass1) > 2 ){
#				browser()
				bin.class.pass1 = colSums(cls.list.pass1.2[-1,]) #apply( cls.list.pass1.2[-1,], MARGIN=2, FUN=sum)
			} else{ bin.class.pass1 = cls.list.pass1.2[2,] }
			# Normalize bin.class.pass1
			if( length(unique(bin.class.pass1)) > 1){
				bin.class.pass1 = REVEALER.normalize.zero.one(bin.class.pass1)
				(bin.class.pass1) #( bin.class.pass1 - min(bin.class.pass1))/(max(bin.class.pass1) - min(bin.class.pass1))
			} else if ( length(unique(bin.class.pass1)) == 1){
				bin.class.pass1 = rep(1, length(cls.list[1,]))
			}
			if( !continuous.data ){
				cls.list.pass1[1,] = ifelse(bin.class.pass1 > 0, "MUT", "WT")
			} else{ cls.list.pass1[1,] = bin.class.pass1 }
		} else {
			bin.class.pass1 <- ifelse(cls.list[1,] == cls.list2[1,1], 1, 0)
		}
		#browser()
		
		if(make.inbetween.heatmaps){
			### Make initial heatmap ###
			cls.phen2.pass1 <- NULL
			if (is.vector(cls.labels)) {
				classes <- unique(cls.list.pass1)
				cls.phen2.pass1 <- classes
				cls.labels.pass1 <- match(cls.list.pass1, cls.phen2.pass1)
			} else {
				#browser()
				for (kk in 1:length(cls.list.pass1[, 1])) {
					classes <- unique(cls.list.pass1[kk,])
#            cls.phen2[[kk]] <- classes
					cls.phen2.pass1 <- c(cls.phen2.pass1, classes)
					cls.labels.pass1[kk,] <- match(cls.list.pass1[kk,], classes)
				}
			}
			#browser()
			cls.labels.pass1 = cls.labels.pass1[1:length(phen.pass1),]
			phen.list.pass1 <- unlist(cls.phen2.pass1)
			colors.list = rep( "gray", length(phen.list.pass1))
			colors.list[phen.list.pass1!="WT"] = cls.phen.colors[1]
#		browser()
			filename <- paste(results.dir, test.file.prefix, file.suffix, ".pre.calculations", sep="")
			pdf(file=paste(filename, ".pdf", sep=""), height = pdf.height, width = pdf.width )
#		browser()
#quartz(height = 11, width = 17)
			if( show.multiple.tissues ){
				#browser()
				#quartz(height = 11, width = 17)
				MSIG.HeatMapPlot.10.show.multiple.tissues(V = m, 
						pathway.mut = bin.class.pass1,
						row.names = model.names,
						col.labels = cls.labels.pass1, 
						col.classes = cls.phen2.pass1, 
						phen.cmap = colors.list, 
						phen.names = phen.pass1,
						col.names = sample.names, 
						main = paste(tissue, "- Initial Heatmap ('Step 0')"), 
						xlab="  ", ylab="  ", row.norm = row.norm,  
						cmap.type = cmap.type, char.rescale = char.rescale,  legend=ifLegend,
						tissue.names = tissue.names,
						tissue.labels = tissue.labels,
						highlight.tissue.name = highlight.tissue.name)
			} else{
				#browser()
				MSIG.HeatMapPlot.10(V = m, 
						pathway.mut = bin.class.pass1,
						row.names = model.names,
						col.labels = cls.labels.pass1, 
						col.classes = cls.phen2.pass1, 
						phen.cmap = colors.list, 
						phen.names = phen.pass1,
						col.names = sample.names, 
						main = paste(tissue, "- Initial Heatmap ('Step 0')"), 
						xlab="  ", ylab="  ", row.norm = row.norm,  
						cmap.type = cmap.type, char.rescale = char.rescale, 
						legend=ifLegend,
						ifScoreBarPlot = FALSE)
			}
			dev.off()
		}
		
		model.descs2.pass1 = vector(length = target.ds$nRow, mode="character")
		
		if( length(unique(bin.class.pass1)) > 1 ){
			
			if( target.ds$nRow > 1 ){
				MI.results = mutual.inf.3.v2(bin.class.pass1, m, 
						target.vector.name="SUMMARY", 
						tissue=tissue)
				MI.list.pass1  = MI.results$MI
				model.descs2.pass1 <- sapply(MI.results$MI, FUN=signif, 3)
				m.order.pass1 = order(MI.list.pass1, 
						decreasing=decreasing.order, na.last=TRUE)
				m2.pass1 <- m[m.order.pass1, ]
				s.order.pass1 <- order(m2.pass1[1,], 
						decreasing = decreasing.order )
				m2.pass1 <- m2.pass1[, s.order.pass1]
			} else{ 
				#browser()
				MI.ref = mutual.inf.2(bin.class.pass1, bin.class.pass1)
				MI.list.pass1 = MI.results = 
						mutual.inf.2(bin.class.pass1, m[1,])/MI.ref
				model.descs2.pass1 <- signif(MI.results, digits=3)
				m2.pass1 <- m ; m.order.pass1 = 1
				s.order.pass1 <- order(m2.pass1[1,], 
						decreasing = decreasing.order )
				m2.pass1 <- t(as.matrix(m[, s.order.pass1]))
				rownames(m2.pass1) = model.names
			}
		} else{ 
			MI.list.pass1 = rep(NA, target.ds$nRow)
			FDR.list.pass1 = rep(NA, target.ds$nRow)
			model.descs2.pass1 = rep(" - (FDR = - )", target.ds$nRow)
			if( target.ds$nRow > 1 ){
				loc <- match(model, model.names)
				s.order.pass1 <- order(m[loc,], decreasing = decreasing.order)
				m2.pass1 <- m[, s.order.pass1]
				correl <- cor(t(m2.pass1))[, loc]
				m.order.pass1 <- order(correl, decreasing=T)
				m2.pass1 <- m2.pass1[m.order.pass1, ]
			} else{ 
				m2.pass1 <- t(as.matrix(m[, s.order.pass1]))
				rownames(m2.pass1) = model.names
				m.order.pass1 = 1 }
		}
		
		MI.list.pass1 = MI.list.pass1[m.order.pass1]
		bin.class.pass1 = bin.class.pass1[s.order.pass1]
		
		model.descs2.pass1.all = model.descs2.pass1
		model.descs2.pass1 = model.descs2.pass1[m.order.pass1]
		sample.names2.pass1 <- colnames(m2.pass1)
		model.names.pass1 <- rownames(m2.pass1)
		print(matrix(c(model.names.pass1, model.descs2.pass1), ncol=2), quote=F)
		if (is.vector(cls.labels)) {
			cls.labels2.pass1 <- cls.labels.pass1[s.order.pass1]
			cls.list2.pass1 <- cls.list.pass1[s.order.pass1]
		} else {
			cls.labels2.pass1 <- cls.labels.pass1[, s.order.pass1]
			cls.list2.pass1 <- cls.list.pass1[, s.order.pass1]
		}
		tissue.labels.pass1 = tissue.labels[s.order.pass1]
		sample.names2 <- colnames(m2.pass1)
		winning.model.ind.pass1 = which(model.names.pass1[1] == rownames(m2.pass1))
		
		MI.list.phen.pass1 = vector(mode="numeric", length=n.phen.pass1)
		phen.descs.pass1 = vector(mode="character", length=n.phen.pass1)
		
		if( length(unique(bin.class.pass1)) > 1){
			MI.signif <- signif(MI.list.pass1[1], digits=3)
			MI.list.phen.pass1[1] = MI.list.pass1[1]
			phen.descs.pass1[1] = model.descs2.pass1[1]
		} else{
			MI.signif <- "-"
			MI.list.phen.pass1[1] = NA
		}
		print(paste(format(phen.pass1[1], width=12), "mutual.inf =", MI.signif
				))
		print(proc.time()-t1)
		print(date())
		phen.descs.pass1[1] <- paste(MI.signif,
				sep="")
		
#	browser()
		if( n.phen.pass1 > 2 ){
			bin.gene.matrix = ifelse(cls.list2.pass1[-1,]=="WT", 0, 1)
			MI.results = mutual.inf.3.v2(
					m2.pass1[winning.model.ind.pass1,],
					bin.gene.matrix)
#			browser()
			MI.list.phen.pass1[-1] = MI.results$MI
			phen.descs.pass1[-1] = sapply(MI.results$MI, FUN=signif, 3)
			
			g.order.pass1 = c(1, order(MI.list.phen.pass1[-1], 
							decreasing=decreasing.order, na.last=TRUE)+1)
			MI.list.phen.pass1 = MI.list.phen.pass1[g.order.pass1]
			phen.descs2.pass1 = phen.descs.pass1[g.order.pass1]
			cls.list2.pass1 = cls.list2.pass1[g.order.pass1,]
			phen.names.pass1 = phen.pass1[g.order.pass1]
		} else{
#		bin.gene.matrix = ifelse(cls.list2.pass1[-1,]=="WT", 0, 1)
#		MI.ref = mutual.inf.2(m2.pass1[winning.model.ind.pass1,],
#				m2.pass1[winning.model.ind.pass1,])
#		MI.results = mutual.inf.2(
#				m2.pass1[winning.model.ind.pass1,],
#				bin.gene.matrix)/MI.ref
#		MI.list.phen.pass1[-1] = MI.results
			phen.descs.pass1[-1] = phen.descs.pass1[1]
			
			g.order.pass1 = 1:2
#	MI.list.phen.pass1 = MI.list.phen.pass1[g.order.pass1]
			phen.descs2.pass1 = phen.descs.pass1
#	cls.list2.pass1 = cls.list2.pass1[g.order.pass1,]
			phen.names.pass1 = phen.pass1
		}
		print(matrix(c(phen.names.pass1, phen.descs2.pass1), ncol=2), quote=F)
		print(proc.time()-t1)
		print(date())
		
		# Recompute cls.list2 as some mutations or copy numbers may have been removed
		# Recompute cls.phen and cls.labels2 as order may have changed
		cls.phen2.pass1 <- NULL
		if (is.vector(cls.labels)) {
			classes <- unique(cls.list2.pass1)
			cls.phen2.pass1 <- classes
			cls.labels2.pass1 <- match(cls.list2.pass1, cls.phen2.pass1)
		} else {
			for (kk in 1:length(cls.list2.pass1[, 1])) {
				classes <- unique(cls.list2.pass1[kk,])
#            cls.phen2[[kk]] <- classes
				cls.phen2.pass1 <- c(cls.phen2.pass1, classes)
				cls.labels2.pass1[kk,] <- match(cls.list2.pass1[kk,], classes)
			}
		}
		cls.labels2.pass1 = cls.labels2.pass1[1:n.phen.pass1,]
		
		phen.list.pass1 <- unlist(cls.phen2.pass1)
		colors.list.pass1 = rep( "gray", length(phen.list.pass1))
		colors.list.pass1[phen.list.pass1!="WT"] = cls.phen.colors[1]
#		colors.list.pass1[phen.list.pass1=="DEL"] = cls.phen.colors[3]
#		colors.list.pass1[phen.list.pass1=="AMP"] = cls.phen.colors[4]
#		colors.list.pass1[phen.list.pass1=="ALT"] = cls.phen.colors[5]
		
		
		filename <- paste(results.dir, test.file.prefix, file.suffix, ".",
				phen.names.pass1[2],
				".initial.calculations", sep="")
#    pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 9)
		pdf(file=paste(filename, ".pdf", sep=""), height = pdf.height, width = pdf.width )
#		browser()
#   windows(width=12, height=8)
		if( show.multiple.tissues ){
			MSIG.HeatMapPlot.10.show.multiple.tissues(V = m2.pass1, 
					pathway.mut = bin.class.pass1,
					row.names = model.names.pass1,
					row.names2 = model.descs2.pass1, 
					col.labels = cls.labels2.pass1, 
					col.classes = cls.phen2.pass1, 
					phen.cmap = colors.list.pass1, 
					phen.names = phen.names.pass1,
					phen.names2 = c(" ", phen.descs2.pass1),
					col.names = sample.names2.pass1, 
					main = paste(tissue, "- Step 1: Known KRAS Pathway Abnormalities (MI)"), 
					xlab="  ", ylab="  ", row.norm = row.norm,  
					cmap.type = cmap.type, char.rescale = char.rescale,  legend=F,
					tissue.names = tissue.names,
					tissue.labels = tissue.labels,
					highlight.tissue.name = highlight.tissue.name,
					ifScoreBarPlot = ifScoreBarPlot,
					MI.list.model = MI.list.pass1,
					MI.list.phen = MI.list.phen.pass1)
		} else{
			MSIG.HeatMapPlot.10(V = m2.pass1, 
					pathway.mut = bin.class.pass1,
					row.names = model.names.pass1,
					row.names2 = model.descs2.pass1, 
					col.labels = cls.labels2.pass1, 
					col.classes = cls.phen2.pass1, 
					phen.cmap = colors.list.pass1, 
					phen.names = phen.names.pass1,
					phen.names2 = phen.descs2.pass1,
					col.names = sample.names2.pass1, 
					main = paste(tissue, "- Step 1: Known KRAS Pathway Abnormalities (MI)"), 
					xlab="  ", ylab="  ", row.norm = row.norm,  
					cmap.type = cmap.type, char.rescale = char.rescale,  
					legend=F,
					ifScoreBarPlot = ifScoreBarPlot,
					MI.list.model = MI.list.pass1,
					MI.list.phen = MI.list.phen.pass1)
		}
		dev.off()
#	stop("Don't do Step 2!")
		### Begin Pass 2 ###
		
		if( #n.phen.pass1 > 2 || 
				!skip.aberration.refinement ){
			print( "--- Begin Step 2 ---")
			if( is.na(MI.list.pass1[1]) || is.na(MI.list.phen.pass1[1]) ){
				dev.off()
				print("Genomic abnormalities specified are not in the cls file of the samples provided.")
				return()
			}
			MI.thresholds = c(0.2, 0.1, 0.08, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01, 0)
			#	MI.threshold = 0.03
			ind.top.MI = vector(mode="integer")
			MI.i = 0
#	FDR.i = 0
			#	browser()
			while( length(ind.top.MI) < 1)
			{
				MI.i = MI.i + 1
#		FDR.i = FDR.i + 1
				if( MI.i > length(MI.thresholds)){
					dev.off()
					print("Selected genomic aberrations do not have 
									positive mutual information with a low enough false 
									discovery rate with the selected pathways")
					return()
				}
				ind.top.MI = which( MI.list.phen.pass1[-1] >= MI.thresholds[MI.i]) +1 #& MI.list.phen.pass1[-1] > 0 
#		) + 1
				
			}
			ind.MI.threshold = c(1, ind.top.MI)
			
			
			n.phen.pass2 = length(ind.MI.threshold)
			cls.list2.pass2 = cls.list2.pass1[ind.MI.threshold,]
			phen.names.pass2 = phen.names.pass1[ind.MI.threshold]
			cls.labels2.pass2 = cls.labels2.pass1[ind.MI.threshold,]
			
			#	browser()
			cls.list2.pass2.2 = ifelse( cls.list2.pass2 == "WT", 0, 1)
			cls.list2.pass2.3 = ifelse( cls.list2.pass2 == "DEL" | cls.list2.pass2 == "AMP", 1, 0)
			if( n.phen.pass2 > 2 ){ 
				pathway.mut.pass2 = colSums(cls.list2.pass2.2[2:n.phen.pass2,]) 
				#apply(cls.list2.pass2.2[2:n.phen.pass2,], MARGIN=2, FUN=sum)
				bin.class.pass2 = pathway.mut.pass2/length(pathway.mut.pass2)
				bin.class.pass2 = ((bin.class.pass2 - min(bin.class.pass2))/
							(max(bin.class.pass2) - min(bin.class.pass2)))
				cls.list2.pass2[1,] = ifelse( bin.class.pass2 > 0, "MUT", "WT")
			} else{
				pathway.mut.pass2 = ifelse( cls.list2.pass2.2[2,] == 1, "MUT", "WT")
				bin.class.pass2 = ifelse( pathway.mut.pass2 == "MUT", 1, 0 ) #+ runif(Ns, min=-.05, max=.05)
				cls.list2.pass2[1,] = pathway.mut.pass2
			}
			
			MI.list.pass2 = vector( length=target.ds$nRow, mode="double" )
			#browser()
			
			if(make.inbetween.heatmaps){
				#### Print Step 2's initial heatmap ###
				cls.phen2.pass1.5 <- NULL
				cls.labels2.pass1.5 = cls.labels2.pass2
				if (is.vector(cls.labels)) {
					classes <- unique(cls.list2.pass1)
					cls.phen2.pass1.5 <- classes
					cls.labels2.pass1.5 <- match(cls.list2.pass1, cls.phen2.pass1.5)
				} else {
					for (kk in 1:length(cls.list2.pass2[, 1])) {
						classes <- unique(cls.list2.pass1[kk,])
#            cls.phen2[[kk]] <- classes
						cls.phen2.pass1.5 <- c(cls.phen2.pass1.5, classes)
						cls.labels2.pass1.5[kk,] <- match(cls.list2.pass2[kk,], classes)
					}
				}
#	cls.labels2.pass1.5 = cls.labels2.pass1.5[1:n.phen.pass2,]
				
				phen.list.pass1.5 <- unlist(cls.phen2.pass1.5)
				colors.list.pass1.5 = rep( "gray", n.phen.pass1)
				colors.list.pass1.5[phen.list.pass1.5!="WT"] = cls.phen.colors[1]
				filename <- paste(results.dir, test.file.prefix, file.suffix, ".inbetween.initial.and.refined", sep="")
				pdf(file=paste(filename, ".pdf", sep=""), height = pdf.height, width = pdf.width )
				if( show.multiple.tissues ){
					MSIG.HeatMapPlot.10.show.multiple.tissues(V = m2.pass1, 
							pathway.mut = bin.class.pass2,
							row.names = model.names.pass1,
							col.labels = cls.labels2.pass1.5, 
							col.classes = cls.phen2.pass1.5, 
							phen.cmap = colors.list.pass1.5, 
							phen.names = phen.names.pass1[ind.MI.threshold],
							col.names = sample.names2.pass1, 
							main = paste(tissue, "- Post-Step 1, Pre-Step 2 (Step 1.5)"), 
							xlab="  ", ylab="  ", row.norm = row.norm,  
							cmap.type = cmap.type, char.rescale = char.rescale,  legend=F,
							tissue.names = tissue.names,
							tissue.labels = tissue.labels,
							highlight.tissue.name = highlight.tissue.name)
				} else{
					MSIG.HeatMapPlot.10(V = m2.pass1, 
							pathway.mut = bin.class.pass2,
							row.names = model.names.pass1,
							col.labels = cls.labels2.pass1.5, 
							col.classes = cls.phen2.pass1.5, 
							phen.cmap = colors.list.pass1.5, 
							phen.names = phen.names.pass1[ind.MI.threshold],
							col.names = sample.names2.pass1, 
							main = paste(tissue, "- Post-Step 1, Pre-Step 2 (Step 1.5)"), 
							xlab="  ", ylab="  ", row.norm = row.norm,  
							cmap.type = cmap.type, char.rescale = char.rescale,  legend=F) 
				}
				dev.off()
			}
			
			print("--- starting step 2! for reals this time (post making in-between heatmaps) ---")
			#	browser()
			#pdf(file=paste(tissue, n.randomizations, "randomizations.Step2", "pdf", sep="."))
#		if(skip.step1){
#			browser()
#			m2.pass1 = m
#			phen.ind = which(phen.names %in% phen.pass1)
#			cls.list2.pass2 = cls.list[c(1, phen.ind+1 ),]
#			cls.labels2.pass2 = cls.labels[c(1, phen.ind+1 ),]
#			phen.names.pass2 = c("SUMMARY", phen.names[phen.ind])
#			n.phen.pass2 =  length(phen.names.pass2)
#			bin.class.pass2 = REVEALER.normalize.zero.one(colSums(ifelse(cls.list2.pass2=="WT", 0, 1)))
#			tissue.labels.pass1 = tissue.labels
#			MI.list.phen.pass2 = vector("numeric", length=n.phen.pass2)
#			phen.descs.pass2 = vector("character", length=n.phen.pass2)
#			
#		}
			model.descs2.pass2 = vector(length = target.ds$nRow, mode="character")
			if( length(unique(bin.class.pass2)) > 1 ){
				
				if( target.ds$nRow > 1 ){
					MI.results = mutual.inf.3.v2(bin.class.pass2, m2.pass1, 
							target.vector.name="SUMMARY", 
					)
					MI.list.pass2  = MI.results$MI
					model.descs2.pass2 = sapply(MI.results$MI, signif, 3)
					m.order.pass2 = order(MI.list.pass2, decreasing=decreasing.order, na.last=TRUE)
					m2.pass2 <- m2.pass1[m.order.pass2, ]
					s.order.pass2 <- order(m2.pass2[1,], decreasing = decreasing.order )
					m2.pass2 <- m2.pass2[, s.order.pass2]
				} else{
					MI.ref = mutual.inf.2(bin.class.pass2, bin.class.pass2)
					MI.list.pass2 = MI.results = 
							mutual.inf.2(bin.class.pass2, 
									m2.pass1[1,])/MI.ref
					#MI.list.pass2  = MI.results$MI
					model.descs2.pass2 = signif(MI.list.pass2, digits=3)
					m.order.pass2 = 1 #order(MI.list.pass2, decreasing=TRUE, na.last=TRUE)
					m2.pass2 <- m2.pass1#[m.order.pass2, ]
					s.order.pass2 <- order(m2.pass2[1,], decreasing = decreasing.order )
					m2.pass2 <- t(as.matrix(m2.pass2[, s.order.pass2]))
					rownames(m2.pass2) = model.names.pass1
				}
			} else{ 
				MI.list.pass2 = rep(NA, target.ds$nRow)
				model.descs2.pass2 = rep(" - ", target.ds$nRow)
				if( target.ds$nRow > 1 ){
					loc <- match(model, model.names)
					s.order.pass2 <- order(m2.pass1[loc,], decreasing = decreasing.order)
					m2.pass2 <- m2.pass1[, s.order.pass2]
					correl <- cor(t(m2.pass2))[, loc]
					m.order.pass2 <- order(correl, decreasing=decreasing.order)
					m2.pass2 <- m2.pass2[m.order.pass2, ]
				} else{
					#loc <- match(model, model.names)
					s.order.pass2 <- order(m2.pass1[1,], decreasing = decreasing.order)
					m2.pass2 <- t(as.matrix(m2.pass1[, s.order.pass2]))
					rownames(m2.pass2) = model.names.pass1
					m.order.pass2 = 1
#		correl <- cor(t(m2.pass2))[, loc]
#		m.order.pass2 <- order(correl, decreasing=T)
#		m2.pass2 <- m2.pass2[m.order.pass2, ]
				}
			}
			
			MI.list.pass2 = MI.list.pass2[m.order.pass2]
			bin.class.pass2 = bin.class.pass2[s.order.pass2]
			tissue.labels = tissue.labels[s.order.pass2]
			model.descs2.pass2 = model.descs2.pass2[m.order.pass2]
			sample.names2.pass2 <- colnames(m2.pass2)
			model.names.pass2 <- rownames(m2.pass2)
			print(matrix(c(model.names.pass2, model.descs2.pass2), ncol=2), quote=F)
#	browser()
			if (is.vector(cls.labels2.pass2)) {
				cls.labels2.pass2 <- cls.labels2.pass2[s.order.pass2]
				cls.list2.pass2 <- cls.list2.pass2[s.order.pass2]
			} else {
				cls.labels2.pass2 <- cls.labels2.pass2[, s.order.pass2]
				cls.list2.pass2 <- cls.list2.pass2[, s.order.pass2]
			}
			
			tissue.labels.pass2 = tissue.labels.pass1[s.order.pass2]
			sample.names2 <- colnames(m2.pass2)
			
			winning.model.ind.pass2 = which(model.names.pass2[1] == rownames(m2.pass2))
			MI.list.phen.pass2 = vector("double", length=n.phen.pass2)
			phen.descs.pass2 = vector("character", length=n.phen.pass2)
			
			if( length(unique(bin.class.pass2)) > 1){
				MI.signif <- signif(MI.list.pass2[1], digits=3)
				MI.list.phen.pass2[1] = MI.list.pass2[1]
			} else{
				MI.signif <- "-"
				MI.list.phen.pass2[1] = NA
			}
			print(paste(format(phen.names.pass2[1], width=12), "mutual.inf =", MI.signif#, "  FDR =", FDR.signif
					))
			print(proc.time()-t1)
			print(date())
			phen.descs.pass2[1] <- paste(MI.signif) #, " (FDR = ", FDR.signif, ")", sep="")
			
			if( n.phen.pass2 == 2 ){
				phen.descs.pass2[2] <- phen.descs.pass2[1] 
#		FDR.list.phen.pass2[2] = FDR.list.phen.pass2[1] 
				MI.list.phen.pass2[2] = MI.list.phen.pass2[1] 
				g.order.pass2 = c(1,2)
				print(paste(format(phen.names.pass2[2], width=12), "mutual.inf =", MI.signif)) #, "  FDR =", FDR.signif))
			} else{
				bin.gene.matrix = ifelse(cls.list2.pass2[-1,]=="WT", 0, 1)
				n.aberrations = rowSums(bin.gene.matrix) # apply(bin.gene.matrix, MARGIN=1, FUN=sum)
				
				MI.results = mutual.inf.3.v2( 
						m2.pass2[winning.model.ind.pass2,],
						bin.gene.matrix,
						target.vector.name=phen.pass2
#				n.randomizations = n.randomizations
				)
				
				MI.list.phen.pass2[-1] = MI.results$MI
				phen.descs.pass2[-1] = sapply(MI.results$MI, signif, 3)
				ind.zeros = which(n.aberrations==0) + 1
				MI.list.phen.pass2[ind.zeros] = NA
#		FDR.list.phen.pass2[ind.zeros] = NA
				phen.descs.pass2[ind.zeros] = " - "
				g.order.pass2 = c(1, order(MI.list.phen.pass2[-1], 
								decreasing=decreasing.order, na.last=TRUE)+1)  # skip PATHWAY.MUT
			}
			#dev.off()
			phen.descs2.pass2 = phen.descs.pass2[g.order.pass2]
			cls.list2.pass2 = cls.list2.pass2[g.order.pass2,]
			phen.names.pass2 = phen.names.pass2[g.order.pass2]
			#	browser()
			# Recompute cls.list2 as some mutations or copy numbers may have been removed
			print(matrix(c(phen.names.pass2, phen.descs2.pass2), ncol=2), quote=F)
			
			
			# Recompute cls.phen and cls.labels2 as order may have changed
			
			cls.phen2.pass2 <- NULL
			if (is.vector(cls.labels)) {
				classes <- unique(as.vector(cls.list2.pass2))
				cls.phen2.pass2 <- classes
				cls.labels2.pass2 <- match(cls.list2.pass2, cls.phen2.pass2)
			} else {
				for (kk in 1:length(cls.list2.pass2[, 1])) {
					classes <- unique(cls.list2.pass2[kk,])
					#            cls.phen2[[kk]] <- classes
					cls.phen2.pass2 <- c(cls.phen2.pass2, classes)
					cls.labels2.pass2[kk,] <- match(cls.list2.pass2[kk,], classes)
				}
			}
			cls.labels2.pass2 = cls.labels2.pass2[1:n.phen.pass2,]
			
			phen.list.pass2 <- unlist(cls.phen2.pass2)
			colors.list.pass2 = rep( "gray", length(phen.list.pass2))
			colors.list.pass2[phen.list.pass2=="MUT"] = cls.phen.colors[1]
#		colors.list.pass2[phen.list.pass2=="DEL"] = cls.phen.colors[3]
#		colors.list.pass2[phen.list.pass2=="AMP"] = cls.phen.colors[4]
#		colors.list.pass2[phen.list.pass2=="ALT"] = cls.phen.colors[5]
			
			filename <- paste(results.dir, test.file.prefix, file.suffix, ".refined.aberrations", sep="")
			pdf(file=paste(filename, ".pdf", sep=""), height = pdf.height, width = pdf.width )
			if( show.multiple.tissues ){
				MSIG.HeatMapPlot.10.show.multiple.tissues(V = m2.pass2, 
						pathway.mut = bin.class.pass2,
						row.names = model.names.pass2,
						row.names2 = model.descs2.pass2, 
						col.labels = cls.labels2.pass2, 
						col.classes = cls.phen2.pass2, 
						phen.cmap = colors.list.pass2, phen.names = phen.names.pass2,
						phen.names2 = c(" ", phen.descs2.pass2),
						col.names = sample.names2.pass2, 
						main = paste(tissue, 
								"- Step 2: only MI >=", MI.thresholds[MI.i],"from Step 1 (MI)"), 
						xlab="  ", ylab="  ", row.norm = row.norm,  
						cmap.type = cmap.type, char.rescale = char.rescale,  legend=F,
						tissue.names = tissue.names,
						tissue.labels = tissue.labels,
						highlight.tissue.name = highlight.tissue.name,
						ifScoreBarPlot = ifScoreBarPlot,
						MI.list.model = MI.list.pass2,
						MI.list.phen = MI.list.phen.pass2)
			} else{
				MSIG.HeatMapPlot.10(V = m2.pass2, 
						pathway.mut = bin.class.pass2,
						row.names = model.names.pass2,
						row.names2 = model.descs2.pass2, 
						col.labels = cls.labels2.pass2, 
						col.classes = cls.phen2.pass2, 
						phen.cmap = colors.list.pass2, phen.names = phen.names.pass2,
						phen.names2 = phen.descs2.pass2,
						col.names = sample.names2.pass2, 
						main = paste(tissue, 
								"- Step 2: only MI >=", MI.thresholds[MI.i],"from Step 1 (MI)"), 
						xlab="  ", ylab="  ", row.norm = row.norm,  
						cmap.type = cmap.type, char.rescale = char.rescale,  
						legend=F,
						ifScoreBarPlot = ifScoreBarPlot,
						MI.list.model = MI.list.pass2,
						MI.list.phen = MI.list.phen.pass2)
			}
			dev.off()
		} 
	} else{ print("'de.novo.search' on -- skipping Steps 1 and 2 and simply 'filling in' from scratch")}
	### 3rd Pass ###	
	if(!align.only.no.discovery #|| de.novo.search
			){
		print( "--- Begin Pass 3 (Iterative Method)---")
		print("2 in explained vector = previous explained vector   1 = new additions")
		if( de.novo.search && !skip.aberration.refnement ){
			model.names.pass2 = rownames(m)
			model.ind = which(model.names.pass2 == model)
			m[c(1,model.ind),] = m[c(model.ind,1),]
			s.order.pass3 = order(m[1,], decreasing=decreasing.order)
			if( target.ds$nRow > 1 ){
				m = m[,s.order.pass3]
				cls.list2.pass3 = cls.list2.pass2 = cls.list[,s.order.pass3]
				cls.labels2.pass3 = cls.labels2.pass2 = cls.labels[,s.order.pass3]
			} else{ m = t(as.matrix(m[,s.order.pass3])); rownames(m) = model.names.pass2
				cls.list2.pass3 = cls.list2.pass2 = t(as.matrix(cls.list[,s.order.pass3]))
				cls.labels2.pass3 = cls.labels2.pass2 = t(as.matrix(cls.labels[,s.order.pass3]))
			}
			m2.pass3 = m2.pass2 = m
			
			#cls.list2.pass3 = cls.list2.pass2 = cls.list
			
			phen.names.pass2 = phen.names
			file.suffix = paste(file.suffix, "_de.novo.search", sep="")
		} else if( skip.aberration.refinement && !de.novo.search ){
			m2.pass3 = m2.pass2 = m2.pass1
			model.names.pass2 = model.names.pass1
			cls.list2.pass2 = cls.list2.pass1
			cls.labels2.pass2 = cls.labels2.pass1
			cls.list2.pass3 = cls.list[, s.order.pass1]
			cls.labels2.pass3 = cls.labels[,s.order.pass1]
			phen.names.pass2 = phen.names.pass1
		}else if ( !skip.aberration.refinement && !de.novo.search){
			m2.pass3 = m2.pass2
			cls.list2.pass3 = cls.list[, s.order.pass1][, s.order.pass2]
			cls.labels2.pass3 = cls.labels[,s.order.pass1][,s.order.pass2]
			
		}
		model.names.pass3 = rownames(m2.pass3)
		sample.names2.pass3 = colnames(m2.pass3)
		n.phen.pass3 = 40
		
		top.genes.ind = NULL
		top.genes.names = NULL
		top.genes.vectors = NULL
		top.genes.MI = NULL
		top.diffs = NULL
		explained.vectors = NULL
#	browser()
		bin.gene.matrix.3 =  ifelse(cls.list2.pass3[-1,]=="WT", 0, 1)
		m2.pass2.1.median = median(m2.pass2[1,])
		mid.point <- which.min(abs(m2.pass2[1,] - quantile(m2.pass2[1,], 0.5)))
		grey.and.black = c("#C0C0C0", "#000000")
		pathway.name = model.names.pass2[1]
		MI.ref = mutual.inf.2(m2.pass3[1,], m2.pass3[1,])
		
		mycol <- vector(length=512, mode = "numeric")
		for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
		for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
		mycol <- rev(mycol)
		
		if( !de.novo.search && skip.aberration.refinement ){
			explained.initial = ifelse(
					bin.class.pass1 == 0, 0, 1) 
			explained.vectors = rbind(explained.vectors, explained.initial)
			explained = explained.initial
			explained.MI.initial = mutual.inf.2(explained, m2.pass2[1,])/MI.ref
			print(paste("explained.MI.initial =", explained.MI.initial))
			print(explained)
		} else if( !de.novo.search && !skip.aberration.refinement ){
			explained.initial = ifelse(
					bin.class.pass2 == 0, 0, 1) 
			explained.vectors = rbind(explained.vectors, explained.initial)
			explained = explained.initial
			explained.MI.initial = mutual.inf.2(explained, m2.pass2[1,])/MI.ref
			print(paste("explained.MI.initial =", explained.MI.initial))
			print(explained)
		}
		
		cex.axis = 1
		ncolors <- length(mycol)
		
		initial.gene.string = ifelse(de.novo.search, phen.pass1[1], phen.names.pass2[2])
		
		pdf(file=paste(results.dir, test.file.prefix, file.suffix, ".",
						initial.gene.string,
						".REVEALER.pdf", sep=""), 
				height=8.5, width=11)
		#par(mar = c(1, 15, 1, 5))
		
		
		## If we had naively searched the space without removing the explained cell lines
		MI.results = mutual.inf.3.v2(m2.pass2[1,], bin.gene.matrix.3)
		
		#browser()
		MI.order = order(MI.results$MI, decreasing=decreasing.order, na.last=TRUE)+1
		gc()
#		browser()
		write.table(matrix(c(c("Genomic Aberrations", phen.names[MI.order]), 
								c(paste("Normalized Mutual Information to", pathway.name), MI.results$MI[MI.order-1])),
						ncol=2), 
				file = paste(results.dir, test.file.prefix, file.suffix, ".Step3_naive.txt", sep=""),
				quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE, append=FALSE)
		top10.names = c( #paste(c(phen.names.pass2[-1], "(from Step 2)"), collapse="  " ), 
				phen.names[MI.order[1:40]])
		top10.MI = c( #signif(explained.MI.initial, digits=4), 
				signif(MI.results$MI[MI.order[1:40]-1], digits=4))
		top10.labels = rbind(#explained+1, 
				bin.gene.matrix.3[MI.order[1:40]-1,]+1)
		num.redundant = length(which(MI.results$MI == top10.MI[1]))
		par(mar = c(2, 16, 2, 5))
		nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(2, 10), FALSE)
		max.v <- max(max(m2.pass2[1,]), -min(m2.pass2[1,]))
		if( match.to.blue){
			V1 <- c( (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,1:mid.point]), 
					+ (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,(mid.point+1):Ns])+ncolors/2 )
		} else{
			V1 <- c( (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,1:mid.point]) + ncolors/2, 
					(ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,(mid.point+1):Ns]))
		}
		image(1:length(m2.pass2[1,]), 1:1, as.matrix(V1), 
				zlim = c(0, ncolors), col=mycol, axes=FALSE, 
				main="Naive REVEALER without exclusion of Step 2 aberrations", sub = "", xlab= "", ylab="")
		axis(2, at=1:1, labels=pathway.name, adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
		V1 <- top10.labels
		V1 <- apply(V1, MARGIN=2, FUN=rev)      
#		max.v <- max(max(V1), -min(V1))
#		V1 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
		image(1:length(m2.pass2[1,]), 1:dim(V1)[1], t(V1), 
				zlim = c(0, length(grey.and.black)), col=grey.and.black, axes=FALSE, 
				main="", #paste("step:", i), 
				sub = "", xlab= "", ylab="")
		axis(2, at=1:dim(V1)[1], labels=rev(top10.names), adj= 0.5, tick=FALSE, las = 1, cex.axis=1, 
				font.axis=1, line=-1)
		axis(4, at=1:dim(V1)[1], labels=rev(top10.MI), adj= 0.5, tick=FALSE, las = 1, cex.axis=1, 
				font.axis=1, line=-1)
		
		if( de.novo.search && u.gene.names.known == "NA" ){
			print("top10.labels[1:5,] :")
			print(top10.labels[1:5,])
#			browser()
			explained.prev = explained.initial = bin.gene.matrix.3[MI.order[1]-1,]
			explained.MI.prev = explained.MI.initial = top10.MI[1]
			explained.vectors = rbind(explained.vectors, explained.initial)
#			browser()
			initial.genes = paste(top10.names[1], "+", num.redundant, 
					ifelse(num.redundant==1, "other", "others"))
		} else if(de.novo.search && u.gene.names.known != "NA"){  # use provided gene sets as initial explained set
			if(length(u.gene.names.known) > 1){
				summary = colSums(bin.gene.matrix.3[ which(phen.names %in% phen.pass1), ])
				summary = REVEALER.normalize.zero.one(summary)
				summary.MI = mutual.inf.2(m2.pass3[1,], summary)/MI.ref
				explained.prev = explained.initial = ifelse(summary > 0, 1, 0)
				phen.pass1.MIs = mutual.inf.3.v2(m2.pass3[1,], 
						bin.gene.matrix.3[ which(phen.names %in% phen.pass1)-
										ifelse(phen.pass1[1]=="SUMMARY", 1, 0), ])$MI
			} else{
				explained.prev = explained.initial = summary = bin.gene.matrix.3[ 
						which(phen.names %in% phen.pass1), ]
			}
			explained.MI.prev = explained.MI.initial = mutual.inf.2(m2.pass3[1,], explained.prev)/MI.ref
			initial.genes = paste(phen.pass1, collapse="  ")
#			browser()
#			cls.phen2.pass3 <- NULL
#			if (is.vector(cls.labels)) {
#				classes <- unique(cls.list2.pass3[ 
#								(which(phen.names %in% phen.pass1)+ifelse(phen.pass1[1]=="SUMMARY", 0, 1)), ])
#				cls.phen2.pass3 <- classes
#				cls.labels.pass3 <- match(cls.list2.pass3[ 
#								(which(phen.names %in% phen.pass1)+ifelse(phen.pass1[1]=="SUMMARY", 0, 1)), ], 
#						cls.phen2.pass3)
#			} else {
#				#browser()
#				for (kk in (which(phen.names %in% phen.pass1) +
#							ifelse(phen.pass1[1]=="SUMMARY", 0, 1) ) ) {
#					classes <- unique(cls.list2.pass3[kk,])
			##            cls.phen2[[kk]] <- classes
#					cls.phen2.pass3 <- c(cls.phen2.pass3, classes)
#					cls.labels2.pass3[kk,] <- match(cls.list2.pass3[kk,], classes)
#				}
#			}
#			
#			browser()
#			phen.list.pass3 <- unlist(cls.list2.pass3[ 
#							(which(phen.names %in% phen.pass1)+1), ])  # Have to add +1 because added a row for "summary"
#			colors.list.pass3 = rep( "gray", length(phen.list.pass3))
#			colors.list.pass3[phen.list.pass3=="MUT"] = cls.phen.colors[1]
#			MSIG.HeatMapPlot.10(V = m2.pass2[1,], 
#					pathway.mut = summary,
#					row.names = pathway.name,
#					row.names2 = paste(1), 
#					col.labels = cls.labels2.pass3[ which(phen.names %in% phen.pass1), ], 
#					col.classes = cls.phen2.pass3[ which(phen.names %in% phen.pass1), ], 
#					phen.cmap = c("gray", "black"),#colors.list.pass3, 
#					phen.names = phen.pass1,
#					phen.names2 = signif(c(summary.MI, phen.pass1.MIs), digits=3),  ## Calculate MI of individual genes with pathway
#					col.names = sample.names2.pass2, 
#					main = paste(tissue, 
#							"Initial signature with provided genes"), 
#					xlab="  ", ylab="  ", 
#					row.norm = row.norm,  
#					cmap.type = cmap.type, char.rescale = char.rescale,  legend=F)
			
			explained.vectors = rbind(explained.vectors, explained.initial)
		} else{
			explained.prev = explained
			explained.MI.prev = explained.MI.initial
#			explained.vectors = rbind(explained.vectors, explained.initial)  # this is taken care of at the beginning of REVEALER
		}
		samples.without.mut = ifelse(
				explained.prev == 1, 0, 1)
#		browser()
		
		# "wo.mut.or.nonMatchingSide" :
		# without mutations or not on target side to match to, side being "red" or "blue",
		# i.e. align to positive ("red") or negative ("blue") pathway expression
		# (samples to find "explanations"/genomic aberrations for)
		# "wo.mut.and.matchingSide" :
		# with mutations and on target side (samples to cut out)
		if (match.to.blue){
			wo.mut.or.nonMatchingSide = (samples.without.mut | m2.pass2[1,] >= m2.pass2.1.median)
#					ifelse(c(rep(1, length=(Ns - mid.point)),
#							samples.without.mut[(mid.point+1):Ns])==1, TRUE, FALSE)
			wo.mut.and.matchingSide = (samples.without.mut & m2.pass2[1,] < m2.pass2.1.median)
#					ifelse(c(rep(0, length=(Ns - mid.point)),
#							samples.without.mut[1:mid.point] )==1, TRUE, FALSE)
		} else{
			wo.mut.or.nonMatchingSide = (samples.without.mut | m2.pass2[1,] <= m2.pass2.1.median)
#					ifelse(c(samples.without.mut[1:mid.point], 
#							rep(1, length=(Ns - mid.point)) )==1, TRUE, FALSE)
			wo.mut.and.matchingSide = (samples.without.mut & m2.pass2[1,] > m2.pass2.1.median)
#					ifelse(c(samples.without.mut[1:mid.point], 
#							rep(0, length=(Ns - mid.point)) )==1, TRUE, FALSE)
		}
		n.iter = 0
		MI.diff.latest = explained.MI.prev
		#browser()
		print(paste("initial # samples:", length(explained)))
		while( ifelse(match.to.blue, -1, 1)*MI.diff.latest > 0  && n.iter <=5 
				){
			#for( i in 1:n.iter){
			n.iter = n.iter +  1
			i = n.iter
			gc()
			print(paste("iteration:", i))
			print(paste("# samples left: ", sum(wo.mut.or.nonMatchingSide), "/", length(wo.mut.or.nonMatchingSide), 
							"(", length(wo.mut.or.nonMatchingSide)-sum(wo.mut.or.nonMatchingSide),"removed)", sep=""))
			MI.results = mutual.inf.3.v2( 
					m2.pass2[1,wo.mut.or.nonMatchingSide], 
					bin.gene.matrix.3[,wo.mut.or.nonMatchingSide] )
			MI.order = order(MI.results$MI, decreasing=decreasing.order, na.last=TRUE)+1
			
			write.table(matrix(c(c("Genomic Aberrations", phen.names[MI.order]), 
									c(paste("Normalized Mutual Information to", pathway.name), MI.results$MI[MI.order-1])),
							ncol=2), 
					file = paste(results.dir, test.file.prefix, file.suffix, ".", initial.gene.string,
							".Step3_iter", i, ".txt", sep=""),
					quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
			top10.names = phen.names[MI.order[1:10]]
			top10.MI = MI.results$MI[MI.order[1:10]-1]
			top10.labels = bin.gene.matrix.3[MI.order[1:10]-1,wo.mut.or.nonMatchingSide] +1	
			
			top.genes.ind = c(top.genes.ind, MI.order[1] )
			num.redundant = length(which(MI.results$MI == MI.results$MI[MI.order[1]-1]))-1
			top.genes.names = c( top.genes.names, paste(phen.names[MI.order[1]], "+", num.redundant, 
							ifelse(num.redundant==1, "other", "others")))
			
			mut = bin.gene.matrix.3[(MI.order[1]-1),]
			explained = ifelse(mut+explained.prev>0, 1, 0)
			explained.MI = mutual.inf.2( m2.pass2[1,wo.mut.or.nonMatchingSide], 
					explained[wo.mut.or.nonMatchingSide])/MI.ref
			MI.diff = MI.diff.latest = explained.MI - explained.MI.prev
			print(paste("Explained.MI = ", explained.MI, 
							"  MI.diff = ", ifelse(MI.diff<0, "-", "+"), 
							signif(abs(MI.diff), digits=4), sep=""))
			explained.vectors = rbind(explained.vectors, explained)
			print(2*explained.prev + mut)
			
			top.diffs = c(top.diffs, MI.diff)
			top.genes.vectors = rbind(top.genes.vectors, mut)
			top.genes.MI = c(top.genes.MI, paste(signif(MI.results$MI[MI.order[1]-1], digits=4), 
							sep="" ))
			
			
			par(mar = c(2, 16, 2, 12))
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(2, 10), FALSE)
			max.v <- max(max(m2.pass2[1,wo.mut.or.nonMatchingSide]), -min(m2.pass2[1,wo.mut.or.nonMatchingSide]))
			if( match.to.blue){
				V1 <- c( (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,1:mid.point]), 
						+ (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,(mid.point+1):Ns])+ncolors/2 )[wo.mut.or.nonMatchingSide]
			} else{
				V1 <- c( (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,1:mid.point]) + ncolors/2, 
						(ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,(mid.point+1):Ns]))[wo.mut.or.nonMatchingSide]
			}
#			V1 <- c( (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,1:mid.point]) + ncolors/2, 
#					(ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,(mid.point+1):Ns]))
			
			image(1:length(m2.pass2[1,wo.mut.or.nonMatchingSide]), 1:1, as.matrix(V1), 
					zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
			axis(2, at=1:1, labels=pathway.name, adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
			V1 <- rbind( explained[wo.mut.or.nonMatchingSide]+1, top10.labels)
			V1 <- apply(V1, MARGIN=2, FUN=rev)      
			image(1:length(m2.pass2[1,wo.mut.or.nonMatchingSide]), 1:dim(V1)[1], t(V1), 
					zlim = c(0, length(grey.and.black)), col=grey.and.black, 
					axes=FALSE, main=paste("iteration:", i), sub = "", xlab= "", ylab="")
			axis(2, at=1:dim(V1)[1], labels=rev(c("explained with top result", top10.names)), 
					adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
			axis(4, at=1:dim(V1)[1], labels=rev(c( paste(signif(explained.MI, digits=4), 
											sep="" ), 
									signif(top10.MI,digits=4))), 
					adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
			
			samples.without.mut[wo.mut.or.nonMatchingSide] = samples.without.mut[
					wo.mut.or.nonMatchingSide] - mut[wo.mut.or.nonMatchingSide]
			if( match.to.blue ){
				wo.mut.or.nonMatchingSide = (samples.without.mut | m2.pass2[1,] >= m2.pass2.1.median)
				wo.mut.and.matchingSide = (samples.without.mut & m2.pass2[1,] < m2.pass2.1.median)
			} else{
				wo.mut.or.nonMatchingSide = (samples.without.mut | m2.pass2[1,] <= m2.pass2.1.median)
				wo.mut.and.matchingSide = (samples.without.mut & m2.pass2[1,] > m2.pass2.1.median)
			}
			explained.MI.prev = explained.MI
			explained.prev = ifelse(mut+explained.prev>0, 1, 0)
			print(proc.time()-t1)
			print(date())
			
		}
		
#		browser() # plotting isn't working with KRAS RNAi data - too many "cumulative" labels?
		
		explained = ifelse( colSums(rbind(explained.initial, 
								top.genes.vectors)) >=1,
				#apply(rbind(ifelse(cls.list2.pass2[1,]=="WT", 0,1), 
				#				top.genes.vectors), MARGIN=2, FUN=sum)>=1, 
				1, 0)
		explained.MI = mutual.inf.2(m2.pass2[1,], explained)/MI.ref
		top.genes.MI = signif(mutual.inf.3.v2(m2.pass2[1,], top.genes.vectors)$MI, digits=4)
		par(mar = c(2, 16, 2, 12))
		nf <- nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(2, 10), respect = FALSE)
		max.v <- max(max(m2.pass2[1,]), -min(m2.pass2[1,]))
		if( match.to.blue){
			V1 <- c( (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,1:mid.point]), 
					+ (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,(mid.point+1):Ns])+ncolors/2 )
		} else{
			V1 <- c( (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,1:mid.point]) + ncolors/2, 
					(ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,(mid.point+1):Ns]))
		}
#		V1 <- c( (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,1:mid.point]) + ncolors/2, 
#				(ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,(mid.point+1):Ns]))
		image(1:length(m2.pass2[1,]), 1:1, as.matrix(V1), 
				zlim = c(0, ncolors), col=mycol, axes=FALSE, 
				main="Final results from REVEALER iterations", sub = "", xlab= "", ylab="")
		axis(2, at=1:1, labels=pathway.name, adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
		V1 <- rbind(explained.initial, top.genes.vectors, explained)+1
		#rbind(ifelse(cls.list2.pass2[1,]=="WT", 0,1), top.genes.vectors, explained)+1
		V1 <- apply(V1, MARGIN=2, FUN=rev)      
		image(1:length(m2.pass2[1,]), 1:dim(V1)[1], t(V1), 
				zlim = c(0, length(grey.and.black)), col=grey.and.black, axes=FALSE, 
				sub = "", xlab= "", ylab="")
#		browser()
		axis(2, at=1:dim(V1)[1], labels=rev(c(paste("seed:", ifelse(de.novo.search, initial.genes, 
												paste(phen.names.pass2[-1], collapse="  ")) ), 
								paste("iter:", ifelse(n.iter > 1, n.iter-1, 1),top.genes.names), 
								"explained")), adj= 0.5, tick=FALSE, 
				las = 1, cex.axis=1, font.axis=1, line=-1)
		axis(4, at=1:dim(V1)[1], labels=rev( c(signif(explained.MI.initial, digits=4), 
								paste(top.genes.MI, sep=""), 
								signif(explained.MI,digits=4))), adj= 0.5, tick=FALSE, 
				las = 1, cex.axis=1, font.axis=1, line=-1)
		### End plotting individual genomic aberrations 
		
		
		explained.vectors.MI = mutual.inf.3.v2(m2.pass2[1,], explained.vectors)$MI
		MI.diffs = explained.vectors.MI[-1] - explained.vectors.MI[1:n.iter]
		
		par(mar = c(2, 16, 2, 12))
		nf <- nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(2, 10), respect = FALSE)
		max.v <- max(max(m2.pass2[1,]), -min(m2.pass2[1,]))
		if( match.to.blue){
			V1 <- c( (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,1:mid.point]), 
					+ (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,(mid.point+1):Ns])+ncolors/2 )
		} else{
			V1 <- c( (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,1:mid.point]) + ncolors/2, 
					(ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,(mid.point+1):Ns]))
		}
#		V1 <- c( (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,1:mid.point]) + ncolors/2, 
#				(ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,(mid.point+1):Ns]))
		image(1:length(m2.pass2[1,]), 1:1, as.matrix(V1), 
				zlim = c(0, ncolors), col=mycol, axes=FALSE, 
				main="Final results from REVEALER iterations", sub = "", xlab= "", ylab="")
		axis(2, at=1:1, labels=pathway.name, adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
		V1 = apply(explained.vectors+1, MARGIN=2, FUN=rev)
		image(1:length(m2.pass2[1,]), 1:dim(V1)[1], t(V1), 
				zlim = c(0, length(grey.and.black)), col=grey.and.black, axes=FALSE, 
#			main=paste("step:", i), 
				sub = "", xlab= "", ylab="")
		if( n.iter > 1){
			left.labels = c("INITIAL cumulative", paste("cumulative, iter: ", 1:(n.iter-1), " ", sep=""),
					paste("FINAL cumulative, iter: ", n.iter, " ",sep=""))
		} else{
			left.labels = c("INITIAL cumulative",
					paste("FINAL cumulative, iter: ", n.iter, " ",sep=""))
		}
		right.labels = c( paste(" ", signif(explained.vectors.MI[1],digits=4), sep=""), 
				paste( " ", signif(explained.vectors.MI[2:(n.iter+1)], digits=4), 
						" (",ifelse(MI.diffs < 0, "-", "+"), 
						signif(abs(MI.diffs), digits=4),")", sep="") 
		)
		axis(2, at=1:dim(V1)[1], labels=rev(left.labels), adj= 0.5, tick=FALSE, las = 1, cex.axis=1, 
				font.axis=1, line=-1)
		axis(4, at=1:dim(V1)[1], labels=rev(right.labels), adj= 0.5, tick=FALSE, 
				las = 1, cex.axis=1, font.axis=1, line=-1)
		
		cls.labels2.pass3 = rbind(cls.labels2.pass2[1,], top.genes.vectors+1)
		cls.list2.pass3 = rbind(cls.list2.pass2[1,], ifelse(top.genes.vectors==0, "WT", "MUT"))
		
		MI.results = mutual.inf.3.v2( explained,
				m2.pass2)  ## must subtract 1 from the indices because bin.gene.matrix.3 
		## doesn't include SUMMARY
		#target.vector.name=phen.pass1[ind.master],
#			n.randomizations = n.randomizations)
#	g.order.pass3.top40
		phen.descs.pass3 = ifelse( is.nan(MI.results$MI), " - ", signif(MI.results$MI, digits=3))
#	phen.descs2.pass3 = phen.descs.pass3[g.order.pass3]
		print(proc.time()-t1)
		print(date())
		
		dev.off()
#else{ print("skipping iterative method!")}
		
	}
	
	if (!is.na(output.dataset)) {
#		V.GCT <- m.all
#		print("Figure out why model.descs2.pass1 does not correspond to the rows of m.all")
#		browser()
#		colnames(V.GCT) <- sample.names2
#		row.names(V.GCT) <- model.names.pass1
		write.gct(gct.data.frame = m, descs = model.descs2.pass1, filename =output.dataset)  
	}
	
}

REVEALER <- function(
		###  "target" = typically continuous readout such as expression, gene set enrichment, drug/RNAi 
		### 		    sensitivity, degree of phosphorylation .... some values that you would like to
		###				"explain" with the features
		target.ds.file,
		target.ds.file.type = c("gct", "cls", "txt.sample.rows", "txt.sample.cols"), 
		## 											"gct" is an expression file
		##											"cls" is a class/phenotype file
		##											"txt.sample.rows" is a tab-delimited file with sample names on the row
		##												(i.e. first column is sample names)
		##											"txt.sample.cols" is a tab-delimited file with sample names on the columns
		##												(i.e. first row is sample names)
		target.names.to.match = c("NA", "ALL"),
		target.remove.na.samples = FALSE, ## target columns are the samples, ie the things that have a value for a target
		target.remove.na.targets = FALSE, ## target rows are the targets, ie the thing(s) you're trying to explain
		target.if.unmatched = "NA",  ## if no samples have the features specified, which target name to sort by
		
		### "feature" = continuous or binary data. Could be genomic such as mutation, copy number, 
		### 			or epigenomic such as methylation, phosphorylation status, etc
		feature.ds.file,
		feature.ds.file.type = c("gct", "cls", "txt.sample.rows", "txt.sample.cols"),
		
		##### Begin UNTESTED arguments #####
		feature.remove.na.samples = FALSE, ## feature columns are samples, ie the things that have the features
		feature.remove.na.features = FALSE, ## feature rows are features, ie the things you're searching for
		## 											"gct" is an expression file
		##											"cls" is a class/phenotype file
		##											"txt.sample.rows" is a tab-delimited file with sample names on the row
		##												(i.e. first column is sample names)
		##											"txt.sample.cols" is a tab-delimited file with sample names on the columns
		##												(i.e. first row is sample names)
		###### End UNTESTED arguments #####
		
		
		feature.seed.names.to.match = "NA",
		feature.annotations.file = NULL, ## Tab-delimited file with feature names from feature.ds.file in the first
		##									column (naming the rows) and all the rest of the columns will be collapsed
		##									to a single annotation of the corresponding gene.
		
#		feature.suffix.or = c("", "_MUT", "mut"),	##  OR suffix. Suffixes one of "", "_MUT" or "mut", or whatever the 
#		## 											user specified to feature names. 
#		## 											Can be helpful if you just want to specify 
#		## 											a gene list but your dataset says "GENE_MUT"
#		## 											Setting to "" suffixes nothing (retains original naming)
#		feature.suffix.and = c("_AMP", "_DEL"),	##  AND suffix. Suffixes all of "_AMP" and "_DEL" (or whatever 
#		##											the user specified) to the feature names. Set to boolean FALSE to 
#		##											retain original naming. Specifically suffixes separately
#		##											from the OR suffix (as in, you will have "GENE_MUT" and "GENE_AMP"
#		##											separately, but no "GENE_MUT_AMP")
#		feature.add.chrom.locs = FALSE, ##  if feature names are in form "GENE_MUT" or "GENE" then the 
#		## 									chromosome location can be looked up using the BioconductoR
#		## 									library "org.Hs.eg.db" If features are not in a known gene id format
#		##									that can be recognized, will be "NA." Can be helpful to keep track of genes
		feature.skip.refinement = FALSE, ## FALSE: You want to use all the features specified in "feature.seed.names.to.match" and
		##									don't want to do any refinement of your features
		##								TRUE: You're not sure if all your features are relevant so you just want 
		##									to pick and use the top ones for feature selection
		feature.remove.from.selection.search.space = NULL, ## Remove from "feature.ds.original" and don't consider them when doing
		##											REVEALER.feature.selector
		match.to.blue = FALSE,   ## Match to the negative side of the target
		de.novo.search = FALSE, ## Ignore feature.seed.names.to.match and seed with the best feature from feature.ds.file
		## ("best" feature meaning has highest normalized mutual information with target)
		align.only.no.discovery = FALSE,
		target.phen = NA,
		target.class = NA,
		
		## Add these:
#		target.remove.na.samples = TRUE,   # Samples that have NA in target will be removed in feature
#		target.remove.na.targets = FALSE,  # Independent of feature
#		feature.remove.na.samples = FALSE, # Samples that have NA in feature will be removed in target
#		feature.remove.na.features = TRUE,  # Independent of target
		
		
		
		
		### Filename and results parameters
		results.dir = NULL,   ## Where do you want to find your pdfs and txt output?
		tissue = "NA",  ## For pdf naming and figure titling
		file.suffix = "",
		
		
		normalize.score = T,
		normalization.type = "zero.one",
		
		
		
		
		#n.random.target.names.to.match = 10,
		
		
		### Iterations of REVEALER, or until MI starts to decrease (or increase
		### if match.to.blue == TRUE)
		n.iter = 5,
		
		### Plotting parameters
		popup.heatmap = FALSE,  # pops up a heatmap window instead of saving to pdf (mac only)
		plot.unproductive.iteration = FALSE,
		pdf.height = 11,
		pdf.width = 17,
		highlight.tissue.name = NULL,
		show.multiple.tissues = FALSE,
		ifLegend = FALSE,
		ifScoreBarPlot = FALSE,
		if.feature.summary.on.top = FALSE,
		user.colors = NA,
#		decreasing.order = T,
		output.dataset = NA,
		char.rescale = 1,
		cmap.type = 3,
		row.norm = T,
		draw.white.lines = FALSE,
		
		#do.mRMR = FALSE,
		#skip.step1 = FALSE,
		
#		continuous.data = FALSE,
		
		### Debugging params
		make.inbetween.heatmaps = FALSE,  ## for debugging and preparing intermediate steps for presentations
		exclude.samples = NA,     ## For cleaning up and debugging datasets
		kde2d.timing = FALSE

)
# Calls MSIG.HeatMapPlot.9 and makes a plot sorted by the highest-scoring
# target.names.to.match and abnormalities (gene mutations or copy number alterations)
# i.e. doesn't require a "model" to score by as OPAM.sort.projection.by.score.2 does.
# However, it *will* use "model" if it cannot calculate p-values on the gene target.names.to.match, which
# happens when every cell line has a genomic aberration.
#
# Runs 3 passes on the data:
# 1st pass: looks at the genes and copy number alterations specified by feature.seed.names.to.match
# 2nd pass: looks at only the top abnormalities (using a p-value cutoff) from the 1st pass, and adjusts 
# the PATHWAY.MUT vector accordingly (only according to the genes, not by copy number data)
# 3rd pass: Takes the winning signature from the 2nd pass and then looks all the genes available
#
# Very similar to OPAM.sort.projection.by.score.4, however this version uses mutual.inf instead of
# roc.area to calculate mutual information scores and p-values for PATHWAY.MUT, the vector of total genomic aberrations
# in all samples
#
# Differs from OPAM.sort.projection.by.score.6 by requiring the gct file of expression in 
# all pathways by the input tissue ("target.ds")
{
	
	REVEALER.load.libraries(c("gtools", "verification", "ROCR", "MASS", "RColorBrewer", "heatmap.plus"))
	
	decreasing.order = ifelse(match.to.blue, FALSE, TRUE)
#	if(match.to.blue){
	file.suffix = paste(file.suffix, ifelse(match.to.blue, ".match.to.blue", ".match.to.red"), sep="")
	
#	}
#	if( feature.add.chrom.locs ) file.suffix = paste(file.suffix, "_with.chrom.locs", sep="")
	
	ds.parser.lookup = c(REVEALER.parse.gct, REVEALER.parse.cls, REVEALER.parse.txt.sample.rows,
			REVEALER.parse.txt.sample.cols)
	names(ds.parser.lookup) = c("gct", "cls", "txt.sample.rows", "txt.sample.cols")
	
	target.parser = eval(substitute(ds.parser.lookup$file.type, list(file.type = as.name(target.ds.file.type))))
	target.ds.original <- target.parser(filename = target.ds.file)  # "ds" for "dataset"
	
	feature.parser = eval(substitute(ds.parser.lookup$file.type, list(file.type = as.name(feature.ds.file.type))))
	feature.ds.original = feature.parser(filename = feature.ds.file)  # "ds" for "dataset"

	
	if(feature.ds.original$nSample != target.ds.original$nSample){
		stop("feature dataset and target dataset do not have the same number of samples")
	}
	
	#browser()
	
	if(!is.na(exclude.samples)){
		target.ds.original = REVEALER.exclude.samples(target.ds.original, exclude.samples)
		feature.ds.original = REVEALER.exclude.samples(feature.ds.original, exclude.samples)		
	}
	if(!is.null(feature.remove.from.selection.search.space)){
		feature.ds.original = REVEALER.exclude.rows(feature.ds.original, 
				feature.remove.from.selection.search.space)
	}
	feature.ds.original$sample.names = target.ds.original$sample.names
	colnames(feature.ds.original$matrix) = target.ds.original$sample.names
	
	
	tissues = REVEALER.check.multiple.tissues(target.ds.original, show.multiple.tissues, tissue)
	results.dir = REVEALER.fix.results.dir(results.dir)
	extracted.file.prefix <- REVEALER.get.extracted.file.prefix(target.ds.file, target.ds.file.type)
	if(target.names.to.match[1] != "ALL" && target.names.to.match!="NA"){
		target.ds = REVEALER.extract.row.names.to.match(target.ds.original, target.names.to.match)
	}else if(target.names.to.match == "ALL"){
		target.ds = target.ds.original
	} else{
		stop("Must provide target names to match")
	}
	
	
#	CLS <- MSIG.ReadPhenFile.2(file = feature.ds) # Read phenotype file (CLS format)
#	cls.labels <- CLS$class.v
#	cls.phen <- CLS$phen
#	cls.list <- gsub("\\s+", "", CLS$class.list, perl=TRUE)  # get rid of trailing whitespace   
	
	
#	if (is.vector(cls.labels)) {
#		n.phen <- 1
#	} else {
#		n.phen <- length(cls.labels[,1])
#	}
#	if (!is.na(user.colors[1])) {
#		c.test <- user.colors
#	} else {
#		if (!is.null(CLS$col.phen)) {
#			c.test <- CLS$col.phen
#		} else {
#			c.test <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
#					brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
#					brewer.pal(n=8, name="BuGn"),
#					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
#					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
#					brewer.pal(n=8, name="BuGn"),
#					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
#					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
#					brewer.pal(n=8, name="BuGn"))
#		}
#	}
	
	
#	if (!is.null(CLS$phen.names)) {
#		phen.names <- CLS$phen.names
#	} else {
#		phen.names <- "NA"
#	}
#	
#	
#	cls.phen.index <- unlist(cls.phen)
#	cls.phen.colors <- c.test[1:length(cls.phen.index)]
	##	print("cls.phen.colors:")
	##	print(cls.phen.colors)
#	
#	n.classes <- vector(length=n.phen, mode="numeric")
#	if (n.phen == 1) {
#		max.classes <- length(cls.phen)
#		n.classes[1] <- max.classes
#	} else {
#		max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
#		n.classes = unlist(lapply(cls.phen, length))
	##		for (i in 1:n.phen) {
	##			n.classes[i] <- length(cls.phen[[i]])
	##		}
#	}
	
	
	
	## Get rid of trailing whitespace
#	browser()
#	cls.list = rbind(rep("WT", length=ifelse(is.vector(cls.list), length(cls.list), length(cls.list[1,]))), 
#			apply(cls.list, MARGIN=2, FUN=function(x) gsub("\\s+$", "", x, perl=TRUE)))
#	cls.labels = rbind(rep(1, length=ifelse(is.vector(cls.list), length(cls.list), length(cls.list[1,]))), 
#			cls.labels)
	##			gsub("\\s+$", "", cls.labels, #unlist(cls.phen2.pass1), 
	##					perl=TRUE))
#	
	##	browser()
	
	## Remove "Commented out" gene names (with # at beginning of name)
	
#	browser()
	
#	if(kde2d.timing){
#		kde2d.t0 <<- kde2d.t1 <<- kde2d.t2 <<- kde2d.t3 <<- after.kde2d <<- 
#				before.kde2d <<- vector(length=(12*feature.ds.original$nRow), mode = "numeric")
#		kde2d.i <<- 1
#	}
	
#	browser()
	if(target.remove.na.samples){
		target.ds = REVEALER.remove.na.sample(target.ds)
		if( !is.na(target.ds$removed.na.sample.index[1] )){
			## Add removed column, its name, and index so it can be reinserted (just in case)
			feature.ds.original$removed.na.sample = feature.ds.original$matrix[,target.ds$removed.na.sample.index]
			feature.ds.original$removed.na.sample.index = target.ds$removed.na.sample.index
			feature.ds.original$removed.na.sample.name = target.ds$removed.na.sample.name
			
			## Adjust matrix/sample.names/nSample
			feature.ds.original$matrix = feature.ds.original$matrix[,-target.ds$removed.na.sample.index]
			feature.ds.original$sample.names = colnames(feature.ds.original$matrix)
			feature.ds.original$nSample = dim(feature.ds.original$matrix)[2]
#	print(feature.ds)
		}
	}
	
	if(target.remove.na.targets){
		target.ds = REVEALER.remove.na.rows(target.ds)
		if(target.ds$nRow == 0 ){
			stop(paste("removing all NA rows from target.ds resulted in removing all target rows",
							"- please revise the data and remove NAs"))
			
		}
	}
	if(feature.remove.na.features){
		feature.ds.original = REVEALER.remove.na.rows(feature.ds.original)
		if(feature.ds.original$nRow == 0 ){
			stop(paste("removing all NA rows from feature.ds resulted in removing all target rows",
							"- please revise the data and remove NAs"))
			
		}
	}
	
	if( !is.null(feature.annotations.file)){
		feature.ds.original = REVEALER.add.feature.annotation(feature.ds.original, feature.seed.names.to.match,
				feature.annotations.file)   ## this function still needs testing!
	} else{
		feature.ds.original$feature.seed.names.to.match = feature.seed.names.to.match
	}
	
	## Remove feature's na cols
	if(feature.remove.na.samples){
		feature.ds.original = REVEALER.remove.na.sample(feature.ds.original)
		if(!is.na(feature.ds.original$removed.na.sample.index[1])){
			
			## Add removed column, its name, and index so it can be reinserted (just in case)
			target.ds$removed.na.sample.index = feature.ds.original$removed.na.sample.index
			target.ds$removed.na.sample = target.ds$matrix[,target.ds$removed.na.sample.index]
			target.ds$removed.na.sample.name = target.ds$sample.names[target.ds$removed.na.sample.index]
			
			## Adjust matrix/sample.names/nSample
			target.ds$matrix = target.ds$matrix[,-target.ds$removed.na.sample.index]
			target.ds$sample.names = target.ds$sample.names[-target.ds$removed.na.sample.index]
			target.ds$nSample = length(target.ds$sample.names)
			
			if(target.ds$nRow == 1 ){
				target.ds$matrix = t(as.matrix(target.ds$matrix))
			}
		}
	}
	
	if(sum(is.na(target.ds$matrix)) > 0 ){
		stop("target dataset has NAs")
	}
	if(sum(is.na(feature.ds.original$matrix)) > 0 ){
		stop("feature dataset has NAs")
	}
	
	
#	kde2d.t0 <<- kde2d.t1 <<- kde2d.t2 <<- kde2d.t3 <<- 
#			after.kde2d <<- 
#			before.kde2d <<- vector(length=(12*feature.ds.original$nRow), mode = "numeric")
#	kde2d.i <<- 1
	
	char.res <-  0.013 * target.ds$nRow + 0.65
	target.ds = REVEALER.normalize.scores(target.ds, normalize.score, normalization.type)
	pdf.options(height=pdf.height, width=pdf.width, colormodel="rgb", bg="transparent")
	if( !de.novo.search || align.only.no.discovery ){
#		browser()
		feature.seed.names.to.match = feature.ds.original$feature.seed.names.to.match
		feature.ds = REVEALER.extract.row.names.to.match(feature.ds.original, feature.seed.names.to.match)
		
		#browser()
		preprocessed.target.and.feature = REVEALER.preprocessing(
				target.ds, 
				feature.ds,
				make.inbetween.heatmaps, 
				show.multiple.tissues,
				pdf.height, 
				pdf.width, 
				results.dir, 
				extracted.file.prefix, 
				paste(file.suffix, ".initial.calculations", sep=""),
				decreasing.order,
				tissue.names,
				tissue.labels,
				if.feature.summary.on.top,
				kde2d.timing = kde2d.timing,
				target.if.unmatched = target.if.unmatched,
				draw.white.lines = draw.white.lines)
		target.ds = preprocessed.target.and.feature$target.ds
		feature.ds = preprocessed.target.and.feature$feature.ds
		
#		browser()
		### Begin Pass 2 ###
		
		if( #n.phen.pass1 > 2 || 
				!feature.skip.refinement ){
			print( "--- Begin Step 2 ---")
#			if( is.na(MI.list.pass1[1]) || is.na(MI.list.phen.pass1[1]) ){
#				dev.off()
#				stop("Genomic abnormalities specified are not in the cls file of the samples provided.")
			##				return()
#			}
			#browser()	
			MI.thresholds = ifelse(match.to.blue, -1, 1) * c(0.2, 0.1, 0.08, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01)
			#	MI.threshold = 0.03
			ind.top.feature.MI = vector(mode="integer")
			MI.i = 0
			
			while( length(ind.top.feature.MI) < 1){
				MI.i = MI.i + 1
				if( MI.i > length(MI.thresholds)){
#					dev.off()
					print(paste("Selected genomic aberrations do not have",
									ifelse(match.to.blue, "negative", "positive"), 
									"mutual information with a low enough false",
									"discovery rate with the selected pathways."))
					print("If you'd like to do REVEALER anyway, set feature.skip.refinement=TRUE")
					print(paste("---- end", tissue, "----"))
					return()
				}
				if(match.to.blue){
					ind.top.feature.MI = which( feature.ds$MI.features <= -MI.thresholds[MI.i])
				} else{
					ind.top.feature.MI = which( feature.ds$MI.features >= MI.thresholds[MI.i]) #& MI.list.phen.pass1[-1] > 0 
				}
			}
			feature.ds$matrix = feature.ds$matrix[ind.top.feature.MI,]
			if(is.vector(feature.ds$matrix)){
				feature.ds$matrix = t(as.matrix(feature.ds$matrix))
				feature.ds$summary = feature.ds$matrix
			} else{
				feature.ds$summary = colSums(feature.ds$matrix)
				
			}
			feature.ds$row.names = feature.ds$row.names[ind.top.feature.MI]
			feature.ds$nRow = length(feature.ds$row.names)
			feature.ds$MI.features = feature.ds$MI.features[ind.top.feature.MI]
			
#			browser()
			preprocessed.target.and.feature = REVEALER.preprocessing(
					target.ds, 
					feature.ds,
					make.inbetween.heatmaps, 
					show.multiple.tissues,
					pdf.height, 
					pdf.width, 
					results.dir, 
					extracted.file.prefix, 
					paste(file.suffix, ".refinement", sep=""),
					decreasing.order,
					tissue.names,
					tissue.labels,
					if.feature.summary.on.top,
					kde2d.timing = kde2d.timing,
					target.if.unmatched = target.if.unmatched,
					draw.white.lines = draw.white.lines)
			target.ds = preprocessed.target.and.feature$target.ds
			feature.ds = preprocessed.target.and.feature$feature.ds
		} 
	} else{ print("'de.novo.search' on -- skipping preprocessing and simply 'filling in' from scratch")}
	### 3rd Pass ###	
	if(!align.only.no.discovery #|| de.novo.search
			){
		print( "--- Begin Pass 3 (Iterative Method)---")
		print("2 in explained vector = previous explained vector   1 = new additions")
		gc()
		#browser()
		initial.feature.string = ifelse(de.novo.search, "de.novo", 
				paste(c(feature.ds$row.names), collapse="."))
#		browser()
		if(!popup.heatmap){
			pdf(file=paste(results.dir, extracted.file.prefix, file.suffix, ".",
							initial.feature.string, ".",
							target.ds$row.names[1],
							".REVEALER.n.iter=", n.iter,
							".pdf", sep=""), 
					height=8.5, width=11)
		}
		
		#browser()
		if(de.novo.search){
			sample.order = order(target.ds$matrix[1,], decreasing = !match.to.blue)
#			if(target.ds$nRow == 1){
			target.ds$matrix = t(as.matrix(target.ds$matrix[1,sample.order]))
#			} else{
#				target.ds$matrix = target.ds$matrix[1,sample.order]
#			}
			target.ds$sample.names = target.ds$sample.names[sample.order]
			
			#sample.order = match(feature.ds.original$sample.names, target.ds$sample.names)
			feature.ds.original$matrix = feature.ds.original$matrix[,sample.order]
			feature.ds.original$sample.names = target.ds$sample.names#feature.ds.original$sample.names[sample.order]
		} else{
			sample.order = match(feature.ds$sample.names, feature.ds.original$sample.names)
			feature.ds.original$matrix = feature.ds.original$matrix[,sample.order]
			feature.ds.original$sample.names = feature.ds.original$sample.names[sample.order]
		}
		target.median = median(target.ds$matrix[1,])
		target.midpoint.index = REVEALER.find.median.index(target.ds$matrix[1,], match.to.blue)
		tissue.string = ifelse(tissue == "NA", "", tissue)
		if( de.novo.search ){
			file.suffix = paste(file.suffix, ".de.novo.search.", sep="")
			de.novo.top.feature = REVEALER.feature.selector(
					target.ds,
#				feature.ds.original,
					feature.ds.original,
					match.to.blue,
					decreasing.order,
					ifNaive = TRUE,
#					iterations.remaining = NULL,   ## Number of iterations left, for recursion. NULL if ifNaive==TRUE
					nTopFeatures = 20,
#					selected.features = NULL,   ## NULL if ifNaive==TRUE
					explained = NULL,
					unexplained.samples.remaining = NULL,
					results.dir = results.dir,
					extracted.file.prefix = extracted.file.prefix,
					file.suffix = file.suffix,
					target.median = target.median,
					target.midpoint.index = target.midpoint.index,
					popup.heatmap = popup.heatmap,
					tissue = tissue.string,
					kde2d.timing = kde2d.timing,
					plot.unproductive.iteration = plot.unproductive.iteration,
					draw.white.lines = draw.white.lines)
			#browser()
			
			explained = ifelse(de.novo.top.feature$explained >0, 1, 0)
			feature.ds.original = REVEALER.exclude.rows(feature.ds.original, 
					de.novo.top.feature$name)
			explained.original = list(explained=explained, MI=de.novo.top.feature$MI.explained)
#			browser()
			seed.ds = list(
					matrix=t(as.matrix(explained)), 
					row.names=de.novo.top.feature$name,
					nRow = 1,
					sample.names = target.ds$sample.names,
					nSample = target.ds$nSample,
					MI.features = de.novo.top.feature$MI.vector,
					MI.summary = de.novo.top.feature$MI.vector,
					summary = t(as.matrix(explained)))
			colnames(seed.ds$matrix) = seed.ds$sample.names
			rownames(seed.ds$matrix) = seed.ds$row.names
		} else{
#			if(!is.na(feature.ds.original$removed.na.sample.index[1])){
#				explained = feature.ds$summary[-feature.ds.original$removed.na.sample.index]
#			} else{
			explained = feature.ds$summary = ifelse(feature.ds$summary>0, 1, 0)
#			}
			explained.original = list(explained=explained, MI=feature.ds$MI.summary)
			feature.ds.original = REVEALER.exclude.rows(feature.ds.original, 
					feature.ds$row.names)
			seed.ds = feature.ds
		}
		#browser()
		REVEALER.feature.selector.heatmap(target.ds,
				target.midpoint.index,
				seed.ds,
				match.to.blue,
				target.main.text = paste(tissue),
				feature.main.text = paste("target=", target.ds$row.names[1], 
						"   seed=", paste(seed.ds$row.names, collapse=",", sep="")),
				initial.heatmap=TRUE,
				draw.white.lines = draw.white.lines
		)
#				nTopFeatures = nTopFeatures,
#				row.label.text.resize = row.label.text.resize,
#				col.label.text.resize = col.label.text.resize,
#				unexplained.samples.remaining = unexplained.samples.remaining
#			iteration#,
#			results.dir,
#			extracted.file.prefix,
#			file.suffix,
#			initial.feature.string
		### --- Write samples in order of target, for printing and reference ---
		write.table(matrix(c(
								c(paste("Rank_from_Left=1_to_Right=", target.ds$nSample, sep=""), 
										1:target.ds$nSample),
								c("Sample_Name", target.ds$sample.names), 
								c(paste("Status_at_", paste(seed.ds$row.names, collapse="_"), sep=""), 
										seed.ds$matrix[1,])),
						ncol=3), 
				file = paste(results.dir, extracted.file.prefix, file.suffix, ".REVEALER.samples.and.seed.",
						target.ds$row.names[1], ".",
						paste(seed.ds$row.names, collapse="_"),
						".txt", sep=""),
				quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE, append=FALSE)
		### --- End writing of samples ---
		
		
		
		
		### ---- Feature Selection ----
#		browser()
		selected.feature = vector(length=n.iter, mode="list")
		this.iter = 0
		MI.difference = ifelse(match.to.blue, -1, 1) ## Difference between previous "explained" vector and adding 
		##												this new feature
#		explained.original = list(explained=explained, MI=
		unexplained.samples.remaining = 
				REVEALER.get.unexplained.samples.remaining(
						explained=explained, 
						target=target.ds$matrix[1,], 
						target.median=target.median, 
						match.to.blue=match.to.blue)
		this.iter.index = this.iter + 1
		while( this.iter < n.iter && 
				ifelse(match.to.blue, MI.difference < 0, MI.difference > 0) &&
				sum(unexplained.samples.remaining[1:target.midpoint.index]) > 0
				){
			
			print(paste("Iteration:", this.iter.index))
#			if(this.iter.index == 10) browser()
			selected.feature[[this.iter.index]] = REVEALER.feature.selector(
					target.ds = target.ds,
					feature.ds = feature.ds.original,
					match.to.blue = match.to.blue,
					decreasing.order = decreasing.order,
					ifNaive = FALSE,
					nTopFeatures = 20,
					explained.samples = explained,
					unexplained.samples.remaining = unexplained.samples.remaining,
					target.median = target.median,
					target.midpoint.index = target.midpoint.index,
					iteration = this.iter.index,
					results.dir = results.dir,
					extracted.file.prefix = extracted.file.prefix,
					file.suffix = file.suffix,
					initial.feature.string = paste(seed.ds$row.names, collapse = "  "),
					popup.heatmap = popup.heatmap,
					tissue = tissue.string,
					kde2d.timing = kde2d.timing,
					draw.white.lines = draw.white.lines)
			
			if(selected.feature[[this.iter.index]][1] == "NULL"){
				break
			}
			
			# Remove the feature we just found from the search space
			feature.ds.original = REVEALER.exclude.rows(feature.ds.original, 
					selected.feature[[this.iter.index]]$name)
			
			MI.difference = selected.feature[[this.iter.index]]$MI.difference
			unexplained.samples.remaining = selected.feature[[this.iter.index]]$unexplained
			explained = selected.feature[[this.iter.index]]$explained
			this.iter = this.iter + 1
			this.iter.index = this.iter + 1
		}
		print("which criteria caused the 'while' loop to end:")
		if( !(this.iter < n.iter)){
			print("completed iterations")
			this.iter.index = n.iter
#			this.iter = n.iter
		} else if( 
				selected.feature[[this.iter.index]] == "NULL" ||
				!(ifelse(match.to.blue, MI.difference < 0, MI.difference > 0))){
			print("MI did not improve with the addition of the best feature found")
#			if(this.iter == 1) {
#				this.iter = this.iter - 1
#			}
#			print(paste("MI.difference = ", MI.difference))
		} else if( !(sum(unexplained.samples.remaining[1:target.midpoint.index]) > 0)){
			print("All explainable samples explained")
		}
#		browser()
#		if( this.iter != n.iter && (ifelse(match.to.blue, MI.difference < 0, MI.difference >0) ||
#				sum(unexplained.samples.remaining) > target.midpoint.index)){
#			## Cut down "selected.feature" to only be the iterations that improved
#			## Mutual information, unless the one iteration performed did not improve MI, then
#			## keep it just for good measure
#			selected.feature = selected.feature[[1:ifelse(this.iter==2, 1, (this.iter-1))]]
#		}
		
#		browser()
		REVEALER.plot.final.results(selected.feature, 
				seed.ds, 
				target.ds, 
				explained.original,
				match.to.blue = match.to.blue,
				target.midpoint.index = target.midpoint.index,
				n.productive.iter = this.iter, #ifelse(this.iter==2, 1, (this.iter-2)),  ## "productive" because they helped find features
				popup.heatmap = popup.heatmap,
				tissue = tissue.string,
				draw.white.lines = draw.white.lines)
		if(!popup.heatmap) dev.off()
		print(proc.time()-t1)
		print(date())
		### --- End Feature Selection ----
	}
	
	if (!is.na(output.dataset)) {
#		V.GCT <- m.all
#		print("Figure out why model.descs2.pass1 does not correspond to the rows of m.all")
#		browser()
#		colnames(V.GCT) <- sample.names2
#		row.names(V.GCT) <- model.names.pass1
		write.gct(gct.data.frame = m, descs = model.descs2.pass1, filename =output.dataset)  
	}
	
}

REVEALER.kde2d.timing <- function(
		###  "target" = typically continuous readout such as expression, gene set enrichment, drug/RNAi 
		### 		    sensitivity, degree of phosphorylation .... some values that you would like to
		###				"explain" with the features
		target.ds.file,
		target.ds.file.type = c("gct", "cls", "txt.sample.rows", "txt.sample.cols"), 
		## 											"gct" is an expression file
		##											"cls" is a class/phenotype file
		##											"txt.sample.rows" is a tab-delimited file with sample names on the row
		##												(i.e. first column is sample names)
		##											"txt.sample.cols" is a tab-delimited file with sample names on the columns
		##												(i.e. first row is sample names)
		target.names.to.match = c("NA", "ALL"),
		target.remove.na.samples = FALSE, ## target columns are the samples, ie the things that have a value for a target
		target.remove.na.targets = FALSE, ## target rows are the targets, ie the thing(s) you're trying to explain
		target.if.unmatched = "NA",  ## if no samples have the features specified, which target name to sort by
		
		### "feature" = continuous or binary data. Could be genomic such as mutation, copy number, 
		### 			or epigenomic such as methylation, phosphorylation status, etc
		feature.ds.file,
		feature.ds.file.type = c("gct", "cls", "txt.sample.rows", "txt.sample.cols"),
		feature.remove.na.samples = FALSE, ## feature columns are samples, ie the things that have the features
		feature.remove.na.features = FALSE, ## feature rows are features, ie the things you're searching for
		## 											"gct" is an expression file
		##											"cls" is a class/phenotype file
		##											"txt.sample.rows" is a tab-delimited file with sample names on the row
		##												(i.e. first column is sample names)
		##											"txt.sample.cols" is a tab-delimited file with sample names on the columns
		##												(i.e. first row is sample names)
		feature.seed.names.to.match = "NA",
		feature.suffix.or = c("", "_MUT", "mut"),	##  OR suffix. Suffixes one of "", "_MUT" or "mut", or whatever the 
		## 											user specified to feature names. 
		## 											Can be helpful if you just want to specify 
		## 											a gene list but your dataset says "GENE_MUT"
		## 											Setting to "" suffixes nothing (retains original naming)
		feature.suffix.and = c("_AMP", "_DEL"),	##  AND suffix. Suffixes all of "_AMP" and "_DEL" (or whatever 
		##											the user specified) to the feature names. Set to boolean FALSE to 
		##											retain original naming. Specifically suffixes separately
		##											from the OR suffix (as in, you will have "GENE_MUT" and "GENE_AMP"
		##											separately, but no "GENE_MUT_AMP")
		feature.add.chrom.locs = FALSE, ##  if feature names are in form "GENE_MUT" or "GENE" then the 
		## 									chromosome location can be looked up using the BioconductoR
		## 									library "org.Hs.eg.db" If features are not in a known gene id format
		##									that can be recognized, will be "NA." Can be helpful to keep track of genes
		feature.skip.refinement = FALSE, ## FALSE: You want to use all the features specified in "feature.seed.names.to.match" and
		##									don't want to do any refinement of your features
		##								TRUE: You're not sure if all your features are relevant so you just want 
		##									to pick and use the top ones for feature selection
		feature.remove.from.selection.search.space = NULL, ## Remove from "feature.ds.original" and don't consider them when doing
		##											REVEALER.feature.selector
		match.to.blue = FALSE,   ## Match to the negative side of the target
		de.novo.search = FALSE, ## Ignore feature.seed.names.to.match and seed with the best feature from feature.ds.file
		## ("best" feature meaning has highest normalized mutual information with target)
		align.only.no.discovery = FALSE,
		target.phen = NA,
		target.class = NA,
		
		## Add these:
#		target.remove.na.samples = TRUE,   # Samples that have NA in target will be removed in feature
#		target.remove.na.targets = FALSE,  # Independent of feature
#		feature.remove.na.samples = FALSE, # Samples that have NA in feature will be removed in target
#		feature.remove.na.features = TRUE,  # Independent of target
		
		
		
		
		### Filename and results parameters
		results.dir = NULL,   ## Where do you want to find your pdfs and txt output?
		tissue = "NA",  ## For pdf naming and figure titling
		file.suffix = "",
		
		
		normalize.score = T,
		normalization.type = "zero.one",
		
		
		
		
		#n.random.target.names.to.match = 10,
		
		
		### Iterations of REVEALER, or until MI starts to decrease (or increase
		### if match.to.blue == TRUE)
		n.iter = 5,
		
		### Plotting parameters
		popup.heatmap = FALSE,  # pops up a heatmap window instead of saving to pdf (mac only) 
		pdf.height = 11,
		pdf.width = 17,
		highlight.tissue.name = NULL,
		show.multiple.tissues = FALSE,
		ifLegend = FALSE,
		ifScoreBarPlot = FALSE,
		if.feature.summary.on.top = FALSE,
		user.colors = NA,
#		decreasing.order = T,
		output.dataset = NA,
		char.rescale = 1,
		cmap.type = 3,
		row.norm = T,
		
		#do.mRMR = FALSE,
		#skip.step1 = FALSE,
		
#		continuous.data = FALSE,
		
		### Debugging params
		make.inbetween.heatmaps = FALSE,  ## for debugging and preparing intermediate steps for presentations
		exclude.samples = NA     ## For cleaning up and debugging datasets

)
# Calls MSIG.HeatMapPlot.9 and makes a plot sorted by the highest-scoring
# target.names.to.match and abnormalities (gene mutations or copy number alterations)
# i.e. doesn't require a "model" to score by as OPAM.sort.projection.by.score.2 does.
# However, it *will* use "model" if it cannot calculate p-values on the gene target.names.to.match, which
# happens when every cell line has a genomic aberration.
#
# Runs 3 passes on the data:
# 1st pass: looks at the genes and copy number alterations specified by feature.seed.names.to.match
# 2nd pass: looks at only the top abnormalities (using a p-value cutoff) from the 1st pass, and adjusts 
# the PATHWAY.MUT vector accordingly (only according to the genes, not by copy number data)
# 3rd pass: Takes the winning signature from the 2nd pass and then looks all the genes available
#
# Very similar to OPAM.sort.projection.by.score.4, however this version uses mutual.inf instead of
# roc.area to calculate mutual information scores and p-values for PATHWAY.MUT, the vector of total genomic aberrations
# in all samples
#
# Differs from OPAM.sort.projection.by.score.6 by requiring the gct file of expression in 
# all pathways by the input tissue ("target.ds")
{
	
	REVEALER.load.libraries(c("gtools", "verification", "ROCR", "MASS", "RColorBrewer", "heatmap.plus"))
	
	decreasing.order = ifelse(match.to.blue, FALSE, TRUE)
	if( feature.add.chrom.locs ) file.suffix = paste(file.suffix, "_with.chrom.locs", sep="")
	
	ds.parser.lookup = c(REVEALER.parse.gct, REVEALER.parse.cls, REVEALER.parse.txt.sample.rows,
			REVEALER.parse.txt.sample.cols)
	names(ds.parser.lookup) = c("gct", "cls", "txt.sample.rows", "txt.sample.cols")
	
	target.parser = eval(substitute(ds.parser.lookup$file.type, list(file.type = as.name(target.ds.file.type))))
	target.ds.original <- target.parser(filename = target.ds.file)  # "ds" for "dataset"
	
	feature.parser = eval(substitute(ds.parser.lookup$file.type, list(file.type = as.name(feature.ds.file.type))))
	feature.ds.original = feature.parser(filename = feature.ds.file)  # "ds" for "dataset"
	
#	browser()
	
	if(!is.na(exclude.samples)){
		target.ds.original = REVEALER.exclude.samples(target.ds.original, exclude.samples)
		feature.ds.original = REVEALER.exclude.samples(feature.ds.original, exclude.samples)		
	}
	if(!is.null(feature.remove.from.selection.search.space)){
		feature.ds.original = REVEALER.exclude.rows(feature.ds.original, 
				feature.remove.from.selection.search.space)
	}
	
	tissues = REVEALER.check.multiple.tissues(target.ds.original, show.multiple.tissues, tissue)
	results.dir = REVEALER.fix.results.dir(results.dir)
	extracted.file.prefix <- REVEALER.get.extracted.file.prefix(target.ds.file, target.ds.file.type)
	if(target.names.to.match[1] != "ALL" && target.names.to.match!="NA"){
		target.ds = REVEALER.extract.row.names.to.match(target.ds.original, target.names.to.match)
	}else if(target.names.to.match == "ALL"){
		target.ds = target.ds.original
	} else{
		stop("Must provide target names to match")
	}
	
	
#	CLS <- MSIG.ReadPhenFile.2(file = feature.ds) # Read phenotype file (CLS format)
#	cls.labels <- CLS$class.v
#	cls.phen <- CLS$phen
#	cls.list <- gsub("\\s+", "", CLS$class.list, perl=TRUE)  # get rid of trailing whitespace   
	
	
#	if (is.vector(cls.labels)) {
#		n.phen <- 1
#	} else {
#		n.phen <- length(cls.labels[,1])
#	}
#	if (!is.na(user.colors[1])) {
#		c.test <- user.colors
#	} else {
#		if (!is.null(CLS$col.phen)) {
#			c.test <- CLS$col.phen
#		} else {
#			c.test <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
#					brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
#					brewer.pal(n=8, name="BuGn"),
#					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
#					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
#					brewer.pal(n=8, name="BuGn"),
#					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
#					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
#					brewer.pal(n=8, name="BuGn"))
#		}
#	}
	
	
#	if (!is.null(CLS$phen.names)) {
#		phen.names <- CLS$phen.names
#	} else {
#		phen.names <- "NA"
#	}
#	
#	
#	cls.phen.index <- unlist(cls.phen)
#	cls.phen.colors <- c.test[1:length(cls.phen.index)]
	##	print("cls.phen.colors:")
	##	print(cls.phen.colors)
#	
#	n.classes <- vector(length=n.phen, mode="numeric")
#	if (n.phen == 1) {
#		max.classes <- length(cls.phen)
#		n.classes[1] <- max.classes
#	} else {
#		max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
#		n.classes = unlist(lapply(cls.phen, length))
	##		for (i in 1:n.phen) {
	##			n.classes[i] <- length(cls.phen[[i]])
	##		}
#	}
	
	
	
	## Get rid of trailing whitespace
#	browser()
#	cls.list = rbind(rep("WT", length=ifelse(is.vector(cls.list), length(cls.list), length(cls.list[1,]))), 
#			apply(cls.list, MARGIN=2, FUN=function(x) gsub("\\s+$", "", x, perl=TRUE)))
#	cls.labels = rbind(rep(1, length=ifelse(is.vector(cls.list), length(cls.list), length(cls.list[1,]))), 
#			cls.labels)
	##			gsub("\\s+$", "", cls.labels, #unlist(cls.phen2.pass1), 
	##					perl=TRUE))
#	
	##	browser()
	
	## Remove "Commented out" gene names (with # at beginning of name)
	
#	browser()
	
	if(target.remove.na.samples){
		target.ds = REVEALER.remove.na.sample(target.ds)
		if( !is.na(target.ds$removed.na.sample.index[1] )){
			feature.ds.original$removed.na.sample = feature.ds.original$matrix[,target.ds$removed.na.sample.index]
			feature.ds.original$removed.na.sample.index = target.ds$removed.na.sample.index
			feature.ds.original$removed.na.sample.name = target.ds$removed.na.sample.name
			feature.ds.original$matrix = feature.ds.original$matrix[,-target.ds$removed.na.sample.index]
			feature.ds.original$sample.names = colnames(feature.ds.original$matrix)
			feature.ds.original$nSample = length(feature.ds$sample.names)#dim(feature.ds.original$matrix)[2]
			
#	print(feature.ds)
		}
	}
	
	if(target.remove.na.targets){
		target.ds = REVEALER.remove.na.rows(target.ds)
	}
	if(feature.remove.na.features){
		feature.ds.original = REVEALER.remove.na.rows(feature.ds.original)
	}
	
#	target.ds = REVEALER.remove.na.rows(target.ds)
	
	feature.ds.original = REVEALER.add.feature.suffixes(feature.ds.original, feature.seed.names.to.match,
			feature.suffix.or, feature.suffix.and, feature.add.chrom.locs)   ## this function still needs testing!
	
	## Remove feature's na cols
	if(feature.remove.na.samples){
		feature.ds.original = REVEALER.remove.na.sample(feature.ds.original)
		if(!is.na(feature.ds.original$removed.na.sample.index[1])){
			
			
			## Add removed column, its name, and index so it can be reinserted (just in case)
			target.ds$removed.na.sample.index = feature.ds.original$removed.na.sample.index
			target.ds$removed.na.sample = target.ds$matrix[,target.ds$removed.na.sample.index]
			target.ds$removed.na.sample.name = target.ds$sample.names[target.ds$removed.na.sample.index]
			
			## Adjust matrix/sample.names
			target.ds$matrix = target.ds$matrix[,-target.ds$removed.na.sample.index]
			target.ds$sample.names = target.ds$sample.names[-target.ds$removed.na.sample.index]
			target.ds$nSample = length(target.ds$sample.names)
			
			if(target.ds$nRow == 1 ){
				target.ds$matrix = t(as.matrix(target.ds$matrix))
			}
		}
	}
	
	feature.seed.names.to.match = feature.ds.original$feature.seed.names.to.match
	feature.ds = REVEALER.extract.row.names.to.match(feature.ds.original, feature.seed.names.to.match)
	kde2d.t0 <<- kde2d.t1 <<- kde2d.t2 <<- kde2d.t3 <<- after.kde2d <<- 
			before.kde2d <<- vector(length=(12*feature.ds.original$nRow), mode = "numeric")
	kde2d.i <<- 1
	
	char.res <-  0.013 * target.ds$nRow + 0.65
	target.ds = REVEALER.normalize.scores(target.ds, normalize.score, normalization.type)
	pdf.options(height=pdf.height, width=pdf.width, colormodel="rgb", bg="transparent")
	if( !de.novo.search || align.only.no.discovery ){
		preprocessed.target.and.feature = REVEALER.preprocessing(
				target.ds, 
				feature.ds,
				make.inbetween.heatmaps, 
				show.multiple.tissues,
				pdf.height, 
				pdf.width, 
				results.dir, 
				extracted.file.prefix, 
				paste(file.suffix, ".initial.calculations", sep=""),
				decreasing.order,
				tissue.names,
				tissue.labels,
				if.feature.summary.on.top)
		target.ds = preprocessed.target.and.feature$target.ds
		feature.ds = preprocessed.target.and.feature$feature.ds
		
#		browser()
		### Begin Pass 2 ###
		
		if( #n.phen.pass1 > 2 || 
				!feature.skip.refinement ){
			print( "--- Begin Step 2 ---")
#			if( is.na(MI.list.pass1[1]) || is.na(MI.list.phen.pass1[1]) ){
#				dev.off()
#				stop("Genomic abnormalities specified are not in the cls file of the samples provided.")
			##				return()
#			}
			#browser()	
			MI.thresholds = ifelse(match.to.blue, -1, 1) * c(0.2, 0.1, 0.08, 0.05, 0.03, 0.025, 0.02, 0.015, 0.01, 0)
			#	MI.threshold = 0.03
			ind.top.feature.MI = vector(mode="integer")
			MI.i = 0
			
			while( length(ind.top.feature.MI) < 1){
				MI.i = MI.i + 1
				if( MI.i > length(MI.thresholds)){
#					dev.off()
					print(paste("Selected genomic aberrations do not have",
									ifelse(match.to.blue, "negative", "positive"), 
									"mutual information with a low enough false",
									"discovery rate with the selected pathways."))
					print("If you'd like to do REVEALER anyway, set feature.skip.refinement=TRUE")
					print(paste("---- end", tissue, "----"))
					return()
				}
				if(match.to.blue){
					ind.top.feature.MI = which( feature.ds$MI.features <= MI.thresholds[MI.i])
				} else{
					ind.top.feature.MI = which( feature.ds$MI.features >= MI.thresholds[MI.i]) #& MI.list.phen.pass1[-1] > 0 
				}
			}
			feature.ds$matrix = feature.ds$matrix[ind.top.feature.MI,]
			if(is.vector(feature.ds$matrix)){
				feature.ds$matrix = t(as.matrix(feature.ds$matrix))
				feature.ds$summary = feature.ds$matrix
			} else{
				feature.ds$summary = colSums(feature.ds$matrix)
				
			}
			feature.ds$row.names = feature.ds$row.names[ind.top.feature.MI]
			feature.ds$nRow = length(feature.ds$row.names)
			feature.ds$MI.features = feature.ds$MI.features[ind.top.feature.MI]
			
#			browser()
			preprocessed.target.and.feature = REVEALER.preprocessing(
					target.ds, 
					feature.ds,
					make.inbetween.heatmaps, 
					show.multiple.tissues,
					pdf.height, 
					pdf.width, 
					results.dir, 
					extracted.file.prefix, 
					paste(file.suffix, ".refinement", sep=""),
					decreasing.order,
					tissue.names,
					tissue.labels,
					if.feature.summary.on.top)
			target.ds = preprocessed.target.and.feature$target.ds
			feature.ds = preprocessed.target.and.feature$feature.ds
		} 
	} else{ print("'de.novo.search' on -- skipping preprocessing and simply 'filling in' from scratch")}
	### 3rd Pass ###	
	if(!align.only.no.discovery #|| de.novo.search
			){
		print( "--- Begin Pass 3 (Iterative Method)---")
		print("2 in explained vector = previous explained vector   1 = new additions")
		gc()
		initial.feature.string = ifelse(de.novo.search, "de.novo", 
				paste(c(feature.ds$row.names), collapse="."))
#		browser()
		if(!popup.heatmap){
			pdf(file=paste(results.dir, extracted.file.prefix, file.suffix, ".",
							initial.feature.string, ".",
							target.ds$row.names[1],
							".REVEALER.n.iter=", n.iter,
							".pdf", sep=""), 
					height=8.5, width=11)
		}
		
		#browser()
		if(de.novo.search){
			sample.order = order(target.ds$matrix[1,], decreasing = !match.to.blue)
#			if(target.ds$nRow == 1){
			target.ds$matrix = t(as.matrix(target.ds$matrix[1,sample.order]))
#			} else{
#				target.ds$matrix = target.ds$matrix[1,sample.order]
#			}
			target.ds$sample.names = target.ds$sample.names[sample.order]
			
			#sample.order = match(feature.ds.original$sample.names, target.ds$sample.names)
			feature.ds.original$matrix = feature.ds.original$matrix[,sample.order]
			feature.ds.original$sample.names = target.ds$sample.names#feature.ds.original$sample.names[sample.order]
		} else{
			sample.order = match(feature.ds$sample.names, feature.ds.original$sample.names)
			feature.ds.original$matrix = feature.ds.original$matrix[,sample.order]
			feature.ds.original$sample.names = feature.ds.original$sample.names[sample.order]
		}
		target.median = median(target.ds$matrix[1,])
		target.midpoint.index = REVEALER.find.median.index(target.ds$matrix[1,], match.to.blue)
		tissue.string = ifelse(tissue == "NA", "", tissue)
		
		de.novo.top.feature = REVEALER.feature.selector(
				target.ds,
#				feature.ds.original,
				feature.ds.original,
				match.to.blue,
				decreasing.order,
				ifNaive = TRUE,
#					iterations.remaining = NULL,   ## Number of iterations left, for recursion. NULL if ifNaive==TRUE
				nTopFeatures = 20,
#					selected.features = NULL,   ## NULL if ifNaive==TRUE
				explained = NULL,
				unexplained.samples.remaining = NULL,
				results.dir = results.dir,
				extracted.file.prefix = extracted.file.prefix,
				file.suffix = file.suffix,
				target.median = target.median,
				target.midpoint.index = target.midpoint.index,
				popup.heatmap = popup.heatmap,
				tissue = tissue.string)
		#browser()
		if( de.novo.search ){
			explained = de.novo.top.feature$explained
			feature.ds.original = REVEALER.exclude.rows(feature.ds.original, 
					de.novo.top.feature$name)
			explained.original = list(explained=explained, MI=de.novo.top.feature$MI.explained)
#			browser()
			seed.ds = list(
					matrix=t(as.matrix(explained)), 
					row.names=de.novo.top.feature$name,
					nRow = 1,
					sample.names = target.ds$sample.names,
					nSample = target.ds$nSample,
					MI.features = de.novo.top.feature$MI.vector)
			colnames(seed.ds$matrix) = seed.ds$sample.names
			rownames(seed.ds$matrix) = seed.ds$row.names
		} else{
#			if(!is.na(feature.ds.original$removed.na.sample.index[1])){
#				explained = feature.ds$summary[-feature.ds.original$removed.na.sample.index]
#			} else{
			explained = feature.ds$summary
#			}
			explained.original = list(explained=explained, MI=feature.ds$MI.summary)
			feature.ds.original = REVEALER.exclude.rows(feature.ds.original, 
					feature.ds$row.names)
			seed.ds = feature.ds
		}
		
		
		### ---- Feature Selection ----
#		browser()
		selected.feature = vector(length=n.iter, mode="list")
		this.iter = 0
		MI.difference = ifelse(match.to.blue, -1, 1) ## Difference between previous "explained" vector and adding 
		##												this new feature
#		explained.original = list(explained=explained, MI=
		unexplained.samples.remaining = 
				REVEALER.get.unexplained.samples.remaining(
						explained=explained, 
						target=target.ds$matrix[1,], 
						target.median=target.median, 
						match.to.blue=match.to.blue)
		while( this.iter < n.iter && 
				ifelse(match.to.blue, MI.difference < 0, MI.difference > 0) &&
				sum(unexplained.samples.remaining[1:target.midpoint.index]) > 0
				){
			this.iter.index = this.iter + 1
			print(paste("Iteration:", this.iter.index))
			selected.feature[[this.iter.index]] = REVEALER.feature.selector(
					target.ds = target.ds,
					feature.ds = feature.ds.original,
					match.to.blue = match.to.blue,
					decreasing.order = decreasing.order,
					ifNaive = FALSE,
					nTopFeatures = 20,
					explained.samples = explained,
					unexplained.samples.remaining = unexplained.samples.remaining,
					target.median = target.median,
					target.midpoint.index = target.midpoint.index,
					iteration = this.iter.index,
					results.dir = results.dir,
					extracted.file.prefix = extracted.file.prefix,
					file.suffix = file.suffix,
					initial.feature.string = paste(seed.ds$row.names, collapse = "  "),
					popup.heatmap = popup.heatmap,
					tissue = tissue.string)
			
			# Remove the feature we just found from the search space
			feature.ds.original = REVEALER.exclude.rows(feature.ds.original, 
					selected.feature[[this.iter.index]]$name)
			
			MI.difference = selected.feature[[this.iter.index]]$MI.difference
			unexplained.samples.remaining = selected.feature[[this.iter.index]]$unexplained
			explained = selected.feature[[this.iter.index]]$explained
			this.iter = this.iter + 1
		}
#		browser()
#		if( this.iter != n.iter && (ifelse(match.to.blue, MI.difference < 0, MI.difference >0) ||
#				sum(unexplained.samples.remaining) > target.midpoint.index)){
#			## Cut down "selected.feature" to only be the iterations that improved
#			## Mutual information, unless the one iteration performed did not improve MI, then
#			## keep it just for good measure
#			selected.feature = selected.feature[[1:ifelse(this.iter==2, 1, (this.iter-1))]]
#		}
		
#		browser()
		REVEALER.plot.final.results(selected.feature, 
				seed.ds, 
				target.ds, 
				explained.original,
				match.to.blue = match.to.blue,
				target.midpoint.index = target.midpoint.index,
				n.productive.iter = this.iter, #ifelse(this.iter==2, 1, (this.iter-2)),  ## "productive" because they helped find features
				popup.heatmap = popup.heatmap,
				tissue = tissue.string)
		if(!popup.heatmap) dev.off()
		print(proc.time()-t1)
		print(date())
		### --- End Feature Selection ----
	}
	
	if (!is.na(output.dataset)) {
#		V.GCT <- m.all
#		print("Figure out why model.descs2.pass1 does not correspond to the rows of m.all")
#		browser()
#		colnames(V.GCT) <- sample.names2
#		row.names(V.GCT) <- model.names.pass1
		write.gct(gct.data.frame = m, descs = model.descs2.pass1, filename =output.dataset)  
	}
	
}


OPAM.sort.projection.by.score.8 <- function(
#		input.ds,
#		signatures = "NA",
		input.all.pathways.ds,
		input.cls,
		tissue = "NA",
		results.dir,
		normalize.score = T,
#		normalization.type = "zero.one",
		model = "NA",
		target.phen = NA,
		target.class = NA,
		user.colors = NA,
		decreasing.order = T,
		output.dataset = NA,
		char.rescale = 1,
		cmap.type = 3,
		row.norm = T,
		u.gene.names.known = "NA",
		n.randomizations = 10
)
# Calls MSIG.HeatMapPlot.9 and makes a plot sorted by the highest-scoring
# signatures and abnormalities (gene mutations or copy number alterations)
# i.e. doesn't require a "model" to score by as OPAM.sort.projection.by.score.2 does.
# However, it *will* use "model" if it cannot calculate p-values on the gene signatures, which
# happens when every cell line has a genomic aberration.
#
# Runs 3 passes on the data:
# 1st pass: looks at the genes and copy number alterations specified by u.gene.names.known
# 2nd pass: looks at only the top abnormalities (using a p-value cutoff) from the 1st pass, and adjusts 
# the PATHWAY.MUT vector accordingly (only according to the genes, not by copy number data)
# 3rd pass: Takes the winning signature from the 2nd pass and then looks all the genes available
#
# Very similar to OPAM.sort.projection.by.score.4, however this version uses mutual.inf instead of
# roc.area to calculate mutual information scores and p-values for PATHWAY.MUT, 
# the vector of total genomic aberrations in all samples
#
# Differs from OPAM.sort.projection.by.score.6 by requiring the gct file of expression in 
# all pathways by the input tissue ("input.all.pathways.ds")
#
# Differs from OPAM.sort.projection.by.score.7 by finding the top enriched pathways that 
# differentiate according to phenotype
# from testing all the pathways in "input.all.pathways.ds." Does not require signatures to be known a priori.
{
	
	library(gtools)
	library(verification)
	library(ROCR)
	library(MASS)
	library(RColorBrewer)
	library(heatmap.plus)
	
#	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
#	m <- data.matrix(dataset$ds)
#	model.names <- dataset$row.names
	##	model.descs <- dataset$descs
#	Ns <- length(m[1,])
#	dim(m)
#	sample.names <- dataset$names
	
	dataset.all <- MSIG.Gct2Frame( filename = input.all.pathways.ds)
	m.all <- data.matrix(dataset.all$ds)#[1:30,]
	#model.names <- dataset.all$row.names#[1:30]
	model.names <- make.unique(dataset.all$descs)
	m.all <- na.omit(t(apply(m.all, MARGIN=1, FUN=normalize)))
	Ns = length(m.all[1,])
	sample.names = dataset.all$names
	
#	if( signatures == "NA" ){
#		stop("Must provide a vector of signature names to evaluate, or specify 'ALL'")
#	}
#	if( signatures == "ALL"){
#		model.names = model.names.all
#		m = m.all
#		model.descs = dataset.all$descs
#	} else{
#		model.names = signatures
#		model.ind = match(signatures, model.names.all)
#		m = m.all[model.ind,]
#		model.descs = dataset.all$descs[model.ind]
	##	browser()
#	}
#	model.names = model.names.all
#	m = m.all
	model.descs = dataset.all$descs#[1:30]
	target.ds$nRow <- length(m.all[,1])
	temp <- strsplit(input.all.pathways.ds, split="/") # Extract test file name
	s <- length(temp[[1]])
	test.file.name <- temp[[1]][s]
	temp <- strsplit(test.file.name, split=".gct")
	test.file.prefix <-  temp[[1]][1]
	char.res <-  0.013 * target.ds$nRow + 0.65
	
	# normalize scores
	
#	if (normalize.score == T) {
#		if (normalization.type == "zero.one") {
#			for (i in 1:target.ds$nRow) {
#				m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
#			}
#		} else if (normalization.type == "z.score") {
#			for (i in 1:target.ds$nRow) {
#				m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
#			}
#		} else if (normalization.type == "r.z.score") {
#			for (i in 1:target.ds$nRow) {
#				m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
#			}
#		}         
#	}
	
	CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
	cls.labels <- CLS$class.v
	cls.phen <- CLS$phen
	cls.list <- CLS$class.list
#	browser()
	if (is.vector(cls.labels)) {
		n.phen <- 1
	} else {
		n.phen <- length(cls.labels[,1])
	}
	if (!is.na(user.colors[1])) {
		c.test <- user.colors
	} else {
		if (!is.null(CLS$col.phen)) {
			c.test <- CLS$col.phen
		} else {
			c.test <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"),
					brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
					brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
					brewer.pal(n=8, name="BuGn"))
		}
	}
#	browser()
	
	if (!is.null(CLS$phen.names)) {
		phen.names <- CLS$phen.names
	} else if( !is.null(CLS$phen.list)){
		phen.names = CLS$phen.list
	} else {
		phen.names <- "NA"
	}
	
	cls.phen.index <- unlist(cls.phen)
	cls.phen.colors <- c.test[1:length(cls.phen.index)]
#	print("cls.phen.colors:")
#	print(cls.phen.colors)
	
	n.classes <- vector(length=n.phen, mode="numeric")
	if (n.phen == 1) {
		max.classes <- length(cls.phen)
		n.classes[1] <- max.classes
	} else {
		max.classes <- max(unlist(lapply(cls.phen, FUN=length)))
		for (i in 1:n.phen) {
			n.classes[i] <- length(cls.phen[[i]])
		}
	}
	print("--- Begin Pass 1 ---")
	phen.pass1 = u.gene.names.known
	n.phen.pass1 = length(u.gene.names.known)
	ind.phen.pass1 = which( phen.names %in% phen.pass1 )
	phen.pass1 = phen.names[ind.phen.pass1]
	
	MI.list.pass1 = vector( length=target.ds$nRow, mode="numeric" )
	FDR.list.pass1 = vector( length=target.ds$nRow, mode="numeric" )
	
	if( !is.vector(cls.labels)){
		cls.list.pass1 = cls.list[ind.phen.pass1,]
		cls.labels.pass1 = cls.labels[ind.phen.pass1,]
	} else{ 
		cls.list.pass1 = cls.list
		cls.labels.pass1 = cls.labels
	}
	cls.list.pass1.2 = t(as.matrix(ifelse(cls.list.pass1 == "WT", 0, 1)))
	
	if (!is.na(target.phen)) {
		if( is.vector(cls.list.pass1.2)){ bin.class.pass1 = cls.list.pass1.2 
		} else { bin.class.pass1 = colSums(cls.list.pass1.2) } #apply( cls.list.pass1.2, MARGIN=2, FUN=sum) }
		# Normalize bin.class.pass1
		if( length(unique(bin.class.pass1)) > 1){
			bin.class.pass1 = ( 
						bin.class.pass1 - min(bin.class.pass1))/(max(bin.class.pass1) - min(bin.class.pass1))
		} else if ( length(unique(bin.class.pass1)) == 1){
			bin.class = rep(1, length(cls.list[1,]))
		}
		if( is.vector( cls.list.pass1) ){
			cls.list.pass1 = ifelse(bin.class.pass1 > 0, "MUT", "WT")
		} else{	cls.list.pass1[1,] = ifelse(bin.class.pass1 > 0, "MUT", "WT") }
	} else {
#		browser()
		bin.class.pass1 <- ifelse(cls.list == cls.phen[1], 0, 1)
	}
	model.descs2.pass1 = vector(length = target.ds$nRow, mode="character")
	pdf(file=paste(tissue, test.file.name, ".Phase1", "pdf", sep="."))
	skipped.indices = 21:(target.ds$nRow)
	if( length(unique(bin.class.pass1)) > 1 ){
		MI.results = mutual.inf.3.v2(bin.class.pass1, m.all) #signature.indices = 1:target.ds$nRow, )
		MI.list.pass1  = MI.results$MI
		for (i in 1:target.ds$nRow) {
			MI.signif <- signif(MI.list.pass1[i], digits=3)
			model.descs2.pass1[i] <- paste(MI.signif, sep="")
		}
		m.order.pass1 = order(MI.list.pass1, decreasing=FALSE, na.last=TRUE)
		m.order.pass1.top10 = m.order.pass1[-skipped.indices]
		m2.pass1 <- m.all[m.order.pass1, ]
		s.order.pass1 <- order(m2.pass1[1,], decreasing = FALSE )
		m2.pass1 <- m2.pass1[-skipped.indices, s.order.pass1]
	} else{ 
		MI.list.pass1 = rep(NA, target.ds$nRow)
		model.descs2.pass1 = rep(" - ", target.ds$nRow)
		loc <- match(model, model.names)
		s.order.pass1 <- order(m.all[loc,], decreasing = decreasing.order)
		m2.pass1 <- m.all[, s.order.pass1]
		correl <- cor(t(m2.pass1))[, loc]
		m.order.pass1 <- order(correl, decreasing=T)
		m2.pass1 <- m2.pass1[m.order.pass1, ]
	}
	
	dev.off()
	
	MI.list.pass1.top10 = MI.list.pass1[m.order.pass1.top10]
	MI.list.pass1 = MI.list.pass1[m.order.pass1]
	
	bin.class.pass1 = bin.class.pass1[s.order.pass1]
#	m2.pass1 <- m2.pass1[m.order.pass1, ]
	model.descs2.pass1.top10 = model.descs2.pass1[m.order.pass1.top10]	
	model.descs2.pass1 = model.descs2.pass1[m.order.pass1]
	
	sample.names2.pass1 <- colnames(m2.pass1)
	model.names.pass1.top10 <- rownames(m2.pass1)
	print(matrix(c(model.names.pass1.top10, model.descs2.pass1.top10), ncol=2), quote=F)
#	browser()
	if (is.vector(cls.labels)) {
		cls.labels2.pass1 <- cls.labels.pass1[s.order.pass1]
		cls.list2.pass1 <- cls.list.pass1[s.order.pass1]
	} else {
		cls.labels2.pass1 <- cls.labels.pass1[, s.order.pass1]
		cls.list2.pass1 <- cls.list.pass1[, s.order.pass1]
	}
	m.all = m.all[, s.order.pass1]
#	browser()
	winning.model.ind.pass1 = which(model.names.pass1.top10[1] == rownames(m.all))
	
#	pathway.name <- "KRAS_ALL_UP"
#	pathway <- m[1,]
#	pathway0 <- ifelse(pathway < median(pathway), 0, 1) # disctretized version
	
#	MI.ref.genes.pass1 <- mutual.inf.2(m[1,], m[1,])$MI
	
#	browser()
#	m.score.pass1 <- m2.pass1[1,]
#	m.score.norm.pass1 <- (m.score.pass1 - min(m.score.pass1))/(max(m.score.pass1) - min(m.score.pass1))
#	m.score.pass1 = ifelse( m.score.pass1 < median(m.score.pass1), -1, 1)   # discretized version
#	MI.ref.genes.pass1 <- mutual.inf.2(m.score.norm.pass1, m.score.norm.pass1)$MI
#	print(paste("MI.ref.genes.pass1 =", MI.ref.genes.pass1))
	MI.list.phen.pass1 = vector(mode="numeric", length=n.phen.pass1)
#	FDR.list.phen.pass1 = vector(mode="numeric", length=n.phen.pass1)
	phen.descs2.pass1 = vector(mode="character", length=n.phen.pass1)
	
	if( length(unique(bin.class.pass1)) > 1){
#		MI.results <-(mutual.inf.3(bin.class.pass1, m.all, 
#							winning.model.ind.pass1, gene.target.name = phen.pass1[1]))#/MI.ref.genes.pass1
		MI.signif <- signif(MI.list.pass1[1], digits=3)
		MI.list.phen.pass1[1] = MI.list.pass1[1]
#		FDR.signif <- signif(FDR.list.pass1[1], digits=3)
#		FDR.list.phen.pass1[1] = FDR.list.pass1[1]
	} else{
		MI.signif <- "-"
		MI.list.phen.pass1[1] = NA
#		FDR.signif <- "- "
#		FDR.list.phen.pass1[1] = NA
	}
	print(paste(format(phen.pass1[1], width=12), "mutual.inf =", MI.signif))
	phen.descs2.pass1[1] <- paste(MI.signif)
#	browser()
	if( n.phen >= 2 ){
		bin.gene.matrix = ifelse(cls.list2.pass1[-1,]=="WT", 0, 1)
		n.aberrations = rowSums(bin.gene.matrix) #apply(bin.gene.matrix, MARGIN=1, FUN=sum)
		u.n.aberrations = unique(n.aberrations[n.aberrations != 0])
		for( i in 1:length(u.n.aberrations)){
			ind.without.SUMMARY = which(n.aberrations == u.n.aberrations[i])
			ind.master = ind.without.SUMMARY + 1
#		browser()
#		bin.gene.matrix.temp = bin.gene.matrix[ind.without.SUMMARY,]
			
			MI.results = mutual.inf.3.v2(bin.gene.matrix[ind.without.SUMMARY,], 
					m.all, winning.model.ind.pass1, gene.target.name=phen.pass1[ind.master],
					n.randomizations = n.randomizations)
			
			MI.list.phen.pass1[ind.master] = MI.results$MI
			#		FDR.list.phen.pass1[ind.master] = MI.results$FDR
			for( j in 1:length(ind.master)){
				phen.descs.pass1[ind.master[j]] = 
						paste( signif(MI.results$MI[j], digits=3), 
								sep="")
			}
		}
		ind.zeros = which(n.aberrations==0) + 1
		MI.list.phen.pass1[ind.zeros] = NA
#		FDR.list.phen.pass1[ind.zeros] = NA
		phen.descs.pass1[ind.zeros] = " - "
	}
	
	phen.names.pass1 = phen.pass1#[g.order.pass1]#[1:n.phen.pass1]
#	browser(text="Figure out how to print phen.descs.pass1 and phen.names.pass1 in a nice table")
	
	# Recompute cls.list2 as some mutations or copy numbers may have been removed
	
	
	# Recompute cls.phen and cls.labels2 as order may have changed
	
	cls.phen2.pass1 <- NULL
	if (is.vector(cls.labels)) {
		classes <- unique(cls.list2.pass1)
		cls.phen2.pass1 <- classes
		cls.labels2.pass1 <- match(cls.list2.pass1, cls.phen2.pass1)
	} else {
		for (kk in 1:length(cls.list2.pass1[, 1])) {
			classes <- unique(cls.list2.pass1[kk,])
#            cls.phen2[[kk]] <- classes
			cls.phen2.pass1 <- c(cls.phen2.pass1, classes)
			cls.labels2.pass1[kk,] <- match(cls.list2.pass1[kk,], classes)
		}
	}
#	cls.labels2.pass1 = cls.labels2.pass1[1:n.phen.pass1,]
	
	
#	browser()
#	correl <- cor(t(m2))[, loc]
#	m.order <- order(correl, decreasing=decreasing.order)
#	correl2 <- correl[m.order]
	
#	model.descs2 <- paste(model.descs[m.order], signif(correl2, digits=3))
	phen.list.pass1 <- unlist(cls.phen2.pass1)
	colors.list.pass1 = rep( "gray", length(phen.list.pass1))
	colors.list.pass1[phen.list.pass1=="MUT"] = cls.phen.colors[1]
	colors.list.pass1[phen.list.pass1=="DEL"] = cls.phen.colors[3]
	colors.list.pass1[phen.list.pass1=="AMP"] = cls.phen.colors[4]
	colors.list.pass1[phen.list.pass1=="ALT"] = cls.phen.colors[5]
#	browser()
#	colors.list.pass1[1,] = grey(bin.class.pass1)
#	print("cls.phen2:")
#	print(unlist(cls.phen2))
#	
#	print("cls.phen:")
#	print(unlist(cls.phen))
#	
#	print("colors.list:")
#	print(colors.list)
	
#	browser()
	
	filename <- paste(results.dir, test.file.prefix, ".Phase1.MI|HXY", sep="")
#    pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 9)
	pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 17 )
#	browser()
#   windows(width=12, height=8)
	MSIG.HeatMapPlot.9(V = m2.pass1, 
#			pathway.mut = bin.class.pass1,
			row.names = model.names.pass1.top10,
			row.names2 = model.descs2.pass1.top10, 
			col.labels = cls.labels2.pass1, 
			col.classes = cls.phen2.pass1, 
			phen.cmap = colors.list.pass1, 
			phen.names = phen.names.pass1,
			phen.names2 = phen.descs2.pass1,
			col.names = sample.names2.pass1, 
			main = paste(tissue, test.file.prefix), 
			xlab="  ", ylab="  ", row.norm = row.norm,  
			cmap.type = cmap.type, char.rescale = char.rescale,  legend=F)
	
	dev.off()
#	browser()
	if (!is.na(output.dataset)) {
#		V.GCT <- m.all
#		colnames(V.GCT) <- sample.names2
#		row.names(V.GCT) <- model.names2
		write.gct(gct.data.frame = m.all, descs = model.descs2.pass1, filename = paste(output.dataset, ".gct", sep=""))
		write.cls.2( class.v = cls.labels2.pass1, phen = cls.phen, filename = paste(output.dataset, ".cls", sep=""))
	}
	
}




MSIG.HeatMapPlot.11 <- function(
		## For Plotting expression heatmap only!! (No phenotypes)
		V, 
		row.names = "NA",
		row.names2 = "NA", 
		col.labels = "NA",
		col.labels2 = "NA", 
		col.classes = "NA", 
		phen.cmap = NULL, 
		col.names = "NA",
		phen.names = "NA", 
		phen.names2 = "NA",
		main = " ", 
		sub = " ", 
		xlab=" ", 
		ylab=" ",
		row.norm = TRUE,
		char.rescale = 0.85,                               
		cmap.type = 1,   # 1 = vintage pinkogram, 2 = scale of blues, 3 = high-resolution pinkogram for scores or probabilities [0, 1], 4 = high-resolution pinkogram for general values, 5 = color map for normalized enrichment scores, 6 = scale of red purples, 7 = scale of Oranges, 8 = scale of Greens, 9 = scale of Blues
		max.v = "NA",
		legend = T)
{
#
# Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
	
	n.rows <- length(V[,1])
	n.cols <- length(V[1,])
	V1 <- matrix(0, nrow=n.rows, ncol=n.cols)
	
#       if ((cmap.type == 5) | (cmap.type == 3)) {
	if (cmap.type == 5) {
		row.norm <- F
	}
	
	if (row.norm == TRUE) {
		row.mean <- apply(V, MARGIN=1, FUN=mean)
		row.sd <- apply(V, MARGIN=1, FUN=sd)
		row.n <- length(V[,1])
		for (i in 1:n.rows) {
			if (row.sd[i] == 0) {
				V1[i,] <- 0
			} else {
				V1[i,] <- (V[i,] - row.mean[i])/(0.333 * row.sd[i])
			}
			V1[i,] <- ifelse(V1[i,] < -4, -4, V1[i,])
			V1[i,] <- ifelse(V1[i,] > 4, 4, V1[i,])
		}
	} else {
		V1 <- V
	}
	
	if (cmap.type == 1) { 
		mycol <- c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA",
				"#FF9DB0", "#FF7080", 
				"#FF5A5A", "#FF4040", "#FF0D1D") # blue-pinkogram colors. This is the 1998-vintage,
		# pre-gene cluster, original pinkogram color map
	} else if (cmap.type == 2) {
		violet.palette <- colorRampPalette(c("#400030", "white"), space = "rgb")
		mycol <- rev(violet.palette(20))
		
#          mycol <- c("#FCFBFD","#F4F2F8","#F8F7FB","#EFEDF5","#E1E1EF","#E8E7F2","#DADAEB","#C6C7E1","#D0D1E6",
#                        "#BCBDDC","#A8A6CF",
#                        "#B2B2D6","#9E9AC8","#8A87BF","#9491C4","#807DBA","#7260AB","#796FB3","#6A51A3","#5C3596",
#                        "#63439D","#54278F","#460D83","#4D1A89","#3F007D")
	} else if (cmap.type == 6) {
		mycol <- c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E",
				"#7A0177", "#49006A")
	} else if (cmap.type == 7) {
		mycol <- c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801",
				"#A63603", "#7F2704")
	} else if (cmap.type == 8) {
		mycol <- c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45",
				"#006D2C", "#00441B")
	} else if (cmap.type == 9) {
		mycol <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5",
				"#08519C", "#08306B")
	} else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 5)) {
		mycol <- vector(length=512, mode = "numeric")
		
		for (k in 1:256) {
			mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
		}
		for (k in 257:512) {
			mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
		}
		mycol <- rev(mycol)
	}
	
	ncolors <- length(mycol)
	
	if (cmap.type == 5) {
		if (max.v == "NA") {
			max.v <- max(max(V1), -min(V1))
		}
		V2 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
		
	} else {
		V2 <- ceiling(ncolors * (V1 - min(V1))/(1.001*(max(V1) - min(V1))))
	}
	
	if (col.labels[1] == "NA") {      
		heatm <- matrix(0, nrow = n.rows, ncol = n.cols)
		heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
		tot.cols <- ncolors
		browser()
		if (legend == T) {
			nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(5, 1), heights = c(10, 1), respect = FALSE)
		} else {
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(8, 1), respect = FALSE)
		}
		par(mar = c(3, 16, 3, 16))
		mycol <- c(mycol, phen.cmap[1:length(col.classes)])
		image(1:n.cols, 1:n.rows, t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=main, sub = sub, xlab= xlab, ylab=ylab)
		n.rows.phen <- 0
	} else {
		tot.cols <- ncolors
		if (is.vector(col.labels)) {
			heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
			heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
			n.rows.phen <- 1
			heatm[n.rows + 1,] <- tot.cols + col.labels
			cols.row <- length(unique(col.labels))
			tot.cols <- tot.cols + cols.row
			phen.cmap <- phen.cmap[1:cols.row]
		} else {
			n.rows.phen <- length(col.labels[,1])
			cols.row <- vector(length=n.rows.phen, mode = "numeric")
			heatm <- matrix(0, nrow = n.rows + n.rows.phen, ncol = n.cols)
			heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
			for (k in seq(n.rows + n.rows.phen, n.rows + 1, -1)) {
				heatm[k,] <- tot.cols + col.labels[n.rows + n.rows.phen - k + 1,]
				cols.row[n.rows + n.rows.phen - k + 1] <- length(unique(col.labels[n.rows + n.rows.phen - k + 1,]))
				tot.cols <- tot.cols + cols.row[n.rows + n.rows.phen - k + 1]
#                 print(c("col:", k, ":", tot.cols + col.labels[n.rows + n.rows.phen - k + 1,], "tot.cols:", tot.cols))
				
			}
			phen.cmap <- phen.cmap[1:sum(unlist(lapply(col.classes, length)))]
		}
		if (legend == T) {
#              nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(10, 2), heights = c(6, 1), respect = FALSE)
			nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow=T), heights = c(8, 4, 1), respect = FALSE)
		} else {
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(5, 1), respect = FALSE)
		}
		par(mar = c(3, 16, 3, 16))
		mycol <- c(mycol, phen.cmap)
		image(1:n.cols, 1:(n.rows + n.rows.phen), t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=main, sub = sub, xlab= xlab, ylab=ylab)
	}
	
# Add lines to separate phenotypes or subgroups
	
	if (col.labels2[1] != "NA") {
		groups <-  split(col.labels2, col.labels2)
		len.vec <- lapply(groups, length)
		plot.div <- c(0.51, cumsum(len.vec) + 0.5)
		for (i in plot.div) {
			lines(c(i, i), c(0, n.rows + n.rows.phen + 0.48), lwd = 2, cex = 0.9, col = "black")
		}
		lines(c(0.51, n.cols + 0.49), c(0.51, 0.51), lwd = 2, cex = 0.9, col = "black")
		lines(c(0.51, n.cols + 0.49), c(n.rows + n.rows.phen + 0.48, n.rows + n.rows.phen + 0.48), lwd = 2,
				cex = 0.9, col = "black")
		lines(c(0.51, n.cols + 0.49), c(n.rows + 0.50, n.rows + 0.50), lwd = 2,
				cex = 0.9, col = "black")
	}
	if (row.names[1] != "NA") {
#		browser()
		numC <- nchar(row.names)
		size.row.char <- char.rescale*25/(n.rows + 20)
		for (i in 1:n.rows) {
			row.names[i] <- substr(row.names[i], 1, 40)
			row.names[i] <- paste(row.names[i], " ", sep="")
		}
		if (phen.names[1] == "NA") {
			head.names <- paste("Class", seq(n.rows.phen, 1, -1))
		} else {
			head.names <- as.character(rev(phen.names))
		}
		row.names <- c(row.names[seq(n.rows, 1, -1)], head.names)
#            print(paste("n.rows:", n.rows))
#            print(paste("Phen names:", phen.names))
#            print(paste("Head names:", head.names))
#            print(paste("Row names:", row.names))
		axis(2, at=1:(n.rows + n.rows.phen), labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char,
				font.axis=2, line=-1)
	}
	
	if (row.names2[1] != "NA") {
#		browser()
		numC <- nchar(row.names2)
		size.row.char <- char.rescale*25/(n.rows + 20)
		for (i in 1:n.rows) {
			row.names2[i] <- substr(row.names2[i], 1, 40)
			row.names2[i] <- paste(" ", row.names2[i], sep="")
			
		}
		for( i in 1:n.rows.phen ){
			phen.names2[i] <- substr(phen.names2[i], 1, 40)
			phen.names2[i] <- paste( " ", phen.names2[i], sep="")
		}
		
		row.names2 <- rev(row.names2)
		phen.names2 <- rev(phen.names2)
		axis(4, at=1:(n.rows + n.rows.phen), labels=c(row.names2, phen.names2), adj= 0.5, tick=FALSE, las = 1, 
				cex.axis=size.row.char, font.axis=2, line=-1)
	}
	
	if (col.names[1] != "NA") {
		size.col.char <- char.rescale*20/(n.cols + 25)
		axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
	}
	
	# Phenotype Legend 
	
#      print("--------------------------------------------------------------------------------------------")
	if (legend == T) {
		leg.txt <- NULL
		p.vec <- NULL
		c.vec <- NULL
		c2.vec <- NULL
		ind <- 1
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		for (i in 1:n.rows.phen) {  
			if (is.vector(col.labels)) {
				phen.v <- as.character(col.classes)
			} else {
				phen.v <- as.character(col.classes[[i]])
			}
			p.name <- paste(as.character(rev(head.names)[i]), ":   ", sep="")
			leg.txt <- c(p.name, phen.v)  
			p.vec <-  rep(22, cols.row[i] + 1)
			c.vec <-  c("#FFFFFF", phen.cmap[ind:(ind + cols.row[i] - 1)])
			c2.vec <- c("#FFFFFF", rep("black", cols.row[i]))
			ind <- ind + cols.row[i]
			offset <- 0.07
			legend(x=0, y= 1 - offset*i, 
					horiz = T, x.intersp = 0.5, legend=leg.txt, bty="n", xjust=0, yjust= 1, pch = p.vec, 
					pt.bg = c.vec, col = c2.vec, cex = 1.20, pt.cex=1.75)
		}
	}
	
	# Color map legend
	
#       print(c("range V=", range(V)))
#       print(c("range V1=", range(V1)))
#       print(c("range V2=", range(V2)))
	
	par(mar = c(2, 12, 2, 12))
	num.v <- 20
	range.v <- range(V2)
	incr <-  (range.v[1] - range.v[2])/(num.v - 1)
	heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow=num.v, ncol=1)
	image(1:num.v, 1:1, heatm.v, zlim = c(0, tot.cols), col=mycol, axes=FALSE,
			main=" ", sub = " ", xlab= ylab, ylab=xlab)
	range.v <- range(V1)
	incr <-  (range.v[1] - range.v[2])/(num.v - 1)
	heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits=2), nrow=num.v, ncol=1)
#          print(c("heatm.v2=", heatm.v2))
	axis(3, at=1:num.v, labels=heatm.v2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.5*char.rescale, font.axis=1)
	
	return()
	
}

MSIG.HeatMapPlot.9 <- function(
		V, 
		row.names = "NA",
		row.names2 = "NA", 
		col.labels = "NA",
		col.labels2 = "NA", 
		col.classes = "NA", 
		phen.cmap = "NA", 
		col.names = "NA",
		phen.names = "NA", 
		phen.names2 = "NA",
		main = " ", 
		sub = " ", 
		xlab=" ", 
		ylab=" ",
		row.norm = TRUE,
		char.rescale = 0.85,                               
		cmap.type = 1,   # 1 = vintage pinkogram, 2 = scale of blues, 3 = high-resolution pinkogram for scores or probabilities [0, 1], 4 = high-resolution pinkogram for general values, 5 = color map for normalized enrichment scores, 6 = scale of red purples, 7 = scale of Oranges, 8 = scale of Greens, 9 = scale of Blues
		max.v = "NA",
		legend = F)
{
#
# Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
# Doesn't plot the spectrum on the bottom
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
	
	n.rows <- length(V[,1])
	n.cols <- length(V[1,])
	V1 <- matrix(0, nrow=n.rows, ncol=n.cols)
	
#       if ((cmap.type == 5) | (cmap.type == 3)) {
	if (cmap.type == 5) {
		row.norm <- F
	}
	
	if (row.norm == TRUE) {
		row.mean <- apply(V, MARGIN=1, FUN=mean)
		row.sd <- apply(V, MARGIN=1, FUN=sd)
		row.n <- length(V[,1])
		for (i in 1:n.rows) {
			if (row.sd[i] == 0) {
				V1[i,] <- 0
			} else {
				V1[i,] <- (V[i,] - row.mean[i])/(0.333 * row.sd[i])
			}
			V1[i,] <- ifelse(V1[i,] < -4, -4, V1[i,])
			V1[i,] <- ifelse(V1[i,] > 4, 4, V1[i,])
		}
	} else {
		V1 <- V
	}
	
	if (cmap.type == 1) { 
		mycol <- c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA",
				"#FF9DB0", "#FF7080", 
				"#FF5A5A", "#FF4040", "#FF0D1D") # blue-pinkogram colors. This is the 1998-vintage,
		# pre-gene cluster, original pinkogram color map
	} else if (cmap.type == 2) {
		violet.palette <- colorRampPalette(c("#400030", "white"), space = "rgb")
		mycol <- rev(violet.palette(20))
		
#          mycol <- c("#FCFBFD","#F4F2F8","#F8F7FB","#EFEDF5","#E1E1EF","#E8E7F2","#DADAEB","#C6C7E1","#D0D1E6",
#                        "#BCBDDC","#A8A6CF",
#                        "#B2B2D6","#9E9AC8","#8A87BF","#9491C4","#807DBA","#7260AB","#796FB3","#6A51A3","#5C3596",
#                        "#63439D","#54278F","#460D83","#4D1A89","#3F007D")
	} else if (cmap.type == 6) {
		mycol <- c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E",
				"#7A0177", "#49006A")
	} else if (cmap.type == 7) {
		mycol <- c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801",
				"#A63603", "#7F2704")
	} else if (cmap.type == 8) {
		mycol <- c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45",
				"#006D2C", "#00441B")
	} else if (cmap.type == 9) {
		mycol <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5",
				"#08519C", "#08306B")
	} else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 5)) {
		mycol <- vector(length=512, mode = "numeric")
		
		for (k in 1:256) {
			mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
		}
		for (k in 257:512) {
			mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
		}
		mycol <- rev(mycol)
	}
	
	ncolors <- length(mycol)
	
	if (cmap.type == 5) {
		if (max.v == "NA") {
			max.v <- max(max(V1), -min(V1))
		}
		V2 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
		
	} else {
		V2 <- ceiling(ncolors * (V1 - min(V1))/(1.001*(max(V1) - min(V1))))
	}
	
	if (col.labels[1] == "NA") {      
		heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
		heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
		tot.cols <- ncolors
		if (legend == T) {
			nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(5, 1), heights = c(10, 1), respect = FALSE)
		} else {
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(8, 1), respect = FALSE)
		}
		par(mar = c(3, 16, 3, 16))
		mycol <- c(mycol, phen.cmap[1:length(col.classes)])
		image(1:n.cols, 1:n.rows, t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=main, sub = sub, xlab= xlab, ylab=ylab)
		n.rows.phen <- 0
	} else {
		tot.cols <- ncolors
		if (is.vector(col.labels)) {
			heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
			heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
			n.rows.phen <- 1
			heatm[n.rows + 1,] <- tot.cols + col.labels
			cols.row <- length(unique(col.labels))
			tot.cols <- tot.cols + cols.row
			phen.cmap <- phen.cmap[1:cols.row]
		} else {
			n.rows.phen <- length(col.labels[,1])
			cols.row <- vector(length=n.rows.phen, mode = "numeric")
			heatm <- matrix(0, nrow = n.rows + n.rows.phen, ncol = n.cols)
			heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
			for (k in seq(n.rows + n.rows.phen, n.rows + 1, -1)) {
				heatm[k,] <- tot.cols + col.labels[n.rows + n.rows.phen - k + 1,]
				cols.row[n.rows + n.rows.phen - k + 1] <- length(unique(col.labels[n.rows + n.rows.phen - k + 1,]))
				tot.cols <- tot.cols + cols.row[n.rows + n.rows.phen - k + 1]
#                 print(c("col:", k, ":", tot.cols + col.labels[n.rows + n.rows.phen - k + 1,], "tot.cols:", tot.cols))
				
			}
			phen.cmap <- phen.cmap[1:sum(unlist(lapply(col.classes, length)))]
		}
		if (legend == T) {
#              nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(10, 2), heights = c(6, 1), respect = FALSE)
			nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow=T), heights = c(8, 4, 1), respect = FALSE)
		} else {
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(5, 1), respect = FALSE)
		}
		par(mar = c(3, 16, 3, 16))
		mycol <- c(mycol, phen.cmap)
		image(1:n.cols, 1:(n.rows + n.rows.phen), t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=main, sub = sub, xlab= xlab, ylab=ylab)
	}
	
# Add lines to separate phenotypes or subgroups
	
	if (col.labels2[1] != "NA") {
		groups <-  split(col.labels2, col.labels2)
		len.vec <- lapply(groups, length)
		plot.div <- c(0.51, cumsum(len.vec) + 0.5)
		for (i in plot.div) {
			lines(c(i, i), c(0, n.rows + n.rows.phen + 0.48), lwd = 2, cex = 0.9, col = "black")
		}
		lines(c(0.51, n.cols + 0.49), c(0.51, 0.51), lwd = 2, cex = 0.9, col = "black")
		lines(c(0.51, n.cols + 0.49), c(n.rows + n.rows.phen + 0.48, n.rows + n.rows.phen + 0.48), lwd = 2,
				cex = 0.9, col = "black")
		lines(c(0.51, n.cols + 0.49), c(n.rows + 0.50, n.rows + 0.50), lwd = 2,
				cex = 0.9, col = "black")
	}
	if (row.names[1] != "NA") {
#		browser()
		numC <- nchar(row.names)
		size.row.char <- char.rescale*25/(n.rows + 20)
		for (i in 1:n.rows) {
			row.names[i] <- substr(row.names[i], 1, 40)
			row.names[i] <- paste(row.names[i], " ", sep="")
		}
		if (phen.names[1] == "NA") {
			head.names <- paste("Class", seq(n.rows.phen, 1, -1))
		} else {
			head.names <- as.character(rev(phen.names))
		}
		row.names <- c(row.names[seq(n.rows, 1, -1)], head.names)
#            print(paste("n.rows:", n.rows))
#            print(paste("Phen names:", phen.names))
#            print(paste("Head names:", head.names))
#            print(paste("Row names:", row.names))
		axis(2, at=1:(n.rows + n.rows.phen), labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char,
				font.axis=2, line=-1)
	}
	
	if (row.names2[1] != "NA") {
#		browser()
		numC <- nchar(row.names2)
		size.row.char <- char.rescale*25/(n.rows + 20)
		for (i in 1:n.rows) {
			row.names2[i] <- substr(row.names2[i], 1, 40)
			row.names2[i] <- paste(" ", row.names2[i], sep="")
			
		}
		for( i in 1:n.rows.phen ){
			phen.names2[i] <- substr(phen.names2[i], 1, 40)
			phen.names2[i] <- paste( " ", phen.names2[i], sep="")
		}
		
		row.names2 <- rev(row.names2)
		phen.names2 <- rev(phen.names2)
		axis(4, at=1:(n.rows + n.rows.phen), labels=c(row.names2, phen.names2), adj= 0.5, tick=FALSE, las = 1, 
				cex.axis=size.row.char, font.axis=2, line=-1)
	}
	
	if (col.names[1] != "NA") {
		size.col.char <- char.rescale*20/(n.cols + 25)
		axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
	}
	
	# Phenotype Legend 
	
#      print("--------------------------------------------------------------------------------------------")
	if (legend == T) {
		leg.txt <- NULL
		p.vec <- NULL
		c.vec <- NULL
		c2.vec <- NULL
		ind <- 1
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		for (i in 1:n.rows.phen) {  
			if (is.vector(col.labels)) {
				phen.v <- as.character(col.classes)
			} else {
				phen.v <- as.character(col.classes[[i]])
			}
			p.name <- paste(as.character(rev(head.names)[i]), ":   ", sep="")
			leg.txt <- c(p.name, phen.v)  
			p.vec <-  rep(22, cols.row[i] + 1)
			c.vec <-  c("#FFFFFF", phen.cmap[ind:(ind + cols.row[i] - 1)])
			c2.vec <- c("#FFFFFF", rep("black", cols.row[i]))
			ind <- ind + cols.row[i]
			offset <- 0.07
			legend(x=0, y= 1 - offset*i, 
					horiz = T, x.intersp = 0.5, legend=leg.txt, bty="n", xjust=0, yjust= 1, pch = p.vec, 
					pt.bg = c.vec, col = c2.vec, cex = 1.20, pt.cex=1.75)
		}
	}
	
	# Color map legend
	
#       print(c("range V=", range(V)))
#       print(c("range V1=", range(V1)))
#       print(c("range V2=", range(V2)))
	
#	par(mar = c(2, 12, 2, 12))
#	num.v <- 20
#	range.v <- range(V2)
#	incr <-  (range.v[1] - range.v[2])/(num.v - 1)
#	heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow=num.v, ncol=1)
#	image(1:num.v, 1:1, heatm.v, zlim = c(0, tot.cols), col=mycol, axes=FALSE,
#			main=" ", sub = " ", xlab= ylab, ylab=xlab)
#	range.v <- range(V1)
#	incr <-  (range.v[1] - range.v[2])/(num.v - 1)
#	heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits=2), nrow=num.v, ncol=1)
	##          print(c("heatm.v2=", heatm.v2))
#	axis(3, at=1:num.v, labels=heatm.v2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.5*char.rescale, font.axis=1)
	
	return()
	
}

MSIG.HeatMapPlot.10 <- function(
		V, 
		pathway.mut,
		row.names = "NA",
		row.names2 = "NA", 
		col.labels = "NA",
		col.labels2 = "NA", 
		col.classes = "NA", 
		phen.cmap = "NA", 
		col.names = "NA",
		phen.names = "NA", 
		phen.names2 = "NA",
		main = " ", 
		sub = " ", 
		xlab=" ", 
		ylab=" ",
		row.norm = TRUE,
		char.rescale = 0.85,                               
		cmap.type = 1,   # 1 = vintage pinkogram, 2 = scale of blues, 3 = high-resolution pinkogram for scores or probabilities [0, 1], 4 = high-resolution pinkogram for general values, 5 = color map for normalized enrichment scores, 6 = scale of red purples, 7 = scale of Oranges, 8 = scale of Greens, 9 = scale of Blues
		max.v = "NA",
		legend = T,
		tissue.names = "NA",
		tissue.labels = NA,
		ifScoreBarPlot = FALSE,
		MI.list.model = NULL,
		MI.list.phen  = NULL)
{
#
# Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
# Doesn't plot the spectrum on the bottom
#
# Plots PATHWAY.MUT as a continuous vector in a greyscale spectrum
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
	
	n.rows <- length(V[,1])
	n.cols <- length(V[1,])
	V1 <- matrix(0, nrow=n.rows, ncol=n.cols)
	set3 = brewer.pal(12, "Set3")
	accent = brewer.pal(8, "Accent")
	
#       if ((cmap.type == 5) | (cmap.type == 3)) {
	if (cmap.type == 5) {
		row.norm <- F
	}
	
	if(ifScoreBarPlot) char.rescale = 1.5*char.rescale
	
	if (row.norm == TRUE) {
		row.mean <- apply(V, MARGIN=1, FUN=mean)
		row.sd <- apply(V, MARGIN=1, FUN=sd)
		row.n <- length(V[,1])
		for (i in 1:n.rows) {
			if (row.sd[i] == 0) {
				V1[i,] <- 0
			} else {
				V1[i,] <- (V[i,] - row.mean[i])/(0.333 * row.sd[i])
			}
			V1[i,] <- ifelse(V1[i,] < -4, -4, V1[i,])
			V1[i,] <- ifelse(V1[i,] > 4, 4, V1[i,])
		}
	} else {
		V1 <- V
	}
	
	if (cmap.type == 1) { 
		mycol <- c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA",
				"#FF9DB0", "#FF7080", 
				"#FF5A5A", "#FF4040", "#FF0D1D") # blue-pinkogram colors. This is the 1998-vintage,
		# pre-gene cluster, original pinkogram color map
	} else if (cmap.type == 2) {
		violet.palette <- colorRampPalette(c("#400030", "white"), space = "rgb")
		mycol <- rev(violet.palette(20))
		
#          mycol <- c("#FCFBFD","#F4F2F8","#F8F7FB","#EFEDF5","#E1E1EF","#E8E7F2","#DADAEB","#C6C7E1","#D0D1E6",
#                        "#BCBDDC","#A8A6CF",
#                        "#B2B2D6","#9E9AC8","#8A87BF","#9491C4","#807DBA","#7260AB","#796FB3","#6A51A3","#5C3596",
#                        "#63439D","#54278F","#460D83","#4D1A89","#3F007D")
	} else if (cmap.type == 6) {
		mycol <- c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E",
				"#7A0177", "#49006A")
	} else if (cmap.type == 7) {
		mycol <- c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801",
				"#A63603", "#7F2704")
	} else if (cmap.type == 8) {
		mycol <- c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45",
				"#006D2C", "#00441B")
	} else if (cmap.type == 9) {
		mycol <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5",
				"#08519C", "#08306B")
	} else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 5)) {
		mycol <- vector(length=512, mode = "numeric")
		
		for (k in 1:256) {
			mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
		}
		for (k in 257:512) {
			mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
		}
		mycol <- rev(mycol)
	}
#	browser()
	ncolors <- length(mycol)
	pathway.mut = (-(pathway.mut*.749 + 0.251 - 1))
#	image(1:n.cols, 1, as.matrix(pathway.mut), col=gray(n.cols:0/n.cols))
	if (cmap.type == 5) {
		if (max.v == "NA") {
			max.v <- max(max(V1), -min(V1))
		}
		V2 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
		
	} else {
		V2 <- ceiling(ncolors * (V1 - min(V1))/(1.001*(max(V1) - min(V1))))
	}
	
	if (col.labels[1] == "NA") {      
		heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
		heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
		tot.cols <- ncolors
		if (legend == T && ifScoreBarPlot == FALSE) {
			print("legend == T && ifScoreBarPlot == FALSE")
			nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(5, 1), heights = c(10, 1), respect = FALSE)
		} else if( legend == FALSE && ifScoreBarPlot == FALSE ){
			print("legend == FALSE && ifScoreBarPlot == FALSE")
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(8, 1), respect = FALSE)
		} else if( legend == TRUE && ifScoreBarPlot == TRUE){
			print("legend == TRUE && ifScoreBarPlot == TRUE")
			nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(5, 1), heights = c(10, 1), respect = FALSE)
		} else if( legend == FALSE && ifScoreBarPlot == TRUE ){
			print("legend == FALSE && ifScoreBarPlot == TRUE")
			nf <- layout(matrix(c(1, 2, 3, 3 ), ncol=2, byrow=TRUE), 
					heights = c(8, 1), widths=c(10,1), respect = FALSE)
		}
		par(mar = c(3, 16, 3, 16))
		mycol <- c(mycol, phen.cmap[1:length(col.classes)])
		image(1:n.cols, 1:n.rows, t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=main, sub = sub, xlab= xlab, ylab=ylab)
		n.rows.phen <- 0
	} else {
		tot.cols <- ncolors
		if (is.vector(col.labels)) {
			heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
			heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
			n.rows.phen <- 1
			heatm[n.rows + 1,] <- tot.cols + col.labels
			cols.row <- length(unique(col.labels))
			tot.cols <- tot.cols + cols.row
			phen.cmap <- phen.cmap[1:cols.row]
			pathway.mut.grey = grey(pathway.mut)
			u.pathway.mut.grey = unique(pathway.mut.grey)
			heatm[n.rows + n.rows.phen,] = match(pathway.mut.grey, u.pathway.mut.grey) + tot.cols
		} else {
			n.rows.phen <- length(col.labels[,1])
			cols.row <- vector(length=n.rows.phen, mode = "numeric")
			heatm <- matrix(0, nrow = n.rows + n.rows.phen, ncol = n.cols)
			heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
#			heatm[n.rows+n.rows.phen,] = t(as.matrix(gray(pathway.mut)))
			for (k in seq(n.rows + n.rows.phen, n.rows + 1, -1)) {
				heatm[k,] <- tot.cols + col.labels[n.rows + n.rows.phen - k + 1,]
				cols.row[n.rows + n.rows.phen - k + 1] <- length(unique(col.labels[n.rows + n.rows.phen - k + 1,]))
				tot.cols <- tot.cols + cols.row[n.rows + n.rows.phen - k + 1]
#                 print(c("col:", k, ":", tot.cols + col.labels[n.rows + n.rows.phen - k + 1,], "tot.cols:", tot.cols))
				
			}
#			browser()
			pathway.mut.grey = grey(pathway.mut)
			u.pathway.mut.grey = unique(pathway.mut.grey)
			heatm[n.rows + n.rows.phen,] = match(pathway.mut.grey, u.pathway.mut.grey) + tot.cols
			tot.cols = tot.cols + length(u.pathway.mut.grey)
			phen.cmap <- phen.cmap[1:sum(unlist(lapply(col.classes, length)))]
		}
#		image(as.matrix(grey(pathway.mut)))
		if (legend == T && ifScoreBarPlot == FALSE
				) {
			print("legend == T && ifScoreBarPlot == FALSE")
#              nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(10, 2), heights = c(6, 1), respect = FALSE)
			nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow=T), heights = c(8, 4, 1), respect = FALSE)
		} else if( legend == FALSE && ifScoreBarPlot == FALSE 
				){
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(5, 1), respect = FALSE)
			print("legend == FALSE && ifScoreBarPlot == FALSE")
#			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(8, 1), respect = FALSE)
		} else if( legend == TRUE && ifScoreBarPlot == TRUE){
			print("legend == TRUE && ifScoreBarPlot == TRUE")
			#			nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(5, 1), heights = c(10, 1), respect = FALSE)
		} else if( legend == FALSE && ifScoreBarPlot == TRUE ){
			print("legend == FALSE && ifScoreBarPlot == TRUE")
			nf <- layout(matrix(c(1, 1, 2, 1, 1, 2, 3, 3, 2), ncol=3, byrow=F), heights = c(5, 5, 2), 
					widths=c(5,5,1), respect = FALSE)
#			layout.show(3)
#			browser()
		}
		if(ifScoreBarPlot){
			par( mar = c(8, 20, 8, 0))
		} else{
			par( mar = c(3, 16, 3, 16))
		}
#		par(mar = c(8, 20, 8, ifelse(ifScoreBarPlot, 0, 20)))
#		browser()
		mycol <- c(mycol, phen.cmap, u.pathway.mut.grey)
#		browser()
		image(1:n.cols, 1:(n.rows + n.rows.phen), t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=main, sub = sub, xlab= xlab, ylab=ylab, cex.main=ifelse(ifScoreBarPlot, 2, 1))
	}
	
# Add lines to separate phenotypes or subgroups
	
	if (col.labels2[1] != "NA") {
		groups <-  split(col.labels2, col.labels2)
		len.vec <- lapply(groups, length)
		plot.div <- c(0.51, cumsum(len.vec) + 0.5)
		for (i in plot.div) {
			lines(c(i, i), c(0, n.rows + n.rows.phen + 0.48), lwd = 2, cex = 0.9, col = "black")
		}
		lines(c(0.51, n.cols + 0.49), c(0.51, 0.51), lwd = 2, cex = 0.9, col = "black")
		lines(c(0.51, n.cols + 0.49), c(n.rows + n.rows.phen + 0.48, n.rows + n.rows.phen + 0.48), lwd = 2,
				cex = 0.9, col = "black")
		lines(c(0.51, n.cols + 0.49), c(n.rows + 0.50, n.rows + 0.50), lwd = 2,
				cex = 0.9, col = "black")
	}
	if (row.names[1] != "NA") {
#		browser()
		numC <- nchar(row.names)
		size.row.char <- ifelse(ifScoreBarPlot, char.rescale*25/(n.rows + 25), char.rescale*25/(n.rows + 25))
		for (i in 1:n.rows) {
			row.names[i] <- substr(row.names[i], 1, 40)
			row.names[i] <- paste(row.names[i], " ", sep="")
		}
		if (phen.names[1] == "NA") {
			head.names <- paste("Class", seq(n.rows.phen, 1, -1))
		} else {
			head.names <- as.character(rev(phen.names))
		}
		row.names <- c(row.names[seq(n.rows, 1, -1)], head.names)
#            print(paste("n.rows:", n.rows))
#            print(paste("Phen names:", phen.names))
#            print(paste("Head names:", head.names))
#            print(paste("Row names:", row.names))
#		browser()
		axis(2, at=1:(n.rows + n.rows.phen), labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char,
				font.axis=2, line=-1)
	}
	
	if (row.names2[1] != "NA") {
#		browser()
		numC <- nchar(row.names2)
		size.row.char <- char.rescale*25/(n.rows + 20)
		for (i in 1:n.rows) {
			row.names2[i] <- substr(row.names2[i], 1, 40)
			row.names2[i] <- paste(" ", row.names2[i], sep="")
			
		}
		for( i in 1:n.rows.phen ){
			phen.names2[i] <- substr(phen.names2[i], 1, 40)
			phen.names2[i] <- paste( " ", phen.names2[i], sep="")
		}
		
		row.names2 <- rev(row.names2)
		phen.names2 <- rev(phen.names2)
		axis(4, at=1:(n.rows + n.rows.phen), labels=c(row.names2, phen.names2), adj= 0.5, tick=FALSE, las = 1, 
				cex.axis=size.row.char, font.axis=2, line=-1)
	}
	
	if (col.names[1] != "NA") {
		size.col.char <- ifelse(ifScoreBarPlot, char.rescale*20/(n.cols + 15), char.rescale*20/(n.cols + 25))
		axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
	}
	
	
	
	# Phenotype Legend 
	
#      print("--------------------------------------------------------------------------------------------")
	if (legend == T) {
		leg.txt <- NULL
		p.vec <- NULL
		c.vec <- NULL
		c2.vec <- NULL
		ind <- 1
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		for (i in 1:n.rows.phen) {  
#			browser()
			if (is.vector(col.labels)) {
				phen.v <- as.character(col.classes)
			} else {
				phen.v <- as.character(col.classes[[i]])
			}
			p.name <- paste(as.character(rev(head.names)[i]), ":   ", sep="")
			leg.txt <- c(p.name, phen.v)  
			p.vec <-  rep(22, cols.row[i] + 1)
			c.vec <-  c("#FFFFFF", phen.cmap[ind:(ind + cols.row[i] - 1)])
			c2.vec <- c("#FFFFFF", rep("black", cols.row[i]))
			ind <- ind + cols.row[i]
			offset <- 0.07
			legend(x=0, y= 1 - offset*i, 
					horiz = T, x.intersp = 0.5, legend=leg.txt, bty="n", xjust=0, yjust= 1, pch = p.vec, 
					pt.bg = c.vec, col = c2.vec, cex = 1.20, pt.cex=1.75)
		}
	}
	
	
	
	
	
	# Color map legend
	
#       print(c("range V=", range(V)))
#       print(c("range V1=", range(V1)))
#       print(c("range V2=", range(V2)))
	
	par(mar = c(2, 12, 6, 12))
	num.v <- 20
	range.v <- range(V2)
	incr <-  (range.v[1] - range.v[2])/(num.v - 1)
	heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow=num.v, ncol=1)
	image(1:num.v, 1:1, heatm.v, zlim = c(0, tot.cols), col=mycol, axes=FALSE,
			main=" ", sub = " ", xlab= ylab, ylab=xlab)
	range.v <- range(V1)
	incr <-  (range.v[1] - range.v[2])/(num.v - 1)
	heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits=2), nrow=num.v, ncol=1)
	#          print(c("heatm.v2=", heatm.v2))
	axis(3, at=1:num.v, labels=heatm.v2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.5*char.rescale, font.axis=1)
	
#	browser()
	## Bar plots of MI values
	if(ifScoreBarPlot){
		par(mar = c(6, 5, 6, 1))
		barplot(rev(c(MI.list.phen, MI.list.model)), xlab="MI scores\nbarplot", ylab="", xlim=c(-1,1), 
				axes=TRUE, horiz=TRUE, axisnames=FALSE, 
#				main="MI scores\nbarplot", font.main=1
		)
#		browser()
#		text("MI scores\nbarplot")
	}
	return()
	
}

MSIG.HeatMapPlot.10.show.multiple.tissues <- function(
		V, 
		pathway.mut,
		row.names = "NA",
		row.names2 = "NA", 
		col.labels = "NA",
		col.labels2 = "NA", 
		col.classes = "NA", 
		phen.cmap = "NA", 
		col.names = "NA",
		phen.names = "NA", 
		phen.names2 = "NA",
		main = " ", 
		sub = " ", 
		xlab=" ", 
		ylab=" ",
		row.norm = TRUE,
		char.rescale = 0.85,                               
		cmap.type = 1,   # 1 = vintage pinkogram, 2 = scale of blues, 3 = high-resolution pinkogram for scores or probabilities [0, 1], 4 = high-resolution pinkogram for general values, 5 = color map for normalized enrichment scores, 6 = scale of red purples, 7 = scale of Oranges, 8 = scale of Greens, 9 = scale of Blues
		max.v = "NA",
		legend = T,
		tissue.names = "NA",
		tissue.labels = NA,
		highlight.tissue.name = NULL,
		ifScoreBarPlot = FALSE,
		MI.list.model = NULL,
		MI.list.phen  = NULL)
{
#
# Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
# Doesn't plot the spectrum on the bottom
#
# Plots PATHWAY.MUT as a continuous vector in a greyscale spectrum
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
	library(RColorBrewer)
	
	n.tissues = length(tissue.names)
	
	n.rows <- length(V[,1])
	n.cols <- length(V[1,])
	V1 <- matrix(0, nrow=n.rows, ncol=n.cols)
	
#       if ((cmap.type == 5) | (cmap.type == 3)) {
	if (cmap.type == 5) {
		row.norm <- F
	}
	
	if (row.norm == TRUE) {
		row.mean <- apply(V, MARGIN=1, FUN=mean)
		row.sd <- apply(V, MARGIN=1, FUN=sd)
		row.n <- length(V[,1])
		for (i in 1:n.rows) {
			if (row.sd[i] == 0) {
				V1[i,] <- 0
			} else {
				V1[i,] <- (V[i,] - row.mean[i])/(0.333 * row.sd[i])
			}
			V1[i,] <- ifelse(V1[i,] < -4, -4, V1[i,])
			V1[i,] <- ifelse(V1[i,] > 4, 4, V1[i,])
		}
	} else {
		V1 <- V
	}
	
	if (cmap.type == 1) { 
		mycol <- c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA",
				"#FF9DB0", "#FF7080", 
				"#FF5A5A", "#FF4040", "#FF0D1D") # blue-pinkogram colors. This is the 1998-vintage,
		# pre-gene cluster, original pinkogram color map
	} else if (cmap.type == 2) {
		violet.palette <- colorRampPalette(c("#400030", "white"), space = "rgb")
		mycol <- rev(violet.palette(20))
		
#          mycol <- c("#FCFBFD","#F4F2F8","#F8F7FB","#EFEDF5","#E1E1EF","#E8E7F2","#DADAEB","#C6C7E1","#D0D1E6",
#                        "#BCBDDC","#A8A6CF",
#                        "#B2B2D6","#9E9AC8","#8A87BF","#9491C4","#807DBA","#7260AB","#796FB3","#6A51A3","#5C3596",
#                        "#63439D","#54278F","#460D83","#4D1A89","#3F007D")
	} else if (cmap.type == 6) {
		mycol <- c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E",
				"#7A0177", "#49006A")
	} else if (cmap.type == 7) {
		mycol <- c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801",
				"#A63603", "#7F2704")
	} else if (cmap.type == 8) {
		mycol <- c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45",
				"#006D2C", "#00441B")
	} else if (cmap.type == 9) {
		mycol <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5",
				"#08519C", "#08306B")
	} else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 5)) {
		mycol <- vector(length=512, mode = "numeric")
		
		for (k in 1:256) {
			mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
		}
		for (k in 257:512) {
			mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
		}
		mycol <- rev(mycol)
	}
#	browser()
	ncolors <- length(mycol)
	pathway.mut = (-(pathway.mut*.749 + 0.251 - 1))
#	image(1:n.cols, 1, as.matrix(pathway.mut), col=gray(n.cols:0/n.cols))
	if (cmap.type == 5) {
		if (max.v == "NA") {
			max.v <- max(max(V1), -min(V1))
		}
		V2 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
		
	} else {
		V2 <- ceiling(ncolors * (V1 - min(V1))/(1.001*(max(V1) - min(V1))))
	}
	
	if (col.labels[1] == "NA") {      
		heatm <- matrix(0, nrow = n.rows + 1, ncol = n.cols)
		heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
		tot.cols <- ncolors
		if (legend == T && ifScoreBarPlot == FALSE) {
			nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(lcm(30), 1), 
					heights = c(10, 1), respect = FALSE)
		} else if( legend == FALSE && ifScoreBarPlot == FALSE ) {
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(lcm(30), 1), respect = FALSE)
		} else if( legend == TRUE && ifScoreBarPlot == TRUE ){
			nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(lcm(30), 1), 
					heights = c(10, 1), respect = FALSE)
		} else if( legend == FALSE && ifScoreBarPlot == FALSE){
			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(lcm(30), 1), respect = FALSE)
		}
		par(mar = c(3, 20, 3, 12))
		mycol <- c(mycol, phen.cmap[1:length(col.classes)])
		image(1:n.cols, 1:n.rows, t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=main, sub = sub, xlab= xlab, ylab=ylab)
		n.rows.phen <- 0
	} else {
		tot.cols <- ncolors
		if (is.vector(col.labels)) {
			heatm <- matrix(0, nrow = n.rows + 2, ncol = n.cols)
			heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
			n.rows.phen <- 1
			heatm[n.rows + 1,] <- tot.cols + col.labels
			cols.row <- length(unique(col.labels))
			tot.cols <- tot.cols + cols.row
			phen.cmap <- phen.cmap[1:cols.row]
			pathway.mut.grey = grey(pathway.mut)
			u.pathway.mut.grey = unique(pathway.mut.grey)
			heatm[n.rows + n.rows.phen,] = match(pathway.mut.grey, u.pathway.mut.grey) + tot.cols
		} else {
			n.rows.phen <- length(col.labels[,1])
			cols.row <- vector(length=n.rows.phen, mode = "numeric")
			heatm <- matrix(0, nrow = n.rows + n.rows.phen, ncol = n.cols)
			heatm[1:n.rows,] <- V2[seq(n.rows, 1, -1),]
#			heatm[n.rows+n.rows.phen,] = t(as.matrix(gray(pathway.mut)))
			for (k in seq(n.rows + n.rows.phen, n.rows + 1, -1)) {
				heatm[k,] <- tot.cols + col.labels[n.rows + n.rows.phen - k + 1,]
				cols.row[n.rows + n.rows.phen - k + 1] <- length(unique(col.labels[n.rows + n.rows.phen - k + 1,]))
				tot.cols <- tot.cols + cols.row[n.rows + n.rows.phen - k + 1]
#                 print(c("col:", k, ":", tot.cols + col.labels[n.rows + n.rows.phen - k + 1,], "tot.cols:", tot.cols))
				
			}
#			browser()
			pathway.mut.grey = grey(pathway.mut)
			u.pathway.mut.grey = unique(pathway.mut.grey)
			heatm[n.rows + n.rows.phen,] = match(pathway.mut.grey, u.pathway.mut.grey) + tot.cols
			tot.cols = tot.cols + length(u.pathway.mut.grey)
			phen.cmap <- phen.cmap[1:sum(unlist(lapply(col.classes, length)))]
		}
#		image(as.matrix(grey(pathway.mut)))
		if (legend == T && ifScoreBarPlot == FALSE) {
#              nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(10, 2), heights = c(6, 1), respect = FALSE)
			nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow=T), heights = c(9, 5, 1, 1), respect = FALSE)
		} else if( legend == FALSE && ifScoreBarPlot == FALSE ) {
			nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow=T), heights = c(6, 1, 1), respect = FALSE)
		} else if( legend == TRUE && ifScoreBarPlot == TRUE ){
			nf <- layout(matrix(c(1, 2, 3, 4, 0, 3), nrow=3, byrow=FALSE), 
					heights = c(9, 5, 1, 1), widths=c(10,1), respect = FALSE)
		} else if( legend == FALSE && ifScoreBarPlot == TRUE){
			nf <- layout(matrix(c(1, 2, 3, 4, 0, 3), ncol=2, byrow=FALSE), heights = c(6, 1, 1), 
					widths=c(10,1), respect = FALSE)
		}
		if(ifScoreBarPlot){
			par( mar = c(8, 20, 8, 0))
		} else{
			par( mar = c(3, 16, 3, 16))
		}
		#browser()
		
		mycol <- c(mycol, phen.cmap, u.pathway.mut.grey)
		if( length(tissue.names) > 1 ){
			#browser()
			tissue.colors = c(brewer.pal(12, "Set3"), brewer.pal(12,"Paired"))[1:n.tissues]
			if( !is.null(highlight.tissue.name)){
#				browser()
				highlight.tissue.index = which(tissue.names %in% highlight.tissue.name)
#			tissue.labels = ifelse(tissue.labels == highlight.tissue.index, tissue.labels, 1)
				tissue.colors[-highlight.tissue.index] = c(brewer.pal(12, "Set3"), 
						brewer.pal(12,"Paired"))[2]  # 2nd color is yellow and is nice to have as the default background
				# as the other colors are darker
				tissue.colors[highlight.tissue.index] = c(brewer.pal(12, "Set3"), 
						brewer.pal(12,"Paired"))[1:(length(highlight.tissue.name)+1)][-2]
			} 
#			print(matrix(c(tissue.names, tissue.colors), ncol=2), quote=FALSE)
			
			#row.names = c(row.names, "Tissue Types")
			mycol <- c(mycol, tissue.colors)
			n.rows.phen = n.rows.phen + 1
			heatm = rbind(heatm, (tissue.labels + tot.cols))
			tot.cols = tot.cols + length(tissue.colors)
			
		}
		#browser()
		
		
		#browser()
		image(1:n.cols, 1:(n.rows + n.rows.phen), t(heatm), zlim = c(0, tot.cols), col=mycol, axes=FALSE,
				main=main, sub = sub, xlab= xlab, ylab=ylab)
	}
	
# Add lines to separate phenotypes or subgroups
	
	if (col.labels2[1] != "NA") {
		groups <-  split(col.labels2, col.labels2)
		len.vec <- lapply(groups, length)
		plot.div <- c(0.51, cumsum(len.vec) + 0.5)
		for (i in plot.div) {
			lines(c(i, i), c(0, n.rows + n.rows.phen + 0.48), lwd = 2, cex = 0.9, col = "black")
		}
		lines(c(0.51, n.cols + 0.49), c(0.51, 0.51), lwd = 2, cex = 0.9, col = "black")
		lines(c(0.51, n.cols + 0.49), c(n.rows + n.rows.phen + 0.48, n.rows + n.rows.phen + 0.48), lwd = 2,
				cex = 0.9, col = "black")
		lines(c(0.51, n.cols + 0.49), c(n.rows + 0.50, n.rows + 0.50), lwd = 2,
				cex = 0.9, col = "black")
	}
	if (row.names[1] != "NA") {
#		browser()
		numC <- nchar(row.names)
		size.row.char <- char.rescale*25/(n.rows + 20)
		for (i in 1:n.rows) {
			row.names[i] <- substr(row.names[i], 1, 40)
			row.names[i] <- paste(row.names[i], " ", sep="")
		}
		if (phen.names[1] == "NA") {
			head.names <- paste("Class", seq(n.rows.phen, 1, -1))
		} else {
			head.names <- as.character(rev(phen.names))
		}
		row.names <- c(row.names[seq(n.rows, 1, -1)], head.names)
		if( length(tissue.names) > 1){ row.names = c(row.names, "Tissue Types")}
#            print(paste("n.rows:", n.rows))
#            print(paste("Phen names:", phen.names))
#            print(paste("Head names:", head.names))
#            print(paste("Row names:", row.names))
#		browser()
		axis(2, at=1:(n.rows + n.rows.phen), labels=row.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char,
				font.axis=2, line=-1)
	}
	
	if (row.names2[1] != "NA") {
#		browser()
		numC <- nchar(row.names2)
		size.row.char <- char.rescale*25/(n.rows + 20)
		for (i in 1:n.rows) {
			row.names2[i] <- substr(row.names2[i], 1, 40)
			row.names2[i] <- paste(" ", row.names2[i], sep="")
			
		}
		for( i in 1:n.rows.phen ){
			phen.names2[i] <- substr(phen.names2[i], 1, 40)
			phen.names2[i] <- paste( " ", phen.names2[i], sep="")
		}
		
		row.names2 <- rev(row.names2)
		phen.names2 <- rev(phen.names2)
		axis(4, at=1:(n.rows + n.rows.phen), labels=c(row.names2, phen.names2), adj= 0.5, tick=FALSE, las = 1, 
				cex.axis=size.row.char, font.axis=2, line=-1)
	}
	
	if (col.names[1] != "NA") {
		size.col.char <- char.rescale*20/(n.cols + 25)
		axis(1, at=1:n.cols, labels=col.names, tick=FALSE, las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
	}
	
	# Phenotype Legend 
	
#      print("--------------------------------------------------------------------------------------------")
	if (legend == T) {
		leg.txt <- NULL
		p.vec <- NULL
		c.vec <- NULL
		c2.vec <- NULL
		ind <- 1
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		for (i in 1:(n.rows.phen-1)) {  
#			browser()
			if (is.vector(col.labels)) {
				phen.v <- as.character(col.classes)
			} else {
				phen.v <- as.character(col.classes[[i]])
			}
			p.name <- paste(as.character(rev(head.names)[i]), ":   ", sep="")
			leg.txt <- c(p.name, phen.v)  
			p.vec <-  rep(22, cols.row[i] + 1)
			c.vec <-  c("#FFFFFF", phen.cmap[ind:(ind + cols.row[i] - 1)])
			c2.vec <- c("#FFFFFF", rep("black", cols.row[i]))
			ind <- ind + cols.row[i]
			offset <- 0.07
			legend(x=0, y= 1 - offset*i, 
					horiz = T, x.intersp = 0.5, legend=leg.txt, bty="n", xjust=0, yjust= 1, pch = p.vec, 
					pt.bg = c.vec, col = c2.vec, cex = 1.20, pt.cex=1.75)
		}
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		legend(x=0, y= 10, horiz = T, x.intersp = 0.5, legend=tissue.names, bty="n", xjust=0, yjust= 1, 
				fill = tissue.colors, cex = 1.20, pt.cex=1.75, ncol=1)
	}
	#browser()
	## Tissue Legend
	if(length(tissue.names)>1){
		#browser()
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		legend(x=0, y= 1, #horiz = T, x.intersp = 0.5, y.intersp=.25, 
				legend=tissue.names, bty="n", xjust=0, yjust= 1, 
				fill = tissue.colors, #cex = 1.20, pt.cex=1.75, 
				ncol=4)
	}
	# Color map legend
	
#       print(c("range V=", range(V)))
#       print(c("range V1=", range(V1)))
#       print(c("range V2=", range(V2)))
	
	par(mar = c(2, 12, 2, 12))
	num.v <- 20
	range.v <- range(V2)
	incr <-  (range.v[1] - range.v[2])/(num.v - 1)
	heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow=num.v, ncol=1)
	image(1:num.v, 1:1, heatm.v, zlim = c(0, tot.cols), col=mycol, axes=FALSE,
			main=" ", sub = " ", xlab= ylab, ylab=xlab)
	range.v <- range(V1)
	incr <-  (range.v[1] - range.v[2])/(num.v - 1)
	heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits=2), nrow=num.v, ncol=1)
#          print(c("heatm.v2=", heatm.v2))
	axis(3, at=1:num.v, labels=heatm.v2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.5*char.rescale, font.axis=1)
	
	
#	if(legend==TRUE){
#		par(mar = c(2, 12, 2, 12))
#		num.v <- 20
#		range.v <- range(V2)
#		incr <-  (range.v[1] - range.v[2])/(num.v - 1)
#		heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow=num.v, ncol=1)
#		image(1:num.v, 1:1, heatm.v, zlim = c(0, tot.cols), col=mycol, axes=FALSE,
#				main=" ", sub = " ", xlab= ylab, ylab=xlab)
#		range.v <- range(V1)
#		incr <-  (range.v[1] - range.v[2])/(num.v - 1)
#		heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits=2), nrow=num.v, ncol=1)
#		#          print(c("heatm.v2=", heatm.v2))
#		axis(3, at=1:num.v, labels=heatm.v2, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.5*char.rescale, font.axis=1)
#	}
	if(ifScoreBarPlot){
		par(mar = c(6, 5, 8, 1))
		barplot(rev(c(MI.list.phen, MI.list.model)), xlab="MI scores\nbarplot", ylab="", xlim=c(-1,1), 
				axes=TRUE, horiz=TRUE, axisnames=FALSE, 
#				main="MI scores\nbarplot", font.main=1
		)
#		browser()
#		text("MI scores\nbarplot")
	}
	return()
	
}

MSIG.Gct2Frame <- function(filename = "NULL") { 
#
# Reads a gene expression dataset in GCT format and converts it into an R data frame
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
	
	ds <- read.delim(filename, header=T, sep="\t", skip=2, row.names=1, 
			blank.lines.skip=T, comment.char="", as.is=T, na.strings = "")
	descs <- ds[,1]
	ds <- ds[-1]
	row.names <- row.names(ds)
	names <- names(ds)
	return(list(ds = ds, row.names = row.names, descs = descs, names = names))
}

Read.GeneSets.db <- function(
		gs.db,
		thres.min = 2,
		thres.max = 2000,
		gene.names = NULL)
{
	
	temp <- readLines(gs.db)
	max.Ng <- length(temp)
	temp.size.G <- vector(length = max.Ng, mode = "numeric") 
	for (i in 1:max.Ng) {
		temp.size.G[i] <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
	}
	max.size.G <- max(temp.size.G)      
	gs <- matrix(rep("null", max.Ng*max.size.G), nrow=max.Ng, ncol= max.size.G)
	temp.names <- vector(length = max.Ng, mode = "character")
	temp.desc <- vector(length = max.Ng, mode = "character")
	gs.count <- 1
	for (i in 1:max.Ng) {
		gene.set.size <- length(unlist(strsplit(temp[[i]], "\t"))) - 2
		gs.line <- noquote(unlist(strsplit(temp[[i]], "\t")))
		gene.set.name <- gs.line[1] 
		gene.set.desc <- gs.line[2] 
		gene.set.tags <- vector(length = gene.set.size, mode = "character")
		for (j in 1:gene.set.size) {
			gene.set.tags[j] <- gs.line[j + 2]
		}
		if (is.null(gene.names)) {
			existing.set <- rep(TRUE, length(gene.set.tags))
		} else {
			existing.set <- is.element(gene.set.tags, gene.names)
		}
		set.size <- length(existing.set[existing.set == T])
		if ((set.size < thres.min) || (set.size > thres.max)) next
		temp.size.G[gs.count] <- set.size
		gs[gs.count,] <- c(gene.set.tags[existing.set], rep("null", max.size.G - temp.size.G[gs.count]))
		temp.names[gs.count] <- gene.set.name
		temp.desc[gs.count] <- gene.set.desc
		gs.count <- gs.count + 1
	}
	Ng <- gs.count - 1
	gs.names <- vector(length = Ng, mode = "character")
	gs.desc <- vector(length = Ng, mode = "character")
	size.G <- vector(length = Ng, mode = "numeric") 
	
	gs.names <- temp.names[1:Ng]
	gs.desc <- temp.desc[1:Ng]
	size.G <- temp.size.G[1:Ng]
	
	return(list(N.gs = Ng, gs = gs, gs.names = gs.names, gs.desc = gs.desc, size.G = size.G, max.N.gs = max.Ng))
}

write.cls.2 <- function (class.v, phen, filename) 
{
	f <- file(filename, "w")
	n <- length(phen)
	l <- length(class.v)
	cat(l, n, "1", "\n", file = f, append = TRUE, sep = " ")
	cat("#", unlist(phen), "\n", file = f, append = TRUE, sep = " ")
	if (is.vector(class.v)) {
		class.v <- phen[class.v]
		cat(class.v, "\n", file = f, append = TRUE, sep = " ")
	} else {
		class.list <- matrix(0, nrow=length(class.v[,1]), ncol=length(class.v[1,]))
		for (i in 1:length(class.v[,1])) {
			class.list[i,] <- unlist(phen[[i]])[class.v[i,]]
			cat(class.list[i,], "\n", file = f, append = TRUE, sep = " ")
		}
	}
	close(f)
}

write.gct <- function(gct.data.frame, descs = "", filename) 
{
	f <- file(filename, "w")
	cat("#1.2", "\n", file = f, append = TRUE, sep = "")
	cat(dim(gct.data.frame)[1], "\t", dim(gct.data.frame)[2], "\n", file = f, append = TRUE, sep = "")
	cat("Name", "\t", file = f, append = TRUE, sep = "")
	cat("Description", file = f, append = TRUE, sep = "")
	
	names <- names(gct.data.frame)
	cat("\t", names[1], file = f, append = TRUE, sep = "")
	
	if (length(names) > 1) {
		for (j in 2:length(names)) {
			cat("\t", names[j], file = f, append = TRUE, sep = "")
		}
	}
	cat("\n", file = f, append = TRUE, sep = "\t")
	
	oldWarn <- options(warn = -1)
	m <- matrix(nrow = dim(gct.data.frame)[1], ncol = dim(gct.data.frame)[2] +  2)
	m[, 1] <- row.names(gct.data.frame)
	if (length(descs) > 1) {
		m[, 2] <- descs
	} else {
		m[, 2] <- row.names(gct.data.frame)
	}
	index <- 3
	for (i in 1:dim(gct.data.frame)[2]) {
		m[, index] <- gct.data.frame[, i]
		index <- index + 1
	}
	write.table(m, file = f, append = TRUE, quote = FALSE, sep = "\t", eol = "\n", col.names = FALSE, row.names = FALSE)
	close(f)
	options(warn = 0)
	
}
MSIG.ReadPhenFile <- function(file = "NULL") {
#
# Reads a matrix of class vectors from a CLS file and defines phenotype and class labels vectors
#  (numeric and character) for the samples in a gene expression file (RES or GCT format)
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
	
	cls.cont <- readLines(file)
	num.lines <- length(cls.cont)
	temp <- unlist(strsplit(cls.cont[[1]], " "))
	if (length(temp) == 3) {
		phen.names <- NULL
		col.phen <- NULL
	} else {
		l.phen.names <- match("phen.names:", temp)
		l.col.phen <- match("col.phen:", temp)
		phen.names <- temp[(l.phen.names + 1):(l.col.phen - 1)]
		col.phen <- temp[(l.col.phen + 1):length(temp)]
	}
	temp <- unlist(strsplit(cls.cont[[2]], " "))
	phen.list <- temp[2:length(temp)]
	
	for (k in 1:(num.lines - 2)) {
		temp <- unlist(strsplit(cls.cont[[k + 2]], " "))
		if (k == 1) {
			len <- length(temp)
			class.list <- matrix(0, nrow = num.lines - 2, ncol = len)
			class.v <- matrix(0, nrow = num.lines - 2, ncol = len)
			phen <- list(NULL)
		}
		class.list[k, ] <- temp
		classes <- unique(temp)
		class.v[k, ] <- match(temp, classes)
		phen[[k]] <- classes
	}
	if (num.lines == 3) {
		class.list <- as.vector(class.list)
		class.v <- as.vector(class.v)
		phen <- unlist(phen)
	}
	return(list(phen.list = phen.list, phen = phen, phen.names = phen.names, col.phen = col.phen,
					class.v = class.v, class.list = class.list))
}

MSIG.ReadPhenFile.2 <- function(file = "NULL") { 
#
# Reads a matrix of class vectors from a CLS file and defines phenotype and class labels vectors
#  (numeric and character) for the samples in a gene expression file (RES or GCT format)
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
	
	cls.cont <- readLines(file)
	num.lines <- length(cls.cont)
	temp <- unlist(strsplit(cls.cont[[1]], " "))
	if (length(temp) == 3) {
		phen.names <- NULL
		col.phen <- NULL
	} else {
#		browser()
		l.phen.names <- match("phen.names:", temp)
		l.col.phen <- match("col.phen:", temp)
		phen.names <- temp[(l.phen.names + 1):(l.col.phen - 1)]
		col.phen <- temp[(l.col.phen + 1):length(temp)]
	}
	temp <- unlist(strsplit(cls.cont[[2]], " "))
	phen.list <- temp[2:length(temp)]
	
	phen <- NULL
	for (k in 1:(num.lines - 2)) {
		temp <- unlist(strsplit(cls.cont[[k + 2]], " "))
		if (k == 1) {
			len <- length(temp)
			class.list <- matrix(0, nrow = num.lines - 2, ncol = len)
			class.v <- matrix(0, nrow = num.lines - 2, ncol = len)
#           phen <- NULL
		}
		class.list[k, ] <- temp
		classes <- unique(temp)
		class.v[k, ] <- match(temp, classes)
#        phen[[k]] <- classes
		phen <- c(phen, classes)
	}
	if (num.lines == 3) {
		class.list <- as.vector(class.list)
		class.v <- as.vector(class.v)
#         phen <- unlist(phen)
	}
	return(list(phen.list = phen.list, phen = phen, phen.names = phen.names, col.phen = col.phen,
					class.v = class.v, class.list = class.list))
}


MSIG.Subset.Dataset.2 <- function(
		input.ds,
		input.cls = NULL,
		column.subset = "ALL",    # subset of column numbers or names (or phenotypes)
		column.sel.type = "samples",  # "samples" or "phenotype"
		row.subset = "ALL",       # subset of row numbers or names
		output.ds,
		output.cls = NULL) {
	
# start of methodology
	
	print(c("Running MSIG.Subset.Dataset... on GCT file:", input.ds))
	print(c("Running MSIG.Subset.Dataset... on CLS file:", input.cls))
	
# Read input datasets
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)
	m <- data.matrix(dataset$ds)
	gs.names <- dataset$row.names
	gs.descs <- dataset$descs
	sample.names <- dataset$names
	
# Read CLS file
	
	if (!is.null(input.cls)) {
		CLS <- MSIG.ReadPhenFile.2(file=input.cls)
		class.labels <- CLS$class.v
		class.phen <- CLS$phen
		class.list <- CLS$class.list 
	}
	
# Select desired column subset
	
	if (column.sel.type == "samples") {
		if (column.subset[1] == "ALL") {
			m2 <- m
			sample.names2 <- sample.names
			if (!is.null(input.cls)) {
				class.labels2 <- class.labels
			}
		} else {
			if (is.numeric(column.subset[1])) {
				m2 <- m[,column.subset]
				sample.names2 <- sample.names[column.subset]
				if (!is.null(input.cls)) {
					if (is.vector(class.labels)) {
						class.labels2 <- class.labels[column.subset]
					} else {
						class.labels2 <- class.labels[, column.subset]
					}
				}
			} else {
				locations <- !is.na(match(sample.names, column.subset))
				sample.names2 <- sample.names[locations]
				m2 <- m[, locations]
				if (!is.null(input.cls)) {
					if (is.vector(class.labels)) {
						class.labels2 <- class.labels[locations]
					} else {
						class.labels2 <- class.labels[, locations]
					}
				}
			}
		}
	} else if (column.sel.type == "phenotype") {
		locations <- !is.na(match(class.list, column.subset))
		sample.names2 <- sample.names[locations]
		m2 <- m[,locations]
		if (!is.null(input.cls)) {
			if (is.vector(class.labels)) {
				class.labels2 <- class.labels[locations]
			} else {
				class.labels2 <- class.labels[, locations]
			}
		}
	}
	
	if (row.subset[1] == "ALL") {
		m3 <- m2
		gs.names2 <- gs.names
		gs.descs2 <- gs.descs
	} else {
		locations <- !is.na(match(gs.names, row.subset))
		m3 <- m2[locations,]
		gs.names2 <- gs.names[locations]
		gs.descs2 <- gs.descs[locations]
	}
	
# Save datasets
	
	V <- data.frame(m3)
	names(V) <- sample.names2
	row.names(V) <- gs.names2
	write.gct(gct.data.frame = V, descs = gs.descs2, filename = output.ds)  
	
	if (!is.null(input.cls)) {
		write.cls.2(class.v = class.labels2, phen = class.phen, filename = output.cls) 
	}
}

OPAM.match.projection.to.pathway  <- function(
		input.ds,
		input.cls          = NA,
		results.dir,
		normalize.score    = F,
		normalization.type = "zero.one",
		pathway,
		max.n              = 10,
		user.colors        = NA,
		decreasing.order   = T,
		sort.columns       = F,
		char.rescale       = 1.25,
		cmap.type          = 3,
		row.norm           = T,
		output.dataset     = NA)
{
	library(gtools)
	library(verification)
	library(ROCR)
	library(MASS)
	library(RColorBrewer)
	library(heatmap.plus)
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	pathway.names <- dataset$row.names
	pathway.descs <- dataset$descs
	Ns <- length(m[1,])
	dim(m)
	sample.names <- dataset$names
	
	n.pathways <- length(m[,1])
	temp <- strsplit(input.ds, split="/") # Extract test file name
	s <- length(temp[[1]])
	test.file.name <- temp[[1]][s]
	temp <- strsplit(test.file.name, split=".gct")
	test.file.prefix <-  temp[[1]][1]
#   char.res <-  0.013 * n.pathways + 0.65
	
	# normalize scores
	
	if (normalize.score == T) {
		if (normalization.type == "zero.one") {
			for (i in 1:n.pathways) {
				m[i,] <- (m[i,] - min(m[i,]))/(max(m[i,]) - min(m[i,]))
			}
		} else if (normalization.type == "z.score") {
			for (i in 1:n.pathways) {
				m[i,] <- (m[i,] - mean(m[i,]))/sd(m[i,])
			}
		} else if (normalization.type == "r.z.score") {
			for (i in 1:n.pathways) {
				m[i,] <- (m[i,] - median(m[i,]))/mad(m[i,])
			}
		}         
	}
	
	
	loc <- match(pathway, pathway.names)
	print(c("loc:", loc))
	if (sort.columns == T) {
		s.order <- order(m[loc,], decreasing = decreasing.order)
		m2 <- m[, s.order]
		sample.names2 <- sample.names[s.order]
	} else {
		m2 <- m
		sample.names2 <- sample.names
	}
	correl <- cor(t(m2))[, loc]
	m.order <- order(correl, decreasing=T)
	correl2 <- correl[m.order]
	m2 <- m2[m.order[1:max.n],]
	pathway.names2 <- pathway.names[m.order]
	pathway.descs2 <- signif(correl2, digits=3)
	
	if (input.cls == "NA") {
		cls.labels2 <- c(rep(0, 10), rep(1, length(sample.names2) - 10))
		cls.phen2 <- c(" ")
		colors.list <- c("white")
		phen.names2 <- "    "
	} else {
		CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
		cls.labels <- CLS$class.v
		cls.phen <- CLS$phen
		cls.list <- CLS$class.list 
		if (!is.null(CLS$phen.names)) {
			phen.names <- CLS$phen.names
		} else {
			phen.names <- "  "
		}
		if (is.vector(cls.labels)) {
			if (sort.columns == T) {
				cls.labels2 <- cls.labels[s.order]
				cls.list2 <- cls.list[s.order]
			} else {
				cls.labels2 <- cls.labels
				cls.list2 <- cls.list
			}
			n.phen <- 1
		} else {
			if (sort.columns == T) {
				cls.labels2 <- cls.labels[, s.order]
				cls.list2 <- cls.list[, s.order]
			} else {
				cls.labels2 <- cls.labels
				cls.list2 <- cls.list
			}
			n.phen <- length(cls.labels2[,1])
		}
		cls.phen2 <- list(NULL)
		if (is.vector(cls.labels2)) {
			classes <- unique(cls.list2)
			cls.phen2 <- classes
			cls.labels2 <- match(cls.list2, cls.phen2)
		} else {
			for (kk in 1:length(cls.list2[, 1])) {
				classes <- unique(cls.list2[kk,])
				cls.phen2[[kk]] <- classes
				cls.labels2[kk,] <- match(cls.list2[kk,], cls.phen2[[kk]])
			}
		}
		phen.names2 <- phen.names
		if (!is.na(user.colors[1])) {
			c.test <- user.colors
		} else {
			if (!is.null(CLS$col.phen)) {
				c.test <- CLS$col.phen
			} else {
				c.test <- c(brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Set1"),
						brewer.pal(n=8, name="Accent"), brewer.pal(n=9, name="Spectral"), brewer.pal(n=8, name="Set3"), 
						brewer.pal(n=8, name="BuGn"),
						brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
						brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
						brewer.pal(n=8, name="BuGn"),
						brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
						brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
						brewer.pal(n=8, name="BuGn"))
			}
		}
	}
	cls.phen.index <- unlist(cls.phen2)
	colors.list <- c.test[1:length(cls.phen.index)]
	
	filename <- paste(results.dir, test.file.prefix, ".SORT.PROJ.TO.", pathway, sep="")
	pdf(file=paste(filename, ".pdf", sep=""), height = 8.5, width = 10.5)
	
	MSIG.HeatMapPlot.7(V = m2, row.names = pathway.names2[1:max.n],
			row.names2 = pathway.descs2[1:max.n], col.labels = cls.labels2, 
			col.classes = cls.phen2, phen.cmap = colors.list, phen.names = phen.names,
			col.names = sample.names2, main = " ", xlab="  ", ylab="  ", row.norm = row.norm,  
			cmap.type = cmap.type, char.rescale = char.rescale,  legend=T)
	dev.off()
	
	if (!is.na(output.dataset)) {
		V.GCT <- m2
		colnames(V.GCT) <- sample.names2
		row.names(V.GCT) <- pathway.names2
		write.gct(gct.data.frame = V.GCT, descs = pathway.descs2, filename =output.dataset)  
	}
	
}

MSIG.Define.Dataset.from.Table2 <- function(
		input.gct,
		table.txt,
		output.gct,
		output.txt = NULL,  # optional version of table with overlap (GCT & TAB) samples
		output.cls,
		prefix_entries = F,
		remove_samples_with_na = TRUE,
		remove_na_of_specific_gene = "")
{
# Read input dataset
	
	library(RColorBrewer)
	
	dataset1 <- MSIG.Gct2Frame(filename = input.gct)
	m <- data.matrix(dataset1$ds)
	gene.names <- dataset1$row.names
	gene.decs  <- dataset1$descs
	sample.names.gct <- dataset1$names
	Ns <- length(sample.names.gct)
	
# Read Table 
	
	tab <- read.delim(table.txt, header=T, row.names = 1, 
			sep="\t", skip=0, blank.lines.skip=T, comment.char="", as.is=T)
	sample.names.tab <- row.names(tab)
	phen.names <- names(tab)
#	browser()
	
	
	if(remove_samples_with_na){
		#browser()
		ind.na = lapply(tab, function(x) which( is.na(x)))
		#browser()
		if(remove_na_of_specific_gene != ""){
			ind.na.u = ind.na[which(phen.names == remove_na_of_specific_gene)][[1]]
		} else{
			ind.na.u = unique(unlist(unique(ind.na)))
		}
		if(sum(unlist(lapply(ind.na.u, length))) > 0){
			sample.names.tab = sample.names.tab[-ind.na.u]
			tab = tab[-ind.na.u,]
		}
	}
	
	overlap <- intersect(sample.names.tab, sample.names.gct)
	print("sample names GCT")
	print(sample.names.gct)
	print("sample names TAB")
	print(sample.names.tab)
	
	
	locs.gct <- match(overlap, sample.names.gct)
	print(match(sample.names.tab, sample.names.gct))
	print(match(sample.names.gct, sample.names.tab))
	
	locs.tab <- match(overlap, sample.names.tab)
	print(locs.tab)
	print(c("GCT matching set (", length(locs.gct), " samples):", sample.names.gct[locs.gct]))
	print(c("TAB matching set (", length(overlap), " samples):", sample.names.tab[locs.tab]))
	print(c("overlap set (", length(overlap), " samples):", overlap))
	
	m2 <- m[, locs.gct]
	sample.names.gct <- sample.names.gct[locs.gct]
	sample.names.tab <- sample.names.tab[locs.tab]
	
	if (!is.null(output.txt)) {
		tab2 <- tab[locs.tab,]
		sample.names.tab2 <- sample.names.tab[locs.tab]
		col.names <- paste(colnames(tab2), collapse = "\t")
		col.names <- paste("SAMPLE", col.names, sep= "\t")
		write(noquote(col.names), file = output.txt, append = F, ncolumns = length(col.names))
		write.table(tab2, file=output.txt, quote=F, col.names = F, row.names = T, append = T, sep="\t")
	}
	#browser()
	cls.table <- t(tab[locs.tab,])
	
	
	
	if (prefix_entries == TRUE) {
		for (i in 1:length(cls.table[,1])) {
#        cls.table[i,] <- paste(row.names(cls.table)[i], cls.table[i,], sep=".")
			cls.table[i,] <- paste(colnames(tab)[i], tab[,i], sep=".")
		}
	}
	## Get rid of leading whitespace
	cls.table = apply(cls.table, MARGIN=2, FUN=function(x) gsub("^\\s+", "", x, perl=TRUE))
	
	
	if (!is.null(output.gct)) {      
		V <- data.frame(m2)
		names(V) <- sample.names.gct
		row.names(V) <- gene.names
		write.gct(gct.data.frame = V, descs = gene.decs, filename = output.gct)
	}
	
	class.phen <- unique(cls.table)
	n <- length(class.phen)
	l <- length(cls.table[1,])
	
	col.list <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
			brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
			brewer.pal(n=8, name="BuGn"),
			brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
			brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
			brewer.pal(n=8, name="BuGn"),
			brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
			brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
			brewer.pal(n=8, name="BuGn"))
	num <- 0
	class.order.list <- NULL
	for (i in 1:length(cls.table[,1])) {
		num <- num + length(unique(cls.table[i,]))
		class.order.list <- c(class.order.list, unique(cls.table[i,]))
		#browser()
	}
	
	phen.names.string <- paste("phen.names:", paste(phen.names, collapse=" "), sep=" ")
	sig.col <- col.list[1:num]
	col.phen.string <- paste("col.phen:", paste(sig.col, collapse=" "), sep=" ")
	
	browser()
	cat(paste(l, num, length(cls.table[, 1]), phen.names.string, col.phen.string, sep=" "), "\n", 
			file = output.cls, append = FALSE, sep = "")
	cat("# ", paste(class.order.list, collapse=" "), "\n", file = output.cls, append = TRUE, sep = "")
	for (i in 1:length(cls.table[,1])) {
		cat(paste(cls.table[i,], collapse=" "), "\n", file = output.cls, append = TRUE, sep = "")
	}
}

MSIG.Define.Dataset.from.Table.2 <- function(
		input.gct,
		table.txt,
		output.gct,
		output.cls,
		prefix_entries = F)
{
# Read input dataset
	
	library(RColorBrewer)
	
	dataset1 <- MSIG.Gct2Frame(filename = input.gct)
	m <- data.matrix(dataset1$ds)
	gene.names <- dataset1$row.names
	gene.decs  <- dataset1$descs
	sample.names.gct <- dataset1$names
	Ns <- length(sample.names.gct)
	
#	browser()
# Read Table 
	
	
	tab <- read.delim(table.txt, header=T, row.names = 1, 
			sep="\t", skip=0, blank.lines.skip=T, comment.char="", as.is=T)
	sample.names.tab <- row.names(tab)
	phen.names <- names(tab)
	overlap <- intersect(sample.names.tab, sample.names.gct)
	if(length(overlap)==0){ return(NULL)}
#	print(overlap)
#	print("sample names GCT")
#	print(sample.names.gct)
#	print("sample names TAB")
#	print(sample.names.tab)
	
	if(length(overlap)==0){ return(NULL)}
	
	locs.gct <- match(overlap, sample.names.gct)
	print(match(sample.names.tab, sample.names.gct))
	print(match(sample.names.gct, sample.names.tab))
	locs.tab <- match(overlap, sample.names.tab)
#	print(locs.tab)
#	print(c("GCT matching set (", length(locs.gct), " samples):", sample.names.gct[locs.gct]))
#	print(c("TAB matching set (", length(overlap), " samples):", sample.names.tab[locs.tab]))
#	print(c("overlap set (", length(overlap), " samples):", overlap))
	
	
	
	m2 <- m[, locs.gct]
	sample.names.gct <- sample.names.gct[locs.gct]
	sample.names.tab <- sample.names.tab[locs.tab]
	cls.table <- t(tab[locs.tab,])
	
	if (prefix_entries == TRUE) {
		for (i in 1:length(cls.table[,1])) {
#        cls.table[i,] <- paste(row.names(cls.table)[i], cls.table[i,], sep=".")
			cls.table[i,] <- paste(colnames(tab)[i], tab[,i], sep=".")
		}
	}
	
	if (!is.null(output.gct)) {      
		V <- data.frame(m2)
		names(V) <- sample.names.gct
		row.names(V) <- gene.names
		write.gct(gct.data.frame = V, descs = gene.decs, filename = output.gct)
	}
	
	class.phen <- unique(cls.table)
	n <- length(class.phen)
	l <- length(cls.table[1,])
	
	col.list <- c(brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
			brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
			brewer.pal(n=8, name="BuGn"),
			brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
			brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
			brewer.pal(n=8, name="BuGn"),
			brewer.pal(n=7, name="Set2"), brewer.pal(n=7, name="Dark2"), brewer.pal(n=7, name="Set1"),
			brewer.pal(n=8, name="Accent"), brewer.pal(n=10, name="Spectral"), brewer.pal(n=8, name="Set3"), 
			brewer.pal(n=8, name="BuGn"))
#	num <- 0
#	class.order.list <- NULL
	class.order.list = apply(cls.table, 1, unique)
	num = sum(unlist(lapply(class.order.list, length)))
	class.order.list = unlist(class.order.list)
	
	
#	for (i in 1:length(cls.table[,1])) {
	##		num <- num + length(unique(cls.table[i,]))
#		class.order.list <- c(class.order.list, unique(cls.table[i,]))
#	}
	
	phen.names.string <- paste("phen.names:", paste(phen.names, collapse=" "), sep=" ")
	sig.col <- col.list[1:num]
	col.phen.string <- paste("col.phen:", paste(sig.col, collapse=" "), sep=" ")
	cat(paste(l, num, length(cls.table[, 1]), phen.names.string, col.phen.string, sep=" "), "\n", 
			file = output.cls, append = FALSE, sep = "")
	cat("# ", paste(class.order.list, collapse=" "), "\n", file = output.cls, append = TRUE, sep = "")
	for (i in 1:length(cls.table[,1])) {
		cat(paste(cls.table[i,], collapse=" "), "\n", file = output.cls, append = TRUE, sep = "")
	}
}


rec.area <- function(
		obs,
		pred,
		metric = "absolute.deviation",       # Either "squared.error" or "absolute.deviation"
#		null.distribution = "gaussian",  
		# Either "gaussian" [null.model = mean(obs)] or "laplacian" [null.model = median(obs)]
		interval = 0.01
){
#	browser()
	error.windows = seq(0, 1, by=interval)
	n.errors = length(error.windows)
	intervals = rep( interval, n.errors )
	n.obs = length(obs)
	n.pred = length(pred)
	if( n.obs != n.pred ){ stop( "The number of observations does not equal the number of predictions." ) }
#	if( null.distribution == "gaussian" ){ null.model = mean(obs) 
#	} else if( null.distribution == "laplacian" ){ null.model = median(obs) }
	
	if( metric == "squared.error" ){
		difference = (obs-pred)^2
		accuracy = unlist(lapply(error.windows, FUN=squared.error, difference, n.obs))
	} else if( metric == "absolute.deviation" ){
		difference = abs(obs-pred)
		accuracy = unlist(lapply(error.windows, FUN=absolute.deviation, difference, n.obs))
	}
#	plot(accuracy, type="l"); par(new=TRUE); plot(error.windows, type="l")
	
	triangle.heights = accuracy - c(0, accuracy[1:(n.errors-1)])
	triangles = triangle.heights*intervals/2
	rectangle.heights = c(0, accuracy[1:(n.errors-1)])
	rectangles = rectangle.heights*intervals
#	A = (cumsum(accuracy)*intervals)[n.errors]
	A = sum( rectangles + triangles)
	
	# Calculate p-value using Kolmogorov-Smirnov Test
#	Dn = max(accuracy-error.windows)
#	i = 1:100
#	x = sqrt( n.obs*n.pred/(n.obs+n.pred) )*Dn
#	p.value = 1 - (sqrt(2*pi)/x)*sum( exp(-(2*i - 1)^2 * pi^2/ (8*x^2)) )
#	browser()
#	pred.scrambled = sample(pred)
#	difference.scrambled = abs(pred.scrambled - obs)
#	accuracy.scrambled = unlist(lapply(error.windows, FUN=squared.error, difference.scrambled, n.obs))
#	triangle.heights = accuracy.scrambled - c(0, accuracy.scrambled[1:(n.errors-1)])
#	triangles = triangle.heights*intervals/2
#	rectangle.heights = c(0, accuracy.scrambled[1:(n.errors-1)])
#	rectangles = rectangle.heights*intervals
#	A.scrambled = sum( rectangles + triangles)
#	T2.scrambled = .5*(sum((accuracy.scrambled-error.windows)^2))
#	p.value.scrambled = cvmts.pval(T2.scrambled, n.errors, n.errors)
	
	# Calculate p-value using Cramer-Von-Mises Criterion
	T2 = .25*(sum((accuracy-error.windows)^2))  # accuracy-error.windows = integral difference between null model and REC
#	browser()
	p.value = cvmts.pval(T2, n.errors, n.errors)
	
	
#	T2.norm = (T2- min(T2))/(max(T2)-min(T2))
#	wilcox.test(T2.norm, error.)
#browser()
#	U =  2*n.errors^2*(T2 + (4*n.errors^2-1)/12*n.errors)
#	p.value.u = cvmts.pval(U, n.errors, n.errors)
	
#	print('calculating REC...')
#	browser()
#	rec.list.ccle[master.ind] <<- A
#	p.value.list.ccle[master.ind] <<- p.value
#	T2.list.ccle[master.ind] <<- T2
#	
#	
#	
#	rec.list.scrambled[master.ind] <<- A.scrambled
#	p.value.list.scrambled[master.ind] <<- p.value.scrambled
#	T2.list.scrambled[master.ind] <<- T2.scrambled
#	master.ind <<- master.ind+1
#	browser()
#	stat = ks.test(accuracy, error.windows, exact=TRUE)
	return( list(A = A, p.value = p.value, T2=T2) )
}

squared.error <- function( error, squared.difference, n ){
	return( length(which( squared.difference <= error ))/n )
}

absolute.deviation <- function( error, absolute.difference, n ) {
	return( length(which( absolute.difference <= error ))/n )
}

#mutual.inf <- function(x, y, n.grid=100) {
#	
#	kde2d.xy <- kde2d(x, y, n = n.grid, h = c(width.SJ(x, method="dpi"), width.SJ(y, method="dpi")))
##	X <- kde2d.xy$x
##	Y <- kde2d.xy$y
#	PXY <- kde2d.xy$z/sum(kde2d.xy$z)
#	
#	PX <- apply(PXY, MARGIN=1, sum)
#	PX <- PX/sum(PX)
#	PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
#	
#	PY <- apply(PXY, MARGIN=2, sum)
#	PY <- PY/sum(PY)
#	PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)
#	
#	MI <- sum(PXY * log2(PXY/(PX*PY)))
##	browser()
#	return(MI)
#}

mutual.inf.2 <- function(x, y, n.grid=20, normalize.by ="HXY", # Whether to normalize by HXY, HX, or HY
		pos.and.neg = T

) {
	# x and y can be binary or continuous
	# If there is not sufficient variation in x and y, 
	# will take the standard deviation as the bandwidth
	# (IQR finds the inter-quartile range of the data vector)
	
	## Add random noise if vectors are constant to prevent errors in bcv 
	## (because the inter-quartile range of a constant vector is 0, and the IQR is 
	## used in calculating the bandwidth)
#	if( length(unique(x)) == 1 ){
	##		browser()
#		x = x + (10e-10 * runif(length(x)))
#	}
#	if( length(unique(y)) == 1 ){
#		y = y + (10e-10 * runif(length(y)))
#	}
	y = unlist(y); y = as.vector(y)
	x = unlist(x); x = as.vector(x)
	if( length(unique(y)) == 1 || length(unique(x)) == 1){ return(0) }
	if( sd(y) == 0){
		y = y + runif(n=length(y), min=mean(y)-0.001, max=mean(y)+0.001)
	}
	if( sd(x) == 0){
		x = x + runif(n=length(x), min=mean(x)-0.001, max=mean(x)+0.001)
	}
#	bandwidth.x = ifelse(IQR(x) == 0, bcv(x, n.grid), width.SJ(x, method="dpi"))
#	bandwidth.y = ifelse(IQR(y) == 0, bcv(y, n.grid), width.SJ(y, method="dpi"))
#	print("---")
#	print(x)
#	print(y)
	
	## Using suppressWarnings because bcv(.) gets mad if the minimum is on one side of the vector
#	before.kde2d[kde2d.i] <<- proc.time()[3]
	kde2d.xy <- kde2d(x, y, n = n.grid, h = c(suppressWarnings(bcv(x)),
		suppressWarnings(bcv(y))) )
#	after.kde2d[kde2d.i] <<- proc.time()[3]
#	kde2d.i <<- kde2d.i + 1
#	X <- kde2d.xy$x
#	Y <- kde2d.xy$y
#	Z = kde2d.xy$z
	PXY <- kde2d.xy$z/sum(kde2d.xy$z)
	
	PX <- rowSums(PXY)#apply(PXY, MARGIN=1, sum)
	PX <- PX/sum(PX)
	HX = -sum(PX * log2(PX))
	PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
	
	PY <- colSums(PXY)#apply(PXY, MARGIN=2, sum)
	PY <- PY/sum(PY)
	HY = -sum( PY * log2(PY))
	PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)
	
#	browser()
	
#	MIXY = PXY * log2(PXY/(PX*PY))
#	
#	if( pos.and.neg ){
#	q1 = MIXY[1:(n.grid/2), 1:(n.grid/2)]
#	q2 = MIXY[1:(n.grid/2), (n.grid/2 + 1):n.grid]
#	q3 = MIXY[(n.grid/2+1):n.grid, 1:(n.grid/2)]
#	q4 = MIXY[(n.grid/2+1):n.grid, (n.grid/2+1):n.grid]
#	
#	# q's divide MIXY into quarters. If the sum of q2 and q3 is greater than the sum of q1 and q4, then
#	# x and y are negatively correlated.
#	# on heatmap:   q2  q4
#	#               q1  q3
#	
	## Ignore NaN's that are a result of underflow (experimentally derived)
#	MI <- ifelse( sum(q1+q4, na.rm=TRUE) < sum(q2+q3, na.rm=TRUE), 
#			-sum(MIXY, na.rm=TRUE), sum(MIXY, na.rm=TRUE))
#} else{ MI = sum(MIXY, na.rm=TRUE)}
#	MI <- ifelse( sum(q1+q4, na.rm=TRUE) < sum(q2+q3, na.rm=TRUE), 
# -sum(q2+q3-q1-q4, na.rm=TRUE), sum(q1+q4-q2-q3, na.rm=TRUE))
	HXY <- - sum(PXY * log2(PXY), na.rm=TRUE)
	
#	HX = -sum( PX * log2(PX) )
#	HY = -sum( PY * log2(PY) )
#	MI.norm = (HX+HY)/HXY
#browser()
#	normalization.factor = 1 #ifelse(normalize.by=="HXY", HXY, ifelse(normalize.by=="HX", HX, HY))
	## browser()
	MI.norm =  ifelse(pos.and.neg, sign(cor(x, y)), 1) * ((HX + HY)/HXY - 1) #MI/normalization.factor
#	browser()
	
	return( MI.norm )#list(MI=MI, HXY=HXY))
}

mutual.inf.2.kde2d.olga <- function(x, y, n.grid=20, normalize.by ="HXY", # Whether to normalize by HXY, HX, or HY
		pos.and.neg = T

) {
	# x and y can be binary or continuous
	# If there is not sufficient variation in x and y, 
	# will take the standard deviation as the bandwidth
	# (IQR finds the inter-quartile range of the data vector)
	
	## Add random noise if vectors are constant to prevent errors in bcv 
	## (because the inter-quartile range of a constant vector is 0, and the IQR is 
	## used in calculating the bandwidth)
#	if( length(unique(x)) == 1 ){
	##		browser()
#		x = x + (10e-10 * runif(length(x)))
#	}
#	if( length(unique(y)) == 1 ){
#		y = y + (10e-10 * runif(length(y)))
#	}
	y = unlist(y); y = as.vector(y)
	x = unlist(x); x = as.vector(x)
	if( length(unique(y)) == 1 || length(unique(x)) == 1){ return(0) }
	if( sd(y) == 0){
		y = y + runif(n=length(y), min=mean(y)-0.001, max=mean(y)+0.001)
	}
	if( sd(x) == 0){
		x = x + runif(n=length(x), min=mean(x)-0.001, max=mean(x)+0.001)
	}
#	bandwidth.x = ifelse(IQR(x) == 0, bcv(x, n.grid), width.SJ(x, method="dpi"))
#	bandwidth.y = ifelse(IQR(y) == 0, bcv(y, n.grid), width.SJ(y, method="dpi"))
#	print("---")
#	print(x)
#	print(y)
	
	## Using suppressWarnings because bcv(.) gets mad if the minimum is on one side of the vector
	before.kde2d[kde2d.i] <<- proc.time()[3]
	kde2d.xy <- kde2d.olga(x, y, n = n.grid, h = c(suppressWarnings(bcv(x)), suppressWarnings(bcv(y))) )
	after.kde2d[kde2d.i] <<- proc.time()[3]
	kde2d.i <<- kde2d.i + 1
#	X <- kde2d.xy$x
#	Y <- kde2d.xy$y
#	Z = kde2d.xy$z
	PXY <- kde2d.xy$z/sum(kde2d.xy$z)
	
	PX <- rowSums(PXY)#apply(PXY, MARGIN=1, sum)
	PX <- PX/sum(PX)
	HX = -sum(PX * log2(PX))
	PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
	
	PY <- colSums(PXY)#apply(PXY, MARGIN=2, sum)
	PY <- PY/sum(PY)
	HY = -sum( PY * log2(PY))
	PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)
	
#	browser()
	
#	MIXY = PXY * log2(PXY/(PX*PY))
#	
#	if( pos.and.neg ){
#	q1 = MIXY[1:(n.grid/2), 1:(n.grid/2)]
#	q2 = MIXY[1:(n.grid/2), (n.grid/2 + 1):n.grid]
#	q3 = MIXY[(n.grid/2+1):n.grid, 1:(n.grid/2)]
#	q4 = MIXY[(n.grid/2+1):n.grid, (n.grid/2+1):n.grid]
#	
#	# q's divide MIXY into quarters. If the sum of q2 and q3 is greater than the sum of q1 and q4, then
#	# x and y are negatively correlated.
#	# on heatmap:   q2  q4
#	#               q1  q3
#	
	## Ignore NaN's that are a result of underflow (experimentally derived)
#	MI <- ifelse( sum(q1+q4, na.rm=TRUE) < sum(q2+q3, na.rm=TRUE), 
#			-sum(MIXY, na.rm=TRUE), sum(MIXY, na.rm=TRUE))
#} else{ MI = sum(MIXY, na.rm=TRUE)}
#	MI <- ifelse( sum(q1+q4, na.rm=TRUE) < sum(q2+q3, na.rm=TRUE), 
# -sum(q2+q3-q1-q4, na.rm=TRUE), sum(q1+q4-q2-q3, na.rm=TRUE))
	HXY <- - sum(PXY * log2(PXY), na.rm=TRUE)
	
#	HX = -sum( PX * log2(PX) )
#	HY = -sum( PY * log2(PY) )
#	MI.norm = (HX+HY)/HXY
#browser()
#	normalization.factor = 1 #ifelse(normalize.by=="HXY", HXY, ifelse(normalize.by=="HX", HX, HY))
	## browser()
	MI.norm =  ifelse(pos.and.neg, sign(cor(x, y)), 1) * ((HX + HY)/HXY - 1) #MI/normalization.factor
#	browser()
	
	return( MI.norm )#list(MI=MI, HXY=HXY))
}

mutual.inf = mutual.inf.2

#mutual.inf.2.single.gene.target <- function( signature, gene.target){
#	return(mutual.inf.2(gene.target, signature))
#}

mutual.inf.2.multiple.gene.targets <- function( signature, gene.targets ){
	return(apply(gene.targets, 
					MARGIN=1, FUN=mutual.inf.2, 
					signature ) )
}

mutual.inf.3 <- function(gene.target, signature.matrix, signature.indices, n.grid=100, gene.target.name = "",
		n.randomizations = 100, tissue = "NA") {
	## How this is different from mutual.inf.2:
	## x is a target pathway and y is a matrix
	## calculates a false discovery rate for the mutual information 
	## of the pathway with a randomized version of x
	## y.ind indicates the row index of the target gene set / signature
	
	# x and the elements of y can be binary or continuous
	# If there is not sufficient variation in x and y, 
	# will take the standard deviation as the bandwidth
	# (IQR finds the inter-quartile range of the data vector)
	
	# "signature.indices" is the indices of signatures that you are interested in using.
	# This code is used in comparing the chosen signatures to SUMMARY, the mutation
	# status of all the cell lines. 
#	browser()
	n.signatures = length(signature.matrix[,1])
	MI.vector = vector(length=n.signatures, mode="double")
#	MI.vector.rand = vector(length=length(signature.matrix[,1]), mode="double")
	gene.target.rand = t(replicate(n.randomizations, sample(gene.target)))
	MI.matrix.rand = matrix(ncol = n.signatures, nrow = n.randomizations)
#	browser()
#	for( i in 1:length(signature.matrix[,1]) ){
	MI.vector = apply(signature.matrix, MARGIN=1, FUN=mutual.inf.2, gene.target)
#		browser()
#		for( j in 1:n.randomizations ){
	MI.matrix.rand = apply(signature.matrix, 
			MARGIN=1, FUN=mutual.inf.2.multiple.gene.targets, 
			gene.target.rand)
	#mutual.inf.2(gene.target.rand[j,], signature.matrix[i,])$MI
#		}
#		browser()
#		MI.vector.rand[i] = mean(temp.MI.rand)
#	}
#	x.rand = sample(x)
#	print("Make plot of densities!! And save the output!")
#	browser()
	
	quartz()
#	if( gene.target.name =="SUMMARY"){
	temp <- density(MI.vector, adjust=1, n = 512, from=min(MI.vector), to=max(MI.vector))
	x <- temp$x
	y <- temp$y/sum(temp$y)
	
	temp.rand <- density(MI.matrix.rand, adjust=1, n = 512, from=min(MI.matrix.rand), to=max(MI.matrix.rand))
	x.rand <- temp.rand$x
	y.rand <- temp.rand$y/sum(temp.rand$y)
#		pdf(file=paste(tissue, gene.target.name, n.randomizations, "pdf", sep=".") )
#		quartz(file=paste(tissue, gene.target.name, n.randomizations, "pdf", sep="."))
	plot(x.rand, y.rand, type="l", lwd=2, xlab="MI", #xlim = c(max(min(x), 10^-5), max(x)), ylim = range(c(y, y.rand)), 
			col="red", 
			ylab = "P(MI)", main=paste(tissue, gene.target.name, n.randomizations, sep="  "))
	points(x, y, type="l", lwd=2, col = "black")
	legend("topright", c("actual gene target vs all gene sets", 
					"randomized gene target vs all gene sets"), col=c("black", "red"), lwd=c(2,2))
	browser()
#		dev.off()
#	}
	
	MI = MI.vector[signature.indices]
	FDR = vector(length=length(MI))
#	browser()
	ranked.MI.vector = rank(-MI.vector)  # take negative so rank 1 corresponds to highest value
	ranked.MI.matrix.rand = rank(MI.matrix.rand)
	median.MI.rand = median(MI.matrix.rand)
	if( gene.target.name[1] == "EGFR_AMP" #|| gene.target.name =="TP53"
			){ browser() }
	for( i in 1:length(MI)){
		if( MI[i] > median.MI.rand ){
			rank.observed = ranked.MI.vector[signature.indices[i]]
			rank.randomized = sum(MI[i] < MI.matrix.rand)
		} else{ 
			
			rank.observed = n.signatures - ranked.MI.vector[signature.indices[i]] + 1
			rank.randomized = sum(MI[i] > MI.matrix.rand)
#			browser()
		}
		
		FDR[i] = (rank.randomized/n.randomizations)/rank.observed
#		if( MI[i] <= median.MI.rand ){ browser() }
	}
	
	
#	MI.rand.ind = which(x.rand >= MI)
	
#	MI.integral = sum(x[MI.ind]*y[MI.ind], na.rm=T)
#	MI.rand.integral = sum(x.rand[MI.ind]*y.rand[MI.ind], na.rm=T)
#	FDR = MI.rand.integral/MI.integral
#	browser()
	
	return(list(MI=MI, FDR=FDR))
}

mutual.inf.3.v2 <- function(target.vector, comparison.matrix, n.grid=100, target.vector.name = "",
		tissue = "NA", normalize.by = "HXY", pos.and.neg=T, print.MI.ref=TRUE,
		normalize.by.MI.ref = TRUE) {
	## How this is different from mutual.inf.2:
	## x is a target pathway and y is a matrix
	## calculates a false discovery rate for the mutual information 
	## of the pathway with a randomized version of x
	## y.ind indicates the row index of the target gene set / signature
	
	# x and the elements of y can be binary or continuous
	# If there is not sufficient variation in x and y, 
	# will take the standard deviation as the bandwidth
	# (IQR finds the inter-quartile range of the data vector)
	
	# "signature.indices" is the indices of signatures that you are interested in using.
	# This code is used in comparing the chosen signatures to SUMMARY, the mutation
	# status of all the cell lines. 
	
	n.signatures = length(comparison.matrix[,1])
	MI.vector = vector(length=n.signatures, mode="double")
	
	MI.ref = mutual.inf.2(target.vector, target.vector, normalize.by=normalize.by)
	if(print.MI.ref) print(paste("MI.ref =", MI.ref))
	MI.vector = apply(comparison.matrix, 
			MARGIN=1, 
			FUN=mutual.inf.2, 
			target.vector, 
			normalize.by=normalize.by, 
			pos.and.neg=pos.and.neg)
	MI = MI.vector/ifelse(normalize.by.MI.ref, MI.ref, 1)
	FDR = rep(1, length=length(MI))
	
	return(list(MI=MI, FDR=FDR))
}



mutual.inf.4 <- function(gene.targets, signature.matrix, signature.index, n.grid=100, gene.target.name = "",
		n.randomizations = 100) {
	## How this is different from mutual.inf.2:
	## x is a target pathway and y is a matrix
	## calculates a false discovery rate for the mutual information 
	## of the pathway with a randomized version of x
	## y.ind indicates the row index of the target gene set / signature
	
	# x and the elements of y can be binary or continuous
	# If there is not sufficient variation in x and y, 
	# will take the standard deviation as the bandwidth
	# (IQR finds the inter-quartile range of the data vector)
	
	# "signature.indices" is the index of the "winning" signature that is used to compare to
	# all the genomic aberrations.
#	browser()
	if( is.vector(gene.targets)){
		gene.targets = t(as.matrix(gene.targets))
	}
	n.gene.targets = length(gene.targets[,1])
	n.signatures = length(signature.matrix[,1])
	n.samples = length(gene.targets[1,])
	MI.matrix = matrix(ncol = n.signatures, nrow = n.gene.targets)
	MI.matrix.rand = matrix(ncol = n.signatures, nrow = n.randomizations)
	gene.target.rand = t(replicate(n.randomizations, sample(gene.targets[1,])))
#	temp.MI.rand = vector(length=n.iter)
#	browser()
#	for( i in 1:length(signature.matrix[,1]) ){
#		browser()
	MI.matrix = apply(signature.matrix, MARGIN=1, 
			FUN=mutual.inf.2.multiple.gene.targets, 
			gene.targets)
	MI.matrix.rand = apply(signature.matrix, 
			MARGIN=1, FUN=mutual.inf.2.multiple.gene.targets, 
			gene.target.rand)
#		MI.matrix.rand[i,] = apply(gene.target.rand, 
#						MARGIN=1, FUN=mutual.inf.2, 
#						signature.matrix[i,]) 
#		browser()
#		MI.vector.rand[i] = mean(temp.MI.rand)
#	}
#	x.rand = sample(x)
	
#	browser()
	quartz()
	temp <- density(MI.matrix, adjust=1, n = 512, from=min(MI.matrix), to=max(MI.matrix))
	x <- temp$x
	y <- temp$y/sum(temp$y)
	
	temp.rand <- density(MI.matrix.rand, adjust=1, n = 512, from=min(MI.matrix), to=max(MI.matrix))
	x.rand <- temp.rand$x
	y.rand <- temp.rand$y/sum(temp.rand$y)
	
#	pdf(file=paste(tissue, paste(gene.target.name, collapse="-"), 
#rownames(signature.matrix)[signature.index], n.randomizations, "pdf", sep="."))
#	quartz(file=paste(tissue, paste(gene.target.name, collapse="-"), 
#rownames(signature.matrix)[signature.index], n.randomizations, "pdf", sep="."))
#if( gene.gene.target.name =="SUMMARY"){
	plot(x.rand, y.rand, type="l", lwd=2, xlab="MI", 
			#xlim = c(max(min(x), 10^-5), max(x)), ylim = range(c(y, y.rand)), 
			col="red", 
			ylab = "P(MI)", main=paste(tissue, 
					paste(gene.target.name, collapse=" "), 
					rownames(signature.matrix)[signature.index], 
					n.randomizations, sep="  "))
	points(x, y, type="l", lwd=2, col = "black")
	legend("topright", c("actual gene target(s) vs all gene sets", 
					"randomized gene target vs all gene sets"), col=c("black", "red"), lwd=c(2,2))
#	if( gene.target.name[1] =="KRAS_AMP") {browser()}
	browser()
#	dev.off()
#}
#
#	browser()
#	MI = ifelse(is.matrix(MI.matrix), MI.matrix[,signature.index], MI.matrix[signature.index])
#	ranked.MI.matrix = ifelse( is.matrix(MI.matrix), apply(-MI.matrix, MARGIN=1, rank), rank(-MI.matrix))
#	MI.vector = ifelse(is.matrix(MI.matrix))
	
	
	if(is.matrix(MI.matrix)){
		MI = MI.matrix[,signature.index]
		ranked.MI.matrix =  apply(-MI.matrix, MARGIN=1, rank)
		FDR = vector(length=n.gene.targets, mode="numeric")
#		browser()
		for( i in 1:n.gene.targets){
			if( MI[i] > median(MI.matrix.rand) ){
				rank.observed = ranked.MI.matrix[signature.index, i]
				rank.randomized = sum(MI[i] < MI.matrix.rand)
			} else{ 
#				browser()
				rank.observed = n.signatures - ranked.MI.matrix[signature.index, i]
				rank.randomized = sum(MI[i] > MI.matrix.rand)
			}
			FDR[i] = (rank.randomized/n.randomizations)/rank.observed
#		browser()
		}
	} else{ 
		MI = MI.matrix[signature.index]
		ranked.MI.matrix = rank(-MI.matrix)
		if( MI > median(MI.matrix.rand)){
			rank.observed = ranked.MI.matrix[signature.index]
			rank.randomized = sum(MI <= MI.matrix.rand)
		} else{
			rank.observed = n.signatures - ranked.MI.matrix[signature.index]
			rank.randomized = sum(MI >= MI.matrix.rand)
		}
		FDR = (rank.randomized/n.randomizations)/rank.observed
	}
#	if( gene.target.name == "EGFR_AMP" #|| gene.target.name =="TP53"
#			){ browser() }
	
#	if( MI > 0 ){
#		MI.ind = which(x >= MI) 
#	} else{ MI.ind = which( x <= MI) }
#	MI.rand.ind = which(x.rand >= MI)
	
#	MI.integral = sum(x[MI.ind]*y[MI.ind], na.rm=T)
#	MI.rand.integral = sum(x.rand[MI.ind]*y.rand[MI.ind], na.rm=T)
#	FDR = MI.rand.integral/MI.integral
#	browser()
	
	return(list(MI=MI, FDR=FDR))
}

mutual.inf.4.v2 <- function(gene.targets, signature.matrix, signature.index, n.grid=100, gene.target.name = "",
		n.randomizations = 100) {
	## How this is different from mutual.inf.2:
	## x is a target pathway and y is a matrix
	## calculates a false discovery rate for the mutual information 
	## of the pathway with a randomized version of x
	## y.ind indicates the row index of the target gene set / signature
	
	# x and the elements of y can be binary or continuous
	# If there is not sufficient variation in x and y, 
	# will take the standard deviation as the bandwidth
	# (IQR finds the inter-quartile range of the data vector)
	
	# "signature.indices" is the index of the "winning" signature that is used to compare to
	# all the genomic aberrations.
#	browser()
	if( is.vector(gene.targets)){
		gene.targets = t(as.matrix(gene.targets))
	}
	n.gene.targets = length(gene.targets[,1])
	n.signatures = length(signature.matrix[,1])
	n.samples = length(gene.targets[1,])
	MI.matrix = matrix(ncol = n.signatures, nrow = n.gene.targets)
	if( n.randomizations > 0 ){
		MI.matrix.rand = matrix(ncol = n.signatures, nrow = n.randomizations)
		gene.target.rand = t(replicate(n.randomizations, sample(gene.targets[1,])))
		MI.matrix.rand = apply(signature.matrix, 
				MARGIN=1, FUN=mutual.inf.2.multiple.gene.targets, 
				gene.target.rand)
		MI.matrix.rand = normalize(MI.matrix.rand)
	}
#	temp.MI.rand = vector(length=n.iter)
#	browser()
#	for( i in 1:length(signature.matrix[,1]) ){
#		browser()
	MI.matrix = apply(signature.matrix, MARGIN=1, 
			FUN=mutual.inf.2.multiple.gene.targets, 
			gene.targets)
	MI.matrix = normalize(MI.matrix)
#		MI.matrix.rand[i,] = apply(gene.target.rand, 
#						MARGIN=1, FUN=mutual.inf.2, 
#						signature.matrix[i,]) 
#		browser()
#		MI.vector.rand[i] = mean(temp.MI.rand)
#	}
#	x.rand = sample(x)
	
#	browser()
#	quartz()
#	temp <- density(MI.matrix, adjust=1, n = 512, from=min(MI.matrix), to=max(MI.matrix))
#	x <- temp$x
#	y <- temp$y/sum(temp$y)
#	
#	temp.rand <- density(MI.matrix.rand, adjust=1, n = 512, from=min(MI.matrix), to=max(MI.matrix))
#	x.rand <- temp.rand$x
#	y.rand <- temp.rand$y/sum(temp.rand$y)
#	
	##	pdf(file=paste(tissue, paste(gene.target.name, collapse="-"), 
#	rownames(signature.matrix)[signature.index], n.randomizations, "pdf", sep="."))
	##	quartz(file=paste(tissue, paste(gene.target.name, collapse="-"), 
#	rownames(signature.matrix)[signature.index], n.randomizations, "pdf", sep="."))
	##if( gene.gene.target.name =="SUMMARY"){
#	plot(x.rand, y.rand, type="l", lwd=2, xlab="MI", #xlim = c(max(min(x), 
#10^-5), max(x)), ylim = range(c(y, y.rand)), 
#			col="red", 
#			ylab = "P(MI)", main=paste(tissue, paste(gene.target.name, collapse=" "), 
#rownames(signature.matrix)[signature.index], n.randomizations, sep="  "))
#	points(x, y, type="l", lwd=2, col = "black")
#	legend("topright", c("actual gene target(s) vs all gene sets", 
#"randomized gene target vs all gene sets"), col=c("black", "red"), lwd=c(2,2))
	##	if( gene.target.name[1] =="KRAS_AMP") {browser()}
#	browser()
#	dev.off()
#}
#
#	browser()
#	MI = ifelse(is.matrix(MI.matrix), MI.matrix[,signature.index], MI.matrix[signature.index])
#	ranked.MI.matrix = ifelse( is.matrix(MI.matrix), apply(-MI.matrix, MARGIN=1, rank), rank(-MI.matrix))
#	MI.vector = ifelse(is.matrix(MI.matrix))
#	browser()
#MI.ref = mutual.inf.2( signature.matrix )
	
	
	
	if(is.matrix(MI.matrix)){
		MI = MI.matrix[,signature.index]
		if( n.randomizations > 0 ){
			ranked.MI.matrix =  apply(-MI.matrix, MARGIN=1, rank)
			FDR = vector(length=n.gene.targets, mode="numeric")
#		browser()
			for( i in 1:n.gene.targets){
				if( MI[i] > median(MI.matrix.rand) ){
					rank.observed = ranked.MI.matrix[signature.index, i]
					rank.randomized = sum(MI[i] < MI.matrix.rand)
				} else{ 
#				browser()
					rank.observed = n.signatures - ranked.MI.matrix[signature.index, i]
					rank.randomized = sum(MI[i] > MI.matrix.rand)
				}
				FDR[i] = (rank.randomized/n.randomizations)/rank.observed
#		browser()
			}
		} else{ FDR = rep(1, length=n.gene.targets) }
	} else{ 
		MI = MI.matrix[signature.index]
		if( n.randomizations > 0 ){
			ranked.MI.matrix = rank(-MI.matrix)
			if( MI > median(MI.matrix.rand)){
				rank.observed = ranked.MI.matrix[signature.index]
				rank.randomized = sum(MI <= MI.matrix.rand)
			} else{
				rank.observed = n.signatures - ranked.MI.matrix[signature.index]
				rank.randomized = sum(MI >= MI.matrix.rand)
			}
			FDR = (rank.randomized/n.randomizations)/rank.observed
		} else{ FDR = rep(1, length=n.gene.targets) }
	}
#	if( gene.target.name == "EGFR_AMP" #|| gene.target.name =="TP53"
#			){ browser() }
	
#	if( MI > 0 ){
#		MI.ind = which(x >= MI) 
#	} else{ MI.ind = which( x <= MI) }
#	MI.rand.ind = which(x.rand >= MI)
	
#	MI.integral = sum(x[MI.ind]*y[MI.ind], na.rm=T)
#	MI.rand.integral = sum(x.rand[MI.ind]*y.rand[MI.ind], na.rm=T)
#	FDR = MI.rand.integral/MI.integral
#	browser()
	
	return(list(MI=MI, FDR=FDR))
}



mise <- function( x ) {
	n.x = length(x)
	r = seq(2,10)
	f = vector(length=n.x, mode=mode(x))
	for( i in 1:n.x ){
		f = f + (x - x[i])^2
	}
	expected.value = sum(f)/n.x
}

amise <- function( v ){
	# Reference: 
	# "Very fast optimal bandwith selection for univariate kernel density estimation"
	# Vikas Chandrakant Raykar and Ramani Duraiswami
	# [CS-TR-4774/UMIACS-TR-2005-73] June 28, 2006
	
	H4 <- function( x ){ x^4 - 6*x^2 + 3 }
	H6 <- function( x ){ x^6 - 15*x^4 + 45*x^2 - 15 } 
	
	N = length(v)
	
	# Step 1 on page 11 of reference
	sigma = mean(v)
	
	# Step 2 on p. 11
	Phi6 = sigma^(-7)*(-15/(16*sqrt(pi)))
	Phi8 = sigma^(-9)*(-105/(32*sqrt(pi)))
	
	# Step 3 on p. 11
	g1 = ( -6/(sqrt(2*pi) * Phi6 * N))^(1/7)
	g2 = ( 30/(sqrt(2*pi) * Phi8 * N))^(1/9)
	
	# Make a matrix Z where Z(i,j) = x_i - x_j
	Z = matrix(v, ncol = N, nrow = N, byrow=TRUE) - matrix(v, ncol = N, nrow = N, byrow=FALSE)
	
	Phi4 <- function( g ) 1/(N*(N-1)*sqrt(2*pi)*g^5) * sum( H4(Z/g) * exp( -(Z^2)/(2*g^2)) )
	Phi4.g1 = Phi4(g1)
	Phi6.g2 = 1/(N*(N-1)*sqrt(2*pi)*g1^5) * sum( H4(Z/g1) * exp( -(Z^2)/(2*g1^2)) )
	
	# Step 4
	Y <- function( h ){
		( (-6*sqrt(2)*Phi4.g1)/Phi6.g2)^(1/7)*h^(5/7)
	}
	
	fxn <- function( h ){
		h - ( 1/ (sqrt(2) * Phi4(Y(h)) * N))^(1/5)
	}
	newtonraphson(fxn, mean(v))
}

#parzen.window <- function(z, h){
#	Sigma = cov(z,z)
#	
#	
#}

write.cls.with.locs <- function( output.cls,
		cls.labels,
		phen.names){
	
	class.order.list = apply(cls.labels, 1, unique)
	num = sum(unlist(lapply(class.order.list, length)))
	class.order.list = unlist(class.order.list)
	
	class.phen <- unique(cls.labels)
	n <- length(class.phen)
	l <- length(cls.list[1,])
	
	phen.names.string <- paste("phen.names:", paste(phen.names, collapse=" "), "col.names:", sep=" ")
	cat(paste(l, num, length(cls.list[, 1]), phen.names.string, sep=" "), "\n", 
			file = output.cls, append = FALSE, sep = "")
	cat("# ", paste(class.order.list, collapse=" "), "\n", file = output.cls, append = TRUE, sep = "")
	for (i in 1:length(cls.list[,1])) {
		cat(paste(cls.list[i,], collapse=" "), "\n", file = output.cls, append = TRUE, sep = "")
	}
	
}

mutual.inf.P <- function(x, y, n.grid=100, use.sign=TRUE) {
	# for definitions of mutual information and the universal metric (NMI) see the 
	# definition of "Mutual Information" in wikipedia and Thomas and Cover's book
	
#   kde2d.xy <- kde2d(x, y, n = n.grid, h = c(width.SJ(x, method="dpi"), width.SJ(y, method="dpi")))
	kde2d.xy <- kde2d(x, y, n = n.grid, h = c(bcv(x), bcv(y)))
	X <- kde2d.xy$x
	Y <- kde2d.xy$y  
	PXY <- kde2d.xy$z/sum(kde2d.xy$z)
	
#   PX <- apply(PXY, MARGIN=1, sum)
	PX <- rowSums(PXY)
	PX <- PX/sum(PX)
	HX <- -sum(PX * log2(PX))
	PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
	
#   PY <- apply(PXY, MARGIN=2, sum)
	PY <- colSums(PXY)
	PY <- PY/sum(PY)
	HY <- -sum(PY * log2(PY))
	PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)
	
	MI <- sum(PXY * log2(PXY/(PX*PY)))
	MI
	HXY <- - sum(PXY * log2(PXY))
	
#	NMI <- sign(cor(x, y)) * ((HX + HY)/HXY - 1)  # use peason correlation the get the sign (directionality)
	NMI <- ifelse(use.sign, sign(cor(x, y)), 1) * ((HX + HY)/HXY - 1)  # no directionality
	
	return(list(MI=MI, HXY=HXY, HX=HX, HY=HY, NMI=NMI))
}


REVEALER.normalize.zero.one <- function( v ){
	if(length(unique(v))==1){ 
		#browser()
		normalized = rep(0, length(v))
	} else {
		#browser()
		normalized = (v - min(v))/(max(v) - min(v))
		
	}
	return(normalized)
}

OPAM.Evaluate.Results <- function(
		input.ds,
		input.cls,
		phenotype = NULL,
		target.class = NULL,
		target.type = "discrete",
		output.txt,
		output.pdf) {
	
	pdf(file=output.pdf, height=8.5, width=11)
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	model.names <- dataset$row.names
	model.descs <- dataset$descs
	Ns <- length(m[1,])
	for (i in 1:length(m[,1])) {
		if (sd(m[i,]) == 0) {
			val <- m[i, 1]
			m[i,] <- m[i,] + runif(n=Ns, min= val - 0.001, max=val + 0.001)  # add small noise to flat profiles
		}
	}
	
	dim(m)
	sample.names <- dataset$names
	CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
	cls.labels <- CLS$class.v
	cls.phen <- CLS$phen
	cls.list <- CLS$class.list 
	
	library(verification)
	
	if (is.null(phenotype)) {
		phen.loc <- 1
	} else {
		phen.loc <- match(phenotype, CLS$phen.names)
	}
	if (is.vector(CLS$class.list)) {
		target.vec <- CLS$class.list
	} else {
		target.vec <- CLS$class.list[phen.loc,]
	}
	if (target.type == "continuous") {
		target <- target.vec
	} else if (target.type == "discrete") {
		target <- ifelse(target.vec == target.class, 1, 0)    
	}
	
	ind <- order(target)
	target <- target[ind]
	target.vec <- target.vec[ind]
	m <- m[, ind]
	sample.names[ind]
	class.v <- CLS$class.v
	if (is.vector(class.v)) {
		class.v <- class.v[ind]
	} else {
		class.v <- class.v[, ind]
	}
	annot <- MI <- AUC <- AUC.pval <- t.stat <- t.pval <- vector(length=dim(m)[1], mode="numeric")
	
	NMI.ref <- mutual.inf.P(x = target, y = target, n.grid=100)$NMI
	
	for (i in 1:dim(m)[1]) {
		feature <- m[i,]
		MI[i] <- signif(mutual.inf.P(target, feature, n.grid=100)$NMI/NMI.ref, 4)
		if (target.type == "continuous") {
			AUC[i] <- AUC.pval[i] <- t.stat[i] <- t.pval[i] <- "-"
		} else if (target.type == "discrete") {
			feature.norm <- (feature - min(feature))/(max(feature) - min(feature))
			perf.auc <- roc.area(target, feature.norm)
			AUC[i] <- ifelse(perf.auc$A < 0.5, -(1 - perf.auc$A), perf.auc$A)
			AUC[i] <- signif(AUC[i], digits=4)
			p.val <- perf.auc$p.value
			p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
			AUC.pval[i] <- signif(p.val, digits=4)
			temp <- split(feature, target)
			x <- temp$'1'
			y <- temp$'0'
			t.stat[i] <- signif(t.test(x=x, y=y)$statistic, digits=4)
			p.val <- t.test(x=x, y=y)$p.value
			p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
			t.pval[i] <- signif(p.val, digits=4)
		}
		annot[i] <- paste(MI[i], "     ", AUC[i], " (", AUC.pval[i], ")    ", t.stat[i], " (", t.pval[i], ")", sep="")
	}
	
	mycol <- vector(length=512, mode = "numeric")
	for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
	for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
	mycol <- rev(mycol)
	cex.axis = 1
	ncolors <- length(mycol)
	
	nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(4, 10), FALSE)
	par(mar = c(1, 15, 5, 15))
	max.v <- max(max(target), -min(target))
	V1 <- target
	image(1:length(target), 1:1, as.matrix(V1), zlim = c(0, 1), 
			col=c("yellow", "purple"), axes=FALSE, main="", sub = "", xlab= "", ylab="")
	axis(2, at=1:1, labels=paste(phenotype, target.class), 
			adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
	axis(4, at=1:1, labels=paste("NMI     AUC (p-val)     t-test (p-val)", sep=""), adj= 0.5, tick=FALSE, 
			las = 1, cex.axis=0.80, font.axis=1, line=-1) 
	par(mar = c(5, 15, 1, 15))
	V1 <- m
	for (i in 1:dim(V1)[1]) V1[i,] <- (V1[i,] - mean(V1[i,]))/sd(V1[i,])
	max.v <- max(max(V1), -min(V1))
	V1 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
	V1 <- apply(V1, MARGIN=2, FUN=rev)
	image(1:dim(V1)[2], 1:dim(V1)[1], t(V1), zlim = c(0, ncolors), 
			col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
	axis(2, at=1:dim(V1)[1], labels=row.names(V1), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
	axis(4, at=1:dim(V1)[1], labels=rev(annot), adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
	axis(1, at=1:dim(V1)[2], labels=sample.names, adj= 0.5, tick=FALSE, las = 3, cex.axis=0.60, font.axis=1, line=-1)
	
	dev.off()
	
	annot2 <- data.frame(cbind(MI, AUC, AUC.pval, t.stat, t.pval))
	row.names(annot2) <- row.names(m)
	write(paste(c("gene set ", noquote(colnames(annot2))), collapse="\t"), 
			file = output.txt, append = F, ncolumns = length(colnames(annot2)))
	write.table(annot2, file=output.txt, append=T, quote=F, sep="\t", eol="\n", col.names=F, row.names=T)
	
}

OPAM.Evaluate.Results.2 <- function(
		input.ds,
		input.cls,
		phenotype = NULL,
		target.class = NULL,
		target.type = "discrete",
		sort.results = T,
		display.top.n = 20,
		output.txt,
		output.graphics.filename,
		signatures.of.interest = NULL,
		is.full = FALSE,
		statistic = NULL,
		weight = NULL,
		output.MI.Ranks.txt = NULL,
		graphics.file.type = "pdf", # pdf, eps, jpg / jpeg, bmp, png, tiff
		gene.of.interest = NULL,
		tissue = "NA",
		highlight.tissue.name = NULL,  # highlight a single tissue in a dataset with multiple tissues
		ifMultipleTissues = FALSE,
		ifTissueGenePlotTitle = FALSE,
		ifUseAllClasses = FALSE,
		ifSortWithinClasses = FALSE,
		ifScoreBarPlot = TRUE
) {
	
	isCLS.input.cls = (length(grep(".cls$", input.cls)) > 0)
	isGCT.input.cls = (length(grep(".gct$", input.cls)) > 0)
	
	graphics.file.suffix = ""
	if(!is.null(highlight.tissue.name)){
		graphics.file.suffix = paste(graphics.file.suffix, "_", 
				paste(highlight.tissue.name, collapse="."),
				".highlighted", sep="")
	}
	if( ifUseAllClasses ){
		graphics.file.suffix = paste(graphics.file.suffix, ".allClasses", sep="")
	}
	if( ifSortWithinClasses){
		graphics.file.suffix = paste(graphics.file.suffix, ".sortedWithinClasses", sep="")
	}
	
	
#	graphics.file.types = c("pdf", "eps", "jpeg", "bmp", "png")
#	browser()
	if(graphics.file.type == "pdf"){
		if( length(grep(".pdf$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".pdf", sep="")
		} else{ 
			output.graphics.filename = gsub(".pdf", 
					paste(graphics.file.suffix, ".pdf", sep=""), output.graphics.filename)}
		pdf(file=output.graphics.filename, height=8.5, width=11)
	} else if(graphics.file.type == "eps"){
		if( length(grep(".eps$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".eps", sep="")
		} else{ 
			output.graphics.filename= gsub(".eps", 
					paste(graphics.file.suffix, ".eps", sep=""), output.graphics.filename)}
		eps(file=paste(output.graphics.filename, ".eps", sep=""), height=8.5, width=11)
	} else if(graphics.file.type == "jpeg"){
		if( length(grep(".jpeg$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".jpeg", sep="")
		} else{ 
			output.graphics.filename = gsub(".jpeg", 
					paste(graphics.file.suffix, ".jpeg", sep=""), output.graphics.filename)}
		jpeg(file=paste(output.graphics.filename, ".jpg", sep=""), height=8.5, width=11)
	} else if(graphics.file.type == "bmp"){
		if( length(grep(".bmp$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".bmp", sep="")
		} else{ 
			output.graphics.filename = gsub(".bmp", 
					paste(graphics.file.suffix, ".bmp", sep=""), output.graphics.filename)}
		bmp(file=paste(output.graphics.filename, ".bmp", sep=""), height=8.5, width=11)
	} else if(graphics.file.type == "png"){
		if( length(grep(".png$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".png", sep="")
		} else{ 
			output.graphics.filename = gsub(".png", 
					paste(graphics.file.suffix, ".png", sep=""), output.graphics.filename)}
		png(file=paste(output.graphics.filename, ".png", sep=""), height=8.5, width=11)
	} else if(graphics.file.type == "tiff"){
		if( length(grep(".tiff$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".tiff", sep="")
		} else{ 
			output.graphics.filename = gsub(".tiff", 
					paste(graphics.file.suffix, ".tiff", sep=""), output.graphics.filename)}
		tiff(file=paste(output.graphics.filename, ".eps", sep=""), height=8.5, width=11)
	} else{
		stop("Please specify heatmap output file type")
	}
	
	
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	nGeneSets <- dim(m)[1]
	model.names <- dataset$row.names
	model.descs <- dataset$descs
	nSample <- length(as.matrix(m)[1,])
	for (i in 1:nGeneSets) {
		if ( !sum(is.na(m[i,])) && sd(as.matrix(m)[i,]) == 0) {
			val <- as.matrix(m)[i, 1]
			m[i,] <- as.matrix(m)[i,] + runif(n=nSample, 
					min= val - 0.001, max=val + 0.001)  # add small noise to flat profiles
		}
	}
	
	dim(m)
	sample.names <- dataset$names
	
	if( ifMultipleTissues ){
		tissue.type <- vector(length=nSample, mode="character")
#		temp = strsplit(sample.names, split="_")
		for (k in 1:nSample) {
			temp <- strsplit(sample.names[k], split="_") 
			tissue.type[k] <- paste(temp[[1]][2:length(temp[[1]])], collapse="_")
		}
		tissue.names = unique(tissue.type)
		tissue.labels = match(tissue.type, tissue.names)
		n.tissues = length(tissue.names)
		
		library(RColorBrewer)
		tissue.colors = c(brewer.pal(12, "Set3"), brewer.pal(12,"Paired"))[
				1:n.tissues]
		if( !is.null(highlight.tissue.name)){
			highlight.tissue.index = which(tissue.names %in% highlight.tissue.name)
#			tissue.labels = ifelse(tissue.labels == highlight.tissue.index, tissue.labels, 1)
			tissue.colors[-highlight.tissue.index] = c(brewer.pal(12, "Set3"), 
					brewer.pal(12,"Paired"))[2]  # 2nd color is yellow and is nice to have as the default background
			# as the other colors are darker
			tissue.colors[highlight.tissue.index] = c(brewer.pal(12, "Set3"), 
					brewer.pal(12,"Paired"))[1:(length(highlight.tissue.name)+1)][-2]
		}
		print(matrix(c(tissue.names, tissue.colors), ncol=2), quote=FALSE)
	}
#	else{
#		tissue.names = tissue
#		tissue.labels = rep(1, nSample)
#	}
	
	if( isCLS.input.cls ){
		CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
		cls.labels <- CLS$class.v
		cls.phen <- CLS$phen
		cls.list <- CLS$class.list 
		
#		browser()
		
		library(verification)
		
		if (is.null(phenotype)) {
			phen.loc <- 1
		} else {
			phen.loc <- match(phenotype, CLS$phen.names)
		}
		if (is.vector(CLS$class.list)) {
			target.vec <- CLS$class.list
		} else {
			target.vec <- CLS$class.list[phen.loc,]
		}
		if( ifUseAllClasses ){
			classes = unique(target.vec)
			target.matrix = t(sapply(classes, FUN=function(x) ifelse(target.vec ==x, 1, 0)))
			target.type = "discrete"
			n.classes = length(classes)
			display.top.n.per.class = ceiling(display.top.n/n.classes)
			display.top.n = display.top.n.per.class*n.classes
			if( length(CLS$col.phen) == 
					n.classes){  # If the CLS file provides enough colors to describe all the classes
				class.colors = CLS$col.phen
			} else {  # Otherwise, choose colors for the classes
				class.colors = c(#brewer.pal(12, "Set3"), brewer.pal(12,"Paired"), 
						brewer.pal(12, "Set3"), brewer.pal(12,"Paired"))[
						1:n.classes]
			}
		} else if (!ifUseAllClasses && target.type == "continuous") {
			target <- target.vec
		} else if (!ifUseAllClasses && target.type == "discrete") {
			target <- ifelse(target.vec == target.class, 1, 0) 
			class.colors = c("yellow", "purple")
			other.classes = paste(CLS$class.list[which(CLS$class.list != target.class)], collapse ="  ")
			classes = c( target.class, other.classes )
		}
	} else if( isGCT.input.cls && !is.null(gene.of.interest) ){
		target.dataset <- MSIG.Gct2Frame(filename = input.cls)
		target = data.matrix(target.dataset$ds)[grep(paste("^", 
								gene.of.interest, sep=""), target.dataset$row.names),]
		target = target.vec = REVEALER.normalize.zero.one(target)
		phenotype = paste(gene.of.interest, "expression")
		target.class = ""
		target.type = "continuous"
	}
	
	if(!ifUseAllClasses){
		ind <- order(target)
		target <- target[ind]
		target.vec <- target.vec[ind]
		m <- as.matrix(m)[, ind]
		sample.names[ind]
		class.v <- CLS$class.v
		if( isCLS.input.cls ){
			if (is.vector(class.v)) {
				class.v <- class.v[ind]
			} else {
				class.v <- class.v[, ind]
			}
		}
		target.matrix = matrix(target, nrow=1)
	} else{ 
		annot.multipleClasses <- MI.multipleClasses <- AUC.multipleClasses <- 
				AUC.pval.multipleClasses <- t.stat.multipleClasses <- 
				t.pval.multipleClasses <- vector(length=display.top.n, 
						mode="numeric")
		m.multipleClasses = matrix(nrow=(display.top.n), ncol=nSample)
		rownames(m.multipleClasses) = rep("", display.top.n)
	}
	
	for( target.row in 1:length(target.matrix[,1])){
		target = target.matrix[target.row,]
		annot <- MI <- AUC <- AUC.pval <- t.stat <- t.pval <- vector(length=nGeneSets, mode="numeric")
		
		NMI.ref <- mutual.inf.2(x = target, y = target, n.grid=100)#$NMI
#	browser()
		for (i in 1:nGeneSets) {
			if (nGeneSets == 1) {
				feature <- m
			} else {
				feature <- m[i,]
			}
			MI[i] <- signif(mutual.inf.P(target, feature, n.grid=100)$NMI/NMI.ref, 4)
			if (target.type == "continuous") {
				AUC[i] <- AUC.pval[i] <- t.stat[i] <- t.pval[i] <- "-"
			} else if (target.type == "discrete") {
				feature.norm <- (feature - min(feature))/(max(feature) - min(feature))
				perf.auc <- roc.area(target, feature.norm)
				AUC[i] <- ifelse(perf.auc$A < 0.5, -(1 - perf.auc$A), perf.auc$A)
				AUC[i] <- signif(AUC[i], digits=4)
				p.val <- perf.auc$p.value
				p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
				AUC.pval[i] <- signif(p.val, digits=4)
				temp <- split(feature, target)
				x <- temp$'1'
				y <- temp$'0'
				t.stat[i] <- signif(t.test(x=x, y=y)$statistic, digits=4)
				p.val <- t.test(x=x, y=y)$p.value
				p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
				t.pval[i] <- signif(p.val, digits=4)
			}
			annot[i] <- paste(MI[i], "     ", AUC[i], " (", AUC.pval[i], ")    ", t.stat[i], " (", t.pval[i], ")", sep="")
		}
		
		if ((nGeneSets > 1) & (sort.results == T)) {
			MI.order <- order(MI, decreasing=T)
			MI <- MI[MI.order]
			AUC <- AUC[MI.order]
			AUC.pval <- AUC.pval[MI.order]
			t.stat <- t.stat[MI.order]
			t.pval <= t.pval[MI.order]
			m <- as.matrix(m)[MI.order,]
			annot <- annot[MI.order]
		}
		
		if(ifUseAllClasses){
			
			from.ind = (display.top.n.per.class*(target.row-1)+1)
			to.ind = (display.top.n.per.class*(target.row))
			annot.multipleClasses[from.ind:to.ind] = annot[1:display.top.n.per.class]
			m.multipleClasses[from.ind:to.ind,] = m[1:display.top.n.per.class,]
			row.names(m.multipleClasses)[from.ind:to.ind] = row.names(m)[1:display.top.n.per.class]
			
			MI.multipleClasses[from.ind:to.ind] = MI[1:display.top.n.per.class]
			AUC.multipleClasses[from.ind:to.ind] = AUC[1:display.top.n.per.class]
			AUC.pval.multipleClasses[from.ind:to.ind] = AUC.pval[1:display.top.n.per.class]
			t.stat.multipleClasses[from.ind:to.ind] = t.stat[1:display.top.n.per.class]
			t.pval.multipleClasses[from.ind:to.ind] = t.pval[1:display.top.n.per.class]
			
			if( length(grep(".txt$", output.txt))!=0 ){
				this.class.output.txt = gsub(".txt", paste(".", classes[target.row], ".txt", sep=""), output.txt)	
			} else{
				this.class.output.txt = paste(".", output.txt, classes[target.row], ".txt", sep="")
			}
		} else{ this.class.output.txt = output.txt }
		annot2 <- data.frame(cbind(MI, AUC, AUC.pval, t.stat, t.pval))
		if( nGeneSets > 1){
			row.names(annot2) <- make.unique(row.names(m))
		} else row.names(annot2) = model.names
		write(paste(c("gene set ", noquote(colnames(annot2))), collapse="\t"), file = this.class.output.txt, append = F, 
				ncol = length(colnames(annot2)))
		write.table(annot2, file=this.class.output.txt, append=T, quote=F, sep="\t", eol="\n", col.names=F, row.names=T)
	}
	
	if(ifUseAllClasses){
		annot = annot.multipleClasses
		m = m.multipleClasses
		if(ifSortWithinClasses){
#				browser()
			for( target.row in 1:length(target.matrix[,1])){
				target = target.matrix[target.row,]
				from.ind = (display.top.n.per.class*(target.row-1)+1)
				to.ind = (display.top.n.per.class*(target.row))
				
				class.sample.order = order(m[from.ind, which(target==1)])
				m[,which(target==1)] = m[,which(target==1)[class.sample.order]]
				sample.names[which(target==1)] = sample.names[which(target==1)][class.sample.order]
			}
		}
		
		
		MI = MI.multipleClasses
		AUC = AUC.multipleClasses
		AUC.pval = AUC.pval.multipleClasses
		t.stat = t.stat.multipleClasses
		t.pval = t.pval.multipleClasses
	} else if( !ifUseAllClasses && ifSortWithinClasses ){
		class.indicators = unique(target)
		for( class.indicator in class.indicators ){
			class.sample.order = order(m[1,which(target==class.indicator)])
			m[,which(target==class.indicator)] = m[,which(target==class.indicator)][class.sample.order]
			sample.names[which(target==class.indicator)] = sample.names[
					which(target==class.indicator)][class.sample.order]
		}
	}
	## Save ranks of the signatures of interest, if provided
	if( length(signatures.of.interest) > 0 ){
#		browser()
#		split.up = strsplit(signatures.of.interest, "_UP")
#		split.dn = strsplit(signatures.of.interest, "_DN")
#		browser()
		combined.signature = substr(signatures.of.interest[1], 1, nchar(signatures.of.interest[1])-3) 
		#split.up[[which(split.up %in% split.dn)]]
		signatures.of.interest = c(signatures.of.interest[1], combined.signature, signatures.of.interest[2])
		signature.ranks = which(rownames(m) %in% signatures.of.interest)
#		browser()
		signature.ranks.text = c(paste(na.omit(signature.ranks), rownames(m)[na.omit(signature.ranks)], 
						MI[na.omit(signature.ranks)]), 
				paste("(out of", length(m[,1]), "signatures)"))
		write(paste(weight, signif(MI[na.omit(signature.ranks)], digits=5), na.omit(signature.ranks), 
						rownames(m)[na.omit(signature.ranks)], sep="\t"), 
				file=output.MI.Ranks.txt, append=TRUE)
#		browser()
	} else{ signature.ranks.text = NULL}
	
	mycol <- vector(length=512, mode = "numeric")
	for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
	for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
	mycol <- rev(mycol)
	cex.axis = 1
	n.pinkogram.colors <- length(mycol)
	
#	browser()
	if(ifMultipleTissues || ifSortWithinClasses ){
		nf <- layout(matrix(c(1, 2, 3), ncol=1, byrow=FALSE), 1, heights = c(2, 10, 2),
				respect = FALSE)
	} else{	nf <- layout(matrix(c(1, 2), ncol=1), 1, c(4, 10), FALSE) }
	par(mar = c(1, 15, 5, 15))
	max.v <- max(max(target), -min(target))
	if( isCLS.input.cls && !ifUseAllClasses ){
		V1 <- target
	} else if( isGCT.input.cls && !ifUseAllClasses ){
		target.mid.point = which.min(abs(target - quantile(target, 0.5)))
		V1 = ceiling(n.pinkogram.colors*target)
#		V1 = ceiling(c( (n.pinkogram.colors/2)*target[1:target.mid.point], 
#				(n.pinkogram.colors/2)*target[(target.mid.point+1):nSample] + n.pinkogram.colors/2))
	}
	if(ifUseAllClasses){
		V1 = colSums(target.matrix*1:n.classes)
	} 
	n.total.colors = ifelse( target.type=="discrete" && !ifMultipleTissues 
					&& !ifUseAllClasses, 1, length(class.colors))
	tissue.colors.plus.mycol = class.colors
	if( ifMultipleTissues && target.type=="continuous" ){
#		browser()
		V1 = t(rbind(V1, (tissue.labels + n.pinkogram.colors)))
		n.total.colors = n.pinkogram.colors + n.tissues
		tissue.colors.plus.mycol = c(mycol, tissue.colors)
		
	} else if( target.type == "continuous"){ n.total.colors = n.pinkogram.colors; tissue.colors.plus.mycol = mycol }
	if (ifTissueGenePlotTitle){
		plot.title = paste("tissue:", tissue, "  gene:", gene.of.interest)
	} else{
		plot.title = strsplit(output.pdf, "/")[[1]]
		if(is.full){
			plot.title = strsplit(plot.title[length(plot.title)], ".HEATMAP.FULL.pdf")[[1]]
		} else {
			plot.title = strsplit(plot.title[length(plot.title)], ".HEATMAP.pdf")[[1]]
		}
	}
	image(1:nSample, 1:(ifelse(ifMultipleTissues, 2, 1)), 
			as.matrix(V1), zlim = c(0, n.total.colors), 
			col=tissue.colors.plus.mycol, axes=FALSE, 
			main=plot.title, sub = "", xlab= "", ylab="")
	axis(2, at=1:1, labels=paste(phenotype, target.class), 
			adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
	axis(4, at=1:1, labels=paste("NMI     AUC (p-val)     t-test (p-val)", sep=""), adj= 0.5, tick=FALSE, 
			las = 1, cex.axis=0.80, font.axis=1, line=-1) 
	par(mar = c(5, 15, 1, 15))
	
	if (display.top.n > nGeneSets) display.top.n <- nGeneSets
	
	if (nGeneSets == 1) {
		V1 <- m
		V1 <- (V1 - mean(V1))/sd(V1)
	} else {
		V1 <- m[1:display.top.n, ]
		for (i in 1:display.top.n) V1[i,] <- (V1[i,] - mean(V1[i,]))/sd(V1[i,])
	}
	
	max.v <- max(max(V1), -min(V1))
	V1 <- ceiling(n.pinkogram.colors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
	
	if (nGeneSets > 1) {
		V1 <- apply(V1, MARGIN=2, FUN=rev)
		image(1:nSample, 1:display.top.n, t(V1), zlim = c(0, n.pinkogram.colors), 
				col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
		axis(2, at=1:display.top.n, labels=row.names(V1), adj= 0.5,
				tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
		axis(4, at=1:display.top.n, labels=rev(annot[1:display.top.n]), adj= 0.5, tick=FALSE, 
				las = 1, cex.axis=0.70, font.axis=1, line=-1)
		if( nSample < 100){
			axis(1, at=1:nSample, labels=sample.names, adj= 0.5, tick=FALSE, las = 3, 
					cex.axis=0.60, font.axis=1, line=-1)
		}
	} else {
		image(1:nSample, 1:1, as.matrix(V1), zlim = c(0, n.pinkogram.colors), 
				col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
		axis(2, at=1:1, labels=model.names, adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
		axis(4, at=1:1, labels=annot[1], adj= 0.5, tick=FALSE, las = 1, cex.axis=0.70, font.axis=1, line=-1)
		if( nSample < 100){
			axis(1, at=1:nSample, labels=sample.names, adj= 0.5, tick=FALSE, las = 3, 
					cex.axis=0.60, font.axis=1, line=-1)
		}
	}
	
	## Tissue Legend
	if(ifMultipleTissues){
		#browser()
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		legend(x=0, y= 1, #horiz = T, x.intersp = 0.5, y.intersp=.25, 
				legend=tissue.names, bty="n", xjust=0, yjust= 1, 
				fill = tissue.colors, cex = 0.8, #pt.cex=1.75, 
				ncol=4)
	}
	if( ifSortWithinClasses ){
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		legend(x=0, y= 1, #horiz = T, x.intersp = 0.5, y.intersp=.25, 
				legend=classes, bty="n", xjust=0, yjust= 1, 
				fill = class.colors, cex = 1.2, #pt.cex=1.75, 
				ncol=4)
	}
	
	if(ifScoreBarPlot){
		par(mar = c(6, 5, 6, 1))
		barplot(rev(MI[display.top.n]), xlab="MI scores\nbarplot", ylab="", xlim=c(-1,1), 
				axes=TRUE, horiz=TRUE, axisnames=FALSE, 
#				main="MI scores\nbarplot", font.main=1
		)
#		browser()
#		text("MI scores\nbarplot")
	}
	dev.off()
	
	
	
}

OPAM.Evaluate.Results.3 <- function(  ## Added score bar plots to right from OPAM.Evaluate.Results.2
		input.ds,
		input.cls,
		phenotype = NULL,
		target.class = NULL,
		target.type = "discrete",
		sort.results = T,
		display.top.n = 20,
		output.txt,
		output.graphics.filename,
		signatures.of.interest = NULL,
		is.full = FALSE,
		statistic = NULL,
		weight = NULL,
		output.MI.Ranks.txt = NULL,
		graphics.file.type = "pdf", # pdf, eps, jpg / jpeg, bmp, png, tiff, window,
		graphic.height = 11,
		graphic.width = 17,
		gene.of.interest = NULL,
		tissue = "NA",
		highlight.tissue.name = NULL,  # highlight a single tissue in a dataset with multiple tissues
		ifMultipleTissues = FALSE,
		ifTissueGenePlotTitle = FALSE,
		ifUseAllClasses = FALSE,
		ifSortWithinClasses = FALSE,
		ifScoreBarPlot = TRUE
#		original.expression.ds = NULL
) {
	
	isCLS.input.cls = (length(grep(".cls$", input.cls)) > 0)
	isGCT.input.cls = (length(grep(".gct$", input.cls)) > 0)
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	nGeneSets <- dim(m)[1]
	model.names <- dataset$row.names
	model.descs <- dataset$descs
	nSample <- dim(m)[2]
	for (i in 1:nGeneSets) {
		if ( !sum(is.na(m[i,])) && sd(as.matrix(m)[i,]) == 0) {
			val <- as.matrix(m)[i, 1]
			m[i,] <- as.matrix(m)[i,] + runif(n=nSample, 
					min= val - 0.001, max=val + 0.001)  # add small noise to flat profiles
		}
	}
#	browser()
	
	dim(m)
	sample.names <- dataset$names
	
	if( ifMultipleTissues ){
		tissue.type <- vector(length=nSample, mode="character")
#		temp = strsplit(sample.names, split="_")
		for (k in 1:nSample) {
			temp <- strsplit(sample.names[k], split="_") 
			tissue.type[k] <- paste(temp[[1]][2:length(temp[[1]])], collapse="_")
		}
		tissue.names = unique(tissue.type)
		tissue.labels = match(tissue.type, tissue.names)
		n.tissues = length(tissue.names)
		
		library(RColorBrewer)
		tissue.colors = c(brewer.pal(12, "Set3"), brewer.pal(12,"Paired"))[
				1:n.tissues]
		if( !is.null(highlight.tissue.name)){
			highlight.tissue.index = which(tissue.names %in% highlight.tissue.name)
#			tissue.labels = ifelse(tissue.labels == highlight.tissue.index, tissue.labels, 1)
			tissue.colors[-highlight.tissue.index] = c(brewer.pal(12, "Set3"), 
					brewer.pal(12,"Paired"))[2]  # 2nd color is yellow and is nice to have as the default background
			# as the other colors are darker
			tissue.colors[highlight.tissue.index] = c(brewer.pal(12, "Set3"), 
					brewer.pal(12,"Paired"))[1:(length(highlight.tissue.name)+1)][-2]
		}
		print(matrix(c(tissue.names, tissue.colors), ncol=2), quote=FALSE)
	}
#	else{
#		tissue.names = tissue
#		tissue.labels = rep(1, nSample)
#	}
	
	if( isCLS.input.cls ){
		CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
		cls.labels <- CLS$class.v
		cls.phen <- CLS$phen
		cls.list <- CLS$class.list 
		
#		browser()
		
		library(verification)
		
		if (is.null(phenotype)) {
			phen.loc <- 1
		} else {
			phen.loc <- match(phenotype, CLS$phen.names)
		}
		if (is.vector(CLS$class.list)) {
			target.vec <- CLS$class.list
		} else {
			target.vec <- CLS$class.list[phen.loc,]
		}
		if( ifUseAllClasses ){
			classes = unique(target.vec)
			target.matrix = t(sapply(classes, FUN=function(x) ifelse(target.vec ==x, 1, 0)))
			target.type = "discrete"
			n.classes = length(classes)
			display.top.n.per.class = ceiling(display.top.n/n.classes)
			display.top.n = display.top.n.per.class*n.classes
			if( length(CLS$col.phen) == 
					n.classes){  # If the CLS file provides enough colors to describe all the classes
				class.colors = CLS$col.phen
			} else {  # Otherwise, choose colors for the classes
				class.colors = c(#brewer.pal(12, "Set3"), brewer.pal(12,"Paired"), 
						brewer.pal(12, "Set3"), brewer.pal(12,"Paired"))[
						1:n.classes]
			}
		} else if (!ifUseAllClasses && target.type == "continuous") {
			target <- target.vec
		} else if (!ifUseAllClasses && target.type == "discrete") {
			target <- ifelse(target.vec == target.class, 1, 0) 
			class.colors = c("yellow", "purple")
			other.classes = paste(CLS$class.list[which(CLS$class.list != target.class)], collapse ="  ")
			classes = c( target.class, other.classes )
		}
	} else if( isGCT.input.cls && !is.null(gene.of.interest) ){
		target.dataset <- MSIG.Gct2Frame(filename = input.cls)
		target = data.matrix(target.dataset$ds)[grep(paste("^", 
								gene.of.interest, sep=""), target.dataset$row.names),]
		target = target.vec = REVEALER.normalize.zero.one(target)
		phenotype = paste(gene.of.interest, "expression")
		target.class = ""
		target.type = "continuous"
	}
	
	if(!ifUseAllClasses){
#		browser()
		display.top.n.per.class = display.top.n
		ind <- order(target)
		target <- target[ind]
		target.vec <- target.vec[ind]
		m <- as.matrix(m)[, ind]
		sample.names = sample.names[ind]
#		nSample = dim(m)[2	]
		class.v <- CLS$class.v
		if( isCLS.input.cls ){
			if (is.vector(class.v)) {
				class.v <- class.v[ind]
			} else {
				class.v <- class.v[, ind]
			}
		}
		target.matrix = matrix(target, nrow=1)
	} else{ 
		annot.multipleClasses <- MI.multipleClasses <- AUC.multipleClasses <- 
				AUC.pval.multipleClasses <- t.stat.multipleClasses <- 
				t.pval.multipleClasses <- vector(length=display.top.n, 
						mode="numeric")
		m.multipleClasses = matrix(nrow=(display.top.n), ncol=nSample)
		rownames(m.multipleClasses) = rep("", display.top.n)
	}
	
	for( target.row in 1:length(target.matrix[,1])){
		target = target.matrix[target.row,]
		annot <- MI <- AUC <- AUC.pval <- t.stat <- t.pval <- vector(length=nGeneSets, mode="numeric")
		
		#NMI.ref <- mutual.inf.2(x = target, y = target, n.grid=100)#$NMI
#	browser()
		if( nGeneSets > 1){	
			MI = signif(REVEALER.mutual.inf.one.vs.many.efficient(target, as.matrix(m))$MI, 4)
		} else{
			target.bandwidth = bcv(target)
			MI = signif(REVEALER.mutual.inf.one.y.bandwidth(as.vector(m), target, target.bandwidth)/
							REVEALER.mutual.inf.one.y.bandwidth(target, target, target.bandwidth))
		}
		for (i in 1:nGeneSets) {
			if (nGeneSets == 1) {
				feature <- m
			} else {
				feature <- m[i,]
			}
#			MI[i] <- signif(mutual.inf.P(target, feature, n.grid=100)$NMI/NMI.ref, 4)
			if (target.type == "continuous") {
				AUC[i] <- AUC.pval[i] <- t.stat[i] <- t.pval[i] <- "-"
			} else if (target.type == "discrete") {
				feature.norm <- (feature - min(feature))/(max(feature) - min(feature))
				perf.auc <- roc.area(target, feature.norm)
				AUC[i] <- ifelse(perf.auc$A < 0.5, -(1 - perf.auc$A), perf.auc$A)
				AUC[i] <- signif(AUC[i], digits=4)
				p.val <- perf.auc$p.value
				p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
				AUC.pval[i] <- signif(p.val, digits=4)
				temp <- split(feature, target)
				x <- temp$'1'
				y <- temp$'0'
				t.stat[i] <- signif(t.test(x=x, y=y)$statistic, digits=4)
				p.val <- t.test(x=x, y=y)$p.value
				p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
				t.pval[i] <- signif(p.val, digits=4)
			}
			#annot[i] <- paste(MI[i], "     ", AUC[i], " (", AUC.pval[i], ")    ", t.stat[i], " (", t.pval[i], ")", sep="")
		}
		annot <- paste(MI, "     ", AUC, " (", AUC.pval, ")    ", t.stat, " (", t.pval, ")", sep="")
		
		if ((nGeneSets > 1) & (sort.results == T)) {
			MI.order <- order(MI, decreasing=T)
			MI <- MI[MI.order]
			AUC <- AUC[MI.order]
			AUC.pval <- AUC.pval[MI.order]
			t.stat <- t.stat[MI.order]
			t.pval <= t.pval[MI.order]
			m <- as.matrix(m)[MI.order,]
			annot <- annot[MI.order]
		}
		
#		if(!is.null(signatures.of.interest)){
#			print("calculating number of genesets with a higher rank")
		##			n.higher.genesets = min(which(row.names(m) %in% signatures.of.interest))-1
		##			n.higher.genesets = ifelse(n.higher.genesets > display.top.n.per.class, 
		##					display.top.n.per.class, n.higher.genesets)
#			n.genesets.to.plot.distribution = display.top.n.per.class
#			browser()
#		}
		
		if(ifUseAllClasses){
			
			from.ind = (display.top.n.per.class*(target.row-1)+1)
			to.ind = (display.top.n.per.class*(target.row))
			annot.multipleClasses[from.ind:to.ind] = annot[1:display.top.n.per.class]
			m.multipleClasses[from.ind:to.ind,] = m[1:display.top.n.per.class,]
			row.names(m.multipleClasses)[from.ind:to.ind] = row.names(m)[1:display.top.n.per.class]
			
			MI.multipleClasses[from.ind:to.ind] = MI[1:display.top.n.per.class]
			AUC.multipleClasses[from.ind:to.ind] = AUC[1:display.top.n.per.class]
			AUC.pval.multipleClasses[from.ind:to.ind] = AUC.pval[1:display.top.n.per.class]
			t.stat.multipleClasses[from.ind:to.ind] = t.stat[1:display.top.n.per.class]
			t.pval.multipleClasses[from.ind:to.ind] = t.pval[1:display.top.n.per.class]
			
			if( length(grep(".txt$", output.txt))!=0 ){
				this.class.output.txt = gsub(".txt", paste(".", classes[target.row], ".txt", sep=""), output.txt)	
			} else{
				this.class.output.txt = paste(".", output.txt, classes[target.row], ".txt", sep="")
			}
		} else{ this.class.output.txt = output.txt }
		annot2 <- data.frame(cbind(MI, AUC, AUC.pval, t.stat, t.pval))
		if( nGeneSets > 1){
			row.names(annot2) <- make.unique(row.names(m))
		} else row.names(annot2) = model.names
		write(paste(c("gene set ", noquote(colnames(annot2))), collapse="\t"), file = this.class.output.txt, append = F, 
				ncol = length(colnames(annot2)))
		write.table(annot2, file=this.class.output.txt, append=T, quote=F, sep="\t", eol="\n", col.names=F, row.names=T)
		
		print(paste("finished:", classes[target.row], "  which is", 
						target.row, "out of", length(target.matrix[,1]), "classes"))
	}
	
	if(ifUseAllClasses){
		annot = annot.multipleClasses
		m = m.multipleClasses
		if(ifSortWithinClasses){
#				browser()
			for( target.row in 1:length(target.matrix[,1])){
				target = target.matrix[target.row,]
				from.ind = (display.top.n.per.class*(target.row-1)+1)
				to.ind = (display.top.n.per.class*(target.row))
				
				class.sample.order = order(m[from.ind, which(target==1)])
				m[,which(target==1)] = m[,which(target==1)[class.sample.order]]
				sample.names[which(target==1)] = sample.names[which(target==1)][class.sample.order]
			}
		}
		
		
		MI = MI.multipleClasses
		AUC = AUC.multipleClasses
		AUC.pval = AUC.pval.multipleClasses
		t.stat = t.stat.multipleClasses
		t.pval = t.pval.multipleClasses
	} else if( !ifUseAllClasses && ifSortWithinClasses ){
		class.indicators = unique(target)
		for( class.indicator in class.indicators ){
			class.sample.order = order(m[1,which(target==class.indicator)])
			m[,which(target==class.indicator)] = m[,which(target==class.indicator)][class.sample.order]
			sample.names[which(target==class.indicator)] = sample.names[
					which(target==class.indicator)][class.sample.order]
		}
	}
	## Save ranks of the signatures of interest, if provided
	if( length(signatures.of.interest) > 0 ){
#		browser()
#		split.up = strsplit(signatures.of.interest, "_UP")
#		split.dn = strsplit(signatures.of.interest, "_DN")
#		browser()
		combined.signature = substr(signatures.of.interest[1], 1, nchar(signatures.of.interest[1])-3) 
		#split.up[[which(split.up %in% split.dn)]]
		signatures.of.interest = c(signatures.of.interest[1], combined.signature, signatures.of.interest[2])
		signature.ranks = which(rownames(m) %in% signatures.of.interest)
#		browser()
		signature.ranks.text = c(paste(na.omit(signature.ranks), rownames(m)[na.omit(signature.ranks)], 
						MI[na.omit(signature.ranks)]), 
				paste("(out of", length(m[,1]), "signatures)"))
		write(paste(weight, signif(MI[na.omit(signature.ranks)], digits=5), na.omit(signature.ranks), 
						rownames(m)[na.omit(signature.ranks)], sep="\t"), 
				file=output.MI.Ranks.txt, append=TRUE)
#		browser()
	} else{ signature.ranks.text = NULL}
	
	mycol <- vector(length=512, mode = "numeric")
	for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
	for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
	mycol <- rev(mycol)
	cex.axis = 1
	n.pinkogram.colors <- length(mycol)
	
#	browser()
	graphics.file.suffix = ""
	if(!is.null(highlight.tissue.name)){
		graphics.file.suffix = paste(graphics.file.suffix, "_", 
				paste(highlight.tissue.name, collapse="."),
				".highlighted", sep="")
	}
	if( ifUseAllClasses ){
		graphics.file.suffix = paste(graphics.file.suffix, ".allClasses", sep="")
	}
	if( ifSortWithinClasses){
		graphics.file.suffix = paste(graphics.file.suffix, ".sortedWithinClasses", sep="")
	}
	if( ifScoreBarPlot ){
		graphics.file.suffix = paste(graphics.file.suffix, ".barPlot", sep="")
	}
	
	
#	graphics.file.types = c("pdf", "eps", "jpeg", "bmp", "png")
#	browser()
	if(graphics.file.type == "pdf"){
		if( length(grep(".pdf$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".pdf", sep="")
		} else{ 
			output.graphics.filename = gsub(".pdf", 
					paste(graphics.file.suffix, ".pdf", sep=""), output.graphics.filename)}
		pdf(file=output.graphics.filename, height=graphic.height, width=graphic.width)
	} else if(graphics.file.type == "eps"){
		if( length(grep(".eps$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".eps", sep="")
		} else{ 
			output.graphics.filename= gsub(".eps", 
					paste(graphics.file.suffix, ".eps", sep=""), output.graphics.filename)}
		eps(file=paste(output.graphics.filename, ".eps", sep=""), height=graphic.height, width=graphic.width)
	} else if(graphics.file.type == "jpeg"){
		if( length(grep(".jpeg$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".jpeg", sep="")
		} else{ 
			output.graphics.filename = gsub(".jpeg", 
					paste(graphics.file.suffix, ".jpeg", sep=""), output.graphics.filename)}
		jpeg(file=paste(output.graphics.filename, ".jpg", sep=""), height=graphic.height, width=graphic.width)
	} else if(graphics.file.type == "bmp"){
		if( length(grep(".bmp$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".bmp", sep="")
		} else{ 
			output.graphics.filename = gsub(".bmp", 
					paste(graphics.file.suffix, ".bmp", sep=""), output.graphics.filename)}
		bmp(file=paste(output.graphics.filename, ".bmp", sep=""), height=graphic.height, width=graphic.width)
	} else if(graphics.file.type == "png"){
		if( length(grep(".png$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".png", sep="")
		} else{ 
			output.graphics.filename = gsub(".png", 
					paste(graphics.file.suffix, ".png", sep=""), output.graphics.filename)}
		png(file=paste(output.graphics.filename, ".png", sep=""), height=graphic.height, width=graphic.width)
	} else if(graphics.file.type == "tiff"){
		if( length(grep(".tiff$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".tiff", sep="")
		} else{ 
			output.graphics.filename = gsub(".tiff", 
					paste(graphics.file.suffix, ".tiff", sep=""), output.graphics.filename)}
		tiff(file=paste(output.graphics.filename, ".eps", sep=""), height=graphic.height, width=graphic.width)
	} else if(graphics.file.type=="window") {
		quartz(height=graphic.height, width=graphic.width)
	}
#	browser()
#	} else{
#		stop("Please specify heatmap output file type")
#	}
	if(ifMultipleTissues || ifSortWithinClasses && !ifScoreBarPlot ){
		print("ifMultipleTissues || ifSortWithinClasses && !ifScoreBarPlot")
		nf <- layout(matrix(c(1, 2, 3), ncol=1), 1, heights = c(2, 10, 2),
				respect = FALSE)
	} else if( !ifMultipleTissues && !ifSortWithinClasses && !ifScoreBarPlot ){	
		print("!ifMultipleTissues && !ifSortWithinClasses && !ifScoreBarPlot")
		nf <- layout(matrix(c(1, 2), ncol=1), 1, c(4, 10), respect=TRUE) 
	} else if(ifMultipleTissues || ifSortWithinClasses && ifScoreBarPlot ){
		print("ifMultipleTissues || ifSortWithinClasses && ifScoreBarPlot")
		nf <- layout(matrix(c(1, 2, 3, 0, 4, 3), ncol=2, byrow=FALSE), heights = c(2, 10, 2),
				width=c(12,1),
				respect = FALSE)
#		layout.show(4)
	} else if( !ifMultipleTissues && !ifSortWithinClasses && ifScoreBarPlot ){	
		print("!ifMultipleTissues && !ifSortWithinClasses && ifScoreBarPlot")
		nf <- layout(matrix(c(1, 2, 0, 3), ncol=2), heights=c(4, 10), widths=c(10,1), respect=FALSE)
#		layout.show(3)
	}
	max.v <- max(max(target), -min(target))
	if( isCLS.input.cls && !ifUseAllClasses ){
		V1 <- target
	} else if( isGCT.input.cls && !ifUseAllClasses ){
		target.mid.point = which.min(abs(target - quantile(target, 0.5)))
		V1 = ceiling(n.pinkogram.colors*target)
#		V1 = ceiling(c( (n.pinkogram.colors/2)*target[1:target.mid.point], 
#				(n.pinkogram.colors/2)*target[(target.mid.point+1):nSample] + n.pinkogram.colors/2))
	}
	if(ifUseAllClasses){
		V1 = colSums(target.matrix*1:n.classes)
	} 
	n.total.colors = ifelse( target.type=="discrete" && !ifMultipleTissues 
					&& !ifUseAllClasses, 1, length(class.colors))
	tissue.colors.plus.mycol = class.colors
	if( ifMultipleTissues && target.type=="continuous" ){
#		browser()
		V1 = t(rbind(V1, (tissue.labels + n.pinkogram.colors)))
		n.total.colors = n.pinkogram.colors + n.tissues
		tissue.colors.plus.mycol = c(mycol, tissue.colors)
		
	} else if( target.type == "continuous"){ n.total.colors = n.pinkogram.colors; tissue.colors.plus.mycol = mycol }
	if (ifTissueGenePlotTitle){
		plot.title = paste("tissue:", tissue, "  gene:", gene.of.interest)
	} else{
		plot.title = strsplit(output.pdf, "/")[[1]]
		if(is.full){
			plot.title = strsplit(plot.title[length(plot.title)], ".HEATMAP.FULL.pdf")[[1]]
		} else {
			plot.title = strsplit(plot.title[length(plot.title)], ".HEATMAP.pdf")[[1]]
		}
	}
	par(mar = c(1, ifelse(ifScoreBarPlot, 28, 15), 3, ifelse(ifScoreBarPlot, 18, 15)))
	image(1:nSample, 1:(ifelse(ifMultipleTissues, 2, 1)), 
			as.matrix(V1), zlim = c(0, n.total.colors), 
			col=tissue.colors.plus.mycol, axes=FALSE, 
			main=plot.title, sub = "", xlab= "", ylab="")
	if(!is.null(phenotype) && !is.null(target.class)){
		axis(2, at=1:1, labels=paste(phenotype, target.class), 
				adj= 0.5, tick=FALSE, las = 1, cex.axis=ifelse(ifScoreBarPlot, 1, 0.70), font.axis=1, line=-1)
	}
	axis(4, at=1:1, labels=paste("NMI     AUC (p-val)     t-test (p-val)", sep=""), adj= 0.5, tick=FALSE, 
			las = 1, cex.axis=ifelse(ifScoreBarPlot, 1, 0.80), font.axis=1, line=-1) 
	par(mar = c(5, ifelse(ifScoreBarPlot, 28, 15), 2, ifelse(ifScoreBarPlot, 18, 15)))
	if (display.top.n > nGeneSets) display.top.n <- nGeneSets
	if (nGeneSets == 1) {
		V1 <- m
		V1 <- (V1 - mean(V1))/sd(V1)
	} else {
		V1 <- m[1:display.top.n, ]
		for (i in 1:display.top.n) V1[i,] <- (V1[i,] - mean(V1[i,]))/sd(V1[i,])
	}
	max.v <- max(max(V1), -min(V1))
	V1 <- ceiling(n.pinkogram.colors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
	if (nGeneSets > 1) {
		V1 <- apply(V1, MARGIN=2, FUN=rev)
		image(1:nSample, 1:display.top.n, t(V1), zlim = c(0, n.pinkogram.colors), 
				col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
		left.labels = row.names(V1)
		heatmap.nrow = display.top.n
#		axis(2, at=1:display.top.n, labels=row.names(V1), adj= 0.5,
#				tick=FALSE, las = 1, cex.axis=ifelse(ifScoreBarPlot, 1, 0.70), font.axis=1, line=-1)
#		axis(4, at=1:display.top.n, labels=rev(annot[1:display.top.n]), adj= 0.5, tick=FALSE, 
#				las = 1, cex.axis=ifelse(ifScoreBarPlot, 1, 0.70), font.axis=1, line=-1)
#		if( nSample < 100){
#			axis(1, at=1:nSample, labels=sample.names, adj= 0.5, tick=FALSE, las = 3, 
#					cex.axis=ifelse(ifScoreBarPlot, 1, 0.60), font.axis=1, line=-1)
#		}
	} else {
		image(1:nSample, 1:1, as.matrix(V1), zlim = c(0, n.pinkogram.colors), 
				col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
		left.labels = model.names
		heatmap.nrow = 1
#		axis(2, at=1:1, labels=model.names, adj= 0.5, tick=FALSE, las = 1, 
#				cex.axis=ifelse(ifScoreBarPlot, 1, 0.70), font.axis=1, line=-1)
#		axis(4, at=1:1, labels=annot[1], adj= 0.5, tick=FALSE, las = 1, 
#				cex.axis=ifelse(ifScoreBarPlot, 1, 0.70), font.axis=1, line=-1)
#		if( nSample < 100){
#			axis(1, at=1:nSample, labels=sample.names, adj= 0.5, tick=FALSE, las = 3, 
#					cex.axis=ifelse(ifScoreBarPlot, 1, 0.60), font.axis=1, line=-1)
#		}
	}
	
	REVEALER.label.heatmap(
			heatmap.nrow = heatmap.nrow,
			heatmap.ncol = nSample,
			right.labels = annot,
			left.labels = left.labels,
			bottom.labels = sample.names
			)
	## Tissue Legend
	if(ifMultipleTissues){
		#browser()
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		legend(x=0, y= 1, #horiz = T, x.intersp = 0.5, y.intersp=.25, 
				legend=tissue.names, bty="n", xjust=0, yjust= 1, 
				fill = tissue.colors, cex = 0.8, #pt.cex=1.75, 
				ncol=4)
	}
	if( ifSortWithinClasses ){
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		legend(x=0, y= 1, #horiz = T, x.intersp = 0.5, y.intersp=.25, 
				legend=classes, bty="n", xjust=0, yjust= 1, 
				fill = class.colors, cex = 1.2, #pt.cex=1.75, 
				ncol=4)
	}
	if(ifScoreBarPlot){
		par(mar = c(3, 1, 0, 1))
		barplot(rev(MI[1:display.top.n]), xlab="MI scores\nbarplot", ylab="", xlim=c(-1,1), 
				axes=TRUE, horiz=TRUE, axisnames=FALSE, 
#				main="MI scores\nbarplot", font.main=1
		)
#		browser()
#		text("MI scores\nbarplot")
	}
	dev.off()
	
	
	
}

OPAM.Evaluate.Results.4.enrichment.plots <- function(  ## Added score bar plots to right from OPAM.Evaluate.Results.2
		input.ds,
		input.cls,
		phenotype = NULL,
		target.class = NULL,
		target.type = "discrete",
		sort.results = T,
		display.top.n = 20,
		output.txt,
		output.graphics.filename,
		signatures.of.interest = NULL,
		is.full = FALSE,
		statistic = NULL,
		weight = NULL,
		output.MI.Ranks.txt = NULL,
		graphics.file.type = "pdf", # pdf, eps, jpg / jpeg, bmp, png, tiff, window,
		graphic.height = 11,
		graphic.width = 17,
		gene.of.interest = NULL,
		tissue = "NA",
		highlight.tissue.name = NULL,  # highlight a single tissue in a dataset with multiple tissues
		ifMultipleTissues = FALSE,
		ifTissueGenePlotTitle = FALSE,
		ifUseAllClasses = FALSE,
		ifSortWithinClasses = FALSE,
		ifScoreBarPlot = FALSE,
		original.expression.ds = NULL,
		gene.set.databases = NULL,
		plot.gene.set.distributions = FALSE,
		rankings.data.dir = ""
) {
	
	isCLS.input.cls = (length(grep(".cls$", input.cls)) > 0)
	isGCT.input.cls = (length(grep(".gct$", input.cls)) > 0)
	
	dataset <- MSIG.Gct2Frame(filename = input.ds)  # Read gene expression dataset (GCT format)
	m <- data.matrix(dataset$ds)
	nGeneSets <- dim(m)[1]
	model.names <- dataset$row.names
	model.descs <- dataset$descs
	nSample <- dim(m)[2]
	for (i in 1:nGeneSets) {
		if ( !sum(is.na(m[i,])) && sd(as.matrix(m)[i,]) == 0) {
			val <- as.matrix(m)[i, 1]
			m[i,] <- as.matrix(m)[i,] + runif(n=nSample, 
					min= val - 0.001, max=val + 0.001)  # add small noise to flat profiles
		}
	}
#	browser()
	
	dim(m)
	sample.names <- dataset$names
	
	if( ifMultipleTissues ){
		tissue.type <- vector(length=nSample, mode="character")
#		temp = strsplit(sample.names, split="_")
		for (k in 1:nSample) {
			temp <- strsplit(sample.names[k], split="_") 
			tissue.type[k] <- paste(temp[[1]][2:length(temp[[1]])], collapse="_")
		}
		tissue.names = unique(tissue.type)
		tissue.labels = match(tissue.type, tissue.names)
		n.tissues = length(tissue.names)
		
		library(RColorBrewer)
		tissue.colors = c(brewer.pal(12, "Set3"), brewer.pal(12,"Paired"))[
				1:n.tissues]
		if( !is.null(highlight.tissue.name)){
			highlight.tissue.index = which(tissue.names %in% highlight.tissue.name)
#			tissue.labels = ifelse(tissue.labels == highlight.tissue.index, tissue.labels, 1)
			tissue.colors[-highlight.tissue.index] = c(brewer.pal(12, "Set3"), 
					brewer.pal(12,"Paired"))[2]  # 2nd color is yellow and is nice to have as the default background
			# as the other colors are darker
			tissue.colors[highlight.tissue.index] = c(brewer.pal(12, "Set3"), 
					brewer.pal(12,"Paired"))[1:(length(highlight.tissue.name)+1)][-2]
		}
		print(matrix(c(tissue.names, tissue.colors), ncol=2), quote=FALSE)
	}
#	else{
#		tissue.names = tissue
#		tissue.labels = rep(1, nSample)
#	}
	
	if( isCLS.input.cls ){
		CLS <- MSIG.ReadPhenFile.2(file = input.cls) # Read phenotype file (CLS format)
		cls.labels <- CLS$class.v
		cls.phen <- CLS$phen
		cls.list <- CLS$class.list 
		
#		browser()
		
		library(verification)
		
		if (is.null(phenotype)) {
			phen.loc <- 1
		} else {
			phen.loc <- match(phenotype, CLS$phen.names)
		}
		if (is.vector(CLS$class.list)) {
			target.vec <- CLS$class.list
		} else {
			target.vec <- CLS$class.list[phen.loc,]
		}
		if( ifUseAllClasses ){
			classes = unique(target.vec)
			target.matrix = t(sapply(classes, FUN=function(x) ifelse(target.vec ==x, 1, 0)))
			target.type = "discrete"
			n.classes = length(classes)
			display.top.n.per.class = ceiling(display.top.n/n.classes)
			display.top.n = display.top.n.per.class*n.classes
			if( length(CLS$col.phen) == 
					n.classes){  # If the CLS file provides enough colors to describe all the classes
				class.colors = CLS$col.phen
			} else {  # Otherwise, choose colors for the classes
				class.colors = c(#brewer.pal(12, "Set3"), brewer.pal(12,"Paired"), 
						brewer.pal(12, "Set3"), brewer.pal(12,"Paired"))[
						1:n.classes]
			}
		} else if (!ifUseAllClasses && target.type == "continuous") {
			target <- target.vec
		} else if (!ifUseAllClasses && target.type == "discrete") {
			target <- ifelse(target.vec == target.class, 1, 0) 
			class.colors = c("yellow", "purple")
			other.classes = paste(CLS$class.list[which(CLS$class.list != target.class)], collapse ="  ")
			classes = c( target.class, other.classes )
		}
	} else if( isGCT.input.cls && !is.null(gene.of.interest) ){
		target.dataset <- MSIG.Gct2Frame(filename = input.cls)
		target = data.matrix(target.dataset$ds)[grep(paste("^", 
								gene.of.interest, sep=""), target.dataset$row.names),]
		target = target.vec = REVEALER.normalize.zero.one(target)
		phenotype = paste(gene.of.interest, "expression")
		target.class = ""
		target.type = "continuous"
	}
	
	if(!ifUseAllClasses){
#		browser()
		display.top.n.per.class = display.top.n
		ind <- order(target)
		target <- target[ind]
		target.vec <- target.vec[ind]
		m <- as.matrix(m)[, ind]
		sample.names = sample.names[ind]
#		nSample = dim(m)[2	]
		class.v <- CLS$class.v
		if( isCLS.input.cls ){
			if (is.vector(class.v)) {
				class.v <- class.v[ind]
			} else {
				class.v <- class.v[, ind]
			}
		}
		target.matrix = matrix(target, nrow=1)
	} else{ 
		annot.multipleClasses <- MI.multipleClasses <- AUC.multipleClasses <- 
				AUC.pval.multipleClasses <- t.stat.multipleClasses <- 
				t.pval.multipleClasses <- vector(length=display.top.n, 
						mode="numeric")
		m.multipleClasses = matrix(nrow=(display.top.n), ncol=nSample)
		rownames(m.multipleClasses) = rep("", display.top.n)
	}
	
	
	########### --- Start MI Calculation loop on all target types in matrix -----	
	for( target.row in 1:length(target.matrix[,1])){
		target = target.matrix[target.row,]
		annot <- MI <- AUC <- AUC.pval <- t.stat <- t.pval <- vector(length=nGeneSets, mode="numeric")
		
		NMI.ref <- mutual.inf.2(x = target, y = target, n.grid=100)#$NMI
#	browser()
		for (i in 1:nGeneSets) {
			if (nGeneSets == 1) {
				feature <- m
			} else {
				feature <- m[i,]
			}
			MI[i] <- signif(mutual.inf.P(target, feature, n.grid=100)$NMI/NMI.ref, 4)
			if (target.type == "continuous") {
				AUC[i] <- AUC.pval[i] <- t.stat[i] <- t.pval[i] <- "-"
			} else if (target.type == "discrete") {
				feature.norm <- (feature - min(feature))/(max(feature) - min(feature))
				perf.auc <- roc.area(target, feature.norm)
				AUC[i] <- ifelse(perf.auc$A < 0.5, -(1 - perf.auc$A), perf.auc$A)
				AUC[i] <- signif(AUC[i], digits=4)
				p.val <- perf.auc$p.value
				p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
				AUC.pval[i] <- signif(p.val, digits=4)
				temp <- split(feature, target)
				x <- temp$'1'
				y <- temp$'0'
				t.stat[i] <- signif(t.test(x=x, y=y)$statistic, digits=4)
				p.val <- t.test(x=x, y=y)$p.value
				p.val <- ifelse(p.val > 0.5, 1 - p.val, p.val)
				t.pval[i] <- signif(p.val, digits=4)
			}
			annot[i] <- paste(MI[i], "     ", AUC[i], " (", AUC.pval[i], ")    ", t.stat[i], " (", t.pval[i], ")", sep="")
		}
		
		if ((nGeneSets > 1) & (sort.results == T)) {
			MI.order <- order(MI, decreasing=T)
			MI <- MI[MI.order]
			AUC <- AUC[MI.order]
			AUC.pval <- AUC.pval[MI.order]
			t.stat <- t.stat[MI.order]
			t.pval <= t.pval[MI.order]
			m <- as.matrix(m)[MI.order,]
			annot <- annot[MI.order]
		}
		
		#### --- Sort within classes ---
		if(ifUseAllClasses && ifSortWithinClasses){
#				browser()
#			for( target.row in 1:length(target.matrix[,1])){
#				target = target.matrix[target.row,]
			from.ind = (display.top.n.per.class*(target.row-1)+1)
			to.ind = (display.top.n.per.class*(target.row))
			
			class.sample.order = order(m[from.ind, which(target==1)])
			m[,which(target==1)] = m[,which(target==1)[class.sample.order]]
			sample.names[which(target==1)] = sample.names[which(target==1)][class.sample.order]		
		} else if( !ifUseAllClasses && ifSortWithinClasses ){
			class.indicators = unique(target)
			for( class.indicator in class.indicators ){
				#browser()
				class.sample.order = order(m[1,which(target==class.indicator)])
				m[,which(target==class.indicator)] = m[,which(target==class.indicator)][,class.sample.order]
				sample.names[which(target==class.indicator)] = sample.names[
						which(target==class.indicator)][class.sample.order]
			}
		}
		#### --- Sort within classes ---
		
		if(!is.null(signatures.of.interest) && plot.gene.set.distributions){
			print("calculating number of genesets with a higher rank")
#			n.higher.genesets = min(which(row.names(m) %in% signatures.of.interest))-1
#			n.higher.genesets = ifelse(n.higher.genesets > display.top.n.per.class, 
#					display.top.n.per.class, n.higher.genesets)
			n.genesets.to.plot.distribution = display.top.n.per.class
#			browser()
#			source("~/workspace/REVEALER/OPAM.library.v6.R")
			OPAM.project.dataset.3(  
					input.ds           = original.expression.ds,
					output.ds          = paste(original.expression.ds, ".PATHWAYS.top.from.MI.gct", sep=""),
					gene.set.databases = gene.set.databases,
					gene.set.selection = rownames(m)[1:n.genesets.to.plot.distribution],
					sample.norm.type   = "rank",  # "rank", "log" or "log.rank"
					output.score.type  = "ES",  # "ES" or "NES"
					weight             = 0.5,
					statistic          = "area.under.RES",
					combine.mode       = "combine.add",  # "combine.off", "combine.replace", "combine.add"
					nperm              = 1,
					correl.type        = "z.score",
					save.gene.set.rankings = TRUE,
					sample.names.for.ordering = colnames(m),
					rankings.data.dir = rankings.data.dir,
					target=target)             # "rank", "z.score", "symm.rank"
#			browser()
			top.genesets = rownames(m)[1:n.genesets.to.plot.distribution]
			write(top.genesets, file=paste(rankings.data.dir, 
							"top_", n.genesets.to.plot.distribution, "_genesets.txt", sep=""))
			genesets.w.files = sub("\\.expression\\.txt", "", 
					dir(path=rankings.data.dir, pattern="expression\\.txt"), perl=TRUE)
			gs.matched = match(top.genesets, genesets.w.files)
			genesets.vs.genes.list.this.target = vector(mode="list", length = n.genesets.to.plot.distribution)
			names(genesets.vs.genes.list.this.target) = top.genesets
			genesets.vs.genes.list.not.target = genesets.vs.genes.list.this.target
			
			
			expression = read.delim(paste(rankings.data.dir, genesets.w.files[1], ".expression.txt", sep=""))
			gene.names = NULL
			for(gs in top.genesets[which(!is.na(gs.matched))]){
				ranks = read.delim(paste(rankings.data.dir, gs, ".gene.ranks.txt", sep=""))
				genesets.vs.genes.list.this.target[[match(gs, top.genesets)]] = rowMeans(ranks[,which(target==1)])
				genesets.vs.genes.list.not.target[[match(gs, top.genesets)]] = rowMeans(ranks[,which(target==0)])
				gene.names = c(gene.names, rownames(ranks))
				expression.with.gene.ranks.heatmap(expression, 
						ranks,
						gene.set.name = gs,
						graphic.dir = rankings.data.dir,
						graphic.filetype = "png",
						normalize.function = "zero.one",
						target = target)
			}
			
			combined.gs = top.genesets[which(is.na(gs.matched))]
			for( gs in combined.gs ){
				dn.ranks = read.delim(paste(rankings.data.dir, gs, "_DN.gene.ranks.txt", sep=""))
				up.ranks = read.delim(paste(rankings.data.dir, gs, "_UP.gene.ranks.txt", sep=""))
				ranks = rbind(dn.ranks, up.ranks)
				genesets.vs.genes.list.this.target[[match(gs, top.genesets)]] = rowMeans(ranks[,which(target==1)])
				genesets.vs.genes.list.not.target[[match(gs, top.genesets)]] = rowMeans(ranks[,which(target==0)])
				gene.names = c(gene.names, rownames(ranks))
				expression.with.gene.ranks.heatmap(expression, 
						ranks,
						gene.set.name = gs,
						graphic.dir = rankings.data.dir,
						graphic.filetype = "png",
						normalize.function = "zero.one",
						target = target)
			}
			
			nGenes = dim(expression)[1]
			
			u.gene.names = unique(gene.names)
			genesets.vs.genes = matrix(nrow = n.genesets.to.plot.distribution, 
					ncol=length(u.gene.names), dimnames=list(top.genesets, u.gene.names))
#			gs.vs.g.this.target = lapply(genesets.vs.genes.list.this.target, c, 1, nGenes)
#			gs.vs.g.not.target = lapply(genesets.vs.genes.list.not.target, c, 1, nGenes)
			a = lapply(genesets.vs.genes.list.this.target, 
					function(x) 1-REVEALER.normalize.zero.one(c(x, 1, nGenes))[1:(length(x)-2)])
			mycol = rev(REVEALER.make.pinkogram.colors())
			ncol = length(mycol)
#			genesets.vs.genes.plot = genesets.vs.genes.norm * ncol
#			genesets.vs.genes.plot = genesets.vs.genes.plot+1
			a.plot = lapply(a, function(x) x*ncol+1)
			for(i in 1:length(genesets.vs.genes.list)){ 
				gs = a.plot[[i]]
				genesets.vs.genes[i,match(names(gs), u.gene.names)] = gs 
			}
			
			#genesets.vs.genes.norm = apply(genesets.vs.genes, 1, function(x) REVEALER.normalize.zero.one(na.omit(x)))
			genesets.vs.genes[is.na(genesets.vs.genes)] = 0
			
			mycol = c("black", mycol)
			library(gplots)
			browser()
			quartz(height=8.5, width=17); heatmap.2(genesets.vs.genes, col=mycol, trace="none", margins=c(5,15))
			#heatmap.2(genesets.vs.genes, col=mycol)
			quartz(width=28, height=8.5)
			par(mar=c(5, 15, 4, 2) + 0.1)
			image(1:length(u.gene.names), 1:n.genesets.to.plot.distribution, genesets.vs.genes.plot, 
					col=mycol, axes=FALSE, xlab="", ylab="")
			size.row.char = ifelse(n.genesets.to.plot.distribution >30, 30/(n.genesets.to.plot.distribution), 1)
			axis(2, at=1:n.genesets.to.plot.distribution, 
					labels=top.genesets, tick=FALSE, las = 1, 
					cex.axis=size.row.char, font.axis=1, line=-1)
			size.col.char <- 300/(length(u.gene.names) + 50)
			axis(1, at=1:length(u.gene.names), 
					labels=u.gene.names, tick=FALSE, las = 2, 
					cex.axis=size.col.char, font.axis=1, line=-1)
		}
		
		if(ifUseAllClasses){
			
			from.ind = (display.top.n.per.class*(target.row-1)+1)
			to.ind = (display.top.n.per.class*(target.row))
			annot.multipleClasses[from.ind:to.ind] = annot[1:display.top.n.per.class]
			m.multipleClasses[from.ind:to.ind,] = m[1:display.top.n.per.class,]
			row.names(m.multipleClasses)[from.ind:to.ind] = row.names(m)[1:display.top.n.per.class]
			
			MI.multipleClasses[from.ind:to.ind] = MI[1:display.top.n.per.class]
			AUC.multipleClasses[from.ind:to.ind] = AUC[1:display.top.n.per.class]
			AUC.pval.multipleClasses[from.ind:to.ind] = AUC.pval[1:display.top.n.per.class]
			t.stat.multipleClasses[from.ind:to.ind] = t.stat[1:display.top.n.per.class]
			t.pval.multipleClasses[from.ind:to.ind] = t.pval[1:display.top.n.per.class]
			
			if( length(grep(".txt$", output.txt))!=0 ){
				this.class.output.txt = gsub(".txt", paste(".", classes[target.row], ".txt", sep=""), output.txt)	
			} else{
				this.class.output.txt = paste(".", output.txt, classes[target.row], ".txt", sep="")
			}
		} else{ this.class.output.txt = output.txt }
		annot2 <- data.frame(cbind(MI, AUC, AUC.pval, t.stat, t.pval))
		if( nGeneSets > 1){
			row.names(annot2) <- make.unique(row.names(m))
		} else row.names(annot2) = model.names
		write(paste(c("gene set ", noquote(colnames(annot2))), collapse="\t"), file = this.class.output.txt, append = F, 
				ncol = length(colnames(annot2)))
		write.table(annot2, file=this.class.output.txt, append=T, quote=F, sep="\t", eol="\n", col.names=F, row.names=T)
		
		print(paste("finished:", classes[target.row], "  which is", 
						target.row, "out of", length(target.matrix[,1]), "classes"))
	}
	################################## --- End calculation of MI for all of target.matrix
	
	if(ifUseAllClasses){
		annot = annot.multipleClasses
		m = m.multipleClasses
#		if(ifSortWithinClasses){
		##				browser()
#			for( target.row in 1:length(target.matrix[,1])){
#				target = target.matrix[target.row,]
#				from.ind = (display.top.n.per.class*(target.row-1)+1)
#				to.ind = (display.top.n.per.class*(target.row))
#				
#				class.sample.order = order(m[from.ind, which(target==1)])
#				m[,which(target==1)] = m[,which(target==1)[class.sample.order]]
#				sample.names[which(target==1)] = sample.names[which(target==1)][class.sample.order]
#			}
#		}
		
		
		MI = MI.multipleClasses
		AUC = AUC.multipleClasses
		AUC.pval = AUC.pval.multipleClasses
		t.stat = t.stat.multipleClasses
		t.pval = t.pval.multipleClasses
	}
#	else if( !ifUseAllClasses && ifSortWithinClasses ){
#		class.indicators = unique(target)
#		for( class.indicator in class.indicators ){
#			class.sample.order = order(m[1,which(target==class.indicator)])
#			m[,which(target==class.indicator)] = m[,which(target==class.indicator)][class.sample.order]
#			sample.names[which(target==class.indicator)] = sample.names[
#					which(target==class.indicator)][class.sample.order]
#		}
#	}
	## Save ranks of the signatures of interest, if provided
	if( length(signatures.of.interest) > 0 ){
#		browser()
#		split.up = strsplit(signatures.of.interest, "_UP")
#		split.dn = strsplit(signatures.of.interest, "_DN")
#		browser()
		combined.signature = substr(signatures.of.interest[1], 1, nchar(signatures.of.interest[1])-3) 
		#split.up[[which(split.up %in% split.dn)]]
		signatures.of.interest = c(signatures.of.interest[1], combined.signature, signatures.of.interest[2])
		signature.ranks = which(rownames(m) %in% signatures.of.interest)
#		browser()
		signature.ranks.text = c(paste(na.omit(signature.ranks), rownames(m)[na.omit(signature.ranks)], 
						MI[na.omit(signature.ranks)]), 
				paste("(out of", length(m[,1]), "signatures)"))
		write(paste(weight, signif(MI[na.omit(signature.ranks)], digits=5), na.omit(signature.ranks), 
						rownames(m)[na.omit(signature.ranks)], sep="\t"), 
				file=output.MI.Ranks.txt, append=TRUE)
#		browser()
	} else{ signature.ranks.text = NULL}
	
	mycol <- vector(length=512, mode = "numeric")
	for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
	for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
	mycol <- rev(mycol)
	cex.axis = 1
	n.pinkogram.colors <- length(mycol)
	
#	browser()
	graphics.file.suffix = ""
	if(!is.null(highlight.tissue.name)){
		graphics.file.suffix = paste(graphics.file.suffix, "_", 
				paste(highlight.tissue.name, collapse="."),
				".highlighted", sep="")
	}
	if( ifUseAllClasses ){
		graphics.file.suffix = paste(graphics.file.suffix, ".allClasses", sep="")
	}
	if( ifSortWithinClasses){
		graphics.file.suffix = paste(graphics.file.suffix, ".sortedWithinClasses", sep="")
	}
	if( ifScoreBarPlot ){
		graphics.file.suffix = paste(graphics.file.suffix, ".barPlot", sep="")
	}
	
	
#	graphics.file.types = c("pdf", "eps", "jpeg", "bmp", "png")
#	browser()
	if(graphics.file.type == "pdf"){
		if( length(grep(".pdf$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".pdf", sep="")
		} else{ 
			output.graphics.filename = gsub(".pdf", 
					paste(graphics.file.suffix, ".pdf", sep=""), output.graphics.filename)}
		pdf(file=output.graphics.filename, height=graphic.height, width=graphic.width)
	} else if(graphics.file.type == "eps"){
		if( length(grep(".eps$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".eps", sep="")
		} else{ 
			output.graphics.filename= gsub(".eps", 
					paste(graphics.file.suffix, ".eps", sep=""), output.graphics.filename)}
		eps(file=paste(output.graphics.filename, ".eps", sep=""), height=graphic.height, width=graphic.width)
	} else if(graphics.file.type == "jpeg"){
		if( length(grep(".jpeg$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".jpeg", sep="")
		} else{ 
			output.graphics.filename = gsub(".jpeg", 
					paste(graphics.file.suffix, ".jpeg", sep=""), output.graphics.filename)}
		jpeg(file=paste(output.graphics.filename, ".jpg", sep=""), height=graphic.height, width=graphic.width)
	} else if(graphics.file.type == "bmp"){
		if( length(grep(".bmp$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".bmp", sep="")
		} else{ 
			output.graphics.filename = gsub(".bmp", 
					paste(graphics.file.suffix, ".bmp", sep=""), output.graphics.filename)}
		bmp(file=paste(output.graphics.filename, ".bmp", sep=""), height=graphic.height, width=graphic.width)
	} else if(graphics.file.type == "png"){
		if( length(grep(".png$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".png", sep="")
		} else{ 
			output.graphics.filename = gsub(".png", 
					paste(graphics.file.suffix, ".png", sep=""), output.graphics.filename)}
		png(file=paste(output.graphics.filename, ".png", sep=""), height=graphic.height, width=graphic.width)
	} else if(graphics.file.type == "tiff"){
		if( length(grep(".tiff$", output.graphics.filename)) == 0){
			output.graphics.filename = paste(output.graphics.filename, 
					graphics.file.suffix, ".tiff", sep="")
		} else{ 
			output.graphics.filename = gsub(".tiff", 
					paste(graphics.file.suffix, ".tiff", sep=""), output.graphics.filename)}
		tiff(file=paste(output.graphics.filename, ".eps", sep=""), height=graphic.height, width=graphic.width)
	} else if(graphics.file.type=="window") {
		quartz(height=graphic.height, width=graphic.width)
	}
#	browser()
#	} else{
#		stop("Please specify heatmap output file type")
#	}
	if(ifMultipleTissues || ifSortWithinClasses && !ifScoreBarPlot ){
#		print("ifMultipleTissues || ifSortWithinClasses && !ifScoreBarPlot")
		nf <- layout(matrix(c(1, 2, 3), ncol=1), 1, heights = c(2, 10, 2),
				respect = FALSE)
	} else if( !ifMultipleTissues && !ifSortWithinClasses && !ifScoreBarPlot ){	
#		print("!ifMultipleTissues && !ifSortWithinClasses && !ifScoreBarPlot")
		nf <- layout(matrix(c(1, 2), ncol=1), 1, c(4, 10), respect=TRUE) 
	} else if(ifMultipleTissues || ifSortWithinClasses && ifScoreBarPlot ){
#		print("ifMultipleTissues || ifSortWithinClasses && ifScoreBarPlot")
		nf <- layout(matrix(c(1, 2, 3, 0, 4, 3), ncol=2, byrow=FALSE), heights = c(2, 10, 2),
				width=c(12,1),
				respect = FALSE)
#		layout.show(4)
	} else if( !ifMultipleTissues && !ifSortWithinClasses && ifScoreBarPlot ){	
#		print("!ifMultipleTissues && !ifSortWithinClasses && ifScoreBarPlot")
		nf <- layout(matrix(c(1, 2, 0, 3), ncol=2), heights=c(4, 10), widths=c(10,1), respect=FALSE)
#		layout.show(3)
	}
	max.v <- max(max(target), -min(target))
	if( isCLS.input.cls && !ifUseAllClasses ){
		V1 <- target
	} else if( isGCT.input.cls && !ifUseAllClasses ){
		target.mid.point = which.min(abs(target - quantile(target, 0.5)))
		V1 = ceiling(n.pinkogram.colors*target)
#		V1 = ceiling(c( (n.pinkogram.colors/2)*target[1:target.mid.point], 
#				(n.pinkogram.colors/2)*target[(target.mid.point+1):nSample] + n.pinkogram.colors/2))
	}
	if(ifUseAllClasses){
		V1 = colSums(target.matrix*1:n.classes)
	} 
	n.total.colors = ifelse( target.type=="discrete" && !ifMultipleTissues 
					&& !ifUseAllClasses, 1, length(class.colors))
	tissue.colors.plus.mycol = class.colors
	if( ifMultipleTissues && target.type=="continuous" ){
#		browser()
		V1 = t(rbind(V1, (tissue.labels + n.pinkogram.colors)))
		n.total.colors = n.pinkogram.colors + n.tissues
		tissue.colors.plus.mycol = c(mycol, tissue.colors)
		
	} else if( target.type == "continuous"){ n.total.colors = n.pinkogram.colors; tissue.colors.plus.mycol = mycol }
	if (ifTissueGenePlotTitle){
		plot.title = paste("tissue:", tissue, "  gene:", gene.of.interest)
	} else{
		plot.title = strsplit(output.pdf, "/")[[1]]
		if(is.full){
			plot.title = strsplit(plot.title[length(plot.title)], ".HEATMAP.FULL.pdf")[[1]]
		} else {
			plot.title = strsplit(plot.title[length(plot.title)], ".HEATMAP.pdf")[[1]]
		}
	}
	par(mar = c(1, ifelse(ifScoreBarPlot, 28, 15), 3, ifelse(ifScoreBarPlot, 18, 15)))
	image(1:nSample, 1:(ifelse(ifMultipleTissues, 2, 1)), 
			as.matrix(V1), zlim = c(0, n.total.colors), 
			col=tissue.colors.plus.mycol, axes=FALSE, 
			main=plot.title, sub = "", xlab= "", ylab="")
	if(!is.null(phenotype) && !is.null(target.class)){
		axis(2, at=1:1, labels=paste(phenotype, target.class), 
				adj= 0.5, tick=FALSE, las = 1, cex.axis=ifelse(ifScoreBarPlot, 1, 0.70), font.axis=1, line=-1)
	}
	axis(4, at=1:1, labels=paste("NMI     AUC (p-val)     t-test (p-val)", sep=""), adj= 0.5, tick=FALSE, 
			las = 1, cex.axis=ifelse(ifScoreBarPlot, 1, 0.80), font.axis=1, line=-1) 
	par(mar = c(ifelse(ifScoreBarPlot, 10, 5), ifelse(ifScoreBarPlot, 28, 15), 2, ifelse(ifScoreBarPlot, 18, 15)))
	if (display.top.n > nGeneSets) display.top.n <- nGeneSets
	if (nGeneSets == 1) {
		V1 <- m
		V1 <- (V1 - mean(V1))/sd(V1)
	} else {
		V1 <- m[1:display.top.n, ]
		for (i in 1:display.top.n) V1[i,] <- (V1[i,] - mean(V1[i,]))/sd(V1[i,])
	}
	max.v <- max(max(V1), -min(V1))
	V1 <- ceiling(n.pinkogram.colors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
	if (nGeneSets > 1) {
		V1 <- apply(V1, MARGIN=2, FUN=rev)
		image(1:nSample, 1:display.top.n, t(V1), zlim = c(0, n.pinkogram.colors), 
				col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
		axis(2, at=1:display.top.n, labels=row.names(V1), adj= 0.5,
				tick=FALSE, las = 1, cex.axis=ifelse(ifScoreBarPlot, 1, 0.70), font.axis=1, line=-1)
		axis(4, at=1:display.top.n, labels=rev(annot[1:display.top.n]), adj= 0.5, tick=FALSE, 
				las = 1, cex.axis=ifelse(ifScoreBarPlot, 1, 0.70), font.axis=1, line=-1)
		if( nSample < 100){
			axis(1, at=1:nSample, labels=sample.names, adj= 0.5, tick=FALSE, las = 3, 
					cex.axis=ifelse(ifScoreBarPlot, 1, 0.60), font.axis=1, line=-1)
		}
	} else {
		image(1:nSample, 1:1, as.matrix(V1), zlim = c(0, n.pinkogram.colors), 
				col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
		axis(2, at=1:1, labels=model.names, adj= 0.5, tick=FALSE, las = 1, 
				cex.axis=ifelse(ifScoreBarPlot, 1, 0.70), font.axis=1, line=-1)
		axis(4, at=1:1, labels=annot[1], adj= 0.5, tick=FALSE, las = 1, 
				cex.axis=ifelse(ifScoreBarPlot, 1, 0.70), font.axis=1, line=-1)
		if( nSample < 100){
			axis(1, at=1:nSample, labels=sample.names, adj= 0.5, tick=FALSE, las = 3, 
					cex.axis=ifelse(ifScoreBarPlot, 1, 0.60), font.axis=1, line=-1)
		}
	}
	## Tissue Legend
	if(ifMultipleTissues){
		#browser()
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		legend(x=0, y= 1, #horiz = T, x.intersp = 0.5, y.intersp=.25, 
				legend=tissue.names, bty="n", xjust=0, yjust= 1, 
				fill = tissue.colors, cex = 0.8, #pt.cex=1.75, 
				ncol=4)
	}
	if( ifSortWithinClasses ){
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		legend(x=0, y= 1, #horiz = T, x.intersp = 0.5, y.intersp=.25, 
				legend=classes, bty="n", xjust=0, yjust= 1, 
				fill = class.colors, cex = 1.2, #pt.cex=1.75, 
				ncol=4)
	}
	if(ifScoreBarPlot){
		par(mar = c(8, 1, 0, 1))
		barplot(rev(MI[1:display.top.n]), xlab="MI scores\nbarplot", ylab="", xlim=c(-1,1), 
				axes=TRUE, horiz=TRUE, axisnames=FALSE, 
#				main="MI scores\nbarplot", font.main=1
		)
#		browser()
#		text("MI scores\nbarplot")
	}
	dev.off()
	
	
	
}

rnorm2d_xfixed = function (
		n, 
		rho.set = c(seq(0, 0.8, 0.05), seq(0.82, 0.98, 0.02)),
		random.seed = NULL, # = 7125809
		x = NULL,
		ifPrint = FALSE#,
)
#		distribution.type = c("gaussian", "gamma")) 
# randomly generate a vector, then generate random vectors with specified correlation
{
	
#	if(distribution.type == "gaussian"){
	distribution.function = rnorm
#	}
	if(!is.null(random.seed)) set.seed(random.seed)  # Ensure that the experiment is repeatable	
	if(is.null(x)) x = distribution.function(n)
	
	correl.set = vector("numeric", length=length(rho.set))
	y = matrix(nrow=n, ncol=length(rho.set))
	for( i in 1:length(rho.set) ){
#		print("---",quote=FALSE)
		rho.signed = rho.set[i]
		correl.sign = ifelse(rho.signed < 0, -1, 1)
		rho = abs(rho.signed)
#		mean = c(0, 0)
		sigma = diag(2)
		sigma[1, 2] = sigma[2, 1] = rho       # set off-diagonal values to rho (correlation)
#		print(sigma)
		L = chol(sigma)
#		print(L)
		xy = matrix(c(x, distribution.function(n)), nrow=n); 
		xy = xy %*% L
		xy[,2] = correl.sign*xy[,2]
		y[,i] = xy[,2]
		correl.set[i] = cor(xy[,1], xy[,2])
		if(ifPrint){
			print(paste("rho:", rho.signed, "  correl:", correl.set[i]))
		}
	}
	return(list(x=x, y=y, rho.set=rho.set, correl.set=correl.set))
}

REVEALER.make.pinkogram.colors = function(){
	mycol <- vector(length=512, mode = "numeric")
	for (k in 1:256) mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
	for (k in 257:512) mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
	mycol <- rev(mycol)
	return(mycol)
}

rgammaexp2d = function(
		n,
		theta.set=c(seq(0.05,1, 0.05), seq(2,20,1)),
		random.seed = NULL,
		x = NULL,
		ifPrint = FALSE
){
	x = runif(n)
	y = runif(n)
#	for( theta in theta.set ){
	gamma.x = matrix(nrow=n, ncol=length(theta.set))
	gamma.y = matrix(nrow=n, ncol=length(theta.set))
	gamma.xy = matrix(n,n)
	for(ind in 1:length(theta.set)){
		gamma.xy = sapply(y, function(Y) x^theta*exp(-x-x*Y)/gamma(theta))
#	}
		gamma.x[,ind] = rowSums(gamma.xy)  ## dimensions: n (rows) by length(theta.set) (columns)
		gamma.y[,ind] = colSums(gamma.xy)  ## dimensions: n (rows) by length(theta.set) (columns)
	}
	return(list(x=gamma.x, y=gamma.y, theta.set=theta.set))
}

normalize.for.pinkogram.plotting = function(v){
	max.v = max(v)
	return(v - (- max.v))/(1.001*(max.v - (- max.v)))
}

expression.with.gene.ranks.heatmap = function(expression.matrix,  ## cols = samples, rows = all genes
		gene.ranks.matrix,  ## cols = samples, rows = gene set genes
		gene.set.name = "",
		graphic.dir = NULL,
		graphic.filetype = c("png", "quartz", "pdf"),
		normalize.function = "z.score",
		plot.only.genes.in.gene.set = TRUE,
		target
## Plot gene ranks' expression heatmaps
## Author: Olga Botvinnnik
## Date: 07/08/2011
){
	library(MASS)
	
	### Make "Pinkogram" colors
	#browser()
	mycol <- REVEALER.make.pinkogram.colors()
	n.pinkogram.colors = length(mycol)
	
	n.samples = length(gene.ranks.matrix[1,])
	n.total.genes = length(expression.matrix[,1])
	n.genes.in.gene.set = length(gene.ranks.matrix[,1])
	
	if( normalize.function == "zero.one"){
		normalize.function.string = "(v - min(v))/(max(v) - min(v))"
		normalized.expression = apply(expression.matrix, 2, function(x) REVEALER.normalize.zero.one(x)
#					REVEALER.normalize.zero.one(rev(x))
		)
	} else if( normalize.function == "z.score"){
		normalize.function.string = "(x - mean(x)/sd(x)"
		normalized.expression = apply(expression.matrix, 2, function(x) (x - mean(x)/sd(x)))
	} else if (normalize.function == "pablo.pinkogram"){
		normalize.function.string = "v - (- max.v))/(1.001*(max.v - (- max.v))"
		normalized.expression = normalize.for.pinkogram.plotting(expression.matrix)
	}
	
	normalized.expression = ceiling(n.pinkogram.colors*normalized.expression) 
	normalized.expression = ifelse(normalized.expression==0, 1,normalized.expression)
	normalized.expression = apply(normalized.expression, 2, rev)
	
	#if(plot.only.genes.in.gene.set){
	filename.suffix = "_gs.genes.only"
	#} else{ filename.suffix = ""}
	#### ------ Make box and whisker plots of rank distribution ----
	#browser()
	pdf(height=8.5, width=11, file=paste(graphic.dir, gene.set.name, # "_", normalize.function, 
					filename.suffix, "_boxplots.pdf", sep=""))
	boxplot(gene.ranks.matrix, ylim=c(length(expression.matrix[,1]),1), show.names=FALSE, main=gene.set.name)
	size.col.char <- 30/(length(gene.ranks.matrix[1,]) + 15)
	axis(1, at=1:length(gene.ranks.matrix[1,]), labels=paste(colnames(gene.ranks.matrix), ""), tick=FALSE, las = 3, 
			cex.axis=size.col.char, font.axis=2, line=-1)
	dev.off()
	#### ------ End box and whisker plots of rank distribution ----
	
#	if( is.null(pdf.filename)){
	if(graphic.filetype == "png"){
		png(height=8.5, width=11, units="in", res=150,
				filename=paste(graphic.dir, gene.set.name, "_", normalize.function, 
						filename.suffix, ".png", sep=""))
	} else if(graphic.filetype == "quartz"){ 
		quartz(height=8.5,width=11)
	} else if(graphic.filetype == "pdf"){ 
		pdf(height=8.5,width=11, file=paste(graphic.dir, gene.set.name, "_", normalize.function, 
						filename.suffix, ".pdf", sep="")) }
#	par(mar=rep(0,4))
	par(mar=c(6,5,4,5))	
	
	if(plot.only.genes.in.gene.set){
		normalized.expression = sapply(1:n.samples, function(sample.ind) 
					normalized.expression[n.total.genes+1-gene.ranks.matrix[,sample.ind], sample.ind] )
		n.total.genes = n.genes.in.gene.set
		if(!missing(target)){
			gs.genes.MI = mutual.inf.3.v2(c(rep(0,24), rep(1,15)), normalized.expression)$MI
			gene.order = order(gs.genes.MI, decreasing=TRUE)
			normalized.expression = normalized.expression[rev(gene.order),]
			n.total.genes = n.genes.in.gene.set + 1
			target.classes = unique(target)
			if(length(target.classes)==2){
				target.for.plotting = ifelse(target == min(target.classes), 1, 2)+n.pinkogram.colors
				mycol = c(mycol, "grey", "black")
			} 
			## Add support for multiple classes - Olga 07/08/2011
#			else{
#				target.for.plotting = 
#			}
			normalized.expression = rbind(normalized.expression, target.for.plotting)
		} else{ gene.order = 1:n.genes.in.gene.set}
	}
#	browser()
	image(1:n.samples, 1:n.total.genes, t(normalized.expression), col=mycol, axes=FALSE, 
			main=paste(gene.set.name,"  normalized by:", normalize.function.string), 
			zlim=c(1, length(mycol)), xlab="", ylab="")
# Make lines to delineate samples
	for(x.ind in 1:(n.samples-1)){ 
		lines(c(x.ind+0.5, x.ind+0.5),c(0,n.total.genes+1), col="white", lwd=0.25)
	}
	if(plot.only.genes.in.gene.set){  # make lines to delineate genes
		for(y.ind in 1:(n.total.genes-1)){ 
			lines(c(0,n.samples+1), c(y.ind+0.5, y.ind+0.5), col="white", lwd=0.25)
		}
	}
	
# Make lines at the gene set ranks of each
	if(!plot.only.genes.in.gene.set){
		for(s.ind in 1:n.samples){ 
			for(rank.ind in n.total.genes+1-gene.ranks.matrix[,s.ind]){  
				# have to gene ranks subtract from total number of genes+1 because want 
				# rank of 1 to appear at the top position, which is the number of genes
				# Also, the +1 is because if the rank is 100 (out of 100 total genes), 
				# then 100-100 = 0, and R indexes from 1
				lines(c(s.ind-0.49, s.ind+0.49), c(rank.ind, rank.ind), lwd=1 )
			}
		}
	}
	
	## Label sample names along x-axis
	size.col.char <- 30/(n.samples + 15)
	axis(1, at=1:n.samples, labels=paste(colnames(gene.ranks.matrix), ""), tick=FALSE, las = 3, 
			cex.axis=size.col.char, font.axis=2, line=-1)
#	browser()
#	target.btwn.phenotypes.ind = max(which)
#	browser()
	# Plot gene names and MI values
	if(plot.only.genes.in.gene.set){
		size.row.char = ifelse(n.total.genes >30, 30/(n.total.genes), 1)
		axis(2, at=1:n.total.genes, 
				labels=rev(paste(c("target", rownames(gene.ranks.matrix)[gene.order]), "")), tick=FALSE, las = 1, 
				cex.axis=size.row.char, font.axis=2, line=-1)	
		axis(4, at=1:n.total.genes, 
				labels=rev(c(" 1", paste("", signif(gs.genes.MI[gene.order], digits=4)))), 
				tick=FALSE, las = 1, 
				cex.axis=size.row.char, font.axis=1, line=-1)
		## Delineate between top MI's (in our experience, >0.05, <-0.05 are the most significant)
		have.genes.above.0.05 = sum(rev(gs.genes.MI[gene.order])>0.05)
		if( have.genes.above.0.05){
			MI.above.0.05.ind = min(which(rev(gs.genes.MI[gene.order])>0.05))-1
			lines(c(0,n.samples+1), c(MI.above.0.05.ind+0.5, MI.above.0.05.ind+0.5), col="black", lwd=1)
#			lines()
			# Lines to delinate target vs non-target
			target.btwn.phenotypes.ind = max(which(target==target[1]))
			lines(rep(target.btwn.phenotypes.ind+0.5, 2), c(MI.above.0.05.ind+0.5,n.total.genes+1), col="black", lwd=1)
		}
		have.genes.below.neg0.05 = sum(rev(gs.genes.MI[gene.order])<(-0.05))
		if( have.genes.below.neg0.05){
			MI.below.neg0.05.ind = max(which(rev(gs.genes.MI[gene.order])<(-0.05)))
			lines(c(0,n.samples+1), c(MI.below.neg0.05.ind+0.5, MI.below.neg0.05.ind+0.5), col="black", lwd=1)
			# Lines to delinate target vs non-target
			target.btwn.phenotypes.ind = max(which(target==target[1]))
			lines(rep(target.btwn.phenotypes.ind+0.5, 2), c(0, MI.below.neg0.05.ind+0.5), col="black", lwd=1)
		}
	}
#	browser()
	
	if(graphic.filetype != "quartz"){
		dev.off()
	}
}

REVEALER.load.libraries = function(library.names){
	### loads multiple libraries in one line. library.names is a character vector
	### of the library names you'd like to load
#	if(is.missing(libary.names)){ 
#		stop("No libary names specified")
#	} else 
	if( length(library.names) > 1 ){
		sapply(library.names, library, character.only = TRUE)
	} else{ library(library.names, character.only = TRUE)}
}

REVEALER.parse.gct = function(filename){
	dataset = MSIG.Gct2Frame(filename)
#	daseset.matrix <- data.matrix(dataset$ds)
#	dataset.row.names <- dataset$row.names
#	dataset.nRows = length(dataset.matrix[1,])
#	dataset.col.names = dataset$names
#	dataset.nCols = length(dataset.col.names)
	return(list(
					matrix = data.matrix(dataset$ds),
					row.names = dataset$row.names,
					nRow = length(dataset$ds[,1]),
					sample.names = dataset$names,
					nSample = length(dataset$names)
			))
}

REVEALER.parse.cls = function(filename){
	## Still needs to be able to read non-numeric CLS files
	dataset = MSIG.ReadPhenFile.2(filename)
	dataset.matrix <- apply(dataset$class.list, 2, as.numeric)
	colnames(dataset.matrix) = colnames(data.frame(dataset$class.list))
	rownames(dataset.matrix) = dataset$phen.names
#	dataset.row.names <- dataset$phen.names
#	dataset.nRows = length(dataset.matrix[1,])
#	dataset.col.names = dataset$names
#	dataset.nCols = length(dataset.col.names)
	return(list(
					matrix = dataset.matrix,
					row.names = rownames(dataset.matrix),
					nRow = length(dataset.matrix[,1]),
					sample.names = colnames(dataset.matrix),
					nSample = length(dataset.matrix[1,])
			))
}

REVEALER.parse.txt.sample.rows = function(filename){
	## ambiguous tab-delimited file: assume samples are rows.
	## Have to correct for this and re-orient the data so samples are columns
	dataset = read.delim(filename)
#	rownames(dataset) = dataset[,1]
	dataset.matrix = t(as.matrix(dataset))
	if( is.null(colnames(dataset.matrix)) ){
		dataset.matrix = t(as.matrix(dataset[,-1]))
		colnames(dataset.matrix) = as.vector(dataset[,1])
		rownames(dataset.matrix) = colnames(dataset)[-1]
	}
#	browser()
#	daseset.matrix <- data.matrix(dataset$ds)
#	dataset.row.names <- dataset$row.names
#	dataset.nRows = length(dataset.matrix[1,])
#	dataset.col.names = dataset$names
#	dataset.nCols = length(dataset.col.names)
#	browser()
	return(list(
					matrix = dataset.matrix,
					row.names = rownames(dataset.matrix),
					nRow = length(dataset.matrix[,1]),
					sample.names = colnames(dataset.matrix),
					nSample = length(dataset.matrix[1,])
			))
}

REVEALER.parse.txt.sample.cols = function(filename){
	## ambiguous tab-delimited file: assume samples are rows.
	## Have to correct for this and re-orient the data so samples are columns
	dataset = read.delim(filename)
#	browser()
#	rownames(dataset) = dataset[,1]
	dataset.matrix = as.matrix(dataset)
#	daseset.matrix <- data.matrix(dataset$ds)
#	dataset.row.names <- dataset$row.names
#	dataset.nRows = length(dataset.matrix[1,])
#	dataset.col.names = dataset$names
#	dataset.nCols = length(dataset.col.names)
	return(list(
					matrix = dataset.matrix,
					row.names = rownames(dataset.matrix),
					nRow = length(dataset.matrix[,1]),
					sample.names = colnames(dataset.matrix),
					nSample = length(dataset.matrix[1,])
			))
}

REVEALER.exclude.samples = function(ds, exclude.sample.names){
	exclude.samples.ind = which(ds$sample.names %in% exclude.sample.names)
	if(length(exclude.samples.ind) > 0 ){
		ds$sample.names = ds$sample.names[-exclude.samples.ind]
		ds$matrix = ds$matrix[, -exclude.samples.ind]
		if(ds$nRow == 1){
			ds.matrix = t(as.matrix(ds.matrix))
		}
		ds$nSample = length(ds$matrix[1,])
	}
	return(ds)
}

REVEALER.exclude.rows = function(ds, exclude.row.names){
#	browser() 
	exclude.row.ind = which(ds$row.names %in% exclude.row.names)
	if(length(exclude.row.ind) > 0 ){
		ds$row.names = ds$row.names[-exclude.row.ind]
		ds$matrix = ds$matrix[-exclude.row.ind,]
		ds$nRow = length(ds$row.names)
		
		if(ds$nRow == 1){
			ds.matrix = t(as.matrix(ds.matrix))
		}
	}
	return(ds)
}

REVEALER.fix.results.dir = function(results.dir){
	if( !is.null(results.dir) ){
		#browser()
		if( length(grep("./$", results.dir)) == 0 ){
			results.dir = paste(results.dir, "/", sep="")
		}
	}
	return(results.dir)
}

REVEALER.extract.row.names.to.match = function(ds.original, row.names.to.match){
	ds = ds.original
	if( is.null(row.names.to.match[1]) || is.na(row.names.to.match[1]) ){
		stop("Must provide a vector of target names to evaluate, or specify 'ALL'")
	}
	## Remove "Commented out" target.names.to.match (with "#" at beginning of name)
	if( length(grep("^#", row.names.to.match)) > 0){
		row.names.to.match = row.names.to.match[-grep("^#", row.names.to.match)]
	}
	if( row.names.to.match[1] == "ALL"){
		ds$row.names = ds.original$row.names
		ds$matrix = ds.original$matrix
#		model.descs = dataset.all$descs
	} else if( row.names.to.match[1] != "ALL"){
		ds$row.names = row.names.to.match
		#browser()
		row.names.ind = match(row.names.to.match, ds.original$row.names)
		ds$matrix= ds.original$matrix[row.names.ind,]
#		model.descs = dataset.all$descs[model.ind]
#	browser()
	}
	if( length(row.names.ind) == 1 ){
		ds$matrix = t(as.matrix(ds$matrix))
		rownames(ds$matrix) = ds$row.names
	}
	ds$nSample = ds.original$nSample
	ds$nRow = length(ds$row.names)
	
	
	return(ds)
}

REVEALER.check.multiple.tissues = function(ds, show.multiple.tissues, tissue){
	tissues = vector(length=2, mode="list")
	if( show.multiple.tissues ){
		## Get tissue labels for each sample for heatmap plotting
		tissue.type <- vector(length=ds$nSample, mode="character")
		for (k in 1:ds$nSample) {
			temp <- strsplit(ds$sample.names[k], split="_") 
			tissue.type[k] <- paste(temp[[1]][2:length(temp[[1]])], collapse="_")
		}
		tissues$names = unique(tissue.type)
		tissues$labels = match(tissue.type, tissue.names)
	} else{
		tissues$names = tissue
		tissues$labels = rep(1, ds$nSample)
	}
	return(tissues)
}

REVEALER.get.extracted.file.prefix = function(filename, filetype){
	temp <- strsplit(filename, split="/") # Extract test file name
	s <- length(temp[[1]])
	test.file.name <- temp[[1]][s]
	temp <- strsplit(test.file.name, split=paste(".", substr(filetype, 1,3)))
	extracted.file.prefix <-  temp[[1]][1]
	return(extracted.file.prefix)
}

REVEALER.remove.na.rows = function(target.ds){
	## Remove NA rows
	na.row.index = unique(unlist(as.vector(apply(target.ds$matrix, MARGIN=2, FUN=function(x) which(is.na(x))))))
	if(length(na.row.index) > 0 ) { 
		target.ds$removed.na.row.index = na.row.index
		target.ds$removed.na.row = target.ds$matrix[na.row.index,]
		target.ds$removed.na.row.nam = target.ds$row.names[na.row.index]
		target.ds$matrix = target.ds$matrix[-na.row.index,] 
		target.ds$row.names = target.ds$row.names[-na.row.index]#rownames(target.ds$matrix)
		target.ds$nRow = length(target.ds$row.names)#dim(target.ds$matrix)[1]
		if( target.ds$nRow == 1 ){
			target.ds$matrix = t(as.matrix(target.ds$matrix))
		}
	}
	return(target.ds)
}

REVEALER.remove.na.sample = function(target.ds){
	## Remove NA rows
	na.col.index = unique(unlist(as.vector(apply(target.ds$matrix, MARGIN=1, FUN=function(x) which(is.na(x))))))
	if(length(na.col.index) > 0 ) { 
#		print(na.col.index)
#		browser()
		target.ds$removed.na.sample.index = na.col.index
		target.ds$removed.na.sample = target.ds$matrix[,na.col.index]
		target.ds$removed.na.sample.name = target.ds$sample.names[na.col.index]
		target.ds$matrix = target.ds$matrix[,-na.col.index]
		target.ds$sample.names = target.ds$sample.names[-na.col.index]#colnames(target.ds$matrix)
		target.ds$nSample = length(target.ds$sample.names)
		
		if(target.ds$nRow == 1 ){
			target.ds$matrix = t(as.matrix(target.ds$matrix))
		}
	} else{ target.ds$removed.na.sample = target.ds$removed.na.sample.index = NA }
	return(target.ds)
}

REVEALER.normalize.scores = function(target.ds, normalize.score, normalization.type){
# normalize scores	
	if (normalize.score == T) {
		if (normalization.type == "zero.one") {
			#browser()
#			for (i in 1:target.ds$nRow) {
			apply(target.ds$matrix, 1, function(m)
						m <- (m - min(m))/(max(m) - min(m)))
#			}
		} else if (normalization.type == "z.score") {
#			for (i in 1:target.ds$nRow) {
			apply(target.ds$matrix, 1, function(m)
						m <- (m - mean(m))/sd(m))
#			}
		} else if (normalization.type == "r.z.score") {
#			for (i in 1:target.ds$nRow) {
			apply(target.ds$matrix, 1, function(m)
						m <- (m - median(m))/mad(m))
#			}
		}
	}
	return(target.ds)
}

REVEALER.add.feature.annotation = function(feature.ds.original, feature.seed.names.to.match,
		feature.annotation.file){
	### This function still needs to be tested - should look up and make sure the names
	### in the annotation file and feature row names match up
	#browser()
	
	if( length(grep("^#", feature.seed.names.to.match)) > 0){
		feature.seed.names.to.match = feature.seed.names.to.match[-grep("^#", feature.seed.names.to.match)]
	}
	seed.ind.in.original = which(feature.ds.original$row.names %in% feature.seed.names.to.match)
	
#	feature.ds = feature.ds.original
	
	feature.annotation.df = read.delim(feature.annotation.file)
#	if(is.null(rownames(feature.annotation.df))){
	feature.ds.original$row.names = paste(feature.ds.original$row.names,
			feature.annotation.df[,-1], sep="_" )
#	}
	
	feature.seed.names.names.to.match = feature.ds.original$row.names[seed.ind.in.original]
	feature.ds.original$feature.seed.names.to.match = feature.seed.names.to.match
	
	return(feature.ds.original)
#	
#	if(feature.suffix.or[1] != ""){
#		feature.seed.names.to.match.or = paste(feature.seed.names.to.match, feature.suffix.or[1], sep="")
#		
#	} else{ feature.seed.names.to.match.or = NULL }
#	if ( feature.suffix.and[1] != FALSE){
#		feature.seed.names.to.match.and = c( feature.seed.names.to.match, 
#				paste(feature.seed.names.to.match, feature.suffix.and, sep=""))
#	} else{feature.seed.names.to.match.and = NULL }
#	
#	if( feature.suffix.and[1] != FALSE || feature.suffix.or[1] != "" ){
#		feature.seed.names.to.match = c(feature.seed.names.to.match.or, feature.seed.names.to.match.and)
#	} #else{ feature.ds.original$row.names = c( feature.seed.names.to.match ) }
#	#browser()
#	#print(feature.seed.names.to.match)
#	
#	## Find chromosomal locations of genes specified
#	## See "if( find.chromosomal.locations)" for more transparent code
#	if( feature.add.chrom.locs ){
#		library(org.Hs.eg.db)
#		
	##			browser()
#		feature.split = strsplit(feature.seed.names.to.match, split="_")
#		feature.noampdel = unlist( lapply(feature.split, function(x) x[1]))
#		feature.egIDs = mget(feature.noampdel, org.Hs.egALIAS2EG, ifnotfound=NA)
#		feature.no.chrom.loc = which(lapply(lapply(feature.egIDs, is.na), sum) > 0)
#		if( length(feature.no.chrom.loc) == 0){
#			feature.egIDs.w.locs = feature.egIDs
#		} else{
#			feature.egIDs.w.locs = feature.egIDs[-feature.no.chrom.loc]
#		}
#		feature.locs.list = lapply(feature.egIDs.w.locs, mget, org.Hs.egMAP)
#		feature.locs = vector(mode="character", length=length(feature.ds$row.names))
#		if( length(feature.no.chrom.loc) == 0){
#			feature.locs = unlist(lapply(feature.locs.list, function(x) paste(unlist(x), collapse="_")))
#		} else{
#			feature.locs[feature.no.chrom.loc] = "NA"
#			feature.locs[-feature.no.chrom.loc] = (
#						unlist(lapply(feature.locs.list, function(x) paste(unlist(x), collapse="_"))) )
#		}
#		feature.w.locs = paste(feature.noampdel, ".", feature.locs, sep="")
#		feature.ampdel.suffix = unlist(lapply(feature.split, function(x) x[2]))
#		feature.seed.names.to.match = paste(feature.w.locs, "_", feature.ampdel.suffix, sep="")
#		
#		print("Found chromosomal locations for all genes in feature.seed.names.to.match!")
#		
#		
#		feature.names.split = strsplit(feature.ds.original$row.names, split="_")
#		feature.names.noampdel = unlist( lapply(feature.names.split, function(x) x[1]))
#		feature.names.egIDs = mget(feature.names.noampdel, org.Hs.egALIAS2EG, ifnotfound=NA)
#		feature.names.no.chrom.loc = which(lapply(lapply(feature.names.egIDs, is.na), sum) > 0)
#		if( length(feature.names.no.chrom.loc) == 0){
#			feature.names.egIDs.w.locs = feature.names.egIDs
#		} else{
#			feature.names.egIDs.w.locs = feature.names.egIDs[-feature.names.no.chrom.loc]
#		}
#		feature.names.locs.list = lapply(feature.names.egIDs.w.locs, FUN=mget, org.Hs.egMAP)
#		feature.names.locs = vector(mode="character", length=length(feature.names))
#		if( length(feature.names.no.chrom.loc) == 0){
#			feature.names.locs = unlist(lapply(feature.names.locs.list, function(x) paste(unlist(x), collapse="_")))
#		} else{
#			feature.names.locs[feature.names.no.chrom.loc] = "NA"
#			feature.names.locs[-feature.names.no.chrom.loc] = (
#						unlist(lapply(feature.names.locs.list, function(x) paste(unlist(x), collapse="_"))) )
#		}
#		feature.names.w.locs = paste(feature.names.noampdel, ".", feature.names.locs, sep="")
#		feature.names.ampdel.suffix = unlist(lapply(feature.names.split, function(x) x[2]))
#		#feature.names.no.suffix = which(is.na(feature.names.ampdel.suffix))
#		#feature.names[feature.names.no.suffix] = feature.names.w.locs[feature.names.no.suffix]
	##			feature.names[-feature.names.no.suffix] = paste(feature.names.w.locs[-feature.names.no.suffix], 
	##					"_", feature.names.ampdel.suffix[-feature.names.no.suffix], sep="")
#		feature.ds.original$row.names = paste(feature.names.w.locs, "_", feature.names.ampdel.suffix, sep="")
#		print("Found chromosomal locations for all genes in cls file!")
#		
#		rm(list=c("feature.names.split", "feature.names.noampdel", "feature.names.egIDs", "feature.names.no.chrom.loc", 
#						"feature.names.egIDs.w.locs", "feature.names.locs.list", "feature.names.locs", "feature.names.w.locs",
#						"feature.names.ampdel.suffix", "feature.split", "feature.noampdel", 
#						"feature.egIDs", "feature.no.chrom.loc", 
#						"feature.egIDs.w.locs", "feature.locs.list", "feature.locs", "feature.w.locs",
#						"feature.ampdel.suffix"))
#	}
#	#print(feature.seed.names.to.match)
#	feature.ds.original$feature.seed.names.to.match = feature.seed.names.to.match
#	return(feature.ds.original)
}

REVEALER.preprocessing = function(
		target.ds, 
		feature.ds, 
		make.inbetween.heatmaps, 
		show.multiple.tissues,
		pdf.height, 
		pdf.width, 
		results.dir, 
		test.file.prefix, 
		file.suffix,
		decreasing.order,
		tissue.names,
		tissue.labels,
		if.feature.summary.on.top,
		popup.heatmap = FALSE,
		kde2d.timing = FALSE,
		target.if.unmatched,
		draw.white.lines = FALSE
)
### Aligns and plots the targets vs the features. Sorts by normalized mutual information.

{
#	phen.names = c("SUMMARY", phen.names)
#	browser()
#	ind.phen.pass1 = which( feature.ds$row.names %in% feature.seed.names.to.match )
#	## I tried to do this as a tryCatch but couldn't get it to work, so I kept it simple
#	if(length(ind.phen.pass1)==0){
#		stop(simpleError(paste("The gene names provided do not match with",
#								"the feature file. Some common mistakes:",
#								"the genes must either have '_MUT', '_AMP',",
#								"or '_DEL' suffixed to them, or", 
#								"'feature.suffix.or' and/or 'feature.suffix.and' set to TRUE")))
#	}
#	phen.pass1 = phen.names[ind.phen.pass1]
#	n.phen.pass1 = length(phen.pass1)
	MI.targets.vs.feature.summary = vector( length=target.ds$nRow, mode="numeric" )
	feature.summary = colSums(feature.ds$matrix)
	
	if(kde2d.timing){
		mutual.inf.one.vs.one = mutual.inf.2.kde2d.olga
		mutual.inf.one.vs.many = REVEALER.mutual.inf.one.vs.many
	} else{
		mutual.inf.one.vs.one = mutual.inf.2
		mutual.inf.one.vs.many = REVEALER.mutual.inf.one.vs.many.efficient#mutual.inf.3.v2
	}
	
	if(make.inbetween.heatmaps){
		filename <- paste(results.dir, test.file.prefix, file.suffix, ".pre.calculations", sep="")
		if( popup.heatmap ){
			quartz(height = pdf.height, width = pdf.width )
		} else{
			pdf(file=paste(filename, ".pdf", sep=""), height = pdf.height, width = pdf.width )
		}
#		browser()
#quartz(height = 11, width = 17)
		if( show.multiple.tissues ){
			print("show.multiple.tissues doesn't work yet!")
			#browser()
			#quartz(height = 11, width = 17)
#			MSIG.HeatMapPlot.10.show.multiple.tissues(V = target.ds$matrix, 
#					pathway.mut = feature.summary,
#					row.names = target.ds$row.names,
#					col.labels = cls.labels.pass1,
#					col.classes = cls.phen2.pass1, 
#					phen.cmap = colors.list, 
#					phen.names = phen.pass1,
#					col.names = sample.names, 
#					main = paste(tissue, "- Initial Heatmap ('Step 0')"), 
#					xlab="  ", ylab="  ", row.norm = row.norm,  
#					cmap.type = cmap.type, char.rescale = char.rescale,  legend=ifLegend,
#					tissue.names = tissue.names,
#					tissue.labels = tissue.labels,
#					highlight.tissue.name = highlight.tissue.name)
		} else{
#			browser()
			REVEALER.preprocessing.heatmap(
					target.ds,
					feature.ds,
					feature.summary,
					if.feature.summary.on.top = if.feature.summary.on.top,
					main = filename,
					draw.white.lines = draw.white.lines#,
#					if.greyscale.feature.summary = FALSE,
#					if.greyscale.feature.ds = FALSE
			)
		}
		dev.off()
	}
	
#	model.descs2.pass1 = vector(length = target.ds$nRow, mode="character")
	
	#### ---- Normalized Mutual Information of targets vs feature.summary ----
	if( length(unique(feature.summary)) > 1 ){
		
		if( target.ds$nRow > 1 ){
#			browser()
			MI.targets.vs.feature.summary = 
#					mutual.inf.3.v2(
#					REVEALER.mutual.inf.one.vs.many(
					mutual.inf.one.vs.many(
							feature.summary, 
							target.ds$matrix, 
					#target.vector.name="SUMMARY", 
					#tissue=tissue
					)$MI
			target.row.order = order(MI.targets.vs.feature.summary,
					decreasing = decreasing.order)
			target.ds$matrix <- target.ds$matrix[target.row.order, ]
			sample.order <- order(target.ds$matrix[1,], 
					decreasing = decreasing.order )
			target.ds$matrix <- target.ds$matrix[, sample.order]
		} else{ 
			#browser()
			MI.summary.vs.summary = 
#					mutual.inf.2(
					mutual.inf.one.vs.one(
#					mutual.inf.2.kde2d.olga(		
							feature.summary, feature.summary)
			MI.targets.vs.feature.summary =
					mutual.inf.one.vs.one(
#					mutual.inf.2(
#					mutual.inf.2.kde2d.olga(
							feature.summary, target.ds$matrix[1,])/MI.summary.vs.summary
#			model.descs2.pass1 <- signif(MI.results, digits=3)
#			m2.pass1 <- m ; 
			target.row.order = 1
			sample.order <- order(target.ds$matrix[1,], 
					decreasing = decreasing.order )
			target.ds$matrix <- t(as.matrix(target.ds$matrix[, sample.order]))
#			rownames(m2.pass1) = model.names
		}
	} else{ 
		MI.targets.vs.feature.summary = rep(NA, target.ds$nRow)
#		FDR.list.pass1 = rep(NA, target.ds$nRow)
#		model.descs2.pass1 = rep(" - (FDR = - )", target.ds$nRow)
		if( target.ds$nRow > 1 ){
			loc <- match(target.if.unmatched, target.ds$row.names)
			sample.order <- order(target.ds$matrix[loc,], decreasing = decreasing.order)
			target.ds$matrix <- target.ds$matrix[, sample.order]
			correl <- cor(t(target.ds$matrix))[, loc]
			target.row.order <- order(correl, decreasing=decreasing.order)
			target.ds$matrix <- target.ds$matrix[target.row.order, ]
		} else{ 
			sample.order <- order(target.ds$matrix[1,], decreasing = decreasing.order)
			target.ds$matrix <- t(as.matrix(target.ds$matrix[, sample.order]))
#			rownames(m2.pass1) = model.names
			target.row.order = 1 }
	}
	
	target.ds$sample.names = target.ds$sample.names[sample.order]
	target.ds$row.names = target.ds$row.names[target.row.order]
	target.ds$MI.vs.feature.summary = MI.targets.vs.feature.summary
	MI.targets.vs.feature.summary = MI.targets.vs.feature.summary[target.row.order]
	feature.summary = feature.summary[sample.order]
#	feature.ds
	
	feature.ds$sample.names = feature.ds$sample.names[sample.order] 
	print(matrix(c(target.ds$row.names, signif(MI.targets.vs.feature.summary, digits=3)), ncol=2), quote=F)
	if(show.multiple.tissues) tissue.labels = tissue.labels[sample.order]
	#### ---- End Normalized Mutual Information of targets vs feature.summary ----	
	
	
	MI.features.vs.best.target = vector(mode="numeric", length=feature.ds$nRow)
	#browser()
	MI.target.vs.target = 
			mutual.inf.one.vs.one(
#			mutual.inf.2(
#			mutual.inf.2.kde2d.olga(
					target.ds$matrix[1,], target.ds$matrix[1,])
	if( length(unique(feature.summary)) > 1){
		
		MI.feature.summary.vs.best.target = 
				mutual.inf.one.vs.one(
#				mutual.inf.2(
#				mutual.inf.2.kde2d.olga(		
						feature.summary, target.ds$matrix[1,])/MI.target.vs.target
	} else{
		MI.feature.summary.vs.best.target = NA
	}
	print(paste(format("feature summary", width=12), "mutual.inf =", MI.feature.summary.vs.best.target
			))
	print(proc.time()-t1)
	print(date())
	
	if( feature.ds$nRow > 1 ){
		feature.ds$matrix = feature.ds$matrix[,sample.order]
		MI.features.vs.best.target = 
				mutual.inf.one.vs.many(
#				mutual.inf.3.v2(
#				REVEALER.mutual.inf.one.vs.many(
						target.ds$matrix[1,],
						feature.ds$matrix)$MI
		feature.row.order = order(MI.features.vs.best.target, 
				decreasing=ifelse(if.feature.summary.on.top, 
						decreasing.order, !decreasing.order), na.last=TRUE)
		MI.features.vs.best.target = MI.features.vs.best.target[feature.row.order]
		feature.ds$matrix = feature.ds$matrix[feature.row.order,]
		feature.ds$row.names = feature.ds$row.names[feature.row.order]
		
	} else{
		feature.ds$matrix = t(as.matrix(feature.ds$matrix[,sample.order], nrow=1))
#		MI.ref = mutual.inf.2(target.ds$matrix[1,], target.ds$matrix[1,])
		MI.features.vs.best.target = 
				mutual.inf.one.vs.one(
#				mutual.inf.2(
#				mutual.inf.2.kde2d.olga(
						target.ds$matrix[1,],
						feature.ds$matrix[1,])/MI.target.vs.target
		feature.row.order = 1
		#non.unique.row.ind = which(sapply(apply(feature.ds$matrix, 1, unique), length) == 1)
	}
	#non.unique.row.ind = which(apply(apply(feature.ds$matrix, 1, unique), 2, length)==1)
	non.unique.row.ind = apply(feature.ds$matrix, 1, unique)
	if(is.list(non.unique.row.ind)){
		non.unique.row.ind = which(sapply(non.unique.row.ind, length)==1)
		if(length(non.unique.row.ind) > 0 ){
			MI.features.vs.best.target[non.unique.row.ind] = NA
		}
	} else if( is.matrix(non.unique.row.ind)){
		non.unique.row.ind = which(apply(non.unique.row.ind, 2, length)==1)
		if(length(non.unique.row.ind) > 0 ){
			MI.features.vs.best.target[non.unique.row.ind] = NA
		}
	}
	feature.ds$MI.features = MI.features.vs.best.target
	feature.ds$summary = feature.summary
	feature.ds$MI.summary = MI.feature.summary.vs.best.target
	print(matrix(c(feature.ds$row.names, 
							signif(MI.features.vs.best.target, digits=3)), ncol=2), quote=F)
	print(proc.time()-t1)
	print(date())
	
	filename <- paste(results.dir, test.file.prefix, file.suffix, ".",
			feature.ds$row.names[1],
#			".initial.calculations", 
			sep="")
#    pdf(file=paste(filename, ".pdf", sep=""), height = 11, width = 9)
	if(popup.heatmap){
		quartz(height = pdf.height, width = pdf.width)
	} else{
		pdf(file=paste(filename, ".pdf", sep=""), height = pdf.height, width = pdf.width )
	}
	#browser()
#   windows(width=12, height=8)
	if( show.multiple.tissues ){
		print("show.multiple.tissues doesn't work yet!")
#		MSIG.HeatMapPlot.10.show.multiple.tissues(V = target.ds$matrix, 
#				pathway.mut = bin.class.pass1,
#				row.names = model.names.pass1,
#				row.names2 = model.descs2.pass1, 
#				col.labels = cls.labels2.pass1, 
#				col.classes = cls.phen2.pass1, 
#				phen.cmap = colors.list.pass1, 
#				phen.names = phen.names.pass1,
#				phen.names2 = c(" ", phen.descs2.pass1),
#				col.names = sample.names2.pass1, 
#				main = paste(tissue, "- Step 1: Known KRAS Pathway Abnormalities (MI)"), 
#				xlab="  ", ylab="  ", row.norm = row.norm,  
#				cmap.type = cmap.type, char.rescale = char.rescale,  legend=F,
#				tissue.names = tissue.names,
#				tissue.labels = tissue.labels,
#				highlight.tissue.name = highlight.tissue.name,
#				ifScoreBarPlot = ifScoreBarPlot,
#				MI.list.model = MI.list.pass1,
#				MI.list.phen = MI.list.phen.pass1)
	} else{
		REVEALER.preprocessing.heatmap(
				target.ds,    ## list with: $matrix, $row.names, $nRow, $nSample
				feature.ds,   ## list with: $matrix, $row.names, $nRow, $nSample
				feature.summary,  ## Sing
				MI.targets.vs.feature.summary,   ## Raw values (not strings)
				MI.features.vs.best.target,  ## Raw values (not strings)
				MI.feature.summary.vs.best.target,
				if.feature.summary.on.top = if.feature.summary.on.top,
				main = filename,
				draw.white.lines = draw.white.lines)
	}
	dev.off()
	return(list(target.ds=target.ds, feature.ds=feature.ds))
}

REVEALER.preprocessing.heatmap <- function(
		target.ds,    ## list with: $matrix, $row.names, $nRow, $nSample
		feature.ds,   ## list with: $matrix, $row.names, $nRow, $nSample
		feature.summary,  ## vector
		MI.targets,   ## Raw values (not strings)
		MI.features,  ## Raw values (not strings)
		MI.feature.summary,  ## Raw values (not strings)
#		if.greyscale.feature.summary = FALSE,   ## greyscale can be helpful to show ranges, otherwise if 
#		## normalized feature summary > 0, values are interpreted as 1
#		## and plotted as black. if =0, then background (grey)
#		if.greyscale.feature.ds = FALSE,
		if.feature.summary.on.top = FALSE,  ## Put feature summary on top of rest of features or on bottom
		## Bottom is more like you're summing the features, and it then
		## aligns right next to the highest scoring target, which can be helpful
		main = " ", 
		sub = " ", 
		xlab=" ", 
		ylab=" ",
		row.norm = TRUE,
		char.rescale = 0.85,                               
		cmap.type = 4,   # 	1 = vintage pinkogram, 
		#					2 = scale of blues, 
		#					3 = high-resolution pinkogram for scores or probabilities [0, 1], 
		#					4 = high-resolution pinkogram for general values, 
		#					5 = color map for normalized enrichment scores, 
		#					6 = scale of red purples, 
		#					7 = scale of Oranges, 
		#					8 = scale of Greens, 
		#					9 = scale of Blues
		max.v = "NA",
		legend = T,
		tissue.names = "NA",
		tissue.labels = NA,
		ifScoreBarPlot = FALSE,
		row.label.text.resize = 1,
		col.label.text.resize = 1,
		draw.white.lines = FALSE
)
{
# Based off of MSIG.HeatMapPlot.10
# Plots a heatmap "pinkogram" of a gene expression matrix including phenotype vector and gene, sample and phenotype labels
# Doesn't plot the spectrum on the bottom
#
# Plots PATHWAY.MUT as a continuous vector in a greyscale spectrum
#
# The Broad Institute
# SOFTWARE COPYRIGHT NOTICE AGREEMENT
# This software and its documentation are copyright 2003 by the
# Broad Institute/Massachusetts Institute of Technology.
# All rights are reserved.
#
# This software is supplied without any warranty or guaranteed support
# whatsoever. Neither the Broad Institute nor MIT can be responsible for
# its use, misuse, or functionality.
	
	
#       if ((cmap.type == 5) | (cmap.type == 3)) {
	if (cmap.type == 5) {
		row.norm <- F
	}
	
	if(ifScoreBarPlot) char.rescale = 1.5*char.rescale
	
	if (row.norm == TRUE) {
		V1 <- matrix(0, nrow=target.ds$nRow, ncol=target.ds$nSample)
		row.mean <- apply(target.ds$matrix, MARGIN=1, FUN=mean)
		row.sd <- apply(target.ds$matrix, MARGIN=1, FUN=sd)
		row.n <- length(target.ds$matrix[,1])
		for (i in 1:target.ds$nRow) {
			if (row.sd[i] == 0) {
				V1[i,] <- 0
			} else {
				V1[i,] <- (target.ds$matrix[i,] - row.mean[i])/(0.333 * row.sd[i])
			}
			V1[i,] <- ifelse(V1[i,] < -4, -4, V1[i,])
			V1[i,] <- ifelse(V1[i,] > 4, 4, V1[i,])
		}
	} else {
		V1 <- target.ds$matrix
	}
	
	if (cmap.type == 1) { 
		mycol <- c("#0000FF", "#4040FF", "#7070FF", "#8888FF", "#A9A9FF", "#D5D5FF", "#EEE5EE", "#FFAADA",
				"#FF9DB0", "#FF7080", 
				"#FF5A5A", "#FF4040", "#FF0D1D") # blue-pinkogram colors. This is the 1998-vintage,
		# pre-gene cluster, original pinkogram color map
	} else if (cmap.type == 2) {
		violet.palette <- colorRampPalette(c("#400030", "white"), space = "rgb")
		mycol <- rev(violet.palette(20))
		
#          mycol <- c("#FCFBFD","#F4F2F8","#F8F7FB","#EFEDF5","#E1E1EF","#E8E7F2","#DADAEB","#C6C7E1","#D0D1E6",
#                        "#BCBDDC","#A8A6CF",
#                        "#B2B2D6","#9E9AC8","#8A87BF","#9491C4","#807DBA","#7260AB","#796FB3","#6A51A3","#5C3596",
#                        "#63439D","#54278F","#460D83","#4D1A89","#3F007D")
	} else if (cmap.type == 6) {
		mycol <- c("#FFF7F3", "#FDE0DD", "#FCC5C0", "#FA9FB5", "#F768A1", "#DD3497", "#AE017E",
				"#7A0177", "#49006A")
	} else if (cmap.type == 7) {
		mycol <- c("#FFF5EB", "#FEE6CE", "#FDD0A2", "#FDAE6B", "#FD8D3C", "#F16913", "#D94801",
				"#A63603", "#7F2704")
	} else if (cmap.type == 8) {
		mycol <- c("#F7FCF5", "#E5F5E0", "#C7E9C0", "#A1D99B", "#74C476", "#41AB5D", "#238B45",
				"#006D2C", "#00441B")
	} else if (cmap.type == 9) {
		mycol <- c("#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", "#4292C6", "#2171B5",
				"#08519C", "#08306B")
	} else if ((cmap.type == 3) | (cmap.type == 4) | (cmap.type == 5)) {
		mycol <- vector(length=512, mode = "numeric")
		
#		for (k in 1:256) {
#			mycol[k] <- rgb(255, k - 1, k - 1, maxColorValue=255)
#		}
#		for (k in 257:512) {
#			mycol[k] <- rgb(511 - (k - 1), 511 - (k - 1), 255, maxColorValue=255)
#		}
		mycol <- REVEALER.make.pinkogram.colors()
	}
#	browser()
	ncolors <- length(mycol)
	
#	image(1:n.cols, 1, as.matrix(feature.summary), col=gray(n.cols:0/n.cols))
	if (cmap.type == 5) {
		if (max.v == "NA") {
			max.v <- max(max(V1), -min(V1))
		}
		V2 <- ceiling(ncolors * (V1 - (- max.v))/(1.001*(max.v - (- max.v))))
		
	} else {
		V2 <- ceiling(ncolors * (V1 - min(V1))/(1.001*(max(V1) - min(V1))))
	}
	
	total.ncolors <- ncolors
	heatm.n.row = target.ds$nRow + feature.ds$nRow + 1
	heatm <- matrix(0, nrow = heatm.n.row, ncol = target.ds$nSample)
	heatm[1:target.ds$nRow,] <- V2[seq(target.ds$nRow, 1, -1),]
	#browser()
#			heatm[target.ds$nRow+feature.ds$nRow,] = t(as.matrix(gray(feature.summary)))
	feature.summary.heatmap = 
			REVEALER.preprocessing.heatmap.adjust.to.greyscale(feature.summary#, if.greyscale.feature.ds
			)
	feature.ds.heatmap = 
			REVEALER.preprocessing.heatmap.adjust.to.greyscale(feature.ds$matrix#, if.greyscale.feature.ds
			)
#	browser()
	if(!if.feature.summary.on.top){  # want feature summary on bottom of feature heatmap
		heatm[(target.ds$nRow+1):heatm.n.row,] = rbind(
				total.ncolors + feature.summary.heatmap$color.labels, 
				total.ncolors + feature.summary.heatmap$n.colors.unique + 
						feature.ds.heatmap$color.labels[seq(feature.ds$nRow, 1, -1),])
		mycol = c(mycol, feature.summary.heatmap$colors.unique, feature.ds.heatmap$colors.unique)
	} else{
		heatm[(target.ds$nRow+1):heatm.n.row,] = rbind(
				total.ncolors + feature.ds.heatmap$color.labels[seq(feature.ds$nRow, 1, -1),],
				total.ncolors + feature.ds.heatmap$n.colors.unique + feature.summary.heatmap$color.labels
		)
		mycol = c(mycol, feature.ds.heatmap$colors.unique, feature.summary.heatmap$colors.unique)
	}
	
	total.ncolors = length(mycol)
	if (legend == T && ifScoreBarPlot == FALSE
			) {
#		print("legend == T && ifScoreBarPlot == FALSE")
#              nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(10, 2), heights = c(6, 1), respect = FALSE)
		nf <- layout(matrix(c(1, 2, 3), 3, 1, byrow=T), heights = c(8, 4, 1), respect = FALSE)
	} else if( legend == FALSE && ifScoreBarPlot == FALSE 
			){
		nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(5, 1), respect = FALSE)
#		print("legend == FALSE && ifScoreBarPlot == FALSE")
#			nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), heights = c(8, 1), respect = FALSE)
	} else if( legend == TRUE && ifScoreBarPlot == TRUE){
#		print("legend == TRUE && ifScoreBarPlot == TRUE")
		#			nf <- layout(matrix(c(1, 2, 3, 0), 2, 2, byrow=T), widths = c(5, 1), heights = c(10, 1), respect = FALSE)
	} else if( legend == FALSE && ifScoreBarPlot == TRUE ){
#		print("legend == FALSE && ifScoreBarPlot == TRUE")
		nf <- layout(matrix(c(1, 1, 2, 1, 1, 2, 3, 3, 2), ncol=3, byrow=F), heights = c(5, 5, 2), 
				widths=c(5,5,1), respect = FALSE)
#			layout.show(3)
#			browser()
	}
	if(ifScoreBarPlot){
		par( mar = c(8, 20, 8, 0))
	} else{
		par( mar = c(3, 16, 3, 16))
	}
#		par(mar = c(8, 20, 8, ifelse(ifScoreBarPlot, 0, 20)))
#		browser()
#	mycol <- c(mycol, phen.cmap, u.feature.summary.grey)
#		browser()
	image(1:target.ds$nSample, 1:(heatm.n.row), 
			t(heatm), zlim = c(0, total.ncolors), col=mycol, axes=FALSE,
			main=main, sub = sub, xlab= xlab, ylab=ylab, 
			cex.main= 1 )#ifelse(ifScoreBarPlot, 2, 1))
	REVEALER.place.normalized.mutual.information.text(target.ds$nSample)
#	}
	
# Add lines to separate phenotypes or subgroups
	if( draw.white.lines){
		REVEALER.add.white.lines.to.heatmap(heatm.n.row, target.ds$nSample)
	}
#	for(x.ind in 1:(target.ds$nSample-1)){ # Make lines to delineate samples
#		lines(c(x.ind+0.5, x.ind+0.5),c(0,heatm.n.row+1), col="white", lwd=0.25)
#	}
#	for(y.ind in 1:(heatm.n.row-1)){ # Make lines to delineate rows
#		lines(c(0,target.ds$nSample+1), c(y.ind+0.5, y.ind+0.5), col="white", lwd=0.25)
#	}
	## Label rows and columns
	## FYI: axis(side, ...), where side = {1=below, 2=left, 3=above and 4=right}
	
	
	if (target.ds$row.names[1] != "NA") {
#		browser()
#		numC <- nchar(target.ds$row.names)
		size.row.char <- ifelse(ifScoreBarPlot, char.rescale*25/(target.ds$nRow + 25), 
				char.rescale*50/(target.ds$nRow + 25))
		heatm.row.names.left = vector(mode="character", length=target.ds$nRow)
		for (i in 1:target.ds$nRow) {
			heatm.row.names.left[i] <- substr(target.ds$row.names[i], 1, 40)
			heatm.row.names.left[i] <- paste(heatm.row.names.left[i], " ", sep="")
		}
		if (feature.ds$row.names[1] == "NA") {
			feature.names.left <- paste("Class", seq(feature.ds$nRow, 1, -1))
		} else {
			feature.names.left <- as.character(paste(feature.ds$row.names, ""))
#			feature.names.left <- as.character(paste(rownames(feature.ds$matrix), ""))
		}
		if(if.feature.summary.on.top){
			feature.names.left = c( "feature summary ", feature.names.left) 
		} else{
			feature.names.left = c(feature.names.left, "feature summary ")
		}
		
		heatm.row.names.left <- c(feature.names.left, heatm.row.names.left)
#		print("feature.ds:")
#		print(feature.ds)
#		print(paste("target.ds$nRow:", target.ds$nRow))
#		print(paste("feature.ds$row.names:", feature.ds$row.names))
#		print(paste("rownames(feature.ds$matrix):", rownames(feature.ds$matrix)))
#		print(paste("feature names.left:", feature.names.left))
#		print(paste("heatm.row.names.left:", heatm.row.names.left))
#		browser()
#		axis(2, at=1:heatm.n.row, labels=heatm.row.names.left, adj= 0.5, 
#				tick=FALSE, las = 1, cex.axis=size.row.char,
#				font.axis=2, line=-1)
	} else{ heatm.row.names.left = "NA"}
	
	if (!missing(MI.targets) && !missing(MI.features)) {  
		## For some reason this evaluates as true if both MI.targets 
		## and MI.features are provided
#		browser()
#		numC <- nchar(row.names2)
		MI.targets.signif = signif(MI.targets, digits=4)
		MI.features.signif = signif(MI.features, digits=4)
		MI.feature.summary.signif = paste("", signif(MI.feature.summary, digits=4), 
				"\n NMI of ",target.ds$row.names[1], 
				"\n with itself is 1", sep="")
		
		heatm.row.names.right = vector(mode="character", length=target.ds$nRow)
		feature.names.right = vector(mode="character", length=feature.ds$nRow)
#		size.row.char <- char.rescale*40/(target.ds$nRow + 20)
		
		for (i in 1:target.ds$nRow) {
			heatm.row.names.right[i] <- substr(MI.targets.signif[i], 1, 40)
			heatm.row.names.right[i] <- paste(" ", heatm.row.names.right[i], sep="")
		}
		heatm.row.names.right[1] = paste(heatm.row.names.right[1], "\n NMI of feature summary\n with itself is 1", sep="")
		for( i in 1:feature.ds$nRow ){
			feature.names.right[i] <- substr(MI.features.signif[i], 1, 40)
			feature.names.right[i] <- paste( " ", feature.names.right[i], sep="")
		}
		
		if(if.feature.summary.on.top){
			heatm.row.names.right = c(MI.feature.summary.signif,
					feature.names.right,
					heatm.row.names.right)
		} else{
			heatm.row.names.right = c(feature.names.right,
					MI.feature.summary.signif,
					heatm.row.names.right)
#			heatm.row.names.right = c(rev(heatm.row.names.right), 
#					MI.feature.summary.signif, 
#					rev(feature.names.right))
		}
		
#		heatm.row.names.right <- rev(heatm.row.names.right)
#		feature.names.right <- rev(feature.names.right)
#		axis(4, at=1:heatm.n.row, 
#				labels=heatm.row.names.right, 
#				adj= 0.5, tick=FALSE, las = 1, 
#				cex.axis=size.row.char, font.axis=2, line=-1)
	} else{ heatm.row.names.right = rep("", heatm.n.row)}
	
	if (target.ds$sample.names[1] != "NA" || feature.ds$sample.names[1] != "NA") {
#		size.col.char <- ifelse(ifScoreBarPlot, char.rescale*20/(target.ds$nSample + 15), 
#				char.rescale*100/(target.ds$nSample + 10
#							))
		if(target.ds$sample.names[1] != "NA"){
			heatm.bottom.labels = target.ds$sample.names
#			axis(1, at=1:target.ds$nSample, 
#					labels=paste(target.ds$sample.names, ""), 
#					tick=FALSE, 
#					las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
		} else{
			heatm.bottom.labels = feature.ds$sample.names
#			axis(1, at=1:feature.ds$nSample, 
#					labels=paste(target.ds$sample.names, ""),  
#					tick=FALSE, 
#					las = 3, cex.axis=size.col.char, font.axis=2, line=-1)
		}
	} else{ heatm.bottom.labels = "NA"}
	
	REVEALER.label.heatmap(
			heatmap.nrow = heatm.n.row, 
			heatmap.ncol = feature.ds$nSample, 
			left.labels = heatm.row.names.left, 
			right.labels = heatm.row.names.right, 
			bottom.labels = heatm.bottom.labels,
			row.label.text.resize = row.label.text.resize,
			col.label.text.resize = col.label.text.resize)
	
	
	# Phenotype Legend 
	
#      print("--------------------------------------------------------------------------------------------")
	if (legend == T) {
		legend.text <- NULL
		phenotype.vec <- NULL
		box.fill <- NULL
		box.border <- NULL
#		ind <- 1
		par(mar = c(0, 0, 0, 0))
		plot(c(0,0), c(1, 1), xlim = c(0, 1), ylim = c(-.050, 1.05), axes=F, type="n", xlab = "", ylab="")
		for (i in 1:feature.ds$nRow) { 
			this.feature.colors = unique(feature.ds.heatmap$raw.colors[i,])
			this.feature.num.colors = length(this.feature.colors)
#			browser()
#			if (is.vector(col.labels)) {
#			if (is.vector(feature.ds.heatmap$color.labels)){
#				phenotypes <- as.character(col.classes)
#			} else {
#				phenotypes <- as.character(col.classes[[i]])
#			}
			phenotypes = ifelse(this.feature.colors == "#000000", "MUT", "WT")
			phenotype.name <- paste(as.character(rev(feature.ds$row.names)[i]), ":   ", sep="")
			legend.text <- c(phenotype.name, phenotypes)  
			phenotype.vec <-  rep(22, this.feature.num.colors + 1)
			box.fill = c("#FFFFFF", this.feature.colors) #c("#FFFFFF", phen.cmap[ind:(ind + cols.row[i] - 1)])
			box.border <- c("#FFFFFF", rep("black", this.feature.num.colors))
#			ind <- ind + this.feature.num.colors
			offset <- 0.07
			legend(x=0, y= 0 + offset*i, 
					horiz = T, x.intersp = 0.5, legend=legend.text, bty="n", xjust=0, yjust= 1, pch = phenotype.vec, 
					pt.bg = box.fill, col = box.border, cex = 1.20, pt.cex=1.75)
		}
		i = i+1
		this.feature.colors = unique(as.vector(feature.summary.heatmap$raw.colors))
		this.feature.num.colors = length(this.feature.colors)
#			browser()
#			if (is.vector(col.labels)) {
#			if (is.vector(feature.ds.heatmap$color.labels)){
#				phen.v <- as.character(col.classes)
#			} else {
#				phen.v <- as.character(col.classes[[i]])
#			}
#		print("this.feature.colors")
#		print(this.feature.colors)
#		print(feature.summary.heatmap$raw.colors)
		phenotypes = ifelse(this.feature.colors == "#000000", "MUT", "WT")
		phenotype.name <- paste("feature summary :   ", sep="")
		legend.text <- c(phenotype.name, phenotypes)  
		phenotype.vec <-  rep(22, this.feature.num.colors + 1)
		box.fill = c("#FFFFFF", this.feature.colors) #c("#FFFFFF", phen.cmap[ind:(ind + cols.row[i] - 1)])
		box.border <- c("#FFFFFF", rep("black", this.feature.num.colors))
#			ind <- ind + this.feature.num.colors
		offset <- 0.07
		legend(x=0, y= 0 + offset*i, 
				horiz = T, x.intersp = 0.5, legend=legend.text, bty="n", xjust=0, 
				yjust= 1, pch = phenotype.vec, 
				pt.bg = box.fill, col = box.border, cex = 1.20, pt.cex=1.75)
	}
	
	
	
#	browser()
	
	# Color map legend
	
#       print(c("range V=", range(V)))
#       print(c("range V1=", range(V1)))
#       print(c("range V2=", range(V2)))
	
	#par(mar = c(2, 12, 6, 12))
	par(mar=rep(2,4))
	num.v <- 20
	range.v <- range(V2)
	incr <-  (range.v[1] - range.v[2])/(num.v - 1)
	heatm.v <- matrix(rev(seq(range.v[2], range.v[1], incr)), nrow=num.v, ncol=1)
	image(1:num.v, 1:1, heatm.v, zlim = c(0, total.ncolors), col=mycol, axes=FALSE,
			main=" ", sub = " ", xlab= ylab, ylab=xlab)
	range.v <- range(V1)
	incr <-  (range.v[1] - range.v[2])/(num.v - 1)
	heatm.v2 <- matrix(signif(rev(seq(range.v[2], range.v[1], incr)), digits=2), nrow=num.v, ncol=1)
	#          print(c("heatm.v2=", heatm.v2))
	axis(3, at=1:num.v, labels=heatm.v2, adj= 0.5, tick=FALSE, las = 1, cex.axis=1.5*char.rescale, font.axis=1)
	
#	browser()
	## Bar plots of MI values
	if(ifScoreBarPlot){
		par(mar = c(6, 5, 6, 1))
		barplot(rev(c(MI.list.phen, MI.list.model)), xlab="MI scores\nbarplot", ylab="", xlim=c(-1,1), 
				axes=TRUE, horiz=TRUE, axisnames=FALSE, 
#				main="MI scores\nbarplot", font.main=1
		)
#		browser()
#		text("MI scores\nbarplot")
	}
	return()
	
}

REVEALER.preprocessing.heatmap.adjust.to.greyscale = function(vector.or.matrix
#	ifGreyscale
){
#	browser()
	
	if(is.vector(vector.or.matrix) || ncol(vector.or.matrix) == 1){
		vector.or.matrix = t(as.matrix(REVEALER.normalize.zero.one(vector.or.matrix)))
	} else{
		vector.or.matrix = t(apply(vector.or.matrix, 1,
						REVEALER.normalize.zero.one))
	}
	
#	if(ifGreyscale){ 
#		vector.or.matrix = -(vector.or.matrix * .749 + 0.251 - 1)
#	} else{
	vector.or.matrix = 1 - ifelse(vector.or.matrix > 0, 1, 0.251)
#	}
	raw.colors = grey(vector.or.matrix)  # 0 = black, 1 = white
	colors.unique = unique(raw.colors)   # #000000 = black, #BFBFBF = 25.1% grey
	color.labels = match(raw.colors, colors.unique)
	n.colors.unique = length(colors.unique)
	
	if(is.matrix(vector.or.matrix)){
		color.labels = matrix(color.labels, nrow=nrow(vector.or.matrix))
		raw.colors = matrix(raw.colors, nrow=nrow(vector.or.matrix))
	}
	
#	browser()
	
	return(list(normalized.data = vector.or.matrix,
					color.labels = color.labels,
					raw.colors = raw.colors, 
					colors.unique = colors.unique , 
					n.colors.unique = n.colors.unique))
}

REVEALER.feature.selector = function(
		target.ds,
		feature.ds,
#		feature.ds.original = NULL,
		match.to.blue,
		decreasing.order,
		ifNaive = FALSE,
#		iterations.remaining = NULL,   ## Number of iterations left, for recursion. NULL if ifNaive==TRUE
		nTopFeatures = 20,
#		selected.features = NULL,   ## NULL if ifNaive==TRUE
		explained.samples = NULL,
		unexplained.samples.remaining = NULL,
		target.median,
		target.midpoint.index,
		iteration = NULL,
		results.dir = NULL,
		extracted.file.prefix = NULL,
		file.suffix = NULL,
		initial.feature.string = NULL,
		popup.heatmap = FALSE,
		row.label.text.resize = 0.75,
		col.label.text.resize = 0.5,
		tissue = "",
		kde2d.timing = FALSE,
		plot.unproductive.iteration = FALSE,
		draw.white.lines = FALSE
)
## The workhorse of REVEALER!
##
## Returns a list with:
## 	$top.hit.vector - vector with the same length as the target, feature with greatest mutual information with target
##	$MI.top.hit.vs.target - double, top-scoring MI
## 	$top.hit.name - string, name of top hit
##	$num.redundant.with.top - integer, number of features from dataset which had the exact same MI score,
##								meaning they had exactly the same aberration pattern, such as a genes from a
##								deleted chromosome
##	$explained - vector, 	0 = "unexplained" ie this sample doesn't have a feature that could 
##								"explain" its association with the target
##							1 = "explained," ie this sample has a feature that explains its
##								association with the target
## 	$unexplained - boolean vector, FALSE = "unexplained" **or** from non-matching side.
##							i.e. if "match.to.blue" then "unexplained" **or** on negative or "blue" side
{
#	browser()
	if(ifNaive){
		iteration = explained.samples = NULL
		unexplained.samples.remaining = rep(TRUE, target.ds$nSample)
	} else{
#		i = length(selected.features) - iterations.remaining + 1
#		print(paste("iteration:", i))
#		unexplained.samples.remaining = 
		print(paste("# samples left: ", sum(unexplained.samples.remaining), 
						"/", length(unexplained.samples.remaining), 
						"  (", length(unexplained.samples.remaining)-sum(unexplained.samples.remaining),
						" removed, and remember we only seek to 'explain' up to the median ('half') of the target",
						" or up to ", target.midpoint.index, " samples)",
						sep=""))
	}
	
	if(kde2d.timing){
		mutual.inf.one.vs.one = mutual.inf.2.kde2d.olga
		mutual.inf.one.vs.many = REVEALER.mutual.inf.one.vs.many
	} else{
		mutual.inf.one.vs.one = mutual.inf.2
		mutual.inf.one.vs.many = REVEALER.mutual.inf.one.vs.many.efficient
#		mutual.inf.one.vs.many = mutual.inf.3.v2
	}
	#browser()
#	top.genes.ind = NULL
#	top.genes.names = NULL
#	top.genes.vectors = NULL
#	top.genes.MI = NULL
#	top.diffs = NULL
#	explained.vectors = NULL
#	browser()
#	bin.gene.matrix.3 =  ifelse(cls.list2.pass3[-1,]=="WT", 0, 1)
#	m2.pass2.1.median = median(m2.pass2[1,])
	
	## target.midpoint.index and target.median are needed for plotting and
	## determining there are any samples remaning on the side of interest ("red" or "blue")
#	target.midpoint.index <- which.min(abs(target.ds$matrix[1,] - quantile(target.ds$matrix[1,], 0.5)))
#	target.median = median(target.ds$matrix[1,])
	
	## Shortened form of target for abbreviation purposes
#	target = target.ds$matrix[1,unexplained.samples.remaining]
	
#	grey.and.black = c("#C0C0C0", "#000000")
#	pathway.name = model.names.pass2[1]
	
	
	## Adjust "target" for removed samples
	target.unexplained = target.ds$matrix[1,unexplained.samples.remaining]
	
	
#	initial.gene.string = ifelse(de.novo.search, phen.pass1[1], phen.names.pass2[2])
	
#	 If we had naively searched the space without removing the explained cell lines
	mutual.inf.index <<- mutual.inf.index + 1
	#mutual.inf.t1[mutual.inf.index] <<- proc.time()[3]
	
	MI.features.vs.target = 
			mutual.inf.one.vs.many(
#			mutual.inf.3.v2(
#			REVEALER.mutual.inf.one.vs.many(
					target.unexplained, 
					feature.ds$matrix[,unexplained.samples.remaining])$MI
	#mutual.inf.t2[mutual.inf.index] <<- proc.time()[3]
	print(proc.time()-t1)
	print(date())
	
	
	feature.ds.MI.sorted = REVEALER.make.top.features.ds(feature.ds, 
			MI.features.vs.target,
			decreasing.order)
	
#	if(!popup.heatmap) dev.off()
	
#browser()
	
	MI.target.vs.target = 
			mutual.inf.one.vs.one(
#			mutual.inf.2(
#			mutual.inf.2.kde2d.olga(
					target.ds$matrix[1,], target.ds$matrix[1,])
	top.hit.ind.in.original = which(feature.ds$row.names ==
					feature.ds.MI.sorted$row.names[1])
	top.hit.vector = feature.ds$matrix[top.hit.ind.in.original,]
	top.hit.name = feature.ds$row.names[top.hit.ind.in.original]
	top.hit.MI = 
			mutual.inf.one.vs.one(
#			mutual.inf.2(
#			mutual.inf.2.kde2d.olga(		
					top.hit.vector,target.ds$matrix[1,])/MI.target.vs.target
	if(is.null(explained.samples)){
		new.explained = ifelse(top.hit.vector > 0, 1, 0)
		MI.difference = new.explained.MI = 
				mutual.inf.one.vs.one(
#				mutual.inf.2(
#				mutual.inf.2.kde2d.olga(		
						new.explained, target.ds$matrix[1,])/MI.target.vs.target
#		MI.difference = new.explained.MI - explained.MI
	} else{
		explained.MI = 
				mutual.inf.one.vs.one(
#				mutual.inf.2(
#				mutual.inf.2.kde2d.olga(		
						explained.samples, target.ds$matrix[1,])/MI.target.vs.target
		new.explained = ifelse(explained.samples + top.hit.vector > 0, 1, 0)
		new.explained.MI = 
				mutual.inf.one.vs.one(
#				mutual.inf.2(
#				mutual.inf.2.kde2d.olga(		
						new.explained, target.ds$matrix[1,])/MI.target.vs.target
		MI.difference = new.explained.MI - explained.MI
	}
#	if(!ifNaive) {browser()}
	if(!ifNaive && 
			!(ifelse(match.to.blue, MI.difference < 0, MI.difference > 0)) &&
			!plot.unproductive.iteration
#			ifelse(match.to.blue, MI.difference >= 0, MI.difference <= 0)
			){
		#browser()
		return("NULL")
	}
	
	write.table(matrix(c(c("Genomic Aberrations", feature.ds.MI.sorted$row.names), 
							c(paste("Normalized Mutual Information to", target.ds$row.names[1]), 
									feature.ds.MI.sorted$MI)),
					ncol=2), 
			file = paste(results.dir, extracted.file.prefix, file.suffix, ".REVEALER.feature.selector.",
					target.ds$row.names[1], ".",
					ifelse(ifNaive || is.null(iteration), "naive", paste("iter", iteration, sep="")),
					".txt", sep=""),
			quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE, append=FALSE)
	gc()
	
	
#	pdf(file=paste(results.dir, extracted.file.prefix, file.suffix, ".",
#					initial.gene.string,
#					".REVEALER.pdf", sep=""), 
#			height=8.5, width=11)
#	browser()
	if(popup.heatmap) quartz(height=8.5, width=11)
	#browser()
	REVEALER.feature.selector.heatmap(target.ds,
			target.midpoint.index,
			feature.ds.MI.sorted,
			match.to.blue,
			target.main.text = paste(tissue),
			feature.main.text = ifelse(ifNaive, 
					"Naive REVEALER: what features match best to this target?",
					paste("REVEALER top hits, iteration:", iteration)),
			nTopFeatures = nTopFeatures,
			row.label.text.resize = row.label.text.resize,
			col.label.text.resize = col.label.text.resize,
			unexplained.samples.remaining = unexplained.samples.remaining,
			draw.white.lines = draw.white.lines
#			iteration#,
#			results.dir,
#			extracted.file.prefix,
#			file.suffix,
#			initial.feature.string
	)
	
	#browser()
	unexplained.samples.remaining = 
			REVEALER.get.unexplained.samples.remaining(
					explained=new.explained, 
					target=target.ds$matrix[1,], 
					target.median=target.median, 
					match.to.blue=match.to.blue)
#	if(match.to.blue){
#		unexplained.samples.remaining = (explained == 0 | target.ds$matrix[1,] <= target.median)
#	} else{
#		unexplained.samples.remaining = (explained == 0 | target.ds$matrix[1,] >= target.median)
#	}
	#browser()
	return( list(vector = top.hit.vector,
					MI.vector = top.hit.MI,
					name = top.hit.name,
					num.redundant = feature.ds.MI.sorted$num.redundant.with.no.1.hit,
					MI.difference = MI.difference,
					explained = new.explained,
					MI.explained = new.explained.MI,
					unexplained = unexplained.samples.remaining,
					iteration = iteration))
	
#	if( !is.null(iterations.remaining) && !ifNaive ){
#		top.hit.ind.in.original = which(feature.ds.original$row.names ==
#						feature.ds.MI.sorted$row.names[1])
#		top.hit.vector = feature.ds.original$matrix[top.hit.ind.in.original,]
#		top.hit.name = feature.ds.original$row.names[top.hit.ind.in.original]
#		top.hit.MI = mutual.inf.2(top.hit.vector,target.ds$matrix[1,])/MI.target.vs.target
#		explained.MI = mutual.inf.2(explained, target.ds$matrix[1,])/MI.target.vs.target
#		new.explained = ifelse(explained + top.hit.vector > 0, 1, 0)
#		new.explained.MI = mutual.inf.2(new.explained, target.ds$matrix[1,])/MI.target.vs.target
#		MI.difference = new.explained.MI - explained.MI
#		### If MI is improving in the direction it should
#		if( (match.to.blue && MI.difference < 0) || ## Matching to "blue" (negative) side => MI should become more negative
#				(!match.to.blue && MI.difference > 0) ## Matching to "red" (positive) side => MI should become more positive
#				){ 
#			if(is.null(selected.features)){
#				selected.features = vector(mode="list", length=iterations.remaining)
#				selected.features[[1]]$original.explained = explained
#				selected.features[[1]]$original.explained.MI = explained.MI
#				selected.features[[1]]$vector = top.hit.vector
#				selected.features[[1]]$MI = top.hit.MI
#				selected.features[[1]]$name = top.hit.name
#				selected.features[[1]]$num.redundant = features.ds.MI.sorted$num.redundant.with.no.1.hit
#				selected.features[[1]]$MI.difference = top.hit.MI - explained.MI
#				selected.features[[1]]$explained = new.explained
#			} else{
	##				i = length(selected.features) - iterations.remaining + 1
#				selected.features[[i]]$vector = top.hit.vector
#				selected.features[[i]]$MI = top.hit.MI
#				selected.features[[i]]$name = top.hit.name
#				selected.features[[i]]$num.redundant = features.ds.MI.sorted$num.redundant.with.no.1.hit
#				selected.features[[i]]$MI.difference = top.hit.MI - explained.MI
#				selected.features[[i]]$explained = new.explained
#			}
#			
#			if( match.to.blue ){
#				unexplained.samples.remaining = (!new.explained | target.ds$matrix[1,] >= target.median) #m2.pass2.1.median)
#				if( (unexplained.samples.remaining == (target.ds$matrix[1,] >= target.median))){
#					continue.iterations = FALSE
#				} else continue.iterations = TRUE
	##			wo.mut.and.matchingSide = (samples.without.mut & m2.pass2[1,] < m2.pass2.1.median)
#			} else{
#				unexplained.samples.remaining = (!new.explained | target.ds$matrix[1,] <= target.median)
#				if( (unexplained.samples.remaining == (target.ds$matrix[1,] <= target.median))){
#					continue.iterations = FALSE
#				} else continue.iterations = TRUE
	##			wo.mut.and.matchingSide = (samples.without.mut & m2.pass2[1,] > m2.pass2.1.median)
#			}
#			
	##			if( continue.iterations )
	##				## Recursion! call REVEALER.feature.selector again
	##				Recall(
	##						target.ds,
	##						feature.ds,
	##						feature.ds.original,
	##						match.to.blue,
	##						decreasing.order,
	##						ifNaive = FALSE,
	##						iterations.remaining = iterations.remaining - 1,
	##						nTopFeatures = 20,
	##						selected.features = selected.features,   ## NULL if ifNaive==TRUE
	##						explained = new.explained,
	##						unexplained.samples.remaining = unexplained.samples.remaining)
#			
#		} #else{ }
#		return(selected.features)
#	} else{
#		return(feature.ds.MI.sorted)
#	}
	
#	while( ifelse(match.to.blue, -1, 1)*MI.diff.latest > 0  && n.iter <=5 
#			){
#		#for( i in 1:n.iter){
#		n.iter = n.iter +  1
#		i = n.iter
#		gc()
#		print(paste("iteration:", i))
#		print(paste("# samples left: ", sum(wo.mut.or.nonMatchingSide), "/", length(wo.mut.or.nonMatchingSide), 
#						"(", length(wo.mut.or.nonMatchingSide)-sum(wo.mut.or.nonMatchingSide),"removed)", sep=""))
#		MI.results = mutual.inf.3.v2( 
#				m2.pass2[1,wo.mut.or.nonMatchingSide], 
#				bin.gene.matrix.3[,wo.mut.or.nonMatchingSide] )
#		MI.order = order(MI.results$MI, decreasing=decreasing.order, na.last=TRUE)+1
#		
#		write.table(matrix(c(c("Genomic Aberrations", phen.names[MI.order]), 
#								c(paste("Normalized Mutual Information to", pathway.name), MI.results$MI[MI.order-1])),
#						ncol=2), 
#				file = paste(results.dir, extracted.file.prefix, file.suffix, ".", initial.gene.string,
#						".Step3_iter", i, ".txt", sep=""),
#				quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
#		top10.names = phen.names[MI.order[1:10]]
#		top10.MI = MI.results$MI[MI.order[1:10]-1]
#		top10.labels = bin.gene.matrix.3[MI.order[1:10]-1,wo.mut.or.nonMatchingSide] +1	
#		
#		top.genes.ind = c(top.genes.ind, MI.order[1] )
#		num.redundant = length(which(MI.results$MI == MI.results$MI[MI.order[1]-1]))-1
#		top.genes.names = c( top.genes.names, paste(phen.names[MI.order[1]], "+", num.redundant, 
#						ifelse(num.redundant==1, "other", "others")))
#		
#		mut = bin.gene.matrix.3[(MI.order[1]-1),]
#		explained = ifelse(mut+explained.prev>0, 1, 0)
#		explained.MI = mutual.inf.2( m2.pass2[1,wo.mut.or.nonMatchingSide], 
#				explained[wo.mut.or.nonMatchingSide])/MI.ref
#		MI.diff = MI.diff.latest = explained.MI - explained.MI.prev
#		print(paste("Explained.MI = ", explained.MI, 
#						"  MI.diff = ", ifelse(MI.diff<0, "-", "+"), 
#						signif(abs(MI.diff), digits=4), sep=""))
#		explained.vectors = rbind(explained.vectors, explained)
#		print(2*explained.prev + mut)
#		
#		top.diffs = c(top.diffs, MI.diff)
#		top.genes.vectors = rbind(top.genes.vectors, mut)
#		top.genes.MI = c(top.genes.MI, paste(signif(MI.results$MI[MI.order[1]-1], digits=4), 
#						sep="" ))
#		
#		
#		par(mar = c(2, 16, 2, 12))
#		nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(2, 10), FALSE)
#		max.v <- max(max(m2.pass2[1,wo.mut.or.nonMatchingSide]), -min(m2.pass2[1,wo.mut.or.nonMatchingSide]))
#		if( match.to.blue){
#			V1 <- c( (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,1:mid.point]), 
#					+ (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,(mid.point+1):Ns])+ncolors/2 )[wo.mut.or.nonMatchingSide]
#		} else{
#			V1 <- c( (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,1:mid.point]) + ncolors/2, 
#					(ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,(mid.point+1):Ns]))[wo.mut.or.nonMatchingSide]
#		}
	##			V1 <- c( (ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,1:mid.point]) + ncolors/2, 
	##					(ncolors/2)*REVEALER.normalize.zero.one(m2.pass2[1,(mid.point+1):Ns]))
#		
#		image(1:length(m2.pass2[1,wo.mut.or.nonMatchingSide]), 1:1, as.matrix(V1), 
#				zlim = c(0, ncolors), col=mycol, axes=FALSE, main="", sub = "", xlab= "", ylab="")
#		axis(2, at=1:1, labels=pathway.name, adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
#		V1 <- rbind( explained[wo.mut.or.nonMatchingSide]+1, top10.labels)
#		V1 <- apply(V1, MARGIN=2, FUN=rev)      
#		image(1:length(m2.pass2[1,wo.mut.or.nonMatchingSide]), 1:dim(V1)[1], t(V1), 
#				zlim = c(0, length(grey.and.black)), col=grey.and.black, 
#				axes=FALSE, main=paste("iteration:", i), sub = "", xlab= "", ylab="")
#		axis(2, at=1:dim(V1)[1], labels=rev(c("explained with top result", top10.names)), 
#				adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
#		axis(4, at=1:dim(V1)[1], labels=rev(c( paste(signif(explained.MI, digits=4), 
#										sep="" ), 
#								signif(top10.MI,digits=4))), 
#				adj= 0.5, tick=FALSE, las = 1, cex.axis=1, font.axis=1, line=-1)
#		
#		samples.without.mut[wo.mut.or.nonMatchingSide] = samples.without.mut[
#				wo.mut.or.nonMatchingSide] - mut[wo.mut.or.nonMatchingSide]
#		if( match.to.blue ){
#			wo.mut.or.nonMatchingSide = (samples.without.mut | m2.pass2[1,] >= m2.pass2.1.median)
#			wo.mut.and.matchingSide = (samples.without.mut & m2.pass2[1,] < m2.pass2.1.median)
#		} else{
#			wo.mut.or.nonMatchingSide = (samples.without.mut | m2.pass2[1,] <= m2.pass2.1.median)
#			wo.mut.and.matchingSide = (samples.without.mut & m2.pass2[1,] > m2.pass2.1.median)
#		}
#		explained.MI.prev = explained.MI
#		explained.prev = ifelse(mut+explained.prev>0, 1, 0)
#		print(proc.time()-t1)
#		print(date())
#		
#	}
	
}

REVEALER.feature.selector.heatmap = function(
		target.ds,
		target.midpoint.index,
		feature.ds,
		match.to.blue,
		ifNaive,
		iteration,
		row.label.text.resize = 0.75,
		col.label.text.resize = 0.75,
		target.main.text = "",
		feature.main.text = "",
		calculate.MI.diff = FALSE,
		nTopFeatures = NA,
		unexplained.samples.remaining = NA,
		initial.heatmap = FALSE,
		draw.white.lines = FALSE
#,
#		results.dir,
#		extracted.file.prefix,
#		file.suffix,
#		initial.feature.string
){
#	browser()
	target.main.text = paste("tissue:", target.main.text)
	mycol = REVEALER.make.pinkogram.colors()
	cex.axis = 1
	ncolors <- length(mycol)
	grey.and.black = c("#C0C0C0", "#000000")
	
	if(is.na(unexplained.samples.remaining[1])){
		plot.these.samples = 1:target.ds$nSample
		nSample = target.ds$nSample
	} else{
		plot.these.samples = unexplained.samples.remaining
		nSample = sum(plot.these.samples)
	}
	
	par(mar = c(0, 16, 2, 10))
	nf <- layout(matrix(c(1, 2), 2, 1, byrow=T), 1, c(2, 10), FALSE)
#	max.v <- max(max(target.ds$matrix[1,]), -min(target.ds$matrix[1,]))
	target.1st.half = target.ds$matrix[1,1:target.midpoint.index]
	target.2nd.half = target.ds$matrix[1,(target.midpoint.index+1):target.ds$nSample]
	target.max = max(target.ds$matrix[1,])
	target.min = min(target.ds$matrix[1,])
	
	## This crazy logic is to check if all of one side is equal to one value, the minimum or the maximum,
	## since the target.ds$matrix[1,] is (presumably) sorted. This avoids NaN values produced
	## from using REVEALER.normalize.zero.one on a vector of constants.
	#browser()
	if( match.to.blue){
		if( length(unique(target.1st.half)) == 1 || length(unique(target.2nd.half)) == 1){ 
			#unique(target.1st.half)[1] == target.min || unique(target.2nd.half)[1] == target.max){
			if( unique(target.2nd.half) == target.max ){
				target.heatmap <- c( (ncolors/2)*REVEALER.normalize.zero.one(target.1st.half), 
						(ncolors/2)*rep(1, target.ds$nSample-target.midpoint.index)+ncolors/2 )
			} else if (unique(target.1st.half) == target.min) {
				target.heatmap <- c( (ncolors/2)*rep(1/n.colors,target.midpoint.index), #rep(0,target.midpoint.index),  
						(ncolors/2)*REVEALER.normalize.zero.one(target.2nd.half)+ncolors/2 )
			}
		} else{
			target.heatmap <- c( (ncolors/2)*REVEALER.normalize.zero.one(target.1st.half), 
					(ncolors/2)*REVEALER.normalize.zero.one(target.2nd.half)+ncolors/2 )
		}
	} else{
		if( length(unique(target.1st.half)) == 1 || length(unique(target.2nd.half)) == 1){ 
			#(unique(target.1st.half)[1] == target.max || unique(target.2nd.half)[1] == target.min)){
			if( length(unique(target.2nd.half)) == 1 && unique(target.2nd.half) == target.min ){
				target.heatmap <- c( (ncolors/2)*REVEALER.normalize.zero.one(target.1st.half) + ncolors/2, 
						(ncolors/2)*rep(1/n.colors, target.ds$nSample-target.midpoint.index)) 
				# rep(0, target.ds$nSample-target.midpoint.index))
			} else if (length(unique(target.1st.half)) == 1 && unique(target.1st.half) == target.max) {
				target.heatmap <- c( (ncolors/2)*rep(1, target.midpoint.index) + ncolors/2, 
						(ncolors/2)*REVEALER.normalize.zero.one(target.2nd.half))
			}
		} else{
			target.heatmap <- c( (ncolors/2)*REVEALER.normalize.zero.one(target.1st.half) + ncolors/2, 
					(ncolors/2)*REVEALER.normalize.zero.one(target.2nd.half))
		}
	}
	names(target.heatmap) = target.ds$sample.names[plot.these.samples]
#	browser()
#	target.main.text = ifelse(ifNaive, 
#			"Naive REVEALER - what features match best to this target?",
#			paste("REVEALER top hits, iteration:", iteration))
	image(1:nSample, 1:1, as.matrix(target.heatmap[plot.these.samples]), 
			zlim = c(0, ncolors), col=mycol, axes=FALSE, 
			main=target.main.text, sub = "", xlab= "", ylab="")
	REVEALER.place.normalized.mutual.information.text(nSample)
#	max.col.label.nchar = max(unlist(sapply(target.ds$row.names[1], function(x) sapply(strsplit(x, "\n"), nchar))))
	max.col.label.nchar = nchar(target.ds$row.names[1])
#	browser()
	
	REVEALER.label.heatmap(
			heatmap.nrow = 1, 
			heatmap.ncol = nSample, 
			left.labels = paste("target:", ifelse(max.col.label.nchar > 15, " \n", ""), target.ds$row.names[1]), 
			right.labels = "1",
			row.label.text.resize = ifelse(max.col.label.nchar > 15, 
					(15/max.col.label.nchar)*row.label.text.resize, 
					row.label.text.resize))
	REVEALER.add.white.lines.to.heatmap(target.ds$nRow, nSample)
	par(mar = c(10, 16, 2, 10))
	if(initial.heatmap){
		feature.heatmap = rbind(feature.ds$matrix, feature.ds$summary)
		nTopFeatures = feature.ds$nRow + 1
	} else{
		if(is.na(nTopFeatures)){
			feature.heatmap <- feature.ds$matrix
			nTopFeatures = feature.ds$nRow
		} else{
			if( nTopFeatures > feature.ds$nRow ){ 
				nTopFeatures = feature.ds$nRow
			}
			feature.heatmap <- feature.ds$matrix[1:nTopFeatures,]
		}
	}
#	browser()
	if(is.vector(feature.heatmap)){
		feature.heatmap = t(as.matrix(feature.heatmap))
#		feature.heatmap
	}	
	feature.heatmap <- apply(feature.heatmap, MARGIN=2, FUN=rev)
	if(is.vector(feature.heatmap)){
		feature.heatmap = t(as.matrix(feature.heatmap))
#		feature.heatmap
	}
	if(nTopFeatures > 1){
		image(1:nSample, 1:nTopFeatures, 
				t(feature.heatmap[,plot.these.samples]), 
				zlim = c(0, 1), col=grey.and.black, axes=FALSE, 
				main=feature.main.text, #paste("step:", i), 
				sub = "", xlab= "", ylab="")
	} else{
		image(1:nSample, nTopFeatures, 
				as.matrix(feature.heatmap[,plot.these.samples]), 
				zlim = c(0, 1), col=grey.and.black, axes=FALSE, 
				main=feature.main.text, #paste("step:", i), 
				sub = "", xlab= "", ylab="")
	}
	if(calculate.MI.diff){
		MI.diff = feature.ds$MI[2:feature.ds$nRow] - feature.ds$MI[1:(feature.ds$nRow-1)]
		MI.diff.sign = ifelse(MI.diff >= 0, "+", "")
		right.labels = signif(feature.ds$MI, digits=4)
		right.labels[-1] = paste(right.labels[-1], " (", MI.diff.sign, 
				signif(MI.diff, 4), ")", sep="")
	} else if(initial.heatmap){
		right.labels = signif(c(feature.ds$MI.features, feature.ds$MI.summary), digits=4)
	} else{
		right.labels = signif(feature.ds$MI, digits=4)
	}
	if(initial.heatmap){
		left.labels = c(feature.ds$row.names, "explained original")
	} else{
		left.labels = feature.ds$row.names[1:ifelse(is.na(nTopFeatures), feature.ds$nRow, nTopFeatures)]
	}
	#browser()
#	if(nTopFeatures < 10){
	
#	max.col.label.nchar = max(sapply(left.labels, nchar))
	max.col.label.nchar = max(unlist(sapply(left.labels, function(x) sapply(strsplit(x, "\n"), nchar))))
#	browser()
	if( nTopFeatures < 5 && length(grep("seed:", left.labels)) > 0 && max.col.label.nchar > 25){
		left.labels = unlist(lapply(sapply(left.labels, strsplit, ": "), paste, collapse=": \n"))
	}
	
#		if(max.col.label.nchar > 25){
#			row.label.text.resize = 0.5*row.label.text.resize
#		}
#	}
#	browser()
	REVEALER.label.heatmap(
			heatmap.nrow = ifelse(is.na(nTopFeatures), feature.ds$nRow, nTopFeatures),
			heatmap.ncol = nSample, 
			left.labels = left.labels,
			right.labels = right.labels[1:ifelse(is.na(nTopFeatures), feature.ds$nRow, nTopFeatures)],
			#signif(feature.ds$MI, digits=4), 
			bottom.labels = target.ds$sample.names[plot.these.samples],
			row.label.text.resize = ifelse(nTopFeatures < 5 &&
							left.labels[1] !="explained original" &&
							max.col.label.nchar > 25, 
					(25/max.col.label.nchar)*row.label.text.resize, 
					row.label.text.resize),
			col.label.text.resize = col.label.text.resize)
	if(draw.white.lines){
		REVEALER.add.white.lines.to.heatmap(feature.ds$nRow, nSample)
	}
}

REVEALER.make.top.features.ds = function(feature.ds, 
		MI.features.vs.target,
		decreasing.order){
	feature.ds.MI.sorted = feature.ds
	feature.order = order(MI.features.vs.target, decreasing=decreasing.order, na.last=TRUE)
	
	feature.ds.MI.sorted$MI = MI.features.vs.target[feature.order]
	feature.ds.MI.sorted$matrix = feature.ds.MI.sorted$matrix[feature.order,]
	feature.ds.MI.sorted$row.names = feature.ds.MI.sorted$row.names[feature.order]
	feature.ds.MI.sorted$num.redundant.with.no.1.hit = 
			length(which(feature.ds.MI.sorted$MI[-1] == feature.ds.MI.sorted$MI[1]))
	
	if( feature.ds.MI.sorted$num.redundant.with.no.1.hit > 0){
#		browser()
		feature.ds.MI.sorted$row.names[1:(feature.ds.MI.sorted$num.redundant.with.no.1.hit+1)] = 
				sort(feature.ds.MI.sorted$row.names[1:(feature.ds.MI.sorted$num.redundant.with.no.1.hit+1)])
	}
#	top.names = c( 
#			feature.ds$row.names[feature.order[1:nTopFeatures]])
#	top.MI = c( 
#			signif(MI.features.vs.target[feature.order[1:nTopFeatures]-1], digits=4))
#	top.labels = rbind( 
#			bin.gene.matrix.3[feature.order[1:nTopFeatures]-1,]+1)
#	num.redundant = length(which(MI.features.vs.target == feature.ds.MI.sorted$MI[1]))
	
	
#	write.table(matrix(c(c("Genomic Aberrations", feature.ds.MI.sorted$row.names), 
#							c(paste("Normalized Mutual Information to", target.ds$row.names[1]), 
#									feature.ds.MI.sorted$MI)),
#					ncol=2), 
#			file = paste(results.dir, test.file.prefix, file.suffix, ".Step3_", 
#					ifelse(ifNaive || is.null(nIter), "naive", paste("iter", nIter, sep="")),
#					".txt", sep=""),
#			quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE, append=FALSE)
	
	return(feature.ds.MI.sorted)
}



REVEALER.get.unexplained.samples.remaining = function(explained, target, target.median, match.to.blue){
	## Returns a boolean vector of which samples are "unexplained", that is,
	## do not have a nonzero value for explained, or are on the opposite "matching"
	## side, ie if you're matching to high numbers => you want to explain the "red"
	## => match.to.blue = FALSE => want either explained = 0 or "blue" values to 
	## make sure you don't get any features that have high mutual information with
	## your opposite side of interest.
	if(match.to.blue){
		unexplained.samples = (explained == 0 | target > target.median)
	} else{
		unexplained.samples = (explained == 0 | target < target.median)
	}
	return(unexplained.samples)
}

REVEALER.find.median.index = function(v, match.to.blue){
	## This function assumes that v is sorted,
	## if match.to.blue:
	##		lowest to highest
	## if !match.to.blue:
	##		highest to lowest
	v.minus.median = abs(v - quantile(v, 0.5))
	min.indices = which(v.minus.median == min(v.minus.median))
	return(ifelse(match.to.blue, min(min.indices), max(min.indices)))
}

REVEALER.add.white.lines.to.heatmap = function(heatmap.nrow, heatmap.ncol){
	if(heatmap.ncol > 300){
		return()
#		line.width = 0.01
	} else if(heatmap.ncol > 200) {
		line.width = 0.05
	} else if (heatmap.ncol > 100){
		line.width = 0.1
	} else { line.width = 0.25}
	
	for(x.ind in 1:(heatmap.ncol-1)){ # Make lines to delineate samples
		lines(c(x.ind+0.5, x.ind+0.5),c(0,heatmap.nrow+1), col="white", lwd=line.width)
	}
	for(y.ind in 1:(heatmap.nrow-1)){ # Make lines to delineate rows
		lines(c(0,heatmap.ncol+1), c(y.ind+0.5, y.ind+0.5), col="white", lwd=line.width)
	}
}

REVEALER.place.normalized.mutual.information.text = function(heatmap.ncol){
	nmi.text = "Normalized\nMutual\nInformation"
#	browser()
	text.position = heatmap.ncol+1#heatmap.ncol/20
	mtext(nmi.text, side=3, at=text.position, #target.ds$nSample+6.5
			padj=.5, adj=0)
}

REVEALER.label.heatmap = function(
		heatmap.nrow, 
		heatmap.ncol, 
		left.labels = NULL, 
		right.labels = NULL, 
		bottom.labels = NULL,
		top.labels = NULL,
		row.label.text.resize = 1,
		col.label.text.resize = 1){
	
	size.row.char <- #ifelse(ifScoreBarPlot, char.rescale*20/(target.ds$nSample + 15), 
			row.label.text.resize*20/(heatmap.nrow + 10)#)
	if(!is.null(left.labels)){
		axis(2, at=1:heatmap.nrow, labels=rev(paste(left.labels, "")), 
				adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char, 
				font.axis=1, line=-1)
	}
	if(!is.null(right.labels)){
		axis(4, at=1:heatmap.nrow, labels=rev(paste("", right.labels)), 
				adj= 0.5, tick=FALSE, las = 1, cex.axis=size.row.char, 
				font.axis=1, line=-1)
	}
	
	size.col.char <- #ifelse(ifScoreBarPlot, char.rescale*20/(target.ds$nSample + 15), 
			col.label.text.resize*50/(ifelse(heatmap.ncol < 50, 50, heatmap.ncol) + 
				ifelse(heatmap.ncol > 200, 15, 10))#)
	if(!is.null(bottom.labels)){
		axis(1, at=1:heatmap.ncol, 
				labels=paste(bottom.labels, ""), 
				tick=FALSE, 
				las = 3, cex.axis=size.col.char, font.axis=1, line=-1)
	}
	if(!is.null(top.labels)){
		axis(1, at=1:heatmap.ncol, 
				labels=paste(top.labels, ""), 
				tick=FALSE, 
				las = 3, cex.axis=size.col.char, font.axis=1, line=-1)
	}
}

REVEALER.plot.final.results = function(
		selected.feature, 
		seed.ds,   ## dataset of the features used to seed the feature selection search
		target.ds, 
		explained.original,
		match.to.blue,
		target.midpoint.index,
		n.productive.iter,
		popup.heatmap,
		row.label.text.resize = 0.75,
		col.label.text.resize = 0.5,
		tissue,
		draw.white.lines = FALSE){
#	browser()
	n.feature = n.productive.iter #length(selected.feature)
#	browser()
	
	
	feature.ds = explained.ds = vector(length=length(target.ds), mode="list")
	names(feature.ds) = names(explained.ds) = names(target.ds)
	
	#browser()
	feature.ds$matrix = seed.ds$matrix
	feature.ds$row.names = paste("seed:", seed.ds$row.names)
	feature.ds$nRow = seed.ds$nRow
	feature.ds$sample.names = explained.ds$sample.names = seed.ds$sample.names
	feature.ds$nSample = explained.ds$nSample = seed.ds$nSample
	feature.ds$MI = seed.ds$MI.features
	if(n.feature > 0){
		feature.ds$matrix = rbind(feature.ds$matrix, 
				t(sapply(selected.feature[1:n.feature], function(x) x$vector)))#,
#			selected.feature[[n.feature]]$explained
		
		feature.ds$row.names = c(feature.ds$row.names, 
				sapply(selected.feature[1:n.feature], function(x) 
							paste("iteration ", x$iteration, ": ", x$name, " ",
									paste(ifelse(x$num.redundant > 0, 
													paste("\n(+", x$num.redundant, 
															" other",
															ifelse(x$num.redundant > 1, "s", ""),
															" with the same pattern)", sep=""), "")),
									sep="")))
		feature.ds$nRow = feature.ds$nRow + n.feature
		feature.ds$MI = c(feature.ds$MI, 
				sapply(selected.feature[1:n.feature], function(x) x$MI.vector))#,
	} 
#else{
#	feature.ds$matrix = feature.ds$matrix#,
	##			selected.feature[[n.feature]]$explained
#	feature.ds$row.names = feature.ds$row.names
#	feature.ds$nRow = feature.ds$nRow
#	feature.ds$MI = feature.ds$MI
#}
#			selected.feature[[n.feature]]$MI.explained)
#	feature.ds$row.names = c(feature.ds$row.names, "explained")
#	feature.ds$nRow = feature.ds$nRow + 1
	
	
	if(popup.heatmap) quartz(height=8.5, width=11)
	REVEALER.feature.selector.heatmap(target.ds,
			target.midpoint.index,
			feature.ds,
			match.to.blue,
			target.main.text = paste(tissue),
			feature.main.text =  "REVEALER results: top hit from each iteration",
			row.label.text.resize = row.label.text.resize,
			col.label.text.resize = col.label.text.resize,
			draw.white.lines = draw.white.lines#,
#			results.dir,
#			extracted.file.prefix,
#			file.suffix,
#			initial.feature.string
	)
	
	explained.ds$matrix = explained.original$explained
	explained.ds$MI = explained.original$MI
	explained.ds$nRow = 1
	explained.ds$row.names = "explained original"
	if(n.feature > 0){
		explained.ds$matrix = rbind(explained.ds$matrix, 
				t(sapply(selected.feature[1:n.feature], function(x) x$explained)))
		explained.ds$row.names = c(explained.ds$row.names,
				sapply(selected.feature[1:n.feature], function(x) 
							paste("explained after iteration", x$iteration)))
		explained.ds$nRow = explained.ds$nRow + ifelse(n.feature > 0, n.feature, 0)
		explained.ds$MI = c(explained.ds$MI, 
				sapply(selected.feature[1:n.feature], function(x) x$MI.explained))
	} 
#else{
#	explained.ds$matrix = explained.ds$matrix
#	explained.ds$row.names = explained.ds$row.names
#	explained.ds$nRow = explained.ds$nRow
#	explained.ds$MI = explained.ds$MI
#}
	
	
	if(popup.heatmap) quartz(height=8.5, width=11)
	REVEALER.feature.selector.heatmap(target.ds,
			target.midpoint.index,
			explained.ds,
			match.to.blue,
			target.main.text = paste(tissue),
			feature.main.text = paste("REVEALER results: explained after each iteration"),
			row.label.text.resize = row.label.text.resize,
			col.label.text.resize = col.label.text.resize,
			calculate.MI.diff = TRUE,
			draw.white.lines = draw.white.lines#,
#			results.dir,
#			extracted.file.prefix,
#			file.suffix,
#			initial.feature.string
	)
}

kde2d.olga = function (x, y, h, n = 25, lims = c(range(x), range(y))) 
{
	nx <- length(x)
	if (length(y) != nx) 
		stop("data vectors must be the same length")
	if (any(!is.finite(x)) || any(!is.finite(y))) 
		stop("missing or infinite values in the data are not allowed")
	if (any(!is.finite(lims))) 
		stop("only finite values are allowed in 'lims'")
	n <- rep(n, length.out = 2L)
	gx <- seq.int(lims[1L], lims[2L], length.out = n[1L])
	gy <- seq.int(lims[3L], lims[4L], length.out = n[2L])
	h <- if (missing(h)) 
				c(bandwidth.nrd(x), bandwidth.nrd(y))
			else rep(h, length.out = 2L)
	h <- h/4
	#browser()
#	kde2d.t0[kde2d.i] <<- proc.time()[3]
	ax <- outer(gx, x, "-")/h[1L]
#	kde2d.t1[kde2d.i] <<- proc.time()[3]
	ay <- outer(gy, y, "-")/h[2L]
#	kde2d.t2[kde2d.i] <<- proc.time()[3]
	z <- tcrossprod(matrix(dnorm(ax), , nx), matrix(dnorm(ay), 
					, nx))/(nx * h[1L] * h[2L])
#	kde2d.t3[kde2d.i] <<- proc.time()[3]
	list(x = gx, y = gy, z = z)
}


REVEALER.mutual.inf.y.bandwidth <- function(x, y,
#		x.bandwidth = NULL,
		y.bandwidth = NULL, 
		n.grid=20#, 
#		normalize.by ="HXY", # Whether to normalize by HXY, HX, or HY
#		pos.and.neg = T


) {
	# x and y can be binary or continuous
	# If there is not sufficient variation in x and y, 
	# will take the standard deviation as the bandwidth
	# (IQR finds the inter-quartile range of the data vector)
	
	## Add random noise if vectors are constant to prevent errors in bcv 
	## (because the inter-quartile range of a constant vector is 0, and the IQR is 
	## used in calculating the bandwidth)
#	y = unlist(y); y = as.vector(y)
#	x = unlist(x); x = as.vector(x)
	if( length(unique(y)) == 1 || length(unique(x)) == 1){ return(NA) }
	if( sd(y) == 0){
		y = y + runif(n=length(y), min=mean(y)-0.001, max=mean(y)+0.001)
	}
	if( sd(x) == 0){
		x = x + runif(n=length(x), min=mean(x)-0.001, max=mean(x)+0.001)
	}
	
	if(is.null(y.bandwidth)){
		y.bandwidth = suppressWarnings(bcv(y))
	}
	
	## Using suppressWarnings because bcv(.) gets mad if the minimum is on one side of the vector
#	before.kde2d[kde2d.i] <<- proc.time()[3]
	kde2d.xy <- kde2d(x, y, n = n.grid, h = c(suppressWarnings(bcv(x)), 
					#					suppressWarnings(bcv(y))
					y.bandwidth
			) )
#	kde2d.xy <- kde2d.olga(x, y, n = n.grid, h = c(suppressWarnings(bcv(x)), 
	##					suppressWarnings(bcv(y))
#					y.bandwidth
#			) )
#	after.kde2d[kde2d.i] <<- proc.time()[3]
#	kde2d.i <<- kde2d.i + 1
	PXY <- kde2d.xy$z/sum(kde2d.xy$z)
	
	PX <- rowSums(PXY)#apply(PXY, MARGIN=1, sum)
	PX <- PX/sum(PX)
	HX = -sum(PX * log2(PX))
	PX <- matrix(PX, nrow=n.grid, ncol=n.grid)
	
	PY <- colSums(PXY)#apply(PXY, MARGIN=2, sum)
	PY <- PY/sum(PY)
	HY = -sum( PY * log2(PY))
	PY <- matrix(PY, byrow = TRUE, nrow=n.grid, ncol=n.grid)
	
	HXY <- - sum(PXY * log2(PXY), na.rm=TRUE)
	MI.norm =  sign(cor(x, y))* ((HX + HY)/HXY - 1) #MI/normalization.factor
	
	return( MI.norm )
}


REVEALER.mutual.inf.one.vs.many <- function(one.vector, 
		many.matrix, 
		n.grid=20, 
		target.vector.name = "",
		tissue = "NA", 
#		normalize.by = "HXY", 
#		pos.and.neg=T, 
#		print.MI.ref=TRUE,
		normalize.by.MI.ref = TRUE,
		n.many) {
	## How this is different from mutual.inf.2:
	## x is a target pathway and y is a matrix
	## calculates a false discovery rate for the mutual information 
	## of the pathway with a randomized version of x
	## y.ind indicates the row index of the target gene set / signature
	
	# x and the elements of y can be binary or continuous
	# If there is not sufficient variation in x and y, 
	# will take the standard deviation as the bandwidth
	# (IQR finds the inter-quartile range of the data vector)
	
	# "signature.indices" is the indices of signatures that you are interested in using.
	# This code is used in comparing the chosen signatures to SUMMARY, the mutation
	# status of all the cell lines. 
	
#	n.signatures = length(comparison.matrix[,1])
#	MI.vector = vector(length=n.many, mode="double")
	
	one.bandwidth = bcv(one.vector)
	
	MI.ref = REVEALER.mutual.inf.x.bandwidth(one.vector, one.vector, 
			y.bandwidth = one.bandwidth
	)
	
#	if(print.MI.ref) print(paste("MI.ref =", MI.ref))
	MI.vector = apply(many.matrix, 
			MARGIN=1, 
			FUN=REVEALER.mutual.inf.x.bandwidth, 
			one.vector,
			y.bandwidth = one.bandwidth
	)
	MI = MI.vector/ifelse(normalize.by.MI.ref, MI.ref, 1)
	FDR = rep(1, length=length(MI))
	
	return(list(MI=MI, FDR=FDR))
}

REVEALER.mutual.inf.one.vs.many.efficient <- function(one, # a vector 
		many, # a matrix 
		n.grid=20#, 
#		target.vector.name = "",
#		tissue = "NA", 
##		normalize.by = "HXY", 
##		pos.and.neg=T, 
##		print.MI.ref=TRUE,
#		normalize.by.MI.ref = TRUE,
#		n.many
) {
	## How this is different from mutual.inf.2:
	## x is a target pathway and y is a matrix
	## calculates a false discovery rate for the mutual information 
	## of the pathway with a randomized version of x
	## y.ind indicates the row index of the target gene set / signature
	
	# x and the elements of y can be binary or continuous
	# If there is not sufficient variation in x and y, 
	# will take the standard deviation as the bandwidth
	# (IQR finds the inter-quartile range of the data vector)
	
	# "signature.indices" is the indices of signatures that you are interested in using.
	# This code is used in comparing the chosen signatures to SUMMARY, the mutation
	# status of all the cell lines. 
	
#	n.signatures = length(comparison.matrix[,1])
#	MI.vector = vector(length=n.many, mode="double")
	
	one.bandwidth = bcv(one)
	
	if( length(unique(one)) == 1){
		return(rep(NA, length=length(many[,1])))
	}
	
	if(sd(one) == 0){
		one = one + runif(n=length(one),
				min=mean(one)-0.001, max=mean(one)+0.001)
	}
	
	normalized.mutual.information = apply(
			many, 1,
			function(x){
				if( length(unique(x)) == 1){
					return(NA)
				}
				
				if(sd(x) == 0){
					x = x + runif(n=length(x),
							min=mean(x)-0.001, max=mean(x)+0.001)
				}
				
				kde2d.xy = kde2d(one, x, n = n.grid,
						h=c(one.bandwidth, suppressWarnings(bcv(x))))
				PXY = kde2d.xy$z/sum(kde2d.xy$z)
				
				PX = rowSums(PXY)
				PX = PX/sum(PX)
				HX = -sum(PX * log2(PX))
				PX = matrix(PX, nrow=n.grid, ncol=n.grid)
				
				PY = colSums(PXY)
				PY = PY/sum(PY)
				HY = -sum(PY * log2(PY))
				PY = matrix(PY, byrow=TRUE, nrow=n.grid, ncol=n.grid)
				
				HXY = -sum(PXY * log2(PXY), na.rm=TRUE)
				return( sign(cor(one, x)) * ((HX + HY)/HXY - 1))
			}
	)
	return(list(MI=normalized.mutual.information, FDR=rep(1,length(many[,1]))))
}

REVEALER.mutual.inf.one.vs.many.efficient.fdr <- function(one, # a vector 
		many, # a matrix 
		n.grid=20,
		n.random.perm = 1000
## With p-value calculation!
) {
	## How this is different from mutual.inf.2:
	## x is a target pathway and y is a matrix
	## calculates a false discovery rate for the mutual information 
	## of the pathway with a randomized version of x
	## y.ind indicates the row index of the target gene set / signature
	
	# x and the elements of y can be binary or continuous
	# If there is not sufficient variation in x and y, 
	# will take the standard deviation as the bandwidth
	# (IQR finds the inter-quartile range of the data vector)
	
	# "signature.indices" is the indices of signatures that you are interested in using.
	# This code is used in comparing the chosen signatures to SUMMARY, the mutation
	# status of all the cell lines. 
	
#	n.signatures = length(comparison.matrix[,1])
#	MI.vector = vector(length=n.many, mode="double")
	if( length(unique(one)) == 1){
		return(rep(NA, length=length(many[,1])))
	}
	one.bandwidth = bcv(one)
	
	if(sd(one) == 0){
		one = one + runif(n=length(one),
				min=mean(one)-0.001, max=mean(one)+0.001)
	}
	
	normalized.mutual.information = apply(
			many, 1,
			function(x){
				if( length(unique(x)) == 1){
					return(NA)
				}
				
				if(sd(x) == 0){
					x = x + runif(n=length(x),
							min=mean(x)-0.001, max=mean(x)+0.001)
				}
				
				kde2d.xy = kde2d(one, x, n = n.grid,
						h=c(one.bandwidth, suppressWarnings(bcv(x))))
				PXY = kde2d.xy$z/sum(kde2d.xy$z)
				
				PX = rowSums(PXY)
				PX = PX/sum(PX)
				HX = -sum(PX * log2(PX))
				PX = matrix(PX, nrow=n.grid, ncol=n.grid)
				
				PY = colSums(PXY)
				PY = PY/sum(PY)
				HY = -sum(PY * log2(PY))
				PY = matrix(PY, byrow=TRUE, nrow=n.grid, ncol=n.grid)
				
				HXY = -sum(PXY * log2(PXY), na.rm=TRUE)
				return( sign(cor(one, x)) * ((HX + HY)/HXY - 1))
			}
	)
	
	false.discovery.rate = sapply(
			1:n.random.perm, 
			function(a){
				one.random = sample(one)
				one.random.bandwidth = bcv(one.random.bandwidth)
				
				apply(many, 1, 
						function(x){
							if(sd(x) == 0){
								x = x + runif(n=length(x),
										min=mean(x)-0.001, max=mean(x)+0.001)
							}
							
							kde2d.xy = kde2d(one.random, x, n = n.grid,
									h=c(one.random.bandwidth, suppressWarnings(bcv(x))))
							PXY = kde2d.xy$z/sum(kde2d.xy$z)
							
							PX = rowSums(PXY)
							PX = PX/sum(PX)
							HX = -sum(PX * log2(PX))
							PX = matrix(PX, nrow=n.grid, ncol=n.grid)
							
							PY = colSums(PXY)
							PY = PY/sum(PY)
							HY = -sum(PY * log2(PY))
							PY = matrix(PY, byrow=TRUE, nrow=n.grid, ncol=n.grid)
							
							HXY = -sum(PXY * log2(PXY), na.rm=TRUE)
							return( sign(cor(one.random, x)) * ((HX + HY)/HXY - 1))							
						}
				)
			}
	)
	return(list(MI=normalized.mutual.information, FDR=false.discovery.rate))
}
