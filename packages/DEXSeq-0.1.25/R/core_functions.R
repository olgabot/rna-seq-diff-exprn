setMethod("estimateSizeFactors", signature(object="ExonCountSet"),
   function( object ){
      cds <- object
      stopifnot( is( cds, "ExonCountSet") )
      geomeans <- exp( rowMeans( log( counts(cds) ) ) )
      sizeFactors(cds) <- 
         apply( counts(cds), 2, function(cnts) 
            median( ( cnts / geomeans )[ geomeans>0 ] ) )
      cds
   }
)

countTableForGene <- function( ecs, geneID, normalized=FALSE, withDispersion=FALSE ) {
   stopifnot( is( ecs, "ExonCountSet") )
   rows <- geneIDs(ecs) == geneID
   if( sum(rows) == 0 )
      stop( "geneID not found" )
	ans <- counts( ecs )[ rows, , drop=FALSE ]
	rownames( ans ) <- exonIDs(ecs)[ rows ]
	attr( ans, "geneID" ) <- geneID
	if( normalized )
	   if(any(is.na(sizeFactors(ecs)))){
   	   stop("Please first call function estimateSizeFactors\n")
	   }else{
	   ans <- t( t(ans) / sizeFactors(ecs) )
	   }
	if( withDispersion )
      attr( ans, "dispersion" ) <- fData(ecs)$dispersion[rows]
   ans
}

modelFrameForGene <- function( ecs, geneID, onlyTestable=FALSE) {
   stopifnot( is( ecs, "ExonCountSet") )
   if( onlyTestable & any(colnames(fData(ecs)) %in% "testable")){
      rows <- geneIDs(ecs) == geneID & fData(ecs)$testable
   }else{
      rows <- geneIDs(ecs) == geneID
   }
   
   numExons <- sum(rows)
   exonCol <- rep( factor( exonIDs(ecs)[rows], levels=exonIDs(ecs)[rows] ), ncol( counts(ecs) ) )
   modelFrame <- data.frame(
      sample = rep( factor( colnames( counts(ecs) ) ), each = numExons ),
      exon = exonCol,
      sizeFactor = rep( sizeFactors(ecs), each = numExons ) )
   for( cn in colnames( design(ecs,drop=FALSE) ) )
      modelFrame[[cn]] <- factor(rep( design(ecs,drop=FALSE)[[cn]], each=numExons ))
   modelFrame$dispersion <- fData(ecs)$dispersion[ rows ][ 
      match( modelFrame$exon, exonIDs(ecs)[rows] ) ]
   modelFrame$count <- as.vector( counts(ecs)[rows,] )
   attr( modelFrame, "geneID" ) <- geneID
   modelFrame
}

profileLogLikelihood <- function( disp, mm, y, muhat )
{
   # calculate the log likelihood:
   if(length(disp) != length(y)){
      disp <- rep(disp, length(y))
   }

   ll <- sum( sapply( 1:length(y), function(i) 
      dnbinom( y[i], mu=muhat[i], size=1/disp[i], log=TRUE ) ) )

   # transform the residuals, i.e., y - muhat, to the linear
   # predictor scale by multiplying them with the derivative
   # of the link function, i.e., by 1/muhat, and add this to the
   # linear predictors, log(muhat), to get the predictors that
   # are used in the IWLS regression
   z <- log(muhat) + ( y - muhat ) / muhat

   # the variance function of the NB is as follows
   v0 <- muhat + disp * muhat^2 

   # transform the variance vector to linear predictor scale by
   # multiplying with the squared derivative of the link to
   # get the (reciprocal) weights for the IWLS
   w <- 1 / ( ( 1 / muhat )^2 * v0 )

   # All we need from the IRLS run is the QR decomposition of
   # its matrix
   qrres <- qr( mm*sqrt(w) )

   # from it, we extract we leverages and calculate the Cox-Reid
   # term:
   cr <- sum( log( abs( diag( qrres$qr )[ 1:qrres$rank ] ) ) )  

   # return the profile log likelihood:
   ll - cr 
}

estimateExonDispersionsForModelFrame <- function( modelFrame, formula=NULL, mm=NULL, muhat=NULL, initialGuess = .01 )
{
   exonNames <- as.character(levels(modelFrame$exon))
   if( length( exonNames ) <= 1 ) {
      ans <- NA_real_
      names(ans) <- exonNames
      return( ans )
   }
   if(is.null(formula)){
	formula <- count ~ sample + condition*exon
   }
   
   if(is.null(mm)){
      mm <- model.matrix( formula, modelFrame )
   }

   if( nrow(mm) <= ncol(mm) )
      stop( "Underdetermined model; cannot estimate dispersions. Maybe replicates have not been properly specified." )

   countsums <- tapply( modelFrame$count, modelFrame$exon, sum )
   y <- modelFrame$count

   if(is.null(muhat)){
      y1 <- pmax(y, 1/6)
      fit <- lm.fit(mm, log(y1) - log(modelFrame$sizeFactor))
      mm <- mm[,!is.na(fit$coefficients)]
      start <- fit$coefficients[!is.na(fit$coefficients)]
      muhat <- fitted.values(glmnb.fit(mm, y, initialGuess, log(modelFrame$sizeFactor), start=start))
   }

   disp <- rep( initialGuess, length(exonNames) )
   names(disp) <- exonNames
      for( exon in exonNames[ countsums > 0 ] )
         disp[exon] <- exp(
            optimize( 
               function(logalpha) {
                  disp[exon] <- exp( logalpha )
                  profileLogLikelihood( disp[as.character(modelFrame$exon)], mm, y, muhat ) },
               log( c( 1e-11, 1e5 ) ),
   	       tol = 0.1,
               maximum=TRUE 
            )$maximum )
   disp[ countsums == 0 ] <- NA
   disp[ disp < 1e-10 ] <- 0
   disp
}

fitDispersionFunction <- function( ecs )
{
   stopifnot(is(ecs, "ExonCountSet"))
   if(all(is.na(fData(ecs)$dispBeforeSharing))){
      stop("no CR dispersion estimations found, please first call estimateDispersions function")
   }
   means <- colMeans( t(counts(ecs))/sizeFactors(ecs) )
   disps <- fData(ecs)$dispBeforeSharing
   coefs <- c( .1, 1 )
   iter <- 0   
   while(TRUE) {
      residuals <- disps / ( coefs[1] + coefs[2] / means )
      good <- which((residuals > 1e-4) & (residuals < 15))
      mm <- model.matrix(disps[good] ~ I(1/means[good]))
      fit <- try(glmgam.fit(mm, disps[good], start=coefs), silent=TRUE)
      if(inherits(fit, "try-error")){
         stop("Failed to fit the dispersion function\n")
      }
      oldcoefs <- coefs
      coefs <- coefficients(fit)
      if(coefs[1] < 0){
         coefs[1] <- 0
         warning("Negative intercept value in the dispersion function, it will be set to 0. Check fit diagnostics plot section from the vignette.")
         break
      }
      if( sum( log( coefs / oldcoefs )^2 ) < .005 )
         break
      iter <- iter + 1
      if( iter > 10 ) {
         warning( "Dispersion fit did not converge." )
         break }
    }
    ecs@dispFitCoefs <- coefs
    fData(ecs)$dispFitted <- ecs@dispFitCoefs[1] + ecs@dispFitCoefs[2] / colMeans( t(counts(ecs))/sizeFactors(ecs) )
    fData(ecs)$dispersion <- pmin(
       pmax( 
          fData(ecs)$dispBeforeSharing, 
          fData(ecs)$dispFitted,
          na.rm = TRUE ),
          1e8 )   # 1e8 as an arbitrary way-too-large value to capture infinities
    return(ecs)
}


setMethod("estimateDispersions", signature(object="ExonCountSet"),
   function( object, formula=count ~ sample + condition*exon, initialGuess=.01, nCores=1, minCount=10, maxExon=70, quiet=FALSE, file="")
   {
      cds <- object
      stopifnot(is(cds, "ExonCountSet"))
      if( all( is.na(sizeFactors(cds)) ) ){
         stop( "Estimate size factors first." )	
      }
      cds@formulas[["formulaDispersion"]] <- deparse(formula)


      if(!interactive() & !quiet & file == ""){  ## if the session is not interactive, quiet is FALSE and there is no file, then quiet does not make any sense
         quiet=TRUE
      }
      ########## DEFINE TESTABLE GENES #############
      # Which exons are actually testable? 
       # take away those exons with counts lower than minCounts
      testable <- rowSums(counts(cds)) > minCount
      if(!all(testable) & nCores<=1){
         warning(sprintf("Exons with less than %d counts will be discarded. For more details read the documentation, parameter minCount", minCount + 1))
      }
      # If a gene contains less than two high count exons, all its exons non-testable
      for( r in split( 1:nrow(cds), geneIDs(cds) ) ) {
         if( sum( testable[r] ) <= 1 | sum( testable[r] ) > maxExon )
            testable[r] <- FALSE
      }
      fData(cds)$testable <- testable

      # take away genes with just one exon
      generle <- rle( as.character( geneIDs(cds) ) )
      fData(cds)$testable[which(geneIDs(cds) %in% generle$values[which( generle$lengths ==1 )])] <- FALSE
      # take away those exons bigger than maxExon (default 70)
#      fData(cds)$testable[which(geneIDs(cds) %in% generle$values[which( generle$lengths > maxExon )])] <- FALSE
      ###
      
      if(max(generle$lengths) > maxExon & nCores<=1){
         warning(sprintf("Genes with more than %d testable exons will be kicked out of the analysis. For more details read the documentation, parameter maxExon", maxExon))
      }

      ##### DOES THE SAME AS A TAPPLY, but without considering the order of the levels ##
#      testable <- fData(cds)$testable
      testablegenes <- as.character(unique(fData(cds)[which(fData(cds)$testable),]$geneID))

      if(!quiet & nCores==1 ) {
         cat( "Dispersion estimation. (Progress report: one dot per 100 genes)\n", file=file, append=TRUE)
      }
         
      i <- 0
      if(nCores > 1){
         if(!quiet){
         cat(sprintf("Estimating Cox-Reid exon dispersion estimates using %d cores. (Progress report: one dot per 100 genes)\n", nCores), file=file, append=TRUE)}
         toapply <- function(x){estimateDispersions(x, formula=formula, initialGuess=initialGuess, nCores=-1, minCount=minCount, maxExon=maxExon, file=file, quiet=quiet)}
         cds <- divideWork(cds, funtoapply=toapply, fattr="dispBeforeSharing", mc.cores=nCores, testablegenes)
      }else{
         modelFrames <- lapply( testablegenes, function(gs) {
            mf <- modelFrameForGene(cds, gs, onlyTestable=TRUE)
            mf$offset <- log(mf$sizeFactor)
            mf })

         names(modelFrames) <- testablegenes
         
         modelmatrices <- lapply( modelFrames, function(mf){
            model.matrix( formula, data = mf )} )

         names(modelmatrices) <- testablegenes

         muhats <- lapply( testablegenes, function( gn ) { 
            y <- modelFrames[[gn]]$count
            y1 <- pmax(y, 1/6)
            mf <- modelFrames[[gn]]
            mm <- modelmatrices[[gn]]
            fit <- lm.fit(mm, log(y1) - mf$offset)
            mm <- mm[,!is.na(fit$coefficients)]
            if( nrow(mm) <= ncol(mm) )
               stop( "Underdetermined model; cannot estimate dispersions. Maybe replicates have not been properly specified." )
            start <- fit$coefficients[!is.na(fit$coefficients)]
            muhat <- try(fitted.values(glmnb.fit(mm, y, initialGuess, mf$offset, start=start)), silent=TRUE)
            muhat
         })

         names(muhats) <- testablegenes
   
         badones <- which( sapply( muhats, inherits, "try-error") )
         if( length(badones) > 0 ) {
            testablegenes <- testablegenes[ ! testablegenes %in% names(badones) ]
    	      warning( paste( "Failed to set up model frames for genes ", 
   	         paste( names(badones), collapse=", " ) ) )
         }

         ###### WORKS FASTER THIS WAY THAN IN EVERY CYCLE ACCESSING fData
         fData(cds)$dispBeforeSharing <- NA_real_
         dispBeforeSharing <- fData(cds)$dispBeforeSharing
         testable <- fData(cds)$testable
         exonids <- fData(cds)$exonID

         for( genename in testablegenes ){
            mf <- modelFrames[[genename]]
            disps <- try( 
               estimateExonDispersionsForModelFrame( mf, 
                  mm=modelmatrices[[genename]], muhat=muhats[[genename]] ),
               silent=TRUE )
            if( inherits( disps, "try-error" ) ) {
               disps <- rep( NA_real_, length( muhats[[genename]] ) )
               warning( sprintf( "Failed to fit dispersion for gene %s", genename ) )
            }   
      
            rows <- as.character(geneIDs(cds)) %in% genename & testable
            stopifnot(all(names(disps)==exonids[rows]))
            dispBeforeSharing[rows] <- disps 
      
            i <- i + 1
            if(!quiet & i %% 100 == 0 ){
               cat( ".", file=file, append=TRUE)}    
         }
         fData(cds)$dispBeforeSharing <- dispBeforeSharing
      }
      if(nCores==1 & !quiet){
         cat("\n", file=file, append=TRUE)
      }
      cds
   }
)

testGeneForDEU <- function (ecs, gene, formula0=NULL, formula1=NULL ){
   stopifnot(is(ecs, "ExonCountSet"))
   if( all( is.na(featureData(ecs)$dispersion ) ) ) {
      stop("No dispersion values found, call function fitDispersionFunction first.")		
   }

   mf <- modelFrameForGene(ecs, gene, onlyTestable=TRUE)
   
   ans <- data.frame( row.names = levels( mf$exon ) )
   ans$deviance = NA_real_ 
   ans$df = NA_integer_
   ans$pvalue = NA_real_

   if(is.null(formula0)){
      formula0 <- count ~ sample + exon + condition}
   if(is.null(formula1)){
      formula1 <- count ~ sample + exon + condition * I(exon == exonID)}
   

   # This makes sure that the formula see the 'exonID' variable used
   # below even if it the formula was supplied as parameter
   environment(formula0) <- environment()
   environment(formula1) <- environment()
   
   mm <- model.matrix(formula0, mf)
   y1 <- pmax(mf$count, 1/6)
   fit <- lm.fit(mm, log(y1) - log(mf$sizeFactor))
   mm <- mm[,!is.na(fit$coefficients)]
   start <- fit$coefficients[!is.na(fit$coefficients)]
   
   fit0 <- try(
      glmnb.fit(mm, mf$count, mf$dispersion, log(mf$sizeFactor), start=start),
   silent=TRUE)
   if( inherits( fit0, "try-error" ) ) {
      warning( sprintf( "Error in fit0 for gene %s: %s", gene, fit0 ) )
      return(ans) }

   for( exonID in unique(as.character(mf$exon)) ) {
      mm <- model.matrix(formula1, mf)
      y1 <- pmax(mf$count, 1/6)
      fit <- lm.fit(mm, log(y1) - log(mf$sizeFactor))
      mm <- mm[,!is.na(fit$coefficients)]
      start <- fit$coefficients[!is.na(fit$coefficients)]

      fit1 <- try(
         glmnb.fit(mm, mf$count, mf$dispersion, log(mf$sizeFactor), start=start),
      silent=TRUE)
      if( inherits( fit1, "try-error" ) ) {
         warning( sprintf( "Error in fit1 for gene %s, exon %s: %s", gene, exonID, fit1 ) )
         next }
      ans[ exonID, "deviance" ] <- deviance( fit0 ) - deviance( fit1 )
      ans[ exonID, "df" ] <- length(fit1$coefficients) - length(fit0$coefficients)
      ans[ exonID, "pvalue"   ] <- 1 - pchisq( ans[exonID,"deviance"], ans[exonID,"df"] )
  }
  ans
}


testForDEU <- function( ecs, formula0=NULL, formula1=NULL, nCores=1, quiet=FALSE, file="")
{
   stopifnot(is(ecs, "ExonCountSet"))
   if( all( is.na(featureData(ecs)$dispersion ) ) ) {
      stop("No dispersion values found, call function fitDispersionFunction first.")		
   }
   
   if(!interactive() & !quiet & file == ""){  ## if the session is not interactive, quiet is FALSE and there is no file, then quiet does not make any sense
      quiet=TRUE
   }

   if(!quiet & nCores==1){
      cat( "Testing for differential exon usage. (Progress report: one dot per 100 genes)\n", file=file, append=TRUE)}

   testablegenes <- as.character(unique(fData(ecs)[which(fData(ecs)$testable),]$geneID))

   if(nCores > 1){
      if(!quiet){
      cat(sprintf("Testing for differential exon usage using %d cores. (Progress report: one dot per 100 genes)\n", nCores), file=file, append=TRUE)}
      toapply <- function(x){testForDEU(x, formula0=formula0, formula1=formula1, nCores=-1, file=file, quiet=quiet)}
      ecs <- divideWork(ecs, funtoapply=toapply, fattr="pvalue", mc.cores=nCores, testablegenes)
   }else{
      i <- 0
      testable <- fData(ecs)$testable
      pvalue <- fData(ecs)$pvalue
      exonids <- fData(ecs)$exonID
      for( genename in testablegenes ) {
         i <- i + 1
         if(!quiet & i %% 100 == 0 ){
            cat( ".", file=file, append=TRUE)}
         rows <- as.character(geneIDs(ecs)) %in% genename & testable
         res <- testGeneForDEU( ecs, genename, formula0, formula1 )
         stopifnot(all(rownames(res) == exonids[rows]))    ### makes sure to do the assignment correctly
         pvalue[rows] <- res$pvalue 
      }
      fData(ecs)$pvalue <- pvalue
   }
   if(nCores==1 & !quiet){
      cat("\n", file=file, append=TRUE)
   }
   fData(ecs)$padjust <- p.adjust(fData(ecs)$pvalue, method="BH")
   ######### STORE FORMULAS IN THE OBJECTs
   if(is.null(formula0)){
      formula0 <- count~sample+exon+condition
   }
   if(is.null(formula1)){
      formula1 <- count~sample+exon+condition*I(exon==exonID)
   }
   ecs@formulas[["formula0"]] <- deparse(formula0)
   ecs@formulas[["formula1"]] <- deparse(formula1)
   ecs
}

estimatelog2FoldChanges <- function(ecs, fitExpToVar="condition", nCores=1, quiet=FALSE, file="")
{
   stopifnot(is(ecs, "ExonCountSet"))
   if(any(is.na(sizeFactors(ecs)))){
      stop("Please estimate sizeFactors first\n")}
   if(!fitExpToVar %in% ecs@designColumns){
      stop("fitExpToVar parameter is not in the design columns, double check ecs@designColumns")}
   if(sum(is.na(featureData(ecs)$dispersion))==nrow(counts(ecs))){
      stop("No dispersion parameters found, first call function estimateDispersions...\n")}

   nms <- sort(levels(design(ecs, drop=FALSE)[[fitExpToVar]]))
   colfoldnames <- sapply(2:length(nms), function(x){
      paste("log2fold(", nms[x], "/", nms[1], ")", sep="")
   })

   if(!interactive() & !quiet & file == ""){  ## if the session is not interactive, quiet is FALSE and there is no file, then quiet does not make any sense
      quiet=TRUE
   }
 
   testablegenes <- as.character(unique(fData(ecs)[which(fData(ecs)$testable),]$geneID))
  
   logfold <- data.frame(matrix(ncol=length(colfoldnames), nrow=nrow(fData(ecs))))
   colnames(logfold) <- colfoldnames
   rownames(logfold) <- rownames(fData(ecs))

   if(!any(colfoldnames %in% colnames(fData(ecs)))){
      fData(ecs) <- cbind(fData(ecs), logfold)
   }

   if(!quiet & nCores==1){
      cat("Calculating fold changes. (Progress report: one dot per 100 genes)\n", file=file,append=TRUE)}
 
   if (nCores > 1){
      if(!quiet){
      cat(sprintf("Calculating fold changes using %d cores. (Progress report: one dot per 100 genes)\n", nCores), file=file, append=TRUE)}
      toapply <- function(x){estimatelog2FoldChanges(x, fitExpToVar=fitExpToVar, nCores=-1, file=file, quiet=quiet)}
      ecs <- divideWork(ecs, toapply, fattr=colfoldnames, mc.cores=nCores, testablegenes)
   }else{
      frm <- as.formula(paste("count ~", fitExpToVar,  "* exon"))
      colstosteal <- paste(fitExpToVar, ":exon", sep="")
      colstopass <- c(2:length(nms))
      i <- 0
      for(gene in testablegenes){
         i <- i + 1
         if(!quiet & i %% 100 == 0 ){
            cat( ".", file=file, append=TRUE)}
         tr <- fitAndArrangeCoefs( ecs, gene, frm=frm)
         if(is.null(tr)){
              warning(sprintf("log fold change calculation failed for gene %s", gene))
              next
         }
         rows <- as.character(geneIDs(ecs)) %in% gene
         tr <- tr[[colstosteal]]
         stopifnot(all(fData(ecs)$exonID[rows] == colnames(tr))) ##### makes sure to do good the assignation
         vals <- log2(exp(t(tr)))[,colstopass]
         logfold[rows, colfoldnames] <- vals
      }
      fData(ecs)[,colfoldnames] <- logfold
   }
   if(nCores==1 & !quiet){
      cat("\n", file=file, append=TRUE)
   }
   ecs
}

divideWork <- function(ecs, funtoapply, fattr, mc.cores, testablegenes)
{
#   if(!suppressMessages(suppressWarnings(require("multicore")))){
 #     stop("multicore package not found...")}
   if(!is.loaded("mc_fork", PACKAGE="multicore")){
     stop("Please load first multicore package or set parameter nCores to 1...")}
   stopifnot(mc.cores>=1)

   subgenes <- split(testablegenes, seq(along=testablegenes) %% mc.cores)
   allecs <- lapply(subgenes, function(x) subsetByGenes(ecs, x) )
   allecs <- multicore::mclapply(allecs, FUN=funtoapply, mc.cores=mc.cores)

   for(j in seq(along=allecs)) {
      rownam <-  rownames(fData(allecs[[j]]))
      fData(ecs)[ rownam, fattr] <-  fData(allecs[[j]])[, fattr]
   }
  
   ecs
}

makeCompleteDEUAnalysis <- function(ecs, formulaDispersion=count ~ sample + condition*exon, minCount=10, maxExon=50, formula0=NULL, formula1=NULL, FDR=0.1, fitExpToVar="condition", nCores=1, path=NULL, color=NULL, color.samples=NULL, quiet=FALSE, file="")
{
   stopifnot(is(ecs, "ExonCountSet"))
   ecs <- estimateSizeFactors( ecs )
   ecs <- estimateDispersions( ecs, formulaDispersion, nCores=nCores, quiet=quiet, file=file, minCount=minCount, maxExon=maxExon)
   ecs <- fitDispersionFunction( ecs )
   ecs <- testForDEU( ecs, formula1=formula1, formula0=formula0, nCores=nCores, quiet=quiet, file=file)
   ecs <- estimatelog2FoldChanges(ecs, fitExpToVar=fitExpToVar, nCores=nCores, quiet=quiet, file=file)
   if(!is.null(path)){
      DEXSeqHTML(ecs, path=path, FDR=0.1, fitExpToVar=fitExpToVar, color=color, color.samples=color.samples)
   }
   ecs
}
