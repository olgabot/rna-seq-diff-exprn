estimateSizeFactorsForMatrix <- function( counts, locfunc = median )
{
   loggeomeans <- rowMeans( log(counts) ) 
   apply( counts, 2, function(cnts) 
      exp( locfunc( ( log(cnts) - loggeomeans )[ is.finite(loggeomeans) ] ) ) )
}


getBaseMeansAndVariances <- function( counts, sizeFactors ) {

   # Divides the counts by sizeFactors and calculates the estimates for
   # base means and variances for each gene.
   
   data.frame(
      baseMean = rowMeans( t( t(counts) / sizeFactors ) ),
      baseVar = rowVars( t( t(counts) / sizeFactors ) ) )
}   

modelMatrixToConditionFactor <- function( modelMatrix ) {

   mmconds <- 1:nrow(modelMatrix)
   for( i in 2:nrow(modelMatrix) )
      for( j in 1:(i-1) )
         if( all( modelMatrix[i,] == modelMatrix[j,] ) ) {
            mmconds[i] = mmconds[j]
            break }
   factor( as.integer( factor( mmconds ) ) )
}


getBaseMeansAndPooledVariances <- function( counts, sizeFactors, conditions ) {

   basecounts <- t( t(counts) / sizeFactors )
   replicated_sample <- conditions %in% names(which(table(conditions)>1))
   df <- sum(replicated_sample) - length( unique( conditions[ replicated_sample ] ) ) 

   data.frame(
      baseMean = rowMeans( basecounts ),
      baseVar =
	 rowSums( 
	    sapply( 
               tapply( 
        	  ( 1:ncol(counts) )[ replicated_sample ], 
        	  factor( conditions[ replicated_sample ] ), 
        	  function(cols) 
        	     rowSums( ( basecounts[,cols] - rowMeans(basecounts[,cols]) )^2 ) ), 
               identity ) ) / df )
}

parametricDispersionFit <- function( means, disps ) 
{
   coefs <- c( .1, 1 )
   iter <- 0   
   while(TRUE) {
      residuals <- disps / ( coefs[1] + coefs[2] / means )
      good <- which( (residuals > 1e-4) & (residuals < 15) )
      fit <- glm( disps[good] ~ I(1/means[good]), 
         family=Gamma(link="identity"), start=coefs )
      oldcoefs <- coefs   
      coefs <- coefficients(fit)
      if( !all( coefs > 0 ) )
         stop( "Parametric dispersion fit failed. Try a local fit and/or a pooled estimation. (See '?estimateDispersions')" )
      if( sum( log( coefs / oldcoefs )^2 ) < 1e-6 )
         break
      iter <- iter + 1
      if( iter > 10 ) {
         warning( "Dispersion fit did not converge." )
         break }
   }

   names( coefs ) <- c( "asymptDisp", "extraPois" )
   ans <- function( q )
      coefs[1] + coefs[2] / q
   attr( ans, "coefficients" ) <- coefs
   ans
}   

estimateAndFitDispersionsFromBaseMeansAndVariances <- function( means, 
   variances, sizeFactors, fitType = c( "parametric", "local" ),
   locfit_extra_args=list(), lp_extra_args=list(), adjustForBias=TRUE ) {
   
   fitType <- match.arg( fitType )

   xim <- mean( 1/sizeFactors )
   dispsAll <- ( variances - xim * means ) / means^2
   
   variances <- variances[ means > 0 ]
   disps <- dispsAll[ means > 0 ]
   means <- means[ means > 0 ]

   if( adjustForBias )
      disps <- adjustScvForBias( disps, length( sizeFactors ) )
   
   if( fitType == "local" ) {
   
      fit <- do.call( "locfit", c( 
         list( 
            variances ~ do.call( "lp", c( list( log(means) ), lp_extra_args ) ),
            family = "gamma" ), 
         locfit_extra_args ) )
      
      rm( means )
      rm( variances )
      
      if( adjustForBias )
         ans <- function( q )
            adjustScvForBias( 
               pmax( ( safepredict( fit, log(q) ) - xim * q ) / q^2, 1e-8 ),
               length(sizeFactors) )
      else
         ans <- function( q )
            pmax( ( safepredict( fit, log(q) ) - xim * q ) / q^2, 1e-8 )
            
      # Note: The 'pmax' construct above serves to limit the overdispersion to a minimum
      # of 10^-8, which should be indistinguishable from 0 but ensures numerical stability.

   } else if( fitType == "parametric" ) {

      ans <- parametricDispersionFit( means, disps )
   
   } else
      stop( "Unknown fitType." )
   
   attr( ans, "fitType" ) <- fitType
   list( disps=dispsAll, dispFunc=ans )
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

     
estimateAndFitDispersionsWithCoxReid <- function( counts, modelFormula, modelFrame,
   sizeFactors, fitType = c( "parametric", "local" ),
   locfit_extra_args=list(), lp_extra_args=list(), initialGuess=.1 ) 
{
   if( as.character(modelFormula[2]) == "count" )
      modelFormula <- modelFormula[-2]   # the '[-2]' removes the lhs, i.e., the response
   mm <- model.matrix( modelFormula, modelFrame )  
   disps <- apply( counts, 1, function( y ) {
      fit <- glm.fit( mm, y, family=MASS::negative.binomial( initialGuess ), offset=log(sizeFactors) )
      if( df.residual(fit) == 0 )
         stop( "No residual degrees of freedom. Most likely the design is lacking sufficient replication." )
      exp(
         optimize( 
            function(logalpha)
               profileLogLikelihood( exp(logalpha), mm, y, fitted.values(fit) ),
            log( c( 1e-11, 1e5 ) ),
            maximum=TRUE 
         )$maximum ) } )

   means <- colMeans( t(counts) / sizeFactors )
   xim <- mean( 1/sizeFactors )

   if( fitType == "local" ) {
         
      fit <- do.call( "locfit", c( 
         list( 
            disps[means>0] ~ do.call( "lp", c( list( log(means[means>0]) ), lp_extra_args ) ),
            family = "gamma" ), 
         locfit_extra_args ) )
      
      rm( means )
      
      ans <- function( q )
         pmax( ( safepredict( fit, log(q) ) - xim * q ) / q^2, 1e-8 )
            
      # Note: The 'pmax' construct above serves to limit the overdispersion to a minimum
      # of 10^-8, which should be indistinguishable from 0 but ensures numerical stability.

   } else if( fitType == "parametric" )

      ans <- parametricDispersionFit( means, disps )
   
   else
      stop( "Unkknown fitType." )
   
   attr( ans, "fitType" ) <- fitType
   list( disps=disps, dispFunc=ans )   
   
}      
   
safepredict <- function( fit, x )
{
   # A wrapper around predict to avoid the issue that predict.locfit cannot
   # propagate NAs and NaNs properly.

   res <- rep.int( NA_real_, length(x) )
   res[ is.finite(x) ] <- predict( fit, x[is.finite(x)] )
   res   
}

nbinomTestForMatrices <- function( countsA, countsB, sizeFactorsA, sizeFactorsB, 
   dispsA, dispsB )
{

   kAs <- rowSums( cbind(countsA) )
   kBs <- rowSums( cbind(countsB) )
   
   mus <- rowMeans( cbind(      
      t( t( countsA ) / sizeFactorsA ),
      t( t( countsB ) / sizeFactorsB ) ) )      

   fullVarsA <- pmax( mus * sum( sizeFactorsA ) + dispsA * mus^2 * sum(sizeFactorsA^2), 
      mus * sum( sizeFactorsA ) * (1+1e-8) )
   fullVarsB <- pmax( mus * sum( sizeFactorsB ) + dispsB * mus^2 * sum(sizeFactorsB^2), 
      mus * sum( sizeFactorsB ) * (1+1e-8) )
   
   sumDispsA <- ( fullVarsA - mus * sum( sizeFactorsA ) ) / ( mus * sum( sizeFactorsA ) )^2
   sumDispsB <- ( fullVarsB - mus * sum( sizeFactorsB ) ) / ( mus * sum( sizeFactorsB ) )^2

   sapply( 1:length(kAs), function(i) {
   
      if( kAs[i] == 0 & kBs[i] == 0 )
         return( NA )
      
      # probability of all possible counts sums with the same total count:
      ks <- 0 : ( kAs[i] + kBs[i] )
      ps <- dnbinom(                   ks, mu = mus[i] * sum( sizeFactorsA ), size = 1/sumDispsA[i] ) * 
            dnbinom( kAs[i] + kBs[i] - ks, mu = mus[i] * sum( sizeFactorsB ), size = 1/sumDispsB[i] )
            
      # probability of observed count sums:
      pobs <- dnbinom( kAs[i], mu = mus[i] * sum( sizeFactorsA ), size = 1/sumDispsA[i] ) * 
              dnbinom( kBs[i], mu = mus[i] * sum( sizeFactorsB ), size = 1/sumDispsB[i] )
              
      stopifnot( pobs == ps[ kAs[i]+1 ] )
      if( kAs[i] * sum( sizeFactorsB ) < kBs[i] * sum( sizeFactorsA ) )
         numer <- ps[ 1 : (kAs[i]+1) ]
      else 
         numer <- ps[ (kAs[i]+1) : length(ps) ]
      min( 1, 2 * sum(numer) / sum(ps) )
   } )
}


# Note: The following function is never called; it is here only for
# documentation purposes, as it has been used to produce the data object
# scvBiasCorrectionFits, which is stored in the file
# inst/scvBiasCorrectionFits.rda, gets loadewd by the line after this
# function and is used by the function adjustScvForBias
#
# To do: the correct place for this type of documentation is a vignette
#

prepareScvBiasCorrectionFits <- function( maxnrepl=15, mu=100000, ngenes=10000,
      true_raw_scv = c( seq( 0, 2, length.out=100 )[-1], seq( 2, 10, length.out=20 )[-1] ) )
   lapply( 2:maxnrepl, function( m ) {
      est_raw_scv <- sapply( true_raw_scv, function( alpha ) {
         k <- matrix( rnbinom( ngenes*m, mu=mu, size=1/alpha ), ncol=m )
         k <- k[ rowSums(k)>0, ]
         mean( rowVars(k) / rowMeans(k)^2 ) } )
      locfit( true_raw_scv ~ lp( est_raw_scv, nn=.2 ) ) } )

load( system.file ( "extra/scvBiasCorrectionFits.rda", package="DESeq" ) )

adjustScvForBias <- function( scv, nsamples ) {
   stopifnot( nsamples > 1 )
   if( nsamples - 1 > length( scvBiasCorrectionFits ) )
      scv
   else
      ifelse( scv > .02,
         pmax( safepredict( scvBiasCorrectionFits[[ nsamples-1 ]], scv ), 1e-8 * scv ),
         scv )   # For scv < .02, our fit is too coarse, but no correction seems necessary anyway
}      

nbkd.sf <- function( r, sf ) {
   fam <- list(
     
      family = sprintf( "nbkd,r=%g", r ),
      link = "log_sf",
      linkfun = function( mu ) log( mu/sf ),
      linkinv = function (eta) pmax(sf*exp(eta), .Machine$double.eps),
      mu.eta = function (eta) pmax(sf*exp(eta), .Machine$double.eps),
      variance = function(mu) mu + mu^2 / r,

      dev.resids = function( y, mu, wt )
          2 * wt * ( ifelse( y > 0, y * log( y / mu ), 0 ) + 
            (r + y) * log( (r+mu)/(r+y) ) ), 
      
      aic = function (y, n, mu, wt, dev) NA,   # what for?
      initialize = expression( {
         n <- rep.int(1, nobs)         # What is n?
         mustart <- y + 0.1
      } ),
      valid.mu <- function(mu) all( mu > 0 ),
      valid.eta <- function(eta) TRUE,
      simulate <- NA
   )
   
   class(fam) <- "family"
   fam }        


fitNbinomGLMsForMatrix <- function( counts, sizeFactors, rawScv, modelFormula, 
   modelFrame, quiet=FALSE, reportLog2=TRUE, glmControl=list() ) 
{
   stopifnot( length(sizeFactors) == ncol(counts) )
   stopifnot( length(rawScv) == nrow(counts) )
   stopifnot( nrow(modelFrame) == ncol(counts) )

   stopifnot( is( modelFormula, "formula" ) )  
   if( as.character( modelFormula[[1]] ) != "~" )  
      stop( "Formula does not have a '~' as top-level operator." )
   if( as.character( modelFormula[[2]] ) != "count" )  
      stop( "Left-hand side of model formula must be 'count'." )
   
   goodRows <- is.finite( rawScv ) & rowSums(counts) > 0 
   
   modelMatrix <- model.matrix( modelFormula[c(1,3)], modelFrame )
   res <- 
   t( sapply( which(goodRows), function(i) {
      if( !quiet & i %% 1000 == 0 )
         cat( '.' ) 
      nbfam <- nbkd.sf( 1 / rawScv[i], sizeFactors )      
      fit <- try( 
         glm.fit( modelMatrix, counts[i,], family=nbfam, control = glmControl ),
         silent=TRUE )
      if( !is( fit, "try-error" ) )
         c( 
            coefficients(fit), 
            deviance = deviance(fit), 
            df.residual = fit$df.residual,
            converged = fit$converged ) 
      else {
         coefs <- rep( NA, ncol(modelMatrix ) )
         names(coefs)  <- colnames(modelMatrix) 
         #warning( as.character( fit ) )
         c( 
            coefs, 
            deviance = NA, 
            df.residual = NA,
            converged = FALSE ) } } ) )
      
   if( !quiet )
      cat( "\n" ) 

   df.residual <- max( na.omit( res[ , "df.residual" ] ) )
   if( !all( na.omit( res[ , "df.residual" ] == df.residual ) ) ) {
      res[ res[ , "df.residual" ] != df.residual, "deviance" ] <- NA
      warning( "Some deviances set to NA due to reduction in degrees of freedom." )
   }

   # Put in the NAs
   res2 <- data.frame(
      row.names = row.names( counts ),
      apply( res, 2, function( col ) {
         a <- rep( NA_real_, nrow(counts) )
         a[goodRows] <- col
         a } ) )
   colnames(res2) <- colnames(res)
      
   if( reportLog2 )
      res2[ , 1:(ncol(res2)-3) ] <- res2[ , 1:(ncol(res2)-3) ] / log(2)

   res2$converged <- as.logical( res2$converged )

   res2 <- res2[ , colnames(res) != "df.residual" ]     
   attr( res2, "df.residual" ) <- df.residual
   res2
}
      
